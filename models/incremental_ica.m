function [A,YS]=incremental_ica(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('wdim',NaN,@isnumeric);
pr.addParamValue('alpha',0.1,@isnumeric);
pr.addParamValue('beta',0.5,@isnumeric);
pr.addParamValue('minibatch_size',100,@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('eigen',{});
pr.parse(varargin{:});
options=pr.Results;

[datadim,ndata]=size(dataset);

% lowpass & whitening

if isempty(options.eigen)
    [wfilters,~,evals]=pca(dataset');
else
    wfilters=options.eigen{1};
    evals=options.eigen{2};
end;

if isnan(options.wdim)
    wdim=nunits;
else
    wdim=options.wdim;
end;
K=evals(1:wdim).^(-1/2);
wfilters=wfilters(:,1:wdim)*diag(K);
wdata=wfilters'*dataset;

if ~isempty(options.params)
    A=options.params;
else
    A=randn(wdim,nunits)*0.1;
    A=A./repmat(sqrt(sum(A.^2,1)),wdim,1);
end;

%%%%%

nsteps=10000;
minibatch_size=options.minibatch_size;
display_steps=1000;
alpha=options.alpha;
beta=options.beta;

dA=zeros(wdim,nunits);

YS=zeros(nunits,nsteps);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

f=@(x)log(cosh(x));
g=@(x)tanh(x);
gd=@(x)1-tanh(x).^2;

% f=@(u)-exp(-u.^2/2);
% g=@(u)u.*exp(-u.^2/2);
% gd=@(u)(1-u.^2).*exp(-u.^2/2);

for T=1:nsteps
    if mod(T,minibatch_size)==1 
        A=A+dA/minibatch_size;
        A=A+beta*(A-A*A'*A);
%         A=A*(A'*A)^(-1/2);
        A=A./repmat(sqrt(sum(A.^2,1)),wdim,1);
        dA(:)=0;
    end;
    if mod(T,display_steps)==1 && options.display 
        fprintf('step #%d\n',T);
        disp_weight(h,A,wfilters); 
    end;
%     if mod(T,nsteps/8)==0
%         alpha=alpha/2;
%         fprintf('alpha=%e\n',alpha);
%     end;

    % inference

    X=wdata(:,mod(T-1,ndata)+1);
    Y=X'*A;
    
    dA=dA-alpha*X*g(X'*A);

    YS(:,T)=Y;
    
    
end;

end


function disp_weight(h,A,wfilters)
    figure(h); clf;
    subplot(1,1,1);
    W=pinv(wfilters)'*A;
    W=W./repmat(sqrt(sum(W.^2,1)),size(W,1),1);
    [datadim,nunits]=size(W);
    pixwid=sqrt(datadim);
    wid=ceil(sqrt(nunits)); ht=ceil(nunits/wid);

    M=ones(wid*(pixwid+1),ht*(pixwid+1))*min(W(:));
    e=zeros(nunits,1);
    for I=1:nunits
        [x y]=ind2sub([wid ht],I);
        x=(x-1)*(pixwid+1)+1; y=(y-1)*(pixwid+1)+1;
        v=W(:,I);
        M(x:x+pixwid-1,y:y+pixwid-1)=reshape(v,pixwid,pixwid);
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;


end
                  
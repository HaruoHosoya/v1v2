function [A,W,F,R,YS]=incremental_smica(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('learnscale',true,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('wdim',NaN,@isnumeric);
pr.addParamValue('alpha',0.1,@isnumeric);
pr.addParamValue('minibatch_size',100,@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('eigen',{});
pr.addParamValue('initscale',1.5,@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

[~,ndata]=size(dataset);

% lowpass & whitening

if isempty(options.eigen)
    [E,~,D]=pca(dataset');
else
    E=options.eigen{1};
    D=options.eigen{2};
end;

if isnan(options.wdim)
    wdim=nunits;
else
    wdim=options.wdim;
end;
F=E(:,1:wdim)*diag(D(1:wdim).^(-1/2));
wdata=F'*dataset;

if ~isempty(options.params)
    R=options.params;
else
    R=randn(wdim,nunits)*0.1;
    R=bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
end;

%%%%%

nsteps=10000;
minibatch_size=options.minibatch_size;
display_steps=1000;
alpha=options.alpha;

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

G=@(x)-log(cosh(x));
g=@(x)-tanh(x);
gd=@(x)-1+tanh(x).^2;
gdd=@(x)2*tanh(x).*(1-tanh(x).^2);

dR=zeros(wdim,nunits);
dnu=zeros(nunits,1);
nu=ones(nunits,1)*options.initscale;
YS=zeros(nunits,nsteps);
J=0;

for T=1:nsteps
    if mod(T,minibatch_size)==1 
        R=R-alpha*dR/minibatch_size;
        if options.learnscale nu=nu-alpha*dnu/minibatch_size; end;
        R=bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
        dR(:)=0; dnu(:)=0;
        J=0;
    end;
    if mod(T,display_steps)==1 && options.display 
        fprintf('step #%d\n',T);
        disp_weight(h,R,F,nu); 
    end;

    % inference

    X=wdata(:,mod(T-1,ndata)+1);
    Y=R'*X;
    Y0=g(Y);
    Y1=gd(Y);
    Y2=gdd(Y);
    Z=R*(nu.*Y0);
    U=R'*Z;
    
    dR=dR+X*(nu.*Y2)'+Z*(nu.*Y0)'+X*(nu.*Y1.*U)';
    dnu=dnu+Y1+Y0.*U;
    
    J=J+sum(nu.*Y1)+1/2*sum(Z.^2);
    YS(:,T)=Y;
    
    
end;

W=F*R;
A=pinv(F)'*pinv(R)';

end


function disp_weight(h,R,F,nu)
    figure(h); clf;
    subplot(1,2,1);
    A=pinv(F)'*pinv(R)';
%     A=F*R;
    A=bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));
    [datadim,nunits]=size(A);
    pixwid=sqrt(datadim);
    wid=ceil(sqrt(nunits)); ht=ceil(nunits/wid);

    M=ones(wid*(pixwid+1),ht*(pixwid+1))*min(A(:));
    e=zeros(nunits,1);
    for I=1:nunits
        [x y]=ind2sub([wid ht],I);
        x=(x-1)*(pixwid+1)+1; y=(y-1)*(pixwid+1)+1;
        v=A(:,I);
        M(x:x+pixwid-1,y:y+pixwid-1)=reshape(v,pixwid,pixwid);
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(1,2,2);
    barh(nu);

end
                  
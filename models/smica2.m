function [A,W,F,R,YS]=smica2(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('wdim',NaN,@isnumeric);
pr.addParamValue('alpha',0.1,@isnumeric);
pr.addParamValue('minibatch_size',100,@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('eigen',{});
pr.addParamValue('initscale',1.5,@isnumeric);
pr.addParamValue('learnscale',true,@islogical);
pr.addParamValue('nowhitening',false,@islogical);
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
if options.nowhitening
    F=E(:,1:wdim);
else
    F=E(:,1:wdim)*diag(D(1:wdim).^(-1/2));
end;
wdata=F'*dataset;

if ~isempty(options.params)
    R=options.params;
else
    R=randn(wdim,nunits)*0.1;
    R=bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
end;

%%%%%

nsteps=10000;
alpha=options.alpha;

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

G=@(x)-log(cosh(x));
g=@(x)-tanh(x);
gd=@(x)-1+tanh(x).^2;
gdd=@(x)2*tanh(x).*(1-tanh(x).^2);

% G=@(x)-sqrt(1+x.^2);
% g=@(x)-x./sqrt(1+x.^2);
% gd=@(x)-1./sqrt(1+x.^2).^3;
% gdd=@(x)3/2*x./sqrt(1+x.^2).^5;

YS=zeros(nunits,nsteps);
nu=ones(nunits,1)*options.initscale;

for T=1:nsteps
    dR=zeros(wdim,nunits);
    dnu=zeros(nunits,1);

    Y=R'*wdata;
    Y0=g(Y);
    Y1=gd(Y);
    Y2=gdd(Y);
    Z=R*bsxfun(@times,Y0,nu);
    U=R'*Z;

    dR=dR+wdata*bsxfun(@times,Y2,nu)'+Z*bsxfun(@times,Y0,nu)'+wdata*bsxfun(@times,Y1.*U,nu)';
    dnu=dnu+sum(Y1,2)+sum(Y0.*U,2);

    J=sum(sum(bsxfun(@times,Y1,nu)))+1/2*sum(sum(Z.^2));
    
    Rold=R;
    nuold=nu;
    Jold=J;
    alpha1=alpha;    
    for L=1:10
        R=R-alpha1*dR/ndata;
        R=bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
        if options.learnscale nu=nu-alpha1*dnu/ndata; end;
    
        Y=R'*wdata;
        Y0=g(Y);
        Y1=gd(Y);
        Z=R*bsxfun(@times,Y0,nu);
        J=sum(sum(bsxfun(@times,Y1,nu)))+1/2*sum(sum(Z.^2));

        if J<Jold break; end;
        R=Rold; J=Jold; nu=nuold;
        alpha1=alpha1/2;
    end;

    fprintf('step #%d  J=%e  alpha=%f\n',T,J/ndata,alpha1);

    if options.display disp_weight(h,R,F,nu); end;
        
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
                  
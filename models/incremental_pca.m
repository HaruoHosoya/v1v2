function [E,D,YS]=incremental_pca(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('alpha',0.5,@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('candid',false,@islogical);
pr.parse(varargin{:});
options=pr.Results;

[datadim,ndata]=size(dataset);

E=randn(datadim,nunits)*0.1;
D=randn(nunits,1)*0.1;
Xm=randn(datadim,1)*0.1;

if ~isempty(options.params)
    E=options.params{1};
    D=options.params{2};
    Xm=options.params{3};
end;

%%%%%

nsteps=10000;
minibatch_size=500;
alpha=options.alpha;

dE=zeros(datadim,nunits);
dD=zeros(nunits,1);
dXm=zeros(datadim,1);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

YS=zeros(nunits,nsteps);
LS=zeros(1,nsteps/minibatch_size);
LS_tmp=zeros(1,minibatch_size);

for T=1:nsteps
    if mod(T,minibatch_size)==1 
        LS(1,ceil(T/minibatch_size))=mean(LS_tmp);
        fprintf('step #%d (exp.var.=%3.1f)\n',T,(1-mean(LS_tmp))*100);
        E=E+alpha*dE/minibatch_size;
        
        if ~options.candid
            % oja & karuhnen
            E=orthonorm(E);
            D=D+alpha*dD/minibatch_size;
            Xm=Xm+alpha*dXm/minibatch_size;
        else
            % candid
            D=sqrt(sum(E.^2,1));            
        end;

        dE(:)=0;
        dD(:)=0;
        dXm(:)=0;
        if options.display disp_weight(h,E,D,LS(1,1:ceil(T/minibatch_size))); end;
    end;
%     if mod(T,nsteps/4)==0
%         alpha=alpha/2;
%         fprintf('alpha=%e\n',alpha);
%     end;

    % inference

    X=dataset(:,mod(T-1,ndata)+1);
%     X0=X;
    X0=X-Xm;
    
    if options.candid
        % candid
        for I=1:nunits
            a=E(:,I)/norm(E(:,I));
            dE(:,I)=dE(:,I)+X*(X'*a)-E(:,I);
            X=X-(X'*a)*a;
        end;
    else
        % oja & karuhnen
        dE=dE+X0*(X0'*E);
        dD=dD+(E'*X0).^2-D;    
        dXm=dXm+X0;
    end;
    
    %
    
    A0=E./repmat(sqrt(sum(E.^2,1)),size(E,1),1);
    Y=A0'*X0;
    YS(:,T)=Y;

    L=sum((X0-A0*Y).^2)/sum(X0.^2);
    LS_tmp(1,mod(T-1,minibatch_size)+1)=L;
    
end;

E=E./repmat(sqrt(sum(E.^2,1)),size(E,1),1);

end

function B=orthonorm(A)
    B=A; 
    for I=1:size(A,2)
       for J=1:I-1
           B(:,I)=B(:,I)-B(:,J)'*A(:,I)*B(:,J);
       end
       B(:,I)=B(:,I)/norm(B(:,I));
    end
end

function disp_weight(h,E,D,LS)
    figure(h); clf;
    subplot(2,3,[1 2 4 5]);
    [datadim,nunits]=size(E);
    pixwid=sqrt(datadim);
    wid=ceil(sqrt(nunits)); ht=ceil(nunits/wid);

    E=E./repmat(sqrt(sum(E.^2,1)),size(E,1),1);
    M=ones(wid*(pixwid+1),ht*(pixwid+1))*min(E(:));
    e=zeros(nunits,1);
    for I=1:nunits
        [x y]=ind2sub([wid ht],I);
        x=(x-1)*(pixwid+1)+1; y=(y-1)*(pixwid+1)+1;
        M(x:x+pixwid-1,y:y+pixwid-1)=reshape(E(:,I),pixwid,pixwid);
        e(I)=norm(E(:,I));
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(2,3,3);
    plot(D);

    subplot(2,3,6);
    plot((1-LS)*100);

end
                  
function [A,B,YS_out]=gaussian_component_model(dataset,bank,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('filtered',false);
pr.parse(varargin{:});
options=pr.Results;

[~,ndata]=size(dataset);
[~,bankdim]=size(bank);

if options.filtered
    XS=dataset;
else
    XS=bank'*dataset;
    % XS(XS<0)=0;
end

if ~isempty(options.params)
    A=options.params{1};
    B=options.params{2};
else
    A=randn(bankdim,nunits)*0.1;   % variance
    A=A./repmat(sqrt(sum(A.^2,1)),bankdim,1);
    B=randn(bankdim,1)*0.01;            % mean
    B(:)=0;
end;

%%%%%

nsteps_out=10000;
nsteps_in=100;
minibatch_size=100;
beta=0.01;
alpha=0.01;
eta=0.01;
% rho=2;
rho=0.1;
lambda=0.1;

dB=zeros(bankdim,1);
dA=zeros(bankdim,nunits);

LS_out=zeros(1,nsteps_out);
LS_mean=[];

YS_out=zeros(nunits,nsteps_out);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps_out=options.nsteps; end;

% g=@(x)log(1+exp(x));
% gd=@(x)1./(1+exp(-x));
% gi=@(x)log(exp(x)-1);

g=@(x)exp(x);
gd=@(x)exp(x);
gi=@(x)log(x);

for T=1:nsteps_out
    if mod(T,minibatch_size)==1 
        L=mean(LS_out(max(1,T-1000):T-1));
        LS_mean=[LS_mean;L];
        fprintf('step #%d (mean log likelihood=%1.5f)\n',T,L);
%         B=B+beta*dB;
        A=A+alpha*dA;
        A=A./repmat(sqrt(sum(A.^2,1)),bankdim,1);
        dB(:)=0;
        dA(:)=0;
        if options.display disp_weight(h,A,bank,LS_mean); end;
    end;
    if mod(T,nsteps_out/8)==0
        beta=beta/2; alpha=alpha/2;
        fprintf('alpha=%e beta=%e\n',beta,alpha);
    end;
    X=XS(:,mod(T-1,ndata)+1);
    
    % inference

    Y=A\gi((X-B).^2)*0.01; 
%     Y=randn(nunits,1)*0.01;
    XS_var=zeros(bankdim,nsteps_in);
    
    YS=zeros(nunits,nsteps_in);
    YS(:,1)=Y;
    LS_in=zeros(1,nsteps_in);

    eta2=eta;
    L=-Inf;
    for S=1:nsteps_in
        X1=X-B;
        X2=A*Y;
        X2g=g(X2);
        X2gd=gd(X2);
%         Z=X2gd.*(X1.^2./(X2g.^2)-1./X2g)/2;
%         Z=X2gd.*(abs(X1)./(X2g.^2)/lambda-1./X2g);
%         dY=A'*Z-sign(Y)/rho;
%         dY=A'*Z-Y/rho;
        dY=A'*(X1.^2-X2)/lambda-sign(Y)/rho;
%         dY=A'*(abs(X1)-X2)/lambda-sign(Y)/rho;
        for I=1:10
            Y1=Y+eta2*dY;
%             L1=sum(-X1.^2./X2g/2-log(X2g)/2)/lambda-sum(abs(Y))/rho;
%             L1=sum(-X1.^2./X2g/2-log(X2g)/2)/lambda-sum(Y.^2)/rho/2;
%             L1=sum(-abs(X1)./X2g/lambda-log(X2g))-sum(Y.^2)/rho/2;
%             L1=sum(-abs(X1)./X2g/lambda-log(X2g))-sum(abs(Y))/rho;
            L1=sum(-(X1.^2-X2).^2/lambda/2)-sum(abs(Y))/rho;
%             L1=sum(-(abs(X1)-X2).^2/lambda/2)-sum(abs(Y))/rho;
            if L1>L Y=Y1; L=L1; break; end;
            eta2=eta2/2;
        end;
        YS(:,S)=Y;
        XS_var(:,S)=X2;
        LS_in(1,S)=L;
        if eta2<eta*1e-4 break; end;
    end;

%     subplot(2,2,1);plot(YS');title('latent');
%     subplot(2,2,2);plot(XS_var');title('s.d.');
%     subplot(2,2,4);plot(LS_in(1,1:S)');title('log likelihood');
%     pause;
    
    % update
%     disp(Y');

    X1=X-B;
    X2=A*Y;
    X2g=g(X2);
    X2gd=gd(X2);
%     Z=X2gd.*(X1.^2./(X2g.^2)-1./X2g);
%     Z=X2gd.*(abs(X1)./(X2g.^2)/lambda-1./X2g);
%     L=sum(-X1.^2./X2g/2-log(X2g)/2)/lambda-sum(abs(Y))/rho;
%     L=sum(-X1.^2./X2g/2-log(X2g)/2)/lambda-sum(Y.^2)/rho/2;
%     L=sum(-abs(X1)./X2g/lambda-log(X2g))-sum(Y.^2)/rho/2;
%     L=sum(-abs(X1)./X2g/lambda-log(X2g))-sum(abs(Y))/rho;
    L=sum(-(X1.^2-X2).^2/lambda/2)-sum(abs(Y))/rho;
%     L=sum(-(abs(X1)-X2).^2/lambda/2)-sum(abs(Y))/rho;
    LS_out(1,T)=L;
    YS_out(:,T)=Y;

%     dB=dB+X1./X2g;
%     dA=dA+Z*Y';
    dA=dA+(X1.^2-X2)*Y';
%     dA=dA+(abs(X1)-X2)*Y';
    
end;

% MLS=zeros(1,nsteps_out);
% for I=1:nsteps_out MLS(I)=mean(LS_out(max(1,I-1000):I));end;
% figure;plot(1:nsteps_out,MLS');

end


function disp_weight(h,A,bank,LS)
    figure(h); clf;
    subplot(2,3,[1 4]);
    [~,nunits]=size(A);
    [sz,~]=size(bank);
    wid=sqrt(sz);

    nsub=5;
    M=ones(nsub*(wid+1),nunits*(wid+1));
    for I=1:nunits
        [~,idx]=sort(A(:,I),'descend');
        for J=1:nsub
            x=(J-1)*(wid+1)+1; y=(I-1)*(wid+1)+1;
%             M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(J))/max(bank(:,idx(J))),wid,wid);
            M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(J)),wid,wid);
        end;
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(2,3,[2 3]);
    imagesc(A);
    colorbar;

    subplot(2,3,[5 6]);
    plot(LS);

end
                  
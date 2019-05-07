function [YS_out,A,W]=sparse_coding(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('dcremoval',true,@islogical);
pr.addParamValue('lastEig',NaN,@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

[datadim,ndata]=size(dataset);

fprintf('Input dimensions: %d\n',datadim);
fprintf('Data points: %d\n',ndata);

if options.dcremoval
    fprintf('Removing DC components...\n');
    dataset=bsxfun(@minus,dataset,mean(dataset, 1));
end;

dataset=bsxfun(@minus,dataset,mean(dataset,2));

if isnan(options.lastEig) options.lastEig=datadim; end;
fprintf('Whitening data with selected %d dimensions...\n', options.lastEig);
[E,~,D]=pca(dataset','NumComponents',options.lastEig);
fprintf('Eigen value max=%f min=%f\n',D(1),D(options.lastEig));

% V=E*(diag(D(1:options.lastEig).^(-1/2)));
% V=E*(diag(D(1:opt.lastEig).^-0.2));

% V=E*(diag(D(1:options.lastEig).^(-1/2)))*E';
% whiteData = V'*dataset;

whiteData=dataset;

inputdim=size(whiteData,1);

if ~isempty(options.params)
    A=options.params;
else
    A=randn(inputdim,nunits)*0.1;   % basis matrix
    A=bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));
end;

%%%%%

nsteps_out=10000;
nsteps_in=100;
minibatch_size=100;
alpha=0.0001;
eta=0.01;
% rho=2;

% rho=0.15^2;
% lambda=0.1;

signal_sd=std(whiteData(:));
rho=signal_sd^2*0.12;
lambda=signal_sd/4;

% rho=1;
% lambda=4.5;

dA=zeros(inputdim,nunits);

LS_out=zeros(1,nsteps_out);
LS_mean=[];

YS_out=zeros(nunits,nsteps_out);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps_out=options.nsteps; end;

for T=1:nsteps_out
    if mod(T,minibatch_size)==1 
        L=mean(LS_out(max(1,T-1000):T-1));
        output_sd=std(YS_out(:,max(1,T-1000):T-1),[],2);
        LS_mean=[LS_mean;L];
        fprintf('step #%d (mean log likelihood=%1.5f, std out=%1.5f)\n',T,L,mean(output_sd));
        A=A+alpha*dA;
%         if ~isnan(output_sd) A=bsxfun(@times,A,output_sd'./signal_sd); end;
        A=bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));
        dA(:)=0;
        if options.display disp_weight(h,A,LS_mean); end;
    end;
    if mod(T,nsteps_out/8)==0
        alpha=alpha/2;
        fprintf('alpha=%e\n',alpha);
    end;
%     X=dataset(:,mod(T-1,ndata)+1);
    X=whiteData(:,mod(T-1,ndata)+1);
    
    % inference

%     Y=A\X*1; 
%     Y=A'*X*1e-5; 
    Y=rand(nunits,1); 
    
    YS=zeros(nunits,nsteps_in);
    YS(:,1)=Y;
    LS_in=zeros(1,nsteps_in);

    eta2=eta;
    L=-Inf;
    for S=1:nsteps_in
        Z=X-A*Y; 
        dY=A'*Z/rho-sign(Y)/lambda;
        for I=1:10
            Y1=Y+eta2*dY;
            Z1=X-A*Y1;
            L1=-sum(Z1.^2)/rho/2-sum(abs(Y1))/lambda;
            if L1>L Y=Y1; L=L1; break; end;
            eta2=eta2/2;
        end;
        YS(:,S)=Y;
        LS_in(1,S)=L;
        if eta2<eta*1e-10 break; end;
    end;

%     subplot(2,2,2);plot(YS(:,1:S)');title('latent');
%     subplot(2,2,4);plot(LS_in(1,1:S)');title('log likelihood');
%     pause;
    
    % update
%     disp(Y');

    Z=X-A*Y;
    L=-sum(Z.^2)/rho/2-sum(abs(Y))/lambda;
    LS_out(1,T)=L;
    YS_out(:,T)=Y;

    dA=dA+Z*Y'/rho;

    
end;

% A=V*A;
W=pinv(A);

% MLS=zeros(1,nsteps_out);
% for I=1:nsteps_out MLS(I)=mean(LS_out(max(1,I-1000):I));end;
% figure;plot(1:nsteps_out,MLS');

end


function disp_weight(h,A,LS)
    figure(h); clf;
    subplot(2,1,1);
    [datadim,nunits]=size(A);
    pixwid=sqrt(datadim);
    wid=ceil(sqrt(nunits)); ht=ceil(nunits/wid);

    M=ones(wid*(pixwid+1),ht*(pixwid+1))*min(A(:));
    for I=1:nunits
        [y x]=ind2sub([wid ht],I);
        x=(x-1)*(pixwid+1)+1; y=(y-1)*(pixwid+1)+1;
        M(x:x+pixwid-1,y:y+pixwid-1)=reshape(A(:,I),pixwid,pixwid);
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(2,1,2);
    plot(LS);
    
    drawnow;

end
                  
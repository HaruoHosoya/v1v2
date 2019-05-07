function sparse_autoencoder(dataset,nunits,nsteps)

[datadim,ndata]=size(dataset);

% nonlin=@(x)abs(x);
% nonlind=@(x)sign(x);
nonlin=@(x)log(cosh(x));
nonlind=@(x)tanh(x);

W=randn(datadim,nunits)*0.01;
W=W./repmat(sqrt(sum(W.^2,1)),datadim,1);
dW=zeros(datadim,nunits);

h=figure;
alpha=0.005;
lambda=5;

% for T=1:nsteps
%     if mod(T,1000)==1 
%         fprintf('step #%d %1.3f\n',T,norm(dW)); 
%         disp_weight(h,W,ES);
%         W=W+alpha*dW;
%         dW=zeros(datadim,nunits);
%         W=W./repmat(sqrt(sum(W.^2,1)),datadim,1);
%     end;
% 
%     X=dataset(:,mod(T,ndata)+1);
%     Y=W'*X;
%     
%     dW=dW+((X-W*Y)*Y'+X*(X-W*Y)'*W)-lambda*(X*nonlind(Y)');
% 
%     e=-1/2*norm(X-W*Y)^2-lambda*nonlin(Y);
% end

    function [e,d]=obj(W)
        e=0; d=zeros(datadim,nunits);
        X=dataset;
        Y=W'*X; Z=W*Y; D=X-Z;
        e=e-1/2*sum(sum(D.^2))-lambda*sum(nonlin(Y(:)));
        d=d+D*Y'+X*D'*W-lambda*(X*nonlind(Y)');
        d=d(:);
        e=-e; d=-d;
    end

options=optimset(...
    'Display','iter',...
    'Algorithm','interior-point',...
    'GradObj','on','DerivativeCheck','off',...
    'LargeScale','on');
% options=optimset(...
%     'Display','off','TolFun',record.tol,'TolX',record.tolx,'MaxIter',record.maxiter,'MaxFunEvals',record.maxfunevals,...
%     'PlotFcns',dispfcn,'OutputFcn',@output,...
%     'Algorithm','interior-point',...
%     'Hessian','lbfgs',...
%     'GradObj','on','DerivativeCheck','off');
[W,fval,exitflag,output]=fminunc(@obj,W,options);

disp_weight(h,W,[]);

end

function disp_weight(h,W,E)
    figure(h); clf;
    [datadim,nunits]=size(W);
    width=floor(sqrt(nunits)); height=ceil(nunits/width);
    imwid=sqrt(datadim);
    
    M=ones(width*(imwid+1),height*(imwid+1));
    for I=1:height
        for J=1:width
            x=(J-1)*(imwid+1)+1; y=(I-1)*(imwid+1)+1;
            M(x:x+imwid-1,y:y+imwid-1)=reshape(W(:,(I-1)*width+J),imwid,imwid);
        end;
    end;

%     subplot(121);
    imagesc(M'); colormap(gray); axis image;
    colorbar;
    
%     subplot(122);
%     plot(1:length(E),E);

end


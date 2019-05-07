function [W,ED,YS]=incremental_tica(dataset,nsub,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('alpha',0.1,@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('filters',[]);
pr.parse(varargin{:});
options=pr.Results;

[datadim,ndata]=size(dataset);

% lowpass & whitening

wdim=nsub;

if isempty(options.filters)
    [E,~,evals]=pca(dataset');
    D=evals(1:wdim).^(-1/2);
    ED=E(:,1:wdim)*diag(D);
else
    ED=options.filters;
end;

wdata=ED'*dataset;

if ~isempty(options.params)
    W=options.params;
else
    W=randn(wdim,nsub)*0.1;
    W=W./repmat(sqrt(sum(W.^2,1)),wdim,1);
end;

%%%%%

nsteps=10000;
minibatch_size=100;
display_steps=1000;
alpha=options.alpha;
beta=0.5;

dW=zeros(wdim,nsub);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

% Sparseness nonlinearity

% G=@(x)sqrt(x+eps);
% g=@(x)1./(2*sqrt(x+eps));


% G=@(x)log(cosh(x));
% g=@(x)tanh(x);
% gd=@(x)1-tanh(x).^2;

% G=@(x)-log(cosh(x));
% g=@(x)-tanh(x);
% gd=@(x)-1+tanh(x).^2;
% 

% G=@(x)-log(cosh(x))+x/2;
% g=@(x)-tanh(x)+1/2;
% gd=@(x)-1+tanh(x).^2;

% G=@(x)tanh(x/2);
% g=@(x)0.5*(1-tanh(x/2).^2);

% G=@(u)-exp(-u.^2/2);
% g=@(u)u.*exp(-u.^2/2);
% gd=@(u)(1-u.^2).*exp(-u.^2/2);

G=@(x)-0.4*x+sqrt(x+0.5);
g=@(x)-0.4+1./(2*sqrt(x+0.5));

% General output nonlinearity

% F=@(x)(x<0).*0+(x>=0).*x;
% f=@(x)(x<0).*0+(x>=0).*1;

F=@(x)(x<0).*0+(x>=0).*x.^2;
f=@(x)(x<0).*0+(x>=0).*2.*x;


% Topographic ICA

% wid=floor(sqrt(nsub));
% nunits=nsub;
% L=zeros(nsub,nunits);
% M=zeros(nsub,nunits);
% H=zeros(nsub,nunits);
% nei=-1:1;
% % nei=0:1;
% for x=1:wid
%     for y=1:wid
%         [ix iy]=ndgrid(x+nei,y+nei);
% %         ix=mod(ix-1,wid)+1; iy=mod(iy-1,wid)+1;
%         idx=ix>=1 & ix<=wid & iy>=1 & iy<=wid; ix(~idx)=[]; iy(~idx)=[];
%         I=sub2ind([wid wid],ix(:),iy(:));
%         J=sub2ind([wid wid],x,y);
% %         H(I,repmat(J,length(ix(:)),1))=1/length(ix(:))/2;
%         H(I,repmat(J,length(ix(:)),1))=1;
%     end;
% end;

%
% ISA
% 

% L=zeros(nsub,nunits);
% M=zeros(nsub,nunits);
% H=zeros(nsub,nunits);
% ncomp=nsub/nunits;
% for I=1:nunits
%     for J=1:ncomp
%         H((I-1)*ncomp+J,I)=1;
%     end;
% end;

%
% Quadratic model
% 

% H=zeros(nsub,nunits);
% L=zeros(nsub,nunits);
% M=zeros(nsub,nunits);
% ncomp=nsub/nunits;
% for I=1:nunits
%     L((I-1)*ncomp+1,I)=1;
%     for J=2:ncomp
%         H((I-1)*ncomp+J,I)=1;
%     end;
% end;

% Half topographic ICA

wid=floor(sqrt(nsub));
nunits=nsub;
L=zeros(nsub,nunits);
M=zeros(nsub,nunits);
H=zeros(nsub,nunits);
nei=-1:1;
% nei=0:1;
for x=1:wid
    for y=1:wid
        [ix iy]=ndgrid(x+nei,y+nei);
        ix=mod(ix-1,wid)+1; iy=mod(iy-1,wid)+1;
%         idx=ix>=1 & ix<=wid & iy>=1 & iy<=wid; ix(~idx)=[]; iy(~idx)=[];
        I=sub2ind([wid wid],ix(:),iy(:));
        J=sub2ind([wid wid],x,y);
        M(I,repmat(J,length(ix(:)),1))=1/length(ix(:));
%         M(I,repmat(J,length(ix(:)),1))=1;
    end;
end;

%
% Half ISA
% 

% L=zeros(nsub,nunits);
% M=zeros(nsub,nunits);
% H=zeros(nsub,nunits);
% ncomp=nsub/nunits;
% for I=1:nunits
%     for J=1:ncomp
%         M((I-1)*ncomp+J,I)=1/ncomp;
%     end;
% end;

% Learning

YS=zeros(nunits,nsteps);

for T=1:nsteps
    if mod(T,minibatch_size)==1 
        W=W+dW/minibatch_size;
        W=W+beta*W*(eye(wdim)-W'*W);
%         A=A*(A'*A)^(-1/2);
        W=W./repmat(sqrt(sum(W.^2,1)),wdim,1);
%         A=orthonorm(A);
        dW(:)=0;
    end;
    if mod(T,display_steps)==1  
        fprintf('step #%d\n',T);
        if options.display
            disp_weight(h,W,ED); 
        end;
    end;
%     if mod(T,nsteps/8)==0
%         alpha=alpha/2;
%         fprintf('alpha=%e\n',alpha);
%     end;

    % inference

    X=wdata(:,mod(T-1,ndata)+1);
    U=W'*X;
    Y=L'*U+1/2*H'*(U.^2)+M'*F(U);
    V1=L*g(Y);
    V2=U.*(H*g(Y));
    V3=f(U).*(M*g(Y));
        
    dW=dW-alpha*X*(V1'+V2'+V3');

    YS(:,T)=U;
    
    
end;

end

function B=orthonorm(A)
    B=A; 
    for I=1:size(A,2)
       for J=1:I-1
           B(:,I)=B(:,I)-B(:,J)'*B(:,I)*B(:,J);
       end
       B(:,I)=B(:,I)/norm(B(:,I));
    end
end


function disp_weight(h,W,ED)
    figure(h); clf;
    subplot(1,1,1);
    A=pinv(ED)'*W;
%     A=ED*W;
    A=A./repmat(sqrt(sum(A.^2,1)),size(A,1),1);
    [datadim,nsub]=size(A);
    pixwid=sqrt(datadim);
    wid=ceil(sqrt(nsub)); ht=ceil(nsub/wid);

    M=ones(wid*(pixwid+1),ht*(pixwid+1))*min(A(:));
    e=zeros(nsub,1);
    for I=1:nsub
        [x y]=ind2sub([wid ht],I);
        x=(x-1)*(pixwid+1)+1; y=(y-1)*(pixwid+1)+1;
        v=A(:,I);
        M(x:x+pixwid-1,y:y+pixwid-1)=reshape(v,pixwid,pixwid);
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;


end
                  
function A=nonnegative_pooling(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('nsteps',50,@isnumeric);
pr.addParamValue('alpha',1,@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

[datadim,ndata]=size(dataset);


alpha=options.alpha;

A=rand(datadim,nunits);

% f=@(x)tanh(x/2);
% g=@(x)0.5*(1-tanh(x/2).^2);

% f=@(x)x.^3/3;
% g=@(x)x.^2;

f=@(x)log(cosh(x));
g=@(x)tanh(x);
gd=@(x)1-tanh(x).^2;
 
% f=@(u)-exp(-u.^2/2);
% g=@(u)u.*exp(-u.^2/2);
% gd=@(u)(1-u.^2).*exp(-u.^2/2);

for r=1:options.nsteps
    L=-sum(sum(f(dataset'*A)));
    fprintf('step #%d  L=%f\n',r,L);
    dA=-dataset*g(dataset'*A);
    alpha1=alpha;
    for u=1:10
        A1=A+alpha1*dA/ndata;
        L1=-sum(sum(f(dataset'*A1)));
        if L1>L A=A1; break; end;
        alpha1=alpha1/2;
    end;
    A=A*(A'*A)^(-1/2);
    s=max(A,[],1)>-min(A,[],1);
    A=bsxfun(@times,A,s*2-1);
%     A(A<0)=0;
    A=bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));    
end;
    
end


function disp_weight(h,A,wfilters)
    figure(h); clf;
    subplot(1,1,1);
%     W=pinv(wfilters*A)';
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
                  
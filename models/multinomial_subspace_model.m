function multinomial_subspace_model(dataset,nunits,nsubunits,nsteps)

[datadim,ndata]=size(dataset);

W=randn(datadim,nsubunits,nunits)*0.1+0.1;
W=W./repmat(sqrt(sum(W.^2,1)),datadim,1);
dW=zeros(datadim,nsubunits,nunits);

h=figure;
alpha=0.005;

for T=1:nsteps
    if mod(T,100)==1 
        fprintf('step #%d %1.3f\n',T,0); 
        disp_weight(h,W);
        W=W+alpha*dW;
        dW=zeros(datadim,nsubunits,nunits);
    end;

    X=dataset(:,mod(T,ndata)+1);
    
    % inference

    p=zeros(nunits,1);
    for I=1:nunits
        W1=W(:,:,I);
        p(I)=p(I)-1/2*(log(det(W1*W1'))-X'/(W1*W1')*X);
    end;
    p=exp(p); p=real(p); p(p<0)=0;
    if all(p==0) k=randi(nunits,1); 
    else p=p/sum(p); k=randp(p);
    end;
    ks(k)=ks(k)+1;

    % update
    
    Wk=W(:,:,k);
    X1=(Wk*Wk')\X;
    Y1=Wk'*X1;

    dW(:,:,k)=dW(:,:,k)-inv(W1)+X1*Y1';
end


disp_weight(h,W,[]);

end

function disp_weight(h,W)
    figure(h); clf;
    [sz,nsub,nunits]=size(W);
    wid=sqrt(sz);

    M=ones(nsub*(wid+1),nunits*(wid+1));
    for A=1:nunits
        for B=1:nsub
            x=(B-1)*(wid+1)+1; y=(A-1)*(wid+1)+1;
            W1=W(:,B,A);
            M(x:x+wid-1,y:y+wid-1)=reshape(W1/max(W1),wid,wid);
        end;
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;


end

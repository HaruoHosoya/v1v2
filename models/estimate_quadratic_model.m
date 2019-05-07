function [A,B,C,f]=estimate_quadratic_model(dataset,outset,nquad)

sigma=0.1;
[datadim,ndata]=size(dataset);

A0=randn(datadim,nquad)*0.1;
B0=randn(datadim,1)*0.1;
C0=randn(1)*0.1;
ps0=fold_params(A0,B0,C0);

lb=ones(length(ps0),1)*-Inf;
ub=ones(length(ps0),1)*Inf;

opts=optimset('Algorithm','interior-point','DerivativeCheck','off','GradObj','on','Display','off','Hessian','lbfgs','TolFun',1e-4);

[ps,L,eflag]=fmincon(@obj,ps0,[],[],[],[],lb,ub,[],opts);
[A,B,C]=unfold_params(ps,nquad);
f=sum((A'*dataset).^2,1)+B'*dataset+C;

    function [L,dL]=obj(ps)        
        [A,B,C]=unfold_params(ps,nquad);
        f=sum((A'*dataset).^2,1)+B'*dataset+C;
        d=outset-f;
        dC=sum(d,2);
        dB=sum(dataset.*repmat(d,datadim,1),2);
        dA=zeros(datadim,nquad);
        for I=1:nquad
            dA(:,I)=2*sum(dataset.*repmat((A(:,I)'*dataset).*d,datadim,1),2);
        end;
        L=1/(2*sigma^2)/ndata*sum(d.^2,2);
        dL=-1/(sigma^2)/ndata*[dA(:);dB;dC];
    end
            


end

function ps=fold_params(A,B,C)
    
ps=[A(:);B;C];

end

function [A,B,C]=unfold_params(ps,nquad)

D=(length(ps)-1)/(nquad+1);
A=reshape(ps(1:D*nquad),D,nquad);
B=ps(D*nquad+1:D*(nquad+1));
C=ps(end);

end



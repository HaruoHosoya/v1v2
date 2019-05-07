function net=v2_anzai_analyze2(net,layer)

nunits=net.structure.layers{layer}.numUnits;
nnodes=net.structure.layers{layer}.width^2;
subrf_dim=num2cell(net.content.layers{layer}.layerProperties.subrf_dim);
[npos,~,nori,nfreq,npha]=subrf_dim{:};

dsubrf=net.content.layers{layer}.unitProperties.double_subrf;
dsubrf=reshape(dsubrf,[nori,npos^2,nori,npos^2,nunits,nnodes]);
    
tsi=zeros(npos^2,npos^2,nunits,nnodes);
tsip=zeros(npos^2,npos^2,nunits,nnodes);

min_resp=1e-3;

for N=1:nnodes
    parfor U=1:nunits
        t=zeros(npos^2,npos^2);
        tp=zeros(npos^2,npos^2);
        fprintf('layer %d, node %d, unit %d\n',layer,N,U);
        for P1=1:npos^2
            clf;J=1;
            for P2=1:npos^2               
                subplot(6,6,J);J=J+1;
                r2=reshape(dsubrf(:,P1,:,P2,U,N),nori,nori);
                [~,idx]=max(r2(:));
                [o1,o2]=ind2sub([nori,nori],idx);
                r_opt=r2(:,o2);
                r_inv=r2(:,mod(o2-1+nori/2,nori)+1);
                if all(r_opt<min_resp) || all(r_inv<min_resp) 
                    t(P1,P2)=NaN;
                    tp(P1,P2)=NaN;
                else
                    mdl=fitlm(r_opt,r_inv,'linear');
                    t(P1,P2)=mdl.Coefficients.Estimate(2);
                    tp(P1,P2)=mdl.Coefficients.pValue(2);
                end;
            end;
        end;
        tsi(:,:,U,N)=t;
        tsip(:,:,U,N)=tp;
    end;
end;

figure;hist3([tsi(:) tsip(:)>0.05],'Edges',{-2:0.1:2 [0 1]})

    
net.content.layers{layer}.unitProperties2.tsi=tsi;
net.content.layers{layer}.unitProperties2.tsip=tsip;



end



function net=v2_anzai_test(net,nori,npha,freqs,npos,patchwid)

oris=0:pi/nori:pi-pi/nori;
phas=0:2*pi/npha:2*pi-2*pi/npha;
nfreq=length(freqs);

wid=net.structure.layers{1}.width;
intv=(wid-patchwid)/(npos-1);
pos=1:intv:1+intv*(npos-1);

[pha,freq,ori]=ndgrid(phas,freqs,oris);

[xi,yi]=meshgrid(1:patchwid,1:patchwid);

dataset=zeros(wid,wid,numel(ori),npos,npos);

for YI=1:npos
    for XI=1:npos
        for I=1:numel(ori)
            img=grating(xi,yi,ori(I),freq(I),pha(I));
            x=pos(XI); y=pos(YI);
            dataset(y:y+patchwid-1,x:x+patchwid-1,I,XI,YI)=img;
        end;
    end;
end;

dataset=reshape(dataset,wid^2,numel(ori)*npos^2);
net=run_net(net,dataset);

for L=4:length(net.content.layers)
    fprintf('single patch test in layer %d\n',L);
    nunits=net.structure.layers{L}.numUnits;
    nnodes=net.structure.layers{L}.width^2;
    
    net.content.layers{L}.layerProperties.subrf_dim=[npos,npos,nori,nfreq,npha];
    net.content.layers{L}.layerProperties.orients=oris';
    net.content.layers{L}.layerProperties.phases=phas';
    net.content.layers{L}.layerProperties.freqs=freqs;
    net.content.layers{L}.layerProperties.positions=pos';
    resp=reshape(net.content.layers{L}.unitProperties.resp,[npha*nfreq*nori*npos^2,nunits,nnodes]);
    resp=bsxfun(@rdivide,resp,max(resp,[],1));
    net.content.layers{L}.unitProperties.subrf=resp;
    [~,idx]=max(reshape(resp,[npha*nfreq,nori*npos*npos,nunits,nnodes]),[],1);
    [optpha,optfreq]=ind2sub([npha nfreq],idx);
    net.content.layers{L}.unitProperties.optpha=optpha;
    net.content.layers{L}.unitProperties.optfreq=optfreq;
end;


% for L=4:length(net.content.layers)
%     fprintf('double patch test in layer %d\n',L);
%     nunits=net.structure.layers{L}.numUnits;
%     nnodes=net.structure.layers{L}.width^2;
%     optpha=net.content.layers{L}.unitProperties.optpha;
%     optfreq=net.content.layers{L}.unitProperties.optfreq;
%     resp=zeros((nori*npos*npos)^2,nunits,nnodes);
%     for N=1:nnodes
%         for U=1:nunits
%             fprintf('layer %d, node %d, unit %d\n',L,N,U);
%             grat=zeros(wid,wid,nori*npos*npos);
%             J=1;
%             for YI=1:npos
%                 for XI=1:npos
%                     for O=1:nori
%                         x=pos(XI); y=pos(YI);
%                         grat(x:x+patchwid-1,y:y+patchwid-1,J)=grating(xi,yi,oris(O),freqs(optfreq(J)),phas(optpha(J)));
%                         J=J+1;
%                     end;
%                 end;
%             end;
%             imgset=zeros(wid,wid,nori*npos*npos,nori*npos*npos);
%             for J1=1:nori*npos^2
%                 for J2=1:nori*npos^2
%                     imgset(:,:,J2,J1)=grat(:,:,J1)+grat(:,:,J2);
%                 end;
%             end;
%             net1=run_net(net,reshape(imgset,wid^2,(nori*npos*npos)^2));
%             resp(:,U,N)=net1.content.layers{L}.unitProperties.resp(:,U,N);            
%         end;
%     end;
%     resp=bsxfun(@rdivide,resp,max(resp,[],1));
%     net.content.layers{L}.unitProperties.double_subrf=resp;
% end;
    
    



end



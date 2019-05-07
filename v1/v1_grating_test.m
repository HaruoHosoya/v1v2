function net=v1_grating_test(net,nori,npha,freqs)

oris=0:pi/nori:pi-pi/nori;
phas=0:2*pi/npha:2*pi-2*pi/npha;
[pha,freq,ori]=ndgrid(phas,freqs,oris);
nfreq=length(freqs);

wid=net.structure.layers{1}.width;
[xi,yi]=ndgrid(1:wid,1:wid);

dataset=zeros(wid^2,numel(ori));

for I=1:numel(ori)
    img=grating(xi,yi,ori(I),freq(I),pha(I));
    dataset(:,I)=img(:);
end;

net=run_net(net,dataset);

for L=1:length(net.content.layers)
    nunits=net.structure.layers{L}.numUnits;
    nnodes=net.structure.layers{L}.width^2;
    net.content.layers{L}.layerProperties.orients=oris';
    net.content.layers{L}.layerProperties.phases=phas';
    net.content.layers{L}.layerProperties.freqs=freqs;
    net.content.layers{L}.unitProperties3.grating=...
        reshape(net.content.layers{L}.unitProperties.resp,[npha,nfreq,nori,nunits,nnodes]);
end;


end

function [netf,netb]=v1_simple_cells(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('fixed',false,@islogical);
pr.addParamValue('epsilon',0.0005,@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));

if options.fixed
    nx=0.4; ny=0.4; grid=datawid/4:datawid/8:datawid-datawid/4;
    freq=[1/4 1/6 1/8]; ori=0:pi/8:pi-pi/8; pha=0:pi/2:2*pi-pi/2;
    exponent=0;
    nunits=length(grid)^2*length(freq)*length(ori)*length(pha);
end;

net=struct;
net=addlayer(net,1,'LGN',datawid,1,datawid);
net=addlayer(net,2,'V1s',1,nunits,1);
net.content.layers{2}.layerProperties.meanrem=false;
net.content.layers{2}.layerProperties.nonlin='linear';
net.content.layers{2}.layerProperties.dcrem=false;
net=calc_coverage(net);

if options.fixed
    [bank,params]=gabor_bank(datawid,nx,ny,grid,freq,ori,pha,exponent);
    A=reshape(bank,datadim,nunits);
    W=A';
    net.content.layers{2}.unitProperties.params=reshape(params,8,nunits);
else
    [A,W]=fastica(dataset,'approach','symm','numOfIC',nunits,'g','tanh','epsilon',options.epsilon,'verbose','on','lastEig',nunits);
end;

net.content.layers{2}.weights=reshape(W',[1,datadim,nunits,1]);

net=v1s_analyze(net);

netf=net;

netb=net;
netb.content.layers{2}.weights=reshape(A,[1,datadim,nunits,1]);


end


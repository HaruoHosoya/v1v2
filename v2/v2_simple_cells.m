function [netf,netb,Y]=v2_simple_cells(dataset,v2s_nunits,varargin)

pr=inputParser;
pr.addParamValue('nonlin','half-rect',@isstr);
pr.addParamValue('reduction',8,@isnumeric);
pr.addParamValue('learning','ica',@isstr);
pr.addParamValue('lastEig',NaN,@isnumeric);
pr.addParamValue('epsilon',0.0005,@isnumeric);
pr.addParamValue('nsteps',10000,@isnumeric);
pr.addParamValue('maxIter',100,@isnumeric);
pr.addParamValue('dcrem',true,@islogical);
pr.addParamValue('energy_model',false,@islogical);
pr.addParamValue('norient',8,@isnumeric);

pr.parse(varargin{:});
options=pr.Results;

[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));
net=struct;

% v1 simple

pos=6:4:26;
oris=0:pi/options.norient:pi-pi/options.norient;
freqs=[1/8 1/6 1/4];
% phases=[0 pi/2];
phases=[0 pi/2 pi pi/2*3];
[v1s_bank,v1s_params]=gabor_bank(datawid,0.4,0.4,pos,freqs,oris,phases,1.15);
[~,npha,nori,nfreq,nx,ny]=size(v1s_bank);
v1s_nunits=npha*nori*nfreq*nx*ny;

net=addlayer(net,1,'LGN',datawid,1,datawid);
net=addlayer(net,2,'V1s',1,v1s_nunits,1);
net.content.layers{2}.layerProperties.meanrem=false;
net.content.layers{2}.layerProperties.nonlin=options.nonlin;
net.content.layers{2}.layerProperties.dcrem=false;

net.content.layers{2}.weights=reshape(v1s_bank,[1,datadim,v1s_nunits,1]);
net.content.layers{2}.unitProperties.params=reshape(v1s_params,[8,v1s_nunits,1]);

% v1 complex


if ~options.energy_model
    v1c_nunits=v1s_nunits;
    v1c_nreduced=floor(v1c_nunits/options.reduction);
    net=addlayer(net,3,'V1c',1,v1c_nunits,1);
    net.content.layers{3}.layerProperties.nonlin='linear';
    v1c_pooling=zeros(1,v1c_nunits,v1c_nunits);
    for I=1:v1c_nunits
        v1c_pooling(:,I,I)=1;
    end;    
    net.content.layers{3}.layerProperties.meanrem=false;
    net.content.layers{3}.layerProperties.dcrem=true;
    net.content.layers{3}.weights=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);
    net.content.layers{3}.mean=zeros(1,v1s_nunits,1);
else
    v1c_nunits=v1s_nunits/npha;
    v1c_nreduced=floor(v1c_nunits/options.reduction);
    net=addlayer(net,3,'V1c',1,v1c_nunits,1);
    net.content.layers{2}.layerProperties.nonlin='square';
    net.content.layers{3}.layerProperties.nonlin='sqrt';
    v1c_pooling=zeros(npha,v1c_nunits,v1c_nunits);
    for I=1:v1c_nunits
        v1c_pooling(:,I,I)=1;
    end;
    net.content.layers{3}.layerProperties.meanrem=false;
    net.content.layers{3}.layerProperties.dcrem=true;
    net.content.layers{3}.weights=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);
    net.content.layers{3}.mean=zeros(1,v1s_nunits,1);
end;    


% v2 simple

net=addlayer(net,4,'V2s',1,v2s_nunits,1);
net.content.layers{4}.layerProperties.meanrem=true;
net.content.layers{4}.layerProperties.nonlin='linear';
net.content.layers{4}.layerProperties.dcrem=false;

net.content.layers{4}.mean=zeros(1,v1c_nunits,1);

% net.content.layers{4}.layerProperties.dcrem=options.dcrem;
% net.content.layers{4}.layerProperties.meanrem=true;
% net.content.layers{4}.layerProperties.nonlin=options.nonlin;
net=calc_coverage(net);

net=run_net(net,dataset);
v1c_outputs=reshape(net.content.layers{3}.unitProperties.resp,[datalen,v1c_nunits,1])';

if isnan(options.lastEig)
    options.lastEig=v1c_nreduced;
end;

switch options.learning
    case 'ica'
        [Y,A,W]=fastica(v1c_outputs,'approach','symm','numOfIC',v2s_nunits,'g','tanh','epsilon',options.epsilon,'verbose','on','lastEig',options.lastEig);
    case 'overcomplete-ica'
        [Y,A,W]=smica(v1c_outputs,'numOfIC',v2s_nunits,'lastEig',options.lastEig,'dcremoval',false,'maxIter',options.maxIter);
    case 'sparse-coding'
        [Y,A,W]=sparse_coding(v1c_outputs,v2s_nunits,'lastEig',options.lastEig,'dcremoval',false,'display',false,'nsteps',options.nsteps);
    case 'pca'
        [A,~,D]=pca(v1c_outputs','NumComponents',options.lastEig,'Algorithm','eig');
        W=A'; Y=W*v1c_outputs;
end;

% s=max(W,[],2)>-min(W,[],2);
% s=max(A,[],1)>-min(A,[],1); s=s';
s=skewness(Y,0,2)>=0;
W=bsxfun(@times,W,s*2-1);
A=bsxfun(@times,A,s'*2-1);
Y=bsxfun(@times,Y,s*2-1);

net.content.layers{4}.mean=reshape(mean(v1c_outputs,2),1,v1c_nunits,1);
net.content.layers{4}.weights=reshape(W',[v1c_nunits,1,v2s_nunits,1]);

net=v1c_analyze(net);
net=v2s_bases_comp(net);

net=strip_resp(net);

netf=net;
netb=net;

netb.content.layers{4}.weights=reshape(A,[v1c_nunits,1,v2s_nunits,1]);
netb=v2s_bases_comp(netb);



end

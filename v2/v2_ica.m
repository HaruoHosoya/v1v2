function [net,netf]=v2_ica(dataset,nunits,varargin)

pr=inputParser;
pr.addParamValue('nonlin','half-rect',@isstr);
pr.addParamValue('learning','overcomplete-ica',@isstr);
pr.addParamValue('lastEig',NaN,@isnumeric);
pr.addParamValue('epsilon',0.0001,@isnumeric);
pr.addParamValue('nsteps',10000,@isnumeric);
pr.addParamValue('dcrem',true,@islogical);
pr.parse(varargin{:});
options=pr.Results;

dataset=bsxfun(@minus,dataset,mean(dataset,1));
[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));

[v1s_bank,v1s_params]=make_bank(datawid);
[~,npha,nori,nfreq,nx,ny]=size(v1s_bank);
v1s_nunits=npha*nori*nfreq*nx*ny;
v1c_nunits=nori*nfreq*nx*ny;

v1c_pooling=zeros(npha,v1c_nunits,v1c_nunits);
for I=1:v1c_nunits
    v1c_pooling(:,I,I)=1;
end;

net=struct;
net=addlayer(net,1,'LGN',datawid,1,datawid);
net=addlayer(net,2,'V1s',1,v1s_nunits,1);
net.content.layers{2}.layerProperties.dcrem=false;
net.content.layers{2}.layerProperties.meanrem=false;
net.content.layers{2}.layerProperties.nonlin='square';
net=addlayer(net,3,'V1c',1,v1c_nunits,1);
net.content.layers{3}.layerProperties.dcrem=false;
net.content.layers{3}.layerProperties.meanrem=false;
net.content.layers{3}.layerProperties.nonlin='sqrt';
net=addlayer(net,4,'V2s',1,nunits,1);
net.content.layers{4}.layerProperties.dcrem=options.dcrem;
net.content.layers{4}.layerProperties.meanrem=true;
net.content.layers{4}.layerProperties.nonlin=options.nonlin;
net=calc_coverage(net);

net.content.layers{2}.filters=reshape(v1s_bank,[1,datadim,v1s_nunits,1]);
net.content.layers{2}.weights=reshape(v1s_bank,[1,datadim,v1s_nunits,1]);
net.content.layers{2}.mean=zeros(1,datadim);
net.content.layers{2}.unitProperties.params=reshape(v1s_params,[8,v1s_nunits,1]);
net.content.layers{3}.filters=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);
net.content.layers{3}.weights=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);
net.content.layers{3}.mean=zeros(v1s_nunits,1);

net=run_net(net,dataset);
v1c_outputs=reshape(net.content.layers{3}.unitProperties.resp,[datalen,v1c_nunits,1])';

if isnan(options.lastEig)
    options.lastEig=nunits;
end;

switch options.learning
    case 'ica'
        [Y,A,W]=fastica(v1c_outputs,'approach','symm','numOfIC',nunits,'g','tanh','epsilon',options.epsilon,'verbose','on','lastEig',options.lastEig);
        A2=A;
    case 'overcomplete-ica'
        [Y,A,W,~,~,A2]=smica(v1c_outputs,'numOfIC',nunits,'lastEig',options.lastEig,'dcremoval',false);
    case 'sparse-coding'
        [Y,A,W]=sparse_coding(v1c_outputs,nunits,'lastEig',options.lastEig,'dcremoval',false,'display',false,'nsteps',options.nsteps);
    case 'pca'
        [A,~,D]=pca(v1c_outputs','NumComponents',options.lastEig,'Algorithm','eig');
        A2=A; W=A'; Y=W*v1c_outputs;
end;

% s=max(W,[],2)>-min(W,[],2);
s=max(A,[],1)>-min(A,[],1); s=s';
W=bsxfun(@times,W,s*2-1);
A=bsxfun(@times,A,s'*2-1);
A2=bsxfun(@times,A2,s'*2-1);
Y=bsxfun(@times,Y,s*2-1);

A=A2;
net.content.layers{4}.filters=reshape(W',[v1c_nunits,1,nunits,1]);
net.content.layers{4}.weights2=reshape(A,[v1c_nunits,1,nunits,1]);
net.content.layers{4}.weights=reshape(A2,[v1c_nunits,1,nunits,1]);
net.content.layers{4}.mean=reshape(mean(v1c_outputs,2),v1c_nunits,1);

net=v1c_analyze(net);
net=v2s_bases_comp(net);

netf=net;

netf.content.layers{4}.weights=reshape(W',[v1c_nunits,1,nunits,1]);
netf=v1c_analyze(netf);
netf=v2s_bases_comp(netf);


end
% 
% function [bank,params]=make_bank(wid)
% 
% npos=6;
% pos=wid/(npos+1):wid/(npos+1):wid-wid/(npos+1);
% nori=4;
% oris=0:pi/nori:pi-pi/nori;
% freqs=[1/8 1/6 1/4];
% % freqs=[1/8];
% phases=[0 pi/2];
% [bank,params]=gabor_bank(wid,0.4,0.4,pos,freqs,oris,phases,1.15);
% bank=bank./repmat(sqrt(sum(bank.^2,1)),[size(bank,1),1,1,1,1,1]);
% 
% end


function [bank,params]=make_bank(wid)

pos=6:4:26;
oris=0:pi/12:pi-pi/12;
% oris=0:pi/4:pi-pi/4;
freqs=[1/8 1/6 1/4];
phases=[0 pi/2];
[bank,params]=gabor_bank(wid,0.4,0.4,pos,freqs,oris,phases,1.15);
% [bank,params]=gabor_bank(wid,0.3,0.3,pos,freqs,oris,phases,1.15);
bank=bank./repmat(sqrt(sum(bank.^2,1)),[size(bank,1),1,1,1,1,1]);

end



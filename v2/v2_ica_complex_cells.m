function net=v2_ica_complex_cells(dataset,nunits,nreduced,nshow,varargin)

pr=inputParser;
pr.addParamValue('pooling','dimred',@isstr);
pr.addParamValue('nonlin','half-rect',@isstr);
pr.addParamValue('showneg',false,@islogical);
pr.addParamValue('nsubshow',floor(nunits/nreduced),@isnumeric);
pr.addParamValue('overcomplete',false,@islogical);
pr.addParamValue('lastEig',NaN,@isnumeric);
pr.addParamValue('epsilon',0.0001,@isnumeric);
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
net.content.layers{2}.layerProperties.nonlin='square';
net=addlayer(net,3,'V1c',1,v1c_nunits,1);
net.content.layers{3}.layerProperties.dcrem=false;
net.content.layers{3}.layerProperties.nonlin='sqrt';
net=addlayer(net,4,'V2s',1,nunits,1);
net.content.layers{4}.layerProperties.dcrem=true;
net.content.layers{4}.layerProperties.nonlin=options.nonlin;
net=addlayer(net,5,'V2c',1,nunits,1);
net.content.layers{5}.layerProperties.dcrem=false;
net.content.layers{5}.layerProperties.nonlin='linear';
net=calc_coverage(net);

net.content.layers{2}.filters=reshape(v1s_bank,[1,datadim,v1s_nunits,1]);
net.content.layers{2}.weights=reshape(v1s_bank,[1,datadim,v1s_nunits,1]);
net.content.layers{2}.unitProperties.params=reshape(v1s_params,[8,v1s_nunits,1]);
net.content.layers{3}.filters=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);
net.content.layers{3}.weights=reshape(v1c_pooling,[v1s_nunits,1,v1c_nunits,1]);

net=run_net(net,dataset);
v1c_outputs=reshape(net.content.layers{3}.unitProperties.resp,[datalen,v1c_nunits,1])';

if isnan(options.lastEig)
    options.lastEig=nunits;
end;

if options.overcomplete
    [Y,A,W]=smica(v1c_outputs,'numOfIC',nunits,'lastEig',options.lastEig,'dcremoval',false);
else
    [Y,A,W]=fastica(v1c_outputs,'approach','symm','numOfIC',nunits,'g','tanh','epsilon',options.epsilon,'verbose','on','lastEig',options.lastEig);
end;

s=max(W,[],2)>-min(W,[],2);
% s=max(A,[],1)>-min(A,[],1); s=s';
W=bsxfun(@times,W,s*2-1);
A=bsxfun(@times,A,s'*2-1);
Y=bsxfun(@times,Y,s*2-1);

net.content.layers{4}.filters=reshape(W',[v1c_nunits,1,nunits,1]);
net.content.layers{4}.weights=reshape(A,[v1c_nunits,1,nunits,1]);

net=run_net(net,dataset);
v2s_outputs=reshape(net.content.layers{4}.unitProperties.resp,[datalen,nunits,1])';

switch options.pooling
    case 'ica'
        [Z,B,V]=fastica(v2s_outputs,'approach','symm','numOfIC',nreduced,'g','tanh','epsilon',options.epsilon,'verbose','on','lastEig',nreduced);
    case 'nonneg'
        % preliminary
        V=nonnegative_pooling(v2s_outputs,nreduced); V=V';
    case 'pca'
        [V,~,~]=pca(v2s_outputs');
%         V=V(:,1:nreduced)';
    case 'whitening'
        [F,~,D]=pca(v2s_outputs');
        V=F(:,1:nreduced)*diag(D(1:nreduced).^(-1/2))*F(:,1:nreduced)';
%         V=V(:,1:nreduced)';
    case 'dimred'
        [F,~,D]=pca(v2s_outputs');
        V=F(:,1:nreduced)*F(:,1:nreduced)';
%         V=V(:,1:nreduced)';
end;

s=max(V,[],2)>-min(V,[],2);
V=bsxfun(@times,V,s*2-1);

net.content.layers{5}.filters=reshape(V,[nunits,1,nunits,1]);
net.content.layers{5}.weights=reshape(V,[nunits,1,nunits,1]);
net=run_net(net,dataset);


% displaying

[V2,idx]=sort(V','descend');
if options.showneg
    idx=idx([1:options.nsubshow end-options.nsubshow+1:end],1:nshow);
else
    idx=idx([1:options.nsubshow],1:nshow);
end;

A2=A(:,idx(:));

figure; 
if options.showneg
    show_v2_filters(v1s_bank,A2,1:numel(idx),0.2,size(idx,1));
else
    show_v2_filters(v1s_bank,A2,1:numel(idx),0.2,size(idx,1));
end;

figure;
if options.showneg
    barh(V2([1:4 end-3:end],1:nshow)');
else
    barh(V2(1:4,1:nshow)');
end;
axis ij;
ylim([0 nshow+1]);
colormap(gray);
colorbar;

end

function [bank,params]=make_bank(wid)

npos=6;
pos=wid/(npos+1):wid/(npos+1):wid-wid/(npos+1);
oris=0:pi/8:pi-pi/8;
% oris=0:pi/4:pi-pi/4;
% freqs=[1/8 1/6 1/4];
freqs=[1/8];
phases=0:pi/2:2*pi-pi/2;
[bank,params]=gabor_bank(wid,0.4,0.4,pos,freqs,oris,phases);
bank=bank./repmat(sqrt(sum(bank.^2,1)),[size(bank,1),1,1,1,1,1]);

end



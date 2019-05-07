function [netf,nete]=v2_complex_cells(net,dataset,nreduced,varargin)

pr=inputParser;
pr.addParamValue('pooling','dimred',@isstr);
pr.addParamValue('nonlin','half-rect',@isstr);
pr.parse(varargin{:});
options=pr.Results;

[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));
nunits=net.structure.layers{4}.numUnits;

fprintf('building a V2 complex cells model with %d units from %d reduced eigenspace\n',nunits,nreduced);
fprintf('  using %s method and %s nonlinearity \n',options.pooling,options.nonlin);

net.content.layers{4}.layerProperties.nonlin=options.nonlin;
net=addlayer(net,5,'V2c',1,nunits,1);
net.content.layers{3}.layerProperties.meanrem=false;
net.content.layers{5}.layerProperties.nonlin='half-rect';
net.content.layers{5}.layerProperties.dcrem=false;
net=calc_coverage(net);

net=run_net(net,dataset);
v2s_outputs=reshape(net.content.layers{4}.unitProperties.resp,[datalen,nunits,1])';

switch options.pooling
    case 'ica'
        [Z,B,V]=fastica(v2s_outputs,'approach','symm','numOfIC',nreduced,'g','tanh','epsilon',0.0005,'verbose','on','lastEig',nreduced);
    case 'nonneg'
        % preliminary
        V=nonnegative_pooling(v2s_outputs,nreduced); V=V';
    case 'pca'
        [V,~,~]=pca(v2s_outputs');
        V=V(:,1:nreduced)';
    case 'dimred+whitening'
        [F,~,D]=pca(v2s_outputs');
        V=F(:,1:nreduced)*diag(D(1:nreduced).^(-1/2))*F(:,1:nreduced)';
%         V=V(:,1:nreduced)';
    case 'dimred'
        [F,~,D]=pca(v2s_outputs');
        V=F(:,1:nreduced)*F(:,1:nreduced)';
%         V=V(:,1:nreduced)';
    case 'cov'
        V=cov(v2s_outputs');
    case 'invcov'
        [F,~,D]=pca(v2s_outputs');
        D2=ones(size(D))-D(nreduced)./D;
        D2(D2<0)=0;
        V=F*(diag(D2))*F';
    case 'anti-whitening'
        [F,~,D]=pca(v2s_outputs');
        V=F(:,1:nreduced)*diag(D(1:nreduced))*F(:,1:nreduced)';
%         V=V(:,1:nhigh)';
end;

net.content.layers{5}.weights=reshape(V,[nunits,1,nunits,1]);

net=v1_grating_test(net,12,8,1./(4:12));
net=v1_gratinganalyze(net,'layers',4:5);

net=v2_anzai_test(net,12,4,[1/4 1/6 1/8],6,12);
net=v2_position_invariance(net,'layer',4,'threshold',0.5);
net=v2_position_invariance(net,'layer',5,'threshold',0.5);
net=v2_anzai_analyze1(net,'layer',4);
net=v2_anzai_analyze1(net,'layer',5);

net=strip_resp(net);

netf=net;

nete=net;
nete.content.layers{5}.weights=reshape(F,[nunits,1,nunits,1]);


% displaying

% [V2,idx]=sort(V','descend');
% if options.showneg
%     idx=idx([1:options.nsubshow end-options.nsubshow+1:end],1:nshow);
% else
%     idx=idx([1:options.nsubshow],1:nshow);
% end;
% 
% A2=A(:,idx(:));
% 
% figure; 
% if options.showneg
%     show_v2_filters(v1s_bank,A2,1:numel(idx),0.2,size(idx,1));
% else
%     show_v2_filters(v1s_bank,A2,1:numel(idx),0.2,size(idx,1));
% end;
% 
% figure;
% if options.showneg
%     barh(V2([1:4 end-3:end],1:nshow)');
% else
%     barh(V2(1:4,1:nshow)');
% end;
% axis ij;
% ylim([0 nshow+1]);
% colormap(gray);
% colorbar;
% 
end


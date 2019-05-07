function [netf,nete]=v1_complex_cells(net,dataset,nreduced,varargin)

pr=inputParser;
pr.addParamValue('pooling','dimred',@isstr);
pr.addParamValue('nonlin','half-rect',@isstr);
pr.parse(varargin{:});
options=pr.Results;

[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));
nunits=net.structure.layers{2}.numUnits;

fprintf('building a V1 complex cells model with %d units from %d reduced eigenspace\n',nunits,nreduced);
fprintf('  using %s method and %s nonlinearity \n',options.pooling,options.nonlin);

net.content.layers{2}.layerProperties.nonlin=options.nonlin;
net=addlayer(net,3,'V1c',1,nunits,1);
net.content.layers{3}.layerProperties.meanrem=false;
net.content.layers{3}.layerProperties.nonlin='half-rect';
net.content.layers{3}.layerProperties.dcrem=false;
net=calc_coverage(net);

net=run_net(net,dataset);
v1s_outputs=reshape(net.content.layers{2}.unitProperties.resp,[datalen,nunits,1])';

[F,~,D]=pca(v1s_outputs');

switch options.pooling
    case 'dimred+whitening'
        D(1:nreduced)=D(1:nreduced).^(-1/2);
        D(nreduced+1:end)=0;
    case 'whitening'
        D=D.^(-1/2);
    case 'dimred'
        D(1:nreduced)=1;
        D(nreduced+1:end)=0;
    case 'cov'
        D=D;
    case 'cov2'
        D=D-min(D);
    case 'invcov'
        D=ones(size(D))-D(nreduced)./D;
   case 'denoise'
        D=ones(size(D))-D(nreduced)./D;
        D(D<0)=0;
    case 'dimred+cov'
        D(1:nreduced)=D(1:nreduced);
        D(nreduced+1:end)=0;
end;

V=F*diag(D)*F';

net.content.layers{3}.weights=reshape(V,[nunits,1,nunits,1]);

net=run_net(net,dataset);

net=v1_grating_test(net,12,8,1./(4:12));
net=v1_gratinganalyze(net);

netf=net;

nete=net;
nete.content.layers{3}.weights=reshape(F,[nunits,1,nunits,1]);
nete.content.layers{3}.diagonal_values=D;

end


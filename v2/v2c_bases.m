function v2c_bases(net,varargin)

pr=inputParser;
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',NaN,@isnumeric);
pr.addParamValue('nshow',NaN,@isnumeric);
pr.addParamValue('min_rel_weight',0.05,@isnumeric);
pr.addParamValue('normalize',false,@islogical);
pr.parse(varargin{:});
pr=pr.Results;

V=net.content.layers{5}.weights;
[v2s_nunits,~,v2c_nunits]=size(V);
[V2,idx]=sort(reshape(V,v2s_nunits,v2c_nunits)','descend');

% units_to_show=idx([1:pr.nshow end-pr.nshow+1:end],pr.units);
% weights=V2([1:pr.nshow end-pr.nshow+1:end],pr.units);

units_to_show=idx([1:pr.nshow],pr.units);
weights=V2([1:pr.nshow],pr.units);

figure('position',[0 0 800 600]);
v2s_bases(net,'nodes',1,'units',units_to_show(:),'disp_width',size(units_to_show,1),'onlyextreme',true,'panel_titles','none','axislabels',false,'axisticks',false,'linewidth',0.2);

figure('position',[0 0 200 650]);
barh(weights','FaceColor',[0.2078 0.1647 0.5294]);
axis ij;
set(gca,'FontName','Times','FontSize',12)
set(gca,'YTick',[]);
set(gca,'Box','off');



end

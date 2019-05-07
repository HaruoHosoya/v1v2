function v1c_orientation_tuning(net,varargin)

pr=inputParser;
pr.addParamValue('layer',2,@isnumeric);
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',NaN,@isnumeric);
pr.addParamValue('min_resp',NaN,@isnumeric);
pr.addParamValue('disp_width',25,@isnumeric);
pr.parse(varargin{:});
pr=pr.Results;

info=net.structure.layers{pr.layer};

if isnan(pr.units) pr.units=1:info.numUnits; end;
if isnan(pr.disp_width) pr.disp_width=info.numUnits; end;

disp_height=ceil(length(pr.units)/pr.disp_width);

resp=net.content.layers{pr.layer}.unitProperties3.grating;
[num_phases num_wavelens num_orients num_units num_nodes]=size(resp);


P=1;
for J=1:length(pr.nodes)
    node=pr.nodes(J);
    for K=1:length(pr.units);
        unit=pr.units(K);
        resp1=resp(:,:,:,unit,node);
        resp1=resp1/max(resp1(:));
        resp1=mean(resp1,1);
        respp=permute(resp1,[2 3 1]);
        [~,mi]=max(respp(:));
        [optla,optt]=ind2sub([num_wavelens num_orients],mi);
        x=squeeze(resp1(1,optla,:));
        subplot(disp_height,pr.disp_width,P); 
        plot(0:(180/num_orients):180,[x; x(1)]);
%         title(sprintf('optimal orient=%1.3f freq=%1.3f f1f0=%1.3f',optt,1/optla,f1f0ratio));
%         ylim([0 1]); 
%         ylabel('response');
        xlim([0 180]); 
%         xlabel('orientation'); 
        set(gca,'XTick',[0 90 180]);
        if mod(P,pr.disp_width)~=1 set(gca,'YTick',[]); end;
        if floor((P-1)/pr.disp_width)+1~=disp_height set(gca,'XTick',[]); end;
        set(gca,'FontName','Times','FontSize',12)
        P=P+1;
    end;
end;

end


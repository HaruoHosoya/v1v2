function v1s_bases( net, varargin )

v0_info=net.structure.layers{1};
v1s_info=net.structure.layers{2};
v1s=net.content.layers{2};

pr=inputParser;
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',NaN,@isnumeric);
pr.addParamValue('disp_width',25,@isnumeric);
pr.addParamValue('fitted',false,@islogical);
pr.addParamValue('normalize',false,@islogical);
pr.parse(varargin{:});
pr=pr.Results;

disp_width=pr.disp_width;


width=v0_info.patchWidth;
basis=reshape(v1s.weights,width,width,size(v1s.weights,3),size(v1s.weights,4));
if pr.fitted params=v1s.unitProperties.params; end;

if isnan(pr.units) pr.units=1:v1s_info.numUnits; end;

disp_height=ceil(length(pr.units)/disp_width);

[ix,iy]=ndgrid(1:width,1:width);

% figure('Position',[0 0 800 600]);

for I=1:length(pr.nodes)
    node=pr.nodes(I);
    entire=zeros(disp_width*(width+1)+1,disp_height*(width+1)+1)*-Inf;
    for J=1:length(pr.units)
        unit=pr.units(J);
        if pr.fitted
            img=gabor(ix,iy,params(:,unit,node));
        else
            img=basis(:,:,unit,node);
        end;
        if pr.normalize
            img=img/max(abs(img(:)));
        end;
        [x y]=ind2sub([disp_width disp_height],J);
        x=(x-1)*(width+1)+2; y=(y-1)*(width+1)+2;
        entire(x:x+width-1,y:y+width-1)=img;
    end;
    subplot(length(pr.nodes),1,I);
    imagesc(entire');
%     text(0,-5,['node ' int2str(node)]);
    axis off; axis image; 
    colormap(gray);
end;

end


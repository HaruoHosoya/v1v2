function v1c_bases(net,varargin)

pr=inputParser;
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',NaN,@isnumeric);
pr.addParamValue('disp_width',NaN,@isnumeric);
pr.addParamValue('fitted',false,@islogical);
pr.addParamValue('min_rel_weight',0.05,@isnumeric);
pr.addParamValue('normalize',false,@islogical);
pr.parse(varargin{:});
pr=pr.Results;

v0_info=net.structure.layers{1};
v1s_info=net.structure.layers{2};
v1c_info=net.structure.layers{3};
v1s=net.content.layers{2};
v1c=net.content.layers{3};

if isnan(pr.disp_width) disp_width=floor(v1s_info.numUnits/v1c_info.numUnits); 
else disp_width=pr.disp_width; 
end;

if isnan(pr.units) pr.units=1:v1c_info.numUnits; end;

disp_height=length(pr.units);

width=v0_info.patchWidth;
basis=reshape(v1s.weights,width,width,size(v1s.weights,3),size(v1s.weights,4));
if pr.fitted params=v1s.unitProperties.params; end;

[ix,iy]=ndgrid(1:width,1:width);

% figure('Position',[0 0 800 600]);


subplot(3,4,[1 2 3 5 6 7 9 10 11]);

for J=1:length(pr.nodes)
    node=pr.nodes(J);
    entire=zeros(disp_width*(width+1)+1,disp_height*(width+1)+1)*-Inf;
    for unit=pr.units
        weights=v1c.weights(:,1,unit,node);
        [~,idx]=sort(weights,'descend');
        for I=1:disp_width
            lower=idx(I);
            if weights(lower)<pr.min_rel_weight continue; end;
            if pr.fitted
                img=gabor(ix,iy,params(:,lower,node));
            else
                img=basis(:,:,lower,node);
            end;
            if pr.normalize
                img=img/max(abs(img(:)));
            end;
            x=(I-1)*(width+1)+2; y=(unit-1)*(width+1)+2;
            entire(x:x+width-1,y:y+width-1)=img;
        end;
    end;
%     subplot(1,length(pr.nodes),J);
    imagesc(entire');
%     text(0,-5,['node ' int2str(node)]);
    axis off; axis image; 
    colormap(gray);
end;

subplot(3,4,[4 8 12]);

weights=permute(v1c.weights(:,1,pr.units,pr.nodes(1)),[1,3,2]);
weights=sort(weights,'descend');
weights=weights(1:disp_width,:);

barh(weights','b');
axis ij;
colormap(gray);
set(gca,'FontName','Times','FontSize',12)
set(gca,'YTick',[]);
set(gca,'Box','off');


end


function v1c_bases_ori(net,varargin)

v0_info=net.structure.layers{1};
v1s_info=net.structure.layers{2};
v1c_info=net.structure.layers{3};

pr=inputParser;
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',1:v1c_info.numUnits,@isnumeric);
pr.addParamValue('min_weight',0,@isnumeric);
pr.addParamValue('disp_width',NaN,@isnumeric);
pr.addParamValue('magmeth','eachmax',@isstr);
pr.addParamValue('fast',true,@islogical);
pr.addParamValue('fontsize',NaN,@isnumeric);
pr.addParamValue('axislabels',true,@islogical);
pr.addParamValue('linewidth',1.2,@isnumeric);
pr.addParamValue('barthickness',1,@isnumeric);

pr.addParamValue('min_freq',0,@isnumeric);
pr.addParamValue('numtoshow',1000,@isnumeric);
pr.addParamValue('onlyextreme',false,@islogical);
pr.addParamValue('axisticks',true,@islogical);

pr.addParamValue('fitted',false,@islogical);
pr.addParamValue('fit_type','',@isstr);

pr.parse(varargin{:});
pr=pr.Results;

if length(pr.nodes)==1
    nodes=repmat(pr.nodes,1,length(pr.units));
    units=pr.units;
else
    nodes=pr.nodes;
    units=pr.units;
end;

nunits=length(units);

min_weight=pr.min_weight;

if isnan(pr.disp_width) disp_width=ceil(sqrt(nunits)); 
else disp_width=pr.disp_width; end;

disp_height=double(ceil(nunits/disp_width));

v1c_scale=1;
show_simple=false;
numtoshow=pr.numtoshow;

axes('Position',[0.05 0.95 1 0.1]); 
colormap(red_blue_colormap); 
axis off;

for K=1:nunits
    node=nodes(K);
    unit=units(K);
    if unit==0 continue; end;
    [disp_x disp_y]=ind2sub([disp_width,disp_height],K);
    [subx suby subw subh]=subpos([0.1 0.1 0.8 0.8],disp_width,disp_height,disp_x,disp_y,0.8);
    axes('Position',[subx 1-suby-subh subw subh]);
    if ~isnan(pr.fontsize) set(gca,'fontsize',pr.fontsize); end;
    axis ij; 
    axis square;
    weights=net.content.layers{3}.weights(:,1,unit,node);
    params=net.content.layers{2}.unitProperties.params(:,:,node);
    CX=params(1,:); CY=params(2,:); TH=params(6,:); W=weights(:);
    SX=params(4,:); SY=params(5,:); FR=params(7,:);
    if pr.onlyextreme [CX,CY,TH,W,SX,SY,FR]=retain_onlyextreme(CX,CY,TH,W,SX,SY,FR); end;
    
    [S_sorted idx]=sort(abs(W(:)),'descend');
    idx=idx(~isnan(S_sorted)); idx=idx(1:min(length(idx),numtoshow)); idx=idx(end:-1:1);
    CX=CX(idx); CY=CY(idx); TH=TH(idx); W=W(idx); SX=SX(idx); SY=SY(idx); FR=FR(idx); 
    mag=1/(max(W(:))+eps);
    W=W*mag;
    if isempty(CX) continue; end;

    grid=net.content.layers{3}.nodeProperties.cover(node,:);
    
    linex=[]; liney=[]; cdata=[];
    for J=1:length(CX)
        if (isnan(CX(J)) ||isnan(FR(J)) || FR(J)<pr.min_freq || abs(W(J))<min_weight) continue; end;
        len=SY(J)*v1c_scale;
        wid=pr.barthickness*v1c_scale;
%         len=SY(J)*v1c_scale;
%         wid=1/FR(J)*v1c_scale/8*pr.barthickness;
%         [xs ys]=oriented_ellipse(CX(J),CY(J),TH(J),len,wid);
        [xs ys]=oriented_rectangle(CX(J),CY(J),TH(J),len,wid);
        linex=[linex xs]; liney=[liney ys];
        cdata=[cdata;ones(length(xs),1)*W(J)];
    end;
    xlim([grid(1) grid(1)+grid(3)]); ylim([grid(2) grid(2)+grid(4)]);
    if ~pr.fast
        for L=1:3:length(linex)
            line([linex(L) linex(L+1)],[liney(L) liney(L+1)],'LineWidth',pr.linewidth,'Color',cdata(L,:));  
        end;
    else
        p=patch(linex,liney,0);
%         set(p,'FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','flat','FaceColor','none','LineWidth',2);
        set(p,'FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','none','FaceColor','flat','LineWidth',2);
        caxis([-1 1]);
    end;
    if disp_x~=1 set(gca,'YTick',[]); end;
    if disp_y~=disp_height set(gca,'XTick',[]); end;
    if ~pr.axisticks
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
%         axis off;
    end;    
    set(gca,'FontName','Times','FontSize',12)
    axis square;
end;


end

function [xs ys]=oriented_rectangle(cx,cy,th,len,wid)

xl=sin(th)*len; xw=cos(th)*wid;
yl=-cos(th)*len; yw=sin(th)*wid;

xs=[xl+xw; xl-xw; -xl-xw; -xl+xw]+cx;
ys=[yl+yw; yl-yw; -yl-yw; -yl+yw]+cy;

end

function [xs ys]=oriented_ellipse(cx,cy,th,len,wid)

ori=th+pi/2;
mat=[cos(ori) -sin(ori); sin(ori) cos(ori)];

a=0:pi/12:2*pi;
ps=mat*[cos(a)*len; sin(a)*wid];
xs=ps(1,:)'+cx;
ys=ps(2,:)'+cy;

end

function [CX,CY,TH,W,SX,SY,FR]=retain_onlyextreme(CX,CY,TH,W,SX,SY,FR)
    maxx=max(CX); XY=CX+maxx*CY;
    XYuniq=unique(XY);
    for I=1:length(XYuniq)
        idx_pos=find(XY==XYuniq(I));
        [val_max,idx_max]=max(W(idx_pos));
        [val_min,idx_min]=min(W(idx_pos));
        retain=[];
        if val_max>0 retain=[retain;idx_pos(idx_max)]; end;
        if val_min<0 retain=[retain;idx_pos(idx_min)]; end;
        W(setdiff(idx_pos,retain))=NaN;
    end;
end

function relbase=fitted_relbase(desc_type,relbase,params)
    CX=relbase(1,:); CY=relbase(2,:); TH=relbase(3,:); W=relbase(4,:);
    SX=relbase(5,:); SY=relbase(6,:); FR=relbase(7,:); F1F0=relbase(8,:);
    newW=v2basis_desc(CX,CY,TH,FR,params,desc_type);
    relbase=[CX; CY; TH; newW; SX; SY; FR; F1F0];
end

function cmap=red_blue_colormap

cmap=[];
for v=-1:0.01:-0.01
    cmap=[cmap; hsv2rgb([0.66 -v 1])];
end;
for v=0:0.01:1
    cmap=[cmap; hsv2rgb([0 v 1])];
end;


end


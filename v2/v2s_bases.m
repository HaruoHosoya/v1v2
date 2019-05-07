function h=v2s_bases(net,varargin)

v0_info=net.structure.layers{1};
v1c_info=net.structure.layers{3};
v2s_info=net.structure.layers{4};

pr=inputParser;
pr.addParamValue('nodes',1,@isnumeric);
pr.addParamValue('units',1:v2s_info.numUnits,@isnumeric);
pr.addParamValue('show_3d',false,@islogical);
pr.addParamValue('min_weight',-Inf,@isnumeric);
pr.addParamValue('mag',0.5,@isnumeric);
pr.addParamValue('disp_width',NaN,@isnumeric);
pr.addParamValue('win_title','',@isstr);
pr.addParamValue('panel_titles','node-unit');
pr.addParamValue('tight',true,@islogical);
pr.addParamValue('origin',NaN);
pr.addParamValue('magmeth','eachmax',@isstr);
pr.addParamValue('fast',true,@islogical);
pr.addParamValue('deg_per_pixel',NaN,@isnumeric);
pr.addParamValue('fontsize',NaN,@isnumeric);
pr.addParamValue('axislabels',true,@islogical);
pr.addParamValue('linewidth',2,@isnumeric);
pr.addParamValue('barthickness',1,@isnumeric);

pr.addParamValue('min_freq',0,@isnumeric);
pr.addParamValue('numtoshow',1000,@isnumeric);
pr.addParamValue('onlyextreme',false,@islogical);
pr.addParamValue('axisticks',true,@islogical);

pr.addParamValue('fitted',false,@islogical);
pr.addParamValue('fit_type','',@isstr);

pr.parse(varargin{:});
pr=pr.Results;

mag=pr.mag;
win_title=pr.win_title;
panel_titles=pr.panel_titles;

if length(pr.nodes)==1
    nodes=repmat(pr.nodes,1,length(pr.units));
    units=pr.units;
else
    nodes=pr.nodes;
    units=pr.units;
end;

nunits=length(units);

min_weight=pr.min_weight;

if isnan(pr.disp_width) disp_width=floor(1.5*sqrt(nunits)); 
else disp_width=pr.disp_width; end;

if isnan(pr.origin) x0=1; y0=1; 
else x0=pr.origin(1); y0=pr.origin(2); end;

if isnan(pr.deg_per_pixel) metric='pixel'; deg_per_pixel=1;
else x0=0; y0=0; metric='deg'; deg_per_pixel=pr.deg_per_pixel; end;


disp_height=double(ceil(nunits/disp_width));

v1c_scale=1;
show_simple=false;
numtoshow=pr.numtoshow;

axes('Position',[0.05 0.95 1 0.1]); 
text(0,0,win_title);
colormap(red_blue_colormap); 
if ~isnan(pr.fontsize) set(gca,'fontsize',pr.fontsize); end;
axis off;

v2s_relbases=net.content.layers{4}.unitProperties2.relbases;

magmeth=pr.magmeth;
if isempty(magmeth)
    % magmeth='eachmax';
    % magmeth='eachabsmax';
    magmeth='allmax';
end;

%h=figure('position',[0 0 1400 800]);
for K=1:nunits
    node=nodes(K);
    unit=units(K);
    if unit==0 continue; end;
    [disp_x disp_y]=ind2sub([disp_width,disp_height],K);
    [subx suby subw subh]=subpos([0.1 0.1 0.8 0.9],disp_width,disp_height,disp_x,disp_y,0.8);
    axes('Position',[subx 1-suby-subh subw subh]);
    if ~isnan(pr.fontsize) set(gca,'fontsize',pr.fontsize); end;
    %axis([1 grid_width 1 grid_width]);
    axis ij; 
    axis square;
    relbase=v2s_relbases(:,:,unit,node);
    if pr.fitted 
        desc_type=net.content.layers{4}.unitProperties.desc_types(1,unit,node);
        params=net.content.layers{4}.unitProperties2.desc_params(:,desc_type,unit,node);
        relbase=fitted_relbase(desc_type,relbase,params); 
    end;
    if pr.onlyextreme relbase=retain_onlyextreme(relbase); end;
    CX=relbase(1,:); CY=relbase(2,:); TH=relbase(3,:); W=relbase(4,:);
    SX=relbase(5,:); SY=relbase(6,:); FR=relbase(7,:); F1F0=relbase(8,:);
    
    [S_sorted idx]=sort(abs(W(:)),'descend');
    idx=idx(~isnan(S_sorted)); idx=idx(1:min(length(idx),numtoshow)); idx=idx(end:-1:1);
    CX=CX(idx); CY=CY(idx); TH=TH(idx); W=W(idx); SX=SX(idx); SY=SY(idx); FR=FR(idx); F1F0=F1F0(idx);
    CX=CX+x0-1; CY=CY+y0-1;
    CX=CX*deg_per_pixel; CY=CY*deg_per_pixel; SX=SX*deg_per_pixel; SY=SY*deg_per_pixel;
    switch(magmeth)
        case 'eachabsmax'
            mag=1/(max(abs(W(:)))+eps);
        case 'eachmax'
            mag=1/(max(W(:))+eps);
        case 'allmax'
            Sall=v2s_relbases(4,:,pr.units,node);
            mag=1/(max(Sall(:))+eps);
    end;
    if isempty(CX) continue; end;

    if pr.tight
        grid=net.content.layers{4}.nodeProperties.cover(node,:);
        grid(1)=grid(1)-1; grid(2)=grid(2)-1;
    else
        grid=[0 0 v0_info.width v0_info.width];
    end;
    grid=grid*deg_per_pixel;
    
    linex=[]; liney=[]; cdata=[];
    for J=1:length(CX)
        if (isnan(CX(J)) ||isnan(FR(J)) || FR(J)<pr.min_freq || W(J)*mag<min_weight || (~show_simple && F1F0(J)>1)) continue; end;
        if (W(J)>0)
            col=hsv2rgb([0 restrict(W(J)*mag,0,1) 1]);
        else
            col=hsv2rgb([0.66 restrict(-W(J)*mag,0,1) 1]);
        end;
        if (pr.show_3d)
            scatter3(CX(J),CY(J),mod(TH(J),pi),50,col,'filled');
            xlim([x0+grid(1) x0+grid(1)+grid(3)]); ylim([y0+grid(2) y0+grid(2)+grid(4)-1]); zlim([0 pi]); 
            hold on;
        else
%             len=SY(J)*v1c_scale;
%             linex=[linex CX(J)+sin(TH(J))*len CX(J)-sin(TH(J))*len nan];
%             liney=[liney CY(J)-cos(TH(J))*len CY(J)+cos(TH(J))*len nan];
%             cdata=[cdata; col; col; nan nan nan];
            len=SY(J)*v1c_scale;
            wid=1/FR(J)*v1c_scale/8*pr.barthickness;
%             [xs ys]=oriented_rectangle(CX(J),CY(J),TH(J),len,wid);
            [xs ys]=oriented_ellipse(CX(J),CY(J),TH(J),len,wid);
            linex=[linex xs]; liney=[liney ys];
%             cdata=cat(1,cdata,reshape(col,1,1,3));
%             cdata1=reshape(col,1,3);cdata=[cdata;repmat(cdata1,length(xs),1)];
            cdata=[cdata;ones(length(xs),1)*W(J)*mag];
        end;
    end;
    if (~pr.show_3d)
        xlim([x0+grid(1) x0+grid(1)+grid(3)]); ylim([y0+grid(2) y0+grid(2)+grid(4)]);
        if ~pr.fast
            for L=1:3:length(linex)
                line([linex(L) linex(L+1)],[liney(L) liney(L+1)],'LineWidth',pr.linewidth,'Color',cdata(L,:));  
            end;
        else
            p=patch(linex,liney,0);
%             set(p,'FaceVertexCData',cdata,'CDataMapping','direct','EdgeColor','flat','FaceColor','none','LineWidth',2);
%             set(p,'CData',cdata,'CDataMapping','direct','FaceColor','flat','EdgeColor','none','LineWidth',2);
%             set(p,'FaceVertexCData',cdata,'CDataMapping','direct','EdgeColor','flat','FaceColor','none','LineWidth',2);
            set(p,'FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','flat','FaceColor','none','LineWidth',pr.linewidth);
            caxis([-1 1]);
        end;
    end;
    if pr.axislabels
        xlabel(['x [' metric ']']);
        ylabel(['y [' metric ']']);
    end;
    if disp_x~=1 set(gca,'YTick',[]); end;
    if disp_y~=disp_height set(gca,'XTick',[]); end;
    set(gca,'FontName','Times','FontSize',12)
    axis square;
    switch panel_titles
        case 'none'
        case 'node-unit'
            title(sprintf('(%d,%d)',node,unit));
        case 'unit'
            title(sprintf('%d',unit));
        otherwise
            if iscell(panel_titles)
                title(texlabel(panel_titles{K}));
            end;
    end;
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

function relbase=retain_onlyextreme(relbase)
    CX=relbase(1,:); CY=relbase(2,:); TH=relbase(3,:); W=relbase(4,:);
    SX=relbase(5,:); SY=relbase(6,:); FR=relbase(7,:); F1F0=relbase(8,:);
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
    relbase=[CX; CY; TH; W; SX; SY; FR; F1F0];
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


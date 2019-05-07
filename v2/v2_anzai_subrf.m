function [] = v2_anzai_subrf( net, varargin )

parser=inputParser;
parser.addParamValue('tol',0,@(x)x>=0);
parser.addParamValue('layer',4,@(x)x>=1);
parser.addParamValue('nodes',1:1,@(x)all(x>=1));
parser.addParamValue('units',1:10,@(x)all(x>=1));
parser.addParamValue('minresp',0.1,@(x)x>=0);
parser.addParamValue('minmaxresp',0.0,@(x)x>=0);
parser.addParamValue('entire',true,@islogical);
parser.addParamValue('mag',1,@(x)x>=0);
parser.addParamValue('disp_width',NaN,@(x)x>=1);
parser.addParamValue('normalize',true,@islogical);
parser.addParamValue('linewidth',0.5,@isnumeric);
parser.parse(varargin{:});
options=parser.Results;

L=options.layer;
v2resp=net.content.layers{L}.unitProperties.subrf;
[n,topUnits,topNodes]=size(v2resp);

if isfield(net.content.layers{L}.unitProperties,'sign')
    resp_sign=net.content.layers{L}.unitProperties.sign;
else
    resp_sign=ones(topUnits,topNodes);
end;
    
if isfield(net.content.layers{L}.layerProperties,'subrf_dim')
    subrf_dim=net.content.layers{L}.layerProperties.subrf_dim;
    ny=subrf_dim(1); nx=subrf_dim(2);
    numOrient=subrf_dim(3); numFreqs=subrf_dim(4); numPhases=subrf_dim(5);
    v2resp=reshape(v2resp,[numPhases*numFreqs numOrient nx ny topUnits topNodes]);
    v2resp=bsxfun(@times,v2resp,reshape(resp_sign,[1,1,1,1,topUnits,topNodes]));
    v2resp=reshape(max(v2resp,[],1),[numOrient nx ny topUnits topNodes]);
%     v2resp=reshape(v2resp,[numPhases numFreqs numOrient nx ny topUnits topNodes]);
%     v2resp=bsxfun(@times,v2resp,reshape(resp_sign,[1,1,1,1,1,topUnits,topNodes]));
%     v2resp=shiftdim(max(v2resp,[],1),1);
%     v2resp2=zeros(numOrient,nx,ny,topUnits,topNodes);
%     for T=1:topNodes
%         for U=1:topUnits
%             [~,optfreqi]=max(max(reshape(v2resp(:,:,:,:,U,T),[numFreqs numOrient*nx*ny]),[],2));
%             v2resp2(:,:,:,U,T)=v2resp(optfreqi,:,:,:,U,T);
%         end;
%     end;
%     v2resp=v2resp2;
else    
    numOrient=12;
    % numOrient=8;
    numLocation=n/numOrient;
    ny=floor(sqrt(numLocation));
    nx=ceil(numLocation/ny);
    v2resp=bsxfun(@times,v2resp,reshape(resp_sign,[1,1,1,topUnits,topNodes]));
    v2resp=reshape(v2resp,[numOrient,ny,nx,topUnits,topNodes]);
end;
v2resp=v2resp*options.mag;

% relcover=net.content.layers{L}.nodeProperties.relcover;

num=length(options.nodes)*length(options.units);
if isnan(options.disp_width) disp_width=floor(sqrt(num));
else disp_width=options.disp_width; end;
disp_height=ceil(num/disp_width);

J=0;
for T=options.nodes
    if options.entire
        cx=1; cy=1; cw=ny; ch=nx;
    else
        [cx,cy,cw,ch]=calc_v1cover(net,L,T);
%         cx=relcover(T,1); cy=relcover(T,2); cw=relcover(T,3); ch=relcover(T,4);
    end;
    [lx,ly]=ndgrid(cx:cx+cw-1,cy:cy+ch-1);
    for U=options.units
        J=J+1;
        subplot(disp_height,disp_width,J); 
        title(sprintf('unit(%d,%d)', T, U));
        r2 = v2resp(:,lx,ly,U,T);
        mr = max(r2(:));
        pos=get(gca,'Position');
        unix=pos(1); uniy=pos(2); uniw=pos(3); unih=pos(4);
        for I=1:cw*ch
          [IX,IY]=ind2sub([cw ch],I);
          X=(IX-1)/cw; 
          Y=1-IY/ch;
          axes('Position',[unix+uniw*X,uniy+unih*Y,uniw/cw,unih/ch]);
          t=(pi/numOrient) * (0:(numOrient-1));
          r=v2resp(:,lx(I),ly(I),U,T);
%           r(r<mr*options.tol | r<options.minresp)=0;
          if options.normalize r=r/mr; end;
          polar(-t(:)'-pi/2,r(:)',[0,1],options.linewidth);
        end;
    end;
end;

set(findobj(gcf, 'type','axes'), 'Visible','off')

end

function [cx,cy,cw,ch]=calc_v1cover(net,layer,node)

relcover=net.content.layers{layer}.nodeProperties.relcover;
cx=relcover(node,1); cy=relcover(node,2); cw=relcover(node,3); ch=relcover(node,4);

if layer==3 return; end;

width=net.structure.layers{layer-1}.width;

node1=sub2ind([width,width],cx,cy);
[cx1,cy1,cw1,ch1]=calc_v1cover(net,layer-1,node1);
node2=sub2ind([width,width],cx+cw-1,cy+ch-1);
[cx2,cy2,cw2,ch2]=calc_v1cover(net,layer-1,node2);

cx=cx1; cy=cy1; cw=cx2+cw2-cx1; ch=cy2+ch2-cy1;

end

function polar(theta,rho,radial_limits,linewidth)

th = 0:pi/20:2*pi;
xunit = cos(th);
yunit = sin(th);

rmax = radial_limits(2);
rmin = radial_limits(1);

patch('xdata',xunit*(rmax-rmin),'ydata',yunit*(rmax-rmin),'facecolor','white','edgecolor','black','LineWidth',linewidth)

% hhh=line(xunit*(rmax-rmin),yunit*(rmax-rmin),'linestyle','-');

% nspokes=2;
% th = (1:nspokes)*pi/nspokes;
% cst = cos(th); snt = sin(th);
% cs = [-cst; cst];
% sn = [-snt; snt];
% line((rmax-rmin)*cs,(rmax-rmin)*sn,'linestyle',':','linewidth',1,'color','black');

% % annotate spokes in degrees
%     rt = 1.1*(rmax-rmin);
%     for i = 1:length(th)
%         %text(rt*cst(i),rt*snt(i),int2str(i*180/nspokes),...
%         %     'horizontalalignment','center',...
%         %     'handlevisibility','off','parent',cax);
%         if i == length(th)
%             loc = int2str(0);
%         else
%             loc = int2str(180+i*180/nspokes);
%         end
%         %text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
%         %     'handlevisibility','off','parent',cax)
%     end


xx=min(rmax,rho-rmin).*cos(theta);
yy=min(rmax,rho-rmin).*sin(theta);

xs=[-xx xx -xx(1)];
ys=[-yy yy -yy(1)];

line(xs,ys,'color','r','LineWidth',linewidth);

% rho1=rho(rho-rmin>0);
% rho2=rho(rho-rmin<=0);
% theta1=theta(rho-rmin>0);
% theta2=theta(rho-rmin<=0);
% 
% xx1 = (rho1 - rmin).*cos(theta1);
% yy1 = (rho1 - rmin).*sin(theta1);
% xx2 = (rho2 - rmin).*cos(theta2);
% yy2 = (rho2 - rmin).*sin(theta2);
% 
% xs1=[-xx1 xx1 -xx1(1)];
% ys1=[-yy1 yy1 yy1(1)];
% 
% % xs1=[-xx1; xx1];
% % xs2=[-xx2; xx2];
% % ys1=[-yy1; yy1];
% % ys2=[-yy2; yy2];
% 
% q1=line(xs1,ys1,'color','r','LineWidth',2.0);
% 


end

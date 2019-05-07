function net = v2_position_invariance(net,varargin)
% subreceptive field analysis

parser=inputParser;
parser.addParamValue('layer',4,@isnumeric);
parser.addParamValue('normalize',true,@islogical);
parser.addParamValue('entire',true,@islogical);
parser.addParamValue('threshold',0.5,@isnumeric);
parser.parse(varargin{:});
options=parser.Results;

L=options.layer;

v2resp=net.content.layers{L}.unitProperties.subrf;
[numStimuli topUnits topNodes]=size(v2resp);

subrf_dim=num2cell(net.content.layers{L}.layerProperties.subrf_dim);
[ny,nx,numOrient,numFreqs,numPhases]=subrf_dim{:};
v2resp=reshape(v2resp,[numPhases*numFreqs numOrient nx ny topUnits topNodes]);
v2resp=reshape(max(v2resp,[],1),[numOrient nx ny topUnits topNodes]);

invariance_index=zeros(1,topUnits,topNodes);

for T=1:topNodes
    if options.entire
        x0=1; y0=1; w=ny; h=nx;
    else
        [x0,y0,w,h]=calc_v1cover(net,L,T);
    end;
    [lx,ly]=ndgrid(x0:x0+w-1,y0:y0+h-1);
    for U=1:topUnits
        [~,idx]=max(flatten(v2resp(:,:,:,U,T)));
        [mo mx my]=ind2sub([numOrient nx ny],idx);
        resp1=shiftdim(v2resp(mo,:,:,U,T),1);
        resp1=resp1/max(resp1(:))>options.threshold;
        minx=find(max(resp1,[],2)==1,1,'first');
        maxx=find(max(resp1,[],2)==1,1,'last');
        miny=find(max(resp1,[],1)==1,1,'first');
        maxy=find(max(resp1,[],1)==1,1,'last');
        invariance_index(1,U,T)=max(maxx-minx+1,maxy-miny+1);
    end;
end;

net.content.layers{L}.unitProperties.invariance_index=invariance_index;

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


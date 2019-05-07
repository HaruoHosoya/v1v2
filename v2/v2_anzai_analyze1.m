function net = v2_anzai_analyze1(net,varargin)
% subreceptive field analysis

parser=inputParser;
parser.addParamValue('layer',4,@isnumeric);
parser.addParamValue('normalize',true,@islogical);
parser.addParamValue('entire',true,@islogical);
parser.parse(varargin{:});
options=parser.Results;

L=options.layer;

v2resp=net.content.layers{L}.unitProperties.subrf;
[numStimuli topUnits topNodes]=size(v2resp);

subrf_dim=num2cell(net.content.layers{L}.layerProperties.subrf_dim);
[ny,nx,numOrient,numFreqs,numPhases]=subrf_dim{:};
v2resp=reshape(v2resp,[numPhases*numFreqs numOrient nx ny topUnits topNodes]);
v2resp=reshape(max(v2resp,[],1),[numOrient nx ny topUnits topNodes]);

for T=1:topNodes
    fprintf('Node %d: ',T);
    if options.entire
        x0=1; y0=1; w=ny; h=nx;
    else
        [x0,y0,w,h]=calc_v1cover(net,L,T);
    end;
    [lx,ly]=ndgrid(x0:x0+w-1,y0:y0+h-1);
    for U=1:topUnits
        mr=max(flatten(v2resp(:,lx(:),ly(:),U,T)));
        fprintf('%d ',U);
        parfor v=1:length(lx(:))
            r = v2resp(:,lx(v),ly(v),U,T);
            if options.normalize r=max(0,r/mr); end;
            [c1,gof1]=vonmisesFit(r(:));
            [c2,gof2]=vonmisesFit2(r(:));
            vmparams(:,v,U,T)=...
                [gof1.rsquare, c1.a1, c1.m1, c1.s1, ...
                gof2.rsquare, c2.a1, c2.m1, c2.s1, c2.a2, c2.m2, c2.s2];
        end;
    end;
    fprintf('\n');
end;

net.content.layers{L}.unitProperties2.vmparams=vmparams;

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


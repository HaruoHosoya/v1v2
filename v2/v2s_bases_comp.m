function net=v2s_bases_comp(net)
% compute v2s bases with respect v1c bases
% data=v2s_bases_comp(data)
%  data : network data

v0_info=net.structure.layers{1};
v1s_info=net.structure.layers{2};
v1c_info=net.structure.layers{3};
v2s_info=net.structure.layers{4};

v0=net.content.layers{1};
v1s=net.content.layers{2};
v1c=net.content.layers{3};
v2s=net.content.layers{4};

v0_num_nodes=v0_info.width^2;
v1s_num_nodes=v1s_info.width^2;
v1c_num_nodes=v1c_info.width^2;
v2s_num_nodes=v2s_info.width^2;

v0_step=(v0_info.width-v0_info.patchWidth)/(v1s_info.width-1);

if isfield(v1c.unitProperties,'optwlen') optwlen=v1c.unitProperties.optwlen;
else optwlen=NaN; end;
if isfield(v1c.unitProperties,'f1f0ratio') f1f0ratio=v1c.unitProperties.f1f0ratio;
else f1f0ratio=NaN; end;

v1c_params=v1c.unitProperties.params;

v2s_relbases=zeros(8,v1c_info.patchWidth^2*v1c_info.numUnits,v2s_info.numUnits,v2s_num_nodes);

for node=1:v2s_num_nodes
    for unit=1:v2s_info.numUnits
        [CX CY CZ W SX SY FR F1F0]=v2s_repr(node,unit);
        v2s_relbases(:,:,unit,node)=[CX CY CZ W SX SY FR F1F0]';
    end;
end;

net.content.layers{4}.unitProperties2.relbases=v2s_relbases;



    function [CX CY TH W SX SY FR F1F0]=v2s_repr(node,unit)
        CX=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        CY=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        TH=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        W=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        SX=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        SY=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        FR=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        F1F0=zeros(v1c_info.numUnits,v1c_info.patchWidth^2)*NaN;
        for N1cp=1:v1c_info.patchWidth^2
            w=v2s.weights(:,N1cp,unit,node);
            [~,U1c_range]=sort(w,'ascend');
            N1c=v1c_patch_to_node(node,N1cp);
%             U1c_range=U1c_range(~isnan(squeeze(optwlen(N1c,U1c_range))));
            cx=v1c_params(1,:,N1c);
            cy=v1c_params(2,:,N1c);
            sx=v1c_params(3,:,N1c);
            sy=v1c_params(4,:,N1c);
            th=v1c_params(5,:,N1c);
            fr=v1c_params(6,:,N1c);
            CX(:,N1cp)=cx(U1c_range);
            CY(:,N1cp)=cy(U1c_range);
            TH(:,N1cp)=th(U1c_range);
            W(:,N1cp)=w(U1c_range);
            SX(:,N1cp)=sx(U1c_range);
            SY(:,N1cp)=sy(U1c_range);
%             if ~isnan(optwlen) WL(:,N1cp)=optwlen(1,U1c_range,N1c); 
%             else WL(:,N1cp)=1; end;
            FR(:,N1cp)=fr(U1c_range);
            if ~isnan(f1f0ratio) F1F0(:,N1cp)=f1f0ratio(1,U1c_range,N1c); 
            else F1F0(:,N1cp)=0; end;
        end;
        CX=CX(:); CY=CY(:); TH=TH(:); W=W(:); SX=SX(:); SY=SY(:); FR=FR(:); F1F0=F1F0(:);
    end
    
    function v1c_node=v1c_patch_to_node(v2_node,v1c_patch_node)
        pos=v2s.nodeProperties.relcover(v2_node,:);
        v1cx=pos(1); v1cy=pos(2);
        [v1cpx v1cpy]=ind2sub([v1c_info.patchWidth v1c_info.patchWidth],v1c_patch_node);
        v1c_node=sub2ind([v1c_info.width v1c_info.width],v1cx+v1cpx-1,v1cy+v1cpy-1);
    end


end

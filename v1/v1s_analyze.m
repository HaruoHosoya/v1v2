function net = v1s_analyze(net)
% net = v1s_analyze(net)
%   net : network data
% fits v1s bases with Gabor function
% writes unit props "params" and "residual"

basis=getV1Bases(net);
v1nodes=net.structure.layers{2}.width^2;
v1units=net.structure.layers{2}.numUnits;

all_params=zeros(8,v1units,v1nodes);
all_residual=zeros(1,v1units,v1nodes);

for N=1:v1nodes
    fprintf('%d: ',N);
    parfor U=1:v1units
        fprintf('%d ',U);
        I=(N-1)*v1units+U;
        [params,residual]=gaborFit(basis(:,:,I));
        if isempty(params)
            fprintf('cannot fit unit (%d,%d); replaced with random basis\n',N,U);
            [params,residual]=gaborFit(rand(size(basis,1),size(basis,2))*0.0001);
        end;
        all_params(:,U,N)=params;
        all_residual(1,U,N)=residual;
    end;
    fprintf('\n');
end;

net.content.layers{2}.unitProperties.params=all_params;
net.content.layers{2}.unitProperties.residual=all_residual;

end


      
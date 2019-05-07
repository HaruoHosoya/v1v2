function [A,YS]=energy_component_model(dataset,bank,nunits,varargin)

pr=inputParser;
pr.addParamValue('display',true,@islogical);
pr.addParamValue('nolearning',false,@islogical);
pr.addParamValue('nsteps',[],@isnumeric);
pr.addParamValue('params',[]);
pr.addParamValue('filtered',false);
pr.parse(varargin{:});
options=pr.Results;

[~,ndata]=size(dataset);
[~,bankdim]=size(bank);

if options.filtered
    XS=dataset;
else
    XS=bank'*dataset;
end

if ~isempty(options.params)
    A=options.params;
else
    A=randn(bankdim,nunits)*0.1;   
    A=A./repmat(sqrt(sum(A.^2,1)),bankdim,1);
end;

%%%%%

nsteps=10000;
minibatch_size=100;
display_steps=1000;

alpha=1;
wdim=nunits;

dA=zeros(bankdim,nunits);

YS=zeros(nunits,nsteps);

if options.display h=figure; end;
if options.nolearning beta=0; alpha=0; end;
if ~isempty(options.nsteps) nsteps=options.nsteps; end;

f=@(x)log(cosh(x));
g=@(x)tanh(x);
gd=@(x)1-tanh(x).^2;

% f=@(u)-exp(-u.^2/2);
% g=@(u)u.*exp(-u.^2/2);
% gd=@(u)(1-u.^2).*exp(-u.^2/2);

for T=1:nsteps
    if mod(T,minibatch_size)==1 
        A=A+dA/minibatch_size;
%         A=A*(A'*A)^(-1/2);
        A=A./repmat(sqrt(sum(A.^2,1)),size(A,1),1)*20;
        dA(:)=0;
    end;
    if mod(T,display_steps)==1 && options.display 
        fprintf('step #%d\n',T);
        disp_weight(h,A,bank,[]); 
    end;
%     if mod(T,nsteps/8)==0
%         alpha=alpha/2;
%         fprintf('alpha=%e\n',alpha);
%     end;

    % inference

    X=XS(:,mod(T-1,ndata)+1).^2;
    Y=X'*A;
    
    dA=dA-alpha*X*g(X'*A);

    YS(:,T)=Y;
    
    
end;

end


function disp_weight(h,A,bank,LS)
    figure(h); clf;
    subplot(2,3,[1 4]);
    [~,nunits]=size(A);
    [sz,~]=size(bank);
    wid=sqrt(sz);

    nsub=5;
    M=ones(nsub*(wid+1),nunits*(wid+1));
    for I=1:nunits
        [~,idx]=sort(A(:,I),'descend');
        for J=1:nsub
            x=(J-1)*(wid+1)+1; y=(I-1)*(wid+1)+1;
%             M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(J))/max(bank(:,idx(J))),wid,wid);
            M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(J)),wid,wid);
        end;
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(2,3,[2 3]);
    imagesc(A);
    colorbar;

    subplot(2,3,[5 6]);
    plot(LS);

end
                  
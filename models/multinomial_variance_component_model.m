function multinomial_variance_component_model(dataset,bank)

% wid=12;
% ndata=10000;
% dataset=create_image_dataset(wid,ndata);

[datadim,ndata]=size(dataset);
wid=sqrt(datadim);
nunits=30;

%%%%%

fprintf('estimating filter bank...\n');
% [bank,~]=fastica(dataset,'approach','symm','numOfIC',wid^2,'g','tanh','epsilon',0.0001,'verbose','on');

%%%%%

bank=gabor_bank(wid,wid/3,wid/4:wid/4:wid-wid/4,[1/8 1/6 1/4],0:pi/8:pi-pi/8,0:pi/2:2*pi-pi/2);

bank=bank./repmat(sqrt(sum(bank.^2,1)),datadim,1);
bankdim=size(bank,2);

XS=bank'*dataset;
meanvar=mean(var(XS,[],2));
W=randn(bankdim,nunits)*0.01+0.5;

%%%%%

nsteps=1000000;
alpha=0.01;
dW=zeros(bankdim,nunits);

h=figure;
ks=zeros(nunits);

for T=1:nsteps
    if mod(T,100)==1 
        fprintf('step #%d %1.3f\n',T,norm(dW)); 
        disp_weight(h,W,bank,ks);
        W=W+alpha*dW;
        dW=zeros(bankdim,nunits);
    end;
    X=XS(:,mod(T,ndata)+1);
    
    % inference

    q=zeros(nunits,1);
    for I=1:nunits
        % gaussian
        q(I)=-1/2*sum(log(W(:,I))+X.^2./W(:,I));
        
        % laplacian
%         q(I)=q(I)-sum(log(W(:,I))+abs(X(:))./W(:,I));
    end;
%     p=exp(q);
    q=real(q); p=exp(q-max(q)); p=real(p); p(p<0)=0;
    if all(p==0) k=randi(nunits,1); 
    else p=p/sum(p); k=randp(p);
    end;
    ks(k)=ks(k)+1;
%     disp(p'); disp(k);

    % update

    % gaussian
%     dW(:,k)=dW(:,k)+1/2.*(X.^2-W(:,k))./W(:,k).^2;
    dW(:,k)=dW(:,k)+1/2.*(X.^2-W(:,k));
    
%     dW(:,k)=dW(:,k)+1/2*(X(:).^2./W(:,k).^2 - 1./sum(W(:,k)));
%     dW(:,k)=dW(:,k)+1/2*(X(:).^2. - W(:,k).^2./sum(W(:,k)));
%     dW(:,k)=dW(:,k)+1/2*(X(:).^2./W(:,k).^2 - 1./sum(W(:,k)));
    
    % laplacian
%     dW(:,k)=dW(:,k)+1/2*(abs(X(:))./W(:,k).^2 - 1./sum(W(:,k)));
%     dW(:,k)=dW(:,k)+1/2*(abs(X(:)) - W(:,k).^2./sum(W(:,k)));

    
end;



end

function bank=gabor_bank(fieldwid,gaborwid,grid,freq,ori,pha)

len=length(grid)^2*length(freq)*length(ori)*length(pha);
bank=zeros(fieldwid,fieldwid,len);
U=1;
[IX IY]=meshgrid(1:fieldwid,1:fieldwid);
for Y=grid
    for X=grid
        for F=freq
            for O=ori
                for P=pha
                    g=gabor(IX,IY,[X Y 1 gaborwid/2 gaborwid O F P]);
                    bank(:,:,U)=g;
                    U=U+1;
                end;
            end;
        end;
    end;
end;

bank=reshape(bank,fieldwid*fieldwid,len);

end


function disp_weight(h,W,bank,ks)
    figure(h); clf;
    subplot(131);
    [~,nunits]=size(W);
    [sz,~]=size(bank);
    wid=sqrt(sz);

    nsub=5;
    M=ones(nsub*(wid+1),nunits*(wid+1));
    for A=1:nunits
        [~,idx]=sort(W(:,A),'descend');
        for B=1:nsub
            x=(B-1)*(wid+1)+1; y=(A-1)*(wid+1)+1;
            M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(B))/max(bank(:,idx(B))),wid,wid);
%             M(x:x+wid-1,y:y+wid-1)=reshape(bank(:,idx(B)),wid,wid);
        end;
    end;

    imagesc(M'); colormap(gray); axis image;
    colorbar;

    subplot(132);
    imagesc(W,[0 3*std(W(:))]);
    colorbar;
    
    subplot(133);
    bar(ks);

end
                  
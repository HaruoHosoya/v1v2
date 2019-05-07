function lgn_cells_example(dataset,nreduced,nshow,varargin)

pr=inputParser;
pr.addParamValue('pooling','dimred',@isstr);
pr.parse(varargin{:});
options=pr.Results;

% dataset=bsxfun(@minus,dataset,mean(dataset,1));
[datadim,datalen]=size(dataset); datawid=floor(sqrt(datadim));

switch options.pooling
    case 'whitening'
        [E,~,D]=pca(dataset');
        V=E(:,1:nreduced)*diag(D(1:nreduced).^(-1/2))*E(:,1:nreduced)';
%         V=V(:,1:nreduced)';
    case 'dimred'
        [E,~,D]=pca(dataset');
        V=E(:,1:nreduced)*E(:,1:nreduced)';
%         V=V(:,1:nreduced)';
end;

figure('position',[0 0 600 600]);
subplot(2,2,1);
images=reshape(E(:,1:nshow),datawid,datawid,1,nshow);
images(datawid+1,1:datawid+1,:,:)=-1;
images(1:datawid+1,datawid+1,:,:)=-1;
montage(images,[-0.3 0.3]);
% colorbar;

subplot(2,2,2);
images=reshape(V(:,1:nshow),datawid,datawid,1,nshow);
images(datawid+1,1:datawid+1,1,:)=-1;
images(1:datawid+1,datawid+1,1,:)=-1;
montage(images,[-0.3 0.3]);
% colorbar;

subplot(2,2,3);
imagesc(reshape(V(:,floor(datadim/2)+floor(datawid/2)),datawid,datawid));
axis image;
axis off;

end


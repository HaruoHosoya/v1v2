function energy_correlation_example(dataset)

[datadim,datalen]=size(dataset);
wid=ceil(sqrt(datadim));

[ix iy]=meshgrid(1:wid,1:wid);

g1=gabor(ix,iy,[wid/2 wid/2 1 wid/12 wid/6 pi/8 8/wid 0]);
g2=gabor(ix,iy,[wid/2 wid/2 1 wid/12 wid/6 pi/8 8/wid pi/2]);
g3=gabor(ix,iy,[wid/2 wid/2 1 wid/12 wid/6 pi/8*5 8/wid pi/2]);

y1=dataset'*g1(:);
y2=dataset'*g2(:);
y3=dataset'*g3(:);

y1=zscore(y1);
y2=zscore(y2);
y3=zscore(y3);

figure;
subplot(1,3,1);
imagesc(g1);
axis image; axis off;

subplot(1,3,2);
imagesc(g2);
axis image; axis off;

subplot(1,3,3);
imagesc(g3);
axis image; axis off;
colormap(gray);

m=10;
figure('Position',[0 0 600 400]);

% subplot(2,2,1);
% scatter(y2,y1,1);
% xlim([-m m]);ylim([-m m]);
% hold on;
% plot([-m m],[0 0],'k'); plot([0 0],[-m m],'k');
% axis square;

subplot(2,2,1);
cnt=hist3([y2 y1],{-m:0.5:m -m:0.5:m});
imagesc(-m:0.5:m,-m:0.5:m,log10(cnt/datalen),[-5 0]);
axis square; axis xy;
colorbar('location','EastOutside');
set(gca,'FontName','Times','FontSize',12);

fprintf('correlation coefficient of y1 vs y2: %1.3f\n',corr(y1.^2,y2.^2));

subplot(2,2,2);
cnt=hist3([y2 y1],{-m/2:0.25:m/2 -m/2:0.25:m/2});
cond=bsxfun(@rdivide,cnt,sum(cnt,1)); cond(isnan(cond(:)))=0;
imagesc(-m/2:0.25:m/2,-m/2:0.25:m/2,(cond),[0 0.2]);
axis square; axis xy;
colorbar('location','EastOutside');
set(gca,'FontName','Times','FontSize',12);


% subplot(2,2,3);
% scatter(y2,y3,1);
% xlim([-m m]);ylim([-m m]);
% hold on;
% plot([-m m],[0 0],'k'); plot([0 0],[-m m],'k');
% axis square;

subplot(2,2,3);
cnt=hist3([y2 y3],{-m:0.5:m -m:0.5:m});
h=imagesc(-m:0.5:m,-m:0.5:m,log10(cnt/datalen),[-5 0]);
axis square; axis xy;
colorbar('location','EastOutside');
set(gca,'FontName','Times','FontSize',12)

subplot(2,2,4);
cnt=hist3([y2 y3],{-m/2:0.25:m/2 -m/2:0.25:m/2});
cond=bsxfun(@rdivide,cnt,sum(cnt,1)); cond(isnan(cond(:)))=0;
imagesc(-m/2:0.25:m/2,-m/2:0.25:m/2,(cond),[0 0.5]);
axis square; axis xy;
colorbar('location','EastOutside');
set(gca,'FontName','Times','FontSize',12);

colormap(jet);

fprintf('correlation coefficient of y3 vs y2: %1.3f\n',corr(y3.^2,y2.^2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset32=h5read('~/datasets/imagenet/imagenet10k/processed2/patch32_100k/patches000001.h5','/data');
dataset16=dataset32(9:24,9:24,:);
dataset32=reshape(dataset32,32*32,size(dataset32,3));
dataset32=bsxfun(@minus,dataset32,mean(dataset32,1));
dataset16=reshape(dataset16,16*16,size(dataset16,3));
dataset16=bsxfun(@minus,dataset16,mean(dataset16,1));

%%

energy_correlation_example(dataset32);

%% LGN-like model

lgn_cells_example(dataset16,30,100);

%%

lgn_cells_example(dataset16,60,100);

%%

lgn_cells_example(dataset16,120,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              V1 model                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load model

load('data/v1_simple_cell_model_filter.mat','netf0');
load('data/v1_simple_cell_model_basis.mat','netb0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ica filters -> half rect -> dimension reduction 1/8

[netf1,nete1]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','half-rect');

%% display V1 simple filters

figure;
v1s_bases(netf0,'nodes',1,'normalize',true,'units',1:49,'disp_width',7);

%% display eigenvectors of V1 simple outputs

figure('position',[0 0 600 600]);
v1c_bases_ori(nete1,'nodes',1,'units',1:49,'disp_width',7,'barthickness',0.4,'axisticks',true,'min_weight',0.1);

%% display V1 complex filters

figure('position',[0 0 600 600]);
v1c_bases_ori(netf1,'nodes',1,'units',1:49,'disp_width',7,'barthickness',0.4,'axisticks',true,'min_weight',0.1);

%% display V1 complex filters (v1 simple combination)

figure('position',[0 0 220 500]);
v1c_bases_ori(netf1,'nodes',1,'units',1:7,'disp_width',1,'barthickness',0.3,'axisticks',false,'min_weight',0.1);
figure('position',[0 0 600 500]);
v1c_bases(netf1,'nodes',1,'units',1:7,'min_rel_weight',0,'disp_width',7);

%% ica filters -> half rect -> dimension reduction 1/4

[netf2,nete2]=v1_complex_cells(netf0,dataset16,48,'pooling','dimred','nonlin','half-rect');

%% display V1 complex filters

figure('position',[0 0 600 250]);
v1c_bases_ori(netf2,'nodes',1,'units',1:21,'disp_width',7,'barthickness',0.4,'axisticks',true,'min_weight',0.1);

%% ica filters -> half rect -> dimension reduction 1/4

[netf3,nete3]=v1_complex_cells(netf0,dataset16,12,'pooling','dimred','nonlin','half-rect');

%% display V1 complex filters

figure('position',[0 0 600 250]);
v1c_bases_ori(netf3,'nodes',1,'units',1:21,'disp_width',7,'barthickness',0.4,'axisticks',true,'min_weight',0.1);

%% ica filters -> squaring -> dimension reduction 1/4

[netf4,nete4]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','abs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure('position',[0 0 400 150]);
v1c_phase_tuning(netf1,'layer',2,'units',1:8,'disp_width',4);

figure('position',[0 0 400 150]);
v1c_phase_tuning(netf1,'layer',3,'units',1:8,'disp_width',4);

figure('position',[0 0 400 150]);
v1c_phase_tuning(netf4,'layer',3,'units',1:8,'disp_width',4);

%%

figure('position',[0 0 600 400]);
v1_f1f0ratios(netf1);

figure('position',[0 0 600 400]);
v1_f1f0ratios(netf4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ica filters -> half rect -> whitening dimension reduction 1/8 

[netf5,nete5]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred+whitening','nonlin','half-rect');

%%

figure('position',[0 0 600 170]);
v1_f1f0ratios(netf5,'layers',3);

%% ica filters -> half rect -> whitening 

[netf6,nete6]=v1_complex_cells(netf0,dataset16,24,'pooling','whitening','nonlin','half-rect');

%% ica filters -> half rect -> covariance

[netf7,nete7]=v1_complex_cells(netf0,dataset16,24,'pooling','cov','nonlin','half-rect');

%%

figure('position',[0 0 600 170]);
v1_f1f0ratios(netf6,'layers',3);
figure('position',[0 0 600 170]);
v1_f1f0ratios(netf7,'layers',3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ica filters -> half rect -> denoising 

[netf8,nete8]=v1_complex_cells(netf0,dataset16,24,'pooling','denoise','nonlin','half-rect');

%%

figure('position',[0 0 600 170]);
v1_f1f0ratios(netf8,'layers',3);

%%

figure('position',[0 0 600 400]);
v1_f1f0ratios_exp();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure('position',[0 0 700 300]); hold on;
d1=nete1.content.layers{3}.diagonal_values./max(nete1.content.layers{3}.diagonal_values);
d2=nete5.content.layers{3}.diagonal_values./max(nete5.content.layers{3}.diagonal_values);
d3=nete6.content.layers{3}.diagonal_values./max(nete6.content.layers{3}.diagonal_values);
d4=nete7.content.layers{3}.diagonal_values./max(nete7.content.layers{3}.diagonal_values);
d5=nete8.content.layers{3}.diagonal_values./max(nete8.content.layers{3}.diagonal_values);

plot([1.5:192.5; 1:192; 1:192; 1:192; 1:192]',[d1+0.004 d2+0.008 d3 d4 d5+0.012]);
xlim([1 193]);
ylim([-0.05 1.05]);
legend({'PCA dim. red.','whitening + dim. red.','whitening','covariance','optimal denoiser'},'location','east');
set(gca,'FontName','Times','FontSize',18)
box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              V2 model                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load model

load('data/v2_simple_cell_model_filter.mat','net2f0');
load('data/v2_simple_cell_model_basis.mat','net2b0');

%% v2 simple basis

figure('position',[0 0 550 600]);
v2s_bases(net2b0,'nodes',1,'units',1:20,'disp_width',4,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false,'linewidth',0.2);

%% v2 simple filters

figure('position',[0 0 550 600]);
v2s_bases(net2f0,'nodes',1,'units',1:20,'disp_width',4,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false,'linewidth',0.2);

%%

net2f0_tmp=run_net(net2f0,dataset32);
r=zscore(net2f0_tmp.content.layers{4}.unitProperties.resp);

%%

figure;
joint_dist_table(r,1:10,1:10,'style','histogram');
figure;
joint_dist_table(r,1:10,11:20,'style','histogram');
figure;
joint_dist_table(r,11:20,11:20,'style','histogram');

%%

m=10;
figure('position',[0 0 600 400]);
subplot(2,2,1);
c=hist3([r(:,13),r(:,3)],{-m:0.5:m,-m:0.5:m});
imagesc(-m:0.5:m,-m:0.5:m,log10(c(2:end-1,2:end-1)/size(r,1)),[-5 0]);
axis square;
colorbar('location','EastOutside');
axis xy;
set(gca,'FontName','Times','FontSize',12);

% m2=10;
% subplot(2,2,2);
% c=hist3([r(:,13),r(:,3)],{-m2:0.25:m2,-m2:0.25:m2});
% c=bsxfun(@rdivide,c,sum(c,1)); 
% imagesc(-m2:0.25:m2,-m2:0.25:m2,c,[0 0.2]);
% axis square; axis xy;
% colorbar('location','EastOutside');
% set(gca,'FontName','Times','FontSize',12);


subplot(2,2,3);
c=hist3([r(:,7),r(:,11)],{-m:0.5:m,-m:0.5:m});
imagesc(-m:0.5:m,-m:0.5:m,log10(c(2:end-1,2:end-1)/size(r,1)),[-5 0]);
axis square;
colorbar('location','EastOutside');
axis xy;
set(gca,'FontName','Times','FontSize',12);

% subplot(2,2,4);
% c=hist3([r(:,7),r(:,11)],{-m2:0.25:m2,-m2:0.25:m2});
% c=bsxfun(@rdivide,c,sum(c,1)); 
% imagesc(-m2:0.25:m2,-m2:0.25:m2,c,[0 0.2]);
% axis square; axis xy;
% colorbar('location','EastOutside');
% set(gca,'FontName','Times','FontSize',12);


colormap(jet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

load('data/v2_complex_cell_model_filter_halfrect.mat','net2f_halfrect');
load('data/v2_complex_cell_model_eigen_halfrect.mat','net2e_halfrect');

load('data/v2_complex_cell_model_filter_abs.mat','net2f_abs');
load('data/v2_complex_cell_model_eigen_abs.mat','net2e_abs');

%%

v2c_bases(net2f_halfrect,'nodes',1,'units',1:5,'nshow',6);

%% 

figure('position',[0 0 600 400]);
v1_f1f0ratios(net2f_halfrect,'layers',4:5);

figure('position',[0 0 400 150]);
v1c_orientation_tuning(net2f_halfrect,'layer',4,'units',1:8,'disp_width',4);

figure('position',[0 0 400 150]);
v1c_orientation_tuning(net2f_halfrect,'layer',5,'units',1:8,'disp_width',4);

%%

figure('Position',[0 0 400 250]);
v2_anzai_subrf(net2f_halfrect,'layer',4,'units',1:6,'disp_width',3,'linewidth',0.1);
figure('Position',[0 0 400 250]);
v2_anzai_subrf(net2f_halfrect,'layer',5,'units',1:6,'disp_width',3,'linewidth',0.1);

%%

figure('Position',[0 0 400 400]);
v2_anzai_stat(net2f_halfrect,'minresp',0.5,'layer',4);
figure('Position',[0 0 400 400]);
v2_anzai_stat(net2f_halfrect,'minresp',0.5,'layer',5);

%%

v2c_bases(net2f_abs,'nodes',1,'units',1:5,'nshow',6);

%% 

figure('position',[0 0 600 400]);
v1_f1f0ratios(net2f_abs,'layers',4:5);

figure('position',[0 0 400 150]);
v1c_orientation_tuning(net2f_abs,'layer',5,'units',1:8,'disp_width',4);

%%

figure('Position',[0 0 400 250]);
v2_anzai_subrf(net2f_abs,'layer',5,'units',1:6,'disp_width',3,'linewidth',0.1);

%%

figure('Position',[0 0 400 400]);
v2_anzai_stat(net2f_abs,'minresp',0.5,'layer',5);

%%

ii1=net2f_halfrect.content.layers{4}.unitProperties.invariance_index;
ii2=net2f_halfrect.content.layers{5}.unitProperties.invariance_index;
[h p]=ttest(ii1,ii2);
fprintf('mean increase of position invariance: %f (p=%1.3e) \n',mean(ii2)/mean(ii1),p);


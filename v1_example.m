%%%%%%%%%%%%%%%%%%% process input data %%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data

dataset32=h5read('~/datasets/imagenet/imagenet10k/processed2/patch32_100k/patches000001.h5','/data');
dataset16=dataset32(9:24,9:24,:);
dataset32=reshape(dataset32,32*32,size(dataset32,3));
dataset32=bsxfun(@minus,dataset32,mean(dataset32,1));
dataset16=reshape(dataset16,16*16,size(dataset16,3));
dataset16=bsxfun(@minus,dataset16,mean(dataset16,1));

%%%%%%%%%%%%%%%%%% construct a simple cel model %%%%%%%%%%%%%%%%%
%% construct a simple cell model

[netf0,netb0]=v1_simple_cells(dataset16,192);

%% save model

save('data/v1_simple_cell_model_filter.mat','netf0');
save('data/v1_simple_cell_model_basis.mat','netb0');

%% load model

load('data/v1_simple_cell_model_filter.mat','netf0');
load('data/v1_simple_cell_model_basis.mat','netb0');

%%%%%%%%%%%%%%%%%% construct a complex cell model %%%%%%%%%%%%%%%
%% ica filters -> half rect -> dimension reduction 1/8 

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','half-rect');

%% ica filters -> half rect -> dimension reduction 1/4 

[netf,nete]=v1_complex_cells(netf0,dataset16,48,'pooling','dimred','nonlin','half-rect');

%% ica filters -> half rect -> dimension reduction 1/16

[netf,nete]=v1_complex_cells(netf0,dataset16,12,'pooling','dimred','nonlin','half-rect');

%% ica filters -> square -> dimension reduction 1/8

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','square');

%% ica filters -> abs -> dimension reduction 1/8

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','abs');

%% ica filters -> half square -> dimension reduction 1/8

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred','nonlin','half-square');

%% ica filters -> half rect -> whitening dimension reduction 1/8 

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','dimred+whitening','nonlin','half-rect');

%% ica filters -> half rect -> whitening 

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','whitening','nonlin','half-rect');

%% ica filters -> half rect -> covariance

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','cov','nonlin','half-rect');

%% ica filters -> half rect -> covariance (minus smallest eigenvalue)

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','cov2','nonlin','half-rect');

%% ica filters -> half rect -> inverse covariance

[netf,nete]=v1_complex_cells(netf0,dataset16,144,'pooling','invcov','nonlin','half-rect');

%% ica filters -> square -> inverse covariance

[netf,nete]=v1_complex_cells(netf0,dataset16,144,'pooling','invcov','nonlin','square');

%% ica filters -> half rect -> denoising

[netf,nete]=v1_complex_cells(netf0,dataset16,24,'pooling','denoise','nonlin','half-rect');

%% ica filters -> square -> denoising

[netf,nete]=v1_complex_cells(netf0,dataset16,144,'pooling','denoise','nonlin','square');

%% ica filters -> half rect -> covariance + dimension reduction 1/8

[netf,nete]=v1_complex_cells(netf0,dataset16,144,'pooling','dimred+cov','nonlin','half-rect');

%%%%%%%%%%%%%%%%%%% display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display V1 simple filters

figure;
v1s_bases(netf0,'nodes',1,'normalize',true,'units',1:192,'disp_width',16);

%% display V1 simple basis

figure;
v1s_bases(netb0,'nodes',1,'normalize',true,'units',1:192,'disp_width',16);

%% display V1 complex filters (v1 simple combination)

figure;
v1c_bases_ori(netf,'nodes',1,'units',1:8,'disp_width',1,'barthickness',0.5,'axisticks',false,'min_weight',0.0);

v1c_bases(netf,'nodes',1,'units',1:8,'min_rel_weight',0,'disp_width',8);

%% display V1 complex filters

figure;
v1c_bases_ori(netf,'nodes',1,'units',1:64,'barthickness',0.5,'axisticks',false);
    
%% display eigenvectors of V1 simple outputs

figure;
v1c_bases_ori(nete,'nodes',1,'units',1:64,'barthickness',0.5,'axisticks',false);

%% display analysis summary

figure;
v1_f1f0ratios(netf);

%% display phase tuning functions

figure;
v1c_phase_tuning(netf,'layer',2,'units',1:8,'disp_width',4);

figure;
v1c_phase_tuning(netf,'layer',3,'units',1:8,'disp_width',4);

%% display phase orientation functions

figure;
v1c_orientation_tuning(netf,'layer',2,'units',1:8,'disp_width',4);

figure;
v1c_orientation_tuning(netf,'layer',3,'units',1:8,'disp_width',4);

%%

plot2pdf(gcf,'~/Desktop/a.pdf','size',[40 60]);

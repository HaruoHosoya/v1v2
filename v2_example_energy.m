%%%%%%%%%%%%%%%%%%% process input data %%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data

IMAGES=load('~/workspace/bv/matlab/sparsenet/IMAGES_RAW.mat','IMAGESr');
dataset32=create_image_dataset(IMAGES.IMAGESr,32,100000);

%% save data

save('~/datasets/olshausen/data_pooling/dataset32.mat','dataset32','-v7.3');

%% load data

load('~/datasets/olshausen/data_pooling/dataset32.mat','dataset32');

%%%%%%%%%%%%%%%%%% construct a simple cel model %%%%%%%%%%%%%%%%%
%% fastica

[net2fa,net2ba]=v2_simple_cells(dataset32,100);

%% score-matching / energy model

[net2f0,net2b0]=v2_simple_cells(dataset32,432,'learning','overcomplete-ica','reduction',12,'nonlin','half-rect','energy_model',true,'norient',12);

%% save model

save('data/v2en_simple_cell_model_filter.mat','net2f0','-v7.3');
save('data/v2en_simple_cell_model_basis.mat','net2b0','-v7.3');

%% load model

load('data/v2en_simple_cell_model_filter.mat','net2f0');
load('data/v2en_simple_cell_model_basis.mat','net2b0');

%%%%%%%%%%%%%%%%%% construct a complex cell model %%%%%%%%%%%%%%%
%% v2 simple -> half rect -> dimension reduction 1/8 

[net2f_halfrect,net2e_halfrect]=v2_complex_cells(net2f0,dataset32,54,'pooling','dimred','nonlin','half-rect');

%% v2 simple -> square -> dimension reduction 1/8 

[net2f_square,net2e_square]=v2_complex_cells(net2f0,dataset32,54,'pooling','dimred','nonlin','square');

%% v2 simple -> half square -> dimension reduction 1/8 

[net2f_halfsq,net2e_halfsq]=v2_complex_cells(net2f0,dataset32,54,'pooling','dimred','nonlin','half-square');

%% save model

save('data/v2en_complex_cell_model_filter_halfrect.mat','net2f_halfrect');
save('data/v2en_complex_cell_model_eigen_halfrect.mat','net2e_halfrect');

%% load model

load('data/v2en_complex_cell_model_filter_halfrect.mat','net2f_halfrect');
load('data/v2en_complex_cell_model_eigen_halfrect.mat','net2e_halfrect');

%% save model

save('data/v2en_complex_cell_model_filter_square.mat','net2f_square');
save('data/v2en_complex_cell_model_eigen_square.mat','net2e_square');

%% load model

load('data/v2en_complex_cell_model_filter_square.mat','net2f_square');
load('data/v2en_complex_cell_model_eigen_square.mat','net2e_square');

%% save model

save('data/v2en_complex_cell_model_filter_halfsq.mat','net2f_halfsq');
save('data/v2en_complex_cell_model_eigen_halfsq.mat','net2e_halfsq');

%% load model

load('data/v2en_complex_cell_model_filter_halfsq.mat','net2f_halfsq');
load('data/v2en_complex_cell_model_eigen_halfsq.mat','net2e_halfsq');

%%%%%%%%%%%%%%%%%% display a model and a property %%%%%%%%%%%%%%%
%% compute response

net2f0_tmp=run_net(net2f0,dataset32);
r=net2f0_tmp.content.layers{4}.unitProperties.resp;

%% v2 simple basis

figure;
v2s_bases(net2b0,'nodes',1,'units',1:100,'disp_width',10,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false);

%% v2 simple basis

figure;
v2s_bases(net2b0,'nodes',1,'units',1:20,'disp_width',5,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false);

%% v2 simple filters

figure;
v2s_bases(net2f0,'nodes',1,'units',1:100,'disp_width',10,'onlyextreme',true,'axisticks',false,'axislabels',false);

%% v2 simple filters

figure;
v2s_bases(net2f0,'nodes',1,'units',1:20,'disp_width',5,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false);

%%

plot2pdf(gcf,'~/Desktop/a.pdf','size',[60 40]);

%%

figure;
v2s_bases(net2f0,'nodes',1,'units',[57 92 17 56],'disp_width',2,'onlyextreme',true,'panel_titles','none','axisticks',false,'axislabels',false);

figure;
joint_dist_table(r,57,92,'style','histogram'); 

figure;
joint_dist_table(r,17,56,'style','histogram'); 


%%

figure;
joint_dist_table(r,1:8,1:8,'style','histogram');

%%

figure;
joint_dist_table(r,1:8,1:8,'style','contour');


%% display analysis summary

% figure;
v1_f1f0ratios(net2f_halfrect,'layers',4:5);
figure;
v1c_orientation_tuning(net2f_halfrect,'layer',4,'units',1:8,'disp_width',4);
figure;
v1c_orientation_tuning(net2f_halfrect,'layer',5,'units',1:8,'disp_width',4);

%%

% figure;
v1_f1f0ratios(net2f_square,'layers',4:5);
figure;
v1c_orientation_tuning(net2f_square,'layer',5,'units',1:8,'disp_width',4);

%%

% figure;
v1_f1f0ratios(net2f_halfsq,'layers',4:5);
figure;
v1c_orientation_tuning(net2f_halfsq,'layer',4,'units',1:8,'disp_width',4);
figure;
v1c_orientation_tuning(net2f_halfsq,'layer',5,'units',1:8,'disp_width',4);

%%

v2c_bases(net2f_halfrect,'nodes',1,'units',1:5,'nshow',4);

%%

v2c_bases(net2f_square,'nodes',1,'units',1:5,'nshow',4);

%%

v2c_bases(net2f_halfsq,'nodes',1,'units',1:5,'nshow',4);

%%

v2_anzai_subrf(net2f_halfrect,'layer',4,'units',1:9);
v2_anzai_subrf(net2f_halfrect,'layer',5,'units',1:9);
v2_anzai_stat(net2f_halfrect,'minresp',0.5,'layer',4);
v2_anzai_stat(net2f_halfrect,'minresp',0.5,'layer',5);

%%

v2_anzai_subrf(net2f_square,'layer',4,'units',1:9);
v2_anzai_subrf(net2f_square,'layer',5,'units',1:9);
v2_anzai_stat(net2f_square,'minresp',0.5,'layer',4);
v2_anzai_stat(net2f_square,'minresp',0.5,'layer',5);

%%

v2_anzai_subrf(net2f_halfsq,'layer',4,'units',1:9);
v2_anzai_subrf(net2f_halfsq,'layer',5,'units',1:9);
v2_anzai_stat(net2f_halfsq,'minresp',0.5,'layer',4)
v2_anzai_stat(net2f_halfsq,'minresp',0.5,'layer',5);

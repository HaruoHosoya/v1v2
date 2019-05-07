addpath(pwd);
addpath([pwd filesep 'models']);
addpath([pwd filesep 'preprocess']);
addpath([pwd filesep 'functions']);
addpath([pwd filesep 'v1']);
addpath([pwd filesep 'v2']);
addpath([pwd filesep 'a1']);
addpath([pwd filesep 'denoise']);
wd=pwd; 
cd('../FastICA_25'); initpath; 
% cd('../analysis2'); initpath;
cd('../scoreMatching'); startup;
cd(wd);



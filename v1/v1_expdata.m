function data=v1_expdata

% De Valois RL, Albrecht DG, Thorell LG: Spatial frequency selectivity of cells in macaque visual cortex. 
% Vision Research 1982, 22:545?559.

data.species='macaque';

data.freq_bw.x=[0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 Inf];

data.freq_bw.complex_parafoveal.y=[0 1 4 6 5 6 4 7 6 4 2 3 4];
data.freq_bw.complex_parafoveal.median=1.60;
data.freq_bw.complex_parafoveal.N=52;

data.freq_bw.complex_foveal.y=[0 11 5 10 13 7 7 6 8 3 3 8];
data.freq_bw.complex_foveal.median=1.45;
data.freq_bw.complex_foveal.N=81;

data.freq_bw.simple_parafoveal.y=[1 7 9 15 10 11 9 4 3 2 1 6];
data.freq_bw.simple_parafoveal.median=1.32;
data.freq_bw.simple_parafoveal.N=78;

data.freq_bw.simple_foveal.y=[2 9 12 20 24 21 12 16 8 8 8 7];
data.freq_bw.simple_foveal.median=1.45;
data.freq_bw.simple_foveal.N=147;

data.peak_freq.x=[-Inf 0.5 0.7 1.0 1.4 2.0 2.8 4.0 5.6 8.0 11.2 16.0 Inf];

data.peak_freq.complex_parafoveal.y=[1 2 1 2 5 15 13 9 3 2];
data.peak_freq.complex_parafoveal.mean=3.2;
data.peak_freq.complex_parafoveal.N=53;

data.peak_freq.complex_foveal.y=[0 0 0 9 9 7 10 23 12 8 3 3];
data.peak_freq.complex_foveal.mean=5.1;
data.peak_freq.complex_foveal.N=84;

data.peak_freq.simple_parafoveal.y=[2 4 10 12 18 7 18 3 4 0 0 0];
data.peak_freq.simple_parafoveal.mean=2.2;
data.peak_freq.simple_parafoveal.N=78;

data.peak_freq.simple_foveal.y=[0 4 4 8 25 33 26 28 12 5 2 1];
data.peak_freq.simple_foveal.mean=3.5;
data.peak_freq.simple_foveal.N=148;

data.freq_peak_bw.x=[-Inf 0.7 1.0 2.0 4.0 8.0 11.0 Inf];

data.freq_peak_bw.y=[2.15 1.8 1.62 1.48 1.3 1.05];
data.freq_peak_bw.y_sigma=[1.09 0.75 0.61 0.64 0.48 0.55];
data.freq_peak_bw.N=[13 45 119 127 46 8];

% De Valois RL, Yund EW, Hepler N: The orientation and direction selectivity of cells in macaque visual cortex. 
% Vision Research 1982, 22:531?544.

data.ori_bw.x=[0 10 20 30 40 50 60 70 80 90 180 Inf];

data.ori_bw.simple_foveal.y=[2 9 12 16 13 8 8 7 3 10 4];
data.ori_bw.simple_foveal.median=42.3;
data.ori_bw.simple_foveal.N=92;

data.ori_bw.simple_parafoveal.y=[0 3 12 13 5 4 4 1 3 1 1];
data.ori_bw.simple_parafoveal.median=47;
data.ori_bw.simple_parafoveal.N=33.8;

data.ori_bw.complex_foveal.y=[0 0 11 7 8 4 2 3 0 5 6];
data.ori_bw.complex_foveal.median=44.5;
data.ori_bw.complex_foveal.N=46;

data.ori_bw.complex_parafoveal.y=[0 2 7 9 3 3 4 0 2 5 2];
data.ori_bw.complex_parafoveal.median=37;
data.ori_bw.complex_parafoveal.N=43.8;

data.ori_tuning.x=[0 45 90 135];
data.ori_tuning.parafoveal=[34 34 33 38];
data.ori_tuning.foveal=[66 49 77 54 66];

% 1. Skottun BC, De Valois RL, Grosof DH, Movshon JA, Albrecht DG, Bonds AB: 
% Classifying simple and complex cells on the basis of response modulation. 
% Vision Research 1991, 31:1078?1086.

data.f1f0.x=[0:0.1:2.1 2.2];
data.f1f0.y=[9 17 39 24 30 29 21 19 17 8 9 23 28 34 31 39 42 35 10 4 6 1 38];
data.f1f0.N=513;

end

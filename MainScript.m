rng(1)
RPath='C:\Program Files\R\R-3.4.2\bin\x64';
RPath2=strcat('"',RPath,'"');
callname= strcat({RPath2},{'\'},{'R CMD BATCH '},{'"'},pwd,'\VariableSelection.R',{'"'});
Rcallname=callname{1};
callname2= strcat({RPath2},{'\'},{'R CMD BATCH '},{'"'},pwd,'\Runwfmm.R',{'"'});
Rcallname2=callname2{1};
addpath(genpath(strcat(pwd,'\Matlab_functions')));
addpath(genpath(strcat(pwd,'\CodeDWT')));


% Part (a) - Wavelet basis and data compression

load('Y_outlier_removed')
%load('Y_simulated.mat')
wavespecs.compress=0.995;
wavespecs.ndim=2;
wavespecs.t=[120, 120];
wavespecs.rectangular=1;      %%% rectangular 2d transform
wavespecs.wavelet='db3';
wavespecs.nlevels=[5,5];
wavespecs.boundary=['per';'sym'];

tic
[D,wavespecs_compressed,D_compressed]=DWT_rows(Y_OutRemoved,wavespecs);
%[D,wavespecs_compressed,D_compressed]=DWT_rows(Y_simulated,wavespecs);
toc

save('DCompressed','wavespecs_compressed','D_compressed')

% Part (b) - Call R-script to run model selection
% It takes about 24 minutes to complete
system(Rcallname)

% Part (c) - Call R-script to run the Bayesian model
% It takes 3-4 hours to complete, you can follow the progress in the log
% file eye_filtered_wave_fmm\eye_filtered_wave_fmm_log.log
system(Rcallname2)



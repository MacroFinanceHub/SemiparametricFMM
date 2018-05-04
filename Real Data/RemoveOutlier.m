mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
addpath(genpath(strcat(newdir,'\CodeDWT')))
addpath(genpath(strcat(newdir,'\Matlab_Functions')))%%%%% Load data
load('RawEyeData.mat');
wavespecs.compress=1;
wavespecs.ndim=2;
wavespecs.t=[120, 120];
wavespecs.rectangular=1;      %%% rectangular 2d transform
wavespecs.wavelet='db3';
wavespecs.nlevels=[5,5];
wavespecs.boundary=['per';'sym'];

%Exclude 14L
index_14L=find(sample_info(:,5)==14 & sample_info(:,6)==1);

Y_2=Y;
Y_2(index_14L,:)=[];

tic
[D,wavespecs]=DWT_rows(Y_2,wavespecs);
toc

outlier_flag=find(mean(abs(D))>100*median(abs(D)));
temp_keep=ones(1,size(D,2));
temp_keep(outlier_flag)=0;
keep_logical=(temp_keep==1);
wavespecs2=wavespecs;

wavespecs2.keep=keep_logical;
groupKj=zeros(size(wavespecs2.Kj,2),1);
k=1;
for j=1:size(wavespecs2.Kj,2)
groupKj(j,1)=sum(wavespecs2.keep(k:(k+wavespecs2.Kj(j)-1)));
k=k+wavespecs2.Kj(j);
end;
wavespecs2.K=sum(groupKj);
wavespecs2.Kj=groupKj;
Y1=IDWT_rows(D,wavespecs2); % back to data space

sample_info2=sample_info;
sample_info2(index_14L,:)=[];
%save('OutlierRemovedData','Y1','sample_info2');
%dlmwrite('sample_info2.txt',sample_info2,'delimiter','\t');

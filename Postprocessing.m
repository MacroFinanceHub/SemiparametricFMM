%%%%%%%%%%%%%%%%%%%%%%%  This file is for postprocessing
rng(1)
maindir=pwd;
addpath(genpath(strcat(maindir,'\CodeDWT')))
addpath(genpath(strcat(maindir,'\Matlab_functions')))

load(strcat(maindir,'\DCompressed.mat'));
load(strcat(maindir,'\eye_filtered_wave_fmm\eye_filtered_wave_fmm.mat'),'model');

K = sum(wavespecs_compressed.Kj);
[N, P] = size(model.X);
H = size(model.m,2);
B = 1000;

filename1 = strcat(maindir,'\eye_filtered_wave_fmm\eye_filtered_wave_fmm_Results_0_wbeta.dat');
MCMC_beta = readbeta(B,P,K,filename1);
save('MCMC_beta.mat','MCMC_beta');

filename2 = strcat(maindir,'\eye_filtered_wave_fmm\eye_filtered_wave_fmm_Results_0_theta.dat');
MCMC_theta = readbeta(B,H+1,K,filename2);

save('MCMC_theta.mat','MCMC_theta');

%%%%%%%%% Generating posterior samples for spline random effects


Z_temp=model.Z;
m=model.m;

Z=cell(length(m),1);
Z{1}=Z_temp(:,1:m(1));
cm=cumsum(m);
for i=2:length(m)
    Z{i}=Z_temp(:,cm(i-1)+1:cm(i));
end

N=size(Z{1},1);
B=size(MCMC_beta,1);


model.H=4;
model.p=size(model.X,2);
H_spline=1;
H_0=model.H-H_spline;


h_0=1:H_0;
h_spline=(H_0+1):model.H; %the indices of random effects levels that are splines, assumed to be last
m_0=sum(model.m(h_0));
m_spline=sum(model.m(h_spline));

MCMC_Ustar=zeros(B,wavespecs_compressed.K*m_spline);


for h=1:H_0
    Z0Z0p{h}=Z{h_0(h)}*Z{h_0(h)}';
end;
for h=1:H_spline;
    ZsZsp{h}=Z{h_spline(h)}*Z{h_spline(h)}';
end;
Zs=Z{h_spline(1)};
for h=2:H_spline;
    Zs=[Zs,Z{h_spline(h)}];
end;

% Generating the posterior samples Ustar
% This step takes about 2.5-3 hours

for b=1:B
    betahat=reshape((MCMC_beta(b,:)),wavespecs_compressed.K,model.p)';
    thetahat=reshape(MCMC_theta(b,:),wavespecs_compressed.K,(model.H+1))';
    
    M=zeros(m_spline,wavespecs_compressed.K);
    Uhat=M;
    U=M;
    sig=M;
    for k=1:wavespecs_compressed.K
      
        Sig0_k=thetahat(end,k)*eye(N);
        for h=h_0
            Sig0_k=Sig0_k+Z0Z0p{h}*thetahat(h,k);  
        end;
        Sig0inv_k=pinv(Sig0_k);
        SigSinv_k=Zs'*Sig0inv_k*Zs;
        Qsjk_inv=diag((uneqkron(model.m(h_spline))*thetahat(h_spline,k)).^(-1));
        V_k=pinv(Qsjk_inv+SigSinv_k);
        Uhat_k=pinv(SigSinv_k)*Zs'*Sig0inv_k*(D_compressed(:,k)-model.X*betahat(:,k));
        Shrink_k=V_k*SigSinv_k;
        M_k=Shrink_k*Uhat_k;
        sig(:,k)=diag(V_k);
        M(:,k)=M_k;
        Uhat(:,k)=Uhat_k;
        [evec,eval]=svd(V_k);
        Vrank=sum(diag(eval)>1e-15);
        sqrtV_k=evec(:,1:Vrank)*eval(1:Vrank,1:Vrank)^(1/2)*evec(:,1:Vrank)';
        U_k=sqrtV_k*normrnd(zeros(m_spline,1),1)+M_k;
        U(:,k)=U_k;
    end;

    MCMC_Ustar(b,:)=reshape(U',1,m_spline*wavespecs_compressed.K);
    if mod(b,10) == 0
        fprintf('Ustar at iteration %d...\n',b);
    end
    
end;


T=wavespecs_compressed.T;
K2=wavespecs_compressed.K;
K1=wavespecs_compressed.K;
P=model.p;


%%%%%%%%% Generating posterior samples for random effects in original space
%%%% This step takes around 9 minutes
MCMC_U=zeros(B,m_spline*T);

for i=1:m_spline
  
    idx_in=((i-1)*K1+1):(i*K1);
    temp_U=MCMC_Ustar(:,idx_in);
    idx_out=((i-1)*T+1):(i*T);
    MCMC_U(:,idx_out)=IDWT_rows(temp_U,wavespecs_compressed);
    %fprintf('U at iteration %d...\n',i)
end

save('MCMC_U_all.mat','MCMC_U','MCMC_Ustar','-v7.3')

clear MCMC_U MCMC_Ustar
%%%%%%%%% Generating posterior samples for fixed effects in original space
%%%% This step takes around 45-50 minutes
MCMC_g=zeros(B,P*T);

for p=1:P
    idx_in=((p-1)*K1+1):(p*K1);
    temp_g=MCMC_beta(:,idx_in);
    idx_out=((p-1)*T+1):(p*T);
    MCMC_g(:,idx_out)=IDWT_rows(temp_g,wavespecs_compressed);
end


save('MCMC_g_low_pass.mat', 'MCMC_g')

%%%%%%%% Compute nonparametric age effect

load('Zg_spline.txt');
load('MCMC_U_all.mat','MCMC_U','MCMC_Ustar');



X_pressure_unique=[7,10,15,20,25,30,35,40,45]';
Z_unique=Zg_spline;
Z_age_only=Z_unique([1,6,11,16,21,26,31,36,41,46,51,56,61,66,71],:);
X_age_unique=5*(4:18)';
X_age_only=[zeros(size(X_age_unique)),X_age_unique,zeros(length(X_age_unique),2)];

a_fhat_MCMC_part1=zeros(1000,length(X_age_unique)*14400);
MCMC_fstar=zeros(B,length(X_age_unique)*K);


%This step takes about 15-20 minutes
for b=1:B 
    a_fhat_linear=X_age_only*reshape(MCMC_g(b,:),T,model.p)';   %%% age linear effect
    a_fhat=a_fhat_linear+Z_age_only*reshape(MCMC_U(b,:),T,m_spline)';  %%% age spline_effect
    a_fhat_MCMC_part1(b,:)=reshape(a_fhat,1,length(X_age_unique)*14400);    
end;
save('nonparametric_age_MCMC.mat','a_fhat_MCMC_part1')


for b=1:B %
    a_fhat_linear=X_age_only*reshape(MCMC_beta(b,:),K,model.p)';   %%% age linear effect
    a_fhat=a_fhat_linear+Z_age_only*reshape(MCMC_Ustar(b,:),K,m_spline)';  %%% age spline_effect
    MCMC_fstar(b,:)=reshape(a_fhat,1,length(X_age_unique)*K);    
end;
save('MCMC_fstar.mat','MCMC_fstar')

%%%%%%%%% Generating variance components mean in original space


[Phi_1,Kj1] = Get_DWT('db3',120,'per       ',0,5);
[Phi_2,Kj2] = Get_DWT('db3',120,'reflection',1,5);
Phi_all = kron(Phi_1,Phi_2); % DWT matrix 
IDWT = kron((Phi_1'*Phi_1)\Phi_1',(Phi_2'*Phi_2)\Phi_2'); %Inverse DWT matrix

H_all=model.H+1;
T=wavespecs_compressed.T;
K2=wavespecs_compressed.K;
MCMC_sig_mean=zeros(1,H_all*T);


Tstar = sum(wavespecs_compressed.Kj_all);
T = 14400;
sig_star = zeros(model.H+1, Tstar); 
%%% Get reordered "keep" vector to match order in IDWT matrix
keep_reorder=wavespecs_compressed.keep;
keep_reorder(wavespecs_compressed.reorder)=keep_reorder;



    theta=reshape(mean(MCMC_theta),wavespecs_compressed.K,H_all)';

    for h=1:(H_all)
    sig_tmp = zeros(Tstar,1);
    sig_tmp(wavespecs_compressed.keep==1) = theta(h,:);
    sig_tmp(wavespecs_compressed.reorder) = sig_tmp;
    sig_star(h,:) = sig_tmp;
    tmp_sig = sig_star(h,keep_reorder==1)';
    sig=zeros(T,1);
        for t=1:T
        sig(t)=(IDWT(t,keep_reorder==1).^2)*tmp_sig;     
        end;
     idx_out=((h-1)*T+1):(h*T);
     MCMC_sig_mean(1,idx_out)=sig';
    end;


save('MCMC_sig_mean.mat','MCMC_sig_mean');

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute degrees of freedom - average
load('MCMC_sig_mean.mat');
load('eye_filtered_wave_fmm\eye_filtered_wave_fmm.mat','model')
B=load('B.txt');
Omega=load('Omega.txt');
%%%%% Computing DF
T=14400;
N=size(model.X,1);

M = length(model.m)-1; % the number of random effects excluding one for spline
Z=model.Z(:,1:sum(model.m(1:M)));
N_z = [0,cumsum(model.m)];

mvar_u = MCMC_sig_mean(1,1:T*M);
mvar_s = MCMC_sig_mean(1,(T*M)+1:T*(M+1));
mvar_e = MCMC_sig_mean(1,T*(M+1)+1:T*(M+2));
mvar_u_r=reshape(mvar_u',T,M);
mvar_s_r=mvar_s';
mvar_e_r=mvar_e';
df_age_mean=zeros(1,T);

 
for t=1:T
    
    random_var = zeros(1, size(Z,2)); 
 
    for h=1:M
       tmp_idx = (N_z(h)+1):N_z(h+1); 
       random_var(tmp_idx) = mvar_u_r(t,h);
    end
    
    Sig_U = diag(random_var);
    W_inv = eye(N)*mvar_e_r(t,1)+Z*Sig_U*Z';
    BWB = B'*(W_inv\B);
    df_age_mean(1,t)=trace((BWB+Omega/mvar_s_r(t))\BWB);
    
end

save 'df_mean.mat' df_age_mean;

clear

load('eye_filtered_wave_fmm\eye_filtered_wave_fmm.mat','model')
load('DCompressed.mat')

load('Zg_spline.txt');
load('Zdg_spline.txt');

load('MCMC_U_all.mat','MCMC_U');
load('MCMC_g_low_pass.mat');

%%%%%% Set up evaluation design matrices 
model.p=size(model.X,2);
X_age_unique=(20:90)';
X_pressure_unique=[7,10,15,20,25,30,35,40,45]';
n_eye = size(model.X,1)/9;
unique_idx = n_eye:n_eye:size(model.X,1);
X_IOP_unique= model.X(unique_idx,3:4);

X_age_all=kron(ones(length(X_pressure_unique),1),X_age_unique);
X_IOP_all=kron(X_IOP_unique,ones(length(X_age_unique),1));



X_expand=[ones(size(X_age_all)),X_age_all,X_IOP_all];
Xd_expand=[zeros(size(X_age_all)),ones(size(X_age_all)),zeros(size(X_IOP_all))];

Z_unique=Zg_spline;
Zd_unique=Zdg_spline; %%% derivatives

Z_all=kron(ones(length(X_IOP_unique),1),Z_unique);
Zd_all=kron(ones(length(X_IOP_unique),1),Zd_unique);

n_age=length(X_age_unique);
n_pressure=size(X_IOP_unique,1);
ngrid=n_age*n_pressure;

H_spline=4;
h_spline=H_spline;
m_spline=sum(model.m(H_spline));

B=size(MCMC_U,1);
T=size(MCMC_g,2)/model.p;

%%%% Note: Compute predicted age x pressure (IOP) surface of MPS, and also compute AUC
%%%% of MPS x IOP .
%%%%
%%%% given curve MPS=f(IOP) for IOP={7, 10, 15, 20, ..., 45}, with
%%%% f(5mm)=0, the AUC using the trapezoidal rule is given by:
%%%%
%%%%    5/2 f(7) + 4 f(10) + 5 f(15) + 5 f(20) + ... + 5 f(40) + 5/2 f(45)
AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
AUC_wgts_all=kron(AUC_wgts,eye(n_age));

fhat.mean=zeros(ngrid,T);
fhat.SS=zeros(ngrid,T);
fhat.sd=zeros(ngrid,T);
fhat.AUC_mean=zeros(n_age,T);
fhat.AUC_SS=zeros(n_age,T);
fhat.AUC_sd=zeros(n_age,T);
fhat_linear=fhat;

Y_mean=repmat(IDWT_rows(mean(D_compressed),wavespecs_compressed),ngrid,1);


for b=1:B % 587.933258 sec
    a_fhat_linear=X_expand*reshape(MCMC_g(b,:),T,model.p)';   %%% age linear effect
    a_fhat=a_fhat_linear+Z_all*reshape(MCMC_U(b,:),T,m_spline)';  %%% age spline_effect
    fhat_linear.mean=fhat_linear.mean+a_fhat_linear/B;
    fhat.mean=fhat.mean+a_fhat/B;
    fhat_linear.SS=fhat_linear.SS+a_fhat_linear.^2;
    fhat.SS=fhat.SS+a_fhat.^2;
    % Now AUC
    fhat_linear.AUC_mean=fhat_linear.AUC_mean+AUC_wgts_all*a_fhat_linear/B;
    fhat.AUC_mean=fhat.AUC_mean+AUC_wgts_all*a_fhat/B;
    fhat_linear.AUC_SS=fhat_linear.AUC_SS+(AUC_wgts_all*a_fhat_linear).^2;
    fhat.AUC_SS=fhat.AUC_SS+(AUC_wgts_all*a_fhat).^2;
end;

fhat_linear.sd=((B-1)^(-1)*(fhat_linear.SS-B*(fhat_linear.mean).^2)).^(1/2);
fhat.sd=((B-1)^(-1)*(fhat.SS-B*(fhat.mean).^2)).^(1/2);
fhat_linear.Q025=fhat_linear.mean-1.96*fhat_linear.sd;
fhat_linear.Q975=fhat_linear.mean+1.96*fhat_linear.sd;
fhat.Q025=fhat.mean-1.96*fhat.sd;
fhat.Q975=fhat.mean+1.96*fhat.sd;
%%% Now AUC
fhat_linear.AUC_sd=((B-1)^(-1)*(fhat_linear.AUC_SS-B*(fhat_linear.AUC_mean).^2)).^(1/2);
fhat.AUC_sd=((B-1)^(-1)*(fhat.AUC_SS-B*(fhat.AUC_mean).^2)).^(1/2);
fhat_linear.AUC_Q025=fhat_linear.AUC_mean-1.96*fhat_linear.AUC_sd;
fhat_linear.AUC_Q975=fhat_linear.AUC_mean+1.96*fhat_linear.AUC_sd;
fhat.AUC_Q025=fhat.AUC_mean-1.96*fhat.AUC_sd;
fhat.AUC_Q975=fhat.AUC_mean+1.96*fhat.AUC_sd;
toc 

%%% Now compute 95% simultaneous credible bands

z_fhat_linear=zeros(B,1);
z_fhat=zeros(B,1);
z_fhat_linear_AUC=zeros(B,1);
z_fhat_AUC=zeros(B,1);


for b=1:B %410.652154 sec
    a_fhat_linear=X_expand*reshape(MCMC_g(b,:),T,model.p)';   %%% age linear effect
    a_fhat=a_fhat_linear+Z_all*reshape(MCMC_U(b,:),T,m_spline)';  %%% age spline_effect
    a_fhat_linear_AUC=AUC_wgts_all*a_fhat_linear;
    a_fhat_AUC=AUC_wgts_all*a_fhat;
    
    z_fhat_linear(b,1)=max(max(abs(a_fhat_linear-fhat_linear.mean)./fhat_linear.sd));
    z_fhat(b,1)=max(max(abs(a_fhat-fhat.mean)./fhat.sd));
    z_fhat_linear_AUC(b,1)=max(max(abs(a_fhat_linear_AUC-fhat_linear.AUC_mean)./fhat_linear.AUC_sd));
    z_fhat_AUC(b,1)=max(max(abs(a_fhat_AUC-fhat.AUC_mean)./fhat.AUC_sd));
end;

c_fhat_linear=quantile(z_fhat_linear,0.95);
c_fhat=quantile(z_fhat,0.95);
c_fhat_linear_AUC=quantile(z_fhat_linear_AUC,0.95);
c_fhat_AUC=quantile(z_fhat_AUC,0.95);

fhat_linear.sQ025=fhat_linear.mean-c_fhat_linear*fhat_linear.sd;
fhat_linear.sQ975=fhat_linear.mean+c_fhat_linear*fhat_linear.sd;
fhat.sQ025=fhat.mean-c_fhat*fhat.sd;
fhat.sQ975=fhat.mean+c_fhat*fhat.sd;
fhat_linear.sAUC_Q025=fhat_linear.AUC_mean-c_fhat_linear_AUC*fhat_linear.AUC_sd;
fhat_linear.sAUC_Q975=fhat_linear.AUC_mean+c_fhat_linear_AUC*fhat_linear.AUC_sd;
fhat.sAUC_Q025=fhat.AUC_mean-c_fhat_AUC*fhat.AUC_sd;
fhat.sAUC_Q975=fhat.AUC_mean+c_fhat_AUC*fhat.AUC_sd;
toc

save 'fhat_summary_all.mat' fhat_linear fhat model

%%%% Now construct derivative curves

fdhat.mean=zeros(ngrid,T);
fdhat.SS=zeros(ngrid,T);
fdhat.sd=zeros(ngrid,T);
fdhat.AUC_mean=zeros(n_age,T);
fdhat.AUC_SS=zeros(n_age,T);
fdhat.AUC_sd=zeros(n_age,T);
fdhat_linear=fdhat;


for b=1:B %587.725975 sec
    a_fdhat_linear=Xd_expand*reshape(MCMC_g(b,:),T,model.p)';   %%% age linear effect
    a_fdhat=a_fdhat_linear+Zd_all*reshape(MCMC_U(b,:),T,m_spline)';  %%% age spline_effect
    fdhat_linear.mean=fdhat_linear.mean+a_fdhat_linear/B;
    fdhat.mean=fdhat.mean+a_fdhat/B;
    fdhat_linear.SS=fdhat_linear.SS+a_fdhat_linear.^2;
    fdhat.SS=fdhat.SS+a_fdhat.^2;
    % Now AUC
    fdhat_linear.AUC_mean=fdhat_linear.AUC_mean+AUC_wgts_all*a_fdhat_linear/B;
    fdhat.AUC_mean=fdhat.AUC_mean+AUC_wgts_all*a_fdhat/B;
    fdhat_linear.AUC_SS=fdhat_linear.AUC_SS+(AUC_wgts_all*a_fdhat_linear).^2;
    fdhat.AUC_SS=fdhat.AUC_SS+(AUC_wgts_all*a_fdhat).^2;
 
end;

fdhat_linear.sd=((B-1)^(-1)*(fdhat_linear.SS-B*(fdhat_linear.mean).^2)).^(1/2);
fdhat.sd=((B-1)^(-1)*(fdhat.SS-B*(fdhat.mean).^2)).^(1/2);
fdhat_linear.Q025=fdhat_linear.mean-1.96*fdhat_linear.sd;
fdhat_linear.Q975=fdhat_linear.mean+1.96*fdhat_linear.sd;
fdhat.Q025=fdhat.mean-1.96*fdhat.sd;
fdhat.Q975=fdhat.mean+1.96*fdhat.sd;
%%% Now AUC
fdhat_linear.AUC_sd=((B-1)^(-1)*(fdhat_linear.AUC_SS-B*(fdhat_linear.AUC_mean).^2)).^(1/2);
fdhat.AUC_sd=((B-1)^(-1)*(fdhat.AUC_SS-B*(fdhat.AUC_mean).^2)).^(1/2);
fdhat_linear.AUC_Q025=fdhat_linear.AUC_mean-1.96*fdhat_linear.AUC_sd;
fdhat_linear.AUC_Q975=fdhat_linear.AUC_mean+1.96*fdhat_linear.AUC_sd;
fdhat.AUC_Q025=fdhat.AUC_mean-1.96*fdhat.AUC_sd;
fdhat.AUC_Q975=fdhat.AUC_mean+1.96*fdhat.AUC_sd;
toc 

%%% Now compute 95% simultaneous credible bands

z_fdhat_linear=zeros(B,1);
z_fdhat=zeros(B,1);
z_fdhat_linear_AUC=zeros(B,1);
z_fdhat_AUC=zeros(B,1);


for b=1:B %Elapsed time is 408.233124 sec 
    a_fdhat_linear=X_expand*reshape(MCMC_g(b,:),T,model.p)';   %%% age linear effect
    a_fdhat=a_fdhat_linear+Z_all*reshape(MCMC_U(b,:),T,m_spline)';  %%% age spline_effect
    a_fdhat_linear_AUC=AUC_wgts_all*a_fdhat_linear;
    a_fdhat_AUC=AUC_wgts_all*a_fdhat;
    
    z_fdhat_linear(b,1)=max(max(abs(a_fdhat_linear-fdhat_linear.mean)./fdhat_linear.sd));
    z_fdhat(b,1)=max(max(abs(a_fdhat-fdhat.mean)./fdhat.sd));
    z_fdhat_linear_AUC(b,1)=max(max(abs(a_fdhat_linear_AUC-fdhat_linear.AUC_mean)./fdhat_linear.AUC_sd));
    z_fdhat_AUC(b,1)=max(max(abs(a_fdhat_AUC-fdhat.AUC_mean)./fdhat.AUC_sd));
    
end;

c_fdhat_linear=quantile(z_fdhat_linear,0.975);
c_fdhat=quantile(z_fdhat,0.975);
c_fdhat_linear_AUC=quantile(z_fdhat_linear_AUC,0.975);
c_fdhat_AUC=quantile(z_fdhat_AUC,0.975);

fdhat_linear.sQ025=fdhat_linear.mean-c_fdhat_linear*fdhat_linear.sd;
fdhat_linear.sQ975=fdhat_linear.mean+c_fdhat_linear*fdhat_linear.sd;
fdhat.sQ025=fdhat.mean-c_fdhat*fdhat.sd;
fdhat.sQ975=fdhat.mean+c_fdhat*fdhat.sd;
fdhat_linear.sAUC_Q025=fdhat_linear.AUC_mean-c_fdhat_linear_AUC*fdhat_linear.AUC_sd;
fdhat_linear.sAUC_Q975=fdhat_linear.AUC_mean+c_fdhat_linear_AUC*fdhat_linear.AUC_sd;
fdhat.sAUC_Q025=fdhat.AUC_mean-c_fdhat_AUC*fdhat.AUC_sd;
fdhat.sAUC_Q975=fdhat.AUC_mean+c_fdhat_AUC*fdhat.AUC_sd;
toc 

save 'fdhat_summary_all.mat' fdhat_linear fdhat model

clear

%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3
%%%%%%%%%%%%%%%%%%%%%%%%

load('DCompressed.mat')
load('Y_simulated.mat');
load('eye_filtered_wave_fmm\eye_filtered_wave_fmm.mat','model');

load('Zg_spline.txt')
load('Zdg_spline.txt')
load('MCMC_U_all.mat','MCMC_U')
load('MCMC_g_low_pass.mat')

%%%%%% Aggregated AUC


%%%%%% Set up evaluation design matrices and etc. 

X_age_unique=(20:90)';
X_pressure_unique=[7,10,15,20,25,30,35,40,45]';
n_eye = size(model.X,1)/9;
unique_idx = n_eye:n_eye:size(model.X,1);
X_IOP_unique= model.X(unique_idx,3:4);

X_age_all=kron(ones(length(X_pressure_unique),1),X_age_unique);
X_IOP_all=kron(X_IOP_unique,ones(length(X_age_unique),1));
X_expand=[ones(size(X_age_all)),X_age_all,X_IOP_all];
Xd_expand=[zeros(size(X_age_all)),ones(size(X_age_all)),zeros(size(X_IOP_all))];

Z_unique=Zg_spline;
Zd_unique=Zdg_spline; %%% derivatives

Z_all=kron(ones(length(X_IOP_unique),1),Z_unique);
Zd_all=kron(ones(length(X_IOP_unique),1),Zd_unique);

n_age=length(X_age_unique);
n_pressure=size(X_IOP_unique,1);
ngrid=n_age*n_pressure;

model.H=4;
H_spline=1;
H_0=model.H-H_spline;
h_spline=(H_0+1):model.H; %the indices of random effects levels that are splines, assumed to be last
m_spline=sum(model.m(h_spline));

B=size(MCMC_U,1);
model.p=size(model.X,2);
p=model.p;


%%%% Now compute feature summaries for circumferenctial location 

circum_indicator=kron(ones(120,1),eye(120))/120;

MCMC_beta_circum=MCMC_g*kron(eye(p),circum_indicator);
MCMC_U_circum=MCMC_U*kron(eye(m_spline),circum_indicator);

MCMC_fhat_linear_circum=MCMC_beta_circum*kron(X_expand',eye(120));
MCMC_fhat_circum=MCMC_fhat_linear_circum+MCMC_U_circum*kron(Z_all',eye(120));

AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
AUC_wgts_all=kron(AUC_wgts,eye(length(X_age_unique)));

MCMC_AUC_circum=MCMC_fhat_circum*kron(AUC_wgts_all',eye(120));

fhat_circum_mean = reshape(mean(MCMC_AUC_circum),120,71)'; 
save 'fhat_circum_mean.mat' fhat_circum_mean


%%%% Now compute feature summaries for each "section" x "region"
region_labels={'PP','MP'};
section_labels={'T','ST','S','SN','N','IN','I','IT'};
region_indicator=kron(eye(2),ones(1,60))/60;
section_indicator=kron(eye(8),ones(1,15))/15;
Trs=size(region_indicator,1)*size(section_indicator,1);

Y_section=Y_simulated*kron(section_indicator',region_indicator');
Y_region=Y_simulated*kron(ones(120,1)/120,region_indicator');
n_eyes=34;
AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
AUC_wgts_all=kron(AUC_wgts,eye(n_eyes));
AUC_raw_section=AUC_wgts_all*Y_section;
AUC_raw_region=AUC_wgts_all*Y_region;

save 'region_section_results.mat' Y_section Y_region AUC_raw_section AUC_raw_region region_labels section_labels;

Y_mean=IDWT_rows(mean(D_compressed),wavespecs_compressed);

Y_mean_section=repmat(Y_mean*kron(section_indicator',region_indicator'),B,1);
Y_mean_region=repmat(Y_mean*kron(ones(120,1)/120,region_indicator'),B,1);

MCMC_beta_section=MCMC_g*kron(eye(p),kron(section_indicator',region_indicator'));
MCMC_U_section=MCMC_U*kron(eye(m_spline),kron(section_indicator',region_indicator'));
MCMC_beta_region=MCMC_g*kron(eye(p),kron(ones(120,1)/120,region_indicator'));
MCMC_U_region=MCMC_U*kron(eye(m_spline),kron(ones(120,1)/120,region_indicator'));

MCMC_fhat_linear_region=MCMC_beta_region*kron(X_expand',eye(2));
MCMC_fhat_linear_section=MCMC_beta_section*kron(X_expand',eye(16));
MCMC_fhat_region=MCMC_fhat_linear_region+MCMC_U_region*kron(Z_all',eye(2));
MCMC_fhat_section=MCMC_fhat_linear_section+MCMC_U_section*kron(Z_all',eye(16));

AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
AUC_wgts_all=kron(AUC_wgts,eye(length(X_age_unique)));

MCMC_AUC_linear_region=MCMC_fhat_linear_region*kron(AUC_wgts_all',eye(2));
MCMC_AUC_linear_section=MCMC_fhat_linear_section*kron(AUC_wgts_all',eye(16));
MCMC_AUC_region=MCMC_fhat_region*kron(AUC_wgts_all',eye(2));
MCMC_AUC_section=MCMC_fhat_section*kron(AUC_wgts_all',eye(16));

fhat_linear_region.mean=reshape(mean(MCMC_fhat_linear_region),2,639)';
fhat_linear_region.sd=reshape(sqrt(var(MCMC_fhat_linear_region)),2,639)';
fhat_linear_region.Q025=reshape(quantile(MCMC_fhat_linear_region,.025),2,639)';
fhat_linear_region.Q975=reshape(quantile(MCMC_fhat_linear_region,.975),2,639)';
[upper_CI, lower_CI] = jointband(MCMC_fhat_linear_region,0.05);
fhat_linear_region.sQ025=reshape(lower_CI,2,639)';
fhat_linear_region.sQ975=reshape(upper_CI,2,639)';

fhat_region.mean=reshape(mean(MCMC_fhat_region),2,639)';
fhat_region.sd=reshape(sqrt(var(MCMC_fhat_region)),2,639)';
fhat_region.Q025=reshape(quantile(MCMC_fhat_region,.025),2,639)';
fhat_region.Q975=reshape(quantile(MCMC_fhat_region,.975),2,639)';
[upper_CI, lower_CI] = jointband(MCMC_fhat_region,0.05);
fhat_region.sQ025=reshape(lower_CI,2,639)';
fhat_region.sQ975=reshape(upper_CI,2,639)';

AUC_linear_region.mean=reshape(mean(MCMC_AUC_linear_region),2,71)';
AUC_linear_region.sd=reshape(sqrt(var(MCMC_AUC_linear_region)),2,71)';
AUC_linear_region.Q025=reshape(quantile(MCMC_AUC_linear_region,.025),2,71)';
AUC_linear_region.Q975=reshape(quantile(MCMC_AUC_linear_region,.975),2,71)';
[upper_CI, lower_CI] = jointband(MCMC_AUC_linear_region,0.05);
AUC_linear_region.sQ025=reshape(lower_CI,2,71)';
AUC_linear_region.sQ975=reshape(upper_CI,2,71)';

AUC_region.mean=reshape(mean(MCMC_AUC_region),2,71)';
AUC_region.sd=reshape(sqrt(var(MCMC_AUC_region)),2,71)';
AUC_region.Q025=reshape(quantile(MCMC_AUC_region,.025),2,71)';
AUC_region.Q975=reshape(quantile(MCMC_AUC_region,.975),2,71)';
[upper_CI, lower_CI] = jointband(MCMC_AUC_region,0.05);
AUC_region.sQ025=reshape(lower_CI,2,71)';
AUC_region.sQ975=reshape(upper_CI,2,71)';

fhat_linear_section.mean=reshape(mean(MCMC_fhat_linear_section),16,639)';
fhat_linear_section.sd=reshape(sqrt(var(MCMC_fhat_linear_section)),16,639)';
fhat_linear_section.Q025=reshape(quantile(MCMC_fhat_linear_section,.025),16,639)';
fhat_linear_section.Q975=reshape(quantile(MCMC_fhat_linear_section,.975),16,639)';
[upper_CI, lower_CI] = jointband(MCMC_fhat_linear_section,0.05);
fhat_linear_section.sQ025=reshape(lower_CI,16,639)';
fhat_linear_section.sQ975=reshape(upper_CI,16,639)';

fhat_section.mean=reshape(mean(MCMC_fhat_section),16,639)';
fhat_section.sd=reshape(sqrt(var(MCMC_fhat_section)),16,639)';
fhat_section.Q025=reshape(quantile(MCMC_fhat_section,.025),16,639)';
fhat_section.Q975=reshape(quantile(MCMC_fhat_section,.975),16,639)';
[upper_CI, lower_CI] = jointband(MCMC_fhat_section,0.05);
fhat_section.sQ025=reshape(lower_CI,16,639)';
fhat_section.sQ975=reshape(upper_CI,16,639)';

AUC_linear_section.mean=reshape(mean(MCMC_AUC_linear_section),16,71)';
AUC_linear_section.sd=reshape(sqrt(var(MCMC_AUC_linear_section)),16,71)';
AUC_linear_section.Q025=reshape(quantile(MCMC_AUC_linear_section,.025),16,71)';
AUC_linear_section.Q975=reshape(quantile(MCMC_AUC_linear_section,.975),16,71)';
[upper_CI, lower_CI] = jointband(MCMC_AUC_linear_section,0.05);
AUC_linear_section.sQ025=reshape(lower_CI,16,71)';
AUC_linear_section.sQ975=reshape(upper_CI,16,71)';

AUC_section.mean=reshape(mean(MCMC_AUC_section),16,71)';
AUC_section.sd=reshape(sqrt(var(MCMC_AUC_section)),16,71)';
AUC_section.Q025=reshape(quantile(MCMC_AUC_section,.025),16,71)';
AUC_section.Q975=reshape(quantile(MCMC_AUC_section,.975),16,71)';
[upper_CI, lower_CI] = jointband(MCMC_AUC_section,0.05);
AUC_section.sQ025=reshape(lower_CI,16,71)';
AUC_section.sQ975=reshape(upper_CI,16,71)';

save 'fhat_summary_all_region_section.mat' fhat_linear_region fhat_region AUC_linear_region AUC_region fhat_linear_section fhat_section AUC_linear_section AUC_section

%%%%%% Now compute section/region summaries of derivatives

MCMC_fdhat_linear_region=MCMC_beta_region*kron(Xd_expand',eye(2));
MCMC_fdhat_linear_section=MCMC_beta_section*kron(Xd_expand',eye(16));
MCMC_fdhat_region=MCMC_fdhat_linear_region+MCMC_U_region*kron(Zd_all',eye(2));
MCMC_fdhat_section=MCMC_fdhat_linear_section+MCMC_U_section*kron(Zd_all',eye(16));

AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
AUC_wgts_all=kron(AUC_wgts,eye(length(X_age_unique)));

MCMC_dAUC_linear_region=MCMC_fdhat_linear_region*kron(AUC_wgts_all',eye(2));
MCMC_dAUC_linear_section=MCMC_fdhat_linear_section*kron(AUC_wgts_all',eye(16));
MCMC_dAUC_region=MCMC_fdhat_region*kron(AUC_wgts_all',eye(2));
MCMC_dAUC_section=MCMC_fdhat_section*kron(AUC_wgts_all',eye(16));

fdhat_linear_region.mean=reshape(mean(MCMC_fdhat_linear_region),2,639)';
fdhat_linear_region.sd=reshape(sqrt(var(MCMC_fdhat_linear_region)),2,639)';
fdhat_linear_region.Q025=reshape(quantile(MCMC_fdhat_linear_region,.025),2,639)';
fdhat_linear_region.Q975=reshape(quantile(MCMC_fdhat_linear_region,.975),2,639)';
[upper_CI, lower_CI] = jointband(MCMC_fdhat_linear_region,0.05);
fdhat_linear_region.sQ025=reshape(lower_CI,2,639)';
fdhat_linear_region.sQ975=reshape(upper_CI,2,639)';

fdhat_region.mean=reshape(mean(MCMC_fdhat_region),2,639)';
fdhat_region.sd=reshape(sqrt(var(MCMC_fdhat_region)),2,639)';
fdhat_region.Q025=reshape(quantile(MCMC_fdhat_region,.025),2,639)';
fdhat_region.Q975=reshape(quantile(MCMC_fdhat_region,.975),2,639)';
[upper_CI, lower_CI] = jointband(MCMC_fdhat_region,0.05);
fdhat_region.sQ025=reshape(lower_CI,2,639)';
fdhat_region.sQ975=reshape(upper_CI,2,639)';

dAUC_linear_region.mean=reshape(mean(MCMC_dAUC_linear_region),2,71)';
dAUC_linear_region.sd=reshape(sqrt(var(MCMC_dAUC_linear_region)),2,71)';
dAUC_linear_region.Q025=reshape(quantile(MCMC_dAUC_linear_region,.025),2,71)';
dAUC_linear_region.Q975=reshape(quantile(MCMC_dAUC_linear_region,.975),2,71)';
[upper_CI, lower_CI] = jointband(MCMC_dAUC_linear_region,0.05);
dAUC_linear_region.sQ025=reshape(lower_CI,2,71)';
dAUC_linear_region.sQ975=reshape(upper_CI,2,71)';

dAUC_region.mean=reshape(mean(MCMC_dAUC_region),2,71)';
dAUC_region.sd=reshape(sqrt(var(MCMC_dAUC_region)),2,71)';
dAUC_region.Q025=reshape(quantile(MCMC_dAUC_region,.025),2,71)';
dAUC_region.Q975=reshape(quantile(MCMC_dAUC_region,.975),2,71)';
[upper_CI, lower_CI] = jointband(MCMC_dAUC_region,0.05);
dAUC_region.sQ025=reshape(lower_CI,2,71)';
dAUC_region.sQ975=reshape(upper_CI,2,71)';

fdhat_linear_section.mean=reshape(mean(MCMC_fdhat_linear_section),16,639)';
fdhat_linear_section.sd=reshape(sqrt(var(MCMC_fdhat_linear_section)),16,639)';
fdhat_linear_section.Q025=reshape(quantile(MCMC_fdhat_linear_section,.025),16,639)';
fdhat_linear_section.Q975=reshape(quantile(MCMC_fdhat_linear_section,.975),16,639)';
[upper_CI, lower_CI] = jointband(MCMC_fdhat_linear_section,0.05);
fdhat_linear_section.sQ025=reshape(lower_CI,16,639)';
fdhat_linear_section.sQ975=reshape(upper_CI,16,639)';

fdhat_section.mean=reshape(mean(MCMC_fdhat_section),16,639)';
fdhat_section.sd=reshape(sqrt(var(MCMC_fdhat_section)),16,639)';
fdhat_section.Q025=reshape(quantile(MCMC_fdhat_section,.025),16,639)';
fdhat_section.Q975=reshape(quantile(MCMC_fdhat_section,.975),16,639)';
[upper_CI, lower_CI] = jointband(MCMC_fdhat_section,0.05);
fdhat_section.sQ025=reshape(lower_CI,16,639)';
fdhat_section.sQ975=reshape(upper_CI,16,639)';

dAUC_linear_section.mean=reshape(mean(MCMC_dAUC_linear_section),16,71)';
dAUC_linear_section.sd=reshape(sqrt(var(MCMC_dAUC_linear_section)),16,71)';
dAUC_linear_section.Q025=reshape(quantile(MCMC_dAUC_linear_section,.025),16,71)';
dAUC_linear_section.Q975=reshape(quantile(MCMC_dAUC_linear_section,.975),16,71)';
[upper_CI, lower_CI] = jointband(MCMC_dAUC_linear_section,0.05);
dAUC_linear_section.sQ025=reshape(lower_CI,16,71)';
dAUC_linear_section.sQ975=reshape(upper_CI,16,71)';

dAUC_section.mean=reshape(mean(MCMC_dAUC_section),16,71)';
dAUC_section.sd=reshape(sqrt(var(MCMC_dAUC_section)),16,71)';
dAUC_section.Q025=reshape(quantile(MCMC_dAUC_section,.025),16,71)';
dAUC_section.Q975=reshape(quantile(MCMC_dAUC_section,.975),16,71)';
[upper_CI, lower_CI] = jointband(MCMC_dAUC_section,0.05);
dAUC_section.sQ025=reshape(lower_CI,16,71)';
dAUC_section.sQ975=reshape(upper_CI,16,71)';

save 'fdhat_summary_all_region_section.mat' fdhat_linear_region fdhat_region dAUC_linear_region dAUC_region fdhat_linear_section fdhat_section dAUC_linear_section dAUC_section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part 4 - Compute induced Covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Getting the transform matrix to project the basis values back to the
%%%%% original space
clear

[Phi_1,Kj1] = Get_DWT('db3',120,'per       ',0,5);
[Phi_2,Kj2] = Get_DWT('db3',120,'reflection',1,5);
Phi_all = kron(Phi_1,Phi_2); % DWT matrix 
IDWT = kron((Phi_1'*Phi_1)\Phi_1',(Phi_2'*Phi_2)\Phi_2');

load('DCompressed.mat')
load('eye_filtered_wave_fmm\eye_filtered_wave_fmm.mat','model');

K = sum(wavespecs_compressed.Kj);
[N, P] = size(model.X);
H = size(model.m,2);
B = 1000;

filename2 = strcat('eye_filtered_wave_fmm\eye_filtered_wave_fmm_Results_0_theta.dat');
MCMC_theta = readbeta(B,H+1,K,filename2);
btheta=reshape(mean(MCMC_theta),K,H+1)'; 
btheta=btheta([1:3,5],:);% taking out VC for spline 

Tstar = sum(wavespecs_compressed.Kj_all);
T = 14400;
sig_star = zeros(H, Tstar); 
for h=1:H
    sig_tmp = zeros(Tstar,1);
    sig_tmp(wavespecs_compressed.keep==1) = btheta(h,:);   %%% insert nonzero wavelets in MCMC order
    sig_tmp(wavespecs_compressed.reorder) = sig_tmp;       %%% switches to tensor order
    sig_star(h,:) = sig_tmp;
end
%%% Get reordered "keep" vector to match order in IDWT matrix
%%%  This is a matrix of indicators of which wavelets were seleected after
%%%  compression in the order that matches the IDWT matrix.
%%%  wavespecs.keep is in the order of coefficients from DWT_rows with
%%%  rectangular 2d transform which is the same as that used for the MCMC
keep_reorder=wavespecs_compressed.keep;                    %%% non-thresholded parameters in MCMC order
keep_reorder(wavespecs_compressed.reorder)=keep_reorder;  %%% non-thresholded parameters in tensor order

%%% Test to be sure zeroed out basis levels match in keep and sig_star
thresh=(sum(sig_star==0)>0);
norm(1-thresh-keep_reorder)

tic;
mvar2 = zeros(T,H);
for h=1:H
    tmp_sig = sig_star(h,keep_reorder==1)';
    for t=1:T
        mvar2(t,h) = (IDWT(t,keep_reorder==1).^2)*tmp_sig;
    end
end
toc;

save('InducedCovariances.mat','mvar2','sig_star')

#Geweke diagnostic statistics for the MCMC samples


rm(list=ls())
require(coda)
require(R.matlab)

#############################################
####### Diagnostics for fixed effects samples 
#############################################

# Fixed Effects on basis space

MCMC_beta_basis=readMat("MCMC_beta.mat")
DataMCMC=MCMC_beta_basis$MCMC.beta
Gew_test_beta_basis=NULL
EffectiveSize_beta_basis=NULL
K=dim(DataMCMC)[2]/4
for (j in 1:4){
  for(i in ((j-1)*K+1):(j*K)){
    aa=geweke.diag(DataMCMC[,i], frac1=0.25, frac2=0.25)
    Gew_test_beta_basis[i]=aa$z
    EffectiveSize_beta_basis[i]=effectiveSize(DataMCMC[,i])
  }
}


# Nonparametric age effect on basis space 
MCMC_beta_nonparametric_basis=readMat("MCMC_fstar.mat")
BasisMCMC_nonparametric=MCMC_beta_nonparametric_basis$MCMC.fstar
Gew_test_nonparametric_basis=NULL
EffectiveSize_nonparametric_basis=NULL
for (j in 1:15){
  for(i in ((j-1)*K+1):(j*K)){
    aa=geweke.diag(BasisMCMC_nonparametric[,i], frac1=0.25, frac2=0.25)
    Gew_test_nonparametric_basis[i]=aa$z
    EffectiveSize_nonparametric_basis[i]=effectiveSize(BasisMCMC_nonparametric[,i])
  }
}



#####################################################
# Diagnostics for variance components 
#####################################################

# Variance components in the basis space

MCMC_theta0=readMat("MCMC_theta.mat")
DataMCMC=MCMC_theta0$MCMC.theta
Gew_test=NULL
EffectiveSize=NULL
for (j in 1:5){
  for(i in ((j-1)*K+1):(j*K)){
    aa=geweke.diag(DataMCMC[,i], frac1=0.25, frac2=0.25)
    Gew_test[i]=aa$z
    EffectiveSize[i]=effectiveSize(DataMCMC[,i])
  }
}


writeMat("ResultsDiag_FixedEffects.mat",Geweke_nonparametric_basis=Gew_test_nonparametric_basis, EffectiveSize_nonparametric_basis=EffectiveSize_nonparametric_basis,Geweke_beta_basis=Gew_test_beta_basis, EffectiveSize_beta_basis=EffectiveSize_beta_basis)
writeMat("ResultsDiag_Theta.mat",Geweke_thetabasis=Gew_test, EffectiveSize_thetabasis=EffectiveSize)




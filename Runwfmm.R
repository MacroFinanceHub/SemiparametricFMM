#setwd("C:/Users/mfmiranda/Documents/MATLAB/EYE/FinalScripts/SimulatedDataReproducibility")


require(R.matlab)
require(nlme)
require(MASS)
library(wfmm)

#########################################################################
#### Create covariates from sample information
#########################################################################

D=readMat("DCompressed.mat")
D_wave=D$D
wavespecs=D$wavespecs
sample_info=read.table("sample_info2.txt")
attributes(sample_info)$names=c("MPS_PP","MPS_MP","age","IOP","subject","eye")

#########################################################################
#### Create covariates from sample information
#########################################################################

#### create IOP1 and IOP2 to represent hyperbola
IOP=sample_info$IOP

IOPinv=IOP^(-1)
IOPstd=(IOP-mean(IOP))/sqrt(var(IOP))
IOPinvstd=(IOPinv-mean(IOPinv))/sqrt(var(IOPinv))
IOP1=sqrt(2)/2*IOPstd-sqrt(2)/2*IOPinvstd
IOP2=sqrt(2)/2*IOPstd+sqrt(2)/2*IOPinvstd
sample_info$IOP1=IOP1
sample_info$IOP2=IOP2

#### creat eyeID and eye indicator(left vs right)
N=dim(sample_info)[1]
sample_info$eyeID=rep(1:(N/9),9)
sample_info$one=rep(1,N)
sample_info$eye = 1*(sample_info$eye==1) # transform eye 1 and 2 to eye 1 and 0 (indicator)

##########################################################################
#### Get Spline Z matrices for age and IOP
##########################################################################

source("GetSpline.R")

#### Spline matrix for age

numIntKnots=5
age_unique = unique(sample_info$age)
a=min(age_unique)
b=max(age_unique)
ngrid=length(a:b)

SplineResult<-GetSpline(sample_info$age,numIntKnots,ngrid,a,b)
Zres=SplineResult$Zres
Z=Zres[1:N,]
Z1=Z*IOP1
Z2=Z*IOP2
sample_info$Z=Z
sample_info$Z1=Z1
sample_info$Z2=Z2

#### Spline matrix for IOP

SplineResult_IOP<-GetSpline(sample_info$IOP,numIntKnots,length(7:45),7,45)
Zres_IOP=SplineResult_IOP$Zres
Z_IOP=Zres_IOP[1:N,]
sample_info$Z_IOP=Z_IOP

B_iop=SplineResult_IOP$B
Omega_iop=SplineResult_IOP$Omega

##########################################################################

##########################################################################
#### This tables are for computing DF for nonparametric effect
#### Spline basis evaluation and penalty matrix for computing DF
##########################################################################

B=SplineResult$B
Omega=SplineResult$Omega

#########################################################################
#### Run FMM
#########################################################################

eye.data.filtered.wave <- D_wave

### Format data and set up 'fmmBase' object

# Format data using prep_fmm_data
eye.data.filtered.wave <- as_fmm_data(data = eye.data.filtered.wave,
                                      dim.names = list(cols = 1:ncol(eye.data.filtered.wave)), 
                                      n.dims = 1)

# Put into an 'fmmBase' object together with variables to call 'fmm' upon 
eye.base.filtered.wave <- fmm_base(Y = eye.data.filtered.wave, 
                                   variables = sample_info,
                                   labels = list(cols = "index"))



### Fit an FMM

# Get minimum/maximum of age for evaluation grid of non-parametric effect
a <- min(eye.base.filtered.wave$variables$age)
b <- max(eye.base.filtered.wave$variables$age)

# Call fmm
# NOTES:
# - If you don't want/need the random effects sampled, set PostProcessSpecs$compute_U & keep_U_samples to 0
# - You can select another directory than the current working directory for the model folder to be written
#   using the 'path' argument
# - Adapt the MCMC settings using MCMCspecs
# - If memory may become an issue, you can use 'ff' functionality implemented the function
#   for batchwise processing of posteriors. Relevant arguments of 'fmm' are 'ff.mode' and 'batchMBs'.
# - During the postprocessing, the function quit with an error -- I'm gonna look into this.

eye.filtered.wave.fmm <-
  fmm(x = eye.base.filtered.wave,
      formula = ~ 1 + age + IOP1 + IOP2 + (1 + IOP1 + IOP2 | eyeID), 
      non.parametric = list(age = list(n.int.knots = 5, n.grid = length(a:b), 
                                       eval.from = a, aeval.to = b),
                            evaluate = FALSE),
      disk.id = "eye.filtered.wave.fmm",
      # path = ".",
      basis_specs = list(transformtype = "none",partitions = wavespecs[[8]]),
      PostProcessSpecs = list(compute_U = 0, keep_U_samples = 0),
      MCMCspecs = list(B = 1000, burnin = 5000, thin = 10),
      compute.summaries = FALSE)







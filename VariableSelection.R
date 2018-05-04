#setwd("C:/Eye_Revision/BIC correction")
#setwd("/Users/miranda/Documents/MATLAB/EYE/PaperScripts/SimulatedData_reproducibility")

require(R.matlab)
require(nlme)
require(MASS)


#########################################################################
#### Create covariates from sample information
#########################################################################

D=readMat("DCompressed.mat")
D_wave=D$D.compressed
wavespecs=D$wavespecs.compressed
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

B=SplineResult$B
Omega=SplineResult$Omega

#### Spline matrix for IOP

SplineResult_IOP<-GetSpline(sample_info$IOP,numIntKnots,length(7:45),7,45)
Zres_IOP=SplineResult_IOP$Zres
Z_IOP=Zres_IOP[1:N,]
sample_info$Z_IOP=Z_IOP

B_iop=SplineResult_IOP$B
Omega_iop=SplineResult_IOP$Omega

##########################################################################
#### Model selection using BIC
##########################################################################

#### Fixed effect selection
main_effects_selection = function(basis_input){
  
  K = dim(basis_input)[2]
  input_data = sample_info
  
  n.model = 12
  BIC.all = matrix(1e20,dim(basis_input)[2],n.model)
  AIC.all = matrix(1e20,dim(basis_input)[2],n.model)
  
  # funtion to compute effective degree of freedom (DF) for non-linear effect
  DF_compute <- function(B,Omega,var_e,var_s){
    BtB = t(B)%*%B
    X = ginv(BtB+Omega*var_e/var_s)%*%BtB
    return(sum(diag(X)))
  }
  
  for (k in 1:K){
    input_data$MPS=basis_input[,k]
    dataFr1=groupedData(MPS~age|one,data=input_data)
    
    m1<-try(lme(MPS~age+IOP,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z_IOP-1))))
    m2<-try(lme(MPS~age+IOP1+IOP2,dataFr1,random=list(one=pdIdent(~Z-1))))
    m3<-try(lme(MPS~age+IOP,dataFr1,random=list(one=pdIdent(~Z-1))))
    m4<-try(lme(MPS~age+IOP,dataFr1,random=list(one=pdIdent(~Z_IOP-1))))
    m5<-try(lme(MPS~age+IOP1+IOP2,dataFr1))
    m6<-try(lme(MPS~age+IOP,dataFr1))
    
    m7<-try(lme(MPS~eye+age+IOP,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z_IOP-1))))
    m8<-try(lme(MPS~eye+age+IOP1+IOP2,dataFr1,random=list(one=pdIdent(~Z-1))))
    m9<-try(lme(MPS~eye+age+IOP,dataFr1,random=list(one=pdIdent(~Z-1))))
    m10<-try(lme(MPS~eye+age+IOP,dataFr1,random=list(one=pdIdent(~Z_IOP-1))))
    m11<-try(lme(MPS~eye+age+IOP1+IOP2,dataFr1))
    m12<-try(lme(MPS~eye+age+IOP,dataFr1))
    
    models = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
    for(m in 1:n.model){
      temp_model = models[[m]]
      if(class(temp_model) == "try-error") {print(paste('basis',k,',m',m,':Convergence fail',sep='')); next;
      }else{
        BIC.tmp = BIC(temp_model)
        AIC.tmp = AIC(temp_model)

        Variances = VarCorr(temp_model)
        n = nobs(temp_model)
        if(m%in%c(1,7)){
          DF.age = DF_compute(B,Omega,as.numeric(Variances[17,1]),as.numeric(Variances[2,1]))
          DF.iop = DF_compute(B_iop,Omega_iop,as.numeric(Variances[17,1]),as.numeric(Variances[10,1]))
          BIC.tmp = BIC.tmp + (DF.age-3)*log(n) + (DF.iop-3)*log(n) 
          AIC.tmp = AIC.tmp + (DF.age-3)*2 + (DF.iop-3)*2
        }else if(m%in%c(2,3,8,9)){
          DF.age = DF_compute(B,Omega,as.numeric(Variances[8,1]),as.numeric(Variances[1,1]))
          BIC.tmp = BIC.tmp + (DF.age-3)*log(n) 
          AIC.tmp = AIC.tmp + (DF.age-3)*2 
        }else if(m%in%c(4,10)){
          DF.iop = DF_compute(B_iop,Omega_iop,as.numeric(Variances[8,1]),as.numeric(Variances[1,1]))
          BIC.tmp = BIC.tmp + (DF.iop-3)*log(n) 
          AIC.tmp = AIC.tmp + (DF.iop-3)*2           
        }
        
        BIC.all[k,m]=BIC.tmp
        AIC.all[k,m]=AIC.tmp
        }
    }
    if(k%%100==0){cat("basis",k,"done","\n")}
    
  }
  
  wgt=apply(basis_input,2,function(x)sum(x^2))
  wgt=wgt/sum(wgt)
  BIC.min.ind = matrix(0,k,n.model)
  AIC.min.ind = matrix(0,k,n.model)
  BIC.min.idx = apply(BIC.all,1,which.min)
  AIC.min.idx = apply(AIC.all,1,which.min)
  
  for(j in 1:dim(basis_input)[2]){
    BIC.min.ind[j,BIC.min.idx[j]]=1
    AIC.min.ind[j,AIC.min.idx[j]]=1
  }
  
  BIC.summary=apply(BIC.min.ind*wgt,2,sum)
  AIC.summary=apply(AIC.min.ind*wgt,2,sum)
  opt_idx = which.max(BIC.summary)
  opt_idx_aic = which.max(AIC.summary)
  
  ## Print Results
  output <- data.frame("Age Effect"=rep(rep(c("Nonparametric","Linear"),each=3),2),
             "IOP Effect"=rep(c("Nonparametric","Hyperbolic","Linear"),4),
             "Eye Effect"=rep(c("No","Yes"),each=6), 
             "BIC Score"=sprintf("%.4f",BIC.summary),
             "AIC Score"=sprintf("%.4f",AIC.summary))
  
  sink("model_selection_output.txt")
  cat("Step 1. Selection of Fixed Effects","\n","\n")
  print(output, row.names = FALSE)
  cat("\n", "Best fixed effect model based on BIC:","\n","\n")
  print(output[opt_idx,], row.names = FALSE)
  cat("\n")
  sink()
  return(list(model=opt_idx,score=BIC.summary,score_aic=AIC.summary,model_aix=opt_idx_aic))
}
main_effect_result = main_effects_selection(D_wave) 

#### Interaction term selection
interaction_selection = function(basis_input,main_model=2){
  
  K = dim(basis_input)[2]
  input_data = sample_info
  
  n.model = 4 
  BIC.all=matrix(1e20,dim(basis_input)[2],n.model)
  AIC.all=matrix(1e20,dim(basis_input)[2],n.model)
  
  # funtion to compute effective degree of freedom (DF) for non-linear effect
  DF_compute <- function(B,Omega,var_e, var_s){
    X = ginv(t(B)%*%B+Omega*var_e/var_s)%*%t(B)%*%B
    return(sum(diag(X)))
  }
  
  for (k in 1:K){
    input_data$MPS=basis_input[,k]
    dataFr1=groupedData(MPS~age|one,data=input_data)
    
    if(main_model==2){
      m1<-try(lme(MPS~age+IOP1+IOP2,dataFr1,random=list(one=pdIdent(~Z-1))))
      m2<-try(lme(MPS~age+IOP1+IOP2+Z1,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z1-1))))
      m3<-try(lme(MPS~age+IOP1+IOP2+Z2,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z2-1))))
      m4<-try(lme(MPS~age+IOP1+IOP2+Z1+Z2,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z1-1),one=pdIdent(~Z2-1))))
    }
    
    if(main_model==8){
      m1<-try(lme(MPS~eye+age+IOP1+IOP2,dataFr1,random=list(one=pdIdent(~Z-1))))
      m2<-try(lme(MPS~eye+age+IOP1+IOP2+Z1,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z1-1))))
      m3<-try(lme(MPS~eye+age+IOP1+IOP2+Z2,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z2-1))))
      m4<-try(lme(MPS~eye+age+IOP1+IOP2+Z1+Z2,dataFr1,random=list(one=pdIdent(~Z-1),one=pdIdent(~Z1-1),one=pdIdent(~Z2-1))))
    }
    
    models = list(m1,m2,m3,m4)
    for(m in 1:n.model){
      temp_model = models[[m]]
      if(class(temp_model) == "try-error") {print(paste('basis',k,',m',m,':Convergence fail',sep='')); next;
      }else{
        
        BIC.tmp = BIC(temp_model)
        AIC.tmp = AIC(temp_model)
        Variances = VarCorr(temp_model)
        n = nobs(temp_model)
        
        if(m==1){
          DF.age = DF_compute(B,Omega,as.numeric(Variances[8,1]),as.numeric(Variances[1,1]))
        }else if(m%in%c(2,3)){
          DF.age = 2*DF_compute(B,Omega,as.numeric(Variances[17,1]),as.numeric(Variances[2,1]))
        }else if(m==4){
          DF.age = 3*DF_compute(B,Omega,as.numeric(Variances[25,1]),as.numeric(Variances[2,1]))
        }
        BIC.all[k,m] = BIC.tmp + (DF.age-3)*log(n)  
        AIC.all[k,m] = AIC.tmp + (DF.age-3)*2 
      }		
    }
    if(k%%100==0){cat("basis",k,"done","\n")}
  }
  
  wgt=apply(basis_input,2,function(x)sum(x^2))
  wgt=wgt/sum(wgt)
  BIC.min.ind = matrix(0,k,n.model)
  AIC.min.ind = matrix(0,k,n.model)
  BIC.min.idx = apply(BIC.all,1,which.min)
  AIC.min.idx = apply(AIC.all,1,which.min)
  for(j in 1:dim(basis_input)[2]){
    BIC.min.ind[j,BIC.min.idx[j]]=1
    AIC.min.ind[j,AIC.min.idx[j]]=1
  }
  BIC.summary=apply(BIC.min.ind*wgt,2,sum)
  AIC.summary=apply(AIC.min.ind*wgt,2,sum)
  opt_idx = which.max(BIC.summary)
  opt_idx_aic = which.max(AIC.summary)	
  
  ## Print Results
  output <- data.frame("Interaction Term"=c("No Interaction","f(age)*IOP1","f(age)*IOP2","f(age)*(IOP1 + IOP2)"),
                       "BIC Score"=sprintf("%.4f",BIC.summary),
                       "AIC Score"=sprintf("%.4f",AIC.summary))
  
  sink("model_selection_output.txt", append = TRUE)
  cat("Step 2. Selection of Interaction Terms","\n","\n")
  print(output, row.names = FALSE)
  cat("\n", "Best Interaction Term based on BIC:","\n","\n")
  print(output[opt_idx,], row.names = FALSE)
  cat("\n")
  sink()
  
  return(list(model=opt_idx,score=BIC.summary,model2=opt_idx_aic,socre_aic=AIC.summary))
}
interaction_result = interaction_selection(D_wave,main_model=main_effect_result$model) ## Model 1 chosen (BIC score:1)

#### Random effect selection
random_effects_selection = function(basis_input,fixed_model=2){
  
  K = dim(basis_input)[2]
  input_data = sample_info
  
  ## funtion to compute effective degree of freedom (DF) for non-linear effect
  
  # Effective DF when no other random effects
  DF_compute <- function(B,Omega,var_e, var_s){
    X = ginv(t(B)%*%B+Omega*var_e/var_s)%*%t(B)%*%B
    return(sum(diag(X)))
  }
  
  # Effective DF when other random effects are included (See Section 2 of the supplmentary material)
  DF_compute2 <- function(B,Omega,W,var_s){
    BWB = t(B)%*%W%*%B
    X = ginv(BWB+Omega/var_s)%*%BWB
    return(sum(diag(X)))
  }
  
  # Set-up for computing W in DF_compute2 above
  Z_subj = matrix(0,N,max(sample_info$subject))
  for(i in 1:N){
    Z_subj[i,sample_info$subject[i]] = 1
  }
  n.subj = dim(Z_subj)[2]
  
  eyeID_unique = unique(sample_info$eyeID)
  Z_eyeID = c()
  for(j in 1:length(eyeID_unique)){
    Z_eyeID = 1*cbind(Z_eyeID,sample_info$eyeID == eyeID_unique[j])
  }
  n.eye = dim(Z_eyeID)[2]
  
  Z_iop = Z_eyeID*IOP
  Z_iop1 = Z_eyeID*IOP1
  Z_iop2 = Z_eyeID*IOP2
  
  ZZt_subj = Z_subj%*%t(Z_subj)
  ZZt_eye = Z_eyeID%*%t(Z_eyeID)
  ZZt_iop = Z_iop%*%t(Z_iop)
  ZZt_iop1 = Z_iop1%*%t(Z_iop1)
  ZZt_iop2 = Z_iop2%*%t(Z_iop2)
  
  ### Random effects selection 
  n.model = 10
  BIC.all = matrix(1e20,dim(basis_input)[2],n.model)
  AIC.all = matrix(1e20,dim(basis_input)[2],n.model)
  
  for (k in 1:K){
    input_data$MPS=basis_input[,k]
    dataFr1=groupedData(MPS~age|one,data=input_data)
    
    if(fixed_model==2){fixed_formula=MPS~age+IOP1+IOP2}
    if(fixed_model==8){fixed_formula=MPS~eye+age+IOP1+IOP2}
    
    m1<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1))))
    m2<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),eyeID=pdDiag(~1))))
    m3<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),eyeID=pdDiag(~IOP))))
    m4<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),eyeID=pdDiag(~IOP1+IOP2-1))))
    ptm <- proc.time()
    m5<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),eyeID=pdDiag(~IOP1+IOP2))))
    proc.time() - ptm
    
    
    m6<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),subject=pdDiag(~1))))
    m7<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),subject=pdDiag(~1),eyeID=pdDiag(~1))))
    m8<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),subject=pdDiag(~1),eyeID=pdDiag(~IOP))))
    m9<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),subject=pdDiag(~1),eyeID=pdDiag(~IOP1+IOP2-1))))
    m10<-try(lme(fixed_formula,dataFr1,random=list(one=pdIdent(~Z-1),subject=pdDiag(~1),eyeID=pdDiag(~IOP1+IOP2))))
    
    models = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
    
    for(m in (1:n.model)){
      temp_model = models[[m]]
      if(class(temp_model) == "try-error") {print(paste('basis',k,',m',m,':Convergence fail',sep='')); next;
      }else{
        
        BIC.tmp = BIC(temp_model)
        AIC.tmp = AIC(temp_model)
        Variances = VarCorr(temp_model)
        n = nobs(temp_model)
        
        if(m==1){
          DF.age = DF_compute(B,Omega,as.numeric(Variances[8,1]),as.numeric(Variances[1,1]))
        }else if(m==2){
          W = as.numeric(Variances[10,1])*ZZt_eye + diag(as.numeric(Variances[11,1]),N) 
        }else if(m==3){
          W = as.numeric(Variances[10,1])*ZZt_eye + as.numeric(Variances[11,1])*ZZt_iop + diag(as.numeric(Variances[12,1]),N)
        }else if(m==4){
          W = as.numeric(Variances[10,1])*ZZt_iop1 + as.numeric(Variances[11,1])*ZZt_iop2 + diag(as.numeric(Variances[12,1]),N)
        }else if(m==5){
          W = as.numeric(Variances[10,1])*ZZt_eye + as.numeric(Variances[11,1])*ZZt_iop1 + 
              as.numeric(Variances[12,1])*ZZt_iop2 + diag(as.numeric(Variances[13,1]),N)
        }else if(m==6){
          W = as.numeric(Variances[10,1])*ZZt_subj + diag(as.numeric(Variances[11,1]),N)
        }else if(m==7){
          W = as.numeric(Variances[10,1])*ZZt_subj + as.numeric(Variances[12,1])*ZZt_eye+diag(as.numeric(Variances[13,1]),N)
        }else if(m==8){
          W = as.numeric(Variances[10,1])*ZZt_subj + as.numeric(Variances[12,1])*ZZt_eye+
              as.numeric(Variances[13,1])*ZZt_iop + diag(as.numeric(Variances[14,1]),N)
        }else if(m==9){
          W = as.numeric(Variances[10,1])*ZZt_subj + as.numeric(Variances[12,1])*ZZt_iop1+
              as.numeric(Variances[13,1])*ZZt_iop2 + diag(as.numeric(Variances[14,1]),N)
          
        }else if(m%in%c(10)){
          W = as.numeric(Variances[10,1])*ZZt_subj + as.numeric(Variances[12,1])*ZZt_iop+
              as.numeric(Variances[13,1])*ZZt_iop1 + as.numeric(Variances[14,1])*ZZt_iop2 + diag(as.numeric(Variances[15,1]),N)
        }
        
        if(m!=1){DF.age = DF_compute2(B,Omega,W,as.numeric(Variances[2,1]))}
        BIC.all[k,m] = BIC.tmp + (DF.age-3)*log(n)  
        AIC.all[k,m] = AIC.tmp + (DF.age-3)*2 
      }	
    }
    if(k%%100==0){cat("basis",k,"done","\n")}
  }
  
  wgt=apply(basis_input,2,function(x)sum(x^2))
  wgt=wgt/sum(wgt)
  BIC.min.ind = matrix(0,k,n.model)
  AIC.min.ind = matrix(0,k,n.model)
  BIC.min.idx = apply(BIC.all,1,which.min)
  AIC.min.idx = apply(AIC.all,1,which.min)
  for(j in 1:dim(basis_input)[2]){
    BIC.min.ind[j,BIC.min.idx[j]]=1
    AIC.min.ind[j,AIC.min.idx[j]]=1
  }
  
  BIC.summary=apply(BIC.min.ind*wgt,2,sum)
  AIC.summary=apply(AIC.min.ind*wgt,2,sum)
  opt_idx = which.max(BIC.summary)
  opt_idx_aic = which.max(AIC.summary)
  
  ## Print Results
  output <- data.frame("Subject Level"=rep(c("No","Intercept"),each=5),
                       "Eye Level"=rep(c("No","Intercept","Linear IOP",
                                     "Hyperbolic IOP",
                                     "Intercept + Hyperbolic IOP"),2),
                       "BIC Score"=sprintf("%.4f",BIC.summary),
                       "AIC Score"=sprintf("%.4f",AIC.summary))
  
  sink("model_selection_output.txt", append = TRUE)
  cat("Step 3. Selection of Random Effects","\n","\n")
  print(output, row.names = FALSE)
  cat("\n", "Best Random effect model based on BIC:","\n","\n")
  print(output[opt_idx,], row.names = FALSE)
  cat("\n")
  sink()
  
  return(list(model=opt_idx,score=BIC.summary,model_aic=opt_idx_aic,score_aic=AIC.summary))
}
random_effect_result=random_effects_selection(D_wave,fixed_model=main_effect_result$model)  


##########################################################################
#### Save design matrices for evaluation of funtion and its derivative
#### These are for plotting nonparametric fit and its derivative
##########################################################################

Zg=Zres[(N+1):(N+ngrid),]
Zdg=Zres[(N+ngrid+1):dim(Zres)[1],]
write.table(Zg,"Zg_spline.txt",row.names=F,col.names=F)
write.table(Zdg,"Zdg_spline.txt",row.names=F,col.names=F)

##########################################################################
#### This tables are for computing DF for nonparametric effect
#### Spline basis evaluation and penalty matrix for computing DF
##########################################################################

B=SplineResult$B
Omega=SplineResult$Omega
write.table(B,"B.txt",row.names=F,col.names=F)
write.table(Omega,"Omega.txt",row.names=F,col.names=F)


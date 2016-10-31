

#*****************************************
#R_code_positive_and_negativeSTAN2
#By Daniel Berman
# This code uses hierarchical Bayesian parameter estimation to find parameters
# alpha, beta, delta, gamma, lambda, and luce, described by Cumulative Prospect Theory,
# developed by Amos Tversky and Daniel Kahneman. It uses Stan, developed by Gelman et al. 
# 
#In order to run the code, you need to input customize the following variables to suit your analysis purposes:
#   setwd
#   subjectsperrun
#   totalsubjects
#   subsetOfSubjects
#   location of data files. Input names in "loaded data" section
#   location you wish to save data. Input names in "Saving Data Values" section
#
#
#
# If testing against known values to verify the parameters for an MCMC (chains, interations, thins, and warm ups)
#   result in accurate estimates by producing a reduced chi-squared, make sure to input file names in 
#   "Known Parameter Values" section. You can run this code without doing this, but comment this out first.
#
#Note to user: This code often produces warnings related to values generated outside of the bounds specified
#               in the stan file. This is nothing to worry about.
#*****************************************
# The first step is to set the working directory, e.g.: 
setwd("C:\\Users\\daniel\\Documents\\R")


# Load the R2WinBUGS-package and define where WinBugs is installed
library(Rcpp)
library(inline)
library(rstan)
library(MASS)
library(ggmcmc)
library(coda)
set_cppo("fast")  # for best running speed


#cpt_model<-stan_model(file="cpt_hierarchical_model_positive_and_negative26.stan")
#save(cpt_model, file="cpt_model.RData")
translatedModel<-stanc(file="CPT_hierarchical_model_pos.stan",model_name='CPTbehavioral')
compiledModel<-stan_model(stanc_ret=translatedModel,verbose=FALSE)



#*******************************Loaded Data**************************************
#********************************************************************************
# Load information about the gamble-pairs (a gamble-pair is two opetions each with two possible outcomes occuring with a probability p and p-1).
# The files containing them should have the following structure:
# value of outcome 1 (column 1), probability of outcome 1 (column 2), value of outcome 2 # (column 3), probability of outcome 2 (column 4) 
# (gambles in rows).
# Add the 2nd subject's data in the 5th through 8th columns and the 3rd subject's in the 9th through 12th columns, and so forth. If the number of trials is not equal,
# put NaN until the total number of rows is equal.
prospects_b <- as.matrix(read.table("modeldatagamble18.txt"), fill=T) #gamble option (a gamble choice is a 1)
prospects_a <- as.matrix(read.table("modeldatasure18.txt"), fill=T) #sure option (a sure choice is a 0)

#Since RStan cannot handle NaNs, it converts NaNs into 2, which is an arbitrary value. It doesn't matter what the value because it is a placeholder, and those points aren't used in subsequent analysis.
prospects_b[as.matrix(is.na(prospects_b))]<-2
prospects_a[as.matrix(is.na(prospects_a))]<-2

# Load data (choice made by the first participant when presented the second gamble-pair is saved in
# column 1 row 2). The second subject's choices are in the 2nd column, the 3rd in the 3rd column.
# items_n tells you how many trials each subject had.
rawdata <- as.matrix(read.table("modeldatachoice18.txt"), fill=T)
items_n <- as.matrix(read.table("realpersonntrial18.txt"))

#Since RStan cannot handle NaNs, it converts NaNs into 2, which is an arbitrary value. It doesn't matter the actual value because those points aren't used in analysis.
rawdata[as.matrix(is.na(rawdata))]<-2


#*******************************Required Input Variables*************************
#********************************************************************************
totalsubjects<-300 #total number of subjects stored in the data sets
subsetOfSubjects<-300 #The total number of subjects being run in the analysis (a subset of totalsubjects). Equal to subjectsperrun*(number of loops for jjj). 
subjectsperrun=1 #number of subjects being run at the same time
numberOfInverseNormalDistributions=0
unnecessaryRows=1
numberofrounds<-totalsubjects/subjectsperrun #number of times you'll have to go through to run all subjects

charactersOfSubjects<-as.character(1:subjectsperrun) #string of 1 through number of subjects for the purpose of searching and selecting out from the covariance matrix
indicies<-paste("[",charactersOfSubjects,"]",sep="") #string of search terms for the covariance matrix, written as [number]


#creates a matrix with the first and last subject number run in each round
testcasenumber.n<-matrix(1,2,numberofrounds) #Intializing with "1", then altering in the next two steps.
testcasenumber.n[1,]<- matrix(seq(from=1, to=totalsubjects, by=subjectsperrun),1,)
testcasenumber.n[2,]<- matrix(seq(from=subjectsperrun, to=totalsubjects, by=subjectsperrun),1,)

#creates an empty 1x7 matrix. Used to store mean, sd, 2.5%, 25%, 50%, 75%, 97.5% bounds
alpha<-matrix(c(0),1,7)
gamma<-matrix(c(0),1,7)
luce<-matrix(c(0),1,7)

#The covariance matrix produced by RStan includes a number of categories that aren't needed in further analysis as a result of having to initialize certain values. 
#The covariance terms for parameters of interest are in a square submatrix between lowercovarmatrix and uppercovarmattrix of covarmatrix
lowercovarmatrix<-unnecessaryRows+numberOfInverseNormalDistributions*subjectsperrun #beginning of covariance terms of interest. The 13 and the 6 have to do with the number of inits declared in the "Beginning Bayesian Analysis" section.
uppercovarmatrix<-lowercovarmatrix+3*subjectsperrun-1 #end of covariance terms of interest
covariancematrix<-matrix(c(0),1,3) #initializes an empty 1x6 matrix used to stored covariance terms


listOfSubjects<-as.character(1:totalsubjects) #turns a list of numbers 1:total subjects into characters
subjectIdentifiers<-paste("Subject",listOfSubjects,sep="") #Adds Subject before the number in listOfSubjects
subjectData<-vector("list", length(subjectIdentifiers))#Creates a vector of lists. This becomes a list of covariance matricies
names(subjectData)<-subjectIdentifiers #names each list after each subject (eg Subject1, Subject 2, ...)
alpha_upper=matrix(1,totalsubjects,1)
alpha_lower=matrix(0,totalsubjects,1)

gamma_upper=matrix(1,totalsubjects,1)
gamma_lower=matrix(0,totalsubjects,1)

luce_upper=matrix(1.5,totalsubjects,1)
luce_lower=matrix(0.5,totalsubjects,1)


#realValues = matrix(nrow=3, ncol=subsetOfSubjects)
realValues = matrix(nrow=3, ncol=subsetOfSubjects)
realValues[1,]=t(as.matrix(read.table("alpha.txt"), fill=T))[1:subsetOfSubjects]
realValues[2,]=t(as.matrix(read.table("gamma.txt"), fill=T))[1:subsetOfSubjects]
realValues[3,]=t(as.matrix(read.table("luce.txt"), fill=T))[1:subsetOfSubjects]

#*******************************Beginning Bayesian Analysis**********************
#********************************************************************************

for(jjj in  1:subsetOfSubjects){
  #for loop that runs through the specified subjects.
  
  
  prospects_a_sub=prospects_a[,(4*jjj-3):(4*jjj)]
  prospects_b_sub=prospects_b[,(4*jjj-3):(4*jjj)]
  items_n_pos=items_n[1,jjj]
  rawdata_sub=rawdata[1:items_n_pos,jjj]
  repeatsize=1
  prospects_a_sub.pos=prospects_a_sub[1:items_n_pos,]
  prospects_b_sub.pos=prospects_b_sub[1:items_n_pos,]

  
  ax=prospects_a_sub.pos[,1]^realValues[1,jjj]
  ay=prospects_a_sub.pos[,3]^realValues[1,jjj]
  wax=prospects_a_sub.pos[,2]^realValues[2,jjj]/(prospects_a_sub.pos[,2]^realValues[2,jjj]+(1-prospects_a_sub.pos[,2])^realValues[2,jjj])^(1/realValues[2,jjj])
  way=prospects_a_sub.pos[,4]^realValues[2,jjj]/(prospects_a_sub.pos[,4]^realValues[2,jjj]+(1-prospects_a_sub.pos[,4])^realValues[2,jjj])^(1/realValues[2,jjj])
  vfa=ax*wax+ay*way
  
  bx=prospects_b_sub.pos[,1]^realValues[1,jjj]
  by=prospects_b_sub.pos[,3]^realValues[1,jjj]
  wbx=prospects_b_sub.pos[,2]^realValues[2,jjj]/(prospects_b_sub.pos[,2]^realValues[2,jjj]+(1-prospects_b_sub.pos[,2])^realValues[2,jjj])^(1/realValues[2,jjj])
  wby=prospects_b_sub.pos[,4]^realValues[2,jjj]/(prospects_b_sub.pos[,4]^realValues[2,jjj]+(1-prospects_b_sub.pos[,4])^realValues[2,jjj])^(1/realValues[2,jjj])
  vfb=bx*wbx+by*wby
  
  vf=1/(1+exp(realValues[3,jjj]*(vfb-vfa)))
  
  vfpos=c(vf)#,ay,bx,by)
  
#   axneg=-realValues[5,jjj]*abs(prospects_a_sub.neg[,1])^realValues[2,jjj]
#   ayneg=-realValues[5,jjj]*abs(prospects_a_sub.neg[,3])^realValues[2,jjj]
#   waxneg=prospects_a_sub.neg[,2]^realValues[3,jjj]/(prospects_a_sub.neg[,2]^realValues[3,jjj]+(1-prospects_a_sub.neg[,2])^realValues[3,jjj])^(1/realValues[3,jjj])
#   wayneg=prospects_a_sub.neg[,4]^realValues[3,jjj]/(prospects_a_sub.neg[,4]^realValues[3,jjj]+(1-prospects_a_sub.neg[,4])^realValues[3,jjj])^(1/realValues[3,jjj])
#   vfaneg=axneg*waxneg+ayneg+wayneg
#   
#   bxneg=-realValues[5,jjj]*abs(prospects_b_sub.neg[,1])^realValues[2,jjj]
#   byneg=-realValues[5,jjj]*abs(prospects_b_sub.neg[,3])^realValues[2,jjj]
#   wbxneg=prospects_b_sub.neg[,2]^realValues[3,jjj]/(prospects_b_sub.neg[,2]^realValues[3,jjj]+(1-prospects_b_sub.neg[,2])^realValues[3,jjj])^(1/realValues[3,jjj])
#   wbyneg=prospects_b_sub.neg[,4]^realValues[3,jjj]/(prospects_b_sub.neg[,4]^realValues[3,jjj]+(1-prospects_b_sub.neg[,4])^realValues[3,jjj])^(1/realValues[3,jjj])
#   vfbneg=bxneg*wbxneg+byneg+wbyneg
#   
#   vfneg=realValues[6,jjj]*(vfbneg-vfaneg)
#   
  #rawdata_sub=c(vfpos)
  N<-dim(prospects_a_sub)[1] #placeholder, dimension of prospects_a Used to set limits on input data in the Stan code
  
  testcasenumber<-testcasenumber.n[1:2,jjj] #the subjects that are going to be run in this iteration of the for loop
  
  subjectsbeingrun<-matrix(seq(to=testcasenumber.n[2,jjj], from=testcasenumber.n[1,jjj]),1,) #creates a matrix from testcasenumber[1,jjj] to testcasenumber[2,jjj] incremented by 1
  
  #creates initial values for parameters with the values and dimensions necessary
  inits = function()
  {
    
    list(alpha=c(.12),
         gamma=c(.1),
         luce=c(.5)
      )
  }
  
  # Define what information should be passed into Stan
  player_dat  = list("prospects_a_pos"=prospects_a_sub.pos,
                     "prospects_b_pos"=prospects_b_sub.pos,
                     "N"=N,
                     "rawdata"=rawdata_sub, 
                     "items_n_pos"=items_n_pos,
                     "alpha_upper"=alpha_upper[jjj],"alpha_lower"=alpha_lower[jjj],
                     "gamma_upper"=gamma_upper[jjj],"gamma_lower"=gamma_lower[jjj],
                     "luce_upper"=luce_upper[jjj], "luce_lower"=luce_lower[jjj]
                      ) 
  
  
  # Define the variables of interest. WinBugs will return these to R when the analysis is finished (and WinBugs is # closed).  
  parameters = c("alpha", "gamma", 'luce')
  
  
  
  # The following command calls RStan with specific options. For a detailed description, see 
  # http://mc-stan.org/rstan.html. 
  # Specific to the following code:
  # file: name of the file where model is writt5en in Stan
  # data: variables created in R that contains data for analysis
  # init: initial values specified for distributions specified in the model
  # chains: the number of chains generated
  # iter: the number of samples sampled in each chain
  # warmup: the number of burn-in samples
  # thin: WinBugs will return each x:th sample to R if n.thin is set to x
  
  hierarchical= sampling(compiledModel,data=player_dat, iter=1000, chains=3, seed=2)#,
  #init=inits,
  #chains=2, iter=20,warmup=1, thin=1)  
  
  #prints the estimated parameters from hierarchical, digits is the number of digits it prints
  print(hierarchical, digits=10)
  
  #mcmc_hierarchical <- mcmc.list(lapply(1:ncol(hierarchical), function(x) mcmc(as.array(hierarchical)[,x,])))
  
  #ggs_hierarchical <- ggs(mcmc_hierarchical)
  #ggs_traceplot(ggs_hierarchical)
  #ggs_density(ggs_hierarchical)
  #*******************************Generates Covariance Matricies*******************
  #********************************************************************************
  #produces the covariance matrix for all returned parameters. There are a number more parameters than those of interest,
  #so the following steps are meant to store only the ones of interest: alpha, beta, delta, gamma, lambda, and luce
  covarmatrix<-cov(as.matrix(hierarchical))
  covariancemat<-covarmatrix[lowercovarmatrix:uppercovarmatrix,lowercovarmatrix:uppercovarmatrix]
  parameterMeans<-colMeans(as.matrix(hierarchical))
  
  #The six parameters of interest, alpha, beta, delta, gamma, lambda, and luce, the following variables are described:                                                                    
  # *Mean - Determines the mean from *HatEstimates for each subject
  # *Sd - Determines the standard deviation from *HatEstimates for each subject
  # *Quantiles - Determines the Quantiles from *HatEstimates for each subject
  # *Temp - Stores all information in a single variable with the format: mean, standard devaition, 2.5%, 25%, 50%, 75%, 97.5%
  parametersTemp<-c('alpha', 'gamma','luce')
  orderOfParameters<-colnames(as.matrix(hierarchical))
  for (i in 1:length(parametersTemp)){
    parameterlocation=grep(parametersTemp[i],orderOfParameters)
    varName=paste(parametersTemp[i],'Temp',sep='')
    quantileTemp=quantile(as.matrix(hierarchical)[,parameterlocation])
    sdTemp=sd(as.matrix(hierarchical)[,parameterlocation])
    names(sdTemp)<-'sd'
    assign(varName,c(parameterMeans[parameterlocation],sdTemp,quantileTemp))
  }
  
  #The following loop is used to generate the covariance matrix and store it in a larger variable covariancematrix.
  #The first line stores the coordinates of a specific variable estimate's covariance matrix elements by:
  # as.matrix(row.names(covariancemat))-takes the relevant names and stores them in a matrix
  #   grep - searches and identifies matches of indicies[iii] within as.martix(row.names(covariancemat)) to find alpha[#], beta[#], etc.
  #     list - converts it to a list
  #       rep - doubles the number indices
  #         expand.grid - creates a list of all combinations of the two lists
  #Then pulls out the elements of the matrix from covariancemat
  #concatenates it with the previous iteration to store in one location
  
  for( iii in 1:subjectsperrun){
    labelsofCovarianceMat<-as.matrix(row.names(covariancemat))
    if (subjectsperrun>1){
      indiciesOfInterest<-grep(pattern=indicies[iii],labelsofCovarianceMat,fixed=TRUE)
      combinationOfInterest<-rep(list(indiciesOfInterest),2)
      coordinatesOfCovarMat<-expand.grid(combinationOfInterest)
      blockedcovarMatrix<-matrix(covariancemat[as.matrix(coordinatesOfCovarMat)],nrow=6)
    }
    if (subjectsperrun==1){
      blockedcovarMatrix=covariancemat
    }
    covariancematrix<-rbind(covariancematrix,blockedcovarMatrix)
  }
  
  #Stores all six variable estimates and basic statistics about them
  alpha<-rbind(alpha,alphaTemp)
  gamma<-rbind(gamma,gammaTemp)
  luce<-rbind(luce,luceTemp)

}

#Creates names for the columns in alpha, beta, delta, gamma, lambda, and luce



#*******************************Chi-Square Goodness of Fit***********************
#********************************************************************************
#This section of the code is for testing the goodness of fit against known parameters if using simulated data
#

#Total number of subjects that were run in the experiment
subsetOfSubjects

charactersOfSubjectsRun<-as.character(1:subsetOfSubjects) #turns a list of numbers describing first subject:last subject into characters
subjectIdentifiersRun<-paste("Subject",charactersOfSubjectsRun,sep="")#Adds Subject before the number in listOfSubjects
invCovarianceMatrix<-vector("list", length(subjectIdentifiersRun))#Creates a vector of lists
invVarianceMatrix<-vector("list", length(subjectIdentifiersRun))#Creates a vector of lists
names(subjectData)<-subjectIdentifiersRun #names each list after each subject (eg Subject1, Subject 2, ...)

#This creates a matrix of the means of estimates. Has the following format:
# .....alpha....
# .....beta.....
# .....delta....
# .....gamma....
# .....lambda...
# .....luce.....
estimates<-matrix(nrow=3, ncol=subsetOfSubjects)
estimates[1,]<-alpha[2:(subsetOfSubjects+1),1]
estimates[2,]<-gamma[2:(subsetOfSubjects+1),1]
estimates[3,]<-luce[2:(subsetOfSubjects+1),1]

colnames(estimates)<-subjectIdentifiers[1:subsetOfSubjects]
rownames(estimates)<-labelsofCovarianceMat

#Creates an the inverse matrix for covariance matrix
for (j in 1:subsetOfSubjects){
  invCovarianceMatrix[[j]]<-solve(covariancematrix[(2+(3*(j-1))):(1+(3*j)),1:3])
  colnames(invCovarianceMatrix[[j]])<-labelsofCovarianceMat
  rownames(invCovarianceMatrix[[j]])<-labelsofCovarianceMat
  invVarianceMatrix[[j]]<-solve(diag(3)*covariancematrix[(2+(3*(j-1))):(1+(3*j)),1:3])
  colnames(invVarianceMatrix[[j]])<-labelsofCovarianceMat
  rownames(invVarianceMatrix[[j]])<-labelsofCovarianceMat
}

#*******************************Known Parameter Values***********************
#Uploads the known values of the parameters used to generate simulated data. User must change file names to suit actual data.
# realValues = matrix(nrow=3, ncol=subsetOfSubjects)
# realValues[1,]=t(as.matrix(read.table("alpha.txt"), fill=T))[1:subsetOfSubjects]
# realValues[2,]=t(as.matrix(read.table("gamma.txt"), fill=T))[1:subsetOfSubjects]
rownames(realValues)<-labelsofCovarianceMat
colnames(realValues)<-subjectIdentifiers[1:subsetOfSubjects]

estimateDifference=realValues-estimates
chiSquareElement<-matrix(nrow=1,ncol=subsetOfSubjects)
#chiSquareElement_Variance<-matrix(nrow=1,ncol=subsetOfSubjects)
chiSquareElement_Variance<-vector("list", length(subjectIdentifiersRun))
for (j in 1:subsetOfSubjects){
  chiSquareElement[,j]<-estimateDifference[,j]%*%invCovarianceMatrix[[j]]%*%estimateDifference[,j]
  chiSquareElement_Variance[[j]]<-estimateDifference[,j]*invVarianceMatrix[[j]]*estimateDifference[,j]
}
reducedChiSquare<-sum(chiSquareElement)/(3*subsetOfSubjects)
reducedChiSquare
chiSquareMatrix_Variance<-matrix(Reduce('+', chiSquareElement_Variance),3,3)
colnames(chiSquareMatrix_Variance)<-labelsofCovarianceMat
rownames(chiSquareMatrix_Variance)<-labelsofCovarianceMat
reducedChiSquare_Variance<-chiSquareMatrix_Variance%*%matrix(1,3,1)/(subsetOfSubjects-1)
print(reducedChiSquare_Variance)


mcmc_hierarchical <- mcmc.list(lapply(1:ncol(hierarchical), function(x) mcmc(as.array(hierarchical)[,x,])))

ggs_hierarchical <- ggs(mcmc_hierarchical)
ggs_traceplot(ggs_hierarchical)
ggs_density(ggs_hierarchical)
ggs_compare_partial(ggs_hierarchical)
ggs_running(ggs_hierarchical)
ggs_autocorrelation(ggs_hierarchical)
ggs_crosscorrelation(ggs_hierarchical)
ggs_caterpillar(ggs_hierarchical)

#*******************************Saving Data Values***********************
#These lines write the covariancematrix, subjectData, alpha, beta, delta, gamma, lambda, and luce statistics values to text files for future reference.

write.table(covariancematrix,"C:\\Users\\Daniel\\Documents\\R\\Rstan_covariancemat_dataset18_pos.xls",sep="\t")
write.table(invCovarianceMatrix,"C:\\Users\\Daniel\\Documents\\R\\Rstan_invcovariancemat_dataset18_pos.xls",sep="\t")
write.table(alpha,"C:\\Users\\Daniel\\Documents\\R\\Rstan_alphaestimate_dataset18_pos.xls",sep="\t")
write.table(gamma,"C:\\Users\\Daniel\\Documents\\R\\Rstan_gammaestimate_dataset18_pos.xls",sep="\t")
write.table(luce,"C:\\Users\\Daniel\\Documents\\R\\Rstan_luceestimate_dataset18_pos.xls",sep="\t")
allChiSquares<-matrix(c(chiSquareMatrix_Variance%*%matrix(1,3,1), sum(chiSquareElement)),4,1)
degreesOfFreedom<-rbind(matrix((subsetOfSubjects-1),3,1), 3*subsetOfSubjects)
allChiSquares<-cbind(allChiSquares,degreesOfFreedom)
rownames(allChiSquares)<-c(labelsofCovarianceMat,'covariance matrix')
colnames(allChiSquares)<-c('Chi-Square', 'd.f.')
write.table(allChiSquares,"C:\\Users\\Daniel\\Documents\\R\\Rstan_AllChiSquaresValues_dataset18_pos.xls",sep="\t")


#*****************************************
#R_code_positive_and_negativeSTAN2
#By Daniel Berman
# This code uses hierarchical Bayesian parameter estimation to find parameters
# alpha, beta, delta, gamma, lambda, and luce, described by Cumulative Prospect Theory,
# developed by Amos Tversky and Daniel Kahneman. It uses Stan, developed by Gelman et al. 
# 
#In order to run the code, you need to input customize the following variables to suit your analysis purposes:
#   testingKnownData-for checking the chisquare. If you are, set to 1, if not, set to 0
#
#   workingDirectory - the working directory where the stan file and the data is located
# 
#   nameOfFilesSave - a basic name describing the data that was run
#
#   firstSubject- the index (in rawdata) of the first subject being run
#
#   totalsubjects - total number of subjects
#
#   subsetOfSubjects - the position of the last subject being run in the data files.
#  
#   prospects_b - option B. There should be four columns for each subject in the data file. Columns are tab separated.
#                 The file should have the following structure:
#                 Column 1 - value of outcome 1 for the first subject
#                 Column 2 - probability of outcome 1  for the first subject
#                 Column 3 - value of outcome 2 for the first subject
#                 Column 4 - probability of outcome 2 for the first subject
#                 Repeat for n subjects. 
#                 If this is a sure option, then the values entered for columns 3 and 4 should be 0.
#                 Example: An option will return 100 10% of the time, and 50 90% of the time.
#                          In the data it will appear as 100  .1  50  .9
#                 Add the 2nd subject's data in the 5th through 8th columns and the 3rd subject's in the 9th through 12th columns, and so forth. 
#                 If the number of trials is not equal across subjects, enter NaN in the text file until the total number of rows is equal.
#
#   prospects_a - See prospects_b
#   
#   rawdata - contains the subjects choices. A 1 indicates prospects_b was chosen, a 0 indicates prospects_a was chosen. Each subject's data is
#             in a column. If the number of trials is not equal across subjects, enter NaN in the text file until the total number of rows is equal.
#
#   items_n - row of numbers indicating how many data points exist for each subject, excluding NaNs.
#             If subjects 1, 2, 3, and 4 had 100, 130, 124, and 80 data points respectively, the text file should contain tab separated numbers:
#             100 130 124 80
#
#   knownAlpha, knownBeta, knownGamma, knownDelta, knownLambda, knownLuce - These are the known values of the parameters alpha, gamma, and luce where alpha=(0,1), beta=(0,1), gamma=(0,1), delta=(0,1), Lambda=(1,5), luce=(.5,1.5)
#                Six text files must be listed in these variables
#
#   seed - user can set the value of the seed or remove it completely in the sampling function
#
# If testing against known values to verify the parameters for an MCMC (chains, interations, thins, and warm ups)
#   result in accurate estimates by producing a reduced chi-squared, make sure to input file names in 
#   "Known Parameter Values" section. You can run this code without doing this, but comment this out first.
#
#Note to user: This code often produces warnings related to values generated outside of the bounds specified
#               in the stan file. This is nothing to worry about.
#*****************************************


# Load the rstan and other package
library(Rcpp)
library(inline)
library(rstan)
library(MASS)
library(ggmcmc)
library(coda)
set_cppo("fast")  # for best running speed


#*******************************Loaded Data**************************************
#********************************************************************************
# The first step is to set the working directory, e.g.: 
testingKnownData=1
workingDirectory="C:\\Users\\bermads1\\Documents\\Daniel\\bayesian\\FINAL BAYESIAN CODE 20150522"
setwd(workingDirectory)

nameOfFilesSave='dataset16_pos'

#Sets the model you are working with and translates and compiles it once. 
translatedModel<-stanc(file="CPT_hierarchical_model_all.stan",model_name='CPTbehavioral')
compiledModel<-stan_model(stanc_ret=translatedModel,verbose=FALSE)

#Subject gamble data
prospects_b <- as.matrix(read.table("modeldatagamble16.txt"), fill=T) #gamble option (a gamble choice is a 1)
prospects_a <- as.matrix(read.table("modeldatasure16.txt"), fill=T) #sure option (a sure choice is a 0)

#Subject choice data
rawdata <- as.matrix(read.table("modeldatachoice16.txt"), fill=T)

#This needs to be a table with the dimensions 3 x totalsubjects
items_n <- as.matrix(read.table("realpersonntrial16.txt"))

#Since RStan cannot handle NaNs, it converts NaNs into 2, which is an arbitrary value. It doesn't matter what the value because it is a placeholder, and those points aren't used in subsequent analysis.
#The items_n should prevent this from even happening. If it doesn't, you likely entered the wrong number in items_n as you'll get an error regarding the bernoulli_logit needing a 0 or 1
prospects_b[as.matrix(is.na(prospects_b))]<-2
prospects_a[as.matrix(is.na(prospects_a))]<-2
rawdata[as.matrix(is.na(rawdata))]<-2
if (testingKnownData==1){
  knownAlpha<-"alpha.txt"
  knownBeta<-"beta.txt"
  knownGamma<-'gamma.txt'
  knownDelta<-"delta.txt"
  knownLambda<-"lambda.txt"
  knownLuce<-'luce.txt'
}
#*******************************Required Input Variables*************************
#********************************************************************************
firstSubject<-2 #The index number of the first subject being run.
subsetOfSubjects<-3 #The total number of subjects being run in the analysis (must be less than or equal to totalsubjects). 
totalsubjects<-subsetOfSubjects-firstSubject+1 #total number of subjects stored in the data sets
numberOfParameters=6 #number of parameters being estimated (3, alpha, gamma, luce)

#These are the known values of the parameters

charactersOfSubjects<-as.character(firstSubject:subsetOfSubjects) #string of 1 through number of subjects for the purpose of searching and selecting out from the covariance matrix
indicies<-paste("[",charactersOfSubjects,"]",sep="") #string of search terms for the covariance matrix, written as [number]

#creates an empty 1x7 matrix. Used to store mean, sd, 0%, 25%, 50%, 75%, 100% bounds
alpha<-matrix(c(0),1,7)
beta<-matrix(c(0),1,7)
gamma<-matrix(c(0),1,7)
delta<-matrix(c(0),1,7)
lambda<-matrix(c(0),1,7)
luce<-matrix(c(0),1,7)

#The covariance matrix produced by RStan includes a number of categories that aren't needed in further analysis as a result of having to initialize certain values. 
#The covariance terms for parameters of interest are in a square submatrix between lowercovarmatrix and uppercovarmattrix of covarmatrix
lowercovarmatrix<-1 #beginning of covariance terms of interest. 
uppercovarmatrix<-numberOfParameters #end of covariance terms of interest
covariancematrix<-matrix(c(0),1,numberOfParameters) #initializes an empty 1x6 matrix used to stored covariance terms


listOfSubjects<-as.character(firstSubject:subsetOfSubjects) #turns a list of numbers 1:total subjects into characters
subjectIdentifiers<-paste("Subject",listOfSubjects,sep="") #Adds Subject before the number in listOfSubjects

#*******************************Beginning Bayesian Analysis**********************
#********************************************************************************

for(jjj in  firstSubject:subsetOfSubjects){
  #for loop that runs through the specified subjects one at a time
  
  #Takes the jjjth subject data
  prospects_a_sub=prospects_a[,(4*jjj-3):(4*jjj)]
  prospects_b_sub=prospects_b[,(4*jjj-3):(4*jjj)]
  prospects_b_sub=prospects_b[,(4*jjj-3):(4*jjj)]
  
#  items_n_mix=items_n[1,jjj]
   items_n_mix=items_n[3,jjj]
  
  rawdata_sub=rawdata[(1):items_n_mix,jjj]
  prospects_a_sub.all=prospects_a_sub[(1):items_n_mix,]
  prospects_b_sub.all=prospects_b_sub[(1):items_n_mix,]
    
  # Define what information should be passed into Stan
  player_dat  = list("prospects_a_mix"=prospects_a_sub.all,
                     "prospects_b_mix"=prospects_b_sub.all,
                     "rawdata"=rawdata_sub, 
                     "items_n_mix"=items_n_mix) 
  
  # The following command calls RStan with specific options. For a detailed description, see 
  # http://mc-stan.org/rstan.html. 
  # Specific to the following code:
  # compiledModel is the model compiled earlier
  # data: variables created in R that contains data for analysis
  # chains: the number of chains generated
  # iter: the number of samples sampled in each chain
  # warmup and thin are set to the default values
  
  hierarchical= sampling(compiledModel,data=player_dat, iter=1000, chains=3, seed=2)#,
  
  #prints the estimated parameters from hierarchical, digits is the number of digits it prints
  print(hierarchical, digits=10)
  
  #converts to a mcmc type for the purpose of graphic. This can be disabled if desired.
  mcmc_hierarchical <- mcmc.list(lapply(1:ncol(hierarchical), function(x) mcmc(as.array(hierarchical)[,x,])))
  ggs_hierarchical <- ggs(mcmc_hierarchical)
  
  ggs_traceplot(ggs_hierarchical) #traceplot of the parameter
  ggs_density(ggs_hierarchical) #density plot of the parameter 
  ggs_autocorrelation(ggs_hierarchical) #autocorrelation plots
  ggs_crosscorrelation(ggs_hierarchical) #Grid visualizing covariance matrix
  #   ggs_compare_partial(ggs_hierarchical)
  #   ggs_running(ggs_hierarchical)
  #   

  #********************************************************************************
  #*******************************Generates Covariance Matricies*******************
  #********************************************************************************
  #produces the covariance matrix for all returned parameters. There are a number more parameters than those of interest,
  #so the following steps are meant to store only the ones of interest: alpha, gamma, and luce
  covarmatrix<-cov(as.matrix(hierarchical))
  covariancemat<-covarmatrix[lowercovarmatrix:uppercovarmatrix,lowercovarmatrix:uppercovarmatrix]
  parameterMeans<-colMeans(as.matrix(hierarchical))#Calculates the mean parameter values (this matches the bayesian one)
  
  #The six parameters of interest, alpha, beta, delta, gamma, lambda, and luce, the following variables are described:                                                                    
  # *Mean - Determines the mean from *HatEstimates for each subject
  # *Sd - Determines the standard deviation from *HatEstimates for each subject
  # *Quantiles - Determines the Quantiles from *HatEstimates for each subject
  # *Temp - Stores all information in a single variable with the format: mean, standard devaition, 0%, 25%, 50%, 75%, 100%
  parametersTemp<-c('alpha','beta', 'gamma', 'delta', 'lambda','luce')
 
  orderOfParameters<-colnames(as.matrix(hierarchical)) #looks at the order of parameters in hierarchical in case it's not the same as the one listed above
  for (i in 1:length(parametersTemp)){
    parameterlocation=grep(parametersTemp[i],orderOfParameters) #finds the location of the parameter names
    varName=paste(parametersTemp[i],'Temp',sep='') #concatenates the variable name with with Temp
    quantileTemp=quantile(as.matrix(hierarchical)[,parameterlocation]) #determines the quantiles for that parameter
    sdTemp=sd(as.matrix(hierarchical)[,parameterlocation]) #determines the standard deviation for that parameter
    names(sdTemp)<-'sd'
    assign(varName,c(parameterMeans[parameterlocation],sdTemp,quantileTemp)) #combines the mean, sd, and 5 quantile values
  }
  
  labelsofCovarianceMat<-as.matrix(row.names(covariancemat))#takes the labels from the covariance matrix produced by the estimates of sampling
  blockedcovarMatrix=covariancemat #a temporary matrix to stored the jjjth subject's covariance matrix
  covariancematrix<-rbind(covariancematrix,blockedcovarMatrix) #concatenates with the rest of the covariance matricies


  #Stores all six variable estimates and basic statistics about them
  alpha<-rbind(alpha,alphaTemp)
  beta<-rbind(beta,betaTemp)
  gamma<-rbind(gamma,gammaTemp)
  delta<-rbind(delta,deltaTemp)
  lambda<-rbind(lambda,lambdaTemp)
  luce<-rbind(luce,luceTemp)
  
}


invCovarianceMatrix<-vector("list", length(subjectIdentifiers))#Creates a vector of lists
invVarianceMatrix<-vector("list", length(subjectIdentifiers))#Creates a vector of lists

#This creates a matrix of the means of estimates. Has the following format:
# .....alpha....
# .....beta.....
# .....delta....
# .....gamma....
# .....lambda...
# .....luce.....
estimates<-matrix(nrow=numberOfParameters, ncol=totalsubjects)
estimates[1,]<-alpha[2:(totalsubjects+1),1]
estimates[2,]<-beta[2:(totalsubjects+1),1]
estimates[3,]<-gamma[2:(totalsubjects+1),1]
estimates[4,]<-delta[2:(totalsubjects+1),1]
estimates[5,]<-lambda[2:(totalsubjects+1),1]
estimates[6,]<-luce[2:(totalsubjects+1),1]

colnames(estimates)<-subjectIdentifiers[1:totalsubjects]
rownames(estimates)<-labelsofCovarianceMat

#Creates an the inverse matrix for covariance matrix and for the variance matrix
for (j in 1:totalsubjects){
  #covariance matrix inverse
  invCovarianceMatrix[[j]]<-solve(covariancematrix[(2+(numberOfParameters*(j-1))):(1+(numberOfParameters*j)),1:numberOfParameters])
  colnames(invCovarianceMatrix[[j]])<-labelsofCovarianceMat
  rownames(invCovarianceMatrix[[j]])<-labelsofCovarianceMat
  #variance matrix inverse
  invVarianceMatrix[[j]]<-solve(diag(numberOfParameters)*covariancematrix[(2+(numberOfParameters*(j-1))):(1+(numberOfParameters*j)),1:numberOfParameters])
  colnames(invVarianceMatrix[[j]])<-labelsofCovarianceMat
  rownames(invVarianceMatrix[[j]])<-labelsofCovarianceMat
}

covariancematrix=covariancematrix[-1,]
alpha=alpha[-1,]
beta=beta[-1,]
gamma=gamma[-1,]
delta=delta[-1,]
lambda=lambda[-1,]
luce=luce[-1,]

timeStamp<-gsub("[[:punct:]]","",Sys.time())#timestamp when files are written to prevent overwriting of files

write.table(covariancematrix,paste(workingDirectory,'Rstan_covariancemat_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(invCovarianceMatrix,paste(workingDirectory,'Rstan_invcovariancemat_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(alpha,paste(workingDirectory,'Rstan_alphaestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(beta,paste(workingDirectory,'Rstan_betaestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(gamma,paste(workingDirectory,'Rstan_gammaestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(delta,paste(workingDirectory,'Rstan_deltaestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(lambda,paste(workingDirectory,'Rstan_lambdaestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
write.table(luce,paste(workingDirectory,'Rstan_luceestimate_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")


#*******************************Chi-Square Goodness of Fit***********************
#********************************************************************************
#This section of the code is for testing the goodness of fit against known parameters if using simulated data
#

#Total number of subjects that were run in the experiment
if (testingKnownData==1){
  realValues = matrix(nrow=numberOfParameters, ncol=totalsubjects)
  realValues[1,]=t(as.matrix(read.table(knownAlpha), fill=T))[firstSubject:subsetOfSubjects]
  realValues[2,]=t(as.matrix(read.table(knownBeta), fill=T))[firstSubject:subsetOfSubjects]
  realValues[3,]=t(as.matrix(read.table(knownGamma), fill=T))[firstSubject:subsetOfSubjects]
  realValues[4,]=t(as.matrix(read.table(knownDelta), fill=T))[firstSubject:subsetOfSubjects]
  realValues[5,]=t(as.matrix(read.table(knownLambda), fill=T))[firstSubject:subsetOfSubjects]
  realValues[6,]=t(as.matrix(read.table(knownLuce), fill=T))[firstSubject:subsetOfSubjects]
  
  rownames(realValues)<-labelsofCovarianceMat
  colnames(realValues)<-subjectIdentifiers[1:totalsubjects]
  
  estimateDifference=realValues-estimates
  
  chiSquareElement<-matrix(nrow=1,ncol=totalsubjects)#storage variable for each element of the chisquare using the covariance matrix
  chiSquareElement_Variance<-vector("list", length(subjectIdentifiers)) #storage variable for each element of teh chisquare using the variance matrix
  
  #Calculates the chisquare values for each subject estimate using the formulat:
  # transpose(theta-theta_hat)*(CovarianceMatrix)^-1*(theta-theta_hat)
  for (j in 1:totalsubjects){
    chiSquareElement[,j]<-estimateDifference[,j]%*%invCovarianceMatrix[[j]]%*%estimateDifference[,j]
    chiSquareElement_Variance[[j]]<-estimateDifference[,j]*invVarianceMatrix[[j]]*estimateDifference[,j]
  }
  
  #sums over all chisquare values and calculates the reduced chisquare, where df=numberOfParameters*totalsubjects
  reducedChiSquare<-sum(chiSquareElement)/(numberOfParameters*totalsubjects)
  reducedChiSquare
  
  #sums over all chisquare values and calculates the reduced chisquare, where df=(totalsubjects-1)
  chiSquareMatrix_Variance<-matrix(Reduce('+', chiSquareElement_Variance),numberOfParameters,numberOfParameters)
  colnames(chiSquareMatrix_Variance)<-labelsofCovarianceMat
  rownames(chiSquareMatrix_Variance)<-labelsofCovarianceMat
  reducedChiSquare_Variance<-chiSquareMatrix_Variance%*%matrix(1,numberOfParameters,1)/(totalsubjects-1)
  print(reducedChiSquare_Variance)
  
  
  allChiSquares<-matrix(c(chiSquareMatrix_Variance%*%matrix(1,numberOfParameters,1), sum(chiSquareElement)),(numberOfParameters+1),1)
  degreesOfFreedom<-rbind(matrix((totalsubjects-1),numberOfParameters,1), numberOfParameters*totalsubjects)
  allChiSquares<-cbind(allChiSquares,degreesOfFreedom)
  rownames(allChiSquares)<-c(labelsofCovarianceMat,'covariance matrix')
  colnames(allChiSquares)<-c('Chi-Square', 'd.f.')
  
  #*******************************Saving Data Values***********************
  #These lines write the covariancematrix, subjectData, alpha, beta, delta, gamma, lambda, and luce statistics values to text files for future reference.
  write.table(allChiSquares,paste(workingDirectory,'Rstan_AllChiSquaresValues_',nameOfFilesSave,'_',timeStamp,'.xls'),sep="\t")
}
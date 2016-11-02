# CumulativeProspectTheory
The files uploaded run Bayesian Parameter Estimation for Cumulative Prospect Theory (CPT). The R code loads in the data and runs RStan to execute the Bayesian parameter estimation. 

There are six files. Two files for data with only positive valued gamble options. Two files for data with only negative valued gample options. Two files for data with positive and negative valued gamble options. According to CPT, there are two sets of parameters: one for positive values and one for negative values. Therefore, separating the estimation into positive, negative, and mixed is necessary, as the absence of positive or negative values means there is no data to estimate parameters.

## How to use the code
There are two ways to use the code. The first is running it on data to determine the parameters for CPT. The other is to validate it. This is the same for all three R files.

### Experimental Use
 Set testingKnownData=0
 
 Set the workingDirectory
 
 Set the nameOfFilesSave to say where to save the data
 
 firstSubject- the index (in rawdata) of the first subject being run. If you are running all the subjects, just set to 1.
 
 subsetOfSubjects - the position of the last subject being run in the data files.
 
 prospects_b - option B. There should be four columns for each subject in the data file. Columns are tab separated.
 
     The file should have the following structure:
     
     Column 1 - value of outcome 1 for the first subject
     
     Column 2 - probability of outcome 1  for the first subject
     
     Column 3 - value of outcome 2 for the first subject
     
     Column 4 - probability of outcome 2 for the first subject
     
     Repeat for n subjects. 
     
     If this is a sure option, then the values entered for columns 3 and 4 should be 0.
     Example: An option will return 100 10% of the time, and 50 90% of the time. In the data it will appear as 100  .1  50  .9
     
     Add the 2nd subject's data in the 5th through 8th columns and the 3rd subject's in the 9th through 12th columns, and so forth. If the number of trials is not equal across subjects, enter NaN in the text file until the total number of rows is equal.
 
 prospects_a - See prospects_b
 
 rawdata - contains the subjects choices. A 1 indicates prospects_b was chosen, a 0 indicates prospects_a was chosen. Each subject's data is
 in a column. If the number of trials is not equal across subjects, enter NaN in the text file until the total number of rows is equal.
 
 items_n - 3xn row of numbers indicating how many data points exist for each subject, excluding NaNs. The first row contains the number of positive only
             data points. The second row contains the number of negative only data points. The third row contains the number of mixed (positive and
             negative only) data points. If you have positive and negative data and you only want one set of parameters, you need only the third row. If 
             htat is the case, you can set the other rows to 0. 
             
             If subjects 1, 2, 3, and 4 had 100, 130, 124, and 0 positive only data points, 0, 10, 20, and 200 negative only data points, the file would
             look like this:
             100 130 124 0
             100 140 144 200
             100 140 144 200

### Testing use.
 If you are testing and want to check how the code works by generating data from known parameters, set the following parameters to the files continaing  them:
 knownAlpha, knownBeta, knownGamma, knownDelta, knownLambda, knownLuce - These are the known values of the parameters alpha, gamma, and luce where alpha=(0,1), beta=(0,1), gamma=(0,1), delta=(0,1), Lambda=(1,5), luce=(.5,1.5)

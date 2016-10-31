# CumulativeProspectTheory
The files uploaded run Bayesian Parameter Estimation for Cumulative Prospect Theory (CPT). The R code loads in the data and runs RStan to execute the Bayesian parameter estimation. 

There are six files. Two files for data with only positive valued gamble options. Two files for data with only negative valued gample options. Two files for data with positive and negative valued gamble options. According to CPT, there are two sets of parameters: one for positive values and one for negative values. Therefore, separating the estimation into positive, negative, and mixed is necessary, as the absence of positive or negative values means there is no data to estimate parameters.

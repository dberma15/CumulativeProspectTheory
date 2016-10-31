functions {
  real probWeight(real p, real param) {
	  return  pow(p, param)/pow(pow(p,param)+pow(1-p,param),inv(param));
  }
}

data{
  int<lower=1> N;
  
  int items_n_neg;    #integer indicating where the positive data starts
  int rawdata[items_n_neg]; #1s and 0s indicating subject's choice. Rawdata equals 1 means the subject choose option B.
  real beta_upper;
  real beta_lower;
  real delta_upper;
  real delta_lower;
  real lambda_upper;
  real lambda_lower;
  real luce_upper;
  real luce_lower;
  #The prospect data should contain four columns:
  #Column 1: Value1
  #Column 2: probability of Value1
  #Column 3: Value2
  #Column 4: probability of Value4
  
  matrix[items_n_neg,4] prospects_a_neg; #Option A, negative data 
  matrix[items_n_neg,4] prospects_b_neg; #Option B, negative data 
}

#Transformed data breaks the prospects into individual vectors so it's easier to manipulate. 
transformed data{

  vector[items_n_neg] x_a_neg;
  vector[items_n_neg] y_a_neg;
  vector[items_n_neg] px_a_neg;
  vector[items_n_neg] py_a_neg;
  
  vector[items_n_neg] x_b_neg;
  vector[items_n_neg] y_b_neg;
  vector[items_n_neg] px_b_neg;
  vector[items_n_neg] py_b_neg;


  #negative data
  for (ii in 1:(items_n_neg)){
    x_a_neg[ii]<-prospects_a_neg[ii,1];
    y_a_neg[ii]<-prospects_a_neg[ii,3];
    px_a_neg[ii]<-prospects_a_neg[ii,2];
    py_a_neg[ii]<-prospects_a_neg[ii,4];
    
    x_b_neg[ii]<-prospects_b_neg[ii,1];
    y_b_neg[ii]<-prospects_b_neg[ii,3];
    px_b_neg[ii]<-prospects_b_neg[ii,2];
    py_b_neg[ii]<-prospects_b_neg[ii,4];
  }

 
}

#parameters of interest
parameters{
  real<lower=beta_lower,upper=beta_upper> beta;
  real<lower=delta_lower,upper=delta_upper> delta;
  real<lower=lambda_lower, upper=lambda_upper> lambda;
  real<lower=luce_lower, upper=luce_upper> luce;
}


model {
#Declaring Variables

  vector[items_n_neg] v_x_a_neg;
  vector[items_n_neg] v_y_a_neg;
  vector[items_n_neg] v_x_b_neg;
  vector[items_n_neg] v_y_b_neg;
  

  vector[items_n_neg] w_x_a_neg;
  vector[items_n_neg] w_x_b_neg;
  vector[items_n_neg] w_y_a_neg;
  vector[items_n_neg] w_y_b_neg;

  vector[items_n_neg] vf_a_neg;
  vector[items_n_neg] vf_b_neg;

  vector[items_n_neg] binval;

  
#Actual code
  
  #priors for the variables of interest
  beta~cauchy(0.5,1);
  delta~cauchy(0.5,1);
  lambda~uniform(1,5);
  luce~uniform(0.5,1.5);
  

  
  #Loops through negative data
  for (ii in 1:(items_n_neg)){
    
    #Applies value function to data
    v_x_a_neg[ii]<- -lambda*pow(fabs(x_a_neg[ii]),beta);
    v_y_a_neg[ii]<- -lambda*pow(fabs(y_a_neg[ii]),beta);
    v_x_b_neg[ii]<- -lambda*pow(fabs(x_b_neg[ii]),beta);
    v_y_b_neg[ii]<- -lambda*pow(fabs(y_b_neg[ii]),beta);
    
    #Applied probability weighting function to data
    w_x_a_neg[ii]<-probWeight(px_a_neg[ii],delta);
    w_x_b_neg[ii]<-probWeight(px_b_neg[ii],delta);
    w_y_a_neg[ii]<-probWeight(py_a_neg[ii],delta);
    w_y_b_neg[ii]<-probWeight(py_b_neg[ii],delta);

    #combines the weighted probabilites with the value function to form weighted prospects
    vf_a_neg[ii]<-(v_y_a_neg[ii]*w_y_a_neg[ii])+(v_x_a_neg[ii]*w_x_a_neg[ii]);
    vf_b_neg[ii]<-(v_y_b_neg[ii]*w_y_b_neg[ii])+(v_x_b_neg[ii]*w_x_b_neg[ii]);
  }


  #These three loops combine the weighted prospects into one variable.



  #Applies the Luce Choice rule to find the probability of choosing option B over A.
  for (i in 1:items_n_neg){
    binval[i]<-(luce*(-vf_a_neg[i]+vf_b_neg[i]));
  }

  #likelihood function. Takes rawdata and applies a bernouli distribution for the probability of picking option B over A. 
  #Rawdata equals 1 means the subject choose option B.
  for (ii in 1:items_n_neg){
    rawdata[ii]~bernoulli_logit(binval[ii]);
  }
}
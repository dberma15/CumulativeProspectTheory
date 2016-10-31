functions {
  real probWeight(real p, real param) {
	return  pow(p, param)/pow((pow(p,param)+pow(1-p,param)),inv(param));
  }
}

data{
  int<lower=1> N;
  
  int items_n_pos;    #integer indicating where the positive data starts
  #int items_n_neg;    #integer indicating where the negative data starts
  int rawdata[items_n_pos]; #1s and 0s indicating subject's choice. Rawdata equals 1 means the subject choose option B.
  real alpha_upper;
  real alpha_lower;
  real gamma_upper;
  real gamma_lower;
  real luce_upper;
  real luce_lower;
  #The prospect data should contain four columns:
  #Column 1: Value1
  #Column 2: probability of Value1
  #Column 3: Value2
  #Column 4: probability of Value4
  matrix[items_n_pos,4] prospects_a_pos;  #Option A, positive data.
  #matrix[items_n_neg-items_n_pos,4] prospects_a_neg; #Option A, negative data 
  matrix[items_n_pos,4] prospects_b_pos;  #Option B, positive data 
  #matrix[items_n_neg-items_n_pos,4] prospects_b_neg; #Option B, negative data 
}

#Transformed data breaks the prospects into individual vectors so it's easier to manipulate. 
transformed data{
  vector[items_n_pos] x_a_pos;
  vector[items_n_pos] y_a_pos;
  vector[items_n_pos] px_a_pos;
  vector[items_n_pos] py_a_pos;

  vector[items_n_pos] x_b_pos;
  vector[items_n_pos] y_b_pos;
  vector[items_n_pos] px_b_pos;
  vector[items_n_pos] py_b_pos;

  #vector[items_n_neg-items_n_pos] x_a_neg;
  #vector[items_n_neg-items_n_pos] y_a_neg;
  #vector[items_n_neg-items_n_pos] px_a_neg;
  #vector[items_n_neg-items_n_pos] py_a_neg;
  
  #vector[items_n_neg-items_n_pos] x_b_neg;
  #vector[items_n_neg-items_n_pos] y_b_neg;
  #vector[items_n_neg-items_n_pos] px_b_neg;
  #vector[items_n_neg-items_n_pos] py_b_neg;


  #positive data
  for (ii in 1:items_n_pos){
    x_a_pos[ii]<-prospects_a_pos[ii,1];
    y_a_pos[ii]<-prospects_a_pos[ii,3];
    px_a_pos[ii]<-prospects_a_pos[ii,2];
    py_a_pos[ii]<-prospects_a_pos[ii,4];
    
    x_b_pos[ii]<-prospects_b_pos[ii,1];
    y_b_pos[ii]<-prospects_b_pos[ii,3];
    px_b_pos[ii]<-prospects_b_pos[ii,2];
    py_b_pos[ii]<-prospects_b_pos[ii,4];
  }

  #negative data
  #for (ii in 1:(items_n_neg-items_n_pos)){
  # x_a_neg[ii]<-prospects_a_neg[ii,1];
  # y_a_neg[ii]<-prospects_a_neg[ii,3];
  #  px_a_neg[ii]<-prospects_a_neg[ii,2];
  #  py_a_neg[ii]<-prospects_a_neg[ii,4];
  #  
  #  x_b_neg[ii]<-prospects_b_neg[ii,1];
  #  y_b_neg[ii]<-prospects_b_neg[ii,3];
  #  px_b_neg[ii]<-prospects_b_neg[ii,2];
  #  py_b_neg[ii]<-prospects_b_neg[ii,4];
  #}

 
}

#parameters of interest
parameters{
  real<lower=alpha_lower,upper=alpha_upper> alpha;
  #real<lower=beta_lower,upper=beta_upper> beta;
  #real<lower=delta_lower,upper=delta_upper> delta;
  real<lower=gamma_lower,upper=gamma_upper> gamma;
  #real<lower=lambda_lower, upper=lambda_upper> lambda;
  real<lower=luce_lower, upper=luce_upper> luce;
}


model {
#Declaring Variables
  vector[items_n_pos] v_x_a_pos;
  vector[items_n_pos] v_y_a_pos;
  vector[items_n_pos] v_x_b_pos;
  vector[items_n_pos] v_y_b_pos;

  #vector[items_n_neg-items_n_pos] v_x_a_neg;
  #vector[items_n_neg-items_n_pos] v_y_a_neg;
  #vector[items_n_neg-items_n_pos] v_x_b_neg;
  #vector[items_n_neg-items_n_pos] v_y_b_neg;
  
  
  vector[items_n_pos] w_x_a_pos;
  vector[items_n_pos] w_x_b_pos;
  vector[items_n_pos] w_y_a_pos;
  vector[items_n_pos] w_y_b_pos;
  
  #vector[items_n_neg-items_n_pos] w_x_a_neg;
  #vector[items_n_neg-items_n_pos] w_x_b_neg;
  #vector[items_n_neg-items_n_pos] w_y_a_neg;
  #vector[items_n_neg-items_n_pos] w_y_b_neg;


  #vector[items_n_pos] vf_a_pos;
  #vector[items_n_pos] vf_b_pos;

  #vector[items_n_neg-items_n_pos] vf_a_neg;
  #vector[items_n_neg-items_n_pos] vf_b_neg;

  #vector[items_n_neg] vf_a;
  #vector[items_n_neg] vf_b;
  vector[items_n_pos] vf_a;
  vector[items_n_pos] vf_b;
  
  #vector[items_n_neg] binval;
  vector[items_n_pos] binval;

  
#Actual code
  
  #priors for the variables of interest
  alpha~cauchy(0.5,1);
  #beta~cauchy(0.5,1);
  gamma~cauchy(0.5,1);
  #delta~cauchy(0.5,1);
  #lambda~uniform(1,5);
  luce~uniform(0.5,1.5);
  

  #Loops through positive data. 
  for (ii in 1:items_n_pos){

    #Applies value function to data
    v_x_a_pos[ii]<-pow(x_a_pos[ii],alpha);
    v_y_a_pos[ii]<-pow(y_a_pos[ii],alpha);
    v_x_b_pos[ii]<-pow(x_b_pos[ii],alpha);
    v_y_b_pos[ii]<-pow(y_b_pos[ii],alpha);
    
    
    #Applied probability weighting function to data
    w_x_a_pos[ii]<-probWeight(px_a_pos[ii],gamma);
    w_x_b_pos[ii]<-probWeight(px_b_pos[ii],gamma);
    w_y_a_pos[ii]<-probWeight(py_a_pos[ii],gamma);
    w_y_b_pos[ii]<-probWeight(py_b_pos[ii],gamma);
    
    #combines the weighted probabilites with the value function to form weighted prospects
    #vf_a_pos[ii]<- v_x_a_pos[ii]*w_x_a_pos[ii]+v_y_a_pos[ii]*w_y_a_pos[ii];
    #vf_b_pos[ii]<- v_x_b_pos[ii]*w_x_b_pos[ii]+v_y_b_pos[ii]*w_y_b_pos[ii];
  }
  
  #Loops through negative data
 # for (ii in 1:(items_n_neg-items_n_pos)){
 #   
 #   #Applies value function to data
 #   v_x_a_neg[ii]<- -lambda*pow(fabs(x_a_neg[ii]),beta);
 #   v_y_a_neg[ii]<- -lambda*pow(fabs(y_a_neg[ii]),beta);
 #   v_x_b_neg[ii]<- -lambda*pow(fabs(x_b_neg[ii]),beta);
 #   v_y_b_neg[ii]<- -lambda*pow(fabs(y_b_neg[ii]),beta);
 #   
 #   #Applied probability weighting function to data
 #   w_x_a_neg[ii]<-probWeight(px_a_neg[ii],delta);
 #   w_x_b_neg[ii]<-probWeight(px_b_neg[ii],delta);
 #   w_y_a_neg[ii]<-probWeight(py_a_neg[ii],delta);
 #   w_y_b_neg[ii]<-probWeight(py_b_neg[ii],delta);
#
 #   #combines the weighted probabilites with the value function to form weighted prospects
 #   vf_a_neg[ii]<- v_x_a_neg[ii]*w_x_a_neg[ii]+v_y_a_neg[ii]*w_y_a_neg[ii];
 #   vf_b_neg[ii]<- v_x_b_neg[ii]*w_x_b_neg[ii]+v_y_b_neg[ii]*w_y_b_neg[ii];
 # }


  #These three loops combine the weighted prospects into one variable.
  for (j in 1:items_n_pos){
    vf_a[j]<-v_x_a_pos[j]*w_x_a_pos[j]+v_y_a_pos[j]*w_y_a_pos[j];
	vf_b[j]<-v_x_b_pos[j]*w_x_b_pos[j]+v_y_b_pos[j]*w_y_b_pos[j];
  }
  #for (j in 1:items_n_pos){
#	vf_a[j+items_n_pos]<-v_x_b_pos[j];
#  }
#  for (j in 1:items_n_pos){
#	vf_a[j+items_n_pos*2]<-v_y_a_pos[j];
 # }
 # for (j in 1:items_n_pos){
#	vf_a[j+items_n_pos*3]<-v_y_b_pos[j];
 # }
 # for (j in 1:(items_n_neg-items_n_pos)){
 #   vf_a[j+items_n_pos]<-vf_a_neg[j];
 #   vf_b[j+items_n_pos]<-vf_b_neg[j];
 # }

  #Applies the Luce Choice rule to find the probability of choosing option B over A.
  #for (i in 1:items_n_neg){
  for (i in 1:items_n_pos){
    binval[i]<-(luce*(vf_b[i]-vf_a[i]));
  }

  #likelihood function. Takes rawdata and applies a bernouli distribution for the probability of picking option B over A. 
  #Rawdata equals 1 means the subject choose option B.
  #for (ii in 1:items_n_neg){
  for (ii in 1:items_n_pos){
    rawdata[ii]~bernoulli_logit(binval[ii]);
  }
}
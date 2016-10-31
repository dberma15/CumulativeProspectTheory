functions {
  real probWeight(real p, real param) {
	  return  pow(p, param)/pow((pow(p,param)+pow(1-p,param)),inv(param));
  }

}

data{
  int items_n_mix;    #integer indicating where the positive data starts
  int rawdata[items_n_mix]; #1s and 0s indicating subject's choice. Rawdata equals 1 means the subject choose option B.
  #The prospect data should contain four columns:
  #Column 1: Value1
  #Column 2: probability of Value1
  #Column 3: Value2
  #Column 4: probability of Value4
  matrix[items_n_mix,4] prospects_a_mix;  #Option A, positive data.
  matrix[items_n_mix,4] prospects_b_mix;  #Option B, positive data 
}

#Transformed data breaks the prospects into individual vectors so it's easier to manipulate. 
transformed data{
  vector[items_n_mix] a_x;
  vector[items_n_mix] a_y;
  vector[items_n_mix] pa_x;
  vector[items_n_mix] pa_y;

  vector[items_n_mix] b_x;
  vector[items_n_mix] b_y;
  vector[items_n_mix] pb_x;
  vector[items_n_mix] pb_y;


  for (ii in 1:items_n_mix){
    a_x[ii]<-prospects_a_mix[ii,1];
    a_y[ii]<-prospects_a_mix[ii,3];
    pa_x[ii]<-prospects_a_mix[ii,2];
    pa_y[ii]<-prospects_a_mix[ii,4];
    
    b_x[ii]<-prospects_b_mix[ii,1];
    b_y[ii]<-prospects_b_mix[ii,3];
    pb_x[ii]<-prospects_b_mix[ii,2];
    pb_y[ii]<-prospects_b_mix[ii,4];
  }
}

#parameters of interest
parameters{
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=1> beta;
  real<lower=0,upper=1> gamma;
  real<lower=0,upper=1> delta;
  real<lower=1, upper=5> lambda;
  real<lower=.5, upper=1.5> luce;
}


model {
#Declaring Variables
  vector[items_n_mix] v_a_x;
  vector[items_n_mix] v_a_y;
  vector[items_n_mix] v_b_x;
  vector[items_n_mix] v_b_y;

  vector[items_n_mix] w_a_x;
  vector[items_n_mix] w_b_x;
  vector[items_n_mix] w_a_y;
  vector[items_n_mix] w_b_y;
  
  vector[items_n_mix] vf_a;
  vector[items_n_mix] vf_b;
  
  vector[items_n_mix] binval;
  
#Actual code
  
  #priors for the variables of interest
  alpha~cauchy(0.5,1);
  beta~cauchy(0.5,1);
  gamma~cauchy(0.5,1);
  delta~cauchy(0.5,1);
  lambda~uniform(1,5);
  luce~uniform(0.5,1.5);
  
  #Loops through positive data. 
  for (ii in 1:items_n_mix){

    #Applies value function to data
    if (a_x[ii]>=0){
      v_a_x[ii]<- pow(fabs(a_x[ii]),alpha);
      w_a_x[ii]<- probWeight(pa_x[ii],gamma);
    }
    if (a_x[ii]<0){
      v_a_x[ii]<- -lambda*pow(fabs(a_x[ii]),beta);
      w_a_x[ii]<- probWeight(pa_x[ii],delta);
    }
    if (a_y[ii]>=0){
      v_a_y[ii]<- pow(fabs(a_y[ii]),alpha);
      w_a_y[ii]<- probWeight(pa_y[ii],gamma);
    }
    if (a_y[ii]<0){
      v_a_y[ii]<- -lambda*pow(fabs(a_y[ii]),beta);
      w_a_y[ii]<- probWeight(pa_y[ii],delta);
    }
    if (b_x[ii]>=0){
      v_b_x[ii]<- pow(fabs(b_x[ii]),alpha);
      w_b_x[ii]<- probWeight(pb_x[ii],gamma);
    }
    if (b_x[ii]<0){
      v_b_x[ii]<- -lambda*pow(fabs(b_x[ii]),beta);
      w_b_x[ii]<- probWeight(pb_x[ii],delta);
    }
    if (b_y[ii]>=0){
      v_b_y[ii]<- pow(fabs(b_y[ii]),alpha);
      w_b_y[ii]<- probWeight(pb_y[ii],gamma);
    }
    if (b_y[ii]<0){
      v_b_y[ii]<- -lambda*pow(fabs(b_y[ii]),beta);
      w_b_y[ii]<- probWeight(pb_y[ii],delta);
    }
  }
  #These three loops combine the weighted prospects into one variable.
  for (j in 1:items_n_mix){
    vf_a[j]<-v_a_x[j]*w_a_x[j]+v_a_y[j]*w_a_y[j];
  	vf_b[j]<-v_b_x[j]*w_b_x[j]+v_b_y[j]*w_b_y[j];
  }
  #Applies the Luce Choice rule to find the probability of choosing option B over A.
  for (i in 1:items_n_mix){
    binval[i]<-(luce*(vf_b[i]-vf_a[i]));
  }

  #likelihood function. Takes rawdata and applies a bernouli distribution for the probability of picking option B over A. 
  #Rawdata equals 1 means the subject choose option B.
  for (ii in 1:items_n_mix){
    rawdata[ii]~bernoulli_logit(binval[ii]);
  }
}
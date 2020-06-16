# MIT License
# Copyright (c) 2020 Abhinav Prakash, Vijay Panchang, Yu Ding, and Lewis Ntaimo
  
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
  
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

############## FUNCTION DEFINITIONS #################################

#function to calculate the log-likelihood for nonstationary GEV model with trend in location parameter
GEV.likelihood = function(mu,muslope,sigma,epsilon,data,time_index){
  t = (1+epsilon*((data-mu-(muslope*time_index))/sigma))^(-1/epsilon)
  pdf = (t^(epsilon+1))*exp(-t)/sigma
  log.likelihood = sum(log(pdf))
  return(log.likelihood)
}

#function to create time index for any baseyear
time.index = function(Year,ReferenceYear){
  t = Year - ReferenceYear
  return(t)
}
#function to calculate nonstationary return period for event z
calc.nstat.ret.period = function(mu_0,mu_1,scale,shape,z,baseyear,analysisyear){
  i = 0
  accept_prob = 1
  prod_accept_prob = 1
  t = 1
  while (accept_prob > 0){
    accept_prob = exp(-((1 + (shape*((z-(mu_0+mu_1*(analysisyear-baseyear)+(mu_1*i)))/scale)))^(-1/shape)))
    prod_accept_prob = prod_accept_prob*accept_prob
    t = t + prod_accept_prob
    i = i+1
  }
  return(t)
}

#function to calculate nonstationary "n" year return event (n>=2)
calc.nstat.ret.event = function(mu_0,mu_1,scale,shape,n,baseyear,analysisyear,epsilon){
  stat_params.default = 1
  if (mu_1<=0){
    z = mu_0  - scale*(1-(-log(1 - (1/n)))^(-shape))/shape
    return(z)
  } else {
    #calculate lower bound of n year event
    lb = mu_0  - scale*(1-(-log(1 - (1/n)))^(-shape))/shape
    t_lb = n
    while (t_lb > n-epsilon){
      t_lb = calc.nstat.ret.period(mu_0,mu_1,scale,shape,lb,baseyear,analysisyear)
      if(t_lb>n-epsilon){lb = 0.9*lb}
    }
    #calculate upper bound of n year event
    ub = mu_0 + mu_1*(analysisyear-baseyear) +mu_1*n - scale*(1-(-log(1 - (1/(n))))^(-shape))/shape 
    t_ub = n
    while (t_ub < n+epsilon){
      t_ub = calc.nstat.ret.period(mu_0,mu_1,scale,shape,ub,baseyear,analysisyear)
      if(t_ub<n+epsilon){ub = 1.5*ub}
    }
    #get n year event iteratively with epsilon level of accuracy.
    t = 0
    while (abs(t - n) > epsilon) { 
      z = (lb+ub)/2
      t = calc.nstat.ret.period(mu_0,mu_1,scale,shape,z,baseyear,analysisyear)
      if (t < n) {
        lb = z
      } else {ub = z}
    }
  return(z)
  }
  
}

#function to plot return event with respect to return period
ret.event.plot.data = function(baseyear_par_data, min_ret_index, max_ret_index,mle_year_index, avg_year_index, analysisyear){
  t = seq(5,100,5)
  epsilon_seq = 0.01*t
  min = NULL
  max = NULL
  mle = NULL
  avg = NULL
  for (i in 1:length(t)){
    min[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[min_ret_index],baseyear_par_data$location_slope[min_ret_index],baseyear_par_data$scale[min_ret_index],baseyear_par_data$shape[min_ret_index],t[i],baseyear_par_data$year[min_ret_index],analysisyear,epsilon_seq[i])
    max[i] = calc.nstat.ret.event(baseyear_par_data$location_intercept[max_ret_index],baseyear_par_data$location_slope[max_ret_index],baseyear_par_data$scale[max_ret_index],baseyear_par_data$shape[max_ret_index],t[i],baseyear_par_data$year[max_ret_index],analysisyear,epsilon_seq[i])
    mle[i] = calc.nstat.ret.event(baseyear_par_data$location_intercept[mle_year_index],baseyear_par_data$location_slope[mle_year_index],baseyear_par_data$scale[mle_year_index],baseyear_par_data$shape[mle_year_index],t[i],baseyear_par_data$year[mle_year_index],analysisyear,epsilon_seq[i])
    avg[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[avg_year_index],baseyear_par_data$location_slope[avg_year_index],baseyear_par_data$scale[avg_year_index],baseyear_par_data$shape[avg_year_index],t[i],baseyear_par_data$year[avg_year_index],analysisyear,epsilon_seq[i]) 
  }
  y_data = data.frame(t,min,max,mle,avg)
  write.table(y_data, file = paste0(output_file,"_reteventplot_data.csv"), row.names=FALSE, na="", col.names=c("year","min","max","mle","avg"), sep=",")
  return(y_data)
}
#function to plot return event trend with respect to analysis year for a given return period
ret.event.trend = function(baseyear_par_data,min_ret_index, max_ret_index,mle_year_index, avg_year_index, return_period,first_analysis_year,last_analysis_year,interval){
  epsilon = 0.01*return_period
  analysisyear = seq(first_analysis_year,last_analysis_year,interval)
  min = NULL
  max = NULL
  mle = NULL
  avg = NULL
  for (i in 1:length(analysisyear)){
    min[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[min_ret_index],baseyear_par_data$location_slope[min_ret_index],baseyear_par_data$scale[min_ret_index],baseyear_par_data$shape[min_ret_index],return_period,baseyear_par_data$year[min_ret_index],analysisyear[i],epsilon)
    max[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[max_ret_index],baseyear_par_data$location_slope[max_ret_index],baseyear_par_data$scale[max_ret_index],baseyear_par_data$shape[max_ret_index],return_period,baseyear_par_data$year[max_ret_index],analysisyear[i],epsilon)
    mle[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[mle_year_index],baseyear_par_data$location_slope[mle_year_index],baseyear_par_data$scale[mle_year_index],baseyear_par_data$shape[mle_year_index],return_period,baseyear_par_data$year[mle_year_index],analysisyear[i],epsilon)
    avg[i] =  calc.nstat.ret.event(baseyear_par_data$location_intercept[avg_year_index],baseyear_par_data$location_slope[avg_year_index],baseyear_par_data$scale[avg_year_index],baseyear_par_data$shape[avg_year_index],return_period,baseyear_par_data$year[avg_year_index],analysisyear[i],epsilon)
  }
  ret_trend_data = data.frame(analysisyear,min,max,mle,avg)
  write.table(ret_trend_data, file = paste0(output_file,"_retevent_trend_data.csv"), row.names=FALSE, na="", col.names=c("analysisyear","min","max","mle","avg"), sep=",")
  return(ret_trend_data)
}

#function to calculate likelihood value
likelihood.val = function(data,time_index,mu0,mu1,scale,shape){
  k =(1 + (shape*(data-(mu0+(mu1*time_index)))/scale))^(-1/shape)
  pdf = (1/scale)*((k)^(shape+1))*exp(-k)
  llf = -sum(log(pdf))
  return(llf)
}

#function to calculate stationary return event
calc.stat.ret.event = function(mu, scale, shape, n){
  z = mu  - scale*(1-(-log(1 - (1/n)))^(-shape))/shape
  return(z)
}

#function to do MCMC sampling
draw.MCMC.samples = function(Data,time_index,InitialValues,LocationInterceptPrior,LocationSlopePrior,ScalePrior,ShapePrior,ProposalWidth,n_samples){
  mu_lb = LocationInterceptPrior[1]
  mu_ub = LocationInterceptPrior[2]
  muslope_lb = LocationSlopePrior[1]
  muslope_ub = LocationSlopePrior[2]
  scale_lb = ScalePrior[1]
  scale_ub = ScalePrior[2]
  shape_lb = ShapePrior[1]
  shape_ub = ShapePrior[2]
  mu = InitialValues[1]
  muslope = InitialValues[2]
  scale = InitialValues[3]
  shape = InitialValues[4]
  mu_proposal_width = ProposalWidth[1]
  muslope_proposal_width = ProposalWidth[2]
  scale_proposal_width = ProposalWidth[3]
  shape_proposal_width = ProposalWidth[4]
  mu_post=NULL
  muslope_post=NULL
  scale_post=NULL
  shape_post=NULL
  rng_seed = 1
  sample = 0
  rng = 0
  set.seed(rng_seed)
  
  #start the algorithm; (NOTE: if likelihood not finite for initial values of parameters, select different values).
  
  if(is.finite(GEV.likelihood(mu,muslope,scale,shape,Data,time_index))==FALSE){
    while (is.finite(GEV.likelihood(mu,muslope,scale,shape,Data,time_index))==FALSE){
      mu = runif(1,mu_lb,mu_ub)
      muslope = runif(1,muslope_lb,muslope_ub)
      scale = runif(1,scale_lb,scale_ub)
      shape = runif(1,shape_lb,shape_ub)
      rng = rng+1
    }
  }
  
  while (sample < n_samples){
    
    #start sampling proposals from jump distribution
    mu_proposed = rnorm(1,mu,mu_proposal_width)
    muslope_proposed = rnorm(1,muslope,muslope_proposal_width)
    scale_proposed = rnorm(1,scale,scale_proposal_width)
    shape_proposed = rnorm(1,shape,shape_proposal_width)
    #counter for random number generated
    rng = rng+1
    
    #check if all the values are finite for acceptance ratio calculation
    
    suppressWarnings(if (is.finite(GEV.likelihood(mu_proposed,muslope_proposed,scale_proposed,shape_proposed,Data,time_index))==TRUE & is.finite(dunif(mu_proposed,mu_lb,mu_ub,log = TRUE)) == TRUE & is.finite(dunif(scale_proposed,scale_lb,scale_ub,log = TRUE)) == TRUE & is.finite(dunif(shape_proposed,shape_lb,shape_ub,log = TRUE)) == TRUE & is.finite(dunif(muslope_proposed,muslope_lb,muslope_ub,log = TRUE))==TRUE) {
      
      log.r = GEV.likelihood(mu_proposed,muslope_proposed,scale_proposed,shape_proposed,Data,time_index) + dunif(mu_proposed,mu_lb,mu_ub,log = TRUE) + dunif(muslope_proposed,muslope_lb,muslope_ub,log = TRUE) + dunif(scale_proposed,scale_lb,scale_ub,log = TRUE) + dunif(shape_proposed,shape_lb,shape_ub,log = TRUE) - GEV.likelihood(mu,muslope,scale,shape,Data,time_index) - dunif(mu,mu_lb,mu_ub,log = TRUE) - dunif(scale,scale_lb,scale_ub,log = TRUE) - dunif(shape,shape_lb,shape_ub,log = TRUE) - dunif(muslope,muslope_lb,muslope_ub,log = TRUE)
      
      
      #collect new sample if criteria satisfied
      if (log.r > log(runif(1,0,1))){
        mu = mu_proposed
        muslope = muslope_proposed
        scale = scale_proposed
        shape = shape_proposed
        
      }
    })
    mu_post = c(mu_post,mu)
    scale_post = c(scale_post,scale)
    shape_post = c(shape_post,shape)
    muslope_post = c(muslope_post,muslope)
    sample = sample + 1
  }
  parameters = data.frame(mu_post,muslope_post,scale_post,shape_post)
  return(parameters)
}
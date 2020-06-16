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

source("function.r")
inputfile = "hackensack.csv"
outputfile = "Results/hackensack_constrained_MCMC_"
outputplot = "Hackensack River"
data_col = 3
year_col = 1
ReferenceYear = 2018
ReturnPeriod = 100
burn_in = 500
CSVFile = read.csv(file = inputfile, header = T)
Data = CSVFile[,data_col]
Year = CSVFile[,year_col]
time_index = time.index(Year,ReferenceYear)
parameters = draw.MCMC.samples(Data,time_index,InitialValues = c(50,1,20,0.8),LocationInterceptPrior = c(0,200),LocationSlopePrior = c(-5,5),ScalePrior = c(-50,50),ShapePrior = c(-1,1),ProposalWidth = c(2,0.15,1,0.05),n_samples=10000)
parameters = parameters[-c(1:burn_in),]
mu_post = parameters[,1]
muslope_post = parameters[,2]
scale_post = parameters[,3]
shape_post = parameters[,4]


### Expected and MAP Values of the parameters and the return event 
bw = c(2,0.02,0.75,0.020,10)
d_mu = density(mu_post, bw = bw[1])
d_muslope = density(muslope_post, bw = bw[2])
d_scale = density(scale_post, bw= bw[3])
d_shape = density(shape_post, bw = bw[4])
mu_exp = mean(mu_post)
mu_map = d_mu$x[which.max(d_mu$y)]
muslope_exp = mean(muslope_post)
muslope_map = d_muslope$x[which.max(d_muslope$y)]
scale_exp = mean(scale_post)
scale_map = d_scale$x[which.max(d_scale$y)]
shape_exp =mean(shape_post)
shape_map = d_shape$x[which.max(d_shape$y)]
NLLH_exp = likelihood.val(Data,time_index,mu_exp,muslope_exp,scale_exp,shape_exp)
NLLH_map = likelihood.val(Data,time_index,mu_map,muslope_map,scale_map,shape_map)
cat("Expected Values of the Parameters:",'\n')
cat("Location Intercept:", mu_exp, '\n')
cat("Location Slope:", muslope_exp, '\n')
cat("Scale:",scale_exp, '\n')
cat("Shape:", shape_exp, '\n')
cat("NLLH:", NLLH_exp, '\n')
cat('\n',"MAP values of the Parameters:",'\n')
cat("Location Intercept:", mu_map, '\n')
cat("Location Slope:", muslope_map, '\n')
cat("Scale:",scale_map, '\n')
cat("Shape:", shape_map, '\n')
cat("NLLH:", NLLH_map, '\n')


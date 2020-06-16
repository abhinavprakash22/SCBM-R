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

source("function.r") #source the file with all the subroutines
input_file = "assunpink.csv" #file name for the dataset
output_file = "assunpink" #prefix for output plots
data_col = 3 #column in the dataset with the values of flow in cubic meter per second
year_col = 1 #column in the daatset with the year information
col_names = TRUE #boolean to specify whether the dataset contains the column names as the first line

CSVFile = read.csv(file = input_file, header = col_names) #read the file
Data = CSVFile[,data_col] #find column with flow data
Year = CSVFile[,year_col] #find column with year data
referenceyear = 2018 #set a reference year (t=0) for all the computations 
time_index = time.index(Year,referenceyear) #generate time index, t = 1,2,3,...
lmfit = lm(Data~time_index) #Use a linear regression model to estimate the initial values for mu0 and mu1
mu0.int = lmfit$coefficients[1]
mu1.int = abs(lmfit$coefficients[2])
scale.int = sd(Data) #Use standard deviation as the initial value for scale parameter
set.seed(2) #setting the random seed
shape.int = runif(1,-1,1) #randomly initialiaze the shape parameter

#define the objective function for optimization
llf = function(x) {
    k =(1 + (x[4]*(Data-(x[1]+(x[2]*time_index)))/x[3]))^(-1/x[4])
    pdf = (1/x[3])*((k)^(x[4]+1))*exp(-k)
    -sum(log(pdf))
}

#evaluate the objective function at the intial values; if the value is not finite, try a different value to initialize shape parameter
result = try(
    expr =  optim(par = c(mu0.int,mu1.int,scale.int,shape.int), method = "L-BFGS-B", fn = llf, lower = c(0,-Inf,0,-Inf)), silent = TRUE
)
while (class(result) == "try-error"){
    k = k+1
    shape.int = runif(1,-1,1)
    result = try(
      expr =  optim(par = c(mu0.int,mu1.int,scale.int,shape.int), method = "L-BFGS-B", fn = llf, lower = c(0,-Inf,0,-Inf)), silent = TRUE
    )
}

#extract the optimal values for the parameters
mu0 = result$par[1]
mu1 = result$par[2]
scale = result$par[3]
shape = result$par[4]
NLLH = result$value

#use the optimal values for the parameters to calculate the 100-year return event
retevent = calc.nstat.ret.event(mu0,mu1,scale,shape,100,baseyear = referenceyear,analysisyear = referenceyear,1)

#Display the results
cat('Maximum Likelihood Estimation (Unconstrained)','\n')
cat('mu0 = ',mu0,'mu1 = ',mu1,'scale = ',scale,'shape = ',shape,'100yr return event = ',retevent,"NLLH =",NLLH,'\n')


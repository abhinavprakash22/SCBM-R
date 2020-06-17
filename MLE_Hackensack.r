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
input_file = "hackensack.csv"
output_file = "hackensack"
data_col = 3
year_col = 1
col_names = T
CSVFile = read.csv(file = input_file, header = col_names)
Data = CSVFile[,data_col]
Year = CSVFile[,year_col]
referenceyear = 2018
time_index = time.index(Year,referenceyear)
lmfit = lm(Data~time_index)
mu0.int = lmfit$coefficients[1]
mu1.int = abs(lmfit$coefficients[2])
scale.int = sd(Data)
set.seed(2)
shape.int = runif(1,-1,1)
llf = function(x) {
    k =(1 + (x[4]*(Data-(x[1]+(x[2]*time_index)))/x[3]))^(-1/x[4])
    pdf = (1/x[3])*((k)^(x[4]+1))*exp(-k)
    -sum(log(pdf))
}
result = try(
    expr =  optim(par = c(mu0.int,mu1.int,scale.int,shape.int), method = "L-BFGS-B", fn = llf, lower = c(0,-Inf,0,-Inf)), silent = TRUE
)
while (class(result) == "try-error"){
    shape.int = runif(1,-1,1)
    result = try(
      expr =  optim(par = c(mu0.int,mu1.int,scale.int,shape.int), method = "L-BFGS-B", fn = llf, lower = c(0,-Inf,0,-Inf)), silent = TRUE
    )
  }
  mu0 = result$par[1]
  mu1 = result$par[2]
  scale = result$par[3]
  shape = result$par[4]
  NLLH = result$value
  
  
cat('Maximum Likelihood Estimation (Unconstrained)','\n')
cat('mu0 = ',mu0,'mu1 = ',mu1,'scale = ',scale,'shape = ',shape,"NLLH =",NLLH,'\n')

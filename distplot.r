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

par2015 = c(55,0.3,15,0.2027349)
quantile = seq(0,230,1) 
mu2015 = par2015[1]
sigma = par2015[3]
epsilon = par2015[4]
mu1965 = mu2015 - 50*par2015[2]
t2015 = (1+epsilon*((quantile-mu2015)/sigma))^(-1/epsilon)
pdf2015 = (t2015^(epsilon+1))*exp(-t2015)/sigma
t1965 = (1+epsilon*((quantile-mu1965)/sigma))^(-1/epsilon)
pdf1965 = (t1965^(epsilon+1))*exp(-t1965)/sigma
pdf = data.frame(pdf1965,pdf2015)
pdf(file = "Results/distplot.pdf")
matplot(quantile,pdf,type='l',lwd=2,main= 'Nonstationary GEV distribution', xlab = NA, ylab = 'Density', cex.main = 1.33, cex.lab = 1.33, cex.axis = 1.33)
legend('topright', legend = c('Year 1965', 'Year 2015'),lty = c(1,2), col = c(1,2), cex = 1.33)
dev.off()
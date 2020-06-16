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

input_file = "hackensack.csv"
CSVFile = read.csv(file = input_file, header = T)
data_col = 3
year_col = 1
Data = CSVFile[,data_col]
Year = CSVFile[,year_col]
lmfit = lm(Data~Year)
cat("Slope: ",lmfit$coefficients[2],'\n')
fit_summary = summary(lmfit)
cat('p-value of slope:', fit_summary$coefficients[2,4],'\n')
cat('R2 of fit:', fit_summary$r.squared,'\n')
pdf("Results/TimeSeriesHackensack.pdf")
plot(Year,Data, col = "red", pch = 19,main ="Hackensack River", xlab = "Year",ylab = NA,cex.lab = 1.33,cex.axis =1.33,cex.main = 1.33)
abline(lm(Data~Year), col = "blue", lwd = 2)
dev.off()

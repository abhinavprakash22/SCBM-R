# Sign Constrained Bayesian Inference for Nonstationary Models of Extreme Events

## Description:
R code for reproducing the results  in 
*Prakash, A., Panchang, V., Ding, Y., & Ntaimo, L. (2020). Sign Constrained Bayesian Inference for Nonstationary Models of Extreme Events. Journal of Waterway, Port, Coastal, and Ocean Engineering, 146(5), 04020029.* available at this [link](https://ascelibrary.org/doi/10.1061/%28ASCE%29WW.1943-5460.0000589).

## Datasets:

1. `assunpink.csv`
2. `hackensack.csv`

## Execution:
To execute the code make sure the relevant data files are in the same directory as the code file. Many of the R files generate plots and store it in the `Results` sub-directory. Make sure you create a `Results` sub-directory before running the code. The code does not create the sub-directory on its own.
There are two ways to execute the code:
1. The code can be executed directly from command line using the command:
`Rscript filename.r`
2. From `R` console or `Rstudio` using the command: `source('filename.r')`

The R file `function.r` contains all the sub routines developed for the implementation and cannot be executed independently. It is called internally by other `R` files.



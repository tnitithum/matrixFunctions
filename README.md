
# matrixFunctions

The goal of matrixFunctions is to supply fundamental matrix functions
for matrix differentiation numerical calculations. The following is a
list of some of the notable features:

-   Implementation of Richardson’s extrapolation method for vector
    derivatives of 3rd and 4th-order, and general matrix derivatives of
    up to 1st-order. For vector derivatives of 1st and 2nd order, please
    refer to `numDeriv` package.
-   Numeric improvements for some matrix functions in `matrixcalc`
    package.
-   Matrix derivatives of Euclidean norm up to 4th order.
-   Matrix derivative of spherical integral of up to 4th order.

Note that while there exist many definitions of Matrix Derivative in
literature, our calculations are based on the following definition:

∂*F*/∂*X* = (∂/∂*X*) ⊗ *F*

for some matrix function *F*, some matrix *X*, and Kronecker product ⊗.

## Installation

The package can be installed by running the following R code.

``` r
library(devtools)

install_github("tnitithum/matrixFunctions")
```

## Example

This is an example of numeric vector differentiation of third order
applied to the Euclidean norm.

``` r
library(matrixFunctions)

f <- function(x) norm(x, type="2")
x <- 1:2
vector_deriv_3rd(f, x)
##             [,1]        [,2]
## [1,] -0.21466254 -0.07155418
## [2,] -0.07155418  0.12521981
## [3,] -0.07155418  0.12521981
## [4,]  0.12521981 -0.10733128
D_norm_3rd(x) # Truth
##             [,1]        [,2]
## [1,] -0.21466253 -0.07155418
## [2,] -0.07155418  0.12521981
## [3,] -0.07155418  0.12521981
## [4,]  0.12521981 -0.10733126
```

Fourth order derivative of the Euclidean norm.

``` r
vector_deriv_4th(f, x)
##               [,1]        [,2]        [,3]        [,4]
## [1,]  5.156024e-07  0.21466240  0.21466240 -0.03577708
## [2,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
## [3,]  2.146624e-01 -0.03577708 -0.03577708 -0.10733140
## [4,] -3.577708e-02 -0.10733140 -0.10733140  0.16099882
matrix_deriv(D_norm_3rd, x, transpose_deriv_operator=TRUE)
##               [,1]        [,2]        [,3]        [,4]
## [1,]  2.665172e-12  0.21466253  0.21466253 -0.03577709
## [2,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [3,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [4,] -3.577709e-02 -0.10733126 -0.10733126  0.16099689
D_norm_4th(x) # Truth
##               [,1]        [,2]        [,3]        [,4]
## [1,] -2.081668e-17  0.21466253  0.21466253 -0.03577709
## [2,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [3,]  2.146625e-01 -0.03577709 -0.03577709 -0.10733126
## [4,] -3.577709e-02 -0.10733126 -0.10733126  0.16099689
```

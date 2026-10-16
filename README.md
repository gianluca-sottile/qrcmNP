
<!-- README.md is generated from README.Rmd. Please edit that file -->


# qrcmNP

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/qrcmNP)](https://CRAN.R-project.org/package=qrcmNP)
[![R-CMD-check](https://github.com/gianluca-sottile/qrcmNP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gianluca-sottile/qrcmNP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

'qrcmNP' extends quantile regression coefficients modeling (QRCM, Frumento and
Bottai, 2016; package 'qrcm') in two directions:

* `piqr()`: **penalized** QRCM for covariate selection using the information on
  all quantiles simultaneously (Sottile, Frumento, Chiodi and Bottai, 2020), with
  lasso and group lasso penalties;
* `niqr()`: **nonlinear** parametric quantile models, nonlinear in the parameters,
  in the covariates, or both (Bottai and Cilluffo, 2020).

## Installation

```r
# CRAN version
install.packages("qrcmNP")
# development version
# install.packages("remotes")
remotes::install_github("gianluca-sottile/qrcmNP")
```

## Penalized QRCM

``` r
library(qrcmNP)
set.seed(1234)
n <- 300
d <- data.frame(x1 = rexp(n), x2 = runif(n, 0, 5), x3 = rnorm(n), x4 = runif(n))
u <- runif(n)
d$y <- 1 + qnorm(u) + (1 + 2 * u) * d$x1 + d$x2
fit <- piqr(y ~ x1 + x2 + x3 + x4, formula.p = ~ slp(p, 2), data = d, nlambda = 30)
best <- gof.piqr(fit, method = "BIC", plot = FALSE)
best
#> 
#> Selection of lambda by BIC 
#> 
#> Selected lambda: 9.051  (position 7 )
#> Degrees of freedom: 8 
#> 
#> Coefficients:
#>              (Intercept)  slp1      slp2    
#> (Intercept)   0.17508      1.57385   0.42889
#> x1            0.32843      0.93906  -0.38235
#> x2            0.89827           NA  -0.07583
#> x3                 NA           NA        NA
#> x4                 NA           NA        NA
#> 
m <- refit_piqr(fit, pos.lambda = best$posMinLambda)
summary(m)
#> 
#> Call: qrcm::iqr(formula = y ~ x1 + x2 + x3 + x4, formula.p = ~slp(p, 
#>     2), data = d, s = structure(c(1, 1, 1, 0, 0, 1, 1, 1, 0, 
#> 0, 1, 1, 1, 0, 0), dim = c(5L, 3L)))
#> 
#> converged: TRUE 
#> n. of iterations: 13 
#> n. of free parameters: 9 
#> n. of observations: 300 
#> 
#> ###################### 
#> ###################### 
#> 
#> Coefficients:
#>              (Intercept)  slp1      slp2    
#> (Intercept)  -0.63063      1.67485   0.18364
#> x1            0.94140      0.81264  -0.16361
#> x2            1.03362     -0.01875  -0.05584
#> x3                 NA           NA        NA
#> x4                 NA           NA        NA
#> 
#> Standard errors:
#>              (Intercept)  slp1     slp2   
#> (Intercept)  0.44791      0.24718  0.24268
#> x1           0.27057      0.15578  0.14520
#> x2           0.10471      0.07298  0.06559
#> x3                NA           NA       NA
#> x4                NA           NA       NA
#> 
#> ###################### 
#> ###################### 
#> 
#> Wald test for x:
#>             chi-square df P(> chi)    
#> (Intercept)      75.15  3 3.37e-16 ***
#> x1              442.89  3  < 2e-16 ***
#> x2              343.53  3  < 2e-16 ***
#> x3                  NA  0       NA    
#> x4                  NA  0       NA    
#> 
#> 
#> Wald test for b(p):
#>             chi-square df P(> chi)    
#> (Intercept)    556.503  3   <2e-16 ***
#> slp1           528.999  3   <2e-16 ***
#> slp2             2.728  3    0.436    
#> 
#> Minimized loss function: 128.1989
#> 
```

## Nonlinear parametric quantile models

``` r
set.seed(1234)
x <- runif(n)
fun <- function(theta, p) cbind(theta[1] + exp(theta[2] * p), theta[3] + theta[4] * p)
b <- fun(c(1, 1, 1, 1), runif(n))
dn <- data.frame(y = b[, 1] + b[, 2] * x, x = x)
m2 <- niqr(y ~ x, data = dn, fun = fun, x0 = rep(0.5, 4))
summary(m2)
#> 
#> Call: niqr(formula = y ~ x, data = dn, fun = fun, x0 = rep(0.5, 4))
#> 
#> n. of observations: 300    n. of iterations: 7 
#> Standard errors: sandwich 
#> 
#>        Estimate Std. Error z value Pr(>|z|)    
#> theta0  1.05523    0.06841  15.426  < 2e-16 ***
#> theta1  1.00759    0.03903  25.816  < 2e-16 ***
#> theta2  1.02334    0.14501   7.057 1.70e-12 ***
#> theta3  1.04105    0.20787   5.008 5.49e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Minimized loss function: 0.1852556 
#> AIC = 63.5767   BIC = 78.3918   AIC* = -1.67269   BIC* = -1.64799 
#> 
```

## References

* Frumento, P. and Bottai, M. (2016). Parametric modeling of quantile regression
  coefficient functions. *Biometrics*, 72(1), 74-84. <https://doi.org/10.1111/biom.12410>
* Sottile, G., Frumento, P., Chiodi, M. and Bottai, M. (2020). A penalized approach to
  covariate selection through quantile regression coefficient models. *Statistical
  Modelling*, 20(4), 369-385. <https://doi.org/10.1177/1471082X19825523>
* Bottai, M. and Cilluffo, G. (2020). Nonlinear parametric quantile models.
  *Statistical Methods in Medical Research*, 29(12), 3757-3769.
  <https://doi.org/10.1177/0962280220941159>

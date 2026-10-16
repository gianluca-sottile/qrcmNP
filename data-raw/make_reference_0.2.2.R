# Reference outputs of qrcmNP 0.2.2 (CRAN), used by regression tests of the
# "legacy" algorithm. Run with qrcmNP 0.2.2 installed (not the dev version).
# Output: tests/testthat/fixtures/ref_0.2.2.rds
library(qrcmNP)
stopifnot(packageVersion("qrcmNP") == "0.2.2")
ref <- list()

## piqr, example 1 of ?piqr ---------------------------------------------------
set.seed(1234)
n <- 300
x1 <- rexp(n); x2 <- runif(n, 0, 5)
b <- function(p){matrix(cbind(1, qnorm(p), slp(p, 2)), nrow = 4, byrow = TRUE)}
theta <- matrix(0, nrow = 3, ncol = 4); theta[, 1] <- 1; theta[1, 2] <- 1; theta[2:3, 3] <- 2
qy <- function(p, theta, b, x){rowSums(x * t(theta %*% b(p)))}
y <- qy(runif(n), theta, b, cbind(1, x1, x2))
d1 <- data.frame(y = y, x1 = x1, x2 = x2)
s <- matrix(1, nrow = 3, ncol = 4); s[1, 3:4] <- 0; s[2:3, 2] <- 0
obj <- piqr(y ~ x1 + x2, formula.p = ~ I(qnorm(p)) + slp(p, 2), s = s, nlambda = 50,
            data = d1, display = FALSE)
ref$piqr1 <- list(data = d1, s = s, lambda = obj$lambda, coefficients = obj$coefficients,
                  obj.function = obj$obj.function, df = obj$df, dl = obj$dl,
                  aic = gof.piqr(obj, method = "AIC", plot = FALSE)$posMinLambda,
                  bic = gof.piqr(obj, method = "BIC", plot = FALSE)$posMinLambda)

## piqr, sparse model, default slp basis ---------------------------------------
set.seed(1234)
n <- 500; q <- 5
X <- matrix(abs(rnorm(n * q)), n, q); colnames(X) <- paste0("X", 1:q)
th <- rbind(c(3, 1.5, 1, 1), c(2, 1, 1, 1), 0, 0, c(1.5, 1, 1, 1), 0)
B <- function(p){cbind(1, slp(p, 3))}
y <- rowSums(cbind(1, X) * (B(runif(n)) %*% t(th)))
d2 <- data.frame(y = y, X)
obj2 <- piqr(y ~ X1 + X2 + X3 + X4 + X5, data = d2, nlambda = 30, display = FALSE)
ref$piqr2 <- list(data = d2, lambda = obj2$lambda, coefficients = obj2$coefficients,
                  obj.function = obj2$obj.function, df = obj2$df,
                  bic = gof.piqr(obj2, method = "BIC", plot = FALSE)$posMinLambda)

## niqr, example of ?niqr ------------------------------------------------------
set.seed(1234)
n <- 300
x <- runif(n)
fun <- function(theta, p){
  beta0 <- theta[1] + exp(theta[2]*p)
  beta1 <- theta[3] + theta[4]*p
  cbind(beta0, beta1)}
beta <- fun(c(1, 1, 1, 1), runif(n))
y <- beta[, 1] + beta[, 2]*x
m <- niqr(fun = fun, x0 = rep(.5, 4), X = cbind(1, x), y = y)
ref$niqr1 <- list(data = data.frame(y = y, x = x), coef = m$x, se = m$se, fx = m$fx,
                  iter = m$iter, code = m$code)

ref$info <- list(version = "0.2.2", R = R.version.string, date = Sys.Date())
saveRDS(ref, "tests/testthat/fixtures/ref_0.2.2.rds", version = 2)

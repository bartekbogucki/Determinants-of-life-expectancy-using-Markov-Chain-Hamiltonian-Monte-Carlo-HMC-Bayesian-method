#Cleaning
rm(list=ls())

#Libraries
if(!require("manipulate")) {install.packages("manipulate"); library(manipulate)}
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require("readxl")) {install.packages("readxl"); library(readxl)}
if(!require("coda")) {install.packages("coda"); library(coda)}
if(!require("moments")) {install.packages("moments"); library(moments)}
if(!require("rstan")) {install.packages("rstan"); library(rstan)}
if(!require("readxl")) {install.packages("readxl"); library(readxl)}
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#working directory
#setwd("D:/sgh/mgr 3/ekonometra Bayesowska/PD1/Ekonometria_Bayesowska_PD1_Bartosz_Bogucki_82801")
#if(!is.null(dev.list())) dev.off()
#cat("\014")

#import
dane <- read_excel("Bartosz.Bogucki.82801.xlsx")

#OLS
ols=lm(Life_expectancy~Education+GII+Ratio+Urban+Births, data=dane)
summary(ols)

y <- as.matrix(dane[, c('Life_expectancy')])
N.data <- length(y)
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(dane[, c('Education', 'GII', 'Ratio', 'Urban', 'Births')]))
Beta.ols.data <- ols$coefficients
v.data <- ols$df.residual
XTX.data <- t(X) %*% X
s2.data <- sum((ols$residuals) ^ 2) / v.data

#A priori knowledge

Beta.prior <- c(67.89, 1.087, -0.136, -0.113, 0.07, -0.527)
sm2.prior <- s2.data
U.prior <- diag(6)
U.prior[1,1]=4.534^2
U.prior[2,2]=0.274^2
U.prior[3,3]=0.041^2
U.prior[4,4]=0.024^2
U.prior[5,5]=0.027^2
U.prior[6,6]=0.16^2
v.prior <- 187
vs2.prior <- v.prior / sm2.prior
k <- ncol(X)
N <- length(y)

#A posteriori

Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data-Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)
(Beta.posterior)

#Densisty 

beta.space <- rbind(seq(from = 60, to=80, length.out=400),
                    seq(from=-2, to=2, length.out=400), 
                    seq(from=-0.3, to=0.3, length.out=400), 
                    seq(from=-0.2, to=0.2, length.out=400), 
                    seq(from=-0.1, to=0.1, length.out=400),
                    seq(from=-1, to=1, length.out=400))

n_eval_points <- length(beta.space[1,])
n_parameters <- length(Beta.posterior)
prior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
posterior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
for(ii in 1:length(Beta.posterior)) {
  prior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space[ii,]), 1, dmvt,
                                      delta = Beta.prior[ii], sigma = as.matrix(U.prior[ii, ii] / sm2.prior), df = v.prior, log = FALSE)
  posterior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space[ii,]), 1, dmvt,
                                          delta = Beta.posterior[ii], sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), df = v.posterior, log = FALSE)
}

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)

dev.new()
layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 5, 0), ncol = 4, byrow = TRUE))
for(ii in 2:length(Beta.posterior)) {
  plot(beta.space[ii,], prior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = grey_area,
       ylim = c(0, max(c(max(prior.marg.dens.beta[ii, ]),max(posterior.marg.dens.beta[ii, ]))) + 1), type = "l", ylab = "g?sto??", main = colnames(dane)[ii + 1],
       xlim = c(min(beta.space[ii,]), max(beta.space[ii,])))
  polygon(c(beta.space[ii,], rev(beta.space[ii,])), c(prior.marg.dens.beta[ii, ], 
                                                      rep(0, length(beta.space[ii,]))), col = grey_area, border = NA)
  abline(v = Beta.prior[ii], col = grey_line, lwd = 3)
  text(Beta.prior[ii], max(prior.marg.dens.beta[ii, ]), paste("E(beta) a priori = ", Beta.prior[ii]), col = grey_line)
  abline(v = Beta.ols.data[ii], col = rgb(0, 0, 0, 1), lwd = 3)
  text(Beta.ols.data[ii], max(posterior.marg.dens.beta[ii, ]), paste("parametr OLS = ", round(Beta.ols.data[ii], 4)), col = rgb(0, 0, 0, 1))
  lines(beta.space[ii,], posterior.marg.dens.beta[ii, ], lwd = 2, col = green_line)
  polygon(c(beta.space[ii,], rev(beta.space[ii,])), c(posterior.marg.dens.beta[ii, ], 
                                                      rep(0, length(beta.space[ii,]))), col = green_area, border = NA)
  abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
  text(Beta.posterior[ii], quantile(posterior.marg.dens.beta[ii, ], 0.96), paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
  
}

#HPDI
x1 <- rep(0,5)
x2 <- rep(0,5)
conf_level = 0.95

red_area <- rgb(255, 100, 123, 80, names = NULL, maxColorValue = 255)
red_line <- rgb(200, 0, 30, 160, names = NULL, maxColorValue = 255)

dev.new()
layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 5, 0), ncol = 4, byrow = TRUE))
for(ii in 2:6) {
  df <- data.frame(1:n_eval_points, posterior.marg.dens.beta[ii, ])
  df <- df[order(df[, 2], decreasing = TRUE), ]
  df <- data.frame(df, cumsum(df[, 2]) < conf_level*max(cumsum(df[,2])))
  df <- df[order(df[, 1], decreasing = FALSE), ]
  credible_set_indicator <- as.vector(as.integer(df[, 3]))
  credible_set_begin <- match(1, credible_set_indicator)
  credible_set_end <- length(credible_set_indicator) - match(1, rev(credible_set_indicator))
  x1[ii-1] <- beta.space[ii,credible_set_begin]
  x2[ii-1] <- beta.space[ii,credible_set_end]
  plot(beta.space[ii,], prior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = grey_area,
       ylim = c(0, max(c(max(prior.marg.dens.beta[ii, ]),max(posterior.marg.dens.beta[ii, ]))) + 1), type = "l", ylab = "g?sto??", main = colnames(dane)[ii + 1],
       xlim = c(min(beta.space[ii,]), max(beta.space[ii,])))
  polygon(c(beta.space[ii,], rev(beta.space[ii,])), c(prior.marg.dens.beta[ii, ], 
                                                      rep(0, length(beta.space[ii,]))), col = grey_area, border = NA)
  abline(v = Beta.prior[ii], col = grey_line, lwd = 3)
  text(Beta.prior[ii], max(prior.marg.dens.beta[ii, ]), paste("E(beta) a priori = ", Beta.prior[ii]), col = grey_line)
  abline(v = Beta.ols.data[ii], col = rgb(0, 0, 0, 1), lwd = 3)
  text(Beta.ols.data[ii], max(posterior.marg.dens.beta[ii, ]), paste("parametr OLS = ", round(Beta.ols.data[ii], 4)), col = rgb(0, 0, 0, 1))
  lines(beta.space[ii,], posterior.marg.dens.beta[ii, ], lwd = 2, col = green_line)
  polygon(c(beta.space[ii,], rev(beta.space[ii,])), c(posterior.marg.dens.beta[ii, ], 
                                                      rep(0, length(beta.space[ii,]))), col = green_area, border = NA)
  abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
  text(Beta.posterior[ii], quantile(posterior.marg.dens.beta[ii, ], 0.96), paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
  
  abline(v = x1[ii-1], col = red_line, lwd = 3)
  abline(v = x2[ii-1], col = red_line, lwd = 3)
  
  posterior.cs <- posterior.marg.dens.beta[ii, ] * credible_set_indicator
  polygon(c(beta.space[ii,], rev(beta.space[ii,])), 
          c(posterior.cs, rep(0, length(beta.space[ii,]))), 
          col = red_area, border = NA)
  text(x2[ii-1], quantile(posterior.marg.dens.beta[ii, ], 0.93), paste("95% przedzia? HPDI: (", round(x1[ii-1], digits = 2), " , ", round(x2[ii-1], digits = 2), ")"), col = red_line)
}

#interval HPDI
HPDI_table <- data.frame(names(Beta.ols.data[2:6]), x1, x2)
colnames(HPDI_table) <- c("zmienna", "lewy brzeg 95% przedzialu HPDI", "prawy brzeg 95% przedzialu HPDI")
(HPDI_table)

#Bayes coef

#Model 1: 
P_y_M1 <- ((det(U.posterior) ^ 0.5) * gamma(v.posterior / 2) * ((vs2.posterior) ^ (- v.posterior / 2))) / ((pi ^ (N.data / 2)) * (det(U.prior) ^ 0.5) * gamma(v.prior / 2) * ((vs2.prior) ^ (- v.prior / 2)))

#P_y_M1 need to be logarithmised:
l_P_y_M1 <- 0.5*log(det(U.posterior)) + log(gamma(v.posterior / 2)) - v.posterior/2*log(vs2.posterior) -
  (N.data/2*log(pi) + 0.5*det(U.prior) + log(gamma(v.prior/2)) - v.prior/2*log(vs2.prior))

#Model 2
l_P_y_M2 <- rep(NA, 5)
for (ii in 2:6) {
  X_2 <- X[, -c(ii)]
  eval(parse(text = paste("OLS_results_2 <- lm(Life_expectancy ~ ", substring(paste(paste("`",colnames(X_2),"`", sep=""), collapse = "+"),4), ", data = dane)", sep = "")))
  Beta.ols.data_2 <- OLS_results_2$coefficients
  v.data_2 <- OLS_results_2$df.residual
  XTX.data_2 <- t(X_2) %*% X_2
  s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2
  
  Beta.prior_2 <- Beta.prior[-c(ii)]
  U.prior_2 <- U.prior[- c(ii), -c(ii)]

  sm2.prior_2 <- sm2.prior
  v.prior_2 <- v.prior
  vs2.prior_2 <- vs2.prior
  
  Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
  U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
  v.posterior_2 <- v.prior_2 + N.data
  vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
  sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)

  l_P_y_M2[ii - 1] <- 0.5*log(det(U.posterior_2)) + log(gamma(v.posterior_2 / 2)) - v.posterior_2/2*log(vs2.posterior_2) -
    (N.data/2*log(pi) + 0.5*det(U.prior_2) + log(gamma(v.prior_2/2)) - v.prior_2/2*log(vs2.prior_2))
}

#BF = P_y_M1 / P_y_M2, therefore l_BF = l_P_y_M1 - l_P_y_M2:
l_BF <- as.vector(l_P_y_M1) - as.vector(l_P_y_M2)

#Thus:
options(scipen = 999)
BF <- round(exp(l_BF),2)
BF_table <- data.frame(names(Beta.ols.data[2:6]), BF)
colnames(BF_table) <- c("zmienna", "czynnik Bayesa (analitycznie)")
(BF_table)
options(scipen = 0)

x1 <- dane$Education
x2 <- dane$GII
x3 <- dane$Ratio
x4 <- dane$Urban
x5 <- dane$Births
y <- dane$Life_expectancy
N <- length(y)
X <- cbind(matrix(1, N, 1), x1, x2, x3, x4, x5)

beta_prior <- c(0, 0, 0, 0, 0, 0)
beta_sd_prior <- c(10, 10, 10, 10, 10, 10)
k <- length(beta_prior)

fit <- stan(file = "HMC.stan",
            data = c("N", "k", "y", "x1", "x2", "x3", "x4", "x5", "beta_prior", "beta_sd_prior"),
            iter = 2000,
            warmup = 500,
            seed = 42,
            chains = 4
)
print(fit)
plot(fit)

simulated.chains.all <- extract(fit,
                                inc_warmup = FALSE,
                                par = c("h", "beta"))

#Visualization of density a posteriori
par(mfrow = c(2,3))
for (ii in  1:6) {
  hist(simulated.chains.all$beta[,ii], main=paste("Histogram of beta",ii), xlab=paste("Beta",ii))
}

#HPDI
chain.b <- mcmc(simulated.chains.all$beta)
hpdi.b <- HPDinterval(chain.b)
rownames(hpdi.b) <- c("(Intercept)", "x1", "x2", "x3", "x4", "x5")
(hpdi.b)

rstan::traceplot(fit, pars = c("beta"), inc_warmup = FALSE, nrow = 6)

simulated.chains.1by1 <- extract(fit,
                                 permuted = FALSE,
                                 inc_warmup = FALSE,
                                 par = c("beta"))

combined.chains <- mcmc.list(mcmc(simulated.chains.1by1[, 1, ]), 
                             mcmc(simulated.chains.1by1[, 2, ]),
                             mcmc(simulated.chains.1by1[, 3, ]),
                             mcmc(simulated.chains.1by1[, 4, ]))

plot(combined.chains)
summary(combined.chains, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))

#Scale reduction factors (Koop, p. 66)
gelman.diag(combined.chains)
gelman.plot(combined.chains)

#Geweke statistics (Koop, p. 68)
geweke.diag(combined.chains, frac1 = 0.1, frac2 = 0.5)
geweke.plot(combined.chains, frac1 = 0.1, frac2 = 0.5, nbins = 40, pvalue = 0.05)

#Heidel
heidel.diag(combined.chains)

#Serial correlation in chains
autocorr(combined.chains, lags = c(0, 1, 5, 10, 50), relative = TRUE)
autocorr.diag(combined.chains)
autocorr.plot(combined.chains, lag.max = 50)
#Cross correlation in chains
crosscorr(combined.chains)
par(mfrow = c(1,1))
crosscorr.plot(combined.chains)

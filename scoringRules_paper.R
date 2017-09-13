## ----setup, echo=FALSE---------------------------------------------------
library(knitr)
thm = knit_theme$get("default")
knit_theme$set(thm)

# for inline R code highlighting
knit_hooks$set(inline = function(x) { 
  if (is.numeric(x)) return(knitr:::format_sci(x, 'latex')) 
  highr:::hi_latex(x) 
}) 
opts_chunk$set(concordance=TRUE)
options(scipen = 1, digits = 3)

## ----eval=FALSE----------------------------------------------------------
## crps(y, ...)
## logs(y, ...)
## 
## ## S3 method for class 'numeric'
## crps(y, family, ...)
## logs(y, family, ...)

## ----echo=FALSE----------------------------------------------------------
rm(list=ls())
library(scoringRules)
set.seed(42)

## ------------------------------------------------------------------------
obs <- rnorm(10)
crps(obs, family = "normal", mean = c(1:10), sd = c(1:10))
logs(obs, family = "normal", mean = c(1:10), sd = c(1:10))

## ----echo=FALSE----------------------------------------------------------
plot_prepared_background <- function() {
	plot(NULL, type = "n", xlim = c(0, 9), ylim = c(0, 5), bty = "n", 
	     xlab = "Observation y", ylab = "Score value")
	z <- seq(0, 9, .01)
	bg <- 15 * dgamma(z, shape = 2, scale = 1.5)
	polygon(c(z, rev(z)), c(rep(0, length(bg)), rev(bg)), 
	        col="gray80", border="gray80")
	legend("top", bty = "n", legend = c("LogS", "CRPS"), 
	       col = c("purple", "darkorange"), lty = c(1,1), 
	       lwd = c(2,2))
}

## ----score-illustration, echo=TRUE,  dev='pdf', fig.width=4, fig.height=4----
plot_prepared_background() # initializing plot

crps_y <- function(y) crps(y, family = "gamma", shape = 2, scale = 1.5)
logs_y <- function(y) logs(y, family = "gamma", shape = 2, scale = 1.5)
plot(crps_y, from = 0, to = 9, col = "darkorange", lwd = 2, add = TRUE)
plot(logs_y, from = 0, to = 9, col = "purple", lwd = 2, add = TRUE)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
## crps_norm(y, mean = 0, sd = 1, location = mean, scale = sd)
## logs_norm(y, mean = 0, sd = 1, location = mean, scale = sd)

## ----eval=FALSE----------------------------------------------------------
## crps_sample(y, dat, method = "edf", w = NULL, bw = NULL,
##             num_int = FALSE, show_messages = TRUE)
## logs_sample(y, dat, bw = NULL, show_messages = TRUE)

## ------------------------------------------------------------------------
# single observation
obs <- rnorm(1)
sample <- rnorm(1e4, mean = 2, sd = 3)
crps_sample(obs, dat = sample)  
logs_sample(obs, dat = sample, show_messages = FALSE)

# multiple observations
obs2 <- rnorm(2)
sample2 <- matrix(rnorm(2e4, mean = 2, sd = 3), nrow = 2)
crps_sample(obs2, dat = sample2)
logs_sample(obs2, dat = sample2, show_messages = FALSE)

## ------------------------------------------------------------------------
ngrid <- seq(from = 50, to = length(sample), by = 50)
crps_approx <- logs_approx <- numeric(length(ngrid))
for (i in seq_along(ngrid)) {
  size <- ngrid[i]
  crps_approx[i] <- crps_sample(obs, dat = sample[1:size])
  logs_approx[i] <- logs_sample(obs, dat = sample[1:size],
                                show_messages = FALSE)
}

## ----echo=FALSE----------------------------------------------------------
crps_true <- crps(obs, family = "normal", mean = 2, sd = 3)
logs_true <- logs(obs, family = "normal", mean = 2, sd = 3)

## ----plot1, echo=FALSE, dev='pdf', fig.width=8, fig.height=4-------------
every2nd <- function(x) x[(1:length(x)) %% 2 == 0]
xax <- every2nd(pretty(1:max(ngrid)))

par(mfrow = c(1,2), mar = c(5, 4, 2, 2) + 0.1)

plot(NULL, type = "n", bty = "n",
     main = "CRPS", xlab = "Sample size", ylab = "Score value",
     xlim = c(0, 1e4), ylim = crps_true + c(-0.1, 0.1), xaxt = "n")
axis(1, at = xax)
abline(h = crps_true, lty = 2)
lines(ngrid, crps_approx, col = "darkorange", lwd = 2)

plot(NULL, type = "n", bty = "n",
     main = "LogS", xlab = "Sample size", ylab = "Score value",
     xlim = c(0, 1e4), ylim = logs_true + c(-0.12, 0.12), xaxt = "n")
axis(1, at = xax)
abline(h = logs_true, lty = 2)
lines(ngrid, logs_approx, col = "purple", lwd = 2)

## ------------------------------------------------------------------------
# From Messner et al. (2016):
library(crch)
data(RainIbk)
RainIbk <- sqrt(RainIbk)
RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
RainIbk <- subset(RainIbk, enssd > 0)

## ----tidy=FALSE----------------------------------------------------------
data_train <- subset(RainIbk, as.Date(rownames(RainIbk)) <= "2004-11-30")
data_eval <- subset(RainIbk, as.Date(rownames(RainIbk)) >= "2005-01-01")

CRCHgauss <- crch(rain ~ ensmean | log(enssd), data_train,
                  dist = "gaussian", left = 0)

## ----echo=FALSE----------------------------------------------------------
CRCHlogis <- crch(rain ~ ensmean | log(enssd), data = data_train, 
left = 0, dist = "logistic")
CRCHstud <- crch(rain ~ ensmean | log(enssd), data = data_train, 
left = 0, dist = "student")

## ----tidy=TRUE-----------------------------------------------------------
gauss_mu <- predict(CRCHgauss, data_eval, type = "location")
gauss_sc <- predict(CRCHgauss, data_eval, type = "scale")

ens_fc <- data_eval[, grep('^rainfc', names(RainIbk))]

## ----echo=FALSE----------------------------------------------------------
logis_mu <- predict(CRCHlogis, data_eval, type = "location")
logis_sc <- predict(CRCHlogis, data_eval, type = "scale")
stud_mu <- predict(CRCHstud, data_eval, type = "location")
stud_sc <- predict(CRCHstud, data_eval, type = "scale")
stud_df <- CRCHstud$df

## ----postprocplot, echo=FALSE, warning=FALSE, message=FALSE, dev='pdf', fig.width=6, fig.asp=1/1----
ID.list <- c(206,953,2564)

m <- matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m,heights = c(0.5,0.5))
par(mar = c(5,4,2,1))   

for(ID in ID.list){
  col.logis <- "blue"
  col.gauss <- "green3"
  col.stud <- "darkorange"
  
  z <- seq(0,10,0.01)
  flogis.plot <- flogis(z, logis_mu[ID], logis_sc[ID], lower = 0, lmass = "cens")
  flogis.p0 <- plogis(0, logis_mu[ID], logis_sc[ID])
  fnorm.plot <- fnorm(z, gauss_mu[ID], gauss_sc[ID], lower = 0, lmass = "cens")
  fnorm.p0 <- pnorm(0, gauss_mu[ID], gauss_sc[ID])
  fstud.plot <- ft(z, stud_df, stud_mu[ID], stud_sc[ID], lower = 0, lmass = "cens")
  fstud.p0 <- pt(-stud_mu[ID] / stud_sc[ID], stud_df)
  
  p0.offset <- 0.2
  yrange <- c(-0.025,0.24)
  if(ID == ID.list[3]){yrange <- c(-0.025,0.47)}
  plot(z, flogis.plot, type = "l", bty = "n", col = col.logis, 
       ylim = yrange, xlim = c(-0.4,10),
       ylab = "Density", xlab = "Precipitation amount in mm", 
       main = rownames(data_eval)[ID])
  segments(0, 0, 0, flogis.p0, col = col.logis, lwd = 3)
  lines(z, fnorm.plot, col = col.gauss)
  segments(-p0.offset, 0, -p0.offset, fnorm.p0, col = col.gauss, lwd = 3)
  lines(z, fstud.plot, col = col.stud)
  segments(-2*p0.offset, 0, -2*p0.offset, fstud.p0, col = col.stud, lwd = 3)
  segments(0, 0, 0, flogis.p0, col = col.logis, lwd = 3) 
  segments(data_eval$rain[ID], 0, data_eval$rain[ID], yrange[2], lty = 2)
  ens.fc <- as.numeric(data_eval[, grep('^rainfc',names(RainIbk))][ID,])
  for(j in 1:length(ens.fc)){segments(ens.fc[j], -0.025, ens.fc[j], -0.005)}
}

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("center", inset = 0,
       legend = c("censored logistic", "censored Gaussian", "censored Student t"), 
       lty = rep(1,3), col = c("blue", "green3", "darkorange"), 
       ncol = 1, bty ="n")

## ----tidy=FALSE----------------------------------------------------------
obs <- data_eval$rain
gauss_crps <- crps(obs, family = "cnorm", location = gauss_mu, 
                   scale = gauss_sc, lower = 0, upper = Inf)
ens_crps <- crps_sample(obs, dat = as.matrix(ens_fc))

## ----echo=FALSE----------------------------------------------------------
logis_crps <- crps(obs, family = "clogis", location = logis_mu, 
                   scale = logis_sc, lower = 0, upper = Inf)
stud_crps <- crps(obs, family = "ct", df = stud_df, location = stud_mu, 
                  scale = stud_sc, lower = 0, upper = Inf)

## ----eval=TRUE, echo=FALSE-----------------------------------------------
df <- data.frame(mean(logis_crps), mean(gauss_crps), mean(stud_crps), 
				 sprintf("%1.3f", mean(ens_crps)))
names(df) <- c("CRCHlogis", "CRCHgauss", "CRCHstud", "Ensemble")
print(df, row.names = FALSE)

## ----eval=TRUE, echo=TRUE, message=FALSE, tidy=FALSE---------------------
data(gdp)

# Get training data (2014Q1 vintage; data until 2013Q4)
# Get evaluation data (2015Q1 vintage; data from 2014)
data_train <- subset(gdp, vint == "2014Q1") # 267 observations
data_eval <- subset(gdp, vint == "2015Q1" & grepl("2014", dt)) # 4 obs.

# Draw predictive model parameters
h <- 4L; n <- 20000L
fc_params <- ar_ms(data_train$val, forecast_periods = h, n_rep = n)

## ----eval=TRUE, echo=TRUE, message=FALSE---------------------------------
# Mixture-of-normals approximation
m <- t(fc_params$fcMeans) # matrix[4 * 20000] of means
s <- t(fc_params$fcSds) # matrix[4 * 20000] of standard deviations
w <- matrix(1/n, nrow = h, ncol = n) # matrix of (equal) mixture weights
# Sample approximation
sample <- matrix(rnorm(h * n, mean = m, sd = s), nrow = h, ncol = n)

## ----mcmcplot, echo=FALSE, warning=FALSE, message=FALSE, dev='pdf', fig.width=6, fig.asp=1/1----
fmix <- function(m, s) {
	function(x) {
		40000 * sapply(x, function(z) mean(dnorm(z, mean = m, sd = s)))
	}
}

# Plot histograms
layout(mat = matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE),
heights = c(0.5, 0.5))
par(mar = c(5,4,2,1))
for (jj in seq_along(data_eval$dt)) {
act <- data_eval$val[jj]
x <- sample[jj, ]

hist(x, main = data_eval$dt[jj], xlab = "", yaxt = "n",
xlim = c(-20, 20), ylim = c(0, 8000))
axis(2, at = c(0, 2000, 4000, 6000, 8000))
segments(act, 0, act, 8000, lty = 2)
plot(fmix(m[jj, ], s[jj, ]), from = min(x), to = max(x), 
     lwd = 2, add = TRUE)
}

## ----eval=TRUE, echo=TRUE, message=FALSE---------------------------------
# Compute scores
obs <- data_eval$val
names(obs) <- data_eval$dt
crps_mpe <- crps(obs, "normal-mixture", m = m, s = s, w = w)
logs_mpe <- logs(obs, "normal-mixture", m = m, s = s, w = w)
crps_ecdf <- crps_sample(obs, sample)
logs_kde <- logs_sample(obs, sample, show_message = FALSE)

# Print results
print(cbind(crps_mpe, crps_ecdf, logs_mpe, logs_kde))

## ----eval=FALSE----------------------------------------------------------
## es_sample(y, dat)
## vs_sample(y, dat, w = NULL, p = 0.5)

## ----echo = FALSE--------------------------------------------------------
names(obs) <- NULL

## ------------------------------------------------------------------------
es_sample(obs, dat = sample)
vs_sample(obs, dat = sample)


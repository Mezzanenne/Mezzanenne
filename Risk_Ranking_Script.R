###########################################################################################
# Generic script used in "Quantitative methods for the prioritization of foods implicated in the transmission of hepatititis E to humans in Italy"
# (Moro, Isopi, Suffredini, Schembri, Scavia) to be decline for each food category
###########################################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(stats)
library(MASS)
library(pracma)
library(fBasics)

###########################################################################################
# alpha = probability of contaminated food (estimated from sampled data)
# lambda = viral contamination rate for positive food samples (estimated from sampled data)
# mu = individual infective dose rate (estimated using dose_0 and data)
# gd_i = average portion of food in grams per day (estimated from survey)
# N_s = susceptible population (estimated from census data and from literature)
##########################################################################################

##########################################################################################
##########################             FUNCTION                ###########################
##########################################################################################
p_inf <- function(alpha, lambda, mu, gd_i){
  F_c <- function(x){(1 - alpha)*Heaviside(x, a = 0) + alpha*(1 - exp(-lambda*x))} # Cumulative distribution of viral concentration as a random variables mixture
  f_s <- function(x){dexp(x, mu)} # Density of exponential distribution representing the individual threshold
  
  H <- function(x){F_c(x)*f_s(x)} 
  p_sc <- integral(H, xmin = 0, xmax = Inf)[1] # P(C < S)
  q_sc = 1-p_sc # Exposure probability p after one single meal
  
  dose_i <- 100 # Mean portion of food i in grams (for example 100 g)
  da_i <- round(365*gd_i/dose_i) # Mean number of portion of food i consumed in a year
  p_i_y <- 1-(p_sc)^da_i # Exposure probability in a year
  return(p_i_y)
}

##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
#############################       RAW MODEL OUTPUT        ##############################
##########################################################################################

p_i<- P_inf(alpha, lambda, mu, gd_i) # Evaluation of exposure probability
no_inf <- N_s*p_inf # Mean number of exposed people

##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
###################      OUTPUT PARAMETER AND VISUALIZATION        #######################
##########################################################################################

#### OUTPUT PARAMETER

var_i <- no_inf*(1-p_i) # Variance calculation
sqrt(var_i) # Standard deviation


#### OUTPUT VISUALIZATION

ggplot(data.frame(x = rnorm(10000, mean = no_inf, sd = sqrt(var_inf))), aes(x)) + 
  geom_density(colour = "darkorchid1", size = 1) + theme_bw() +
  labs(x = "number of exposed people per year")



##########################################################################################
####################       CONFIDENCE INTERVAL CALCULATION         #######################
##########################################################################################
# n = Sample size (depending on data)
# x = Sample vector

alph <- .05 # Confidence level 95%

chi.up <- qchisq(alph/2, 2*n)
chi.low <- qchisq(1-alph/2, 2*n)
num <- 2*sum(x)
low <- num/chi.up
up <- num/chi.low

##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
################    BEGIN: UNCERTAINTY AND SENSITIVITY ANALYSIS        ###################
##########################################################################################
source("Param_data.R") #Parameter 10,000 samples 

out <- c() # Output sample vector

for (i in 1:10000){
  out[i] <- p_inf(alpha, lambda_v[i], mu_v[i], gd_v[i])
}

df <- data.frame(out, lambda_v, mu_v, gd_v)


#### DATA MOMENTS
mean(out)
sd(out)
skewness(out)

sd(lambda_v)
sd(mu_v)
sd(gd_v)


#### CORRELATION TEST
cor.test(lambda_v, mu_v)
cor.test(lambda_v, gd_v)
cor.test(mu_v, gd_v)

par <- data.frame(lambda = lambda_v, mu = mu_v, conc_out = gd_v)
pairs(~lambda+mu+conc_out, data = par)


#### VISUALIZATIONS

## LAMBDA
df %>% ggplot(aes(x = lambda_v, fill = lambda_v)) +
  geom_histogram(alpha = 0.7, position ="identity", aes(y = ..density..), color="black", fill = "cornflowerblue") +
  theme_bw() + xlab("lambda samples")
df %>% ggplot(aes(lambda_v)) + stat_ecdf(col = "cornflowerblue", size = 1) + theme_bw() + xlab("lambda") + ylab("cumulative density")
df %>% ggplot(aes(x = lambda_v, y = out)) + geom_point(shape = 1, size = 3) + labs(x = "lambda", y = "output (p_PL)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=0, hjust=0.5)) + theme_bw()


## GD
df %>% ggplot(aes(x = gd_v, fill = gd_v)) +
  geom_histogram(alpha = 0.7, position ="identity", aes(y = ..density..), color="black", fill = "cornflowerblue") +
  xlab("conc_out samples") + theme_bw()
df %>% ggplot(aes(gd_v)) + stat_ecdf(col = "cornflowerblue", size = 1) + theme_bw() + xlab("conc_out") + ylab("cumulative density")
df %>% ggplot(aes(x = gd_v, y = out)) + geom_point(shape = 1, size = 3)+ #VS OUT
  labs(x = "conc_out", y = "output (p_PL)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=0, hjust=0.5)) + theme_bw()

## MU
df %>% ggplot(aes(x = mu_v, fill = mu_v)) +
  geom_histogram(alpha = 0.7, position ="identity", aes(y = ..density..), color="black", fill = "cornflowerblue") +
  theme_bw() + xlab("mu samples")
df %>% ggplot(aes(mu_v)) + stat_ecdf(col = "cornflowerblue", size = 1) + theme_bw() + xlab("mu") + ylab("cumulative density")
df %>% ggplot(aes(x = mu_v, y = out)) + geom_point(shape = 1, size = 3) +
  labs(x = "mu", y = "output (p_PL)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=0, hjust=0.5)) + theme_bw()

#OUT
df %>% ggplot(aes(x = out, fill = out)) +
  geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black", fill = "cornflowerblue") + theme_bw() + xlab("output samples") +
  geom_density(alpha = 0.4, fill = "magenta")
df %>%
  ggplot(aes(x = out)) + stat_ecdf(col = "cornflowerblue") +
  xlab("output (probability of successful exposition per year") + ylab("cumulative probability") + theme_bw()



#### REGRESSION ANALYSIS ####
library(lm.beta) # package containing the function lm.beta to get standardized coefficients
mod <- lm(out ~ -1 + lambda_v + mu_v + gft_v)
# summary(mod)
lm.beta(mod) -> mod.beta
summary(mod.beta)

#### PARAMETER IMPACT ####
mod.beta$standardized.coefficients[1] %>% abs() -> lam
mod.beta$standardized.coefficients[2] %>% abs() -> m
mod.beta$standardized.coefficients[3] %>% abs() -> gf

# Based on ranking of parameter, if mu > gf > lam
(m - gf)/gf
(gf - lam)/lam
(m - lam)/lam


##########################################################################################
################      END: UNCERTAINTY AND SENSITIVITY ANALYSIS        ###################
##########################################################################################

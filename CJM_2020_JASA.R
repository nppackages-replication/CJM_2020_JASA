######################################################################
# Empirical Illustration
# M.D. Cattaneo, M. Jansson, X. Ma
# 24-SEP-2020
######################################################################

# install.packages("ggplot2")
# install.packages("lpdensity")
# install.packages("rddensity")
# install.packages("rdd")

# NOTE: if you are using RDDENSITY version 2020 or newer, the option 
# masspoints=FALSE may be needed to replicate the results in the monograph.
# For example:
#    out = rddensity(X)
# should be replaced by:
#    out = rddensity(X, massPoints=FALSE, regularize=FALSE)

######################################################################
######################################################################
# Headstart Data
######################################################################
######################################################################
rm(list=ls(all=TRUE))
library(foreign); library(MASS); library(Hmisc)
library(rddensity); library(lpdensity); library(ggplot2); library(rdd)

data <- read.csv("headstart.csv", header=T)
data$povrate60 <- data$povrate60 - 59.198
X <- data[!is.na(data[, "povrate60"]), "povrate60"]

######################################################################
# Table: Manipulation Testing
######################################################################
Result <- matrix(NA, ncol=8, nrow=7)
colnames(Result) <- c("hl", "hr", "nhl", "nhr", "t", "pval", "binl", "binr")

j <- 1
model <- rddensity(X, p=1, bwselect="each")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

j <- j + 1
model <- rddensity(X, p=2, bwselect="each")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

j <- j + 1
model <- rddensity(X, p=3, bwselect="each")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

j <- j + 1
model <- rddensity(X, p=1, bwselect="diff")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

j <- j + 1
model <- rddensity(X, p=2, bwselect="diff")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

j <- j + 1
model <- rddensity(X, p=3, bwselect="diff")
Result[j, 1] <- model$h$left    ; Result[j, 2] <- model$h$right
Result[j, 3] <- model$N$eff_left; Result[j, 4] <- model$N$eff_right
Result[j, 5] <- model$test$t_jk ; Result[j, 6] <- model$test$p_jk

McCrary <- DCdensity(runvar=X, cutpoint=0, verbose=TRUE, plot=TRUE, ext.out=TRUE)
Result[7, 7] <- sum(McCrary$data$cellmp < 0); Result[7, 8] <- sum(McCrary$data$cellmp >= 0)
Result[7, 1] <- Result[7, 2] <- McCrary$bw
Result[7, 3] <- sum(McCrary$data$cellmp <  0 & McCrary$data$cellmp >= -1 * McCrary$bw)
Result[7, 4] <- sum(McCrary$data$cellmp >= 0 & McCrary$data$cellmp <=      McCrary$bw)
Result[7, 5] <- McCrary$z; Result[7, 6] <- McCrary$p

round(Result, 3)

######################################################################
# Figure: RD plot, p = 2, different bandwidth on each sides
######################################################################
RDEst <- rddensity(X, bwselect="each") 

rdplotdensity(RDEst, X, plotRange=c(-40, 20), plotN=100, 
                       lcol=c("blue", "red"), xlabel="Running Variable", 
                       histBreaks=seq(-40, 20, length.out=31))$Estplot + 
  geom_vline(xintercept=0, linetype="dashed") + ylim(0, 0.032) + theme(legend.position="none")
#test the performance of JEGN using simulated data
rm(list=ls())

library("JEGN")
source("evaluation_metric.R")
source("generate_data.R")

p = 100
K = 3
G = 4
model = "ER"
umin = 0.5
umax = 1
n = 25
tau = 0.5


rtimes = 5

lambda_list = exp(seq(log(2), log(0.1), length.out = 10));
alpha_list = c(0.15, 0.25, 0.35, 0.45)

Pre = matrix(0, length(lambda_list), length(alpha_list))
Rec = matrix(0, length(lambda_list), length(alpha_list))

for (r in 1:rtimes){
  cat("r=",r)
  # generare simulated data
  dat = generate_data(p, n, tau, K, G, model)

  for (i_lam in 1:length(lambda_list)){
    cat("i_lam=",i_lam)
    for (i_alpha in 1:length(alpha_list)){
      # run JEGN
      result = JEGN.admm(dat$Sigma, lambda_list[i_lam], alpha_list[i_alpha])
      # compute precision and recall
      per = evaluation_metric(dat$Omega, result$Omega)
      Pre[i_lam, i_alpha] = Pre[i_lam, i_alpha] + per$Pre
      Rec[i_lam, i_alpha] = Rec[i_lam, i_alpha] + per$Rec
    }
  }
}

# compute the average of precision and recall
Pre = Pre/rtimes
Rec  = Rec/rtimes


# plot the precision-recall curve
type_list = c("o","o","o","o")
pch_list = c(0,1, 2,3)
plot(Rec[,1], Pre[,1], col="blue", type = type_list[1],
     pch=pch_list[1], lwd=2, xlim = c(0,1),  ylim = c(0,1),  xlab = "Recall", ylab = "Precision")
title(main = "Performance of JEGN on ER network")
for (i in 1:4){
  points(Rec[,i], Pre[,i],
         col="blue",type = type_list[i], pch=pch_list[i], lwd=2)
}

method = c(expression(paste("JEGN ", alpha == 0.15)), expression(paste("JEGN ", alpha == 0.25)),
           expression(paste("JEGN ", alpha == 0.35)), expression(paste("JEGN ", alpha == 0.45)))
pch=  c(0,1, 2,3)
col =  rep("blue",4)
legend(0.5,1.05, legend=  method, ncol=1, col = col, pch=pch, lwd=2,bty = "n")



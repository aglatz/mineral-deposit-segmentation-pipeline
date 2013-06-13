library("R.methodsS3", lib.loc="~/R")
library("R.oo", lib.loc="~/R")
library("R.matlab", lib.loc="~/R")
library("quantreg", lib.loc="~/R/x86_64-redhat-linux-gnu-library/2.15/")

data<-R.matlab::readMat("forR.mat");
hu<-data.frame(vavg=c(log10((data$V+data$V.ref)/2)),
               vdiff=c(2*(data$V-data$V.ref)/(data$V+data$V.ref)),
			   J=c(data$J));

plot(hu$vavg, hu$vdiff)

taus <- c(.05, .25, .5, .75, .95);
for (idx in 1:length(taus)) {
	fit <- rq(vdiff ~ vavg, tau=taus[idx], data=hu)
	summary(fit)
	abline(fit)
}
par(ask=TRUE);

plot(hu$vavg, hu$J)

taus <- c(.05, .25, .5, .75, .95);
for (idx in 1:length(taus)) {
	fit <- rq(J ~ vavg, tau=taus[idx], data=hu)
	summary(fit)
	abline(fit)
}

rm(list=ls())

uu = seq(0, 1, length=1001)


### Select one group
## group 1
xx = seq(0, 10, length=1001)
truth = list(d = function(x) {dgamma(x, 20, sc=0.2)},
			 q = function(u) {qgamma(u, 20, sc=0.2)},
			 t = "gamma(sh=20, sc=0.2)")

pseu = list()
pseu[[1]] = list(d = function(x) {dgamma(x, 2, sc=1)},
				 q = function(u) {qgamma(u, 2, sc=1)},
				 t = "gamma(sh=2, sc=1)")

pseu[[2]] = list(d = function(x) {dgamma(x, 4, sc=1)},
				 q = function(u) {qgamma(u, 4, sc=1)},
				 t = "gamma(sh=4, sc=1)")

pseu[[3]] = list(d = function(x) {dgamma(x, 10, sc=0.4)},
				 q = function(u) {qgamma(u, 10, sc=0.4)},
				 t = "gamma(sh=10, sc=0.4)")

## group 2
xx = seq(-5, 16, length=1001)
truth = list(d = function(x) {dgamma(x, 20, sc=0.2)},
			 q = function(u) {qgamma(u, 20, sc=0.2)},
			 t = "gamma(sh=20, sc=0.2)")

pseu = list()
pseu[[1]] = list(d = function(x) {dt((x-4)/2.5, df=2)/2.5},
				 q = function(u) {qt((u), df=2)*2.5 + 4},
				 t = "t(loc=4, sc=2.5, df=2)")

pseu[[2]] = list(d = function(x) {dt((x-8)/2.5, df=2)/2.5},
				 q = function(u) {qt((u), df=2)*2.5 + 8},
				 t = "t(loc=8, sc=2.5, df=2)")

pseu[[3]] = list(d = function(x) {dt((x-4)/1.2, df=2)/1.2},
				 q = function(u) {qt((u), df=2)*1.2 + 4},
				 t = "t(loc=4, sc=1.2, df=2)")

## group 3
xx = seq(0, 1, length=1001)
truth = list(d = function(x) {dbeta(x, 0.8, 0.2)},
			 q = function(u) {qbeta(u, 0.8, 0.2)},
			 t = "beta(0.8, 0.2)")

pseu = list()

pseu[[1]] = list(d = function(x) {dbeta(x, 2, 2)},
				 q = function(u) {qbeta(u, 2, 2)},
				 t = "beta(2, 2)")

pseu[[2]] = list(d = function(x) {dbeta(x, 1, 1)},
				 q = function(u) {qbeta(u, 1, 1)},
				 t = "unif(0, 1)")

pseu[[3]] = list(d = function(x) {dbeta(x, 0.5, 0.5)},
				 q = function(u) {qbeta(u, 0.5, 0.5)},
				 t = "beta(0.5, 0.5)")


## group 4
xx = seq(-10, 10, length=1001)
truth = list(d = function(x) {dt(x, df=1)},
			 q = function(u) {qt(u, df=1)},
			 t = "t(loc=0, sc=1, df=1)")

pseu = list()
pseu[[1]] = list(d = function(x) {dnorm(x)},
				 q = function(u) {qnorm(u)},
				 t = "normal(0, 1)")

pseu[[2]] = list(d = function(x) {dnorm(x, 0, 3)},
				 q = function(u) {qnorm(u, 0, 3)},
				 t = "normal(0, 3^2)")

pseu[[3]] = list(d = function(x) {dt((x-4)/2, df=3)/2},
				 q = function(u) {qt((u), df=3)*2 + 4},
				 t = "t(loc=4, sc=2, df=3)")


## group 5
xx = seq(-10, 10, length=1001)
truth = list(d = function(x) {dt(x, df=5)},
             q = function(u) {qt(u, df=5)},
             t = "t(df=5)")

pseu = list()
pseu[[1]] = list(d = function(x) {dnorm(x, 0, 10)},
                 q = function(u) {qnorm(u, 0, 10)},
                 t = "normal(0, 10)")

pseu[[2]] = list(d = function(x) {dnorm(x, 0, 1)},
                 q = function(u) {qnorm(u, 0, 1)},
                 t = "normal(0, 1)")

normAppx <- optim(c(0), fn = \(par) -truth$d(par), hessian = TRUE)

pseu[[3]] = list(d = function(x) {dnorm(x, 0, solve(normAppx$hessian) * 3)},
                 q = function(u) {qnorm(u, 0, solve(normAppx$hessian) * 3)},
                 t = paste("normal(0,",round(solve(normAppx$hessian) * 3, 2),")"))


## group 6
xx = seq(-10, 10, length=1001)
truth = list(d = function(x) {dt(x, df=5)},
             q = function(u) {qt(u, df=5)},
             t = "t(loc=0, sc=1, df=5)")

pseu = list()
pseu[[1]] = list(d = function(x) {dnorm(x, 0, 20)},
                 q = function(u) {qnorm(u, 0, 20)},
                 t = "normal(0, 20)")

pseu[[2]] = list(d = function(x) {dt(x, df = 2)},
                 q = function(u) {qt(u, df = 2)},
                 t = "t(loc=0, sc=1, df=2)")

normAppx <- optim(c(0), fn = \(par) -truth$d(par), hessian = TRUE)

pseu[[3]] = list(d = function(x) {dnorm(x, 0, solve(normAppx$hessian) * 3)},
                 q = function(u) {qnorm(u, 0, solve(normAppx$hessian) * 3)},
                 t = paste("normal(0,",round(solve(normAppx$hessian) * 3, 2),")"))

pseu[[4]] = list(d = function(x) {dt(x, df = 4)},
                 q = function(u) {qt(u, df = 4)},
                 t = "t(loc=0, sc=1, df=4)")




### plot
(n_pseu = length(pseu))

m = matrix(c(1,2,3), nrow=1)
layout(mat=m, widths=c(1, 1, 0.7))
par(mar=c(4,1,1,1))

plot(xx, truth$d(xx), type="l", lwd=3,
	ylim=c(0,.6),
	xlab=expression(theta), ylab="", axes=FALSE)
axis(side=1)
for(i in 1:n_pseu) {
lines(xx, pseu[[i]]$d(xx), col=i+1, lwd=2, lty=i+1)	
}

plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=3, 
	ylim=c(0,10),
	xlab=expression(psi), ylab="", axes=FALSE, bty="L")
axis(side=1)

for(i in 1:n_pseu) {
	lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=i+1, lwd=2, lty=i+1)
}

plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left", bty="n", inset=c(0.0, 0.0),
		legend=c(truth$t, sapply(pseu, function(x) x$t)),
		col=1:(n_pseu+1), lwd=2, lty=1:(n_pseu+1), cex = 0.7)







### for poster

xx = seq(-5, 12, length=1001)
truth = list(d = function(x) {dgamma(x, 6, sc=0.5)},
             q = function(u) {qgamma(u, 6, sc=0.5)},
             t = "gamma(sh=6, sc=0.5)")

pseu = list()

pseu[[1]] = list(d = function(x) {dt((x-3)/1.3, df=1)/1.3},
                 q = function(u) {qt((u), df=2)*1.3 + 4},
                 t = "Cauchy(loc=3, sc=1.3)")

pseu[[2]] = list(d = function(x) {dt((x-3)/3, df=1)/3},
                 q = function(u) {qt((u), df=2)*3 + 3},
                 t = "Cauchy(loc=3, sc=3)")

pseu[[3]] = list(d = function(x) {dt((x-7)/3, df=1)/3},
                 q = function(u) {qt((u), df=2)*3 + 7},
                 t = "Cauchy(loc=7, sc=3)")




(n_pseu = length(pseu))
cols = c("dodgerblue3", "goldenrod3", "firebrick3")

pdf(file="transformDemo.pdf", width=8, height=3.75)
par(mar=c(3.5,1,1,1), mfrow=c(1,2))

plot(xx, truth$d(xx), type="l", lwd=6,
     ylim=c(0,.39),
     xlab="x", ylab="", axes=FALSE, cex.lab=1.5, line=2.5)
axis(side=1)
for(i in n_pseu:1) {
  lines(xx, pseu[[i]]$d(xx), col=cols[i], lwd=5, lty=1)	
}
text(0.4, 0.33, expression(g(x)))
text(5.5, 0.2, expression(hat(g)(x)), col=cols[1])
text(5.8, 0.14, expression(hat(g)(x)), col=cols[2])
text(8.8, 0.14, expression(hat(g)(x)), col=cols[3])

par(mar=c(3.5,1,1,1))
plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=5, 
     ylim=c(0,11),
     xlab="u", ylab="", 
     axes=FALSE, cex.lab=1.5, line=2.5)
axis(side=1)
mtext(expression(h(hat(G)^-1*(u))), side=2, line=-1.0, las=3)

for(i in n_pseu:1) {
  lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=cols[i], lwd=5, lty=1)
}

legend("topright", bty="n", inset=c(0.0, 0.0),
       legend=c(truth$t, sapply(pseu, function(x) x$t)), 
       col=c("black", cols), lwd=5, lty=1, cex=1.2)
dev.off()


pdf(file="genericDensity.pdf", width=5, height=4)
xx = seq(-3, 9, length=1001)
par(mar=c(1,1,1,1))
plot(xx, truth$d(xx), type="l", lwd=6,
     xlab="", ylab="", axes=FALSE)
dev.off()

## a function that plots what the transformation of a pseduo target looks like
# author: Sam Johnson

# the function
plot_of_transform <- function(target, pseu, support = c(-10,10)) {
  
  xx <- seq(support[1], support[2], length.out = 1001)
  uu <- seq(0, 1, length=1001)
  
  (n_pseu = length(pseu))
  
  m = matrix(c(1,2,3), nrow=1)
  graphics::layout(mat=m, widths=c(1, 1, 0.7))
  par(mar=c(4,1,1,1))
  
  plot(xx, target$d(xx), type="l", lwd=3,
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
    lines(uu, target$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=i+1, lwd=2, lty=i+1)
  }
  
  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  legend("left", bty="n", inset=c(0.0, 0.0),
         legend=c(target$t, sapply(pseu, function(x) x$t)),
         col=1:(n_pseu+1), lwd=2, lty=1:(n_pseu+1), cex = 1.2)
  
}


# an example
# xx = seq(0, 10, length=1001)
# truth = list(d = function(x) {dinvgamma(x, 2, 1)},
#              q = function(u) {qinvgamma(u, 2, 1)},
#              t = "invgamma(2,1)")
# 
# pseu = list()
# pseu[[1]] = pseudo_t_list(.34, .42, degf = 1, lb = 0)
# pseu[[2]] = pseudo_t_list(.34, .42, degf = 1)
# 
# plot_of_transform(target = truth, pseu = pseu, support = c(0,7))


# ### for poster
# xx = seq(-5, 12, length=1001)
# truth = list(d = function(x) {dgamma(x, 6, sc=0.5)},
#              q = function(u) {qgamma(u, 6, sc=0.5)},
#              t = "Gamma(sh=6, sc=0.5)")
# 
# pseu = list()
# 
# pseu[[1]] = list(d = function(x) {dt((x-3)/1.3, df=1)/1.3},
#                 q = function(u) {qt((u), df=2)*1.3 + 4},
#                 t = "Cauchy(loc=3, sc=1.3)")
# 
# pseu[[2]] = list(d = function(x) {dt((x-3)/3, df=1)/3},
#                 q = function(u) {qt((u), df=2)*3 + 3},
#                 t = "Cauchy(loc=3, sc=3)")
# 
# pseu[[3]] = list(d = function(x) {dt((x-7)/3, df=1)/3},
#                 q = function(u) {qt((u), df=2)*3 + 7},
#                 t = "Cauchy(loc=7, sc=3)")
# 
# 
# (n_pseu = length(pseu))
# cols = c("dodgerblue3", "goldenrod3", "firebrick3")
# 
# pdf(file="transformDemo.pdf", width=8, height=3.75)
# par(mar=c(3.5,1,1,1), mfrow=c(1,2))
# 
# plot(xx, truth$d(xx), type="l", lwd=6,
#     ylim=c(0,.39),
#     xlab=expression(theta), ylab="", axes=FALSE, cex.lab=1.5, line=2.5)
# axis(side=1)
# for(i in n_pseu:1) {
#  lines(xx, pseu[[i]]$d(xx), col=cols[i], lwd=5, lty=1)
# }
# text(0.4, 0.33, expression(f(theta)))
# text(5.5, 0.2, expression(hat(pi)(theta)), col=cols[1])
# text(5.8, 0.14, expression(hat(pi)(theta)), col=cols[2])
# text(8.8, 0.14, expression(hat(pi)(theta)), col=cols[3])
# 
# par(mar=c(3.5,1,1,1))
# uu <- seq(0,1,length.out=1001)
# plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=5,
#     ylim=c(0,11),
#     xlab=expression(psi), ylab="",
#     axes=FALSE, cex.lab=1.5, line=2.5)
# axis(side=1)
# mtext(expression(h(hat(Pi)^-1*(psi))), side=2, line=-1.0, las=3)
# 
# for(i in n_pseu:1) {
#  lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=cols[i], lwd=5, lty=1)
# }
# 
# legend("topright", bty="n", inset=c(0.0, 0.0),
#       legend=c(truth$t, sapply(pseu, function(x) x$t)),
#       col=c("black", cols), lwd=5, lty=1, cex=1.2)
# dev.off()


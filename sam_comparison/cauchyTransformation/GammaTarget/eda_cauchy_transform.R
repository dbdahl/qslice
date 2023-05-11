### exploring using cauchy as a psuedo target for standard normal
scales <- c(0.5,0.75,0.875,1,1.25,1.35,1.5,2)#c(0.5,0.75,0.78125,0.8125,0.875,1,1.25,1.5,2)
color <- RColorBrewer::brewer.pal(length(scales),'Dark2')

curve(dnorm(x), from = -4, to = 4)
for( i in 1:length(scales)) {
  curve(dcauchy(x,scale = scales[i]), add = TRUE, col = color[i], lty = i)
}
legend(x = -3, y = 0.3,
       legend = as.character(scales),
       col = color, 
       lty = 1:length(scales))

# seeing how each it looks after transformation
uu = seq(0, 1, length=1001)
xx = seq(-5, 5, length=1001)
truth = list(d = function(x) {dnorm(x, 0, 1)},
             q = function(u) {qnorm(u, 0, 1)},
             t = "norm(0, 1)")

makePseudo <- function(cauchy_scale, logInd = FALSE) {
  log_pdf <- function(x) dcauchy(x, location = 0, scale = cauchy_scale, log = logInd)
  inv_cdf <- function(u) qcauchy(u, location = 0, scale = cauchy_scale)
  
    list(
      d = log_pdf,
      q = inv_cdf,
      t = paste0("Cauchy(loc=0, sc=",cauchy_scale,")")
    )
}


pseu <- lapply(scales, FUN = \(par) makePseudo(par, logInd = FALSE))

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
       col=1:(n_pseu+1), lwd=2, lty=1:(n_pseu+1), cex = 1)

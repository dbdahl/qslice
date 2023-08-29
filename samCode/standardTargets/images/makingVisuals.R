# Script to make vizuals for the paper
# author: Sam Johnson

library(tidyverse)
library(patchwork)
library(latex2exp)


filesToSource <- Sys.glob('../../utilityFunctions/*.R')
discard <- grepl('*num_of_lines*',filesToSource)
sapply(filesToSource[!discard], source)

get_data <- function(filePath) {
  temp_data <- paste0(Sys.glob(paste0(filePath,'/output/*')),'/*.csv')
  temp_data <- lapply(temp_data, read_in_data)
  df <- combine_best_func(temp_data)
  df
}


theme_set(
  theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = '',
          plot.title = element_text(hjust = 0.5, size = 30),
          axis.title.x = element_text(size = 27), 
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 15)
          # axis.line.y.left = element_line(color = 'black'),
          # axis.line.x.bottom = element_line(color = 'black'),
          # axis.ticks = element_line(color = 'black'),
          # axis.ticks.length=unit(.25, "cm"),
          # plot.title = element_text(hjust = 0.5)
          
  )
)


#### plots of the three standard cuves
pdf(file = 'plotOfStandardTargets.pdf')
par(mfrow = c(1,3))
curve(dnorm(x), -4, 4, bty = 'l', ylab = 'density', main = expression('N(0,1)'), cex = 15)
curve(dgamma(x,2.5,1), 0, 8, bty = 'l', ylab = 'density', main = expression('Ga(2.5,1)'), cex = 15)
curve(dinvgamma(x,2,1), 0, 7, bty = 'l', ylab = 'density', main = expression('IG(2,1)'), cex = 15)
dev.off()

##### 
norm_curve <- data.frame(x = seq(-4,4,0.001)) |> 
  ggplot(aes(x = x)) +
  stat_function(fun = dnorm) +
  labs(
    title = 'Norm(0,1)',
    y = '',
    x = TeX('$\\theta$')
  )

pdf(file = 'normCurve.pdf')
print(norm_curve)
dev.off()

gamma_curve <- data.frame(x = seq(0,8,0.001)) |> 
  ggplot(aes(x = x)) +
  stat_function(fun = dgamma, args = list(shape = 2.5, scale = 1)) +
  labs(
    title = 'Ga(2.5,1)',
    y = '',
    x = TeX('$\\theta$')
  )

pdf(file = 'gammaCurve.pdf')
print(gamma_curve)
dev.off()

invgamma_curve <- data.frame(x = seq(0,8,0.001)) |> 
  ggplot(aes(x = x)) +
  stat_function(fun = dinvgamma, args = list(shape = 2, scale = 1)) + 
  labs(
    title = 'IG(2,1)',
    y = '',
    x = TeX('$\\theta$')
  )

pdf(file = 'inverseGammaCurve.pdf')
print(invgamma_curve)
dev.off()

pdf(file = 'plotOfStandardTargetsggplot.pdf')
norm_curve + gamma_curve + invgamma_curve
dev.off()

#####

theme_set(
  theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = '',
          plot.title = element_text(size = 15),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 8),
          plot.margin = margin(l = -12)
    )
)

# normal
normal_data <- get_data('../normTarget/')

norm_plot <- normal_data |> 
  ggplot(aes(x = sampPsec, y = par, fill = method)) +
  geom_boxplot() + 
  labs(
    title = 'Normal Target',
    x = '',
    y = ''
  ) +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = 'K'))

# gamma
gamma_data <- get_data('../gammaTarget/')

gamma_plot <- gamma_data |> 
  ggplot(aes(x = sampPsec, y = par, fill = method)) +
  geom_boxplot() + 
  labs(
    title = 'Gamma Target',
    x = '',
    y = ''
  ) +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = 'K'))

# inverse gamma
invgamma_data <- get_data('../invGammaTarget/')

invgamma_plot <- invgamma_data |> 
  ggplot(aes(x = sampPsec, y = par, fill = method)) +
  geom_boxplot() +
  theme(
    legend.position = 'bottom',
    legend.key.size = unit(.5, 'cm'),
    axis.title.x = element_text(size = 11)
  ) + 
  labs(
    title = 'Inverse Gamma Target',
    x = 'Effective Samples Per CPU Second',
    y = '',
    fill = 'Method'
  ) +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = 'K'))


# all plots
pdf(file = 'plotOfResults.pdf')
norm_plot / gamma_plot / invgamma_plot
dev.off()


## to plot what the tranfomation looks like
xx = seq(-5, 12, length=1001)
truth = list(d = function(x) {dgamma(x, 6, sc=0.5)},
             q = function(u) {qgamma(u, 6, sc=0.5)},
             t = "Gamma(sh=6, sc=0.5)")

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
     xlab=expression(theta), ylab="", axes=FALSE, cex.lab=1.5, line=2.5)
axis(side=1)
for(i in n_pseu:1) {
  lines(xx, pseu[[i]]$d(xx), col=cols[i], lwd=5, lty=1)
}
text(0.4, 0.33, expression(f(theta)))
text(5.5, 0.2, expression(hat(pi)(theta)), col=cols[1])
text(5.8, 0.14, expression(hat(pi)(theta)), col=cols[2])
text(8.8, 0.14, expression(hat(pi)(theta)), col=cols[3])

par(mar=c(3.5,1,1,1))
uu <- seq(0,1,length.out=1001)
plot(uu, rep(1.0, length(uu)), type="l", col=1, lwd=5,
     ylim=c(0,11),
     xlab=expression(psi), ylab="",
     axes=FALSE, cex.lab=1.5, line=2.5)
axis(side=1)
mtext(expression(h(hat(Pi)^-1*(psi))), side=2, line=-1.0, las=3)

for(i in n_pseu:1) {
  lines(uu, truth$d( pseu[[i]]$q(uu) ) / pseu[[i]]$d( pseu[[i]]$q(uu) ), col=cols[i], lwd=5, lty=1)
}

legend("topright", bty="n", inset=c(0.0, 0.0),
       legend=c(truth$t, sapply(pseu, function(x) x$t)),
       col=c("black", cols), lwd=5, lty=1, cex=1.2)
dev.off()
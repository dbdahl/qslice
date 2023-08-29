# a function that finds the pseduo target based on a laplace approximation
# author: Sam Johnson


if(!exists('pseduo_t_list')) {
  
  pseudo_t_list <- function(loc, sc, degf, lb = -Inf, ub = Inf, log_p = FALSE, name = NULL) {
    
    if(!is.null(name)) {
      t <- paste0("t(loc = ", round(loc,2), ", sc = ", round(sc,2), ", degf = ", round(degf), "), ", name)
    } else {
      t <- paste0("t(loc = ", round(loc,2), ", sc = ", round(sc,2), ", degf = ", round(degf), ")")
    }
    
    plb <- pt((lb - loc)/sc, df = degf)
    pub <- pt((ub - loc)/sc, df = degf)
    normc <- pub - plb
    
    logsc <- log(sc)
    
    list(d = function(x) {dt((x - loc)/sc, df = degf) / sc},
         ld = function(x) {dt((x - loc)/sc, df = degf, log=TRUE) - logsc},
         # dld = function(x) {-2*(x-loc)/sc^2 / (1 + ((x-loc)/sc)^2)},
         q = function(u, log.p = FALSE) {qt(plb + u*normc, log.p = log.p, df = degf)*sc + loc},
         p = function(x) {(pt((x - loc)/sc, df = degf) - plb) / normc},
         t = t,
         loc = loc, sc = sc, degf = degf)
  }
  
}

# this function finds the second derivative of function f at point x
second_derivative <- function( x, h = 1e-5, f ) {
  
  num <- f(x + h) - 2*f(x) + f(x - h)
  denom <- h^2
  
  num/denom
}

lapproxt <- function(f, init, sc_adj = 1.0, lb = -Inf, ub = Inf, ...) {
  
  fit <- optim(par = init, fn = f, control = list(fnscale = -1), method = 'BFGS')
  loc <- fit$par
  hessian <- second_derivative( x = loc, h = 1e-5, f = f )
  sc <- sc_adj / sqrt(-hessian)
  out <- pseudo_t_list(loc = loc, sc = sc, degf = 1, lb = lb, ub = ub, name = 'Laplace')
  
  out
}

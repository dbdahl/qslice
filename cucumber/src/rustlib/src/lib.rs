mod registration;
use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use roxido::*;

#[roxido]
fn slice_sampler(
    current: Rval,
    target: Rval,
    w: Rval,
    max_number_of_steps: Rval,
    on_log_scale: Rval,
) -> Rval {
    let mut rng = Pcg64Mcg::from_seed(r::random_bytes::<16>());
    let uniform = Uniform::from(0.0..1.0);
    let mut u = || uniform.sample(&mut rng);
    let mut evaluation_counter = 0;
    let mut f = |x: f64| {
        evaluation_counter += 1;
        let mut pc = Pc::new();
        target
            .call1(Rval::new(x, &mut pc), &mut pc)
            .unwrap()
            .as_f64()
    };
    let w = w.as_f64();
    if w <= 0.0 {
        panic!("Width w must be strictly positive.");
    }
    let x = current.as_f64();
    let fx = f(x);
    // Step 1 (slicing)
    let y = if on_log_scale.is_true() {
        (u() * fx.exp()).ln()
    } else {
        u() * fx
    };
    // Step 2 (stepping out)
    let mut l = x - u() * w;
    let mut r = l + w;
    let max = max_number_of_steps.as_f64();
    if !Rval::is_finite(max) {
        while y < f(l) {
            l -= w
        }
        while y < f(r) {
            r += w
        }
    } else if max > 0.0 {
        let mut j = (u() * max).floor() as u32;
        let mut k = max as u32 - 1 - j;
        while j > 0 && y < f(l) {
            l -= w;
            j -= 1;
        }
        while k > 0 && y < f(r) {
            r += w;
            k -= 1;
        }
    }
    // Step 3 (shrinkage)
    loop {
        let x1 = l + u() * (r - l);
        let fx1 = f(x1);
        if y < fx1 {
            let result = Rval::new_list(2, &mut pc);
            let names = Rval::new_vector_character(2, &mut pc);
            names.set_character_element(0, "x");
            names.set_character_element(1, "nEvaluations");
            result.names_gets(names);
            result.set_list_element(0, Rval::new(x1, &mut pc));
            drop(f);
            result.set_list_element(1, Rval::new(evaluation_counter, &mut pc));
            return result;
        }
        if x1 < x {
            l = x1;
        } else {
            r = x1;
        }
    }
}

#[roxido]
fn myrnorm(n: Rval, mean: Rval, sd: Rval) -> Rval {
    unsafe {
        use rbindings::*;
        use std::convert::TryFrom;
        let (mean, sd) = (Rf_asReal(mean.0), Rf_asReal(sd.0));
        let length = isize::try_from(Rf_asInteger(n.0)).unwrap();
        let vec = Rf_protect(Rf_allocVector(REALSXP, length));
        let slice = Rval(vec).slice_double().unwrap();
        GetRNGstate();
        for x in slice {
            *x = Rf_rnorm(mean, sd);
        }
        PutRNGstate();
        Rf_unprotect(1);
        Rval(vec)
    }
}

#[roxido]
fn convolve2(a: Rval, b: Rval) -> Rval {
    let (a, xa) = a.coerce_double(&mut pc).unwrap();
    let (b, xb) = b.coerce_double(&mut pc).unwrap();
    let (ab, xab) = Rval::new_vector_double(a.len() + b.len() - 1, &mut pc);
    for xabi in xab.iter_mut() {
        *xabi = 0.0
    }
    for (i, xai) in xa.iter().enumerate() {
        for (j, xbj) in xb.iter().enumerate() {
            xab[i + j] += xai * xbj;
        }
    }
    ab
}

#[roxido]
fn zero(f: Rval, guesses: Rval, stol: Rval, rho: Rval) -> Rval {
    let slice = guesses.slice_double().unwrap();
    let (mut x0, mut x1, tol) = (slice[0], slice[1], stol.as_f64());
    if tol <= 0.0 {
        panic!("non-positive tol value");
    }
    let symbol = Rval::new_symbol("x", &mut pc);
    let feval = |x: f64| {
        let mut pc = Pc::new();
        symbol.assign(Rval::new(x, &mut pc), rho);
        f.eval(rho, &mut pc).unwrap().as_f64()
    };
    let mut f0 = feval(x0);
    if f0 == 0.0 {
        return Rval::new(x0, &mut pc);
    }
    let f1 = feval(x1);
    if f1 == 0.0 {
        return Rval::new(x1, &mut pc);
    }
    if f0 * f1 > 0.0 {
        panic!("x[0] and x[1] have the same sign");
    }
    loop {
        let xc = 0.5 * (x0 + x1);
        if (x0 - x1).abs() < tol {
            return Rval::new(xc, &mut pc);
        }
        let fc = feval(xc);
        if fc == 0.0 {
            return Rval::new(xc, &mut pc);
        }
        if f0 * fc > 0.0 {
            x0 = xc;
            f0 = fc;
        } else {
            x1 = xc;
        }
    }
}

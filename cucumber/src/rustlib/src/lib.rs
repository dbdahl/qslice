mod registration;
use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use roxido::*;

/// The slice sampler using the stepping out and shrinkage procedures
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
    // Step 1 (slicing)
    let y = {
        let u: f64 = u();
        let fx = f(x);
        if on_log_scale.is_true() {
            u.ln() + fx
        } else {
            u * fx
        }
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

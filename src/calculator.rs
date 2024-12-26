#![allow(non_snake_case)] // pH_min, pKa, Kb etc.

use std::thread;

use crate::types::AcidBase;

/// Compute the pH_min of a `sol`ution, with solvent self-ionization constant `Ki` (= Kw = 14 for water)
pub fn compute_pH(sol: Vec<AcidBase>, Ki: f64) -> f64 
{
    // Multithreading for even better performance: One thread deals with pH from 0-6.999, another deals with pH from 7-14        
    let acid_thread = thread::spawn({
        let sol = sol.clone();
        let Ki = Ki.clone();

        move || { _compute_pH_partial(true, sol, Ki) }
    });
    let base_thread = thread::spawn(move || _compute_pH_partial(false, sol, Ki));

    let lowest_acidic = acid_thread.join().unwrap();
    let lowest_basic = base_thread.join().unwrap();

    if lowest_acidic.1.abs() < lowest_basic.1.abs()
    {
        lowest_acidic.0
    }
    else 
    {
        lowest_basic.0
    }
}

/// Compute pH in a specific range (acidic or basic)
fn _compute_pH_partial(acidic_env: bool, sol: Vec<AcidBase>, Ki: f64) -> (f64, f64)
{
    let (mut pH_min, pH_max): (f64, f64) = if acidic_env
    {
        (0.0, 6.9999)
    }
    else 
    {
        (7.0, 14.0)
    };

    let mut most_accur_pH: (f64, f64) = (0.0, std::f64::INFINITY);
    
    while pH_min < pH_max 
    {
        let mut rhs = 10_f64.powf(-Ki + pH_min);

        for species in &sol 
        {
            let mut numer: f64 = 0.0; //times coeff
            let mut denom: f64 = 0.0; 

            if species.is_acidic 
            {
                /* E.g. for a triprotic acid (H3A):
                *
                * H3A = H+ + H2A(-) (Ka1)
                * H2A(-) = H+ + HA(2-) (Ka2)
                * HA(2-) = H+ + A(3-) (Ka3)
                *
                *                 numer = C0 * ([H+]^2 * Ka1 + 2[H+] * Ka1 * Ka2 + 3Ka1 * Ka2 * Ka3)
                * RHS = Kw/[H+] + __________________________________________________________________
                *                 denom = [H+]^3 + [H+]^2 * Ka1 + [H+] * Ka1 * Ka2 + Ka1 * Ka2 * Ka3
                */

                for i in 0..=species.dissoc_consts.len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_min * (species.dissoc_consts.len() as f64 - i as f64)) * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Ka| accum * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                    }

                    denom += 10_f64.powf(-pH_min * (species.dissoc_consts.len() as f64 - i as f64)) * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Ka| accum * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                }

                rhs += species.conc.clamp(0.0, std::f64::INFINITY) * numer / denom;
            } 
            else 
            {
                // Convert Ka to Kb for bases and sort
                /* E.g. for a diprotic base:
                                 numer = C0 * ([H+]Ka1 + 2 * [H+]^2)
                [H+] = Kw/[H+] - ____________________________________
                                 denom = [H+]^2 + [H+]Ka1 + Ka1 * Ka2
                */

                for i in 0..=species.dissoc_consts.len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_min * i as f64) * 
                        match species.dissoc_consts[0..(species.dissoc_consts.len() - i)]
                            .iter()
                            .copied()
                            .reduce(|accum, Kb| accum * 10_f64.powf(-Ki) / Kb)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                    }

                    denom += 10_f64.powf(-pH_min * (species.dissoc_consts.len() as f64 - i as f64)) * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Kb| accum * 10_f64.powf(-Ki) / Kb)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                }

                rhs -= species.conc.clamp(0.0, std::f64::INFINITY) * numer / denom;
            }
        }

        // Replace with more accurate pH_min value (smaller difference between LHS ([H+]) and RHS)
        if (10_f64.powf(-pH_min) - rhs).abs() < most_accur_pH.1 
        {
            most_accur_pH.0 = pH_min;
            most_accur_pH.1 = 10_f64.powf(-pH_min) - rhs;
        }

        pH_min += 0.0001;
    }

    most_accur_pH
}
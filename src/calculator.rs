#![allow(non_snake_case)] // pH_min, pKa, Kb etc.

use crate::types::AcidBase;
use std::thread;

/// Compute the pH_min of a `sol`ution, with solvent self-ionization constant `pKi` (= pKw = 14 for water)
pub fn compute_pH(sol: Vec<AcidBase>, pKi: f64) -> f64 
{
    // Multithreading for even better performance 
    // One thread deals with pH from 0-7, another deals with pH from 7-14
    // Of the 2 threads, the one with the more accurate pH value (smaller difference between RHS and LHS) will be returned

    let acid_thread = thread::spawn({
        let sol = sol.clone();
        let Ki = 10_f64.powf(-pKi);

        move || { _compute_pH_partial(true, sol, Ki) }
    });
    let base_thread = thread::spawn(move || _compute_pH_partial(false, sol, 10_f64.powf(-pKi)));

    let closest_acidic_pH = acid_thread.join().unwrap();
    let closest_basic_pH = base_thread.join().unwrap();

    if closest_acidic_pH.1.abs() < closest_basic_pH.1.abs()
    {
        closest_acidic_pH.0
    }
    else 
    {
        closest_basic_pH.0
    }
}

/// Compute pH in a specific range (acidic or basic)
fn _compute_pH_partial(acidic_env: bool, sol: Vec<AcidBase>, Ki: f64) -> (f64, f64)
{
    let (mut pH_start, pH_end): (f64, f64) = if acidic_env
    {
        (0.0, 7.0)
    }
    else 
    {
        (7.0, 14.0)
    };

    let mut most_accur_pH: (f64, f64) = (0.0, std::f64::INFINITY);
    
    while pH_start < pH_end 
    {
        let mut rhs = Ki / 10_f64.powf(-pH_start);

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

                for i in 0..=species.dissoc_consts_acid.len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_start * (species.dissoc_consts_acid.len() as f64 - i as f64)) * 
                        match species.dissoc_consts_acid[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                    }

                    denom += 10_f64.powf(-pH_start * (species.dissoc_consts_acid.len() as f64 - i as f64)) * 
                        match species.dissoc_consts_acid[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                }

                rhs += species.conc.clamp(0.0, std::f64::INFINITY) * numer / denom;
            } 
            else 
            {
                // Ka of conjugate acid
                /* E.g. for a diprotic base:
                                 numer = C0 * ([H+]Ka1 + 2 * [H+]^2)
                [H+] = Kw/[H+] - ____________________________________
                                 denom = [H+]^2 + [H+]Ka1 + Ka1 * Ka2
                */

                for i in 0..=species.dissoc_consts_acid.len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_start * i as f64) * 
                        match species.dissoc_consts_acid[0..(species.dissoc_consts_acid.len() - i)]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                    }

                    denom += 10_f64.powf(-pH_start * (species.dissoc_consts_acid.len() as f64 - i as f64)) * 
                        match species.dissoc_consts_acid[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                }

                rhs -= species.conc.clamp(0.0, std::f64::INFINITY) * numer / denom;
            }
        }

        // Replace with more accurate pH_min value (smaller difference between LHS ([H+]) and RHS)
        if (10_f64.powf(-pH_start) - rhs).abs() < most_accur_pH.1 
        {
            most_accur_pH.0 = pH_start;
            most_accur_pH.1 = 10_f64.powf(-pH_start) - rhs;
        }

        pH_start += 0.0001;
    }

    most_accur_pH
}
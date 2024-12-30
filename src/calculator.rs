#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

use crate::types::AcidBase;
use std::thread;

/// Compute the pH of a `solution`, with solvent self-ionization constant `pKi` (= pKw = 14 for water)
/// Precision is limited to maximum 5 decimal places
pub fn compute_pH(solution: &[AcidBase], pKi: f64, precision: u8) -> f64 
{
    assert!(precision > 0 && precision < 5, "Precision must range from 0 to 5 d.p.");

    // Multithreading for even better performance 
    // One thread deals with pH from 0-7, another deals with pH from 7-14
    // Of the 2 threads, the one with the more accurate pH value (smaller difference between RHS and LHS) will be returned

    let acid_thread = thread::spawn({
        let sol = solution.to_owned();
        move || compute_pH_partial(true, &sol, 10_f64.powf(-pKi), precision)
    });
    let base_thread = thread::spawn({
        let sol = solution.to_owned();
        move || compute_pH_partial(false, &sol, 10_f64.powf(-pKi), precision)
    });

    // Calculate step-by-step if multithreading fails
    let closest_acidic_pH = match acid_thread.join()
    {
        Ok(result) => result,
        Err(err) =>
        {
            eprintln!("{err:?}");
            compute_pH_partial(true, &solution, 10_f64.powf(-pKi), precision)
        }
    };
    let closest_basic_pH = match base_thread.join()
    {
        Ok(result) => result,
        Err(err) =>
        {
            eprintln!("{err:?}");
            compute_pH_partial(false, &solution, 10_f64.powf(-pKi), precision)
        }
    };

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
fn compute_pH_partial(acidic_env: bool, sol: &[AcidBase], Ki: f64, dp: u8) -> (f64, f64)
{
    let (mut pH_start, pH_end): (f64, f64) = if acidic_env
    {
        (0.0, 7.0)
    }
    else 
    {
        (7.0, 14.0)
    };

    let mut most_accur_pH: (f64, f64) = (0.0, f64::INFINITY);
    
    while pH_start < pH_end 
    {
        let mut rhs = Ki / 10_f64.powf(-pH_start);

        for species in sol 
        {
            let mut numer: f64 = 0.0; //times coeff
            let mut denom: f64 = 0.0; 

            if species.is_acidic() 
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

                for i in 0..=species.dissoc_consts_acid().len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_start * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
                        species.dissoc_consts_acid()[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                            .unwrap_or(1.0);
                    }

                    denom += 10_f64.powf(-pH_start * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
                        species.dissoc_consts_acid()[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                            .unwrap_or(1.0);
                }

                rhs += species.conc() * numer / denom;
            } 
            else 
            {
                // Ka of conjugate acid of the base
                /* E.g. for a diprotic base:
                                 numer = C0 * ([H+]Ka1 + 2 * [H+]^2)
                [H+] = Kw/[H+] - ____________________________________
                                 denom = [H+]^2 + [H+]Ka1 + Ka1 * Ka2
                */

                for i in 0..=species.dissoc_consts_acid().len() 
                {
                    if i > 0
                    {
                        numer += i as f64 * 10_f64.powf(-pH_start * i as f64) * 
                        species.dissoc_consts_acid()[0..(species.dissoc_consts_acid().len() - i)]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                            .unwrap_or(1.0);
                    }

                    denom += 10_f64.powf(-pH_start * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
                        species.dissoc_consts_acid()[0..i]
                            .iter()
                            .copied()
                            .reduce(|accumulator, Ka| accumulator * Ka)
                            .unwrap_or(1.0);
                }

                rhs -= species.conc() * numer / denom;
            }
        }

        // Replace with more accurate pH_min value (smaller difference between LHS ([H+]) and RHS)
        if (10_f64.powf(-pH_start) - rhs).abs() < most_accur_pH.1 
        {
            most_accur_pH.0 = pH_start;
            most_accur_pH.1 = 10_f64.powf(-pH_start) - rhs;
        }

        pH_start += 10_f64.powi(-(dp as i32));
    }

    most_accur_pH
}
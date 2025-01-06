#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

use crate::types::{
    AcidBase, 
    SolProperties
};

/// Compute the pH of a `solution`, with solvent `properties`
pub fn compute_pH(solution: &[AcidBase], properties: &SolProperties) -> f64 
{
    let mut lower_bound_pH: f64 = properties.min_pH();
    let mut upper_bound_pH: f64 = properties.max_pH();
    let mut mean_pH: f64 = (lower_bound_pH + upper_bound_pH) / 2.0;

    for _ in 0..1_000
    {
        let left_value = calculate_diff(solution, 10_f64.powf(-properties.pKi()), lower_bound_pH);
        let right_value = calculate_diff(solution, 10_f64.powf(-properties.pKi()), upper_bound_pH);
        let center_value = calculate_diff(solution, 10_f64.powf(-properties.pKi()), mean_pH);

        if center_value == 0.0
        {
            return mean_pH;
        }
        else if right_value > 0.0 && left_value > 0.0
        {
            return properties.max_pH();
        }
        else if left_value < 0.0 && right_value < 0.0
        {
            return properties.min_pH();
        }
        else if center_value * left_value >= 0.0
        {
            lower_bound_pH = mean_pH;
        }
        else if center_value * right_value >= 0.0
        {
            upper_bound_pH = mean_pH;
        }
        else 
        {
            panic!("Calculation error with lower_bound_pH = {lower_bound_pH}, upper_bound_pH = {upper_bound_pH}");
        }

        mean_pH = (lower_bound_pH + upper_bound_pH) / 2.0;
    }

    mean_pH
}

/// Calculate: LHS (= [H+]) - RHS
fn calculate_diff(sol: &[AcidBase], Ki: f64, pH: f64) -> f64
{
    let mut rhs = Ki / 10_f64.powf(-pH);

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
                    numer += i as f64 * 10_f64.powf(-pH * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
                    species.dissoc_consts_acid()[0..i]
                        .iter()
                        .copied()
                        .reduce(|accumulator, Ka| accumulator * Ka)
                        .unwrap_or(1.0);
                }

                denom += 10_f64.powf(-pH * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
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
                    numer += i as f64 * 10_f64.powf(-pH * i as f64) * 
                    species.dissoc_consts_acid()[0..(species.dissoc_consts_acid().len() - i)]
                        .iter()
                        .copied()
                        .reduce(|accumulator, Ka| accumulator * Ka)
                        .unwrap_or(1.0);
                }

                denom += 10_f64.powf(-pH * (species.dissoc_consts_acid().len() as f64 - i as f64)) * 
                    species.dissoc_consts_acid()[0..i]
                        .iter()
                        .copied()
                        .reduce(|accumulator, Ka| accumulator * Ka)
                        .unwrap_or(1.0);
            }

            rhs -= species.conc() * numer / denom;
        }
    }

    10_f64.powf(-pH) - rhs
}
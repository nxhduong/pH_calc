#![allow(non_snake_case)]

use crate::types::AcidBase;

pub fn calculate_pH(sol: Vec<AcidBase>, kw: f64) -> f64 
{
    let mut most_accur_pH: (f64, f64) = (0.0, std::f64::INFINITY);
    let mut pH = 0.0;
    
    while pH < 14.0 
    {
        let mut rhs = 10_f64.powf(-kw + pH);

        for species in &sol 
        {
            let mut numer: f64 = 0.0; //times coeff
            let mut denom: f64 = 0.0; 

            if species.is_acid 
            {
                for i in 1..species.dissoc_consts.len()
                {
                    // TODO: numer += species
                }

                for i in 0..species.dissoc_consts.len() 
                {
                    denom += 10_f64.powf(pH).powf(species.dissoc_consts.len() as f64 - i as f64) * 
                    species.dissoc_consts[0..i]
                    .iter()
                    .copied()
                    .reduce(|accum, Ka| accum * Ka)
                    .unwrap();
                }

                rhs += species.conc * numer / denom;
            } 
            else 
            {
                for i in 0..species.dissoc_consts.len() 
                {
                    denom += 10_f64.powf(pH).powf(species.dissoc_consts.len() as f64 - i as f64) * 
                    species.dissoc_consts[0..i]
                    .iter()
                    .copied()
                    .reduce(|accum, Kb| accum * 10_f64.powf(-14.0) / Kb)
                    .unwrap();
                }
                rhs -= species.conc * numer / denom;
            }
        }

        if 10_f64.powf(-pH) - rhs < most_accur_pH.1 
        {
            most_accur_pH.0 = pH;
            most_accur_pH.1 = 10_f64.powf(-pH) - rhs;
        }

        pH += 0.0001;
    }

    pH
}
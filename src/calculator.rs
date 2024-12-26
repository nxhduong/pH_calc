#![allow(non_snake_case)] // pH, pKa, Kb etc.

use crate::types::AcidBase;

pub fn compute_pH(sol: Vec<AcidBase>, Kw: f64) -> f64 
{
    let mut most_accur_pH: (f64, f64) = (0.0, std::f64::INFINITY);
    let mut pH = 0.0;
    
    while pH < 14.0 
    {
        let mut rhs = 10_f64.powf(-Kw + pH);

        for species in &sol 
        {
            let mut numer: f64 = 0.0; //times coeff
            let mut denom: f64 = 0.0; 

            if species.is_acid 
            {
                for i in 1..species.dissoc_consts.len()
                {
                    numer += i as f64 * 10_f64.powf(-pH).powf(species.dissoc_consts.len() as f64 - i as f64) * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Ka| accum * Ka)
                        {
                            Some(product) => product,
                            None => 1.0
                        };
                }

                for i in 0..=species.dissoc_consts.len() 
                {
                    denom += 10_f64.powf(-pH).powf(species.dissoc_consts.len() as f64 - i as f64) * 
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
                /* Convert Ka to Kb for bases
                for i in 0..=species.dissoc_consts.len() 
                {
                    denom += 10_f64.powf(pH).powf(species.dissoc_consts.len() as f64 - i as f64) * 
                        species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Kb| accum * 10_f64.powf(-Kw) / Kb)
                            .unwrap();
                }
                rhs -= species.conc.clamp(0.0, std::f64::INFINITY) * numer / denom;*/
            }
        }

        // Replace with more accurate pH value (smaller difference between LHS ([H+]) and RHS)
        if (10_f64.powf(-pH) - rhs).abs() < most_accur_pH.1 
        {
            most_accur_pH.0 = pH;
            most_accur_pH.1 = 10_f64.powf(-pH) - rhs;
        }

        pH += 0.001;
    }

    most_accur_pH.0
}
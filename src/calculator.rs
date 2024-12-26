#![allow(non_snake_case)] // pH, pKa, Kb etc.

use crate::types::AcidBase;
use decimal_rs::*;

/// Compute the pH of a `sol`ution, with solvent autodissociation constant `Ki` (= Kw = 14 for water)
pub fn compute_pH(sol: Vec<AcidBase>, Ki: Decimal) -> Decimal 
{
    let INF = Decimal::from(50);
    let TEN = Decimal::from(10);
    let mut most_accur_pH: (Decimal, Decimal) = (Decimal::ZERO, INF);
    let mut pH = Decimal::ZERO;
    
    while pH < Decimal::from(14)
    {
        let mut rhs = TEN.checked_pow(&(-Ki + pH)).unwrap();

        for species in &sol 
        {
            let mut numer = Decimal::ZERO;
            let mut denom = Decimal::ZERO; 

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

                for i in 1..=species.dissoc_consts.len()
                {
                    numer += Decimal::from(i) * 
                        TEN.checked_pow(&(-pH * (Decimal::from(species.dissoc_consts.len()) - Decimal::from(i)))).unwrap() * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Ka| accum * Ka)
                        {
                            Some(product) => product,
                            None => Decimal::ONE
                        };
                }

                for i in 0..=species.dissoc_consts.len() 
                {
                    denom += TEN.checked_pow(&(-pH * (Decimal::from(species.dissoc_consts.len()) - Decimal::from(i)))).unwrap() * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Ka| accum * Ka)
                        {
                            Some(product) => product,
                            None => Decimal::ONE
                        };
                }

                rhs += species.conc.clamp(Decimal::ZERO, INF) * numer / denom;
            } 
            else 
            {
                // Convert Ka to Kb for bases and sort
                /* E.g. for a diprotic base:
                                 numer = C0 * ([H+]Ka1 + 2 * [H+]^2)
                [H+] = Kw/[H+] - ____________________________________
                                 denom = [H+]^2 + [H+]Ka1 + Ka1 * Ka2
                */

                for i in 1..=species.dissoc_consts.len()
                {
                    numer += Decimal::from(i) * TEN.checked_pow(&(-pH * Decimal::from(i))).unwrap() * 
                        match species.dissoc_consts[0..(species.dissoc_consts.len() - i)]
                            .iter()
                            .copied()
                            .reduce(|accum, Kb| accum * TEN.checked_pow(&-Ki).unwrap() / Kb)
                        {
                            Some(product) => product,
                            None => Decimal::ONE
                        };
                }

                for i in 0..=species.dissoc_consts.len() 
                {
                    denom += TEN.checked_pow(&(-pH * (Decimal::from(species.dissoc_consts.len()) - Decimal::from(i)))).unwrap() * 
                        match species.dissoc_consts[0..i]
                            .iter()
                            .copied()
                            .reduce(|accum, Kb| accum * TEN.checked_pow(&-Ki).unwrap() / Kb)
                        {
                            Some(product) => product,
                            None => Decimal::ONE
                        };
                }
                
                rhs -= species.conc.clamp(Decimal::ZERO, INF) * numer / denom;
            }
        }

        // Replace with more accurate pH value (smaller difference between LHS ([H+]) and RHS)
        if (TEN.checked_pow(&-pH).unwrap() - rhs).abs() < most_accur_pH.1 
        {
            most_accur_pH.0 = pH;
            most_accur_pH.1 = TEN.checked_pow(&-pH).unwrap() - rhs;
        }

        pH += "0.001".parse::<Decimal>().unwrap();
    }

    most_accur_pH.0
}
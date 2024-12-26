#![allow(non_snake_case)]

/// Represents an acid or base in a solution
/// `is_acidic` is true if the species is an acid, false if it is a base
/// `conc` is the concentration of the species in the solution
/// `dissoc_consts` is a vector of dissociation constants for the species (Ka for acids, Kb for bases)
#[derive(Clone)]
pub struct AcidBase 
{
    pub is_acidic: bool,
    pub conc: f64,
    pub dissoc_consts: Vec<f64>
}

impl AcidBase 
{
    // Recommended way of instantiating new species
    pub fn new(is_acid: bool, conc: f64, mut pK_values: Vec<f64>) -> Self 
    {
        if conc < 0.0 
        {
            panic!("Concentration must not be negative")
        }

        // Sort x in 10^x (not sorting 10^x)
        pK_values.sort_by(|val1, val2| val1.partial_cmp(val2).unwrap());

        AcidBase 
        { 
            is_acidic: is_acid, 
            conc, 
            dissoc_consts: pK_values
                .iter()
                .map(|pK| 10_f64.powf(-pK))
                .collect() 
        }
    }
}
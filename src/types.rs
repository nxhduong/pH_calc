#![allow(non_snake_case)]

pub struct AcidBase 
{
    pub is_acid: bool,
    pub conc: f64,
    pub dissoc_consts: Vec<f64>
}

impl AcidBase 
{
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
            is_acid, 
            conc, 
            dissoc_consts: pK_values
                .iter()
                .map(|pK| 10_f64.powf(-pK))
                .collect() 
        }
    }
}
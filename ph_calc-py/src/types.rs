#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

use pyo3::pyclass;

/// Represents an acid or base in a solution
#[pyclass]
#[derive(Clone)]
pub struct AcidBase
{
    /// `is_acid` is true if the species is an acid, false if it is a base
    is_acid: bool,

    /// `conc` is the concentration of the species in the solution
    conc: f64,

    /// `dissoc_consts_acid` is a collection of dissociation constants (Ka) for the acid (or conjugate acid of the base)
    dissoc_consts_acid: Vec<f64>
}

impl AcidBase
{
    /// Recommended way of instantiating new acids/bases
    /// - `is_acidic` is true if the species is an acid, false if it is a base
    /// - `conc` is the concentration of the species in the solution
    /// - `pKa_values` is a collection of dissociation constants (pKa) for the acid (or conjugate acid of the base)
    /// # Panics
    /// - Concentration < 0
    /// - No pKa supplied
    pub fn new(is_acidic: bool, conc: f64, pKa_values: &mut [f64]) -> Self 
    {
        // Prevent dissoc_consts_acid[0..0] and unrealistic numbers
        assert!(
            !(conc < 0.0 || pKa_values.is_empty()), 
            "Concentration must not be negative and pKa must be supplied (For strong acids, pKa = 3 should be enough)."
        );

        // Sort x in 10^x (not sorting 10^x)
        pKa_values.sort_by(|val1, val2| val1.partial_cmp(val2).unwrap());

        Self 
        {
            is_acid: is_acidic, 
            conc, 
            dissoc_consts_acid: pKa_values
                .iter()
                .map(|pK| 10_f64.powf(-pK))
                .collect() 
        }
    }

    pub const fn is_acidic(&self) -> bool
    {
        self.is_acid
    }

    pub const fn conc(&self) -> f64
    {
        self.conc
    }

    pub const fn dissoc_consts_acid(&self) -> &Vec<f64>
    {
        &self.dissoc_consts_acid
    }
}

/// Properties of the solvent
#[pyclass]
#[derive(Clone)]
pub struct SolProperties
{
    /// Self-ionization constant of solvent (`14` for water)
    pKi: f64,

    /// Minimum possible pH value of the solution
    min_pH: f64,

    /// Maximum possible pH value of the solution
    max_pH: f64
}

impl SolProperties
{
    /// Default properties of water
    pub const fn default_water() -> Self
    {
        Self
        {
            pKi: 14.0,
            min_pH: -2.0,
            max_pH: 16.0
        }
    }

    /// Instantiates new `SolProperties`
    /// # Panics 
    /// if `min_pH` is not lower than `max_pH`
    pub fn new(pKi: f64, min_pH: f64, max_pH: f64) -> Self
    {
        assert!(min_pH < max_pH, "Minimum pH must be lower than maximum pH");

        Self
        {
            pKi, // pKw
            min_pH,
            max_pH
        }
    }

    pub const fn pKi(&self) -> f64
    {
        self.pKi
    }

    pub const fn min_pH(&self) -> f64
    {
        self.min_pH
    }

    pub const fn max_pH(&self) -> f64
    {
        self.max_pH
    }
}
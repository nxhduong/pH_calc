#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

#![warn(
    clippy::all,
    clippy::restriction,
    clippy::pedantic,
    clippy::nursery,
    clippy::cargo,
)]

pub mod calculator;
pub mod types;

#[cfg(feature = "pyo3")]
use {
    calculator::compute_pH,
    types::*,
    pyo3::{
        pyfunction,
        pymodule,
        types::{
            PyModule,
            PyModuleMethods
        },
        PyResult,
        Bound,
        wrap_pyfunction
    }
};

/// Python version
/// Compute the pH of a `solution`, with solvent `properties`
#[cfg(feature = "pyo3")]
#[cfg_attr(feature = "pyo3", pyfunction)]
pub fn calculate_pH(solution: Vec<AcidBase>, properties: SolProperties) -> f64 
{
    compute_pH(solution.as_slice(), &properties)
}

/// Python module
#[cfg(feature = "pyo3")]
#[cfg_attr(feature = "pyo3", pymodule)]
fn pH_calc(py_mod: &Bound<'_, PyModule>) -> PyResult<()> {
    py_mod.add_function(wrap_pyfunction!(calculate_pH, py_mod)?)?;
    py_mod.add_class::<AcidBase>()?;
    py_mod.add_class::<SolProperties>()?;
    Ok(())
}
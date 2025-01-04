#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

pub mod calculator;
pub mod types;

use pyo3::{
    *,
    types::{
        PyModule,
        PyModuleMethods
    }
};
use calculator::compute_pH;
use types::*;

/// Compute the pH of a `solution`, with solvent `properties`
#[pyfunction]
pub fn calculate_pH(solution: Vec<AcidBase>, properties: SolProperties) -> f64 
{
    compute_pH(solution.as_slice(), &properties)
}

/// Python module
#[pymodule]
fn pH_calc(py_mod: &Bound<'_, PyModule>) -> PyResult<()> {
    py_mod.add_function(wrap_pyfunction!(calculate_pH, py_mod)?)?;
    Ok(())
}
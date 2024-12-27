# pH_calc
A simple library crate written in Rust that can calculate the pH of a solution from the pKa(b)s and concentrations of species in that solution.
## How it works
(WIP)
## Usage
```rust 
pub fn compute_pH(solution: &[AcidBase], pKi: f64) -> f64
```
Where:
- `solution`: Species in a solution
    - `AcidBase`: comprises the following properties:
        ```rust
        is_acid: bool, // Whether it is an acid or a base
        conc: f64, // Concentration
        dissoc_consts_acid: Vec<f64> // Ka of the acid/conjugate acid of the base
        ```
        To instantiate new acids and bases, use the contructor:
        ```rust
        pub fn new(is_acidic: bool, conc: f64, mut pKa_values: Vec<f64>) -> Self
        ```
- `pKi`: Self-ionization constant of solvent (`14` for water)

Returns:
- pH (set to 4 decimal places)

Example code snippet:
```rust
// Calculate pH of 0.1M phosphoric acid in water
compute_pH(&[
        AcidBase::new(
            true,
            0.1,
            &mut [2.12, 7.21, 12.67]
        )
    ], 14.0)
// Output: ~1.62
```
For amphoteric species, treat them as separate acids and bases.
## TODOs
- Include activity of ions and volume in calculation
- Include automatic detection of amphoteric species.
## Contribution
All contributions are welcomed.
## License
Please see `LICENSE` for more information.
## Contact
- My GitHub: [github.com/nxhduong](https://github.com/nxhduong)
- My email: duong70g@gmail.com.
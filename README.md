# pH_calc
A simple library crate written in Rust that can calculate the pH of a solution, given the pKa(b)s and concentrations of species in that solution.
## How it works

## How to use
```rust 
fn compute_pH(sol: Vec<AcidBase>, Ki: f64) -> f64
```
Where:
- `sol`: Species in a solution
    - `AcidBase`: comprises the following properties:
        ```rust
        pub is_acidic: bool, // Whether it is an acid or a base
        pub conc: f64, // Concentration
        pub dissoc_consts: Vec<f64> // Ka for acids, Kb for bases
        ```
        It is recommended to use the contructor
        ```rust
        pub fn new(is_acid: bool, conc: f64, mut pK_values: Vec<f64>) -> Self
        ```
        which takes equilibrium constants as pKa/pKb and does range check on concentration.
- `Ki`: Self-ionization constant of solvent (`14` for water)

Returns:
- pH (set to 4 decimal places)

Example snippet:
```rust
// Calculate pH of 0.1M phosphoric acid in water
compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0)
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
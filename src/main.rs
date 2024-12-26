pub mod calculator;
pub mod types;

use calculator::compute_pH;
use types::AcidBase;

fn main() 
{
    println!("{}", compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0));
}
pub mod calculator;
pub mod types;

use std::time::Instant;
use calculator::compute_pH;
use types::AcidBase;

fn main() 
{
    // Example usage
    let now = Instant::now();
    println!("{}", compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0));
    println!("{:2?}", now.elapsed());
}
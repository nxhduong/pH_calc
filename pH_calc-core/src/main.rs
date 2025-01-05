use std::time::Instant;
use types::{AcidBase, SolProperties};
use calculator::compute_pH;

pub mod calculator;
pub mod types;

fn main()
{
    let now = Instant::now();

    println!("{}", compute_pH(&[
        AcidBase::new(
            true,
            0.02,
            &mut [2.12, 7.21, 12.67]
        ),
        AcidBase::new(
            true,
            0.01,
            &mut [2.0, 2.7, 6.16, 10.26]
        ),
        AcidBase::new(
            false,
            0.25,
            &mut [9.24, 12.4, 13.3]
        ),
        AcidBase::new(
            true,
            0.001,
            &mut [-3.0, 1.99]
        ),
        AcidBase::new(
            false,
            0.05,
            &mut [4.76]
        )
    ], &SolProperties::default_water()));

    println!("{:2?}", now.elapsed());
}
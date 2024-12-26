use calculator::compute_pH;
use decimal_rs::*;
use types::AcidBase;
use std::time::Instant;

pub mod calculator;
pub mod types;

fn main() 
{
    let now = Instant::now();
    println!("{}", compute_pH(vec![
        AcidBase::new(
            true,
            "0.1".parse::<Decimal>().unwrap(),
            vec!["2.12".parse::<Decimal>().unwrap(), "7.21".parse::<Decimal>().unwrap(), "12.67".parse::<Decimal>().unwrap()]
        )
    ], Decimal::from(14)));
    println!("Elapsed: {:.2?}", now.elapsed());
}
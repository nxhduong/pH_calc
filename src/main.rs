pub mod calculator;
pub mod types;

use calculator::calculate_pH;
use types::AcidBase;

fn main() 
{
    debug_assert_eq!(calculate_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0).round(), 4.3);
}
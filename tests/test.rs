use pH_calc::{
    calculator::compute_pH, 
    types::AcidBase
};
use std::time::Instant;

// Show output (execution time) by running cargo test --nocapture

#[test]
fn triprotic_acid()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0)).round(), 16.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn sulfuric_acid()
{
    let now = Instant::now();

    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![-3.0, 1.99]
        )
    ], 14.0).round(), 1.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn monoprotic_acid()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(vec![
        AcidBase::new(
            true,
            0.02,
            vec![4.76]
        )
    ], 14.0)).round(), 32.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn simple_buffer()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![4.21]
        ),
        AcidBase::new(
            false,
            0.1,
            vec![9.79]
        )
    ], 14.0)).round(), 42.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn monoprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(vec![
        AcidBase::new(
            false,
            0.03,
            vec![4.76]
        )
    ], 14.0).round(), 109.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn diprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(vec![
        AcidBase::new(
            false,
            0.25,
            vec![7.65, 3.67]
        )
    ], 14.0).round(), 119.0);

    println!("{:2?}", now.elapsed());
}

#[test]
#[ignore = "test not completed"]
fn mc_ilvaine_buffer() 
{
    let now = Instant::now();

    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![3.86]
        ),
        AcidBase::new(
            true,
            0.1,
            vec![6.0]
        ),
        AcidBase::new(
            true,
            0.1,
            vec![8.2]
        ),
        AcidBase::new(
            false,
            0.1,
            vec![10.4]
        ),
        AcidBase::new(
            false,
            0.1,
            vec![12.5]
        )
    ], 14.0), 6.0);

    println!("{:2?}", now.elapsed());
}

#[test]
fn super_solution()
{
    let now = Instant::now();

    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.02,
            vec![2.12, 7.21, 12.67]
        ),
        AcidBase::new(
            true,
            0.01,
            vec![2.0, 2.7, 6.16, 10.26]
        ),
        AcidBase::new(
            false,
            0.25,
            vec![4.76, 1.6, 0.7]
        ),
        AcidBase::new(
            true,
            0.001,
            vec![-3.0, 1.99]
        ),
        AcidBase::new(
            false,
            0.05,
            vec![9.24]
        )
    ], 14.0), 7.0);

    println!("{:2?}", now.elapsed());
}
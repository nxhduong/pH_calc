use pH_calc::{
    calculator::compute_pH, 
    types::AcidBase
};
use std::time::Instant;

// Show output (execution time) by running cargo test -- --nocapture

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

    println!("Triprotic acid: {:2?}", now.elapsed());
}

#[test]
fn diprotic_acid()
{
    let now = Instant::now();

    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![-3.0, 1.99]
        )
    ], 14.0).round(), 1.0);

    println!("Diprotic acid: {:2?}", now.elapsed());
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

    println!(" Monoprotic acid: {:2?}", now.elapsed());
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
            vec![4.21]
        )
    ], 14.0)).round(), 42.0);

    println!("Simple buffer: {:2?}", now.elapsed());
}

#[test]
fn monoprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(vec![
        AcidBase::new(
            false,
            0.03,
            vec![9.24]
        )
    ], 14.0).round(), 110.0);

    println!("Monoprotic base: {:2?}", now.elapsed());
}

#[test]
fn diprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(vec![
        AcidBase::new(
            false,
            0.25,
            vec![6.35, 10.33]
        )
    ], 14.0).round(), 120.0);

    println!("Diprotic base: {:2?}", now.elapsed());
}

#[test]
#[ignore = "test not completed"] // TODO
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

    assert_eq!((compute_pH(vec![
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
            vec![9.24, 12.4, 13.3]
        ),
        AcidBase::new(
            true,
            0.001,
            vec![-3.0, 1.99]
        ),
        AcidBase::new(
            false,
            0.05,
            vec![4.76]
        )
    ], 14.0)).round(), 13.0);

    println!("Super solution: {:2?}", now.elapsed());
}
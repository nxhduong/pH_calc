use pH_calc::{
    calculator::compute_pH, 
    types::AcidBase
};

// Due to technical limitations, pH values are rounded to nearest integers

#[test]
fn triprotic_acid()
{
    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0).round(), 2.0);
}

#[test]
fn sulfuric_acid()
{
    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![-3.0, 1.99]
        )
    ], 14.0), 1.0);
}

#[test]
fn monoprotic_acid()
{
    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![4.76]
        )
    ], 14.0), 3.0);
}

#[test]
fn simple_buffer()
{
    assert_eq!(compute_pH(vec![
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
    ], 14.0), 4.0);
}

#[test]
#[ignore = "not completed"]
fn mc_ilvaine_buffer() 
{
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
}

#[test]
fn super_solution()
{
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
    ], 14.0), 7.0);
}
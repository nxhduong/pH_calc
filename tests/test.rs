use pH_calc::*;

#[test]
fn triprotic_acid()
{
    assert_eq!(compute_pH(vec![
        AcidBase::new(
            true,
            0.1,
            vec![2.12, 7.21, 12.67]
        )
    ], 14.0), 1.6);
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
    ], 14.0), 2.9);
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
    ], 14.0), 4.3);
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
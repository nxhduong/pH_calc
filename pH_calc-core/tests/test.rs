use pH_calc_core::{
    calculator::compute_pH, 
    types::{
        AcidBase,
        SolProperties
    }
};
use std::time::Instant;

// Note: tests are slower than running the function in main
// Show output (execution time) by running cargo test -- --nocapture

#[test]
fn water()
{
    let now = Instant::now();

    assert_eq!(compute_pH(&[], &SolProperties::default_water()), 7.0);

    println!("Water: {:2?}", now.elapsed());
}

#[test]
fn strong_acid()
{
    let now = Instant::now();

    assert_eq!((100.0 * compute_pH(&[AcidBase::new(
        true,
        18.0,
        &mut [-3.0, 1.99]
    )], &SolProperties::default_water())).round(), -125.0);

    println!("Strong acid: {:2?}", now.elapsed());
}

#[test]
fn strong_base()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(&[AcidBase::new(
        false,
        84.0,
        &mut [13.5]
    )], &SolProperties::default_water())).round(), 147.0);

    println!("Strong base: {:2?}", now.elapsed());
}

#[test]
fn triprotic_acid()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(&[
        AcidBase::new(
            true,
            0.1,
            &mut [2.12, 7.21, 12.67]
        )
    ], &SolProperties::default_water())).round(), 16.0);

    println!("Triprotic acid: {:2?}", now.elapsed());
}

#[test]
fn diprotic_acid()
{
    let now = Instant::now();

    assert_eq!(compute_pH(&[
        AcidBase::new(
            true,
            0.1,
            &mut [-3.0, 1.99]
        )
    ], &SolProperties::default_water()).round(), 1.0);

    println!("Diprotic acid: {:2?}", now.elapsed());
}

#[test]
fn monoprotic_acid()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(&[
        AcidBase::new(
            true,
            0.02,
            &mut [4.76]
        )
    ], &SolProperties::default_water())).round(), 32.0);

    println!(" Monoprotic acid: {:2?}", now.elapsed());
}

#[test]
fn simple_buffer()
{
    let now = Instant::now();

    assert_eq!((10.0 * compute_pH(&[
        AcidBase::new(
            true,
            0.1,
            &mut [4.21]
        ),
        AcidBase::new(
            false,
            0.1,
            &mut [4.21]
        )
    ], &SolProperties::default_water())).round(), 42.0);

    println!("Simple buffer: {:2?}", now.elapsed());
}

#[test]
fn monoprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(&[
        AcidBase::new(
            false,
            0.03,
            &mut [9.24]
        )
    ], &SolProperties::default_water()).round(), 110.0);

    println!("Monoprotic base: {:2?}", now.elapsed());
}

#[test]
fn diprotic_base()
{
    let now = Instant::now();

    assert_eq!(10.0 * compute_pH(&[
        AcidBase::new(
            false,
            0.25,
            &mut [6.35, 10.33]
        )
    ], &SolProperties::default_water()).round(), 120.0);

    println!("Diprotic base: {:2?}", now.elapsed());
}

#[test]
fn mc_ilvaine_buffer() 
{
    let mixing_table = [
        (22.0, 0.4, 19.6),
        (24.0, 1.24, 18.76),
        (26.0, 2.18, 17.82),
        (28.0, 3.17, 16.83),
        (30.0, 4.11, 15.89),
        (32.0, 4.94, 15.06),
        (34.0, 5.7, 14.3),
        (36.0, 6.44, 13.56),
        (38.0, 7.1, 12.9),
        (40.0, 7.71, 12.29),
        (42.0, 8.28, 11.72),
        (44.0, 8.82, 11.18),
        (46.0, 9.35, 10.65),
        (48.0, 9.86, 10.14),
        (50.0, 10.3, 9.7),
        (52.0, 10.72, 9.28),
        (54.0, 11.15, 8.85),
        (56.0, 11.6, 8.4),
        (58.0, 12.09, 7.91),
        (60.0, 12.63, 7.37),
        (62.0, 13.22, 6.78),
        (64.0, 13.85, 6.15),
        (66.0, 14.55, 5.45),
        (68.0, 15.45, 4.55),
        (70.0, 16.47, 3.53),
        (72.0, 17.39, 2.61),
        (74.0, 18.17, 1.83),
        (76.0, 18.73, 1.27),
        (78.0, 19.15, 0.85),
        (80.0, 19.45, 0.55)        
    ];
    let now = Instant::now();

    for (res, vol_phosphate, vol_citric) in mixing_table
    {
        assert_eq!((10.0 * compute_pH(&[
            AcidBase::new(
                true,
                0.2 * (vol_phosphate) / (vol_citric + vol_phosphate),
                &mut [12.67]
            ),
            AcidBase::new(
                false,
                0.2 * (vol_phosphate) / (vol_citric + vol_phosphate),
                &mut [11.88, 6.79]
            ),
            AcidBase::new(
                true,
                0.1 * (vol_citric) / (vol_citric + vol_phosphate),
                &mut [3.13, 4.76, 6.39, 14.4]
            ),
        ], &SolProperties::default_water())).round(), res);
    }

    println!("{:2?}", now.elapsed());
}

#[test]
fn super_solution()
{
    let now = Instant::now();

    assert_eq!(compute_pH(&[
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
    ], &SolProperties::default_water()).round(), 13.0);

    println!("Super solution: {:2?}", now.elapsed());
}
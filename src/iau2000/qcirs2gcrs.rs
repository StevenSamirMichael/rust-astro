use nalgebra as na;

use super::ierstable::IERSTable;
use crate::astrotime::{AstroTime, Scale};
use crate::frametransform::{qroty, qrotz};

type Quat = na::UnitQuaternion<f64>;
type Delaunay = na::SVector<f64, 14>;

use std::f64::consts::PI;

// Read-only tables used for CIRS to GCRS rotation...
lazy_static::lazy_static! {
    static ref TABLE5A: IERSTable =
        IERSTable::from_file("tab5.2a.txt").unwrap();
    static ref TABLE5B: IERSTable =
        IERSTable::from_file("tab5.2b.txt").unwrap();
    static ref TABLE5D: IERSTable =
        IERSTable::from_file("tab5.2d.txt").unwrap();
}

pub fn qcirs2gcrs(tm: &AstroTime) -> Quat {
    let t_tt = (tm.to_mjd(Scale::TT) - 51544.5) / 36525.0;
    const ASEC2RAD: f64 = PI / 180.0 / 3600.0;

    let mut delaunay = Delaunay::zeros();

    // Arguments for lunisolar nutation
    // Equation 5.43 in IERS technical note 36

    // Mean anomaly of the Moon
    delaunay[0] = ASEC2RAD
        * (3600.0 * 134.96340251
            + t_tt * (1717915923.2178 + t_tt * (31.8792 + t_tt * (0.051635 + t_tt * -0.00024470))));

    // Mean anomaly of the sun
    delaunay[1] = ASEC2RAD
        * (3600.0 * 357.52910918
            + t_tt * (129596581.0481 + t_tt * (-0.5532 + t_tt * (0.000136 + t_tt * -0.00001149))));

    // F = L-Omega
    delaunay[2] = ASEC2RAD
        * (3600.0 * 93.27209062
            + t_tt
                * (1739527262.8478 + t_tt * (-12.7512 + t_tt * (-0.001037 + t_tt * 0.00000417))));

    // D = Mean elongation of the Moon from the Sun
    delaunay[3] = ASEC2RAD
        * (3600.0 * 297.85019547
            + t_tt
                * (1602961601.2090 + t_tt * (-6.37006 + t_tt * (0.006593 + t_tt * -0.00003169))));

    // Omega = mean longitude of ascending node of the Moon
    delaunay[4] = ASEC2RAD
        * (3600.0 * 125.04455501
            + t_tt * (-6962890.5431 + t_tt * (7.4722 + t_tt * (0.007702 + t_tt * -0.00005939))));

    // Planetary nutation
    // Equations 5.44 in IERS technical note 36
    delaunay[5] = 4.402608842 + 2608.7903141574 * t_tt;
    delaunay[6] = 3.176146697 + 1021.3285546211 * t_tt;
    delaunay[7] = 1.753470314 + 628.3075849991 * t_tt;
    delaunay[8] = 6.203480913 + 334.0612426700 * t_tt;
    delaunay[9] = 0.599546497 + 52.9690962641 * t_tt;
    delaunay[10] = 0.874016757 + 21.3299104960 * t_tt;
    delaunay[11] = 5.481293872 + 7.4781598567 * t_tt;
    delaunay[12] = 5.311886287 + 3.8133035638 * t_tt;
    delaunay[13] = t_tt * (0.02438175 + t_tt * 0.00000538691);

    // Polynomial part of X & Y, values in arcseconds
    // Equations 5.16 in IERS technical note 36
    let x0 = -0.016617
        + t_tt
            * (2004.191898
                + t_tt
                    * (-0.4297829
                        + t_tt * (-0.19861834 + t_tt * (0.000007578 + t_tt * 0.0000059285))));
    let y0 = -0.006951
        + t_tt
            * (-0.025896
                + t_tt
                    * (-22.4072747
                        + t_tt * (0.00190059 + t_tt * (0.001112526 + t_tt * 0.0000001358))));

    // Polynomial part of CIO locator s, values in microarcseconds
    // Described in table 5.2d of IERS technical note 36
    let s0 = 94.0
        + t_tt * (3808.65 + t_tt * (-122.68 + t_tt * (-72574.11 + t_tt * (27.98 + t_tt * 15.62))));

    let xsums = TABLE5A.compute(t_tt, &delaunay);
    let ysums = TABLE5B.compute(t_tt, &delaunay);
    let ssums = TABLE5D.compute(t_tt, &delaunay);
    let x = (x0 + xsums * 1.0e-6) * ASEC2RAD;
    let y = (y0 + ysums * 1.0e-6) * ASEC2RAD;
    let s = (s0 + ssums) * 1.0e-6 * ASEC2RAD - x * y / 2.0;

    // Compute expression for the celestial motion of the
    // celestial intermediate pole (CIP)
    // Equations 5.6 & 5.7 of IERS technical note 36
    let e = f64::atan2(y, x);
    let d = f64::atan(f64::sqrt((x * x + y * y) / (1.0 - x * x - y * y)));
    qrotz(-e) * qroty(-d) * qrotz(e + s)
}

use super::satrec::SatRec;
use super::initl::{initl, InitlStruct};

fn getgravconst(whichconst: &str, satrec: &mut SatRec) {
    match whichconst {
        // -- wgs-72 low precision str#3 constants --
        "wgs72old" => {
            satrec.mus = 398600.79964; // in km3 / s2
            satrec.radiusearthkm = 6378.135; // km
            satrec.xke = 0.0743669161; // reciprocal of tumin
            satrec.tumin = 1.0 / satrec.xke;
            satrec.j2 = 0.001082616;
            satrec.j3 = -0.00000253881;
            satrec.j4 = -0.00000165597;
            satrec.j3oj2 = satrec.j3 / satrec.j2;
        }
        // ------------ wgs-72 constants ------------
        "wgs72" => {
            satrec.mus = 398600.8; // in km3 / s2
            satrec.radiusearthkm = 6378.135; // km
            satrec.xke = 60.0
                / f64::sqrt(
                    satrec.radiusearthkm * satrec.radiusearthkm * satrec.radiusearthkm / satrec.mus,
                );
            satrec.tumin = 1.0 / satrec.xke;
            satrec.j2 = 0.001082616;
            satrec.j3 = -0.00000253881;
            satrec.j4 = -0.00000165597;
            satrec.j3oj2 = satrec.j3 / satrec.j2;
        }
        "wgs84" => {
            // ------------ wgs-84 constants ------------
            satrec.mus = 398600.5; // in km3 / s2
            satrec.radiusearthkm = 6378.137; // km
            satrec.xke = 60.0
                / f64::sqrt(
                    satrec.radiusearthkm * satrec.radiusearthkm * satrec.radiusearthkm / satrec.mus,
                );
            satrec.tumin = 1.0 / satrec.xke;
            satrec.j2 = 0.00108262998905;
            satrec.j3 = -0.00000253215306;
            satrec.j4 = -0.00000161098761;
            satrec.j3oj2 = satrec.j3 / satrec.j2;
        }
        _ => (),
    }
} // getgravconst

pub fn sgp4init(
    whichconst: &str,
    opsmode: &str,
    satn: &str,
    epoch: f64,
    xbstar: f64,
    xndot: f64,
    xnddot: f64,
    xecco: f64,
    xargpo: f64,
    xinclo: f64,
    xmo: f64,
    xno_kozai: f64,
    xnodeo: f64,
) -> SatRec {
    /* --------------------- local variables ------------------------ */
    /* c++ comment out
    double ao, ainv, con42, cosio, sinio, cosio2, eccsq,
    omeosq, posq, rp, rteosq,
    cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq,
    cc2, cc3, coef, coef1, cosio4, day, dndt,
    em, emsq, eeta, etasq, gam, argpm, nodem,
    inclm, mm, nm, perige, pinvsq, psisq, qzms24,
    rtemsq, s1, s2, s3, s4, s5, s6,
    s7, sfour, ss1, ss2, ss3, ss4, ss5,
    ss6, ss7, sz1, sz2, sz3, sz11, sz12,
    sz13, sz21, sz22, sz23, sz31, sz32, sz33,
    tc, temp, temp1, temp2, temp3, tsi, xpidot,
    xhdot1, z1, z2, z3, z11, z12, z13,
    z21, z22, z23, z31, z32, z33,
    qzms2t, ss, x2o3, r[3], v[3],
    delmotemp, qzms2ttemp, qzms24temp;
    */

    /* js additions */
    let ainv = 0.0;
    let ao = 0.0;
    let con42 = 0.0;
    let cosio = 0.0;
    let cosio2 = 0.0;
    let eccsq = 0.0;
    let omeosq = 0.0;
    let posq = 0.0;
    let rp = 0.0;
    let rteosq = 0.0;
    let sinio = 0.0;
    /*--------*/
    

    /* ------------------------ initialization --------------------- */
    // sgp4fix divisor for divide by zero check on inclination
    // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    const temp4: f64 = 1.5e-12;
    let mut satrec = SatRec::new();

    /* ------------------------ earth constants ----------------------- */
    // sgp4fix identify constants and allow alternate values
    // this is now the only call for the constants
    /* c++ comment otu
    getgravconst(whichconst, satrec.tumin, satrec.mus, satrec.radiusearthkm, satrec.xke,
        satrec.j2, satrec.j3, satrec.j4, satrec.j3oj2);
    */
    getgravconst(whichconst, &mut satrec);

    //-------------------------------------------------------------------------
    
    satrec.error = 0;
    satrec.operationmode = opsmode.chars().nth(0).unwrap();
    // new alpha5 or 9-digit number

    /* c++ comment out
        #ifdef _MSC_VER
    strcpy_s(satrec.satnum, 6 * sizeof(char), satn);
        #else
    strcpy(satrec.satnum, satn);
        #endif
    */
    satrec.satnum = 0;

    // sgp4fix - note the following variables are also passed directly via satrec.
    // it is possible to streamline the sgp4init call by deleting the "x"
    // variables, but the user would need to set the satrec.* values first. we
    // include the additional assignments in case twoline2rv is not used.
    satrec.bstar = xbstar;
    // sgp4fix allow additional parameters in the struct
    satrec.ndot = xndot;
    satrec.nddot = xnddot;
    satrec.ecco = xecco;
    satrec.argpo = xargpo;
    satrec.inclo = xinclo;
    satrec.mo = xmo;
    // sgp4fix rename variables to clarify which mean motion is intended
    satrec.no_kozai = xno_kozai;
    satrec.nodeo = xnodeo;

  
    /* ------------------------ earth constants ----------------------- */
    // sgp4fix identify constants and allow alternate values no longer needed
    // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
    let ss = 78.0 / satrec.radiusearthkm + 1.0;
    // sgp4fix use multiply for speed instead of pow
    let qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
    let qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
    let x2o3 = 2.0 / 3.0;

    satrec.init = 'y';
    satrec.t = 0.0;

    // sgp4fix remove satn as it is not needed in initl
    let ir =
        initl(satrec.xke, satrec.j2, satrec.ecco, epoch, satrec.inclo, satrec.no_kozai, satrec.operationmode,
            satrec.method, ainv, ao, satrec.con41, con42, cosio, cosio2, eccsq, omeosq,
            posq, rp, rteosq, sinio, satrec.gsto, satrec.no_unkozai)
    satrec.method = ir.method;
    let ainv = ir.ainv;
    let ao = ir.ao;
    satrec.con41 = ir.con41;
    let con42 = ir.con42;
    let cosio = ir.cosio;
    let cosio2 = ir.cosio2;
    let eccsq = ir.eccsq;
    let omeosq = ir.omeosq;
    let posq = ir.posq;
    let rp = ir.rp;
    let rteosq = ir.rteosq;
    let sinio = ir.sinio;
    let satrec.gsto = ir.gsto;
    let satrec.no_unkozai = ir.no_unkozai;


    satrec.a = pow(satrec.no_unkozai * satrec.tumin, (-2.0 / 3.0));
    satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
    satrec.error = 0;

    // sgp4fix remove this check as it is unnecessary
    // the mrt check in sgp4 handles decaying satellite cases even if the starting
    // condition is below the surface of te earth
    //     if (rp < 1.0)
    //       {
    //         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
    //         satrec.error = 5;
    //       }

    if ((omeosq >= 0.0) || (satrec.no_unkozai >= 0.0)) {
        satrec.isimp = 0;
        if (rp < (220.0 / satrec.radiusearthkm + 1.0))
            satrec.isimp = 1;
        sfour = ss;
        qzms24 = qzms2t;
        perige = (rp - 1.0) * satrec.radiusearthkm;

        /* - for perigees below 156 km, s and qoms2t are altered - */
        if (perige < 156.0) {
            sfour = perige - 78.0;
            if (perige < 98.0)
                sfour = 20.0;
            // sgp4fix use multiply for speed instead of pow
            qzms24temp = (120.0 - sfour) / satrec.radiusearthkm;
            qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
            sfour = sfour / satrec.radiusearthkm + 1.0;
        }
        pinvsq = 1.0 / posq;

        tsi = 1.0 / (ao - sfour);
        satrec.eta = ao * satrec.ecco * tsi;
        etasq = satrec.eta * satrec.eta;
        eeta = satrec.ecco * satrec.eta;
        psisq = fabs(1.0 - etasq);
        coef = qzms24 * pow(tsi, 4.0);
        coef1 = coef / pow(psisq, 3.5);
        cc2 = coef1 * satrec.no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
            (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 *
            (8.0 + 3.0 * etasq * (8.0 + etasq)));
        satrec.cc1 = satrec.bstar * cc2;
        cc3 = 0.0;
        if (satrec.ecco > 1.0e-4)
            cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no_unkozai * sinio / satrec.ecco;
        satrec.x1mth2 = 1.0 - cosio2;
        satrec.cc4 = 2.0 * satrec.no_unkozai * coef1 * ao * omeosq *
            (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
                (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (ao * psisq) *
                (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                    (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                    (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
        satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
            (etasq + eeta) + eeta * etasq);
        cosio4 = cosio2 * cosio2;
        temp1 = 1.5 * satrec.j2 * pinvsq * satrec.no_unkozai;
        temp2 = 0.5 * temp1 * satrec.j2 * pinvsq;
        temp3 = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no_unkozai;
        satrec.mdot = satrec.no_unkozai + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 *
            temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
        satrec.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
            (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
            temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
        xhdot1 = -temp1 * cosio;
        satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
            2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
        xpidot = satrec.argpdot + satrec.nodedot;
        satrec.omgcof = satrec.bstar * cc3 * cos(satrec.argpo);
        satrec.xmcof = 0.0;
        if (satrec.ecco > 1.0e-4)
            satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
        satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
        satrec.t2cof = 1.5 * satrec.cc1;
        // sgp4fix for divide by zero with xinco = 180 deg
        if (fabs(cosio + 1.0) > 1.5e-12) {
            satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
        }
        else {
            satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
        }
        satrec.aycof = -0.5 * satrec.j3oj2 * sinio;
        // sgp4fix use multiply for speed instead of pow
        delmotemp = 1.0 + satrec.eta * cos(satrec.mo);
        satrec.delmo = delmotemp * delmotemp * delmotemp;
        satrec.sinmao = f64::sin(satrec.mo);
        satrec.x7thm1 = 7.0 * cosio2 - 1.0;

        /* --------------- deep space initialization ------------- */
        if ((2 * pi / satrec.no_unkozai) >= 225.0) {
            satrec.method = 'd';
            satrec.isimp = 1;
            tc = 0.0;
            inclm = satrec.inclo;

            /*
            let snodm = 0
            let cnodm = 0
            let sinim = 0
            let cosim = 0
            let sinomm = 0
            let cosomm = 0
            let day = 0
            let em = 0
            let emsq = 0
            let gam = 0
            let rtemsq = 0
            let s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0
            let ss1 = 0, ss2 = 0, ss3 = 0, ss4 = 0, ss5 = 0, ss6 = 0, ss7 = 0
            let sz1 = 0, sz2 = 0, sz3 = 0
            let sz11 = 0, sz12 = 0, sz13 = 0, sz21 = 0, sz22 = 0
            let sz23 = 0, sz31 = 0, sz32 = 0, sz33 = 0
            let nm = 0, z1 = 0, z2 = 0, z3 = 0
            let z11 = 0, z12 = 0, z13 = 0, z21 = 0, z22 = 0, z23 = 0
            let z31 = 0, z32 = 0, z33 = 0
            */
            let is = dscom(
                epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo,
                satrec.no_unkozai, snodm, cnodm, sinim, cosim, sinomm, cosomm,
                day, satrec.e3, satrec.ee2, em, emsq, gam,
                satrec.peo, satrec.pgho, satrec.pho, satrec.pinco,
                satrec.plo, rtemsq, satrec.se2, satrec.se3,
                satrec.sgh2, satrec.sgh3, satrec.sgh4,
                satrec.sh2, satrec.sh3, satrec.si2, satrec.si3,
                satrec.sl2, satrec.sl3, satrec.sl4, s1, s2, s3, s4, s5,
                s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3,
                sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33,
                satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
                satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
                satrec.xl3, satrec.xl4, nm, z1, z2, z3, z11,
                z12, z13, z21, z22, z23, z31, z32, z33,
                satrec.zmol, satrec.zmos
            );
            snodm = is.snodm
            cnodm = is.cnodm
            sinim = is.sinim
            cosim = is.cosim
            sinomm = is.sinomm
            cosomm = is.cosomm
            day = is.day
            satrec.e3 = is.e3
            satrec.ee2 = is.ee2
            em = is.em
            emsq = is.emsq
            gam = is.gam
            satrec.peo = is.peo
            satrec.pgho = is.pgho
            satrec.pho = is.pho
            satrec.pinco = is.pinco
            satrec.plo = is.plo
            rtemsq = is.rtemsq
            satrec.se2 = is.se2
            satrec.se3 = is.se3
            satrec.sgh2 = is.sgh2
            satrec.sgh3 = is.sgh3
            satrec.sgh4 = is.sgh4
            satrec.sh2 = is.sh2
            satrec.sh3 = is.sh3
            satrec.si2 = is.si2
            satrec.si3 = is.si3
            satrec.sl2 = is.sl2
            satrec.sl3 = is.sl3
            satrec.sl4 = is.sl4
            s1 = is.s1
            s2 = is.s2
            s3 = is.s3
            s4 = is.s4
            s5 = is.s5
            s6 = is.s6
            s7 = is.s7
            ss1 = is.ss1
            ss2 = is.ss2
            ss3 = is.ss3
            ss4 = is.ss4
            ss5 = is.ss5
            ss6 = is.ss6
            ss7 = is.ss7
            sz1 = is.sz1
            sz2 = is.sz2
            sz3 = is.sz3
            sz11 = is.sz11
            sz12 = is.sz12
            sz13 = is.sz13
            sz21 = is.sz21
            sz22 = is.sz22
            sz23 = is.sz23
            sz31 = is.sz31
            sz32 = is.sz32
            sz33 = is.sz33
            satrec.xgh2 = is.xgh2
            satrec.xgh3 = is.xgh3
            satrec.xgh4 = is.xgh4
            satrec.xh2 = is.xh2
            satrec.xh3 = is.xh3
            satrec.xi2 = is.xi2
            satrec.xi3 = is.xi3
            satrec.xl2 = is.xl2
            satrec.xl3 = is.xl3
            satrec.xl4 = is.xl4
            nm = is.nm
            z1 = is.z1
            z2 = is.z2
            z3 = is.z3
            z11 = is.z11
            z12 = is.z12
            z13 = is.z13
            z21 = is.z21
            z22 = is.z22
            z23 = is.z23
            z31 = is.z31
            z32 = is.z32
            z33 = is.z33
            satrec.zmol = is.zmol
            satrec.zmos = is.zmos

            let dpr = dpper(
                satrec.e3, satrec.ee2, satrec.peo, satrec.pgho,
                satrec.pho, satrec.pinco, satrec.plo, satrec.se2,
                satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
                satrec.sh2, satrec.sh3, satrec.si2, satrec.si3,
                satrec.sl2, satrec.sl3, satrec.sl4, satrec.t,
                satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
                satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
                satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos, inclm, satrec.init,
                satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
                satrec.operationmode
            );
            satrec.ecco = dpr.ep
            satrec.inclo = dpr.inclp
            satrec.nodeo = dpr.nodep
            satrec.argpo = dpr.argpp
            satrec.mo = dpr.mp


            argpm = 0.0;
            nodem = 0.0;
            mm = 0.0;

            let dndt = 0
            let dsr = dsinit(
                satrec.xke,
                cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
                satrec.gsto, satrec.mo, satrec.mdot, satrec.no_unkozai, satrec.nodeo,
                satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
                satrec.irez, satrec.atime,
                satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222,
                satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
                satrec.d5421, satrec.d5433, satrec.dedt, satrec.didt,
                satrec.dmdt, dndt, satrec.dnodt, satrec.domdt,
                satrec.del1, satrec.del2, satrec.del3, satrec.xfact,
                satrec.xlamo, satrec.xli, satrec.xni
            )
            satrec.irez = dsr.irez
            satrec.atime = dsr.atime
            satrec.d2201 = dsr.d2201
            satrec.d2211 = dsr.d2211
            satrec.d3210 = dsr.d3210
            satrec.d3222 = dsr.d3222
            satrec.d4410 = dsr.d4410
            satrec.d4422 = dsr.d4422
            satrec.d5220 = dsr.d5220
            satrec.d5232 = dsr.d5232
            satrec.d5421 = dsr.d5421
            satrec.d5433 = dsr.d5433
            satrec.dedt = dsr.dedt
            satrec.didt = dsr.didt
            satrec.dmdt = dsr.dmdt
            dndt = dsr.dndt
            satrec.dnodt = dsr.dnodt
            satrec.domdt = dsr.domdt
            satrec.del1 = dsr.del1
            satrec.del2 = dsr.del2
            satrec.del3 = dsr.del3
            satrec.xfact = dsr.xfact
            satrec.xlamo = dsr.xlamo
            satrec.xli = dsr.xli
            satrec.xni = dsr.xni
        }

        /* ----------- set variables if not deep space ----------- */
        if (satrec.isimp != 1) {
            cc1sq = satrec.cc1 * satrec.cc1;
            satrec.d2 = 4.0 * ao * tsi * cc1sq;
            temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
            satrec.d3 = (17.0 * ao + sfour) * temp;
            satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
                satrec.cc1;
            satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
            satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
                (12.0 * satrec.d2 + 10.0 * cc1sq));
            satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
                12.0 * satrec.cc1 * satrec.d3 +
                6.0 * satrec.d2 * satrec.d2 +
                15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
        }
    } // if omeosq = 0 ...

    /* finally propogate to zero epoch to initialize all others. */
    // sgp4fix take out check to let satellites process until they are actually below earth surface
    //       if(satrec.error == 0)
    sgp4(satrec, 0.0);

    satrec.init = 'n';

    //#include "debug6.cpp"
    //sgp4fix return boolean. satrec.error contains any error codes
    return true;
    
    satrec
} // sgp4init

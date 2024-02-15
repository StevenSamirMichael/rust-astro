# %%
import satkit as sk
import pytest
import numpy as np
import math as m


class TestJPLEphem:
    def test_jplephem_testvecs(self):
        """
        Test JPL ephemeris against test vectors provided by JPL
        """

        # File contains test calculation vectors provided by NASA
        fname = "../satkit-testvecs/jplephem/testpo.440"

        # Read in the test vectors
        with open(fname, "r") as fd:
            lines = fd.readlines()

        # Function to convert integer index to solar system body
        def int_to_ss(ix: int):
            if ix == 0:
                return sk.solarsystem.Mercury
            elif ix == 1:
                return sk.solarsystem.Venus
            elif ix == 2:
                return sk.solarsystem.EMB
            elif ix == 3:
                return sk.solarsystem.Mars
            elif ix == 4:
                return sk.solarsystem.Jupiter
            elif ix == 5:
                return sk.solarsystem.Saturn
            elif ix == 6:
                return sk.solarsystem.Uranus
            elif ix == 7:
                return sk.solarsystem.Neptune
            elif ix == 8:
                return sk.solarsystem.Pluto
            elif ix == 9:
                return sk.solarsystem.Moon
            elif ix == 10:
                return sk.solarsystem.Sun

        # Go through the test vectors
        # each test vecxtor is a line in the file
        for line in lines[14:]:
            s = line.split()
            assert len(s) >= 7
            # get the fields in the test vector
            jd = float(s[2])
            tar = int(s[3])
            src = int(s[4])
            coord = int(s[5])
            truth = float(s[6])
            time = sk.time.from_jd(jd, sk.timescale.TT)
            # Don't handle any of the exotic test vectors, just do sun, moon,
            # and planetary ephemerides
            if tar <= 10 and src <= 10 and coord <= 6:
                sksrc = int_to_ss(src - 1)
                sktar = int_to_ss(tar - 1)
                tpos, tvel = sk.jplephem.geocentric_state(sktar, time)
                spos, svel = sk.jplephem.geocentric_state(sksrc, time)

                # In test vectors, index 3 is not EMB, but Earth
                # (not obvious...)
                if tar == 3:
                    _mpos, mvel = sk.jplephem.geocentric_state(
                        sk.solarsystem.Moon, time
                    )
                    tvel = tvel - mvel / (1.0 + sk.consts.earth_moon_mass_ratio)
                    tpos = np.array([0, 0, 0])
                if src == 3:
                    spos = np.array([0, 0, 0])
                    _mpos, mvel = sk.jplephem.geocentric_state(
                        sk.solarsystem.Moon, time
                    )
                    svel = svel - mvel / (1.0 + sk.consts.earth_moon_mass_ratio)
                if src == 10:
                    embpos, embvel = sk.jplephem.geocentric_state(
                        sk.solarsystem.EMB, time
                    )
                    svel = svel + (
                        embvel - svel / (1.0 + sk.consts.earth_moon_mass_ratio)
                    )
                if tar == 10:
                    embpos, embvel = sk.jplephem.geocentric_state(
                        sk.solarsystem.EMB, time
                    )
                    tvel = tvel + (
                        embvel - tvel / (1.0 + sk.consts.earth_moon_mass_ratio)
                    )
                # Position test
                if coord <= 3:
                    calc = (tpos - spos)[coord - 1] / sk.consts.au
                    assert calc == pytest.approx(truth, rel=1e-12)
                # Velocity test
                else:
                    calc = (tvel - svel)[coord - 4] / sk.consts.au * 86400.0
                    assert calc == pytest.approx(truth, rel=1e-12)


class TestGravity:
    def test_gravity(self):
        """
        Reference gravity computations, using
        JGM3 model, with 16 terms, found at:
        http://icgem.gfz-potsdam.de/calcstat/
        Outputs from above web page listed
        in function below as reference values
        """

        reference_gravitation = 9.822206169031
        reference_gravity = 9.803696372738
        # Gravity deflection from normal along east-west and
        # north-south direction, in arcseconds
        reference_ew_deflection_asec = -1.283542043355
        reference_ns_deflection_asec = -1.311709802440

        latitude_deg = 42.4473
        longitude_deg = -71.2272
        altitude = 0

        itrf = sk.itrfcoord(
            latitude_deg=latitude_deg, longitude_deg=longitude_deg, altitude=altitude
        )
        gravitation = sk.gravity(itrf, order=16)
        # Add centrifugal force @ Earth surface
        centrifugal = (
            np.array([itrf.vector[0], itrf.vector[1], 0]) * sk.consts.omega_earth**2
        )

        assert np.linalg.norm(gravitation) == pytest.approx(
            reference_gravitation, rel=1e-9
        )
        gravity = gravitation + centrifugal
        assert np.linalg.norm(gravity) == pytest.approx(reference_gravity, rel=1e-9)

        # Rotate gravity into East-North-Up frame to check deflections
        gravity_enu = itrf.qenu2itrf.conj * gravity
        ew_deflection = (
            -m.atan2(gravity_enu[0], -gravity_enu[2]) * 180.0 / m.pi * 3600.0
        )
        ns_deflection = (
            -m.atan2(gravity_enu[1], -gravity_enu[2]) * 180.0 / m.pi * 3600.0
        )
        assert ns_deflection == pytest.approx(reference_ns_deflection_asec, rel=2e-6)
        assert ew_deflection == pytest.approx(reference_ew_deflection_asec, rel=2e-6)


class TestFrameTransform:
    def test_itrf2gcrf(self):
        """
        IAU-2000 Reduction,
        Vallado Example 3-14
        """

        pITRF = np.array([-1033.479383, 7901.2952754, 6380.3565958]) * 1e3
        vITRF = np.array([-3.225636520, -2.872451450, 5.531924446]) * 1.0e3
        tm = sk.time(2004, 4, 6, 7, 51, 28.386009)

        # The example looks up dut1 from start of day, along with x and y polar motion
        # Intermediate check by getting values and comparing against
        # values used in example
        tm2 = sk.time(2004, 4, 6, 0, 0, 0)
        dut1, xp, yp, lod, dX, dY = sk.frametransform.earth_orientation_params(tm2)
        assert dut1 == pytest.approx(-0.4399619, rel=1e-3)
        assert xp == pytest.approx(-0.140682, rel=1e-3)
        assert yp == pytest.approx(0.333309, rel=1e-3)
        jd_tt = tm.to_jd(sk.timescale.TT)
        assert jd_tt == pytest.approx(2453101.828154745)
        t_tt = (tm.to_jd(sk.timescale.TT) - 2451545.0) / 36525.0

        assert t_tt == pytest.approx(0.0426236319, rel=1e-8)
        # Check transform to terrestial intermediate frame
        # with value from example
        pTIRS = sk.frametransform.qitrf2tirs(tm) * pITRF
        assert pTIRS[0] == pytest.approx(-1033475.0312, rel=1e-7)
        assert pTIRS[1] == pytest.approx(7901305.5856, rel=1e-7)
        assert pTIRS[2] == pytest.approx(6380344.5327, rel=1e-7)

        # Check transfomr to celestial intermediate frame
        # with value from example
        pCIRS = sk.quaternion.rotz(sk.frametransform.earth_rotation_angle(tm)) * pTIRS
        assert pCIRS[0] == pytest.approx(5100018.4047, rel=1e-7)
        assert pCIRS[1] == pytest.approx(6122786.3648, rel=1e-7)
        assert pCIRS[2] == pytest.approx(6380344.6237, rel=1e-7)

        # Check transform to geocentric celestial reference frame
        # with value from example
        pGCRF = sk.frametransform.qcirs2gcrf(tm) * pCIRS
        assert pGCRF[0] == pytest.approx(5102508.959, rel=1e-7)
        assert pGCRF[1] == pytest.approx(6123011.403, rel=1e-7)
        assert pGCRF[2] == pytest.approx(6378136.925, rel=1e-7)

        # Now, test the whole transform at once
        pGCRF = sk.frametransform.qitrf2gcrf(tm) * pITRF
        assert pGCRF[0] == pytest.approx(5102508.959)
        assert pGCRF[1] == pytest.approx(6123011.403)
        assert pGCRF[2] == pytest.approx(6378136.925)

    def test_gmst(self):
        """
        Test GMST : vallado example 3-5
        """

        tm = sk.time(1992, 8, 20, 12, 14, 0)

        # Spooof UTC as UT1 value (as is done in example from Vallado)
        tdiff = tm.to_mjd(sk.timescale.UT1) - tm.to_mjd(sk.timescale.UTC)
        tm = tm - sk.duration.from_days(tdiff)
        gmst = sk.frametransform.gmst(tm)
        truth = -207.4212121875 * m.pi / 180
        assert gmst == pytest.approx(truth)


class TestITRFCoord:
    def test_geodetic(self):
        """
        Test geodetic conversions
        """
        latitude_deg = 42.46
        longitude_deg = -71.1516
        altitude = 1000
        itrf = sk.itrfcoord(
            latitude_deg=latitude_deg, longitude_deg=longitude_deg, altitude=altitude
        )
        assert itrf.latitude_deg == pytest.approx(latitude_deg)
        assert itrf.longitude_deg == pytest.approx(longitude_deg)
        assert itrf.altitude == pytest.approx(altitude)

    def test_geodetic2(self):
        """
        Vallado example 3.3
        """
        itrf = sk.itrfcoord(6524.834 * 1e3, 6862.875 * 1e3, 6448.296 * 1e3)
        assert itrf.latitude_deg == pytest.approx(34.352496)
        assert itrf.longitude_deg == pytest.approx(46.4464)


class TestMoon:
    def test_moonpos(self):
        """
        Vallado example 5-3 for
        computing position of the moon
        """
        t0 = sk.time(1994, 4, 28)
        t1 = sk.time.from_mjd(t0.to_mjd(sk.timescale.UTC), sk.timescale.TDB)
        p = sk.moon.pos_gcrf(t1)
        ref_pos = np.array([-134240.626e3, -311571.590e3, -126693.785e3])
        assert p == pytest.approx(ref_pos)


class TestSGP4:
    def test_sgp4_multiple(self):
        """
        Check propagating multiple TLEs at once
        """

        lines = [
            "0 STARLINK-3118",
            "1 49140U 21082L   24030.39663557  .00000076  00000-0  14180-4 0  9995",
            "2 49140  70.0008  34.1139 0002663 260.3521  99.7337 14.98327656131736",
            "0 STARLINK-3093",
            "1 49141U 21082M   24030.50141584 -.00000431  00000-0 -28322-4 0  9990",
            "2 49141  70.0000  73.8654 0002647 256.8611 103.2253 14.98324813131968",
            "0 STARLINK-3042",
            "1 49142U 21082N   24030.19218442  .00000448  00000-0  45331-4 0  9999",
            "2 49142  70.0005  34.6319 0002749 265.6056  94.4790 14.98327526131704",
            "0 STARLINK-3109",
            "1 49143U 21082P   24030.20076173 -.00000320  00000-0 -19071-4 0  9998",
            "2 49143  70.0002  54.6139 0002526 255.5608 104.5271 14.98327699131201",
        ]
        tles = sk.TLE.from_lines(lines)
        print(tles)
        tm = [
            sk.time(2024, 1, 15) + sk.duration.from_seconds(x * 10) for x in range(100)
        ]
        [p, v] = sk.sgp4(tles, tm)
        [p2, v2] = sk.sgp4(tles[2], tm)
        # Verify that propagating multiple TLEs matches propagation of a single TLE
        assert p2 == pytest.approx(np.squeeze(p[2, :, :]))
        assert v2 == pytest.approx(np.squeeze(v[2, :, :]))

    def test_sgp4_vallado(self):
        """
        SGP4 Test Vectors from vallado
        """
        basedir = "../satkit-testvecs/sgp4"
        tlefile = basedir + "/SGP4-VER.TLE"
        with open(tlefile, "r") as fh:
            lines = fh.readlines()

        lines = list(filter(lambda x: x[0] != "#", lines))

        tles = sk.TLE.from_lines(lines)
        for tle in tles:
            fname = f"{basedir}/{tle.satnum:05}.e"
            with open(fname, "r") as fh:
                testvecs = fh.readlines()
            for testvec in testvecs:
                stringvals = testvec.split()

                # Valid lines are all floats of length 7
                if len(stringvals) != 7:
                    continue
                try:
                    vals = [float(s) for s in stringvals]
                except ValueError:
                    continue
                time = tle.epoch + sk.duration.from_seconds(vals[0])
                try:
                    [p, v] = sk.sgp4(
                        tle,
                        time,
                        opsmode=sk.opsmode.afspc,
                        gravconst=sk.gravconst.wgs72,
                    )
                    ptest = np.array([vals[1], vals[2], vals[3]]) * 1e3
                    vtest = np.array([vals[4], vals[5], vals[6]]) * 1e3
                    print(f"p = {p} ;: ptest = {ptest}")
                    print(f"v = {v} :: vtest = {vtest}")
                    assert p == pytest.approx(ptest, rel=1e-4)
                    assert v == pytest.approx(vtest, rel=1e-2)
                except RuntimeError:
                    print("Caught runtime error; this is expected in test vectors")


# %%

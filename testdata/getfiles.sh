#!/bin/sh
curl http://celestrak.org/software/vallado/cpp.zip --output vallado.zip
unzip ./vallado.zip
mv cpp vallado

# Get JPL Ephemeris test data
curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/testpo.440 --output testpo.440

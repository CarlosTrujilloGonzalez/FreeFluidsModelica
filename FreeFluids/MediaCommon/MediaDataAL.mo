within FreeFluids.MediaCommon;

  package MediaDataAL
   "MediaDataAL.mo by Carlos Trujillo
      This file is part of the Free Fluids application
      Copyright (C) 2008-2021  Carlos Trujillo Gonzalez
        
      This program is free software; you can redistribute it and/or
      modify it under the terms of the GNU General Public License version 3
      as published by the Free Software Foundation
        
      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.
        
      You should have received a copy of the GNU General Public License
      along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA."
 
    constant  FreeFluids.MediaCommon.DataRecord Acetone(name="Acetone", CAS = "67-64-1", family = 12, MW = 5.807900e+01, molarMass=0.058079, Tc = 5.081000e+02, criticalPressure = 4.700000e+06, Vc = 2.090000e-04, Zc = 2.320000e-01, w = 3.040000e-01, mu=2.88, Tb = 3.291500e+02, IsothComp = 0.0, lnuA = 0.0111638, lnuB = -5.723473,
     Cp0Corr = 5, Cp0Coef = {4.000000e+00, 5.330000e+05, 3.320000e+02, 4.221700e+01, -2.769400e+01, 4.681000e+07, 3.210000e+02, 1.320000e+02, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 2.000000e+02, Cp0LimS = 1.500000e+03,
     VpCorr = 20, VpCoef = {7.277713e+01, -5.752936e+03, -7.680083e+00, 6.830760e-06, 2.000000e+00, 0.0}, VpLimI = 1.784500e+02, VpLimS = 5.082000e+02,
     HvCorr = 91, HvCoef = {5.136940e+05, -2.047660e+00, 6.554810e+00, -7.019220e+00, 2.833360e+00, 5.081000e+02}, HvLimI = 2.031500e+02, HvLimS = 5.081000e+02,
     lDensCorr = 46, lDensCoef = {2.778900e+02, 5.816900e+02, -1.877150e+00, 1.405200e+02, 7.669280e+01, 5.081000e+02}, lDensLimI = 1.780000e+02, lDensLimS = 5.081000e+02,
     lCpCorr = 19, lCpCoef = {5.763940e+01, 2.655400e+03, -2.221130e+03, 1.445760e+03, 2.770930e+02, 5.081000e+02}, lCpLimI = 1.780000e+02, lCpLimS = 4.980000e+02,
     lViscCorr = 30, lViscCoef = {-1.406400e+01, 1.000700e+03, 4.534900e-01, 3.945600e-07, 2.000000e+00, 0.0}, lViscLimI = 1.900000e+02, lViscLimS = 3.294400e+02,
     lThCondCorr = 51, lThCondCoef = {1.013000e-02, -9.532000e+01, -2.115100e-01, -5.261600e-03, 2.304300e-06, 0.0}, lThCondLimI = 1.784500e+02, lThCondLimS = 3.431500e+02,
     lSurfTensCorr = 61, lSurfTensCoef = {6.220000e-02, 1.124000e+00, 0.0, 0.0, 0.0, 5.082000e+02}, lSurfTensLimI = 1.784500e+02, lSurfTensLimS = 5.082000e+02,
     lBulkModRCorr = 150, lBulkModRCoef = {-1.802980e+01, 7.398040e-02, -1.133730e-04, 8.759310e-08, -2.476600e-11, 0.000000e+00}, lBulkModRLimI = 1.830000e+02, lBulkModRLimS = 4.880000e+02,
     gSatDensCorr = 101, gSatDensCoef = {2.778900e+02, -3.218990e+00, -3.461690e+00, -1.167140e+01, -5.816410e+01, 5.081000e+02}, gSatDensLimI = 2.030000e+02, gSatDensLimS = 5.070000e+02,
     gViscCorr = 110, gViscCoef = {3.101200e-08, 9.761600e-01, 2.304200e+01, 1.483400e+01, 0.0, 0.0}, gViscLimI = 1.784500e+02, gViscLimS = 1.000000e+03,
     gThCondCorr = 120, gThCondCoef = {-2.688200e+01, 9.036000e-01, -1.209500e+08, -6.087900e+08, 0.0, 0.0}, gThCondLimI = 2.731500e+02, gThCondLimS = 1.000000e+03);
    
    constant  FreeFluids.MediaCommon.DataRecord Air(name="Air",MW = 2.896000e+01, molarMass=0.02896, Tc = 1.324400e+02, criticalPressure = 3.790000e+06, Vc = 9.150000e-05, Zc = 3.147819e-01, w = 3.130000e-01, Tb = 7.867001e+01, IsothComp = 0.0, 
    Cp0Corr = 200, Cp0Coef = {9.999320e+02, 3.305700e+02, 3.098800e+03, 2.704480e+02, 1.498290e+03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 5.000000e+01, Cp0LimS = 1.500000e+03, 
    VpCorr = 20, VpCoef = {2.166200e+01, -6.923900e+02, -3.920000e-01, 4.757400e-03, 1.000000e+00, 0.0}, VpLimI = 5.915000e+01, VpLimS = 1.324500e+02,
    HvCorr = 90, HvCoef = {8.474000e+06, 3.822000e-01, 0.0, 0.0, 0.0, 1.324500e+02}, HvLimI = 5.915000e+01, HvLimS = 1.324500e+02,
    lDensCorr = 41, lDensCoef = {2.673100e+00, 2.563700e-01, 1.325100e+02, 2.678800e-01, 0.0, 0.0}, lDensLimI = 5.915000e+01, lDensLimS = 1.325000e+02,
    lCpCorr = 16, lCpCoef = {-2.144600e+05, 9.185100e+03, -1.061200e+02, 4.161600e-01, 0.0, 0.0}, lCpLimI = 7.500000e+01, lCpLimS = 1.150000e+02,
    lTfromHsatCorr = 140, lTfromHsatCoef = {1.078810e+02, 4.513260e-04, -1.068610e-09, -3.516190e-15, -1.128870e-20, 0.0}, lTfromHsatLimI = 7.500000e+01, lTfromHsatLimS = 1.150000e+02,
    lViscCorr = 30, lViscCoef = {-7.233600e+01, 8.134800e+02, 1.268700e+01, -3.306200e-04, 2.000000e+00, 0.0}, lViscLimI = 5.915000e+01, lViscLimS = 1.300000e+02,
    lThCondCorr = 51, lThCondCoef = {-2.119900e-01, -1.631100e+01, -2.305700e-01, -7.619700e-03, 2.501800e-06, 0.0}, lThCondLimI = 7.500000e+01, lThCondLimS = 1.250000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {2.919000e-02, 1.156560e+00, 6.889000e-02, 1.791800e-01, -1.456400e-01, 1.325300e+02}, lSurfTensLimI = 6.305000e+01, lSurfTensLimS = 0.0,
    gSatDensCorr = 101, gSatDensCoef = {3.165030e+02, -1.602760e+00, -4.575450e+00, -6.682720e+00, -3.960370e+01, 1.325000e+02}, gSatDensLimI = 5.315000e+01, gSatDensLimS = 1.325000e+02,
    gViscCorr = 110, gViscCoef = {1.425000e-06, 5.039000e-01, 1.083000e+02, 0.0, 0.0, 0.0}, gViscLimI = 8.000000e+01, gViscLimS = 2.000000e+03,
    gThCondCorr = 120, gThCondCoef = {3.141700e-04, 7.786000e-01, -7.116000e-01, 2.121700e+03, 0.0, 0.0}, gThCondLimI = 7.000000e+01, gThCondLimS = 2.000000e+03);

  constant FreeFluids.MediaCommon.DataRecord Ammonia(
    name = "Ammonia", description = "", CAS = "7664-41-7", family = 15, MW = 1.703100e+01, molarMass=0.017031, Tc = 4.056000e+02, criticalPressure = 1.135000e+07, Vc = 7.250000e-05, Zc = 2.440000e-01, w = 2.560000e-01, Tb = 2.396500e+02, mu = 1.500000e+00, lnuA = 1.177011e-02, lnuB = -5.143356e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 4.830000e+06, 1.727000e+03, 1.385000e+00, 8.872000e+00, -2.077800e+08, 8.280000e+02, 2.000000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 3.000000e+03,
    VpCorr = 20, VpCoef = {9.048300e+01, -4.669700e+03, -1.160700e+01, 1.719400e-02, 1.000000e+00, 0.000000e+00}, VpLimI = 1.954100e+02, VpLimS = 4.056500e+02,
    HvCorr = 90, HvCoef = {2.454200e+07, -1.317800e+00, 4.719400e+00, -5.480800e+00, 2.419600e+00, 4.054000e+02}, HvLimI = 1.954100e+02, HvLimS = 4.031500e+02,
    lDensCorr = 41, lDensCoef = {3.538300e+00, 2.544300e-01, 4.056500e+02, 2.888000e-01, 0.000000e+00, 0.000000e+00}, lDensLimI = 1.954100e+02, lDensLimS = 4.056500e+02,
    lCpCorr = 19, lCpCoef = {2.434490e+02, 4.120130e+03, -2.778720e+03, 1.006620e+04, -1.205040e+04, 4.056000e+02}, lCpLimI = 1.960000e+02, lCpLimS = 3.880000e+02,
    lViscCorr = 30, lViscCoef = {-6.743000e+00, 5.983000e+02, -7.341000e-01, -3.690000e-27, 1.000000e+01, 0.000000e+00}, lViscLimI = 1.954100e+02, lViscLimS = 3.931500e+02,
    lThCondCorr = 50, lThCondCoef = {1.169000e+00, -2.314000e-03, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 1.954100e+02, lThCondLimS = 4.000500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {1.017500e-01, 1.217030e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 4.055000e+02}, lSurfTensLimI = 1.953500e+02, lSurfTensLimS = 3.980000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {1.057600e+00, -5.410970e-02, 2.279270e-04, -3.210800e-07, 1.595860e-10, 0.000000e+00}, lBulkModRLimI = 1.960000e+02, lBulkModRLimS = 3.830000e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.349100e+02, -3.178170e+00, -3.268290e+00, -1.100020e+01, -5.214070e+01, 4.056000e+02}, gSatDensLimI = 1.960000e+02, gSatDensLimS = 3.980000e+02,
    gViscCorr = 110, gViscCoef = {4.185500e-08, 9.806000e-01, 3.080000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.954100e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {9.660800e-06, 1.379900e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.000000e+02, gThCondLimS = 9.000000e+02); 

  constant FreeFluids.MediaCommon.DataRecord Benzene(
    name = "Benzene", description = "", CAS = "71-43-2", family = 5, MW = 7.811400e+01, molarMass = 7.811400e-02, Tc = 5.620100e+02, criticalPressure = 4.894000e+06, Vc = 2.562788e-04, Zc = 2.684000e-01, w = 2.120000e-01, Tb = 3.533000e+02, lnuA = 1.048913e-02, lnuB = -5.561958e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 5.884000e+06, 8.240000e+02, 2.945500e+01, 1.024600e+01, -2.075000e+07, 9.900000e+01, 2.020000e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 2.000000e+02, Cp0LimS = 1.500000e+03,
    VpCorr = 26, VpCoef = {4.894000e+06, 5.620200e+02, 5.561906e-03, 3.700000e-02, -8.662137e-02, 5.050000e-01, -6.964183e+00, 1.014000e+00, 1.124929e+00, 1.469000e+00, -3.961859e+00, 3.711000e+00, -1.310688e+01, 1.264700e+01}, VpLimI = 2.786740e+02, VpLimS = 5.620200e+02,
    HvCorr = 90, HvCoef = {4.881000e+07, 6.106600e-01, -2.588200e-01, 3.223800e-02, 2.247500e-02, 5.620500e+02}, HvLimI = 2.731000e+02, HvLimS = 5.620200e+02,
    lDensCorr = 241, lDensCoef = {3.902000e+03, 5.620200e+02, 2.852588e+00, 4.070000e-01, -5.596548e-01, 5.650000e-01, 1.487205e+01, 4.029000e+00, -6.642959e+01, 5.699000e+00, 1.158133e+03, 9.989000e+00, -3.128774e+03, 1.229900e+01}, lDensLimI = 2.786740e+02, lDensLimS = 5.620200e+02,
    lCpCorr = 19, lCpCoef = {5.889760e+01, 2.197430e+03, -6.094240e+02, -2.549920e+03, 2.653050e+03, 5.620200e+02}, lCpLimI = 2.786800e+02, lCpLimS = 5.581500e+02,
    lViscCorr = 30, lViscCoef = {7.511700e+00, 2.946800e+02, -2.794000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lViscLimI = 2.786800e+02, lViscLimS = 5.450000e+02,
    lThCondCorr = 50, lThCondCoef = {2.344400e-01, -3.057200e-04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 2.786800e+02, lThCondLimS = 4.131000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {7.178000e-02, 1.235930e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.620100e+02}, lSurfTensLimI = 2.786500e+02, lSurfTensLimS = 5.620100e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-2.359460e+01, 9.419370e-02, -1.398300e-04, 1.006980e-07, -2.659620e-11, 0.000000e+00}, lBulkModRLimI = 2.790000e+02, lBulkModRLimS = 5.380000e+02,
    gSatDensCorr = 102, gSatDensCoef = {3.902000e+03, 5.620200e+02, 1.642774e-02, 6.700000e-02, -2.618439e+00, 3.870000e-01, -2.590438e+00, 8.650000e-01, 2.355352e+01, 2.692000e+00, -2.727966e+01, 2.792000e+00, -4.521840e+02, 1.561100e+01}, gSatDensLimI = 2.786740e+02, gSatDensLimS = 5.620200e+02,
    gViscCorr = 110, gViscCoef = {3.136600e-08, 9.675000e-01, 8.028500e+00, -3.562900e+01, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.731000e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {4.954900e-06, 1.451900e+00, 1.541400e+02, 2.620200e+04, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.500000e+02, gThCondLimS = 1.000000e+03); 

  constant FreeFluids.MediaCommon.DataRecord Butane_n(
    name = "n-Butane", description = "", CAS = "106-97-8", family = 1, MW = 5.812000e+01, molarMass=0.05812, Tc = 4.2511e+02, criticalPressure = 3.796000e+06, Vc = 2.549220e-04, Zc = 2.737700e-01, w = 2.010000e-01, Tb = 2.726500e+02, mu = 3.600000e-02, lnuA = 1.378235e-02, lnuB = -5.310408e+00,
    Cp0Corr = 7, Cp0Coef = {4.246805e+00, 0.000000e+00, 0.000000e+00, 5.549133e+00, 3.294041e+02, 1.146490e+01, 1.420174e+03, 7.599876e+00, 2.113089e+03, 9.660333e+00, 4.240857e+03, 0.000000e+00, 0.000000e+00}, Cp0LimI = 1.349000e+02, Cp0LimS = 1.500000e+03,
    VpCorr = 25, VpCoef = {3.796000e+06, -7.035310e+00, 1.500440e+00, -2.589980e+00, -1.537080e+00, 4.251250e+02}, VpLimI = 1.381500e+02, VpLimS = 4.251200e+02,
    HvCorr = 91, HvCoef = {5.506930e+05, 6.402530e-02, 1.184440e+00, -1.670520e+00, 8.036560e-01, 4.251250e+02}, HvLimI = 1.381500e+02, HvLimS = 4.251200e+02,
    lDensCorr = 46, lDensCoef = {2.279910e+02, 3.391410e+02, 5.986920e+02, -8.419580e+02, 5.336690e+02, 4.251250e+02}, lDensLimI = 1.381500e+02, lDensLimS = 4.251200e+02,
    lCpCorr = 19, lCpCoef = {7.330570e+01, 2.744200e+03, -1.983770e+03, 7.346490e+01, 1.392990e+03, 4.251250e+02}, lCpLimI = 1.381500e+02, lCpLimS = 4.181500e+02,
    lViscCorr = 30, lViscCoef = {-7.247100e+00, 5.348200e+02, -5.746900e-01, -4.662500e-27, 1.000000e+01, 0.000000e+00}, lViscLimI = 1.348600e+02, lViscLimS = 4.200000e+02,
    lThCondCorr = 50, lThCondCoef = {2.734900e-01, -7.126700e-04, 5.155500e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 1.348600e+02, lThCondLimS = 4.000000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.203000e-02, 1.219610e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 4.251300e+02}, lSurfTensLimI = 1.347500e+02, lSurfTensLimS = 4.251200e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.354610e+01, 5.939020e-02, -8.424990e-05, 5.623750e-08, -7.183300e-12, 0.000000e+00}, lBulkModRLimI = 1.381500e+02, lBulkModRLimS = 4.131500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.279910e+02, -2.756030e+00, -3.695520e+00, -9.457420e+00, -5.489780e+01, 4.251250e+02}, gSatDensLimI = 1.431500e+02, gSatDensLimS = 4.246500e+02,
    gViscCorr = 111, gViscCoef = {2.688000e-07, 2.513000e-08, -2.326000e-12, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.348600e+02, gViscLimS = 1.200000e+03,
    gThCondCorr = 120, gThCondCoef = {5.109400e-02, 4.525300e-01, 5.455500e+03, 1.979800e+06, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.250000e+02, gThCondLimS = 1.000000e+03); 
  
    constant  FreeFluids.MediaCommon.DataRecord Butanol_n(
    name = "n-Butanol", description = "", CAS = "71-36-3", family = 7, MW = 7.414000e+01, molarMass=0.07414, Tc = 5.629900e+02, criticalPressure = 4.414000e+06, Vc = 2.740000e-04, Zc = 2.590000e-01, w = 5.930000e-01, Tb = 3.911500e+02, mu = 1.800000e+00, lnuA = 1.056184e-02, lnuB = -4.885318e+00,
      Cp0Corr = 200, Cp0Coef = {6.854860e+02, 3.886530e+03, 1.595820e+03, 2.981980e+03, 6.516090e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 2.700000e+03,
      VpCorr = 20, VpCoef = {9.434240e+01, -9.140942e+03, -1.000380e+01, 1.706963e-06, 2.000000e+00, 0.000000e+00}, VpLimI = 1.845100e+02, VpLimS = 5.630000e+02,
      HvCorr = 91, HvCoef = {1.269800e+06, 2.033450e+00, -4.340620e+00, 4.646490e+00, -1.798610e+00, 5.630000e+02}, HvLimI = 1.850000e+02, HvLimS = 5.630000e+02,
      lDensCorr = 46, lDensCoef = {2.700000e+02, 7.772534e+02, -4.468420e+02, 5.788813e+02, -1.729538e+02, 5.630500e+02}, lDensLimI = 1.840000e+02, lDensLimS = 5.630500e+02,
      lCpCorr = 17, lCpCoef = {-1.201990e+02, 3.573830e+01, -2.491730e-01, 7.427300e-04, -7.179160e-07, 0.000000e+00}, lCpLimI = 1.900000e+02, lCpLimS = 4.550000e+02,
      lViscCorr = 37, lViscCoef = {4.356550e+00, 1.010950e+00, 7.422630e+02, 2.691100e+01, 2.140000e-06, 0.000000e+00}, lViscLimI = 1.900000e+02, lViscLimS = 3.920000e+02,
      lThCondCorr = 50, lThCondCoef = {2.136000e-01, -2.034000e-04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 1.838500e+02, lThCondLimS = 3.919000e+02,
      lSurfTensCorr = 61, lSurfTensCoef = {4.839000e-02, 9.106300e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.630500e+02}, lSurfTensLimI = 1.837500e+02, lSurfTensLimS = 4.130000e+02,
      lBulkModRCorr = 150, lBulkModRCoef = {-2.059890e+01, 9.864100e-02, -1.827360e-04, 1.684960e-07, -5.860350e-11, 0.000000e+00}, lBulkModRLimI = 2.331500e+02, lBulkModRLimS = 5.481500e+02,
      gSatDensCorr = 101, gSatDensCoef = {2.705840e+02, -3.600480e+00, -2.231070e+00, -2.070890e+01, -8.293410e+01, 5.630000e+02}, gSatDensLimI = 2.330000e+02, gSatDensLimS = 5.630000e+02,
      gViscCorr = 111, gViscCoef = {-1.178700e-06, 2.894000e-08, -5.708000e-12, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.840000e+02, gViscLimS = 1.000000e+03,
      gThCondCorr = 121, gThCondCoef = {2.093200e-02, -6.300000e-05, 1.789000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 3.700000e+02, gThCondLimS = 8.000000e+02); 
 
    constant  FreeFluids.MediaCommon.DataRecord CO2(name="Carbon dioxide", CAS = "124-38-9", MW = 4.401000e+01, molarMass=0.04401, Tc = 3.041300e+02, criticalPressure = 7.375000e+06, Vc = 9.400000e-05, Zc = 2.741458e-01, w = 2.390000e-01, Tb = 1.947000e+02, lnuA = 7.007000e-03, lnuB = -5.774613e+00,
    Cp0Corr = 5, Cp0Coef = {3.500000e+00, 1.447000e+06, 1.029000e+03, 1.713000e+01, -2.154200e+01, 4.795000e+08, 1.185000e+03, 5.700000e+01, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 5.000000e+01, Cp0LimS = 4.000000e+03,
    VpCorr = 25, VpCoef = {7.374990e+06, -6.956260e+00, 1.196950e+00, -3.126140e+00, 2.994480e+00, 3.041500e+02}, VpLimI = 2.170000e+02, VpLimS = 3.041500e+02,
    HvCorr = 90, HvCoef = {2.109200e+07, 3.536600e-01, -4.613400e-01, 4.355400e-01, 3.767100e-02, 3.042100e+02}, HvLimI = 2.165800e+02, HvLimS = 3.042100e+02,
    lDensCorr = 46, lDensCoef = {4.680000e+02, 8.978727e+02, 1.700410e+02, 1.690516e+02, 3.792180e+01, 3.041300e+02}, lDensLimI = 2.165800e+02, lDensLimS = 3.041300e+02,
    lCpCorr = 19, lCpCoef = {1.065700e+02, 1.092550e+03, 6.512080e+03, -3.150620e+04, 5.170820e+04, 3.041400e+02}, lCpLimI = 2.165800e+02, lCpLimS = 3.037500e+02,
    lViscCorr = 30, lViscCoef = {1.877500e+01, -4.029200e+02, -4.685400e+00, -6.917100e-26, 1.000000e+01, 0.0}, lViscLimI = 2.165800e+02, lViscLimS = 3.031500e+02,
    lThCondCorr = 51, lThCondCoef = {-2.497500e-01, -5.510600e+01, 4.173500e-01, -5.106700e-03, 2.015700e-06, 0.0}, lThCondLimI = 2.165800e+02, lThCondLimS = 3.000000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {8.167000e-02, 1.273390e+00, 0.0, 0.0, 0.0, 3.041300e+02}, lSurfTensLimI = 2.165500e+02, lSurfTensLimS = 3.042000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-7.684150e+01, 2.737650e-01, -3.747250e-04, 2.343830e-07, -5.521780e-11, 0.000000e+00}, lBulkModRLimI = 2.230000e+02, lBulkModRLimS = 2.930000e+02,
    gSatDensCorr = 101, gSatDensCoef = {4.681910e+02, -2.745970e+00, -3.755080e+00, -9.631140e+00, -7.911270e+01, 3.041400e+02}, gSatDensLimI = 2.165800e+02, gSatDensLimS = 3.041400e+02,
    gViscCorr = 110, gViscCoef = {2.246400e-06, 4.549500e-01, 2.926400e+02, 1.669100e+03, 0.0, 0.0}, gViscLimI = 1.700000e+02, gViscLimS = 1.900000e+03,
    gThCondCorr = 120, gThCondCoef = {5.804000e+00, -4.452200e-01, 7.941300e+02, 2.139600e+06, 0.0, 0.0}, gThCondLimI = 1.800000e+02, gThCondLimS = 1.500000e+03);

  constant FreeFluids.MediaCommon.DataRecord D4(
    name = "D4", description = "", CAS = "556-67-2", family = 0, MW = 2.966160e+02, molarMass = 2.966160e-01, Tc = 5.865000e+02, criticalPressure = 1.347215e+06, Vc = 9.587728e-04, Zc = 2.649000e-01, w = 5.890000e-01, Tb = 4.485500e+02, lnuA = 9.458775e-03, lnuB = -4.576029e+00,
    Cp0Corr = 7, Cp0Coef = {4.000000e+00, 0.000000e+00, 0.000000e+00, 2.927570e-01, 4.000000e+01, 3.824560e+01, 2.000000e+02, 5.897500e+01, 1.800000e+03, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 2.902500e+02, Cp0LimS = 1.200000e+03,
    VpCorr = 26, VpCoef = {1.347215e+06, 5.865000e+02, 5.345899e-01, 6.400000e-02, 4.012050e+00, 3.670000e-01, -1.028924e+00, 1.060000e-01, -5.674843e+00, 5.000000e-01, -6.250798e+00, 1.166667e+00, -8.407869e+00, 3.666667e+00}, VpLimI = 2.902500e+02, VpLimS = 5.864990e+02,
    HvCorr = 91, HvCoef = {1.920220e+05, -1.765560e+00, 6.441200e+00, -7.005400e+00, 2.685930e+00, 5.865000e+02}, HvLimI = 2.930000e+02, HvLimS = 5.830000e+02,
    lDensCorr = 241, lDensCoef = {1.043000e+03, 5.865000e+02, 1.699011e+03, 8.900000e-02, -5.407312e+03, 9.300000e-02, 3.713385e+03, 9.500000e-02, -1.408805e+02, 3.575000e-01, 1.384301e+02, 3.640000e-01, 5.512142e-01, 2.166667e+00}, lDensLimI = 2.902500e+02, lDensLimS = 5.864990e+02,
    lCpCorr = 19, lCpCoef = {1.800120e+01, 2.154080e+03, -6.701070e+02, -1.959250e+03, 2.387750e+03, 5.865000e+02}, lCpLimI = 2.930000e+02, lCpLimS = 5.830000e+02,
    lViscCorr = 30, lViscCoef = {-5.130200e+01, 3.232370e+03, 6.130820e+00, -6.678660e-06, 2.000000e+00, 0.000000e+00}, lViscLimI = 2.902500e+02, lViscLimS = 5.865000e+02,
    lThCondCorr = 50, lThCondCoef = {1.688240e-01, 1.612380e-04, -2.274540e-06, 4.780350e-09, -3.373370e-12, 0.000000e+00}, lThCondLimI = 2.902500e+02, lThCondLimS = 5.865000e+02,
    lSurfTensCorr = 64, lSurfTensCoef = {5.865000e+02, 4.246000e-02, 1.207000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lSurfTensLimI = 2.902500e+02, lSurfTensLimS = 5.865000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-5.271820e+00, 7.591170e-03, 1.362770e-05, -1.872630e-08, 7.650020e-12, 0.000000e+00}, lBulkModRLimI = 2.930000e+02, lBulkModRLimS = 5.830000e+02,
    gSatDensCorr = 102, gSatDensCoef = {1.043000e+03, 5.865000e+02, 2.484257e-01, 6.000000e-02, -2.976476e-01, 7.900000e-02, -3.189144e+00, 3.985000e-01, -3.700580e+00, 1.000000e+00, 3.110260e-01, 1.666667e+00, -1.012631e+01, 3.500000e+00}, gSatDensLimI = 2.902500e+02, gSatDensLimS = 5.864990e+02,
    gViscCorr = 111, gViscCoef = {3.280800e-07, 1.753100e-08, 2.927800e-12, -2.287900e-15, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.000000e+02, gViscLimS = 1.500000e+03,
    gThCondCorr = 121, gThCondCoef = {-4.260700e-04, 8.620600e-06, 7.163800e-08, -2.859000e-11, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.000000e+02, gThCondLimS = 1.500000e+03); 
    
  constant FreeFluids.MediaCommon.DataRecord D5(
    name = "D5", description = "", CAS = "541-02-6", family = 0, MW = 3.707700e+02, molarMass = 3.707700e-01, Tc = 6.191500e+02, criticalPressure = 1.160000e+06, Vc = 1.216000e-03, Zc = 2.740000e-01, w = 6.660000e-01, Tb = 4.861500e+02, lnuA = 9.621587e-03, lnuB = -4.516810e+00,
    Cp0Corr = 2, Cp0Coef = {0.000000e+00, 1.861484e+00, -1.403388e-03, 4.999957e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 0.000000e+00, Cp0LimS = 0.000000e+00,
    VpCorr = 26, VpCoef = {1.077700e+06, 6.183000e+02, -9.256000e+00, 1.000000e+00, 3.987000e+00, 1.500000e+00, -1.102000e+01, 2.240000e+00, -1.928600e+01, 3.480000e+00, 1.652400e+01, 2.860000e+00, -8.140000e+00, 1.160000e+01}, VpLimI = 2.000000e+02, VpLimS = 6.183000e+02,
    HvCorr = 91, HvCoef = {1.401280e+05, -3.862430e+00, 1.283490e+01, -1.412780e+01, 5.487790e+00, 6.183000e+02}, HvLimI = 2.330000e+02, HvLimS = 6.170000e+02,
    lDensCorr = 241, lDensCoef = {8.100000e+02, 6.183000e+02, 1.093800e+00, 2.500000e-01, 5.254000e+00, 7.900000e-01, -1.231000e+01, 1.330000e+00, 1.936400e+01, 1.900000e+00, -1.581000e+01, 2.520000e+00, 5.983000e+00, 3.220000e+00}, lDensLimI = 2.000000e+02, lDensLimS = 6.183000e+02,
    lCpCorr = 19, lCpCoef = {1.201210e+01, 2.191390e+03, -1.173770e+03, -1.679440e+02, 2.523120e+02, 6.183000e+02}, lCpLimI = 2.330000e+02, lCpLimS = 6.170000e+02,
    lViscCorr = 30, lViscCoef = {-1.310610e+02, 6.710780e+03, 1.836140e+01, -1.829610e-05, 2.000000e+00, 0.000000e+00}, lViscLimI = 2.730000e+02, lViscLimS = 5.730000e+02,
    lThCondCorr = 50, lThCondCoef = {1.184800e-01, 6.984290e-04, -4.435170e-06, 8.666910e-09, -5.922520e-12, 0.000000e+00}, lThCondLimI = 2.000000e+02, lThCondLimS = 6.150000e+02,
    lSurfTensCorr = 64, lSurfTensCoef = {6.191500e+02, 4.408000e-02, 1.357000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lSurfTensLimI = 2.260000e+02, lSurfTensLimS = 6.191500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {7.410940e-01, -2.079020e-02, 6.326090e-05, -5.593050e-08, 1.783120e-11, 0.000000e+00}, lBulkModRLimI = 2.330000e+02, lBulkModRLimS = 6.170000e+02,
    gSatDensCorr = 102, gSatDensCoef = {8.100000e+02, 6.183000e+02, -9.160000e-01, 2.300000e-01, -5.911000e+00, 6.800000e-01, -1.861700e+01, 2.240000e+00, -7.429000e+01, 5.100000e+00, -1.544000e+02, 1.070000e+01, -2.841000e+02, 1.890000e+01}, gSatDensLimI = 2.000000e+02, gSatDensLimS = 6.183000e+02,
    gViscCorr = 111, gViscCoef = {3.539300e-07, 1.566200e-08, 3.395700e-12, -2.323000e-15, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.000000e+02, gViscLimS = 1.500000e+03,
    gThCondCorr = 121, gThCondCoef = {-1.165000e-03, 1.088800e-05, 6.130600e-08, -2.366700e-11, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.000000e+02, gThCondLimS = 1.500000e+03); 

  constant FreeFluids.MediaCommon.DataRecord DecanoicAcid(
    name = "Decanoic Acid", description = "", CAS = "334-48-5", family = 13, MW = 1.722700e+02, molarMass = 1.722700e-01, Tc = 7.220900e+02, criticalPressure = 2.250000e+06, Vc = 6.080000e-04, Zc = 2.280000e-01, w = 8.060000e-01, Tb = 5.431500e+02, lnuA = 8.937020e-03, lnuB = -3.497030e+00,
    Cp0Corr = 200, Cp0Coef = {2.021890e+03, 1.946620e+03, 1.236670e+03, -1.107520e+05, 2.177990e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 3.040000e+02, Cp0LimS = 1.000000e+03,
    VpCorr = 20, VpCoef = {1.233600e+02, -1.468000e+04, -1.347400e+01, 1.949100e-18, 6.000000e+00, 0.000000e+00}, VpLimI = 3.045500e+02, VpLimS = 7.221000e+02,
    HvCorr = 91, HvCoef = {6.522450e+05, 3.911090e-01, 2.502280e-01, -2.917100e-01, 1.243850e-01, 7.221000e+02}, HvLimI = 3.050000e+02, HvLimS = 7.210000e+02,
    lDensCorr = 48, lDensCoef = {2.871000e-01, 2.682000e-01, 7.130000e+02, 2.686000e-01, 0.000000e+00, 0.000000e+00}, lDensLimI = 3.047500e+02, lDensLimS = 7.130000e+02,
    lCpCorr = 19, lCpCoef = {8.149990e-02, 3.845440e+03, -3.081470e+03, -7.208930e+00, 4.665580e+00, 7.221000e+02}, lCpLimI = 3.050000e+02, lCpLimS = 5.730000e+02,
    lViscCorr = 30, lViscCoef = {-1.230500e+01, 2.324100e+03, -5.549400e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lViscLimI = 3.045500e+02, lViscLimS = 5.431500e+02,
    lThCondCorr = 50, lThCondCoef = {2.060000e-01, -2.000000e-04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 3.047500e+02, lThCondLimS = 5.431500e+02,
    lSurfTensCorr = 60, lSurfTensCoef = {7.226100e-02, -2.426570e-04, 5.142380e-07, -7.110580e-10, 3.771840e-13, 0.000000e+00}, lSurfTensLimI = 3.040000e+02, lSurfTensLimS = 7.020000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-2.968950e+01, 1.481140e-01, -2.777680e-04, 2.466390e-07, -8.211930e-11, 0.000000e+00}, lBulkModRLimI = 3.050000e+02, lBulkModRLimS = 5.730000e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.833390e+02, -4.827200e+00, -1.174380e-01, -2.696370e+01, -1.000000e+02, 7.221000e+02}, gSatDensLimI = 3.050000e+02, gSatDensLimS = 5.730000e+02,
    gViscCorr = 110, gViscCoef = {7.174800e-08, 7.982000e-01, 1.093800e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 3.045500e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {1.704700e-04, 9.313000e-01, 7.576700e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 5.431500e+02, gThCondLimS = 1.000000e+03);
  
  constant FreeFluids.MediaCommon.DataRecord Dichlorodifluormethane(
    name = "Dichlorodifluormethane", description = "", CAS = "75-71-8", family = 0, MW = 1.209130e+02, molarMass=0.120913, Tc = 3.849400e+02, criticalPressure = 4.070961e+06, Vc = 2.170000e-04, Zc = 2.760054e-01, w = 1.796000e-01, Tb = 2.433600e+02, lnuA = 5.679794e-03, lnuB = -5.477153e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 5.420000e+05, 5.350000e+02, 1.003300e+01, -1.990000e-01, -6.290000e+06, 1.800000e+02, 3.800000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 9.000000e+01, Cp0LimS = 1.500000e+03,
    VpCorr = 25, VpCoef = {4.132030e+06, -7.016570e+00, 1.732240e+00, -2.979090e+00, -3.772320e-01, 3.849500e+02}, VpLimI = 1.550000e+02, VpLimS = 3.849500e+02,
    HvCorr = 91, HvCoef = {1.711340e+05, -1.989440e+00, 6.716260e+00, -7.562490e+00, 3.148360e+00, 3.849500e+02}, HvLimI = 1.231500e+02, HvLimS = 3.849500e+02,
    lDensCorr = 46, lDensCoef = {5.650000e+02, 9.148406e+02, 9.949276e+02, -1.324313e+03, 9.610344e+02, 3.851200e+02}, lDensLimI = 0.000000e+00, lDensLimS = 3.851200e+02,
    lCpCorr = 19, lCpCoef = {3.225770e+01, 9.246770e+02, -3.141500e+02, -4.214460e+02, 8.135080e+02, 3.849500e+02}, lCpLimI = 1.431500e+02, lCpLimS = 3.741500e+02,
    lViscCorr = 30, lViscCoef = {-9.000410e+00, 4.672610e+02, -1.522500e-01, 3.889360e-07, 2.000000e+00, 0.000000e+00}, lViscLimI = 1.930000e+02, lViscLimS = 3.290000e+02,
    lThCondCorr = 50, lThCondCoef = {2.016000e-01, -4.522000e-04, 9.000000e-08, -2.270000e-10, 2.120000e-13, 0.000000e+00}, lThCondLimI = 1.500000e+02, lThCondLimS = 3.631500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.898000e-02, 1.288620e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.851200e+02}, lSurfTensLimI = 1.151500e+02, lSurfTensLimS = 3.849500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.605010e+01, 3.124450e-02, -2.153910e-05, 7.423470e-09, -8.516500e-13, 0.000000e+00}, lBulkModRLimI = 1.531500e+02, lBulkModRLimS = 3.741500e+02,
    gSatDensCorr = 101, gSatDensCoef = {5.572030e+02, -2.746820e+00, -3.360380e+00, -9.954450e+00, -5.079580e+01, 3.849500e+02}, gSatDensLimI = 1.531500e+02, gSatDensLimS = 3.741500e+02,
    gViscCorr = 111, gViscCoef = {1.392000e-06, 3.765200e-08, -1.344000e-12, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.500000e+02, gViscLimS = 5.750000e+02,
    gThCondCorr = 121, gThCondCoef = {-1.873000e-03, 2.900000e-05, 3.409000e-08, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 0.000000e+00, gThCondLimS = 0.000000e+00); 
    
    constant  FreeFluids.MediaCommon.DataRecord EG(
     name="Ethylene glycol", CAS = "107-21-1", family = 8, MW = 6.207000e+01, molarMass=0.06207, Tc = 7.199900e+02, criticalPressure = 8.200000e+06, Vc = 1.706000e-04, Zc = 2.340000e-01, w = 5.070000e-01, Tb = 4.707500e+02, mu=2.2, lnuA = 8.525900e-03, lnuB = -6.581590e+00, 
     Cp0Corr = 5, Cp0Coef = {4.000000e+00, 1.700000e+05, 1.380000e+02, 1.926800e+01, 1.388600e+01, -1.217000e+07, 1.440000e+02, 9.400000e+01, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 2.980000e+02, Cp0LimS = 1.000000e+03,
     VpCorr = 20, VpCoef = {8.409000e+01, -1.041100e+04, -8.197600e+00, 1.653600e-18, 6.000000e+00, 0.0}, VpLimI = 2.601500e+02, VpLimS = 7.200000e+02, 
     HvCorr = 90, HvCoef = {8.351800e+07, 4.262500e-01, 0.0, 0.0, 0.0, 7.200000e+02}, HvLimI = 2.601500e+02, HvLimS = 7.200000e+02, 
     lDensCorr = 41, lDensCoef = {1.315000e+00, 2.512500e-01, 7.200000e+02, 2.186800e-01, 0.0, 0.0}, lDensLimI = 2.601500e+02, lDensLimS = 7.200000e+02, 
     lCpCorr = 17, lCpCoef = {1.219460e+03, 2.362250e+00, 7.808290e-03, -8.048960e-06, 0.0, 0.0}, lCpLimI = 2.600000e+02, lCpLimS = 6.330000e+02, 
     lTfromHsatCorr = 140, lTfromHsatCoef = {1.241890e+01, 7.278410e-04, -4.521010e-10, 2.319050e-16, -5.194780e-23, 0.0}, lTfromHsatLimI = 2.600000e+02, lTfromHsatLimS = 4.930000e+02, 
     lViscCorr = 30, lViscCoef = {-2.051500e+01, 2.468500e+03, 1.243500e+00, 2.499800e+12, -5.000000e+00, 0.0}, lViscLimI = 2.601500e+02, lViscLimS = 5.760000e+02, 
     lThCondCorr = 50, lThCondCoef = {1.125000e-01, 6.626000e-04, -8.800000e-08, -2.300000e-09, 1.597000e-12, 0.0}, lThCondLimI = 2.600000e+02, lThCondLimS = 4.700000e+02, 
     lSurfTensCorr = 61, lSurfTensCoef = {7.130000e-02, 7.416200e-01, 0.0, 0.0, 0.0, 7.191500e+02}, lSurfTensLimI = 2.600500e+02, lSurfTensLimS = 4.700000e+02, 
     lBulkModRCorr = 150, lBulkModRCoef = {-2.964600e+01, 1.147840e-01, -1.718190e-04, 1.233410e-07, -3.337110e-11, 0.000000e+00}, lBulkModRLimI = 2.730000e+02, lBulkModRLimS = 6.730000e+02, 
     gSatDensCorr = 101, gSatDensCoef = {3.638340e+02, -4.215380e+00, -1.016890e+00, -1.923510e+01, -7.625050e+01, 7.200000e+02}, gSatDensLimI = 2.600000e+02, gSatDensLimS = 7.030000e+02, 
     gViscCorr = 110, gViscCoef = {8.670600e-08, 8.392300e-01, 7.551200e+01, 0.0, 0.0, 0.0}, gViscLimI = 2.601500e+02, gViscLimS = 1.000000e+03, 
     gThCondCorr = 120, gThCondCoef = {-8.962900e+06, -3.125700e-01, 2.531300e+09, -1.295500e+13, 0.0, 0.0}, gThCondLimI = 2.791300e+02, gThCondLimS = 1.000000e+03);

  constant FreeFluids.MediaCommon.DataRecord Ethane(
    name = "Ethane", description = "", CAS = "74-84-0", family = 1, MW = 3.007000e+01, molarMass = 3.007000e-02, Tc = 3.052900e+02, criticalPressure = 4.872200e+06, Vc = 1.460000e-04, Zc = 2.800000e-01, w = 9.800000e-02, Tb = 1.845500e+02, mu = 0.000000e+00, lnuA = 1.605574e-02, lnuB = -5.754008e+00,
    Cp0Corr = 7, Cp0Coef = {4.003039e+00, 0.000000e+00, 0.000000e+00, 1.117433e+00, 4.302308e+02, 3.467773e+00, 1.224316e+03, 6.941945e+00, 2.014121e+03, 5.970851e+00, 4.268344e+03, 0.000000e+00, 0.000000e+00}, Cp0LimI = 9.000000e+01, Cp0LimS = 1.500000e+03,
    VpCorr = 25, VpCoef = {4.883870e+06, -6.388380e+00, 1.201830e+00, -1.637770e+00, -1.332580e+00, 3.054000e+02}, VpLimI = 9.100000e+01, VpLimS = 3.053000e+02,
    HvCorr = 91, HvCoef = {6.162920e+05, -1.880250e-01, 1.548050e+00, -1.913280e+00, 9.081470e-01, 3.053220e+02}, HvLimI = 9.115000e+01, HvLimS = 3.053000e+02,
    lDensCorr = 46, lDensCoef = {2.059590e+02, 3.031970e+02, 4.465950e+02, -5.838940e+02, 3.751560e+02, 3.053220e+02}, lDensLimI = 9.100000e+01, lDensLimS = 3.053000e+02,
    lCpCorr = 19, lCpCoef = {1.313540e+02, 2.506770e+03, -1.519650e+03, 1.087850e+03, 4.547660e+02, 3.053220e+02}, lCpLimI = 9.115000e+01, lCpLimS = 2.931500e+02,
    lViscCorr = 30, lViscCoef = {-2.305190e+01, 5.428650e+02, 2.351430e+00, -2.513680e-05, 2.000000e+00, 0.000000e+00}, lViscLimI = 9.200000e+01, lViscLimS = 2.920000e+02,
    lThCondCorr = 51, lThCondCoef = {-7.387600e-02, -9.678700e+00, -6.740500e-01, -3.407000e-03, -2.202300e-06, 0.000000e+00}, lThCondLimI = 9.035000e+01, lThCondLimS = 3.000000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {4.863000e-02, 1.198280e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.053200e+02}, lSurfTensLimI = 9.025000e+01, lSurfTensLimS = 0.000000e+00,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.818910e+01, 9.752430e-02, -1.832640e-04, 1.628330e-07, -4.386380e-11, 0.000000e+00}, lBulkModRLimI = 9.115000e+01, lBulkModRLimS = 2.931500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.059590e+02, -2.679700e+00, -3.107080e+00, -8.996600e+00, -4.330590e+01, 3.053220e+02}, gSatDensLimI = 1.181500e+02, gSatDensLimS = 3.053000e+02,
    gViscCorr = 110, gViscCoef = {5.245200e-07, 5.890600e-01, 1.888000e+02, -2.953800e+03, 0.000000e+00, 0.000000e+00}, gViscLimI = 9.035000e+01, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {7.386900e-05, 1.168900e+00, 5.007300e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 1.845500e+02, gThCondLimS = 1.000000e+03); 

  constant FreeFluids.MediaCommon.DataRecord Ethanol(
    name = "Ethanol", description = "", CAS = "64-17-5", family = 7, MW = 4.606844e+01, molarMass = 4.606844e-02, Tc = 5.139100e+02, criticalPressure = 6.268000e+06, Vc = 1.686341e-04, Zc = 2.470000e-01, w = 6.440000e-01, Tb = 3.513900e+02, mu = 1.690900e+00, lnuA = 1.021159e-02, lnuB = -4.904172e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 1.240000e+05, 2.450000e+02, 5.053900e+01, -4.946900e+01, 2.204400e+08, 5.600000e+02, 7.800000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 3.000000e+03,
    VpCorr = 24, VpCoef = {6.132000e+06, -8.685870e+00, 1.178310e+00, -4.876200e+00, 1.588000e+00, 5.139200e+02}, VpLimI = 1.590500e+02, VpLimS = 5.139200e+02,
    HvCorr = 90, HvCoef = {6.389900e+07, 1.278200e+00, -2.673000e+00, 2.797300e+00, -1.020900e+00, 5.139200e+02}, HvLimI = 1.590500e+02, HvLimS = 5.156500e+02,
    lDensCorr = 41, lDensCoef = {1.628800e+00, 2.746900e-01, 5.140000e+02, 2.317800e-01, 0.000000e+00, 0.000000e+00}, lDensLimI = 1.590000e+02, lDensLimS = 5.140000e+02,
    lCpCorr = 19, lCpCoef = {1.688760e+02, 2.850480e+03, 2.161060e+03, -1.618130e+04, 1.543090e+04, 5.147100e+02}, lCpLimI = 1.830000e+02, lCpLimS = 5.090000e+02,
    lViscCorr = 31, lViscCoef = {-6.210000e+00, 1.614000e+03, 6.180000e-03, -1.132000e-05, 0.000000e+00, 0.000000e+00}, lViscLimI = 1.680000e+02, lViscLimS = 5.160000e+02,
    lThCondCorr = 51, lThCondCoef = {1.024700e-01, -1.203900e+02, -4.848700e-01, -7.170600e-03, 3.461000e-06, 0.000000e+00}, lThCondLimI = 1.590500e+02, lThCondLimS = 3.531500e+02,
    lSurfTensCorr = 64, lSurfTensCoef = {5.139000e+02, 5.000000e-02, 9.520000e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lSurfTensLimI = 1.591000e+02, lSurfTensLimS = 5.139000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-2.099050e+01, 9.952410e-02, -1.827080e-04, 1.645080e-07, -5.525560e-11, 0.000000e+00}, lBulkModRLimI = 2.230000e+02, lBulkModRLimS = 4.980000e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.758680e+02, -4.060050e+00, -1.793080e+00, -2.039900e+01, -6.789950e+01, 5.162000e+02}, gSatDensLimI = 2.230000e+02, gSatDensLimS = 5.030000e+02,
    gViscCorr = 110, gViscCoef = {1.061300e-07, 8.066000e-01, 5.270000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.000000e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {-1.001000e-02, 6.492500e-01, -7.360500e+03, -2.552500e+05, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.731500e+02, gThCondLimS = 1.000000e+03);   

  constant FreeFluids.MediaCommon.DataRecord Heptane_n(
    name = "n-Heptane", description = "", CAS = "142-82-5", family = 1, MW = 1.002020e+02, molarMass = 1.002020e-01, Tc = 5.401200e+02, criticalPressure = 2.736000e+06, Vc = 4.320000e-04, Zc = 2.630000e-01, w = 3.490000e-01, Tb = 3.715300e+02, mu = 7.000000e-02, lnuA = 1.283770e-02, lnuB = -5.027860e+00,
    Cp0Corr = 8, Cp0Coef = {3.000000e+00, 1.372660e+01, 1.697890e+02, 3.047070e+01, 8.361950e+02, 4.355610e+01, 1.760460e+03, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 2.000000e+02, Cp0LimS = 1.500000e+03,
    VpCorr = 25, VpCoef = {2.736000e+06, -7.730760e+00, 1.529310e+00, -3.892540e+00, -2.066240e+00, 5.401300e+02}, VpLimI = 1.826500e+02, VpLimS = 5.401300e+02,
    HvCorr = 91, HvCoef = {3.748340e+05, -1.899840e+00, 6.898500e+00, -8.032640e+00, 3.383370e+00, 5.401300e+02}, HvLimI = 1.826500e+02, HvLimS = 5.401300e+02,
    lDensCorr = 46, lDensCoef = {2.319490e+02, 2.347830e+02, 1.319930e+03, -2.004350e+03, 1.153020e+03, 5.401300e+02}, lDensLimI = 1.826500e+02, lDensLimS = 5.401300e+02,
    lCpCorr = 19, lCpCoef = {5.606640e+01, 2.878390e+03, -7.214010e+02, -4.270480e+03, 4.647340e+03, 5.401300e+02}, lCpLimI = 1.826500e+02, lCpLimS = 5.301500e+02,
    lViscCorr = 30, lViscCoef = {-9.462200e+00, 8.770700e+02, -2.344500e-01, 1.402200e+22, -1.000000e+01, 0.000000e+00}, lViscLimI = 1.801500e+02, lViscLimS = 4.321600e+02,
    lThCondCorr = 51, lThCondCoef = {8.365700e-02, 4.911100e+01, -3.453600e+00, 7.798900e-03, -2.511200e-05, 0.000000e+00}, lThCondLimI = 1.825700e+02, lThCondLimS = 3.931500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.422000e-02, 1.252870e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.412300e+02}, lSurfTensLimI = 1.825500e+02, lSurfTensLimS = 5.402600e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-5.171930e+00, 4.401100e-03, 4.462820e-05, -7.282700e-08, 3.819090e-11, 0.000000e+00}, lBulkModRLimI = 1.826500e+02, lBulkModRLimS = 5.301500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.319490e+02, -2.974040e+00, -3.896820e+00, -1.178640e+01, -6.663970e+01, 5.401300e+02}, gSatDensLimI = 2.026500e+02, gSatDensLimS = 5.401300e+02,
    gViscCorr = 110, gViscCoef = {6.672000e-08, 8.283700e-01, 8.575200e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.825700e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {-7.002800e-02, 3.806800e-01, -7.049900e+03, -2.400500e+06, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.500000e+02, gThCondLimS = 1.000000e+03);

  constant FreeFluids.MediaCommon.DataRecord Hexane_n(
    name = "n-Hexane", description = "", CAS = "110-54-3", family = 1, MW = 8.617536e+01, molarMass = 8.617536e-02, Tc = 5.078100e+02, criticalPressure = 3.034000e+06, Vc = 3.695660e-04, Zc = 2.655600e+04, w = 2.990000e-01, Tb = 3.418800e+02, mu = 0.000000e+00, lnuA = 1.279508e-02, lnuB = -5.000528e+00,
    Cp0Corr = 8, Cp0Coef = {3.000000e+00, 1.169770e+01, 1.823260e+02, 2.681420e+01, 8.592070e+02, 3.861640e+01, 1.826590e+03, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 0.000000e+00, Cp0LimS = 0.000000e+00,
    VpCorr = 24, VpCoef = {3.034000e+06, -7.602190e+00, 2.051100e+00, -2.916960e+00, -2.437570e+00, 5.078200e+02}, VpLimI = 1.781500e+02, VpLimS = 5.078200e+02,
    HvCorr = 91, HvCoef = {7.554760e+05, 2.489390e+00, -5.063050e+00, 4.814630e+00, -1.756290e+00, 5.078200e+02}, HvLimI = 1.801500e+02, HvLimS = 5.078000e+02,
    lDensCorr = 46, lDensCoef = {2.341760e+02, 4.129590e+02, 4.794270e+02, -7.868220e+02, 5.742370e+02, 5.076000e+02}, lDensLimI = 1.830000e+02, lDensLimS = 5.030000e+02,
    lCpCorr = 19, lCpCoef = {6.040190e+01, 2.793810e+03, -9.045530e+02, -3.291250e+03, 3.529790e+03, 5.078200e+02}, lCpLimI = 1.801500e+02, lCpLimS = 5.011500e+02,
    lViscCorr = 37, lViscCoef = {6.950800e-01, 5.366500e-01, 7.042950e+02, 1.674700e+01, 5.657000e-05, 0.000000e+00}, lViscLimI = 1.780000e+02, lViscLimS = 4.060000e+02,
    lThCondCorr = 50, lThCondCoef = {2.331000e-01, -4.813000e-04, 7.350000e-07, -1.833000e-09, 1.680000e-12, 0.000000e+00}, lThCondLimI = 1.780000e+02, lThCondLimS = 4.000000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.500000e-02, 1.267690e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.077900e+02}, lSurfTensLimI = 1.777500e+02, lSurfTensLimS = 5.074300e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.145850e+01, 4.667580e-02, -5.865510e-05, 3.661190e-08, -4.452470e-12, 0.000000e+00}, lBulkModRLimI = 1.781500e+02, lBulkModRLimS = 5.001500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.331800e+02, -2.923360e+00, -3.859460e+00, -1.077030e+01, -6.294740e+01, 5.078200e+02}, gSatDensLimI = 1.781500e+02, gSatDensLimS = 5.071500e+02,
    gViscCorr = 111, gViscCoef = {-7.450000e-07, 2.552200e-08, -4.815000e-12, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.780000e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 121, gThCondCoef = {-3.277000e-03, 2.700000e-05, 9.650000e-08, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.900000e+02, gThCondLimS = 1.000000e+03); 
  
    constant  FreeFluids.MediaCommon.DataRecord Isobutane(
      name = "Isobutane", MW = 5.812300e+01, molarMass=0.058123, Tc = 4.077900e+02, criticalPressure = 3.640000e+06, Vc = 2.590000e-04, Zc = 2.780000e-01, w = 1.840000e-01, Tb=261.43, lnuA = 1.434751e-02, lnuB = -5.505380e+00,
      Cp0Corr = 8, Cp0Coef = {3.067140e+00, 8.975750e+00, 4.382700e+02, 5.251560e+00, 1.980180e+02, 2.514230e+01, 1.905020e+03, 1.613880e+01, 8.937650e+02, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 1.000000e+02, Cp0LimS = 1.500000e+03,
      VpCorr = 20, VpCoef = {1.084300e+02, -5.039900e+03, -1.501200e+01, 2.272500e-02, 1.000000e+00, 0.0}, VpLimI = 1.135400e+02, VpLimS = 4.078000e+02,
      HvCorr = 91, HvCoef = {5.293190e+05, 2.749540e-02, 1.509910e+00, -2.234420e+00, 1.089210e+00, 4.078000e+02}, HvLimI = 1.136500e+02, HvLimS = 4.071500e+02,
      lDensCorr = 46, lDensCoef = {2.260000e+02, 3.836237e+02, 3.636942e+02, -4.839167e+02, 3.536428e+02, 4.078100e+02}, lDensLimI = 1.135400e+02, lDensLimS = 4.078100e+02,
      lCpCorr = 19, lCpCoef = {7.678130e+01, 2.635270e+03, -2.142400e+03, 1.297250e+03, -4.830580e+02, 4.078000e+02}, lCpLimI = 1.181500e+02, lCpLimS = 4.071500e+02,
      lViscCorr = 30, lViscCoef = {-3.169540e+01, 1.151380e+03, 3.595160e+00, -1.529250e-05, 2.000000e+00, 0.0}, lViscLimI = 1.136500e+02, lViscLimS = 4.001500e+02,
      lThCondCorr = 50, lThCondCoef = {2.045500e-01, -3.658900e-04, 0.0, 0.0, 0.0, 0.0}, lThCondLimI = 1.135400e+02, lThCondLimS = 4.000000e+02,
      lSurfTensCorr = 61, lSurfTensCoef = {5.217000e-02, 1.272300e+00, 0.0, 0.0, 0.0, 4.081300e+02}, lSurfTensLimI = 1.135400e+02, lSurfTensLimS = 4.081400e+02,
       lBulkModRCorr = 150, lBulkModRCoef = {-7.795790e+00, 1.698280e-02, 3.362320e-05, -8.837230e-08, 5.914560e-11, 0.000000e+00}, lBulkModRLimI = 1.180000e+02, lBulkModRLimS = 3.980000e+02,
      gSatDensCorr = 101, gSatDensCoef = {2.244130e+02, -2.202930e+00, -5.627970e+00, -5.965620e+00, -6.205980e+01, 4.078000e+02}, gSatDensLimI = 1.136500e+02, gSatDensLimS = 4.071500e+02,
      gViscCorr = 110, gViscCoef = {1.087100e-07, 7.813500e-01, 7.063900e+01, 0.0, 0.0, 0.0}, gViscLimI = 1.500000e+02, gViscLimS = 1.000000e+03,
      gThCondCorr = 120, gThCondCoef = {8.977200e-02, 1.850100e-01, 6.392300e+02, 1.114700e+06, 0.0, 0.0}, gThCondLimI = 2.614300e+02, gThCondLimS = 1.000000e+03); 
    
annotation(
      Documentation(info = "<html>
  <body>
  <p>Contains data of the individual substances, mainly parameters for temperature dependent correlations for physical properties. The records can be used both with the ideal gas medium packages, and with the temperature dependent liquid packages.</p>
  <p><strong>There are data records for the following substances:</strong></p>
  <p>Acetone</p>
  <p>Air</p>
  <p>Ammonia</p>
  <p>Benzene</p>
  <p>Butane_n</p>
  <p>Butanol_n</p>
  <p>CO2 (Carbon dioxide)</p>
  <p>D4 (Octamethylcyclotetrasiloxane)</p>
  <p>D5 (Decamethylcyclpentasiloxane)</p>
  <p>DecanoicAcid</p>
  <p>Dichlorodifluormethane (CCl2F2)</p>
  <p>EG (Ethylene glycol)</p>
  <p>Ethane</p> 
  <p>Ethanol</p>
  <p>Heptane_n</p> 
  <p>Hexane_n</p>  
  <p>Isobutane</p>
  </body>
  </html>"));
  end MediaDataAL;

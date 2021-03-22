within FreeFluids.MediaCommon;

  package MediaDataMZ
   "MediaDataMZ.mo by Carlos Trujillo
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
    
    constant  FreeFluids.MediaCommon.DataRecord MarlothermSH(name="Marlotherm SH", MW = 2.720000e+02, molarMass=0.272, Tc = 6.329900e+02, criticalPressure = 1.000000e+00,
     VpCorr = 20, VpCoef = {6.896190e+01, -1.493480e+04, -4.751630e+00, -3.950780e-05, 1.759670e+00, 0.0}, VpLimI = 2.730000e+02, VpLimS = 6.330000e+02,
     lDensCorr = 42, lDensCoef = {4.542129e+01, 1.901230e-01, 1.205922e+03, 4.249540e-01, 0.0, 0.0}, lDensLimI = 2.730000e+02, lDensLimS = 6.330000e+02,
     lCpCorr = 17, lCpCoef = {4.787269e+02, 3.651780e+00, 3.299056e-06, 1.077924e-07, 0.0, 0.0}, lCpLimI = 2.730000e+02, lCpLimS = 6.330000e+02,
     lTfromHsatCorr = 140, lTfromHsatCoef = {5.146940e+01, 1.023560e-03, -8.673710e-10, 5.749250e-16, -1.635000e-22, 0.0}, lTfromHsatLimI = 2.730000e+02, lTfromHsatLimS = 6.330000e+02,
     lViscCorr = 30, lViscCoef = {-6.474040e+00, 1.828941e+03, -7.254549e-01, 1.663671e+25, -1.017037e+01, 0.0}, lViscLimI = 2.730000e+02, lViscLimS = 6.330000e+02,
     lThCondCorr = 50, lThCondCoef = {1.662208e-01, -1.112293e-04, -4.492504e-08, 3.211225e-11, 0.0, 0.0}, lThCondLimI = 2.730000e+02, lThCondLimS = 6.330000e+02);

  constant FreeFluids.MediaCommon.DataRecord Methane(
    name = "Methane", description = "", CAS = "74-82-8", family = 1, MW = 1.604280e+01, molarMass = 1.604280e-02, Tc = 1.905590e+02, criticalPressure = 4.599200e+06, Vc = 9.860000e-05, Zc = 2.860000e-01, w = 1.142000e-02, Tb = 1.116700e+02, mu = 0.000000e+00, lnuA = 2.026200e-02, lnuB = -5.677460e+00,
    Cp0Corr = 7, Cp0Coef = {4.001600e+00, 0.000000e+00, 0.000000e+00, 8.449000e-03, 6.480000e+02, 4.694200e+00, 1.957000e+03, 3.486500e+00, 3.895000e+03, 1.657200e+00, 5.705000e+03, 1.411500e+00, 1.508000e+04}, Cp0LimI = 5.000000e+01, Cp0LimS = 1.500000e+03,
    VpCorr = 24, VpCoef = {4.599200e+06, -6.017680e+00, 1.212630e+00, -8.586310e-01, -1.259300e+00, 1.905640e+02}, VpLimI = 9.115000e+01, VpLimS = 1.905600e+02,
    HvCorr = 91, HvCoef = {6.948980e+05, 7.101700e-01, -1.107900e+00, 9.402850e-01, -1.720100e-01, 1.905640e+02}, HvLimI = 9.115000e+01, HvLimS = 1.905600e+02,
    lDensCorr = 46, lDensCoef = {1.627060e+02, 2.922040e+02, 1.435220e+01, 1.104240e+02, -2.770250e+01, 1.905640e+02}, lDensLimI = 9.115000e+01, lDensLimS = 1.905600e+02,
    lCpCorr = 19, lCpCoef = {2.339550e+02, 3.148000e+03, -1.541750e+03, 3.363190e+03, -2.382490e+03, 1.905640e+02}, lCpLimI = 9.115000e+01, lCpLimS = 1.831500e+02,
    lViscCorr = 30, lViscCoef = {-6.157200e+00, 1.781500e+02, -9.523900e-01, -9.060600e-24, 1.000000e+01, 0.000000e+00}, lViscLimI = 9.069000e+01, lViscLimS = 1.880000e+02,
    lThCondCorr = 50, lThCondCoef = {4.176800e-01, -2.452800e-03, 3.558800e-06, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 9.069000e+01, lThCondLimS = 1.800000e+02,
    lSurfTensCorr = 62, lSurfTensCoef = {-4.070300e-02, -6.952300e+00, -2.311400e+00, -5.556500e-03, 5.759500e-06, 0.000000e+00}, lSurfTensLimI = 9.000000e+01, lSurfTensLimS = 1.931500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-2.033980e+01, 1.519000e-01, -4.212730e-04, 5.686820e-07, -2.703610e-10, 0.000000e+00}, lBulkModRLimI = 9.115000e+01, lBulkModRLimS = 1.831500e+02,
    gSatDensCorr = 101, gSatDensCoef = {1.627060e+02, -2.571110e+00, -2.735020e+00, -8.379570e+00, -3.349500e+01, 1.905640e+02}, gSatDensLimI = 9.115000e+01, gSatDensLimS = 1.905600e+02,
    gViscCorr = 110, gViscCoef = {5.254600e-07, 5.900600e-01, 1.056700e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 9.069000e+01, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {8.398300e-06, 1.426800e+00, -4.965400e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 1.116300e+02, gThCondLimS = 1.400000e+03); 

  constant FreeFluids.MediaCommon.DataRecord Methanol(
    name = "Methanol", description = "", CAS = "67-56-1", family = 7, MW = 3.205000e+01, molarMass = 3.205000e-02, Tc = 5.124990e+02, criticalPressure = 8.215850e+06, Vc = 0.000000e+00, Zc = 2.263000e-01, w = 5.625000e-01, Tb = 3.376900e+02, mu = 1.700000e+00, lnuA = 1.002214e-02, lnuB = -5.574599e+00,
    Cp0Corr = 200, Cp0Coef = {1.157930e+03, 2.662120e+03, 1.476110e+03, 1.021430e+03, 6.774020e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 2.000000e+02, Cp0LimS = 1.500000e+03,
    VpCorr = 26, VpCoef = {8.215850e+06, 5.125000e+02, -5.952507e-02, 1.110000e-01, 3.179978e+00, 6.666667e-01, -9.939771e+00, 8.333333e-01, -4.161042e-01, 1.500000e+00, -2.616089e+00, 2.166667e+00, -9.537880e-01, 5.500000e+00}, VpLimI = 1.756100e+02, VpLimS = 5.125000e+02,
    HvCorr = 90, HvCoef = {5.239000e+07, 3.682000e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.126400e+02}, HvLimI = 1.754700e+02, HvLimS = 5.126400e+02,
    lDensCorr = 241, lDensCoef = {8.520024e+03, 5.125000e+02, -5.161133e+03, 9.900000e-02, 5.793274e+03, 1.000000e-01, -8.508568e+02, 1.150000e-01, 2.212096e+02, 1.370000e-01, 8.350964e-01, 4.000000e+00, -3.224150e+02, 2.500000e+01}, lDensLimI = 1.756100e+02, lDensLimS = 5.125000e+02,
    lCpCorr = 15, lCpCoef = {1.123722e+02, -4.331599e-01, 1.112193e-03, -2.264336e-08, 0.000000e+00, 0.000000e+00}, lCpLimI = 1.754700e+02, lCpLimS = 4.930000e+02,
    lViscCorr = 31, lViscCoef = {-3.935000e+01, 4.826000e+03, 1.091000e-01, -1.127000e-04, 0.000000e+00, 0.000000e+00}, lViscLimI = 2.330000e+02, lViscLimS = 5.120000e+02,
    lThCondCorr = 51, lThCondCoef = {-5.681700e-02, 1.315600e+01, -1.221400e+00, -2.828200e-04, -1.012900e-06, 0.000000e+00}, lThCondLimI = 1.500000e+02, lThCondLimS = 4.300000e+02,
    lSurfTensCorr = 64, lSurfTensCoef = {5.133800e+02, 2.242100e-01, 1.335500e+00, -2.140800e-01, 1.677000e+00, 8.323300e-02}, lSurfTensLimI = 1.756100e+02, lSurfTensLimS = 5.133800e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-5.258950e+01, 2.471690e-01, -4.445610e-04, 3.695940e-07, -1.153720e-10, 0.000000e+00}, lBulkModRLimI = 1.761500e+02, lBulkModRLimS = 4.781500e+02,
    gSatDensCorr = 102, gSatDensCoef = {8.520024e+03, 5.125000e+02, -3.567397e+00, 7.600000e-02, 6.410400e+00, 1.400000e-01, -7.619523e+00, 3.620000e-01, -4.479522e+00, 1.833333e+00, -2.939546e+00, 7.000000e+00, 5.838421e+00, 1.416667e+01}, gSatDensLimI = 1.756100e+02, gSatDensLimS = 5.125000e+02,
    gViscCorr = 110, gViscCoef = {3.065400e-07, 6.965800e-01, 2.048700e+02, 2.430400e+01, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.400000e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {7.836800e-07, 1.756900e+00, 1.081200e+02, -2.110100e+04, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.730000e+02, gThCondLimS = 6.843700e+02); 
 
    constant  FreeFluids.MediaCommon.DataRecord N2(name="Nitrogen", description = "", CAS = "7727-37-9", family = 0,  MW = 2.801340e+01, molarMass=0.0280134, Tc = 1.261900e+02, criticalPressure = 3.390000e+06, Vc = 9.000000e-05, Zc = 2.907460e-01, w = 3.900000e-02, Tb = 7.736000e+01, lnuA = 9.703363e-03, lnuB = -4.978240e+00,
         Cp0Corr = 200, Cp0Coef = {1.039150e+03, 3.296540e+02, 1.713400e+03, -3.422570e+01, 2.446260e+03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 5.000000e+01, Cp0LimS = 5.000000e+03,
      VpCorr = 20, VpCoef = {5.828200e+01, -1.084100e+03, -8.314400e+00, 4.412700e-02, 1.000000e+00, 0.0}, VpLimI = 6.315000e+01, VpLimS = 1.262000e+02,
      HvCorr = 91, HvCoef = {3.971640e+05, 2.372650e+00, -4.840680e+00, 4.484800e+00, -1.556780e+00, 1.262100e+02}, HvLimI = 6.315000e+01, HvLimS = 1.262100e+02, lDensCorr = 47, lDensCoef = {656.007, -2.2046, 5.80173, -5.2628, 1.77008, 126.21}, lDensLimI = 63.15, lDensLimS = 126.21,
      lCpCorr = 16, lCpCoef = {-3.340000e+04, 3.507000e+03, -4.670000e+01, 2.127000e-01, 0.000000e+00, 0.000000e+00}, lCpLimI = 6.400000e+01, lCpLimS = 1.200000e+02,
      lViscCorr = 30, lViscCoef = {3.435800e+00, -2.470600e+01, -2.674800e+00, -4.160300e-05, 2.000000e+00, 0.0}, lViscLimI = 6.315000e+01, lViscLimS = 1.250000e+02,
      lThCondCorr = 51, lThCondCoef = {-2.174300e-01, 1.038300e+01, -1.063100e+00, 3.624500e-04, -2.326500e-05, 0.0}, lThCondLimI = 6.000000e+01, lThCondLimS = 1.240000e+02,
      lSurfTensCorr = 61, lSurfTensCoef = {2.921000e-02, 1.255370e+00, 0.0, 0.0, 0.0, 1.261900e+02}, lSurfTensLimI = 6.305000e+01, lSurfTensLimS = 0.0,
      lBulkModRCorr = 150, lBulkModRCoef = {-4.291280e+00, -6.425800e-03, 5.805370e-05, -7.454690e-08, 3.227420e-11, 0.000000e+00}, lBulkModRLimI = 6.315000e+01, lBulkModRLimS = 1.161500e+02,
      gSatDensCorr = 101, gSatDensCoef = {3.112600e+02, -2.503230e+00, -2.984780e+00, -8.539160e+00, -3.696620e+01, 1.262100e+02}, gSatDensLimI = 6.315000e+01, gSatDensLimS = 1.262100e+02,
      gViscCorr = 110, gViscCoef = {6.559200e-07, 6.081000e-01, 5.471400e+01, 0.0, 0.0, 0.0}, gViscLimI = 6.315000e+01, gViscLimS = 1.970000e+03,
      gThCondCorr = 120, gThCondCoef = {3.314300e-04, 7.722000e-01, 1.632300e+01, 3.737200e+02, 0.0, 0.0}, gThCondLimI = 6.315000e+01, gThCondLimS = 2.000000e+03);
    
    constant  FreeFluids.MediaCommon.DataRecord O2(name="Oxygen", description = "", CAS = "7782-44-7", family = 0, MW = 3.200000e+01, molarMass=0.032, Tc = 1.545800e+02, criticalPressure = 5.043000e+06, Vc = 7.300000e-05, Zc = 2.864150e-01, w = 2.500000e-02, Tb = 8.720000e+01, lnuA = 7.206318e-03, lnuB = -5.283682e+00,
    Cp0Corr = 200, Cp0Coef = {9.097840e+02, 4.331090e+02, 1.282710e+03, -2.590740e+02, 1.913590e+03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 5.000000e+01, Cp0LimS = 3.500000e+03, 
    VpCorr = 25, VpCoef = {5.043000e+06, -6.033340e+00, 1.170520e+00, -1.038500e+00, -1.206520e+00, 1.545900e+02}, VpLimI = 5.015000e+01, VpLimS = 1.545900e+02,
    HvCorr = 91, HvCoef = {2.394190e+05, -5.658560e-01, 2.362100e+00, -2.740890e+00, 1.279100e+00, 1.545900e+02}, HvLimI = 5.435000e+01, HvLimS = 1.545500e+02,
    lDensCorr = 46, lDensCoef = {4.383560e+02, 4.723270e+02, 1.412190e+03, -1.798840e+03, 1.017990e+03, 1.545900e+02}, lDensLimI = 5.015000e+01, lDensLimS = 1.545900e+02,
    lCpCorr = 17, lCpCoef = {1.069570e+04, -4.191670e+02, 7.180700e+00, -5.396550e-02, 1.512440e-04, 0.000000e+00}, lCpLimI = 6.315000e+01, lCpLimS = 1.431500e+02,
    lViscCorr = 30, lViscCoef = {-4.147600e+00, 9.404000e+01, -1.207000e+00, 0.0, 0.0, 0.0}, lViscLimI = 5.436000e+01, lViscLimS = 1.500000e+02,
    lThCondCorr = 50, lThCondCoef = {2.741000e-01, -1.380000e-03, 0.0, 0.0, 0.0, 0.0}, lThCondLimI = 6.000000e+01, lThCondLimS = 1.500000e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {3.797000e-02, 1.210360e+00, 0.0, 0.0, 0.0, 1.546000e+02}, lSurfTensLimI = 5.435000e+01, lSurfTensLimS = 0.0,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.883080e+01, 5.068800e-02, -4.980310e-05, 2.385930e-08, -3.929240e-12, 0.000000e+00}, lBulkModRLimI = 6.315000e+01, lBulkModRLimS = 1.431500e+02,
    gSatDensCorr = 101, gSatDensCoef = {4.383560e+02, -2.460970e+00, -3.278140e+00, -7.420250e+00, -4.001920e+01, 1.545900e+02}, gSatDensLimI = 5.015000e+01, gSatDensLimS = 1.545900e+02,
    gViscCorr = 110, gViscCoef = {8.013400e-07, 6.032100e-01, 5.609000e+01, 1.584900e+03, 0.0, 0.0}, gViscLimI = 5.435000e+01, gViscLimS = 1.950000e+03,
    gThCondCorr = 120, gThCondCoef = {4.499400e-04, 7.456000e-01, 5.669900e+01, 0.0, 0.0, 0.0}, gThCondLimI = 8.000000e+01, gThCondLimS = 2.000000e+03);

  constant FreeFluids.MediaCommon.DataRecord Pentane_n(
    name = "Pentane_n", description = "", CAS = "109-66-0", family = 1, MW = 7.215000e+01, molarMass = 7.215000e-02, Tc = 4.696900e+02, criticalPressure = 3.370000e+06, Vc = 3.110000e-04, Zc = 2.630000e-01, w = 2.510000e-01, Tb = 3.091500e+02, mu = 0.000000e+00, lnuA = 1.342900e-02, lnuB = -5.241821e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 7.220000e+05, 2.510000e+02, 5.015700e+01, 2.998000e+00, -8.770000e+06, 1.760000e+02, 1.230000e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 2.000000e+02, Cp0LimS = 1.500000e+03,
    VpCorr = 25, VpCoef = {3.370000e+06, -7.266130e+00, 1.490250e+00, -2.978460e+00, -1.735520e+00, 4.697000e+02}, VpLimI = 1.441500e+02, VpLimS = 4.697000e+02,
    HvCorr = 91, HvCoef = {7.628270e+05, 2.443590e+00, -5.313080e+00, 5.482680e+00, -2.159360e+00, 4.697000e+02}, HvLimI = 1.461500e+02, HvLimS = 4.651500e+02,
    lDensCorr = 46, lDensCoef = {2.319940e+02, 3.963940e+02, 4.305560e+02, -6.045970e+02, 4.254710e+02, 4.697000e+02}, lDensLimI = 1.441500e+02, lDensLimS = 4.697000e+02,
    lCpCorr = 19, lCpCoef = {5.426700e+01, 2.987200e+03, -2.546090e+03, 1.624410e+02, 1.771460e+03, 4.697000e+02}, lCpLimI = 1.441500e+02, lCpLimS = 4.581500e+02,
    lViscCorr = 30, lViscCoef = {-2.038300e+01, 1.050400e+03, 1.487400e+00, -2.016700e-27, 1.000000e+01, 0.000000e+00}, lViscLimI = 1.434200e+02, lViscLimS = 4.696500e+02,
    lThCondCorr = 50, lThCondCoef = {2.537000e-01, -5.760000e-04, 3.440000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 1.434200e+02, lThCondLimS = 4.696500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.209000e-02, 1.205400e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 4.696600e+02}, lSurfTensLimI = 1.434200e+02, lSurfTensLimS = 4.696600e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-7.767190e+00, 2.086840e-02, 8.012280e-06, -3.875020e-08, 2.772460e-11, 0.000000e+00}, lBulkModRLimI = 1.441500e+02, lBulkModRLimS = 4.581500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.319940e+02, -2.901920e+00, -3.579250e+00, -1.049650e+01, -5.837450e+01, 4.697000e+02}, gSatDensLimI = 1.631500e+02, gSatDensLimS = 4.651500e+02,
    gViscCorr = 110, gViscCoef = {6.341200e-08, 8.475800e-01, 4.171800e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.434200e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {-6.844000e+02, 7.640000e-01, -1.055000e+09, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.731500e+02, gThCondLimS = 1.000000e+03); 

  constant FreeFluids.MediaCommon.DataRecord Propane(
    name = "Propane", description = "", CAS = "74-98-6", family = 1, MW = 4.410000e+01, molarMass=0.0441, Tc = 3.698000e+02, criticalPressure = 4.245518e+06, Vc = 2.030000e-04, Zc = 2.810000e-01, w = 1.520000e-01, Tb = 2.310500e+02, mu = 8.700000e-02, lnuA = 1.449900e-02, lnuB = -5.381233e+00,
    Cp0Corr = 7, Cp0Coef = {4.000000e+00, 0.000000e+00, 0.000000e+00, 3.043000e+00, 3.930000e+02, 5.874000e+00, 1.237000e+03, 9.337000e+00, 1.984000e+03, 7.922000e+00, 4.351000e+03, 0.000000e+00, 0.000000e+00}, Cp0LimI = 8.000000e+01, Cp0LimS = 1.500000e+03,
    VpCorr = 20, VpCoef = {5.907800e+01, -3.492600e+03, -6.066900e+00, 1.091900e-05, 2.000000e+00, 0.000000e+00}, VpLimI = 8.547000e+01, VpLimS = 3.698300e+02,
    HvCorr = 91, HvCoef = {5.041650e+05, -9.494470e-01, 3.787670e+00, -4.335740e+00, 1.837010e+00, 3.698000e+02}, HvLimI = 9.815000e+01, HvLimS = 3.696500e+02,
    lDensCorr = 46, lDensCoef = {2.172410e+02, 4.227030e+02, 1.497270e+02, -1.899320e+02, 2.108340e+02, 3.698000e+02}, lDensLimI = 9.815000e+01, lDensLimS = 3.698000e+02,
    lCpCorr = 19, lCpCoef = {9.074170e+01, 2.657210e+03, -2.568980e+03, 2.548870e+03, -8.598870e+02, 3.698000e+02}, lCpLimI = 9.815000e+01, lCpLimS = 3.631500e+02,
    lViscCorr = 37, lViscCoef = {2.563440e+00, 1.613700e-01, 3.725330e+02, 3.803300e+01, 1.751000e-05, 0.000000e+00}, lViscLimI = 9.000000e+01, lViscLimS = 3.600000e+02,
    lThCondCorr = 50, lThCondCoef = {2.661000e-01, -6.336000e-04, 5.700000e-08, 6.550000e-10, -7.010000e-13, 0.000000e+00}, lThCondLimI = 8.547000e+01, lThCondLimS = 3.698200e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {5.094000e-02, 1.220510e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.698200e+02}, lSurfTensLimI = 8.545000e+01, lSurfTensLimS = 3.698000e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.262750e+01, 5.325320e-02, -6.448980e-05, 2.776990e-08, 8.934950e-12, 0.000000e+00}, lBulkModRLimI = 9.815000e+01, lBulkModRLimS = 3.431500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.172410e+02, -2.688110e+00, -3.419080e+00, -9.319170e+00, -4.978980e+01, 3.698000e+02}, gSatDensLimI = 1.281500e+02, gSatDensLimS = 3.698000e+02,
    gViscCorr = 110, gViscCoef = {4.742200e-08, 9.041600e-01, -4.748400e+00, 4.785700e+02, 0.000000e+00, 0.000000e+00}, gViscLimI = 8.547000e+01, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {-1.120000e+00, 1.097200e-01, -9.834600e+03, -7.535800e+06, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.311100e+02, gThCondLimS = 1.000000e+03); 

  constant FreeFluids.MediaCommon.DataRecord R134A(
    name = "R134A", description = "1,1,1,2-Tetrafluorethane", CAS = "811-97-2", family = 17, MW = 1.020310e+02, molarMass = 1.020310e-01, Tc = 3.7400e+02, criticalPressure=4.056000e+06, Vc = 1.988000e-04, Zc = 2.591786e-01, w = 3.270000e-01, Tb = 2.470500e+02, lnuA = 5.787034e-03, lnuB = -5.099898e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 1.340000e+05, 2.300000e+02, 2.188600e+01, -3.694000e+00, -2.840000e+06, 1.560000e+02, 7.000000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 3.000000e+03,
    VpCorr = 25, VpCoef = {4.056000e+06, -7.540170e+00, 1.307800e+00, -3.313210e+00, -2.633880e+00, 3.741800e+02}, VpLimI = 1.699000e+02, VpLimS = 3.741000e+02,
    HvCorr = 91, HvCoef = {2.580810e+05, -1.410140e+00, 5.186860e+00, -5.842850e+00, 2.412300e+00, 3.741800e+02}, HvLimI = 1.701500e+02, HvLimS = 3.731500e+02,
    lDensCorr = 46, lDensCoef = {5.132340e+02, 9.384090e+02, 7.342780e+02, -9.055000e+02, 7.203790e+02, 3.741800e+02}, lDensLimI = 1.699000e+02, lDensLimS = 3.741000e+02,
    lCpCorr = 19, lCpCoef = {4.350290e+01, 1.340790e+03, -5.439980e+02, -2.457510e+02, 8.908640e+02, 3.741800e+02}, lCpLimI = 1.701500e+02, lCpLimS = 3.531500e+02,
    lViscCorr = 30, lViscCoef = {-5.340380e+01, 2.002990e+03, 7.034950e+00, -2.152880e-05, 2.000000e+00, 0.000000e+00}, lViscLimI = 2.230000e+02, lViscLimS = 3.730000e+02,
    lThCondCorr = 50, lThCondCoef = {3.089520e-01, -1.404510e-03, 3.228500e-06, -3.792470e-09, 5.639790e-13, 0.000000e+00}, lThCondLimI = 1.998000e+02, lThCondLimS = 3.664800e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-9.212050e+00, 1.267010e-02, -1.290250e-06, -2.455060e-09, 9.853880e-13, 0.000000e+00}, lBulkModRLimI = 1.701500e+02, lBulkModRLimS = 3.531500e+02,
    gSatDensCorr = 101, gSatDensCoef = {5.132340e+02, -3.159000e+00, -3.112120e+00, -1.286680e+01, -5.977260e+01, 3.741800e+02}, gSatDensLimI = 1.701500e+02, gSatDensLimS = 3.731500e+02); 

  constant FreeFluids.MediaCommon.DataRecord R410A(
    name = "R410A", description = "", CAS = "R410A.PPF", family = 17, MW = 7.258540e+01, molarMass = 7.258540e-02, Tc = 3.444800e+02, criticalPressure = 4.901200e+06, Vc = 1.581280e-04, Zc = 2.705800e-01, w = 2.960000e-01, Tb = 2.217100e+02, lnuA = 5.986500e-03, lnuB = -5.335073e+00,
    Cp0Corr = 7, Cp0Coef = {0.000000e+00, 2.874900e+00, 1.000000e-01, 2.062300e+00, 6.970009e+02, 5.975100e+00, 1.723000e+03, 1.561200e+00, 3.875006e+03, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = -6.000000e+01, Cp0LimS = 8.000000e+02,
    VpCorr = 25, VpCoef = {4.901200e+06, -7.415000e+00, 1.292920e+00, -2.758960e+00, -2.623060e+00, 3.444940e+02}, VpLimI = 1.731500e+02, VpLimS = 3.444900e+02,
    HvCorr = 91, HvCoef = {4.979860e+05, 1.650880e+00, -3.286540e+00, 3.331220e+00, -1.251180e+00, 3.444940e+02}, HvLimI = 2.031500e+02, HvLimS = 3.444900e+02,
    lDensCorr = 46, lDensCoef = {4.590290e+02, 8.789240e+02, 1.328030e+03, -2.079870e+03, 1.391320e+03, 3.444940e+02}, lDensLimI = 1.731500e+02, lDensLimS = 3.444900e+02,
    lCpCorr = 19, lCpCoef = {6.657140e+01, 1.230010e+03, -5.262890e+01, -9.520020e+02, 2.004840e+03, 3.444940e+02}, lCpLimI = 1.731500e+02, lCpLimS = 3.381500e+02,
    lViscCorr = 30, lViscCoef = {-1.410150e+01, 6.132140e+02, 7.546200e-01, -1.444710e-05, 2.000000e+00, 0.000000e+00}, lViscLimI = 2.231500e+02, lViscLimS = 3.431500e+02,
    lThCondCorr = 50, lThCondCoef = {5.428380e-01, -3.733200e-03, 1.125690e-05, -1.290000e-08, -3.467140e-24, 0.000000e+00}, lThCondLimI = 2.231500e+02, lThCondLimS = 3.431500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-2.057400e+01, 4.908220e-02, -4.428180e-05, 1.989800e-08, -3.343700e-12, 0.000000e+00}, lBulkModRLimI = 1.731500e+02, lBulkModRLimS = 3.291500e+02,
    gSatDensCorr = 101, gSatDensCoef = {4.590290e+02, -2.688940e+00, -4.205650e+00, -1.081900e+01, -5.909060e+01, 3.444940e+02}, gSatDensLimI = 1.731500e+02, gSatDensLimS = 3.431500e+02,
    gViscCorr = 111, gViscCoef = {8.286300e-07, 3.980000e-08, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.431500e+02, gViscLimS = 3.931500e+02,
    gThCondCorr = 121, gThCondCoef = {-8.700000e-03, 7.410000e-05, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 2.431500e+02, gThCondLimS = 3.931500e+02); 

    constant  FreeFluids.MediaCommon.DataRecord ShellS2(name="Shell S2", MW = 3.600000e+02, molarMass=0.36, Tc = 6.129900e+02,
    VpCorr = 23, VpCoef = {2.168690e+01, 6.552780e+03, 1.675470e+01, 0.0, 0.0, 0.0}, VpLimI = 2.730000e+02, VpLimS = 6.280000e+02,
    lDensCorr = 40, lDensCoef = {7.374610e+02, 2.323960e+00, -1.018510e-02, 1.507730e-05, -8.163430e-09, 0.0}, lDensLimI = 2.730000e+02, lDensLimS = 6.130000e+02,
    lCpCorr = 17, lCpCoef = {9.102910e+02, 2.713260e+00, 3.245940e-03, -4.885660e-06, 2.689070e-09, 0.0}, lCpLimI = 2.730000e+02, lCpLimS = 6.130000e+02,
    lTfromHsatCorr = 140, lTfromHsatCoef = {2.466690e+01, 8.381320e-04, -5.270540e-10, 2.754040e-16, -6.340110e-23, 0.0}, lTfromHsatLimI = 2.730000e+02, lTfromHsatLimS = 6.130000e+02,
    lViscCorr = 30, lViscCoef = {-1.581700e+02, 1.000000e+04, 2.153310e+01, -1.153330e-05, 2.000000e+00, 0.0}, lViscLimI = 2.730000e+02, lViscLimS = 6.130000e+02,
    lThCondCorr = 50, lThCondCoef = {1.682120e-01, -1.705060e-04, 2.637920e-07, -2.901120e-10, 1.063780e-13, 0.0}, lThCondLimI = 2.730000e+02, lThCondLimS = 6.130000e+02);

  constant FreeFluids.MediaCommon.DataRecord Toluene(
    name = "Toluene", description = "", CAS = "108-88-3", family = 5, MW = 9.214020e+01, molarMass=0.0921402, Tc = 5.917400e+02, criticalPressure = 4.108000e+06, Vc = 3.160000e-04, Zc = 2.638410e-01, w = 2.630000e-01, Tb = 3.835500e+02, mu = 4.000000e-01, lnuA = 1.076583e-02, lnuB = -5.464864e+00,
    Cp0Corr = 5, Cp0Coef = {4.000000e+00, 8.700000e+04, 2.040000e+02, 3.409500e+01, 1.774300e+01, -3.574000e+07, 2.090000e+02, 7.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 1.500000e+03,
    VpCorr = 24, VpCoef = {4.108000e+06, -7.346970e+00, 1.540940e+00, -3.123090e+00, -2.137540e+00, 5.917500e+02}, VpLimI = 2.031500e+02, VpLimS = 5.917500e+02,
    HvCorr = 90, HvCoef = {5.375200e+07, 5.034100e-01, 2.475500e-01, -7.289800e-01, 3.779400e-01, 5.917500e+02}, HvLimI = 1.781800e+02, HvLimS = 5.850000e+02,
    lDensCorr = 41, lDensCoef = {8.792000e-01, 2.713600e-01, 5.917500e+02, 2.924100e-01, 0.000000e+00, 0.000000e+00}, lDensLimI = 1.781800e+02, lDensLimS = 5.917500e+02,
    lCpCorr = 16, lCpCoef = {1.401400e+05, -1.523000e+02, 6.950000e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lCpLimI = 1.781800e+02, lCpLimS = 5.000000e+02,
    lViscCorr = 37, lViscCoef = {2.954800e+00, 3.000000e-05, 1.048107e+03, 1.371440e+02, 3.970000e-06, 0.000000e+00}, lViscLimI = 1.781800e+02, lViscLimS = 5.880000e+02,
    lThCondCorr = 50, lThCondCoef = {2.038000e-01, -2.353000e-04, -1.900000e-08, 1.000000e-11, 1.300000e-14, 0.000000e+00}, lThCondLimI = 1.781800e+02, lThCondLimS = 4.748500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {6.676000e-02, 1.243860e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.917500e+02}, lSurfTensLimI = 1.781500e+02, lSurfTensLimS = 5.917500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-1.261520e+01, 3.954580e-02, -3.672380e-05, 1.489570e-08, 2.204820e-13, 0.000000e+00}, lBulkModRLimI = 2.331500e+02, lBulkModRLimS = 5.731500e+02,
    gSatDensCorr = 101, gSatDensCoef = {2.915830e+02, -2.743060e+00, -4.283430e+00, -9.579830e+00, -6.301270e+01, 5.917500e+02}, gSatDensLimI = 1.931500e+02, gSatDensLimS = 5.831500e+02,
    gViscCorr = 111, gViscCoef = {-7.109000e-07, 2.788500e-08, -7.300000e-12, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 1.781800e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {2.392000e-05, 1.269400e+00, 5.370000e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gThCondLimI = 3.837800e+02, gThCondLimS = 1.000000e+03); 

  
    constant  FreeFluids.MediaCommon.DataRecord Water(name="Water", CAS = "7732-18-5", family = 6, MW = 1.801500e+01, molarMass=0.018015, Tc = 6.470950e+02, criticalPressure = 2.206400e+07, Vc = 0.0000559, Zc = 2.294000e-01, w = 3.440000e-01, mu=1.8, Tb = 3.731500e+02, IsothComp = 4.900000e-10, lnuA = 7.847711e-03, lnuB = -5.115707e+00, 
    Cp0Corr = 200, Cp0Coef = {1.850940e+03, 1.651260e+03, 2.890990e+03, 5.967120e+02, 1.248300e+03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 5.000000e+01, Cp0LimS = 4.500000e+03, 
    VpCorr = 25, VpCoef = {2.206400e+07, -7.773360e+00, 1.463190e+00, -2.746930e+00, -1.349700e+00, 6.470960e+02}, VpLimI = 2.731500e+02, VpLimS = 6.470960e+02,
    HvCorr = 90, HvCoef = {5.205300e+07, 3.199000e-01, -2.120000e-01, 2.580000e-01, 0.0, 6.473500e+02}, HvLimI = 2.731600e+02, HvLimS = 6.473500e+02,
    lDensCorr = 45, lDensCoef = {3.251621e+01, -3.213004e+00, 7.924110e+00, -7.359898e+00, 2.703522e+00, 6.470960e+02}, lDensLimI = 2.531000e+02, lDensLimS = 6.470960e+02,
    lCpCorr = 19, lCpCoef = {2.487530e+02, 3.588150e+03, -5.485740e+01, -6.383800e+02, 2.312590e+03, 6.470960e+02}, lCpLimI = 2.731500e+02, lCpLimS = 6.371500e+02,
    lTfromHsatCorr = 0, lTfromHsatCoef = {5.580150e+02, 1.852530e-04, -6.677220e-11, -3.939070e-17, -9.552730e-24, 0.000000e+00}, lTfromHsatLimI = 2.731500e+02, lTfromHsatLimS = 6.371500e+02,
    lViscCorr = 30, lViscCoef = {-5.196400e+01, 3.670600e+03, 5.733100e+00, -5.349500e-29, 1.000000e+01, 0.0}, lViscLimI = 2.731500e+02, lViscLimS = 6.431500e+02,
    lThCondCorr = 50, lThCondCoef = {-4.320000e-01, 5.725500e-03, -8.078000e-06, 1.861000e-09, 0.0, 0.0}, lThCondLimI = 2.731600e+02, lThCondLimS = 6.331500e+02,
    lSurfTensCorr = 61, lSurfTensCoef = {1.855000e-01, 2.717000e+00, -3.554000e+00, 2.047000e+00, 0.0, 6.473500e+02}, lSurfTensLimI = 2.731600e+02, lSurfTensLimS = 6.473500e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-5.511280e+01, 2.352350e-01, -3.880190e-04, 2.953210e-07, -8.462820e-11, 0.000000e+00}, lBulkModRLimI = 2.730000e+02, lBulkModRLimS = 5.830000e+02,
    gSatDensCorr = 101, gSatDensCoef = {3.222720e+02, -3.490130e+00, -3.010570e+00, -1.273460e+01, -5.467520e+01, 6.470960e+02}, gSatDensLimI = 2.731500e+02, gSatDensLimS = 6.470960e+02,
    gViscCorr = 110, gViscCoef = {7.002327e-08, 9.345760e-01, 1.956338e+02, -1.304599e+04, 0.0, 0.0}, gViscLimI = 2.731600e+02, gViscLimS = 1.073150e+03, gThCondCorr = 120, gThCondCoef = {6.598600e-06, 1.394700e+00, 5.947800e+01, -1.548400e+04, 0.0, 0.0}, gThCondLimI = 2.731600e+02, gThCondLimS = 1.073150e+03);

  constant FreeFluids.MediaCommon.DataRecord Xylene_m(
    name = "Xylene_m", description = "", CAS = "108-38-3", family = 5, MW = 1.061670e+02, molarMass = 1.061670e-01, Tc = 6.168900e+02, criticalPressure = 3.534600e+06, Vc = 3.750000e-04, Zc = 2.588436e-01, w = 3.250000e-01, Tb = 4.122700e+02, mu = 3.000000e-01, lnuA = 1.086806e-02, lnuB = -5.270932e+00,
    Cp0Corr = 7, Cp0Coef = {2.169909e+00, 0.000000e+00, 0.000000e+00, 4.443120e+00, 1.600000e+02, 2.862794e+00, 1.900000e+02, 2.483298e+01, 1.333000e+03, 1.626077e+01, 3.496000e+03, 0.000000e+00, 0.000000e+00}, Cp0LimI = 5.000000e+01, Cp0LimS = 2.000000e+03,
    VpCorr = 20, VpCoef = {8.509900e+01, -7.615900e+03, -9.307200e+00, 5.564300e-06, 2.000000e+00, 0.000000e+00}, VpLimI = 2.253000e+02, VpLimS = 6.170000e+02,
    HvCorr = 91, HvCoef = {3.515380e+05, -2.297970e+00, 7.660880e+00, -8.661040e+00, 3.611550e+00, 6.168900e+02}, HvLimI = 2.253000e+02, HvLimS = 6.166500e+02,
    lDensCorr = 46, lDensCoef = {2.830000e+02, 6.700403e+02, -3.007884e+02, 6.136168e+02, -1.893359e+02, 6.170500e+02}, lDensLimI = 225.3e+00, lDensLimS = 6.170500e+02,
    lCpCorr = 19, lCpCoef = {3.921800e+01, 2.688330e+03, -2.003000e+03, -6.686670e+02, 1.211730e+03, 6.168900e+02}, lCpLimI = 2.253000e+02, lCpLimS = 6.166500e+02,
    lViscCorr = 31, lViscCoef = {-3.820000e+00, 1.027000e+03, -6.380000e-04, 4.520000e-07, 0.000000e+00, 0.000000e+00}, lViscLimI = 2.260000e+02, lViscLimS = 6.130000e+02,
    lThCondCorr = 50, lThCondCoef = {1.643000e-01, -1.466000e-05, -2.387000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lThCondLimI = 2.250000e+02, lThCondLimS = 6.030000e+02,
    lSurfTensCorr = 64, lSurfTensCoef = {6.168900e+02, 6.445000e-02, 1.256000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, lSurfTensLimI = 2.253000e+02, lSurfTensLimS = 6.168900e+02,
    lBulkModRCorr = 150, lBulkModRCoef = {-8.532470e+00, 2.063220e-02, 0.000000e+00, -1.936320e-08, 1.308310e-11, 0.000000e+00}, lBulkModRLimI = 2.253000e+02, lBulkModRLimS = 6.080000e+02,
    gSatDensCorr = 102, gSatDensCoef = {2.665000e+03, 6.168900e+02, -7.808901e-02, 5.400000e-02, -3.249336e+00, 4.300000e-01}, gSatDensLimI = 2.253000e+02, gSatDensLimS = 6.168900e+02,
    gViscCorr = 110, gViscCoef = {6.829300e-07, 5.219900e-01, 3.241700e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00}, gViscLimI = 2.253000e+02, gViscLimS = 1.000000e+03,
    gThCondCorr = 120, gThCondCoef = {3.059300e-09, 2.418200e+00, -5.692800e+02, 1.210600e+05, 0.000000e+00, 0.000000e+00}, gThCondLimI = 3.200000e+02, gThCondLimS = 1.000000e+03); 
    annotation(
      Documentation(info = "<html>
  <body>
  <p>Contains data of the individual substances, mainly parameters for temperature dependent correlations for physical properties. The records can be used both with the ideal gas medium packages, and with the temperature dependent liquid packages.</p>

  <p><strong>There are data records for the following substances:</strong></p>
  <p>MarlothermSH</p>
  <p>Methane</p> 
  <p>Methanol</p> 
  <p>N2 (Nitrogen)</p>
  <p>O2 (Oxygen)</p>
  <p>Pentane_n</p>
  <p>Propane</p>
  <p>R134A (1,1,1,2-Tetrafluorethane)</p>
  <p>R410A</p> 
  <p>ShellS2</p> 
  <p>Toluene</p>
  <p>Water</p> 
  <p>Xylene-m</p>
  </body>
  </html>"));
  end MediaDataMZ;

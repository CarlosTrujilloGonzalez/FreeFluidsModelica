within FreeFluids;

package MediaCommon "MediaCommon.mo by Carlos Trujillo
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
  record DataRecord
    String name = "", description = "", CAS = "";
    Integer family = 0 "6=water, 7=alcohol, 8=polyol, 17=haloalkane, 18=haloalkene";
    Real MW = 0.0, molarMass = 0.0, Tc = 0.0, criticalPressure = 0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 6.667e-10, lnuA = 0.0, lnuB = 0.0;
    Integer Cp0Corr = 0, VpCorr = 0, BtCorr = 0, HvCorr = 0, lDensCorr = 0, lCpCorr = 0, lTfromHsatCorr = 0, lViscCorr = 0, lThCondCorr = 0, lSurfTensCorr = 0, lBulkModRCorr = 0, gSatDensCorr = 0, gViscCorr = 0, gThCondCorr = 0;
    Real Cp0Coef[:] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, VpCoef[:] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, BtCoef[6] = {0, 0, 0, 0, 0, 0}, HvCoef[6] = {0, 0, 0, 0, 0, 0}, lDensCoef[:] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, lCpCoef[6] = {0, 0, 0, 0, 0, 0}, lTfromHsatCoef[6] = {0, 0, 0, 0, 0, 0}, lViscCoef[6] = {0, 0, 0, 0, 0, 0}, lThCondCoef[6] = {0, 0, 0, 0, 0, 0}, lSurfTensCoef[6] = {0, 0, 0, 0, 0, 0}, lBulkModRCoef[6] = {0, 0, 0, 0, 0, 0}, gSatDensCoef[:] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, gViscCoef[6] = {0, 0, 0, 0, 0, 0}, gThCondCoef[6] = {0, 0, 0, 0, 0, 0};
    Real Cp0LimI = 0, VpLimI = 0, BtLimI = 0, HvLimI = 0, lDensLimI = 0, lCpLimI = 0, lTfromHsatLimI = 0, lViscLimI = 0, lThCondLimI = 0, lSurfTensLimI = 0, lBulkModRLimI = 0, gSatDensLimI = 0, gViscLimI = 0, gThCondLimI = 0;
    Real Cp0LimS = 0, VpLimS = 0, BtLimS = 0, HvLimS = 0, lDensLimS = 0, lCpLimS = 0, lTfromHsatLimS = 0, lViscLimS = 0, lThCondLimS = 0, lSurfTensLimS = 0, lBulkModRLimS = 0, gSatDensLimS = 0, gViscLimS = 0, gThCondLimS = 0;
  end DataRecord;

  constant DataRecord MediaDataTemplate(name = "", description = "", CAS = "", family = 0, MW = 0.0, molarMass = 0.0, Tc = 0.0, criticalPressure = 0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 0.0, lnuA = 0.0, lnuB = 0.0, Cp0Corr = 0, Cp0Coef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 0.0, Cp0LimS = 0.0, VpCorr = 0, VpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, VpLimI = 0.0, VpLimS = 0.0, BtCorr = 0, BtCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, BtLimI = 0.0, BtLimS = 0.0, HvCorr = 0, HvCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, HvLimI = 0.0, HvLimS = 0.0, lDensCorr = 0, lDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lDensLimI = 0.0, lDensLimS = 0.0, lCpCorr = 0, lCpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lCpLimI = 0.0, lCpLimS = 0.0, lTfromHsatCorr = 0, lTfromHsatCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lTfromHsatLimI = 0.0, lTfromHsatLimS = 0.0, lViscCorr = 0, lViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lViscLimI = 0.0, lViscLimS = 0.0, lThCondCorr = 0, lThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lThCondLimI = 0.0, lThCondLimS = 0.0, lSurfTensCorr = 0, lSurfTensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lSurfTensLimI = 0.0, lSurfTensLimS = 0.0, lBulkModRCorr = 0, lBulkModRCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lBulkModRLimI = 0.0, lBulkModRLimS = 0.0, gSatDensCorr = 0, gSatDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gSatDensLimI = 0.0, gSatDensLimS = 0.0, gViscCorr = 0, gViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gViscLimI = 0.0, gViscLimS = 0.0, gThCondCorr = 0, gThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gThCondLimI = 0.0, gThCondLimS = 0.0);

  package Types
    type InputChoice = enumeration(pT "(p,T) as inputs", ph "(p,h) as inputs", ps "(p,s) as inputs", dT "(d,T) as inputs");
    type ReferenceState = enumeration(ASHRAE "0 enthalpy and entropy for saturated liquid at -40ºC", NBP "0 enthalpy and entropy for saturated liquid at normal boiling point", IIR "H=200 KJ/kg, and S=1.0 KJ/(kg·K) for saturated liquid at 0ºC", User "User defined reference_T and reference_P", None "No reference is used");
    type CorrelationEquation = enumeration(FF_DIPPR100, FF_Polynomial, FF_Polynomial2, FF_DIPPR100Ld, FF_expDIPPR100, FF_DIPPR101, FF_DIPPR101Vp, FF_DIPPR101Lv, FF_logDIPPR101, FF_DIPPR102, FF_DIPPR103, FF_DIPPR104, FF_DIPPR105, FF_DIPPR106, FF_DIPPR106Hv, FF_DIPPR106Ld, FF_DIPPR106SurfT, FF_DIPPR107, FF_DIPPR107Cp, FF_DIPPR114, FF_DIPPR115, FF_DIPPR116, FF_DIPPR116Ld, FF_Wilhoit, FF_Cooper, FF_Jaechske, FF_ChemSep16, FF_Antoine1, FF_Antoine2, FF_Wagner25, FF_Wagner36, FF_PPDS9, FF_PPDS10, FF_PCWIN, FF_Rackett, FF_ExtAndrade1, FF_ExtAndrade2, FF_ChericVisc, FF_WagnerGd, FF_Tait, FF_ExtWagner, FF_PPDS15);
  end Types;

  package Functions
  "Contains the functions that could be common to several media packages"
    partial function CorrelationSolver "Anderson-Bjork modification of Regula Falsi method"
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real y;
      input Real max=1000;
      input Real min=50;
      input Real tol=1.0e-7;
      output Real x;
    protected
      replaceable function f
      end f;
      Real xMin,xMax,yMin,yMax,fx;
      Integer i,side;
    algorithm
      xMin:=min;
      xMax:=max;
      yMin:=f(cor,coef,MW,xMin);
      yMax:=f(cor,coef,MW,xMax);
      side:=0;
      for i in 1:20 loop
        x:=xMin+(xMax-xMin)*(y-yMin)/(yMax-yMin);
        if x<min then
          x:=min;
        end if;
        fx:=f(cor,coef,MW,x);
        if (abs(y-fx)<abs(y*tol)) then break;    
        elseif ((y-fx)*(y-yMax)>0)then//the new point is at the same side than yMax
          if(side==1) then
            if(((y-fx)/(y-yMax))<0.5) then
              yMin:=yMin+(y-yMin)*(y-fx)/(y-yMax);
            else
              yMin:=yMin+0.5*(y-yMin);
            end if;
          end if;
          xMax:=x;
          yMax:=fx;
          side:=1;   
        else
          if(side==-1) then
            if(((y-fx)/(y-yMin))<0.5) then
              yMax:=yMax+(y-yMax)*(y-fx)/(y-yMin);
            else
              yMax:=yMax+0.5*(y-yMax);
            end if;
          end if;
          xMin:=x;
          yMin:=fx;
          side:=-1;
        end if;
      end for;
      assert(abs(y-fx)<abs(y*tol), String(x)+" as root for y="+String(y)+" has a relative error of: "+String((y-fx)/y), AssertionLevel.warning);
      annotation(
        Documentation(info = "<html>
        <body>
        <p>It uses the Anderson-Bjork modification of Regula Falsi method for the inverse solving of a correlation. The correlation (f) is defined as a replaceable function. You must extend the CorrelationSolver redeclaring f as the function you want to use. This function must return a Real value, and must accept the calling parameters that will be used by the solver.</p>
        </body>
        </html>"));
    end CorrelationSolver;
    
    function PhysPropCorr
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real y;
    protected
      Real Tm;
      Integer j;
    algorithm
      if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 or cor == 40 or cor == 49 or cor == 50 or cor == 60 or cor == 70 or cor == 71 or cor == 80 or cor == 82 or cor == 100 or cor == 111 or cor == 121 or cor == 140 or cor == 150 then
        y := coef[1] + x * (coef[2] + x * (coef[3] + x * (coef[4] + x * coef[5]))) "FF_DIPPR100.ok";
      elseif cor == 63 then
        Tm := x - 273.15;
        y := coef[1] + Tm * (coef[2] + Tm * (coef[3] + Tm * (coef[4] + Tm * coef[5]))) "FF_DIPPR100 in deg C";
      elseif cor == 10 then
        y := coef[1] + x * (coef[2] + x * (coef[3] + x * (coef[4] + x * (coef[5] + x * coef[6])))) "FF_Polynomial";
      elseif cor == 130 then
        y := coef[1] + coef[2] * x ^ 0.125 + coef[3] * x ^ 0.25 + coef[4] * x ^ 0.5 + coef[5] * x "FF_Polynomial2.ok";
      elseif cor == 20 or cor == 21 or cor == 30 then
        y := exp(coef[1] + coef[2] / x + coef[3] * log(x) + coef[4] * x ^ coef[5]) "FF_DIPPR101.ok";
      elseif cor == 81 or cor == 110 or cor == 120 then
        y := coef[1] * x ^ coef[2] / (1 + coef[3] / x + coef[4] / x ^ 2) "FF_DIPPR102.ok";
      elseif cor == 41 or cor == 42 then
        y := coef[1] / coef[2] ^ (1 + (1 - x / coef[3]) ^ coef[4]) "FF_DIPPR105";
      elseif cor == 45 or cor == 47 or cor == 61 or cor == 90 or cor == 91 then
        Tm := x / coef[6];
        y := coef[1] * (1 - Tm) ^ (coef[2] + coef[3] * Tm + coef[4] * Tm ^ 2 + coef[5] * Tm ^ 3) "FF_DIPPR106.ok";
      elseif cor == 3 or cor == 4 or cor == 200 then
        y := coef[1] + coef[2] * (coef[3] / x / sinh(coef[3] / x)) ^ 2 + coef[4] * (coef[5] / x / cosh(coef[5] / x)) ^ 2 "FF_DIPPR107.ok";
      elseif cor == 5 then
        Tm := (x - coef[8]) / (x + coef[7]);
        y := R * (coef[1] + coef[2] / x ^ 2 * exp(-coef[3] / x) + coef[4] * Tm ^ 2 + (coef[5] - coef[6] / (x - coef[8]) ^ 2) * Tm ^ 8) "FF_Wilhoit.ok";
      elseif cor == 7 then
        y := coef[1] + coef[2] * x ^ coef[3] "FF_Cooper.ok";
        for j in {4, 6, 8, 10} loop
          if coef[j] > 0 then
            y := y + coef[j] * (coef[j + 1] / x) ^ 2 * exp(coef[j + 1] / x) / (exp(coef[j + 1] / x) - 1) ^ 2;
          end if;
        end for;
        y := y * R;
      elseif cor == 8 then
        y := 1 + coef[1] "FF_Jaechske.ok";
        for j in {2, 6} loop
          if coef[j] > 0 then
            y := y + coef[j] * (coef[j + 1] / x / sinh(coef[j + 1] / x)) ^ 2 + coef[j + 2] * (coef[j + 3] / x / cosh(coef[j + 3] / x)) ^ 2;
          end if;
        end for;
        y := y * R;
      elseif cor == 9 or cor == 18 or cor == 51 or cor == 62 then
        y := coef[1] + exp(coef[2] / x + coef[3] + coef[4] * x + coef[5] * x ^ 2) "FF_ChemSep16";
      elseif cor == 36 then
        y := 10 ^ (coef[1] - coef[2] / (x + coef[3])) "FF_Antoine1";
      elseif cor == 22 then
        Tm := x - 273.15;
        y := 10 ^ (coef[1] - coef[2] / (Tm + coef[3])) "FF_Antoine1. In deg C";
      elseif cor == 23 then
        y := exp(coef[1] - coef[2] / (x + coef[3])) "FF_Antoine2";
      elseif cor == 24 then
        Tm := 1 - x / coef[6];
        y := coef[1] * exp((coef[2] * Tm + coef[3] * Tm ^ 1.5 + coef[4] * Tm ^ 2.5 + coef[5] * Tm ^ 5.0) / (1 - Tm)) "FF_Wagner25";
      elseif cor == 25 then
        Tm := 1 - x / coef[6];
        y := coef[1] * exp((coef[2] * Tm + coef[3] * Tm ^ 1.5 + coef[4] * Tm ^ 3 + coef[5] * Tm ^ 6) / (1 - Tm)) "FF_Wagner36.ok";
      elseif cor == 26 or cor== 102 then
        Tm := 1 - x / coef[2];
        y := coef[1] * exp(coef[2] * (coef[3] * Tm ^ coef[4] + coef[5] * Tm ^ coef[6] + coef[7] * Tm ^ coef[8] + coef[9] * Tm ^ coef[10] + coef[11] * Tm ^ coef[12] + coef[13] * Tm ^ coef[14]) / x) "FF_ExtWagner";
      elseif cor== 103 then
        Tm := 1 - x / coef[2];
        y := coef[1] * exp(coef[3] * Tm ^ coef[4] + coef[5] * Tm ^ coef[6] + coef[7] * Tm ^ coef[8] + coef[9] * Tm ^ coef[10] + coef[11] * Tm ^ coef[12] + coef[13] * Tm ^ coef[14]) "FF_ExtWagner2";    
      elseif cor == 31 or cor == 33 then
        y := exp(coef[1] + coef[2] / x + coef[3] * x + coef[4] * x ^ 2 + coef[5] * x ^ 3) "FF_ExtAndrade1";
      elseif cor == 32 then
        y := 10 ^ (coef[1] + coef[2] / x + coef[3] * x + coef[4] * x ^ 2 + coef[5] * x ^ 3) "FF_ExtAndrade2";
      elseif cor == 37 then
        y := coef[5] * exp(coef[1] * ((coef[3] - x) / (x - coef[4])) ^ 0.33333 + coef[2] * ((coef[3] - x) / (x - coef[4])) ^ 1.33333) "FF_PPDS9";
      elseif cor == 19 or cor == 151 then
        Tm := 1 - x / coef[6];
        y := coef[1] / Tm + coef[2] + coef[3] * Tm + coef[4] * Tm * Tm + coef[5] * Tm * Tm * Tm "FF_PPDS15.ok";
      elseif cor == 43 then
        y := coef[1] + coef[2] * (1 - x / coef[5]) + coef[3] * log(1 - x / coef[5]) + coef[4] * (1 - x / coef[5]) ^ 3 "FF_PCWIN";
      elseif cor == 48 then
        y := coef[1] / coef[2] ^ ((1 - x / coef[3]) ^ coef[4]) "FF_Rackett";
      elseif cor == 44 or cor == 46 then
        Tm := 1 - x / coef[6];
        y := coef[1] + coef[2] * Tm ^ 0.35 + coef[3] * Tm ^ 0.666667 + coef[4] * Tm + coef[5] * Tm ^ 1.333333 "FF_DIPPR116";
      elseif cor == 101 then
        Tm := 1 - x / coef[6];
        y := coef[1] * exp(coef[2] * Tm ^ 0.4 + coef[3] * Tm + coef[4] * Tm ^ 2.1 + coef[5] * Tm ^ 5.6) "FF_WagnerGd.ok";
      elseif cor == 240 then
        Tm := x - 273.15;
        y := (coef[1] + coef[2] * Tm + coef[3] * Tm ^ 2) * (1 - 0.0894 * log(1 + 1e5 / (coef[4] * exp(-coef[5] * Tm)))) "FF_Tait at 1 bar. T in centigrades";
      elseif cor == 241 then  
        Tm := 1-x/coef[2];
        y := coef[1]*(1+coef[3]*Tm^coef[4]+coef[5]*Tm^coef[6]+coef[7]*Tm^coef[8]+coef[9]*Tm^coef[10]+
                    coef[11]*Tm^coef[12]+coef[13]*Tm^coef[14])"FF_ExtLogWagner";
      end if;
      if cor == 1 or cor == 21 or cor == 48 then
        y := y * 1e3;
      elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
        y := y / MW * 1e3;
      elseif cor == 3 or cor == 10 then
        y := y / MW * 1e3 * 4.1868;
      elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 or cor == 80 or cor == 81 or cor == 90 then
        y := y / MW;
      elseif cor == 41 or cor == 44 or cor == 45 or cor == 70 then
        y := y * MW;
      elseif cor == 22 then
        y := y * 133.32239;
      elseif cor == 31 or cor == 32 or cor == 34 or cor == 36 or cor == 63 then
        y := y * 1e-3;
      elseif cor == 33 or cor == 102 or cor== 103 or cor == 241 then
        y := y * MW * 1e-3;
      elseif cor == 49 then
        y := MW * 1e3 / y;
      elseif cor == 240 then
        y := 1 / y;
      end if;
      annotation(
        Inline = true,
        smoothOrder = 2,
        Documentation(info = "<html>
        <body>
        <p>Provides the calculations of physical properties as function of one independent variable (normally T). It needs as input: the number of the correlation to use, the coefficients for the correlation, the molecular weight of the substance, and the value of the independent variable. It returns the value of the physical property. It has an equivalent written in C inside the FFphysprop.c file in the Resources directory.</p>
        </body>
        </html>"));
    
    end PhysPropCorr;

    function PhysPropCorrInv "Performs the inverse solving of the PhysPropCorr function"
      extends CorrelationSolver(redeclare function f=PhysPropCorr);
    end PhysPropCorrInv;
  
    function SpecificEnthalpyCorr
    "Provides the calculation of the specific enthalpy as function of T. You can use correlations using ideal gas constant pressure heat capacity, or liquid heat capacity"
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real H "Enthalpy J/kg";
    protected
      Real y, y2, y4, h, a7_a6, a7_a6_2, a7_a6_4, x1, z, w, s;
      Real Tm, Tm2;
      Integer j;
    algorithm
      if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 then
        H := x * (coef[1] + x * (coef[2] / 2 + x * (coef[3] / 3 + x * (coef[4] / 4 + x * coef[5] / 5)))) "DIPPR 100";
      elseif cor == 10 then
        H := x * (coef[1] + x * (coef[2] / 2 + x * (coef[3] / 3 + x * (coef[4] / 4 + x * coef[5] / 5)))) "Polynomial";
      elseif cor == 3 or cor == 4 or cor == 200 then
        H := coef[1] * x + coef[2] * coef[3] * (1 / tanh(coef[3] / x)) - coef[4] * coef[5] * tanh(coef[5] / x) "DIPPR107.ok";
      elseif cor == 5 then
        a7_a6 := coef[8] / coef[7];
        a7_a6_2 := a7_a6 * a7_a6;
        a7_a6_4 := a7_a6_2 * a7_a6_2;
        x1 := (coef[5] * coef[8] * coef[8] - coef[6]) / (coef[7] * coef[7]);
        y := if x <= coef[8] then 0 else (x - coef[8]) / (x + coef[7]);
        y2 := y * y;
        y4 := y2 * y2;
        h := if x <= coef[8] then 0 else (coef[7] + coef[8]) * ((2 * coef[4] + 8 * coef[5]) * log(1 - y) + (coef[4] * (1 + 1 / (1 - y)) + coef[5] * (7 + 1 / (1 - y))) * y + coef[5] * (3 * y2 + 5 * y * y2 / 3 + y4 + 0.6 * y4 * y + y4 * y2 / 3) + (coef[5] - coef[6] / (coef[7] + coef[8]) ^ 2) * y4 * y2 * y / 7);
        H := R * x * (coef[1] + coef[2] * exp(-coef[3] / x) / (coef[3] * x)) + R * h "Wilhoit.ok";
      elseif cor == 7 then
        H := coef[1] * x + coef[2] / (coef[3] + 1) * x ^ (coef[3] + 1);
        for j in {4, 6, 8, 10} loop
          if coef[j] > 0 then
            H := H + coef[j] * coef[j + 1] / (exp(coef[j + 1] / x) - 1);
          end if;
        end for;
        H := H * R "Cooper.ok";
      elseif cor == 8 then
        H := (1 + coef[1]) * x;
        for j in {2, 6} loop
          if coef[j] > 0 then
            H := H + 2 * coef[j] * coef[j + 1] / (exp(2 * coef[j + 1] / x) - 1) + 2 * coef[j + 2] * coef[j + 3] / (exp(2 * coef[j + 3] / x) + 1);
          end if;
        end for;
        H := H * R "Jaeschke.ok";
      elseif cor == 19 then
        Tm := 1 - x / coef[6];
        Tm2 := Tm * Tm;
        H := -coef[6] * (coef[1] * log(Tm) + coef[2] * Tm + 0.5 * coef[3] * Tm2 + coef[4] * Tm2 * Tm / 3 + 0.25 * coef[5] * Tm2 * Tm2) "PPDS15.ok";
      end if;
      if cor == 1 then
        H := H * 1e3;
      elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
        H := H * 1e3 / MW;
      elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 then
        H := H / MW;
      elseif cor == 3 or cor == 10 then
        H := H * 4.1868 * 1e3 / MW;
      end if;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end SpecificEnthalpyCorr;

    function SpecificEnthalpyCorrInv
    "Performs the inverse solving of the SpecificEnthalpyCorr function"
      extends CorrelationSolver(redeclare function f=SpecificEnthalpyCorr);
    end SpecificEnthalpyCorrInv;
  
    function SpecificEntropyCorr
    "Provides the calculation of the specific entropy as function of T. You can use correlations using ideal gas constant pressure heat capacity, or liquid heat capacity"
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real S "Entropy J/(kg·K)";
    protected
      Real y, y2, y4, h, a7_a6, a7_a6_2, a7_a6_4, x1, z, w, s;
      Real Tm, Tm2;
      Integer j;
    algorithm
      if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 then
        S := coef[1] * log(x) + x * (coef[2] + x * (coef[3] / 2 + x * (coef[4] / 3 + x * coef[5] / 4))) "DIPPR100";
      elseif cor == 10 then
        S := coef[1] * log(x) + x * (coef[2] + x * (coef[3] / 2 + x * (coef[4] / 3 + x * (coef[5] / 4 + x * coef[6] / 5)))) "Polynomial";
      elseif cor == 3 or cor == 4 or cor == 200 then
        S := coef[1] * log(x) + coef[2] * (coef[3] / x / tanh(coef[3] / x) - log(sinh(coef[3] / x))) - coef[4] * (coef[5] / x * tanh(coef[5] / x) - log(cosh(coef[5] / x))) "DIPPR107.ok";
      elseif cor == 5 then
        a7_a6 := coef[8] / coef[7];
        a7_a6_2 := a7_a6 * a7_a6;
        a7_a6_4 := a7_a6_2 * a7_a6_2;
        x1 := (coef[5] * coef[8] * coef[8] - coef[6]) / (coef[7] * coef[7]);
        y := if x <= coef[8] then 0 else (x - coef[8]) / (x + coef[7]);
        y2 := y * y;
        y4 := y2 * y2;
        if x <= coef[8] then
          s := 0;
        else
          z := x * (coef[8] + coef[7]) / ((x + coef[7]) * coef[8]);
          w := 0;
          for j in 1:7 loop
            w := w + (x1 * (-a7_a6) ^ (6 - j) - coef[5]) * y ^ j / j;
          end for;
          s := (coef[4] + (coef[5] * coef[8] * coef[8] - coef[6]) * a7_a6_4 / (coef[7] * coef[7])) * a7_a6_2 * log(z) + (coef[4] + coef[5]) * log((x + coef[7]) / (coef[7] + coef[8])) - (coef[4] * (coef[7] + coef[8]) / coef[7] + coef[6] * y4 * y2 / (7. * coef[8] * (coef[7] + coef[8]))) * y + w;
        end if;
        S := R * (coef[1] * log(x) + coef[2] * (1 + coef[3] / x) * exp(-coef[3] / x) / (coef[3] * coef[3]) + s) "Wilhoit.ok";
      elseif cor == 7 then
        S := if coef[2] > 0 then coef[1] * log(x) + coef[2] / coef[3] * x ^ coef[3] else coef[1] * log(x);
        for j in {4, 6, 8, 10} loop
          if coef[j] > 0 then
            S := S + coef[j] * coef[j + 1] * (exp(coef[j + 1] / x) / (x * (exp(coef[j + 1] / x) - 1)) - log(exp(coef[j + 1] / x) - 1) / coef[j + 1]) "Cooper.ok";
          end if;
        end for;
        S := S * R;
      elseif cor == 8 then
        S := (1 + coef[1]) * log(x);
        for j in {2, 6} loop
          if coef[j] > 0 then
            S := S + coef[j] * (coef[j + 1] / x / tanh(coef[j + 1] / x) - log(sinh(coef[j + 1] / x))) - coef[j + 2] * (coef[j + 3] / x * tanh(coef[j + 3] / x) - log(cosh(coef[j + 3] / x))) "Jaeschke.ok";
          end if;
        end for;
        S := S * R;
      elseif cor == 19 then
        Tm2 := coef[6] * coef[6];
        S := (-(6 * Tm2 * x * (coef[3] + 2 * coef[4] + 3 * coef[5]) - 3 * coef[6] * x * x * (3 * coef[5] + coef[4]) + 2 * x * x * x * coef[5]) / (6 * Tm2 * coef[6])) + (coef[1] + coef[2] + coef[3] + coef[4] + coef[5]) * log(x) - coef[1] * log(abs(x - coef[6]));
      end if;
//for (j=1;j<9;j=j+4) if (coef[j]>0) *S = *S+coef[j]*(coef[j+1]/x/tanh(coef[j+1]/x)-log(sinh(coef[j+1]/x)))-coef[j+2]*(coef[j+3]/x*tanh(coef[j+3]/x)-log(cosh(coef[j+3]/x)));
      if cor == 1 then
        S := S * 1e3;
      elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
        S := S * 1e3 / MW;
      elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 then
        S := S / MW;
      elseif cor == 3 or cor == 10 then
        S := S * 4.1868 * 1e3 / MW;
      end if;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end SpecificEntropyCorr;

    function SpecificEntropyCorrInv
    "Performs the inverse solving of the SpecificEntropyCorr function"
      extends CorrelationSolver(redeclare function f=SpecificEntropyCorr);
    end SpecificEntropyCorrInv;
    
      function SpecificEnthalpyCorr2
      "Receives a DataRecord record as input instead of correlation parameters. Calls the SpecificEnthalpyCorr to performs the calculation of the ideal gas specific enthalpy"
        input FreeFluids.MediaCommon.DataRecord data;
        input Real T;
        output Real h;
      
    algorithm
        h:=SpecificEnthalpyCorr(data.Cp0Corr, data.Cp0Coef, data.MW, T);
      
    end SpecificEnthalpyCorr2;
  
    function SpecificEntropyCorr2
    "Receives a DataRecord record as input instead of correlation parameters. Calls the SpecificEntropyCorr to performs the calculation of the ideal gas specific entropy"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real s;
    algorithm
      s:=SpecificEntropyCorr(data.Cp0Corr, data.Cp0Coef, data.MW, T);
    end SpecificEntropyCorr2;

    function Cp0Corr "Calculates ideal gas heat capacity from a given DataRecord and a given temperature."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real y;
    algorithm
      y:=PhysPropCorr(data.Cp0Corr, data.Cp0Coef, data.MW, T);
    end Cp0Corr;

    function liqViscPcorLucas "Lucas liquid viscosity pressure correction."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T, p, Pref, rpVisc "temperature, pressure, reference pressure and viscosity at reference pressure";
      output Real visc;
    protected
      Real Tr, Tr2, Tr3, Tr4, A, C, D, dPr;
    algorithm
      Tr := T / data.Tc;
      Tr2 := Tr * Tr;
      Tr3 := Tr2 * Tr;
      Tr4 := Tr3 * Tr;
      dPr := (p - Pref) / data.criticalPressure;
      A := 0.9991 - 4.674e-4 / (1.0523 * Tr ^ (-0.03877) - 1.0513);
      C := (-0.07921) + 2.1616 * Tr - 13.404 * Tr2 + 44.1706 * Tr3 - 84.8291 * Tr4 + 96.1209 * Tr3 * Tr2 - 59.8127 * Tr3 * Tr3 + 15.6719 * Tr4 * Tr3;
      D := 0.3257 / (1.0039 - Tr ^ 2.573) ^ 0.2906 - 0.2086;
      visc := rpVisc * (1 + D * (dPr / 2.118) ^ A / (1 + C * data.w * dPr));
    end liqViscPcorLucas;

    function waterViscosity
      input Real T, rho;
      output Real visc;
    protected
      Real Tr, Tr2, Tr3, Ti, Ti2, Ti3, Ti4, Ti5, rhoR, ri, ri2, ri3, ri4, ri5, ri6, mu0, mu1, part[6];
    algorithm
      Tr := T / 647.096;
      Tr2 := Tr * Tr;
      Tr3 := Tr2 * Tr;
      rhoR := rho / 322.0;
      mu0 := 1.67752 + 2.20462 / Tr + 0.6366564 / Tr2 - 0.241605 / Tr3;
      mu0 := 100 * Tr ^ 0.5 / mu0;
      Ti := 1 / Tr - 1;
      Ti2 := Ti * Ti;
      Ti3 := Ti2 * Ti;
      Ti4 := Ti3 * Ti;
      Ti5 := Ti4 * Ti;
      ri := rhoR - 1;
      ri2 := ri * ri;
      ri3 := ri2 * ri;
      ri4 := ri3 * ri;
      ri5 := ri4 * ri;
      ri6 := ri5 * ri;
      part[1] := 5.20094e-1 + 2.22531e-1 * ri - 2.81378e-1 * ri2 + 1.61913e-1 * ri3 - 3.25372e-2 * ri4;
      part[2] := (8.50895e-2 + 9.99115e-1 * ri - 9.06851e-1 * ri2 + 2.57399e-1 * ri3) * Ti;
      part[3] := ((-1.08374) + 1.88797 * ri - 7.72479e-1 * ri2) * Ti2;
      part[4] := ((-2.89555e-1) + 1.26613 * ri - 4.89837e-1 * ri2 + 0 + 6.98452e-2 * ri4 - 4.35673e-3 * ri6) * Ti3;
      part[5] := ((-2.57040e-1 * ri2) + 8.72102e-3 * ri5) * Ti4;
      part[6] := (1.20573e-1 * ri - 5.93264e-4 * ri6) * Ti5;
      mu1 := exp(rhoR * (part[1] + part[2] + part[3] + part[4] + part[5] + part[6]));
      visc := mu0 * mu1 * 1e-6;
    end waterViscosity;

    function gasViscCorr "Calculates low pressure gas viscosity from a given DataRecord and a given temperature. Faster than gasViscCorr"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real eta;
    algorithm
      eta := if data.gViscCorr == 111 then data.gViscCoef[1] + T * (data.gViscCoef[2] + T * (data.gViscCoef[3] + T * (data.gViscCoef[4] + T * data.gViscCoef[4]))) else data.gViscCoef[1] * T ^ data.gViscCoef[2] / (1 + data.gViscCoef[3] / T + data.gViscCoef[4] / T ^ 2);
    end gasViscCorr;

    function gasViscLowPressureChung "Dynamic viscosity of a low pressure gas according to Chung"
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.Temperature T "Gas temperature";
      output SI.DynamicViscosity eta "Dynamic viscosity";
    protected
      Real Ta, sigma, omega, muR, muR4, Fc, k;
    algorithm
      if data.family == 6 then
        k := 0.076 "water";
      elseif data.family == 7 then
        k := 0.0682 + 4.74 / data.MW "alcohol";
      elseif data.family == 8 then
        k := 0.0682 + 2 * 4.74 / data.MW "polyol";
      else
        k := 0.0;
      end if;
      Ta := 1.2593 * T / data.Tc;
      sigma := 0.809 * data.Vc ^ 0.33333;
      omega := 1.16145 * Ta ^ (-0.14874) + 0.52487 * exp(-0.7732 * Ta) + 2.16178 * exp(-2.43787 * Ta);
      muR := if data.mu < 999.0 and data.mu > 0.0 then 131.3 * data.mu / (data.Vc * 1e6 * data.Tc) ^ 0.5 else 0.0;
      muR4 := muR * muR * muR * muR;
      Fc := 1 - 0.2756 * data.w + 0.059035 * muR4 + k;
      eta := 26.692 * Fc * (data.MW * T) ^ 0.5 * 1e-11 / (sigma * sigma * omega);
    end gasViscLowPressureChung;
    
    function gasViscTVcpChung "Dynamic viscosity prediction with density correction. Chung method. Fails because the molar volume V is not supplied "
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.Temperature T "Gas temperature";
      input SI.DynamicViscosity eta0In=0 "low density viscosity, if known";
      output SI.DynamicViscosity eta "Dynamic viscosity";
    protected
      Real Ta, sigma, omega, muR, muR4, Fc, k, eta0;
      Real E[10],y,G1,G2,eta2,eta1;
      Real coef[10,4]={{6.324,50.412,-51.680,1189.0},{1.210e-3,-1.154e-3,-6.257e-3,0.03728},{5.283,254.209,-168.48,3898},{6.623,38.096,-8.464,31.42},
                            {19.745,7.63,-14.354,31.53},{-1.9,-12.537,4.985,-18.15},{24.275,3.45,-11.291,69.35},{0.7972,1.117,0.01235,-4.117},
                            {-0.2382,0.0677,-0.8163,4.025},{0.06863,0.3479,0.5926,-0.727}};
    algorithm
      Ta := 1.2593 * T / data.Tc;
      if data.family == 6 then
          k := 0.076 "water";
        elseif data.family == 7 then
          k := 0.0682 + 4.74 / data.MW "alcohol";
        elseif data.family == 8 then
          k := 0.0682 + 2 * 4.74 / data.MW "polyol";
        else
          k := 0.0;
      end if;
      muR := if data.mu < 999.0 and data.mu > 0.0 then 131.3 * data.mu / (data.Vc * 1e6 * data.Tc) ^ 0.5 else 0.0;
        muR4 := muR * muR * muR * muR;
      if eta0In==0 then 
        sigma := 0.809 * data.Vc ^ 0.33333;
        omega := 1.16145 * Ta ^ (-0.14874) + 0.52487 * exp(-0.7732 * Ta) + 2.16178 * exp(-2.43787 * Ta);
        Fc := 1 - 0.2756 * data.w + 0.059035 * muR4 + k;
        eta0 := 26.692 * Fc * (data.MW * T) ^ 0.5 * 1e-11 / (sigma * sigma * omega);
      else
        eta0:=eta0In;
      end if;
      for i in 1:10 loop
        E[i]:=coef[i,1]+coef[i,2]*data.w+coef[i,3]*muR4+coef[i,4]*k;
      end for;
      y:=data.Vc/(6* V);
      G1:=(1-0.5*y)/(1-y)^3;
      G2:=(E[0]*((1-exp(-E[3]*y))/y)+E[1]*G1*exp(E[4]*y)+E[2]*G1)/(E[0]*E[3]+E[1]+E[2]);
      eta2:=E[6]*y*y*G2*exp(E[7]+E[8]/Ta+E[9]/(Ta*Ta));
      eta1:=pTa^0.5*(Fc*(1/G2+E[5]*y))/omega+eta2;
      eta:= if(eta1>1) then eta0*eta1 else eta0;
    end gasViscTVcpChung;

    function gasViscPcorLucas
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T, p, lpVisc;
      output Real visc;
    protected
      Real Z1, mur, Q, xi, Tr, Fp0, Fq0, Z2, Pr, a, b, c, d, e, f, Y, Fp, Fq;
    algorithm
      Tr := T / data.Tc;
      Pr := p / data.criticalPressure;
      if data.mu < 999 and data.mu > 0 then
        mur := 52.46 * data.mu * data.mu * data.criticalPressure * 1e-5 / (data.Tc * data.Tc);
      else
        mur := 0;
      end if;
      if mur < 0.022 then
        Fp0 := 1;
      elseif mur < 0.075 then
        Fp0 := 1 + 30.55 * (0.292 - data.Zc) ^ 1.72;
      else
        Fp0 := 1 + 30.55 * (0.292 - data.Zc) ^ 1.72 * abs(0.96 + 0.1 * (Tr - 0.7));
      end if;
      if data.MW < 4.04 then
        Q := if data.MW < 3 then 0.76 elseif data.MW < 4.02 then 1.38 else 0;
        Fq0 := 1.22 * Q ^ 0.15;
      else
        Fq0 := 1;
      end if;
      xi := 0.176 * (data.Tc / (data.MW * data.MW * data.MW * data.criticalPressure * data.criticalPressure * data.criticalPressure * data.criticalPressure * 1e-20)) ^ 0.1666667 * 1e-7;
      Z1 := lpVisc * xi;
      if Tr <= 1 then
        Tr := 1;
      end if;
      a := 1.245e-3 * exp(5.1726 * Tr ^ (-0.3286)) / Tr;
      b := a * (1.6553 * Tr - 1.2723);
      c := 0.4489 * exp(3.0578 * Tr ^ (-37.7332)) / Tr;
      d := 1.7368 * exp(2.231 * Tr ^ (-7.6351)) / Tr;
      e := 1.3088;
      f := 0.9425 * exp(-0.1853 * Tr ^ 0.4489);
      Z2 := Z1 * (1 + a * Pr ^ e / (b * Pr ^ f + 1 / (1 + c * Pr ^ d)));
      Y := Z2 / Z1;
      Fp := (1 + (Fp0 - 1) * Y ^ (-3)) / Fp0;
      Fq := (1 + (Fq0 - 1) * (1 / Y - 0.007 * log(Y) ^ 4)) / Fq0;
      visc := Z2 * Fp * Fq / xi;
    end gasViscPcorLucas;

    function gasThCondCorr "Calculates low pressure gas thermal conductivity from a given correlation data and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gThCondCoef directly, coming from an iteration."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real lambda;
    algorithm
      lambda := if data.gThCondCorr == 121 then data.gThCondCoef[1] + T * (data.gThCondCoef[2] + T * (data.gThCondCoef[3] + T * (data.gThCondCoef[4] + T * data.gThCondCoef[4]))) else data.gThCondCoef[1] * T ^ data.gThCondCoef[2] / (1 + data.gThCondCoef[3] / T + data.gThCondCoef[4] / T ^ 2);
    end gasThCondCorr;

    function gasThCondLowPressureChung "Thermal conductivity of a low pressure gas according to Chung"
      input FreeFluids.MediaCommon.DataRecord data;
      input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp0;
      input SI.DynamicViscosity eta "low pressure dynamic viscosity";
      input SI.Temperature T "Gas temperature";
      output Modelica.Media.Interfaces.Types.ThermalConductivity lambda;
    protected
      Real R, alpha, beta, Tr, Z, psi;
    algorithm
      R := Modelica.Constants.R * 1000 / data.MW;
      alpha := (Cp0 - R) / R - 1.5;
      beta := 0.7862 - 0.7109 * data.w + 1.3168 * data.w * data.w;
      Tr := T / data.Tc;
      Z := 2.0 + 10.5 * Tr * Tr;
      psi := 1 + alpha * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * Z) / (0.6366 + beta * Z + 1.061 * alpha * beta));
      lambda := 3.75 * psi * eta * R;
    end gasThCondLowPressureChung;

    function gasThCondTVcorChung
    "Density correction for low pressure gas thermal conductivity, according to Chung"
      input SI.Temperature T;
      input SI.Density d;
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.ThermalConductivity ldThCond;
      output Modelica.Media.Interfaces.Types.ThermalConductivity thCond;
    protected
      Real Tr,y,B[7],G1,muR,muR4,k,G2;
      Real coef[7,4]={{2.4166,7.4824e-1,-9.1858e-1,1.2172e2},{-5.0924e-1,-1.5094,-4.9991e1,6.9983e1},{6.6107,5.6207,6.4760e1,2.7039e1},{1.4543e1,-8.9139,-5.6379,7.4344e1},
                            {7.9274e-1,8.2019e-1,-6.9369e-1,6.3173},{-5.8634,1.2801e1,9.5893,6.5529e1},{9.1089e1,1.2811e2,-5.4217e1,5.2381e2}};
    algorithm
      Tr:=T/data.Tc;
      y:=data.Vc*d/(6*data.molarMass);
      G1:=(1-0.5*y)/(1-y)^3;
      muR:= if data.mu>0 then 131.3*data.mu/(data.Vc*1e6*data.Tc)^0.5 else 0;
      muR4:=muR*muR*muR*muR;
      if data.family == 6 then
        k := 0.076 "water";
      elseif data.family == 7 then
        k := 0.0682 + 4.74 / data.MW "alcohol";
      elseif data.family == 8 then
        k := 0.0682 + 2 * 4.74 / data.MW "polyol";
      else
        k := 0.0;
      end if;
      for i in 1:7 loop
        B[i]:=coef[i,1]+coef[i,2]*data.w+coef[i,3]*muR4+coef[i,4]*k;
      end for;
      G2:=(B[1]*(1-exp(-B[4]*y))/y+B[2]*G1*exp(B[5]*y)+B[3]*G1)/(B[1]*B[4]+B[2]+B[3]);
      thCond:=ldThCond*(1/G2+B[6]*y)+3.586e-3*(data.Tc*1000/data.MW)^0.5*B[7]*y*y*Tr^0.5*G2/(data.Vc*1e6)^0.6667;
      if (thCond<ldThCond) then
        thCond:=ldThCond;
      end if;
    end gasThCondTVcorChung;
  
    function gasMixViscosityWilke
      input Modelica.Units.SI.Temperature T;
      input Modelica.Units.SI.MoleFraction X[:];
      input Modelica.Units.SI.DynamicViscosity etaX[:];
      input Modelica.Units.SI.MolarMass MW[:];
      output Modelica.Units.SI.DynamicViscosity eta;
    protected
      Real phi[size(MW, 1), size(MW, 1)];
      Real sigma;
    algorithm
      for i in 1:size(X, 1) loop
        for j in 1:size(X, 1) loop
          phi[i, j] := (etaX[i] / etaX[j]) ^ 0.5 * (MW[i] / MW[j]) ^ 0.25;
          phi[i, j] := (1 + phi[i, j]) ^ 2.0 / (8 * (1 + MW[i] / MW[j])) ^ 0.5;
        end for;
      end for;
      eta := 0;
      for i in 1:size(X, 1) loop
        sigma := 0;
        for j in 1:size(X, 1) loop
          sigma := sigma + X[j] * phi[i, j];
        end for;
        eta := eta + X[i] * etaX[i] / sigma;
      end for;
    end gasMixViscosityWilke;

    function gasMixThCondMason
      input SI.Temperature T;
      input SI.MoleFraction X[:];
      input SI.ThermalConductivity lambdaX[size(X, 1)];
      input SI.MolarMass MW[size(X, 1)];
      input SI.Temperature Tc[size(X, 1)];
      input SI.AbsolutePressure Pc[size(X, 1)];
      output SI.ThermalConductivity lambda;
    protected
      Real Tr[size(X, 1)];
      Real Gamma[size(X, 1)];
      Real A[size(X, 1), size(X, 1)];
      Real ratio;
      Real sigma;
    algorithm
      for i in 1:size(X, 1) loop
        Tr[i] := T / Tc[i];
        Gamma[i] := 210 * (Tc[i] * MW[i] * MW[i] * MW[i] / (Pc[i] * Pc[i] * Pc[i] * Pc[i] * 1e-20)) ^ 0.1666667;
      end for;
      lambda := 0;
      for i in 1:size(X, 1) loop
        sigma := 0;
        for j in 1:size(X, 1) loop
          ratio := Gamma[j] * (exp(0.0464 * Tr[i]) - exp(-0.2412 * Tr[i])) / (Gamma[i] * (exp(0.0464 * Tr[j]) - exp(-0.2412 * Tr[j])));
          A[i, j] := (1 + ratio ^ 0.5 * (MW[i] / MW[j]) ^ 0.25) ^ 2 / (8 * (1 + MW[i] / MW[j])) ^ 0.5;
          sigma := sigma + X[j] * A[i, j];
        end for;
        lambda := lambda + X[i] * lambdaX[i] / sigma;
      end for;
    end gasMixThCondMason;

    function liqThCondLatini
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.Temperature T;
      output SI.ThermalConductivity lambda;
    protected
      Real Tr = T / data.Tc, A, alpha, beta, gamma;
    algorithm
      if data.family == 1 then
        A := 0.00350 "Alkane";
        alpha := 1.2;
        beta := 0.5;
        gamma := 0.167;
      elseif data.family == 2 or data.family == 3 then
        A := 0.0361 "Alkene, Alkyne";
        alpha := 1.2;
        beta := 1.0;
        gamma := 0.167;
      elseif data.family == 4 then
        A := 0.031 "Cycloalkane";
        alpha := 1.2;
        beta := 1.0;
        gamma := 0.167;
      elseif data.family == 5 or data.family == 6 then
        A := 0.0346 "Aromatic, Water";
        alpha := 1.2;
        beta := 1.0;
        gamma := 0.167;
      elseif data.family == 7 or data.family == 8 then
        A := 0.00339 "Alcohol, Polyol";
        alpha := 1.2;
        beta := 0.5;
        gamma := 0.167;
      elseif data.family == 10 then
        A := 0.0385 "Ether";
        alpha := 1.2;
        beta := 1.0;
        gamma := 0.167;
      elseif data.family == 12 then
        A := 0.00383 "Ketone";
        alpha := 1.2;
        beta := 0.5;
        gamma := 0.167;
      elseif data.family == 13 then
        A := 0.00319 "Acid";
        alpha := 1.2;
        beta := 0.5;
        gamma := 0.167;
      elseif data.family == 14 then
        A := 0.0415 "Ester";
        alpha := 1.2;
        beta := 1.0;
        gamma := 0.167;
      elseif data.family == 17 or data.family == 18 then
        A := 0.562 "haloalkanes, haloalkenes";
        alpha := 0.0;
        beta := 0.5;
        gamma := -0.167;
      else
        A := 0.494;
        alpha := 0.0;
        beta := 0.5;
        gamma := -0.167;
      end if;
      lambda := A * data.Tb ^ alpha * (1 - Tr) ^ 0.38 / (data.MW ^ beta * data.Tc ^ gamma * Tr ^ (1 / 6));
    end liqThCondLatini;

    function liqSurfTensSastriRao
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.Temperature T;
      output SI.SurfaceTension sigma;
    protected
      Real k, x, y, z, m;
    algorithm
      if data.family == 7 or data.family == 8 then
        k := 2.28 "alcohol, polyol";
        x := 0.25;
        y := 0.175;
        z := 0;
        m := 0.8;
      elseif data.family == 13 then
        k := 0.125 "acid";
        x := 0.5;
        y := -1.5;
        z := 1.85;
        m := 1.222;
      else
        k := 0.158;
        x := 0.5;
        y := -1.5;
        z := 1.85;
        m := 1.222;
      end if;
      sigma := k * (data.criticalPressure * 1e-5) ^ x * data.Tb ^ y * data.Tc ^ z * ((1 - T / data.Tc) / (1 - data.Tb / data.Tc)) ^ m * 1e-3;
    end liqSurfTensSastriRao;
    
  end Functions;

  annotation(
    Documentation(info = "<html><head></head><body>
    <p>Provides the general data and functions that are common to several media packages. </p>
    <p><b>DataRecord and associated packages</b></p>
    <p></p>
    <p>The DataRecord record in the MediaCommon package defines the container for the individual substance data. This data will be used later by the media packages written directly in Modelica language: LMedia and TMedia. It contains basic data regarding the substance: name, description, CAS, chemical family,...and correlations for several physical properties, normally as function of temperature. Each correlation has: equation to use for its calculation, coefficients, and limits of usage.</p>
    <p>The data for each individual fluid is inside the subpackages: MediaDataAL, and MediaDataMZ. Each data is defined as a constant DataRecord, with the name of the substance.</p>
  <p>You can create a new record, in those subpackages, copying the MediaDataTemplate (it is inside the MediaCommon package) and filling it manually. Nevertheless the faster and more convenient way is to create the record from the FreeFluidsGui program. You need to select the substance from the database, select the correlations you want to be included in the record, and export it in Modelica format. You can put the file in any place, better with the .txt extension. Later you edit the file, copy its content, and paste it inside the MediaData package. You still need to declare the substance inside the Media packages, filling the name and the origin of the data, that will be the record you just copied.</p>
  <p>When using the data in the TMedia package, you will need at least the following correlations: saturated density, vapor pressure, liquid heat capacity. Plus vaporization enthalpy and saturated gas density, if you want to work also with the gas phase. Transport properties correlations as for your needs, and reduced liquid bulk modulus if you want to work at high pressures (till 200 bars). If you want to force the liquid state you can set Tc at a high value, and Pc at a low value (look as example to the MarlothermSH medium), but is better to use LMedia package.</p>
    <p></p>
    <p></p>
    
    </body></html>"));
end MediaCommon;
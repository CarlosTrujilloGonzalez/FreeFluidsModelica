within FreeFluids;

package MediaCommon "MediaCommon.mo by Carlos Trujillo
      This file is part of the Free Fluids application
      Copyright (C) 2008-2019  Carlos Trujillo Gonzalez
        
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
    Real MW = 0.0, molarMass = 0.0, Tc = 0.0, criticalPressure=0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 6.667e-10, lnuA = 0.0, lnuB = 0.0;
    Integer Cp0Corr = 0, VpCorr = 0, BtCorr = 0, HvCorr = 0, lDensCorr = 0, lCpCorr = 0, lTfromHsatCorr = 0, lViscCorr = 0, lThCondCorr = 0, lSurfTensCorr = 0, lBulkModRCorr = 0, gSatDensCorr = 0, gViscCorr = 0, gThCondCorr = 0;
    Real Cp0Coef[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, VpCoef[6] = {0, 0, 0, 0, 0, 0}, BtCoef[6] = {0, 0, 0, 0, 0, 0}, HvCoef[6] = {0, 0, 0, 0, 0, 0}, lDensCoef[6] = {0, 0, 0, 0, 0, 0}, lCpCoef[6] = {0, 0, 0, 0, 0, 0}, lTfromHsatCoef[6] = {0, 0, 0, 0, 0, 0}, lViscCoef[6] = {0, 0, 0, 0, 0, 0}, lThCondCoef[6] = {0, 0, 0, 0, 0, 0}, lSurfTensCoef[6] = {0, 0, 0, 0, 0, 0}, lBulkModRCoef[6] = {0, 0, 0, 0, 0, 0}, gSatDensCoef[6] = {0, 0, 0, 0, 0, 0}, gViscCoef[6] = {0, 0, 0, 0, 0, 0}, gThCondCoef[6] = {0, 0, 0, 0, 0, 0};
    Real Cp0LimI = 0, VpLimI = 0, BtLimI = 0, HvLimI = 0, lDensLimI = 0, lCpLimI = 0, lTfromHsatLimI = 0, lViscLimI = 0, lThCondLimI = 0, lSurfTensLimI = 0, lBulkModRLimI = 0, gSatDensLimI = 0, gViscLimI = 0, gThCondLimI = 0;
    Real Cp0LimS = 0, VpLimS = 0, BtLimS = 0, HvLimS = 0, lDensLimS = 0, lCpLimS = 0, lTfromHsatLimS = 0, lViscLimS = 0, lThCondLimS = 0, lSurfTensLimS = 0, lBulkModRLimS = 0, gSatDensLimS = 0, gViscLimS = 0, gThCondLimS = 0;
  end DataRecord;

  record HelmholtzDerivatives
    Real a, av, avv, at, att, avt;
    Real a, av, avv, at, att, avt;
    Real a, av, avv, at, att, avt;
  end HelmholtzDerivatives;

  record IdealThermo
    Real cp, h, s;
    Real cp, h, s;
    Real cp, h, s;
  end IdealThermo;

  constant DataRecord MediaDataTemplate(name = "", description = "", CAS = "", family = 0, MW = 0.0, molarMass = 0.0, Tc = 0.0, Pc = 0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 0.0, lnuA = 0.0, lnuB = 0.0, Cp0Corr = 0, Cp0Coef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 0.0, Cp0LimS = 0.0, VpCorr = 0, VpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, VpLimI = 0.0, VpLimS = 0.0, BtCorr = 0, BtCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, BtLimI = 0.0, BtLimS = 0.0, HvCorr = 0, HvCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, HvLimI = 0.0, HvLimS = 0.0, lDensCorr = 0, lDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lDensLimI = 0.0, lDensLimS = 0.0, lCpCorr = 0, lCpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lCpLimI = 0.0, lCpLimS = 0.0, lTfromHsatCorr = 0, lTfromHsatCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lTfromHsatLimI = 0.0, lTfromHsatLimS = 0.0, lViscCorr = 0, lViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lViscLimI = 0.0, lViscLimS = 0.0, lThCondCorr = 0, lThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lThCondLimI = 0.0, lThCondLimS = 0.0, lSurfTensCorr = 0, lSurfTensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lSurfTensLimI = 0.0, lSurfTensLimS = 0.0, lBulkModRCorr = 0, lBulkModRCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lBulkModRLimI = 0.0, lBulkModRLimS = 0.0, gSatDensCorr = 0, gSatDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gSatDensLimI = 0.0, gSatDensLimS = 0.0, gViscCorr = 0, gViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gViscLimI = 0.0, gViscLimS = 0.0, gThCondCorr = 0, gThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gThCondLimI = 0.0, gThCondLimS = 0.0);

  package Types
    type InputChoice = enumeration(pT "(p,T) as inputs", ph "(p,h) as inputs", ps "(p,s) as inputs", dT "(d,T) as inputs");
    type ReferenceState = enumeration(ASHRAE "0 enthalpy and entropy for saturated liquid at -40ºC", NBP "0 enthalpy and entropy for saturated liquid at normal boiling point", IIR "H=200 KJ/kg, and S=1.0 KJ/(kg·K) for saturated liquid at 0ºC", User "User defined reference_T and reference_P");
    type CorrelationEquation = enumeration(FF_DIPPR100, FF_Polynomial, FF_Polynomial2, FF_DIPPR100Ld, FF_expDIPPR100, FF_DIPPR101, FF_DIPPR101Vp, FF_DIPPR101Lv, FF_logDIPPR101, FF_DIPPR102, FF_DIPPR103, FF_DIPPR104, FF_DIPPR105, FF_DIPPR106, FF_DIPPR106Hv, FF_DIPPR106Ld, FF_DIPPR106SurfT, FF_DIPPR107, FF_DIPPR107Cp, FF_DIPPR114, FF_DIPPR115, FF_DIPPR116, FF_DIPPR116Ld, FF_Wilhoit, FF_Cooper, FF_Jaechske, FF_ChemSep16, FF_Antoine1, FF_Antoine2, FF_Wagner25, FF_Wagner36, FF_PPDS9, FF_PPDS10, FF_PCWIN, FF_Rackett, FF_ExtAndrade1, FF_ExtAndrade2, FF_ChericVisc, FF_WagnerGd, FF_Tait, FF_ExtWagner, FF_PPDS15);
  end Types;

  package Functions
    function PhysPropCorrCalc
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real y;
    protected
      Real Tm;
      Integer j;
    algorithm
    /*
      if cor == 22 or cor == 63 then
        x := x - 273.15 "First we convert the input variable if necessary";
      else
        x := x;
      end if;*/
      if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 or cor == 40 or cor == 49 or cor == 50 or cor == 60 or cor == 63 or cor == 70 or cor == 71 or cor == 80 or cor == 82 or cor == 100 or cor == 111 or cor == 121 or cor == 140 or cor == 150 then
        y := coef[1] + x * (coef[2] + x * (coef[3] + x * (coef[4] + x * coef[5]))) "FF_DIPPR100.ok";
      elseif cor == 10 then
        y := coef[1] + x * (coef[2] + x * (coef[3] + x * (coef[4] + x * (coef[5] + x * coef[6])))) "FF_Polynomial";
      elseif cor == 130 then
        y := coef[1] + coef[2] * x ^ 0.125 + coef[3] * x ^ 0.25 + coef[4] * x ^ 0.5 + coef[5] * x "FF_Polynomial2.ok";
      elseif cor == 20 or cor == 21 or cor == 30 then
        y := exp(coef[1] + coef[2] / x + coef[3] * log(x) + coef[4] * x ^ coef[5]) "FF_DIPPR101.ok";
      elseif cor == 81 or cor == 110 or cor == 120 then
        y := coef[1] * x^ coef[2] / (1 + coef[3] / x + coef[4] / x ^ 2) "FF_DIPPR102.ok";
      elseif cor == 41 or cor == 42 then
        y := coef[1] / coef[2] ^ (1 + (1 - x / coef[3]) ^ coef[4]) "FF_DIPPR105";
      elseif cor == 45 or cor == 47 or cor == 61 or cor == 90 or cor == 91 then
        Tm := x / coef[6];
        y := coef[1] * (1 - Tm) ^ (coef[2] + coef[3] * Tm + coef[4] * Tm ^ 2 + coef[5] * Tm ^ 3) "FF_DIPPR106.ok";
      elseif cor == 3 or cor == 4 or cor == 200 then
        y := coef[1] + coef[2] * (coef[3] / x / sinh(coef[3] / x)) ^ 2 + coef[4] * (coef[5] / x / cosh(coef[5] / x)) ^ 2 "FF_DIPPR107.ok";
      elseif cor == 5 then
        Tm := (x - coef[8]) / (x + coef[7]);
        y := R * (coef[1] + coef[2] / x^ 2 * exp(-coef[3] / x) + coef[4] * Tm ^ 2 + (coef[5] - coef[6] / (x - coef[8]) ^ 2) * Tm^ 8) "FF_Wilhoit.ok";
      elseif cor == 7 then
        y := coef[1] + coef[2] * x^coef[3] "FF_Cooper.ok";
        for j in {4, 6, 8, 10} loop
          if coef[j] > 0 then
            y := y + coef[j] * (coef[j + 1] / x)^ 2 * exp(coef[j + 1] / x) / (exp(coef[j + 1] / x) - 1) ^ 2;
          end if;
        end for;
        y := y * R;
      elseif cor == 8 then
        y := 1 + coef[1] "FF_Jaechske.ok";
        for j in {2, 6} loop
          if coef[j] > 0 then
            y := y + coef[j] * (coef[j + 1] / x / sinh(coef[j + 1] / x))^2 + coef[j + 2] * (coef[j + 3] / x / cosh(coef[j + 3] / x))^2;
          end if;
        end for;
        y := y * R;
      elseif cor == 9 or cor == 18 or cor == 51 or cor == 62 then
        y := coef[1] + exp(coef[2] / x + coef[3] + coef[4] * x + coef[5] * x ^ 2) "FF_ChemSep16";
      elseif cor == 22 or cor == 36 then
        y:=10^(coef[1]-coef[2]/(x+coef[3]))"FF_Antoine1";
      elseif cor == 23 then
        y:=exp(coef[1]-coef[2]/(x+coef[3]))"FF_Antoine2";
      elseif cor == 24 then
        Tm:=1-x/coef[6];
        y:=coef[1]*exp((coef[2]*Tm+coef[3]*Tm^1.5+coef[4]*Tm^2.5+coef[5]*Tm^5.0)/(1-Tm))"FF_Wagner25";
      elseif cor == 25 then
        Tm:=1-x/coef[6];
        y:=coef[1]*exp((coef[2]*Tm+coef[3]*Tm^1.5+coef[4]*Tm^3+coef[5]*Tm^6)/(1-Tm))"FF_Wagner36.ok";
      elseif cor == 26 then
        Tm:=1-x/coef[2];
        y:=coef[1]*exp(coef[2]*(coef[3]*Tm^coef[4]+coef[5]*Tm^coef[6]+coef[7]*Tm^coef[8]+coef[9]*Tm^coef[10]+
           coef[11]*Tm^coef[12]+coef[13]*Tm^coef[14])/x)"FF_ExtWagner";
      elseif cor == 31 or cor == 33 then
        y:=exp(coef[1]+coef[2]/x+coef[3]*x+coef[4]*x^2+coef[5]*x^3)"FF_ExtAndrade1";
      elseif cor == 32 then
        y:=10^(coef[1]+coef[2]/x+coef[3]*x+coef[4]*x^2+coef[5]*x^3)"FF_ExtAndrade2";
      elseif cor == 37 then
        y:=coef[5]*exp(coef[1]*((coef[3]-x)/(x-coef[4]))^0.33333+coef[2]*((coef[3]-x)/(x-coef[4]))^1.33333)"FF_PPDS9";
      elseif cor == 19 or cor == 151 then
        Tm:=1-x/coef[6];
        y:=coef[1]/Tm+coef[2]+coef[3]*Tm+coef[4]*Tm*Tm+coef[5]*Tm*Tm*Tm "FF_PPDS15.ok";
      elseif cor == 43 then
        y:=coef[1]+coef[2]*(1-x/coef[5])+coef[3]*log(1-x/coef[5])+coef[4]*(1-x/coef[5])^3 "FF_PCWIN";
      elseif cor == 48 then
        y:=coef[1]/coef[2]^((1-x/coef[3])^coef[4]) "FF_Rackett";
      elseif cor == 44 or cor == 46 then
        Tm:=1-x/coef[6];
        y:=coef[1]+coef[2]*Tm^0.35+coef[3]*Tm^0.666667+coef[4]*Tm+coef[5]*Tm^1.333333 "FF_DIPPR116";
      elseif cor == 101 then
        Tm:=1-x/coef[6];
        y:=coef[1]*exp(coef[2]*Tm^0.4+coef[3]*Tm+coef[4]*Tm^2.1+coef[5]*Tm^5.6) "FF_WagnerGd.ok";
      elseif cor == 240 then
        Tm:=x-273.15;
        y:=(coef[1]+coef[2]*Tm+coef[3]*Tm^2)*(1-0.0894*log(1+1e5/(coef[4]*exp(-coef[5]*Tm)))) "FF_Tait at 1 bar. T in centigrades";
      end if;

      if cor == 1 or cor == 21 or cor == 48 then
        y := y * 1e3;
      elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
        y:=y/MW * 1e3;
      elseif cor == 3 or cor == 10 then
        y:=y/MW *1e3*4.1868;
      elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 or cor == 80 or cor == 81 or cor == 90 then
        y:=y/ MW;
      elseif cor == 41 or cor == 44 or cor == 45 or cor == 70 then
        y:=y* MW;
      elseif cor == 22 then
        y:=y*133.32239;
      elseif cor == 31 or cor == 32 or cor == 34 or cor == 36 or cor == 63 then
        y:=y*1e-3;
      elseif cor == 33 then
        y:=y* MW*1e-3;
      elseif cor == 49 then
        y:=MW*1e3/ y;
      elseif cor == 240 then
        y:=1/ y;
      end if;
/*         
        switch (cor){
        case 1://DIPPR 100 Cp0 in KJ/kgr·K
        case 2://DIPPR 100 Cp0 in J/mol·K
        case 6://DIPPR 100 Cp0 in J/kgr·K
        case 15://DIPPR 100 Liquid Cp in J/mol·K
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 17://DIPPR 100 Liquid Cp in J/kgr·K
        case 40://DIPPR 100 Liquid density in kgr/m3
        case 49://DIPPR 100 Liquid density in cm3/mol
        case 50://DIPPR 100 Liquid thermal conductivity W/(m·K)
        case 60://DIPPR 100 Liquid surface tension N/m
        case 63://DIPPR 100 Liquid surface tension dyna/cm
        case 70://DIPPR 100 Solid density in Kmol/m3
        case 71://DIPPR 100 Solid density in kgr/m3
        case 80://DIPPR 100 Solid Cp in J/Kmol·K
        case 82://DIPPR 100 Solid Cp in J/(kgr·K)
        case 100://DIPPR 100 Vapor density in kgr/m3
        case 111://DIPPR100 Gas viscosity in Pa·s
        case 121://DIPPR100 Gas thermal conductivity in W/(m·K)
        case 140://DIPPR100 T in K from liquid enthalpy in J/kgr
    	case 150://DIPPR100 Isothermal compressibility adimensional
            eq=FF_DIPPR100;
            break;
        case 10://Polynomial Cp0 en cal/(mol·K)
            eq=FF_Polynomial;
            //printf("Equation:%i\n",eq);
            break;
        case 130://Polynomial2 boil temp in K
            eq=FF_Polynomial2;
            break;
        case 20://DIPPR101 Vp Pa
        case 21://DIPPR101 Vp KPa
        case 30://DIPPR101 Liquid viscosity Pa·s
            eq=FF_DIPPR101;
            break;
        case 81://DIPPR 102 Solid Cp in J/Kmol·K
        case 110://DIPPR 102 Gas viscosity in Pa·s
        case 120://DIPPR 102 Gas thermal conductivity in W/(m·K)
            eq=FF_DIPPR102;
            break;
        case 41://DIPPR105 Liquid density in mol/dm3
        case 42://DIPPR105 Liquid density in kgr/m3
            eq=FF_DIPPR105;
            break;
        case 45://DIPPR106 Liquid density in mol/dm3
        case 47://DIPPR106 Liquid density in kg/m3
        case 61://DIPPR 106 Liquid surface tension N/m
        case 90://DIPPR 106 HvSat J/Kmol
        case 91://DIPPR 106 HvSat J/kgr
            eq=FF_DIPPR106;
            break;
        case 3://DIPPR107 Cp in cal/mol·K
        case 4://DIPPR107 Cp in J/Kmol·K
    	case 200://DIPPR107 Cp in J/kg·K
            eq=FF_DIPPR107;
            break;
        case 5://Wilhoit Cp0 J/mol·K
            eq=FF_Wilhoit;
            break;
        case 7://Cooper Cp0 J/mol·K
            eq=FF_Cooper;
            break;
        case 8://Jaechske Cp0 J/mol·K
            eq=FF_Jaechske;
            break;
        case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        case 51://ChemSep nº16 Liquid thermal conductivity W/(m·K)
        case 62://chemSep nº16 Liquid surface tension N/m
            eq=FF_ChemSep16;
            break;
        case 22: //Antoine base 10 in C and mmHg
        case 36: //Antoine base 10 in K and cP
            eq=FF_Antoine1;
            break;
        case 23://Antoine base e in K and Pa
            eq=FF_Antoine2;
            break;
        case 24://Wagner 25 in K and Pa
            eq=FF_Wagner25;
            break;
        case 25://Wagner 36 in K and Pa
            eq=FF_Wagner36;
            break;
        case 26://Wagner 36 in K and Pa
            eq=FF_ExtWagner;
            break;
        case 31://Extended Andrade 1 in cP
        case 33://Extended Andrade 1 in cP/MW
            eq=FF_ExtAndrade1;
            break;
        case 32://Extended Andrade 2 in cP
            eq=FF_ExtAndrade2;
            break;
        case 34://Cheric viscosity in cP
            eq=FF_ChericVisc;
            break;
        case 37://PPDS9 in Pa·s
            eq=FF_PPDS9;
            break;
        case 19://PPDS15 Cp liquid in J/kg·K
        case 151://PPDS15 liquid isothermal compressibility adimensional
            eq=FF_PPDS15;
            break;
        case 43://PCWIN liquid density in kgr/m3
            eq=FF_PCWIN;
            break;
        case 48://Rackett liquid density in kgr/dm3
            eq=FF_Rackett;
            break;
        case 44://DIPPR116 liquid density in mol/dm3
        case 46://DIPPR116 liquid density in kg/m3
            eq=FF_DIPPR116;
            break;
        case 101:
            eq=FF_WagnerGd;
            break;
        case 240://Tait equation for polymer density in m3/kg at 1 bar
            eq=FF_Tait;
            break;
        }
    
        FF_CorrelationResult(&eq,coef,x,y);//We make the calculations
    
        switch (cor){//Now we convert to SI units
        case 1://DIPPR 100 Cp in KJ/kgr·K
        case 21://DIPPR101 Vp KPa
        case 48://Rackett liquid density in kgr/dm3
            *y=*y*1e3;
            break;
        case 2://DIPPR 100 Cp in J/mol·K
        case 5://Wilhoit Cp0 J/mol·K
        case 7://Cooper Cp0 J/mol·K
        case 8://Jaechske Cp0 J/mol·K
        case 15://DIPPR 100 Liquid Cp in J/mol·K
            *y=*y/ MW*1e3;
            break;
        case 3://DIPPR107 Cp in cal/mol·K
        case 10://Polynomial Cp0 en cal/(mol·K)
            *y=*y/ MW*1e3*4.1868;
            break;
        case 4://DIPPR107 Cp in J/Kmol·K
        case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        case 80://DIPPR 100 Solid Cp in J/Kmol·K
        case 81://DIPPR 102 Solid Cp in J/Kmol·K
        case 90://DIPPR 106 HvSat J/Kmol
            *y=*y/ MW;
            break;
        case 41://DIPPR105 Liquid density in mol/dm3
        case 44://DIPPR116 liquid density in mol/dm3
        case 45://DIPPR106 liquid density in mol/dm3
        case 70://DIPPR 100 Solid density in Kmol/m3
            *y=*y* MW;
            break;
        case 22://Antoine base 10 in C and mmHg
            *y=*y*133.32239;
            break;
        case 31://Extended Andrade 1 in cP
        case 32://Extended Andrade 2 in cP
        case 34://Cheric viscosity in cP
        case 36://Antoine1 viscosity in cP
        case 63://DIPPR 100 Liquid surface tension dyna/cm
            *y=*y*1e-3;
            break;
        case 33://Extended Andrade 1 in cP/MW
            *y=*y* MW*1e-3;
            break;
        case 49://DIPPR 100 Liquid density in cm3/mol
            *y=MW*1e3/ *y;
            break;
        case 240://Tait equation for polymer density in m3/kg at 1 bar
            *y=1/ *y;
        }
        //Last we reconvert the input variable if necessary
        if ((cor==22)||(cor==63)||(cor==240)) x=x+273.15;*/
    end PhysPropCorrCalc;

    function SpecificEnthalpyCorrCalc
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real H "Enthalpy J/kg";
    protected
      Real y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
      Real Tm,Tm2;
      Integer j;
    algorithm
        if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 then
          H:=x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*(coef[4]/4+x*coef[5]/5)))) "DIPPR 100";
          
        elseif cor == 10 then 
          H:=x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*(coef[4]/4+x*coef[5]/5)))) "Polynomial";
          
        elseif cor == 3 or cor == 4 or cor == 200 then
          H:=coef[1]*x+coef[2]*coef[3]*(1/tanh(coef[3]/x))-coef[4]*coef[5]*tanh(coef[5]/x) "DIPPR107.ok";
        
        elseif cor == 5 then
          a7_a6:=coef[8]/coef[7];
          a7_a6_2:=a7_a6*a7_a6;
          a7_a6_4:=a7_a6_2*a7_a6_2;
          x1:=(coef[5]*coef[8]*coef[8] - coef[6])/(coef[7]*coef[7]);
          y:= if (x<=coef[8]) then 0 else (x-coef[8])/(x+coef[7]);
          y2:=y*y;
          y4:=y2*y2;
          h:= if (x<=coef[8]) then 0 else (coef[7]+coef[8])*((2*coef[4]+8*coef[5])*log(1-y)+ (coef[4]*(1+1/(1-y))+coef[5]*(7+1/(1-y)))*y + coef[5]*(3*y2+5*y*y2/3+y4 + 0.6*y4*y+y4*y2/3)+ (coef[5]-coef[6]/(coef[7]+coef[8])^2)*y4*y2*y/7);
          H:= R*x*(coef[1]+coef[2]*exp(-coef[3]/x)/(coef[3]*x))+R*h "Wilhoit.ok";
          
        elseif cor == 7 then
          H:=coef[1]*x+coef[2]/(coef[3]+1)*x^(coef[3]+1);
            for j in {4, 6, 8, 10} loop
              if (coef[j]>0) then
                H :=  H+coef[j]*coef[j+1] / (exp(coef[j+1]/x)-1);
              end if;
            end for;
            H:=H * R "Cooper.ok";
        
        elseif cor == 8 then
          H:=(1+coef[1])*x;
            for j in {2, 6} loop
              if (coef[j]>0) then 
                H := H+2*coef[j]*coef[j+1] / (exp(2*coef[j+1]/x)-1)+ 2*coef[j+2]*coef[j+3] / (exp(2*coef[j+3]/x)+1);
              end if;
            end for;
            H := H*R "Jaeschke.ok";
            
        elseif cor == 19 then
          Tm:=1-x/coef[6];
          Tm2:=Tm*Tm;
          H:=-coef[6]*(coef[1]*log(Tm)+coef[2]*Tm+0.5*coef[3]*Tm2+coef[4]*Tm2*Tm/3+0.25*coef[5]*Tm2*Tm2) "PPDS15.ok";
        end if;
        
        if cor == 1 then
          H:=H*1e3;
        elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
          H:=H*1e3/ MW;
        elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 then
          H:=H/ MW;
        elseif cor == 3 or cor == 10 then
          H:=H*4.1868*1e3/ MW;
        end if;
    /*    
        switch (cor){
        case 1://DIPPR 100 Cp0 in KJ/kgr·K
        case 2://DIPPR 100 Cp0in J/mol·K
        case 6://DIPPR 100 Cp0 in J/kgr·K
        case 15://DIPPR 100 Liquid Cp in J/mol·K
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 17://DIPPR 100 Liquid Cp in J/kgr·K
            *H=x*(coef[0]+x*(coef[1]/2+x*(coef[2]/3+x*(coef[3]/4+x*coef[4]/5))));//This is the integration from Cp0
            break;
        case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
            *H=x*(coef[0]+x*(coef[1]/2+x*(coef[2]/3+x*(coef[3]/4+x*(coef[4]/5+x*coef[5]/6)))));//This is the integration from Cp0
            break;
        case 3://DIPPR 107 correlation in calories/mol·K
        case 4://DIPPR 107 correlation in J/Kmol·K
    	case 200://DIPPR107 Cp in J/kg·K
            *H=coef[0]*x+coef[1]*coef[2]*(1/tanh(coef[2]/x))-coef[3]*coef[4]*tanh(coef[4]/x);
            break;*/
        /*
        case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
            int j=20;
            double interval;
            double T,Cp1,Cp2;
            for (i=0;i<*nPoints;i++){
                interval=x/i;
                Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
                for (j=1;j<21;j++){
                    T=*refT+j*interval;
                    Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                    th0->H=th0->H+(Cp1+Cp2)/2*interval;
                    th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                    Cp1=Cp2;
                    }
            }
            break;
        */
    /*    case 5:{//Wilhoit equation J/mol·K (8 coefficients)
            int j;
            double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
            a7_a6=coef[7]/coef[6];
            a7_a6_2=a7_a6*a7_a6;
            a7_a6_4=a7_a6_2*a7_a6_2;
            x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);
            if (x<=coef[7]) y=0;
            else y=(x-coef[7])/(x+coef[6]);
            y2=y*y;
            y4=y2*y2;
            if (x<=coef[7]) h=0;
            else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                    coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
            *H= R*x*(coef[0]+coef[1]*exp(-coef[2]/x)/(coef[2]*x))+R*h;
        }
            break;
        case 7:{//Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
            int j;*/
            /*
            th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
            for (i=3;i<13;i=i+2){
                if (coef>0) th0->Cp=th0->Cp+coef*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
            }
            th0->Cp=th0->Cp*R;*/
    /*        *H=coef[0]*x+coef[1]/(coef[2]+1)*pow(x,(coef[2]+1));
            for (j=3;j<13;j=j+2) if (coef[j]>0) *H =*H+coef[j]*coef[j+1] / (exp(coef[j+1]/x)-1);
            *H=*H * R;
        }
            break;
        case 8:{//Jaeschke and Schley equation (9 coefficients). Used by GERG2004
            int j;
            *H=(1+coef[0])*x;
            for (j=1;j<9;j=j+4) if (coef[j]>0) *H = *H+2*coef[j]*coef[j+1] / (exp(2*coef[j+1]/x)-1)+
                    +2*coef[j+2]*coef[j+3] / (exp(2*coef[j+3]/x)+1);
            *H = *H*R;
        }
            break;
        case 19:{//a/Tm+b+c*Tm+d*Tm^2+e*Tm^3 PPDS15 Cpl in J/(kgr·K)
            double Tm,Tm2;
            Tm=1-x/coef[5];
            Tm2=Tm*Tm;
            *H=-coef[5]*(coef[0]*log(Tm)+coef[1]*Tm+0.5*coef[2]*Tm2+coef[3]*Tm2*Tm/3+0.25*coef[4]*Tm2*Tm2);
          }
            break;
        }
        //now is necessary to pass to J/kg·K
        switch (cor){
        case 1://DIPPR 100 Cp0 in KJ/kgr·K
            *H=*H*1e3;
            break;
        case 2://DIPPR 100 Cp0in J/mol·K
        case 5://Wilhoit equation J/mol·K (8 coefficients)
        case 7://Cooper J/mol·K (11 coefficients)
        case 8://Jaeschke and Schley equation J/mol·K (9 coefficients)
        case 15://DIPPR 100 Liquid Cp in J/mol·K
            *H=*H*1e3/ MW;
            break;
        case 4://DIPPR 107 correlation in J/Kmol·K
        case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
            *H=*H/ MW;
            break;
        case 3://DIPPR 107 correlation in calories/mol·K
        case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
            *H=*H*4.1868*1e3/ MW;
            break;
        }*/
    end SpecificEnthalpyCorrCalc;
  
    function SpecificEntropyCorrCalc
      input Integer cor;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real S "Entropy J/(kg·K)";
    protected
      Real y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
      Real Tm,Tm2;
      Integer j;
    algorithm
        if cor == 1 or cor == 2 or cor == 6 or cor == 15 or cor == 16 or cor == 17 then
        S:=coef[1]*log(x)+x*(coef[2]+x*(coef[3]/2+x*(coef[4]/3+x*coef[5]/4))) "DIPPR100";
        
        elseif cor == 10 then 
        S:=coef[1]*log(x)+x*(coef[2]+x*(coef[3]/2+x*(coef[4]/3+x*(coef[5]/4+x*coef[6]/5)))) "Polynomial";
        
        elseif cor == 3 or cor == 4 or cor == 200 then
        S:=coef[1]*log(x)+coef[2]*(coef[3]/x/tanh(coef[3]/x)-log(sinh(coef[3]/x)))-coef[4]*(coef[5]/x*tanh(coef[5]/x)-log(cosh(coef[5]/x)))  "DIPPR107.ok";
        
        elseif cor == 5 then
          a7_a6:=coef[8]/coef[7];
          a7_a6_2:=a7_a6*a7_a6;
          a7_a6_4:=a7_a6_2*a7_a6_2;
          x1:=(coef[5]*coef[8]*coef[8] - coef[6])/(coef[7]*coef[7]);
          y:= if x<=coef[8] then 0 else (x-coef[8])/(x+coef[7]);
          y2:=y*y;
          y4:=y2*y2;
          if x<=coef[8] then
            s:=0;
          else
            z := x*(coef[8] + coef[7])/((x + coef[7])*coef[8]);
            w:=0;
            for j in 1:7 loop
              w:=w+(x1*(-a7_a6)^(6-j) - coef[5])*y^j/j;
            end for;
            s:=(coef[4] + ((coef[5]*coef[8]*coef[8]-coef[6])*a7_a6_4/(coef[7]*coef[7])))*a7_a6_2*log(z)+ (coef[4] + coef[5])*log((x + coef[7])/(coef[7] + coef[8]))- (coef[4]*(coef[7] + coef[8])/coef[7] + coef[6]*y4*y2/(7.*coef[8]*(coef[7] + coef[8])))*y+w;
          end if;
          S:=R*(coef[1]*log(x)+coef[2]*(1+coef[3]/x)*exp(-coef[3]/x)/(coef[3]*coef[3])+s) "Wilhoit.ok";
          
        elseif cor == 7 then
          S:= if (coef[2]>0) then coef[1]*log(x)+coef[2]/coef[3]*x^coef[3] else coef[1]*log(x);
          for j in {4, 6, 8, 10} loop
            if (coef[j]>0) then
              S := S+coef[j]*coef[j+1]*(exp(coef[j+1]/x)/(x*(exp(coef[j+1]/x)-1))- log(exp(coef[j+1]/x) -1)/coef[j+1]) "Cooper.ok";
            end if;
          end for;
          S:=S * R;
          
        elseif cor == 8 then
          S:=(1+coef[1])*log(x);
          for j in {2, 6} loop
            if (coef[j]>0) then
              S := S+coef[j]*(coef[j+1]/x/tanh(coef[j+1]/x)-log(sinh(coef[j+1]/x)))-coef[j+2]*(coef[j+3]/x*tanh(coef[j+3]/x)-log(cosh(coef[j+3]/x))) "Jaeschke.ok";
            end if;
          end for;
            
            //for (j=1;j<9;j=j+4) if (coef[j]>0) *S = *S+coef[j]*(coef[j+1]/x/tanh(coef[j+1]/x)-log(sinh(coef[j+1]/x)))-coef[j+2]*(coef[j+3]/x*tanh(coef[j+3]/x)-log(cosh(coef[j+3]/x)));
            S := S * R;
            
        elseif cor == 19 then
            Tm2:=coef[6]*coef[6];
            S:=-(6*Tm2*x*(coef[3]+2*coef[4]+3*coef[5])-3*coef[6]*x*x*(3*coef[5]+coef[4])+2*x*x*x*coef[5])/(6*Tm2*coef[6])+
                    (coef[1]+coef[2]+coef[3]+coef[4]+coef[5])*log(x)-coef[1]*log(abs(x-coef[6]));
        end if;
        
        if cor == 1 then
          S:=S*1e3;
        elseif cor == 2 or cor == 5 or cor == 7 or cor == 8 or cor == 15 then
          S:=S*1e3/ MW;
        elseif cor == 4 or cor == 9 or cor == 16 or cor == 18 then
          S:=S/ MW;
        elseif cor == 3 or cor == 10 then
          S:=S*4.1868*1e3/ MW;
        end if;
    
    /*    switch (cor){
        case 1://DIPPR 100 Cp0 in KJ/kgr·K
        case 2://DIPPR 100 Cp0in J/mol·K
        case 6://DIPPR 100 Cp0 in J/kgr·K
        case 15://DIPPR 100 Liquid Cp in J/mol·K
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 17://DIPPR 100 Liquid Cp in J/kgr·K
            *S=coef[0]*log(x)+x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*coef[4]/4)));//This is the integration from Cp0
            break;
        case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
            *S=coef[0]*log(x)+x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*(coef[4]/4+x*coef[5]/5))));
            break;
        case 3://DIPPR 107 correlation in calories/mol·K
        case 4://DIPPR 107 correlation in J/Kmol·K
    	case 200://DIPPR107 Cp in J/kg·K
            *S=coef[0]*log(x)+coef[1]*(coef[2]/x/tanh(coef[2]/x)-log(sinh(coef[2]/x)))-coef[3]*(coef[4]/x*tanh(coef[4]/x)-log(cosh(coef[4]/x)));
            break;*/
        /*
        case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
            int j=20;
            double interval;
            double T,Cp1,Cp2;
            for (i=0;i<*nPoints;i++){
                interval=x/i;
                Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
                for (j=1;j<21;j++){
                    T=*refT+j*interval;
                    Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                    th0->H=th0->H+(Cp1+Cp2)/2*interval;
                    th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                    Cp1=Cp2;
                    }
            }
            break;
        */
    /*    case 5:{//Wilhoit equation J/mol·K (8 coefficients)
            int j;
            double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
            a7_a6=coef[7]/coef[6];
            a7_a6_2=a7_a6*a7_a6;
            a7_a6_4=a7_a6_2*a7_a6_2;
            x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);
            if (x<=coef[7]) y=0;
            else y=(x-coef[7])/(x+coef[6]);
            y2=y*y;
            y4=y2*y2;
            if (x<=coef[7]) s=0;
            else{
                z = x*(coef[7] + coef[6])/((x + coef[6])*coef[7]);
                w=0;
                for (j=1;j<8;j++) w=w+(x1*pow(-a7_a6,6-j) - coef[4])*pow(y,j)/j;
                s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                    (coef[3] + coef[4])*log((x + coef[6])/(coef[6] + coef[7]))-
                    (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
                }
            *S=R*(coef[0]*log(x)+coef[1]*(1+coef[2]/x)*exp(-coef[2]/x)/(coef[2]*coef[2])+s);
        }
            break;
        case 7:{//Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
            int j;*/
            /*
            th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
            for (i=3;i<13;i=i+2){
                if (coef>0) th0->Cp=th0->Cp+coef*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
            }
            th0->Cp=th0->Cp*R;*/
    /*        if (coef[1]>0) *S=coef[0]*log(x)+coef[1]/coef[2]*pow(x,coef[2]);
            else *S=coef[0]*log(x);
            //printf("T:%f S0(0):%f\n",x,S*R/th0->MW);
            for (j=3;j<13;j=j+2){
                if (coef[j]>0) *S = *S+coef[j]*coef[j+1]*(exp(coef[j+1]/x)/(x*(exp(coef[j+1]/x)-1))-
                    log(exp(coef[j+1]/x)-1)/coef[j+1]);
                //printf("S0(%i):%f\n",i,S*R/th0->MW);
            }
            *S=*S * R;
        }
            break;
        case 8:{//Jaeschke and Schley equation (9 coefficients). Used by GERG2004
            int j;
            *S=(1+coef[0])*log(x);
            for (j=1;j<9;j=j+4) if (coef[j]>0) *S = *S+coef[j]*(coef[j+1]/x/tanh(coef[j+1]/x)-log(sinh(coef[j+1]/x)))-coef[j+2]*(coef[j+3]/x*tanh(coef[j+3]/x)-log(cosh(coef[j+3]/x)));
            *S = *S * R;
        }
            break;
        case 19:{//PPDS15 equation a/Tm+b+c*Tm+d*Tm^2+e*Tm^3 Cpl in J/(kgr·K)
            double Tc2;
            Tc2=coef[5]*coef[5];
            *S=-(6*Tc2*x*(coef[2]+2*coef[3]+3*coef[4])-3*coef[5]*x*x*(3*coef[4]+coef[3])+2*x*x*x*coef[4])/(6*Tc2*coef[5])+
                    (coef[0]+coef[1]+coef[2]+coef[3]+coef[4])*log(x)-coef[0]*log(abs(x-coef[5]));
          }
            break;
    		
        }
        //now is necessary to pass to J/kg·K
        switch (cor){
        case 1://DIPPR 100 Cp0 in KJ/kgr·K
            *S=*S*1e3;
            break;
        case 2://DIPPR 100 Cp0in J/mol·K
        case 5://Wilhoit equation J/mol·K (8 coefficients)
        case 7://Cooper J/mol·K (11 coefficients)
        case 8://Jaeschke and Schley equation J/mol·K (9 coefficients)
        case 15://DIPPR 100 Liquid Cp in J/mol·K
            *S=*S*1e3/ MW;
            break;
        case 4://DIPPR 107 correlation in J/Kmol·K
        case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
        case 16://DIPPR 100 Liquid Cp in J/Kmol·K
        case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
            *S=*S/ MW;
            break;
        case 3://DIPPR 107 correlation in calories/mol·K
        case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
            *S=*S*4.1868*1e3/ MW;
            break;
        }*/
    end SpecificEntropyCorrCalc;
  
    function PhysPropCorr "Calculates a physical property from a given correlation data and a given temperature, or pressure, using an external function. Should be the general form of calculation, but is not working in OpenModelica 1.14 when iterating in an array of DataRecord"
      input Integer corr;
      input Real coef[:];
      input Real MW;
      input Real x;
      output Real y;
    
      external "C" FF_PhysPropCorr(corr, coef, MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end PhysPropCorr;

    function Cp0Corr "Calculates ideal gas heat capacity from a given DataRecord and a given temperature, using an external function. The intermediate pass of the protected variable is necessary for OpenModelica 1.14"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real y;
    protected
      Real coef[13] = data.Cp0Coef;
    
      external "C" FF_PhysPropCorr(data.Cp0Corr, coef, data.MW, T, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end Cp0Corr;

    function FF_Cp0Corr "Calculates ideal gas heat capacity from a given DataRecord and a given temperature, using DIPPR107 correlation in J/(kg·K)"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real y;
    algorithm
      y := data.Cp0Coef[1] + data.Cp0Coef[2] * (data.Cp0Coef[3] / (T * sinh(data.Cp0Coef[3] / T))) ^ 2 + data.Cp0Coef[4] * (data.Cp0Coef[5] / (T * cosh(data.Cp0Coef[5] / T))) ^ 2;
    end FF_Cp0Corr;

    function SpecificEnthalpyCorr2 "Calculates specific enthalpy from a given DataRecord at a given temperature, using an external function.  The intermediate pass of the protected variable is necessary for OpenModelica 1.14"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real h;
    protected
      Real coef[13] = data.Cp0Coef;
    
      external "C" FF_SpecificEnthalpyCorr(data.Cp0Corr, coef, data.MW, T, h) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end SpecificEnthalpyCorr2;

    function SpecificEntropyCorr2 "Calculates specific enthalpy from a given DataRecord at a given temperature. As per specific enthalpy "
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real s;
    protected
      Real coef[13] = data.Cp0Coef;
    
      external "C" FF_SpecificEntropyCorr(data.Cp0Corr, coef, data.MW, T, s) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end SpecificEntropyCorr2;

    function liqViscPcorLucas "Lucas liquid viscosity pressure correction."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T, p, Pref, rpVisc;
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

    function gasViscCorr "Calculates low pressure gas viscosity from a given DataRecord and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gViscCoef directly, coming from an iteration."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real eta;
    protected
      Real coef[6] = data.gViscCoef;
    
      external "C" FF_PhysPropCorr(data.gViscCorr, coef, data.MW, T, eta) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end gasViscCorr;

    function FF_gasViscCorr "Calculates low pressure gas viscosity from a given DataRecord and a given temperature. Faster than gasViscCorr"
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real eta;
    algorithm
      eta := if data.gViscCorr == 111 then data.gViscCoef[1] + T * (data.gViscCoef[2] + T * (data.gViscCoef[3] + T * (data.gViscCoef[4] + T * data.gViscCoef[4]))) else data.gViscCoef[1] * T ^ data.gViscCoef[2] / (1 + data.gViscCoef[3] / T + data.gViscCoef[4] / T ^ 2);
    end FF_gasViscCorr;

    function gasViscLowPressureChung "Dynamic viscosity of a low pressure gas according to Chung"
      input FreeFluids.MediaCommon.DataRecord data;
      input SI.Temp_K T "Gas temperature";
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
    protected
      Real coef[6] = data.gThCondCoef;
    
      external "C" FF_PhysPropCorr(data.gThCondCorr, coef, data.MW, T, lambda) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end gasThCondCorr;

    function FF_gasThCondCorr "Calculates low pressure gas thermal conductivity from a given correlation data and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gThCondCoef directly, coming from an iteration."
      input FreeFluids.MediaCommon.DataRecord data;
      input Real T;
      output Real lambda;
    algorithm
      lambda := if data.gThCondCorr == 121 then data.gThCondCoef[1] + T * (data.gThCondCoef[2] + T * (data.gThCondCoef[3] + T * (data.gThCondCoef[4] + T * data.gThCondCoef[4]))) else data.gThCondCoef[1] * T ^ data.gThCondCoef[2] / (1 + data.gThCondCoef[3] / T + data.gThCondCoef[4] / T ^ 2);
    end FF_gasThCondCorr;

    function gasThCondLowPressureChung "Thermal conductivity of a low pressure gas according to Chung"
      input FreeFluids.MediaCommon.DataRecord data;
      input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp0;
      input SI.DynamicViscosity eta "low pressure dynamic viscosity";
      input SI.Temp_K T "Gas temperature";
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

    function gasMixViscosityWilke
      input SI.Temp_K T;
      input SI.MoleFraction X[:];
      input SI.DynamicViscosity etaX[:];
      input SI.MolarMass MW[:];
      output SI.DynamicViscosity eta;
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
      input SI.Temp_K T;
      input SI.MoleFraction X[:];
      input SI.ThermalConductivity lambdaX[size(X, 1)];
      input SI.MolarMass MW[size(X, 1)];
      input SI.Temp_K Tc[size(X, 1)];
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
      input SI.Temp_K T;
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
      input SI.Temp_K T;
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

end MediaCommon;

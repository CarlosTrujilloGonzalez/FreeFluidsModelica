within FreeFluids;

package ExternalMix "ExternalMix by Carlos Trujillo
  This file is part of the Free Fluids application
  Copyright (C) 2008-2024  Carlos Trujillo Gonzalez
    
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
  //***** PACKAGE ExternalMixMedium*****
  //*********************************

  package ExternalMixMedium
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(onePhase = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, reference_T = 298.15, reference_p = 101325.0, mediumName = "None",nX=subsNum);
    constant String subsNames = "None";
    constant Integer subsNum = Modelica.Utilities.Strings.count(subsNames,",")+1;
    constant String resDir = Modelica.Utilities.Files.loadResource("modelica://FreeFluids/Resources") "resources directory";
    constant String eosType = "PR" "alternatives are SRK and PCSAFT";
    constant String mixRule = "LCVM" "alternatives are VdW, VdWnoInt, Reid, MKP, HV, MHV1, MHV2, UMR, PSRK, BL, IndAssoc";
    constant String activityModel = "UNIFACdort" "alternatives are: None, UNIFACstd, UNIFACpsrk, UNIQUAC, NRTL";
    constant String viscMixRule="Teja" "alternatives are: Grunberg, Andrade, McAllister3, McAllister4";
  
    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        AbsolutePressure p(displayUnit = "bar") "Pressure in Pa";
        Temperature T(displayUnit = "K") "Kelvin temperature of medium";
        SpecificEnthalpy h "Overall specific enthalpy";
        SpecificEntropy s "Overall specific entropy";
        Density d(displayUnit = "kg/m3") "Overall density in kgr/m3";
        MassFraction gf "gas mass fraction";
        MassFraction z[subsNum] "global mass fraction of each substance";
        MassFraction x[subsNum] "liquid mass fraction of each substance";
        MassFraction y[subsNum] "gas mass fraction of each substance";
        MolarMass MW;
        MolarMass lMW;
        MolarMass gMW;
        Density ld(displayUnit = "kg/m3") "Liquid phase density in kgr/m3";
        Density gd(displayUnit = "kg/m3") "Gas phase density in kgr/m3";
        SpecificEnthalpy lh "liquid phase specific enthalpy in J/kgr";
        SpecificEnthalpy gh "gas phase specific enthalpy in J/kgr";
        SpecificEntropy ls "liquid phase specific entropy";
        SpecificEntropy gs "gas phase specific entropy";
        SpecificHeatCapacity lCv "liquid phase specific heats at constant volume";
        SpecificHeatCapacity lCp "liquid phase specific heat at constant pressure";
        SpecificHeatCapacity gCv "gas phase specific heats at constant volume";
        SpecificHeatCapacity gCp "gas phase specific heat at constant pressure";
        Real lDvp "liquid phase derivative of 1/density regarding pressure";
        Real lDvT "liquid phase derivative of 1/density regarding temperature";
        Real gDvp "gas phase derivative of 1/density regarding pressure";
        Real gDvT "gas phase derivative of 1/density regarding temperature";
    end ThermodynamicState;
  
    pure function molarMasses "return molar masses in kg/mol of the medium substances"
      output Real MM[subsNum];
    
      external "C" FF_molarMassesM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, MM) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end molarMasses;
  
    pure function pressureEOS_dTX "Return pressure given temperature and density, by EOS"
      input Density d;
      input Temperature T;
      input Real z[subsNum];
      output AbsolutePressure p;
      output MolarMass MW;
      external "C" FF_pressure_dTXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, d, T, z, p, MW) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end pressureEOS_dTX;
  
    pure function density_pTX "Return density at give temperature and pressure, by EOS"
      input AbsolutePressure p;
      input Temperature T;
      input Real z[subsNum];
      input String var;
      //input String var "b"=both phases, "g"= only gas, "l"=only liquid";
      output Density ld;
      output Density gd;
      output MolarMass MW;
      external "C" FF_density_pTXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, var, p, T, z, ld, gd, MW) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end density_pTX;
  
    pure function thermo_dTX "Return thermodynamic properties at given temperature, density and composition, by EOS"
      input Density d;
      input Temperature T;
      input Real z[subsNum];
      output Real p, h, s, Cv, Cp, Dvp, DvT;
      external "C" FF_thermo_dTXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, d, T, z, p, h, s, Cv, Cp, Dvp, DvT)  annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end thermo_dTX;
  
    pure function activityCoefficients_TX "Return activity coefficients and excess Gibbs energy at given temperature and composition"
      input Temperature T;
      input Real z[subsNum];
      output Real actCoef[subsNum];
      output Modelica.Units.SI.MolarEnergy gE;
      external "C" FF_activity_TXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, T, z, actCoef, gE) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end activityCoefficients_TX;
  
    pure function fugacityCoefficients_pTX "Return fugacity coefficients, calculated from liquid and gas extremes, and Wilson k value, at given pressure, temperature and composition"
      input AbsolutePressure p;
      input Temperature T;
      input Real z[subsNum];
      output Real fugCoefL[subsNum];
      output Real fugCoefG[subsNum];
      output Real kWilson[subsNum];
      external "C" FF_fugacity_pTXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, z, fugCoefL, fugCoefG, kWilson) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end fugacityCoefficients_pTX;
  
    pure function bubble_TX "Return bubble pressure and gas composition at equilibrium"
      input Temperature T;
      input Real x[subsNum];
      output AbsolutePressure p;
      output Real y[subsNum];
      external "C" FF_bubble_TXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, T, x, p, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end bubble_TX;
  
    pure function bubble_pX "Return bubble pressure and gas composition at equilibrium"
      input AbsolutePressure p;
      input Real x[subsNum];
      output Temperature T;
      output Real y[subsNum] "gas phase composition at equilibrium";
      external "C" FF_bubble_pXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, x, T, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end bubble_pX;
  
    pure function dew_TX "Return dew pressure and liquid composition at equilibrium"
      input Temperature T;
      input Real y[subsNum];
      output Real p;
      output Real x[subsNum];
      external "C" FF_dew_TXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, T, y, p, x) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end dew_TX;
  
    pure function dew_pX "Return dew temperature and liquid composition at equilibrium"
      input AbsolutePressure p;
      input Real y[subsNum];
      output Temperature T;
      output Real x[subsNum];
      external "C" FF_dew_pXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, y, T, x) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end dew_pX;
  
    pure function twoPhasesFlash_pTX "Return liquid and gas composition at equilibrium, plus the gas fraction"
      input AbsolutePressure p;
      input Temperature T;
      input Real z[subsNum];
      output Real x[subsNum];
      output Real y[subsNum];
      output Real gf;
      external "C" FF_TwoPhasesFlashPTXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, z, x, y, gf) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end twoPhasesFlash_pTX;
    
    pure function twoPhasesFlash_phsX "Return liquid and gas composition at equilibrium, plus temperature and gas fraction"
      input String opt "indicates if e is enthalpy (h) or entropy(s)";
      input AbsolutePressure p;
      input Real e;
      input Real z[subsNum];
      output Temperature T;
      output Real x[subsNum];
      output Real y[subsNum];
      output Real gf;
      external "C" FF_TwoPhasesFlashP_HS_XM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, opt, p, e, T, z, x, y, gf) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFmodelicaMedium.c\"");
    end twoPhasesFlash_phsX;
  
    pure function twoPhasesFlash_pThetaX "Return liquid and gas composition at equilibrium, plus temperature and gas fraction"
      input String opt "indicates if theta is in molar or mas basis";
      input AbsolutePressure p;
      input SpecificEnthalpy theta;
      input Real z[subsNum];
      output Temperature T;
      output Real x[subsNum];
      output Real y[subsNum];
  
      external "C" FF_TwoPhasesFlashPThetaXM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, opt, p, theta, T, z, x, y) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFmodelicaMedium.c\"");
    end twoPhasesFlash_pThetaX;
  
    pure function liquidDynamicViscosityAux "Return liquid dynamic viscosity from p,T,X"
      input Real p;
      input Real T;
      input Real x[:];
      output DynamicViscosity eta "Dynamic viscosity";
      external "C" FF_mixLiquidViscosityM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, x, eta) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFmodelicaMedium.c\"");
    end liquidDynamicViscosityAux;
  
    pure function gasDynamicViscosityAux "Return gas dynamic viscosity from a p,T,X"
      input Real p;
      input Real T;
      input Real y[:];
      output DynamicViscosity eta "Dynamic viscosity";
      external "C" FF_mixGasViscosityM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, y, eta) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFmodelicaMedium.c\"");
    end gasDynamicViscosityAux;
   
    pure function liquidThermalConductivityAux "Return liquid thermal conductivity from a p,T,X"
      input Real p;
      input Real T;
      input Real x[:];
      output ThermalConductivity lambda "Thermal conductivity";
      external "C" FF_mixLiquidThCondM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, x, lambda) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFmodelicaMedium.c\"");
    end liquidThermalConductivityAux;
  
    pure function gasThermalConductivityAux "Return liquid thermal conductivity from a p,T,X"
      input Real p;
      input Real T;
      input Real y[:];
      output ThermalConductivity lambda "Thermal conductivity";
      external "C" FF_mixGasThCondM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, y, lambda) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFmodelicaMedium.c\"");
    end gasThermalConductivityAux;
  
    pure function surfaceTensionAux "Return liquid thermal conductivity from a p,T,X"
      input Real p;
      input Real T;
      input Real x[:];
      output SurfaceTension sigma "surface tension";
      external "C" FF_mixSurfaceTensionM(mediumName, subsNum, subsNames, resDir, eosType, mixRule, activityModel, viscMixRule, p, T, x, sigma) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFmodelicaMedium.c\"");
    end surfaceTensionAux;
  
    function massToMoleFractions "Return mole fractions from mass fractions X"
      extends Modelica.Icons.Function;
      input SI.MassFraction X[subsNum] "Mass fractions of mixture";
      input SI.MolarMass MMX[subsNum] "Molar masses of components";
      output SI.MoleFraction moleFractions[subsNum] "Mole fractions of mixture";
    protected
      Real invMMX[subsNum] "Inverses of molar weights";
      SI.MolarMass Mmix "Molar mass of mixture";
      Real nMoles;
    algorithm
      nMoles:=0;
//for i in 1:size(X, 1) loop
      for i in 1:subsNum loop
        invMMX[i] := 1/MMX[i];
        nMoles:=nMoles+X[i]/MMX[i];
      end for;
//Mmix := 1/(X*invMMX);
      if (nMoles>0)then
        Mmix:=1/nMoles;
        for i in 1:size(X, 1) loop
          moleFractions[i] := Mmix*X[i]/MMX[i];
        end for;
      else
        Mmix:=0;
        for i in 1:size(X, 1) loop
          moleFractions[i] := 0;
        end for;
      end if;
    end massToMoleFractions;
  
    function moleToMassFractions "Return mole fractions from mass fractions X"
      extends Modelica.Icons.Function;
      input SI.MoleFraction moleFractions[subsNum] "Mole fractions of mixture";
      input MolarMass MMX[subsNum] "Molar masses of components";
      output SI.MassFraction X[subsNum] "Mass fractions of gas mixture";
    protected
      SI.MolarMass Mmix "Molar mass of mixture";
    algorithm
      Mmix := moleFractions*MMX;
      for i in 1:size(moleFractions, 1) loop
        X[i] := moleFractions[i]*MMX[i]/Mmix;
      end for;
    end moleToMassFractions;
  
    function setLiquidState_pTX "Return ThermodynamicState record as function of p,T and composition X with only liquid"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[nX]  "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    algorithm
      state.phase := 0 "phase has not been checked";
      state.gf := 0;
      state.p := p;
      state.T := T;
      state.z := X;
      state.x := X;
      state.y := fill(0.0, subsNum);
      (state.ld, state.gd, state.lMW) := density_pTX(p, T, X, "l");
      state.d := state.ld;
      state.gd := 0;
      state.MW := state.lMW;
      state.gMW := 0;
      (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
      state.h := state.lh;
      state.gh := 0;
      state.s := state.ls;
      state.gs := 0;
      state.gCv := 0;
      state.gCp := 0;
      state.gDvp := 0;
      state.gDvT := 0;
    end setLiquidState_pTX;
  
    function setGasState_pTX "Return ThermodynamicState record as function of p,T and composition X with only liquid"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[nX] "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    algorithm
      state.phase := 0 "phase has not been checked";
      state.gf := 1;
      state.p := p;
      state.T := T;
      state.z := X;
      state.y := X;
      state.x := fill(0.0, subsNum);
      (state.ld, state.gd, state.gMW) := density_pTX(p, T, X, "g");
      state.d := state.gd;
      state.ld := 0;
      state.MW := state.gMW;
      state.lMW := 0;
      (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
      state.h := state.gh;
      state.lh := 0;
      state.s := state.gs;
      state.ls := 0;
      state.lCv := 0;
      state.lCp := 0;
      state.lDvp := 0;
      state.lDvT := 0;
    end setGasState_pTX;
  
    function setBubbleState_TX "Return bubble state ThermodynamicState record as function of T and composition X"
      extends Modelica.Icons.Function;
      input Temperature T "Temperature";
      input MassFraction X[nX] "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    protected
      Real p, dl, dg;
    algorithm
      state.phase := 1 "there is just 1 liquid phase";
      state.gf := 0;
      state.T := T;
      state.z := X;
      state.x := X;
      (state.p, state.y) := bubble_TX(T, X);
      //p := state.p;
      (state.ld, dg, state.lMW) := density_pTX(state.p, T, X, "l");
      state.d := state.ld;
      state.MW := state.lMW;
      (p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
      state.h := state.lh;
      state.s := state.ls;
      (dl, state.gd, state.gMW) := density_pTX(state.p, T, state.y, "g");
      (p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
    end setBubbleState_TX;
  
    function setBubbleState_pX "Return bubble state ThermodynamicState record as function of T and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "pressure";
      input MassFraction X[nX] "liquid phase Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    protected
      Real p2, dl, dg;
    algorithm
      state.phase := 1 "there is just 1 liquid phase";
      state.gf := 0;
      state.p := p;
      state.z := X;
      state.x := X;
      (state.T, state.y) := bubble_pX(p, X);
      (state.ld, dg, state.lMW) := density_pTX(p, state.T, X, "l");
      state.d := state.ld;
      state.MW := state.lMW;
      (p2, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
      state.h := state.lh;
      state.s := state.ls;
      (dl, state.gd, state.gMW) := density_pTX(p, state.T, state.y, "g");
      (p2, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
      state.p := p;
    end setBubbleState_pX;
  
    function setDewState_TX "Return dew state ThermodynamicState record as function of T and composition X"
      extends Modelica.Icons.Function;
      input Temperature T "Temperature";
      input MassFraction y[nX] "Mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    protected
      Real p, dl, dg;
    algorithm
      state.phase := 1 "there is just 1 gas phase";
      state.gf := 1;
      state.T := T;
      state.z := y;
      state.y := y;
      (state.p, state.x) := dew_TX(T, y);
      //p := state.p;
      (state.ld, dg, state.lMW) := density_pTX(state.p, T, state.x, "l");
      (p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
      (dl, state.gd, state.gMW) := density_pTX(p, T, state.y, "g");
      state.d := state.gd;
      state.MW := state.gMW;
      (p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
      state.h := state.gh;
      state.s := state.gs;
    end setDewState_TX;
  
    function setDewState_pX "Return dew state ThermodynamicState record as function of p and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "pressure";
      input MassFraction y[nX] "gas mass fractions";
      output ThermodynamicState state "Thermodynamic state record";
    protected
      Real p2,dl, dg;
    algorithm
      state.phase := 1 "there is just 1 gas phase";
      state.gf := 1;
      state.p := p;
      state.z := y;
      state.y := y;
      (state.T, state.x) := dew_pX(p, y);
      (state.ld, dg, state.lMW) := density_pTX(p, state.T, state.x, "l");
      (p2, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
      (dl, state.gd, state.gMW) := density_pTX(p, state.T, state.y, "g");
      state.d := state.gd;
      state.MW := state.gMW;
      (p2, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
      state.h := state.gh;
      state.s := state.gs;
      state.p := p;
    end setDewState_pX;
  
    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X"
        extends Modelica.Icons.Function;
  
      protected
        Real dl, dg;
  
      algorithm
        state.T := T;
        state.p := p;
        state.z := X;
        (state.x, state.y, state.gf) := twoPhasesFlash_pTX(p, T, X);
        if (state.gf == 0) then
          state.phase := 1 "only liquid phase";
          (state.ld, dg, state.lMW) := density_pTX(p, T, X, "l");
          state.d := state.ld;
          state.gd := 0;
          state.MW := state.lMW;
          state.gMW := 0;
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
          state.h := state.lh;
          state.gh := 0;
          state.s := state.ls;
          state.gs := 0;
          state.gCv := 0;
          state.gCp := 0;
          state.gDvp := 0;
          state.gDvT := 0;
        elseif (state.gf == 1) then
          state.phase := 1 "only gas phase";
          (dl, state.gd, state.gMW) := density_pTX(p, T, X, "g");
          state.d := state.gd;
          state.ld := 0;
          state.MW := state.gMW;
          state.lMW := 0;
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
          state.h := state.gh;
          state.lh := 0;
          state.s := state.gs;
          state.ls := 0;
          state.lCv := 0;
          state.lCp := 0;
          state.lDvp := 0;
          state.lDvT := 0;
        else
          state.phase := 2 "liquid and gas phases";
          (state.ld, dg, state.lMW) := density_pTX(p, T, state.x, "l");
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
          (dl, state.gd, state.gMW) := density_pTX(p, T, state.y, "g");
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
          state.MW:=1/((1-state.gf)/state.lMW+state.gf/state.gMW);
//state.MW := state.lMW*(1 - state.gf) + state.gMW*state.gf;
          state.d := state.ld*(1 - state.gf) + state.gd*state.gf;
          state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
          state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        end if;
    end setState_pTX;
  
    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,h and composition X"
        extends Modelica.Icons.Function;
  
      protected
        Real dl, dg;
  
      algorithm
        state.p := p;
        state.h := h;
        state.z := X;
        (state.T,state.x, state.y, state.gf) := twoPhasesFlash_phsX("h",p, h, X);
        if (state.gf == 0) then
          state.phase := 1 "only liquid phase";
          (state.ld, dg, state.lMW) := density_pTX(p, state.T, X, "l");
          state.d := state.ld;
          state.gd := 0;
          state.MW := state.lMW;
          state.gMW := 0;
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
//state.h := state.lh;
          state.gh := 0;
          state.s := state.ls;
          state.gs := 0;
          state.gCv := 0;
          state.gCp := 0;
          state.gDvp := 0;
          state.gDvT := 0;
        elseif (state.gf == 1) then
          state.phase := 1 "only gas phase";
          (dl, state.gd, state.gMW) := density_pTX(p, state.T, X, "g");
          state.d := state.gd;
          state.ld := 0;
          state.MW := state.gMW;
          state.lMW := 0;
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
//state.h := state.gh;
          state.lh := 0;
          state.s := state.gs;
          state.ls := 0;
          state.lCv := 0;
          state.lCp := 0;
          state.lDvp := 0;
          state.lDvT := 0;
        else
          state.phase := 2 "liquid and gas phases";
          (state.ld, dg, state.lMW) := density_pTX(p, state.T, state.x, "l");
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
          (dl, state.gd, state.gMW) := density_pTX(p, state.T, state.y, "g");
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
          state.MW:=1/((1-state.gf)/state.lMW+state.gf/state.gMW);
//state.MW := state.lMW*(1 - state.gf) + state.gMW*state.gf;
          state.d := state.ld*(1 - state.gf) + state.gd*state.gf;
          state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
          state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        end if;
    end setState_phX;
  
    redeclare function extends setState_psX "Return ThermodynamicState record as function of p,s and composition X"
        extends Modelica.Icons.Function;
  
      protected
        Real dl, dg;
  
      algorithm
        state.p := p;
        state.s := s;
        state.z := X;
        (state.T,state.x, state.y, state.gf) := twoPhasesFlash_phsX("s",p, s, X);
        if (state.gf == 0) then
          state.phase := 1 "only liquid phase";
          (state.ld, dg, state.lMW) := density_pTX(p, state.T, X, "l");
          state.d := state.ld;
          state.gd := 0;
          state.MW := state.lMW;
          state.gMW := 0;
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
//state.h := state.lh;
          state.gs := 0;
          state.h := state.lh;
          state.gh := 0;
          state.gCv := 0;
          state.gCp := 0;
          state.gDvp := 0;
          state.gDvT := 0;
        elseif (state.gf == 1) then
          state.phase := 1 "only gas phase";
          (dl, state.gd, state.gMW) := density_pTX(p, state.T, X, "g");
          state.d := state.gd;
          state.ld := 0;
          state.MW := state.gMW;
          state.lMW := 0;
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
//state.h := state.gh;
          state.ls := 0;
          state.h := state.gh;
          state.lh := 0;
          state.lCv := 0;
          state.lCp := 0;
          state.lDvp := 0;
          state.lDvT := 0;
        else
          state.phase := 2 "liquid and gas phases";
          (state.ld, dg, state.lMW) := density_pTX(p, state.T, state.x, "l");
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
          (dl, state.gd, state.gMW) := density_pTX(p, state.T, state.y, "g");
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
          state.MW:=1/((1-state.gf)/state.lMW+state.gf/state.gMW);
//state.MW := state.lMW*(1 - state.gf) + state.gMW*state.gf;
          state.d := state.ld*(1 - state.gf) + state.gd*state.gf;
          state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
          state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        end if;
    end setState_psX;
  
    function setState_pThetaX "Return saturated ThermodynamicState record as function of p,gas fraction and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input MassFraction theta "gas mass fraction";
      input MassFraction X[:]=reference_X "Mass fractions of components";
      output ThermodynamicState state "Thermodynamic state record";
  
      protected
        Real dl, dg;
      algorithm
        if (theta==1.0) then state:=setDewState_pX(p,X);
        elseif (theta==0.0) then state:=setBubbleState_pX(p,X);
        else
          state.p := p;
          state.z := X;
          state.gf:=theta;
          (state.T, state.x, state.y):=twoPhasesFlash_pThetaX("mass",p,theta,X);
          state.phase := 2 "liquid and gas phases";
          (state.ld, dg, state.lMW) := density_pTX(p, state.T, state.x, "l");
          (state.p, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := thermo_dTX(state.ld, state.T, state.x);
          (dl, state.gd, state.gMW) := density_pTX(p, state.T, state.y, "g");
          (state.p, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT) := thermo_dTX(state.gd, state.T, state.y);
          state.MW:=1/((1-state.gf)/state.lMW+state.gf/state.gMW);
//state.MW := state.lMW*(1 - state.gf) + state.gMW*state.gf;
          state.d := state.ld*(1 - state.gf) + state.gd*state.gf;
          state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
          state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        end if;
    end setState_pThetaX;
  
    redeclare function extends pressure "Return pressure"
        extends Modelica.Icons.Function;
      algorithm
        p := state.p;
    end pressure;
  
    redeclare function extends temperature "Return temperature"
        extends Modelica.Icons.Function;
      algorithm
        T := state.T;
    end temperature;
  
    redeclare function extends density "Return density"
        extends Modelica.Icons.Function;
      algorithm
        d := state.d;
    end density;
  
    redeclare function extends specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
      algorithm
        h := state.h;
    end specificEnthalpy;
  
    redeclare function extends specificInternalEnergy "Return specific internal energy"
        extends Modelica.Icons.Function;
      algorithm
        u := state.h - state.p/state.d;
    end specificInternalEnergy;
  
    redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
      algorithm
        s := state.s;
    end specificEntropy;
  
    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
      algorithm
        g := state.h - state.T*state.s;
    end specificGibbsEnergy;
  
    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
      algorithm
        f := state.h - state.p/state.d - state.T*state.s;
    end specificHelmholtzEnergy;
  
    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
      algorithm
        cp := if state.gf == 0 then state.lCp else if state.gf == 1 then state.gCp else 0.0;
    end specificHeatCapacityCp;
    
    redeclare function extends specificHeatCapacityCv "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
      algorithm
        cv := if state.gf == 0 then state.lCv else if state.gf == 1 then state.gCv else 0.0;
    end specificHeatCapacityCv;
  
    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
      algorithm
        a := if state.gf == 0.0 then (-state.lCp/(state.lDvp*state.lCv))^0.5/state.ld else if state.gf == 1.0 then (-state.gCp/(state.gDvp*state.gCv))^0.5/state.gd else 0.0;
    end velocityOfSound;
    
    redeclare function extends isobaricExpansionCoefficient "Returns the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
        algorithm
          beta := if state.gf == 0 then state.lDvT*state.ld else if state.gf == 1 then state.gDvT*state.gd else Modelica.Constants.small;
    end isobaricExpansionCoefficient;
    
    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
    algorithm
        kappa := if state.gf ==0 then -state.lDvp*state.ld else if state.gf ==1 then -state.gDvp*state.gd else Modelica.Constants.inf;
    end isothermalCompressibility;
  
    function specificVaporizationHeat "Return bubble state ThermodynamicState record as function of T and composition X"
      extends Modelica.Icons.Function;
      input ThermodynamicState state;
      output SpecificEnergy hv;
    algorithm
      hv := state.gh - state.lh;
    end specificVaporizationHeat;
  
    redeclare function molarMass "Return the molar mass of the medium"
      extends Modelica.Icons.Function;
      input ThermodynamicState state;
      output Real MM;
    algorithm
      MM := state.MW;
    end molarMass;
  
    function dynamicViscosity "Return dynamic viscosity from a ThermodynamicState record"
      extends Modelica.Icons.Function;
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    algorithm
      if (state.gf==0) then eta:=liquidDynamicViscosityAux(state.p,state.T,state.x);
      elseif (state.gf==1) then eta:=gasDynamicViscosityAux(state.p,state.T,state.y);
      else eta:=0;
      end if; 
    annotation(
        Documentation(info = "<html><head></head><body>Uses the Teja-Rice method for the liquid phase and the Lucas one for the gas phase. The viscosity is set to 0 in biphasic systems.</body></html>"));
    end dynamicViscosity;
    
    function liquidDynamicViscosity "Return liquid dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    algorithm
      eta:=liquidDynamicViscosityAux(state.p,state.T,state.x);
    annotation(
        Documentation(info = "<html><head></head><body>Is useful in obtaining the viscosity of the liquid phase in a biphasic system</body></html>"));
    end liquidDynamicViscosity;
  
    function gasDynamicViscosity "Return gas dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    algorithm
      eta:=gasDynamicViscosityAux(state.p,state.T,state.y);
    annotation(
        Documentation(info = "<html><head></head><body>Is useful in obtaining the viscosity of the gas phase in a biphasic system</body></html>"));
    end gasDynamicViscosity;
  
    function thermalConductivity "Return thermal conductivity from a ThermodynamicState record"
      extends Modelica.Icons.Function;
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity lambda "Thermal conductivity";
    algorithm
      if (state.gf==0) then lambda:=liquidThermalConductivityAux(state.p,state.T,state.x);
      elseif (state.gf==1) then lambda:=gasThermalConductivityAux(state.p,state.T,state.y);
      else lambda:=0;
      end if; 
    annotation(
        Documentation(info = "<html><head></head><body>Uses the Li method for the liquid phase and the Mason one for the gas phase. The viscosity is set to 0 in biphasic systems.</body></html>"));
    end thermalConductivity;
  
    function liquidThermalConductivity "Return liquid thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    algorithm
      lambda := liquidThermalConductivityAux(state.p, state.T, state.x);
      annotation(
        Documentation(info = "<html><head></head><body>Is useful in obtaining the thermal conductivity of the liquid phase in a biphasic system</body></html>"));
    end liquidThermalConductivity;
    
    function gasThermalConductivity "Return liquid thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    algorithm
      lambda := gasThermalConductivityAux(state.p, state.T, state.y);
      annotation(
        Documentation(info = "<html><head></head><body>Is useful in obtaining the thermal conductivity of the gas phase in a biphasic system</body></html>"));
    end gasThermalConductivity;
    
    function surfaceTension "Return liquid surface tension from a ThermodynamicState record"
      extends Modelica.Icons.Function;
      input ThermodynamicState state "Thermodynamic state record";
      output SurfaceTension sigma "Thermal conductivity";
    algorithm
      if (state.gf==0) then sigma := surfaceTensionAux(state.p, state.T, state.x);
      else sigma:=0;
      end if;
      annotation(
        Documentation(info = "<html><head></head><body>Uses the Winterfeld method.</body></html>"));
    end surfaceTension; 
    
    annotation(
      Documentation(info = "<html><head></head><body><br></body></html>"));
  end ExternalMixMedium;
  
      annotation(
    Documentation(info = "<html><head></head><body><div><b>Introduction</b></div>The medium is designed for substance mixtures in liquid, gas or biphasic (VL) state. It uses external C functions for the calculations, and it follows, when possible, the masterlines of Modelica. Media.Interfaces.PartialTwoPhaseMedium.<div>The implementation is based in a static object created in the C code, as it has been impossible to make it work, at least in OpenModelica, using the Modelica external object feature.</div><div>The data for the pure, or pseudo pure, substances are in the Resources/Fluids directory.</div><div>At the moment the implementation is still in development and few testing has been done, nevertheless bubble and dew calculations and the pTX, phX and psX flashes are already working.</div><div>For the thermodynamic models, you can use cubic (Peng-Robinson or SRK) &nbsp;or PCSAFT EOS. For PCSAFT the Wolbach and Sandler mixing rule is always used, but you can choose to use induced association by selecting induced association as mixrule. For cubics you can choose between Van der Waals, Panagiotopoulos-Reid, Mathias-Klotz-Prausnitz, HV, MHV1, MHV2, LCVM, UMR and PSRK mixing rules. &nbsp;When using gE mixing rules you can choose between UNIFAC (standard, PSRK, or Dortmund, implementations), UNIQUAC, and NRTL, for the activity model.</div><div><br></div><div><b>Media configuration</b></div><div>It is very easy. You extend the the ExternalMixMedium package and reconfigure it as follows:</div><div>mediumName: A string containing the medium name. In a simulation with more than one mediums you can't repeat the name.</div><div>subsNames: Is a string containing the names of the substances in the mix separated by ,. These names must match the name of files in the Resources/Fluids folder, but without the .sd extension. For the creation of these files with the aid of the FreeFluidsGui program, look at the ExternalPure package information.</div><div>eosType: A string. You can use \"PR\", \"SRK\" or \"PCSAFT\".</div><div>mixRule: A string. If a cubic EOS has been selected, you can use \"VdW\", \"VdWnoInt\", \"Reid\", \"MKP\", \"HV\", \"MHV1\", \"MHV2\", \"LCVM\", \"UMR\", \"PSRK\". Where VdWnoInt means Van der Waals without using BIPs. If PCSAFT is selected as EOS, the mixing rule used will be always that of Berthelot-Lorentz modified by Wolbach and Sandler, but by selecting \"IndAssoc\" as mixing rule you will allow the use of induced association for not self associating polar molecules.</div><div>activityModel: A string. If a gE mixing rule is used, you can use \"UNIFACstd\", \"UNIFACpsrk\" \"UNIFACdort\", \"NRTL\" or \"UNIQUAC\"</div><div>viscMixRule: A string where you can use: \"Teja\", \"Grunberg, \"Andrade\", \"McAllister3\" and \"McAllister4\"</div><div><br></div><div><b>EOS</b></div><div>Only Peng--Robinson and Soave-Redlich-Kwong cubic EOS are supported, but with many different alpha term calculations. You can use Van der Waals mixing rule, or gE ones, using an activity model.</div><div>For SAFT only PCSAFT is supported. Association and polar terms are implemented, the last one according to Gross and Vrabeck &nbsp;model. The Berthelot-Lorentz, as modified by Wolbach and Sandler, mixing rule is used. By making mixRue=\"IndAssoc\" you will allow the use of induced association between self associating and polar molecules, according to the methodology of Kleiner et alt.</div><div><br></div><div><b>Binary interaction parameters for EOS or activity</b></div><div>BIPs are read from text files (space delimited) placed in the Resources/Interactions folder. A lot of BIPs are supplied. Those for NRTL and UNIQUAC come mainly from the PychemQt open source program. For Van der Waals mixing rule, there are separate files for Peng-Robinson and SRK.</div><div>It is recommended that you create your own BIP files, containing only the needed parameters, because if a big file is used you loose speed and, as the program search all the file, you can finish using undesired BIPs. You can rename the supplied files and create new ones containing just the needed BIPs.</div><div>CAS numbers are used in order to seek the needed line inside the files.</div><div>For Van der Waals mixing rule there are 4 parameters, the first 3 are for the calculation of Kij in the form Kij=a+b*T+c/T. the last one is for Lij. As Kij=Kji and Lij=Lji just one entry for a substances pair is needed.</div><div>For PCSAFT we have the same situation, just without the Lij parameter. For this EOS it is very important to check that the BIPs are for the variation of the EOS used. They are not the same for a non-associating EOS than for a associating one, or a polar one.</div><div>For NRTL and UNIQUAC, prior to the BIPs there is the information of the equation that will be used with these parameters. The alternatives are Pol1K, Pol1J, Pol1C, Pol2K, Pol2J, Pol2C, Pol3K, Pol3J, Pol3K. &nbsp;the last letter identifies the units of energy used (K, J/mol or cal/mol). Pol1 will use the equation a+b*T+c*T^2, Pol2 a+b*T+c/T^2, and Pol3 a+b/T+c*T. As the BIPs are not simetrical the line contains both ij and ji parameters. In NRTL the fith BIP is for alphaij. There is a limitation in the type of BIPs used for activity models: all the BIPs used in a mixture must use the same equation. I hope to achive more freedom in the future.</div><div>&nbsp;</div><div><b>Transport properties</b></div><div>&nbsp;At the momment only viscosity and thermal conductivity calculations are supported. No BIPs are used except for liquid viscosity.</div><div>For liquid viscosity you can choose between Grunberg-Nissan, Teja-Rice, Andrade, McAllister 3 body and McAllister 4 body methods. But take into account that McAllister models are only for binary mixtures. BIPs are read from the corresponding file in the Resources/Interactions folder. If no BIPs are found in the file a less acurate value will be calculated.&nbsp;</div><div>Grunberg-Nissan and Teja-Rice use just one BIP, that is calculated according to a+b/T, being a and b the two first numbers following the CAS numbers of the substances in the row. Nevertheless you must fill with 0 the following four numbers.</div><div>Andrade and McAllister3 use two BIPs, calculated as a+b/T and c+d/T, being c and d the next numbers in the row.</div><div>McAllister4 uses a+b/T, c+d/T and e+f/T</div><div>For gas viscosity the Lucas method is used.</div><div>For liquid thermal conductivity the Li method.</div><div>For gas thermal conductivity that of Mason and Saxena</div><div><br></div><div><br></div><div><br></div></body></html>"));
end ExternalMix;

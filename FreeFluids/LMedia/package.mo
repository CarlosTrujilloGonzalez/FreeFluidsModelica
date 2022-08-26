within FreeFluids;

package LMedia "LMedia.mo by Carlos Trujillo
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
  partial package LMedium
    extends Modelica.Media.Interfaces.PartialPureSubstance(singleState = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, reference_T = 298.15, redeclare record FluidConstants = FreeFluids.MediaCommon.DataRecord);
    constant FluidConstants[nS] fluidConstants "Constant data for the fluid";
    import PhysPropCorr = FreeFluids.MediaCommon.Functions.PhysPropCorr;
    import SpecificEnthalpyCorr = FreeFluids.MediaCommon.Functions.SpecificEnthalpyCorr;
    import SpecificEnthalpyCorrInv = FreeFluids.MediaCommon.Functions.SpecificEnthalpyCorrInv;
    import SpecificEntropyCorr = FreeFluids.MediaCommon.Functions.SpecificEntropyCorr;
    import SpecificEntropyCorrInv = FreeFluids.MediaCommon.Functions.SpecificEntropyCorrInv;
    constant String refState = "None" "Enthalpy/entropy reference state. Alternatives: ASHRAE, IIR, NBP, User(value at reference_T as 0), otherwise no reference";
    //constant FreeFluids.MediaCommon.Types.ReferenceState refState=FreeFluids.MediaCommon.Types.ReferenceState.None  "Enthalpy/entropy reference state. Alternatives: ASHRAE, IIR, NBP, User(value at reference_T as 0), otherwise no reference";
    constant String inputChoice = "ph" "Allows to choose the input choise to use for the construction of a BaseProperties object. Alternative are: pT, ph, dT";
    constant Boolean highPressure = false;
    /*constant SpecificEnthalpy UserRefH = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T);
    constant SpecificEnthalpy ASHRAErefH = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15);
    constant SpecificEnthalpy NBPrefH = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb);
    constant SpecificEnthalpy IIRrefH = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5;
    constant SpecificEntropy UserRefS = SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T);
    constant SpecificEntropy ASHRAErefS = SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15);
    constant SpecificEntropy NBPrefS = SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb);
    constant SpecificEntropy IIRrefS = SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3;
    constant SpecificEnthalpy critical_h = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tc - 0.01);
    constant SpecificEnthalpy critical_h0 = SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, fluidConstants[1].Tc-0.01);*/

    //functions duplicated here because OpenModelica gets confused with the passed record
    //-----------------------------------------------------------------------------------
    function liqViscPcorLucas "Lucas liquid viscosity pressure correction."
      input Real T, p, Pref, rpVisc "temperature, pressure, reference pressure and viscosity at reference pressure";
      output Real visc;
    protected
      Real Tr, Tr2, Tr3, Tr4, A, C, D, dPr;
    algorithm
      Tr := T / fluidConstants[1].Tc;
      Tr2 := Tr * Tr;
      Tr3 := Tr2 * Tr;
      Tr4 := Tr3 * Tr;
      dPr := (p - Pref) / fluidConstants[1].criticalPressure;
      A := 0.9991 - 4.674e-4 / (1.0523 * Tr ^ (-0.03877) - 1.0513);
      C := (-0.07921) + 2.1616 * Tr - 13.404 * Tr2 + 44.1706 * Tr3 - 84.8291 * Tr4 + 96.1209 * Tr3 * Tr2 - 59.8127 * Tr3 * Tr3 + 15.6719 * Tr4 * Tr3;
      D := 0.3257 / (1.0039 - Tr ^ 2.573) ^ 0.2906 - 0.2086;
      visc := rpVisc * (1 + D * (dPr / 2.118) ^ A / (1 + C * fluidConstants[1].w * dPr));
    end liqViscPcorLucas;
  
  function liqThCondLatini
    input SI.Temperature T;
    output SI.ThermalConductivity lambda;
  protected
    Real Tr = T / fluidConstants[1].Tc, A, alpha, beta, gamma;
  algorithm
    if fluidConstants[1].family == 1 then
      A := 0.00350 "Alkane";
      alpha := 1.2;
      beta := 0.5;
      gamma := 0.167;
    elseif fluidConstants[1].family == 2 or fluidConstants[1].family == 3 then
      A := 0.0361 "Alkene, Alkyne";
      alpha := 1.2;
      beta := 1.0;
      gamma := 0.167;
    elseif fluidConstants[1].family == 4 then
      A := 0.031 "Cycloalkane";
      alpha := 1.2;
      beta := 1.0;
      gamma := 0.167;
    elseif fluidConstants[1].family == 5 or fluidConstants[1].family == 6 then
      A := 0.0346 "Aromatic, Water";
      alpha := 1.2;
      beta := 1.0;
      gamma := 0.167;
    elseif fluidConstants[1].family == 7 or fluidConstants[1].family == 8 then
      A := 0.00339 "Alcohol, Polyol";
      alpha := 1.2;
      beta := 0.5;
      gamma := 0.167;
    elseif fluidConstants[1].family == 10 then
      A := 0.0385 "Ether";
      alpha := 1.2;
      beta := 1.0;
      gamma := 0.167;
    elseif fluidConstants[1].family == 12 then
      A := 0.00383 "Ketone";
      alpha := 1.2;
      beta := 0.5;
      gamma := 0.167;
    elseif fluidConstants[1].family == 13 then
      A := 0.00319 "Acid";
      alpha := 1.2;
      beta := 0.5;
      gamma := 0.167;
    elseif fluidConstants[1].family == 14 then
      A := 0.0415 "Ester";
      alpha := 1.2;
      beta := 1.0;
      gamma := 0.167;
    elseif fluidConstants[1].family == 17 or fluidConstants[1].family == 18 then
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
    lambda := A * fluidConstants[1].Tb ^ alpha * (1 - Tr) ^ 0.38 / (fluidConstants[1].MW ^ beta * fluidConstants[1].Tc ^ gamma * Tr ^ (1 / 6));
  end liqThCondLatini;
  
  function liqSurfTensSastriRao
    input SI.Temperature T;
    output SI.SurfaceTension sigma;
  protected
    Real k, x, y, z, m;
  algorithm
    if fluidConstants[1].family == 7 or fluidConstants[1].family == 8 then
      k := 2.28 "alcohol, polyol";
      x := 0.25;
      y := 0.175;
      z := 0;
      m := 0.8;
    elseif fluidConstants[1].family == 13 then
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
    sigma := k * (fluidConstants[1].criticalPressure * 1e-5) ^ x * fluidConstants[1].Tb ^ y * fluidConstants[1].Tc ^ z * ((1 - T / fluidConstants[1].Tc) / (1 - fluidConstants[1].Tb / fluidConstants[1].Tc)) ^ m * 1e-3;
  end liqSurfTensSastriRao;
  
    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates and localInputChoice == "ph" then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates and (localInputChoice == "ph" or localInputChoice == "pT") then StateSelect.prefer else StateSelect.default), T(stateSelect = if preferredMediumStates and (localInputChoice == "pT" or localInputChoice == "dT") then StateSelect.prefer else StateSelect.default), d(stateSelect = if preferredMediumStates and localInputChoice == "dT" then StateSelect.prefer else StateSelect.default))
        constant String localInputChoice = inputChoice;
      equation
        MM = fluidConstants[1].MW / 1000 "in kg/mol";
        R_s = Modelica.Constants.R / MM;
        if localInputChoice == "ph" then
//state = setState_phX(p, h);
          T = temperature_ph(p, h);
          d = density_ph(p, h);
        elseif localInputChoice == "pT" then
//state = setState_pTX(p, T);
          d = density_pT(p, T);
          h = specificEnthalpy_pT(p, T);
        elseif localInputChoice == "dT" then
//state = setState_dTX(d, T);
          p = pressure_dT(d, T);
          h = specificEnthalpy_dT(d, T);
        else
          assert(false, "Invalid choice for BaseProperties inputChoice");
        end if;
//phase = state.phase;
        u = h - p / d;
        state.p = p;
        state.T = T;
        state.h = h;
        state.d = d;
    end BaseProperties;

    //Thermodynamic state definition and constructors
    //-----------------------------------------------

    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        AbsolutePressure p "Pressure in Pa";
        Temperature T "Kelvin temperature";
        Density d(displayUnit = "kg/m3") "density in kg/m3";
        SpecificEnthalpy h "specific enthalpy";
    end ThermodynamicState;

    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        AbsolutePressure Vp "saturated pressure at T";
        Density Dls "saturated liquid density at T";
        SpecificEnthalpy Href;

      algorithm
        state.p := p;
        state.T := T;
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, T) "saturated liquid density";
        Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
        state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "liquid density pressure correction";
        state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy at T.";
        if highPressure == true then
          state.h := if T <= 0.45 * fluidConstants[1].Tc then state.h + state.p / state.d - Vp / Dls else if T < 0.85 * fluidConstants[1].Tc then state.h + (0.85 * fluidConstants[1].Tc - T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else state.h "pressure correction for liquid enthalpy";
        end if;
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        state.h := state.h - Href "referenced enthalpy";
    end setState_pTX;

    redeclare function extends setState_dTX "Return ThermodynamicState record as function of T,d and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Density Dls;
        Density Dgs;
        SpecificEnthalpy Hls "enthalpy of the saturated liquid at given T";
        SpecificEnthalpy Href "reference enthalpy";
        AbsolutePressure Vp "vapor pressure at given T";
      algorithm
        state.T := T;
        state.d := d;
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, T) "saturated liquid density at T";
        Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy";
        Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
        state.p := if fluidConstants[1].lnuA > 0 then Vp + (exp((d - Dls) * fluidConstants[1].lnuA) - 1) * 1000 * Modelica.Constants.R * T * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) / (fluidConstants[1].lnuA * fluidConstants[1].MW) else (1 - Dls / d) / fluidConstants[1].IsothComp + Vp;
        state.h := if T <= 0.45 * fluidConstants[1].Tc then Hls + (state.p - Vp) / state.d else if T < 0.85 * fluidConstants[1].Tc then Hls + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p - Vp) / state.d else Hls;

        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        state.h := state.h - Href "referenced enthalpy";
    end setState_dTX;

    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,H and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEnthalpy Hls "saturated liquid enthalpy at given pressure";
        SpecificEnthalpy Hgs "saturated gas enthalpy at given pressure";
        Density Dls "saturated liquid density at given pressure";
        Density Dgs "saturated gas density at given pressure";
        Density Dgi "ideal gas density at p, state.T";
        SpecificEnthalpy Hraw "Input enthalpy without reference correction";
        AbsolutePressure Vp;
        SpecificEnthalpy Href;

      algorithm
        state.p := p;
        state.h := h;
        Tb:=FreeFluids.MediaCommon.Functions.PhysPropCorrInv(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, p, fluidConstants[1].VpLimS);
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Hraw := h + Href;
        Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific enthalpy at Tb";
        state.T := SpecificEnthalpyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Hraw, fluidConstants[1].lCpLimS - 0.01);
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
        Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
        state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
        if highPressure == true then
          Hraw := if state.T <= 0.45 * fluidConstants[1].Tc then Hraw - state.p / state.d + Vp / Dls else if state.T < 0.85 * fluidConstants[1].Tc then Hraw - (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else Hraw "pressure correction for liquid enthalpy";
          state.T := if fluidConstants[1].lTfromHsatCorr > 0 then PhysPropCorr(fluidConstants[1].lTfromHsatCorr, fluidConstants[1].lTfromHsatCoef, fluidConstants[1].MW, Hraw) else SpecificEnthalpyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Hraw, fluidConstants[1].lCpLimS - 0.01);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
          Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
        end if;
    end setState_phX;

    redeclare function extends setState_psX "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEntropy Sl;
        Density Dls;
        SpecificEntropy Sraw "Input entropy without reference correction";
        SpecificEntropy Sref;
        SpecificEnthalpy Href;
        AbsolutePressure Vp;

      algorithm
        state.p := p;
                Tb:=FreeFluids.MediaCommon.Functions.PhysPropCorrInv(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, p, fluidConstants[1].VpLimS);
        Sref := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Sraw := s + Sref;
        Sl := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific entropy at T";
        state.T := SpecificEntropyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Sraw, fluidConstants[1].Tc - 0.01);
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
        Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
        if highPressure == true then
          Sraw := if state.T <= 0.45 * fluidConstants[1].Tc then Sraw else Sraw + (state.T - 0.45 * fluidConstants[1].Tc) * (state.p - Vp) / (0.4 * fluidConstants[1].Tc * state.d * state.T) "pressure correction for liquid entropy";
          state.T := SpecificEntropyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Sraw, fluidConstants[1].Tc - 0.01);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
          Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
            state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
        end if;
        state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid enthalpy";
        if highPressure == true then
          state.h := if state.T <= 0.45 * fluidConstants[1].Tc then state.h + state.p / state.d - Vp / Dls else if state.T < 0.85 * fluidConstants[1].Tc then state.h + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else state.h "pressure correction for liquid enthalpy";
        end if;
      state.h := state.h - Href "referenced enthalpy";
    end setState_psX;

    //Special thermodynamic states constructors
    //------------------------------------------
    //Properties calculation from thermodynamic state
    //-----------------------------------------------

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

    redeclare function extends density_derp_T "Return density derivative w.r.t. pressure at const temperature. OK"
        extends Modelica.Icons.Function;

      algorithm
        ddpT := state.d * isothermalCompressibility(state);
    end density_derp_T;

    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure. OK"
        extends Modelica.Icons.Function;

      algorithm
        ddTp := -state.d * isobaricExpansionCoefficient(state);
    end density_derT_p;

    redeclare function extends density_derp_h "Return density derivative w.r.t. pressure at const specific enthalpy. Numerical derivative in two phases differs from HelmholtzMedia, but is similar to ExternalMedia"
        extends Modelica.Icons.Function;

      protected
        Real dvt "(dV/dT)P";
        Real dvp "(dV/dT)P";

      algorithm
        dvt := isobaricExpansionCoefficient(state) / state.d;
        dvp := -isothermalCompressibility(state) / state.d;
        ddph := -state.d * state.d * (dvp + (state.T * dvt * dvt - dvt / state.d) / specificHeatCapacityCp(state));
    end density_derp_h;

    redeclare function extends density_derh_p "Return density derivative w.r.t. specific enthalpy at constant pressure. OK"
        extends Modelica.Icons.Function;

      algorithm
        ddhp := -state.d * isobaricExpansionCoefficient(state) / specificHeatCapacityCp(state);
    end density_derh_p;

    redeclare function extends specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        h := state.h;
    end specificEnthalpy;

    redeclare function extends specificInternalEnergy "Return specific internal energy"
        extends Modelica.Icons.Function;

      algorithm
        u := state.h - state.p / state.d;
    end specificInternalEnergy;

    redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;

      protected
        SpecificEntropy Sl "entropy of the liquid phase";
        AbsolutePressure Vp;
        SpecificEntropy Sref;
        Temperature Tb;

      algorithm
        s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid entropy at T";
        if highPressure == true then
          Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
          s := if state.T <= 0.45 * fluidConstants[1].Tc then s else s - (state.T - 0.45 * fluidConstants[1].Tc) * (state.p - Vp) / (0.4 * fluidConstants[1].Tc * state.d * state.T) "pressure correction for liquid entropy";
        end if;
        Sref := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        s := s - Sref "referenced entropy";
    end specificEntropy;

    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
        extends Modelica.Icons.Function;

      algorithm
        g := state.h - state.T * specificEntropy(state);
    end specificGibbsEnergy;

    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;

      algorithm
        f := state.h - state.p / state.d - state.T * specificEntropy(state);
    end specificHelmholtzEnergy;

    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure. OK"
        extends Modelica.Icons.Function;

      protected
        Real ds;
        Real Tb;

      algorithm
        cp := PhysPropCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T);
        if highPressure == true then
          ds := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T);
          cp := if fluidConstants[1].family == 17 then cp * exp(-5.0 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5)) else cp * exp(-2.8 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5));
        end if;

    end specificHeatCapacityCp;

    redeclare function extends specificHeatCapacityCv "Return heat capacity at constant volume"
        extends Modelica.Icons.Function;

      protected
        Real beta, kappa;

      algorithm
        beta := isobaricExpansionCoefficient(state);
        kappa := isothermalCompressibility(state);
        cv := specificHeatCapacityCp(state) - state.T * beta * beta / (kappa * state.d);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificHeatCapacityCv;

    redeclare function extends isentropicExponent "Return isentropic exponent"
        extends Modelica.Icons.Function;

      protected
        Real cp, cv;

      algorithm
        cp := specificHeatCapacityCp(state);
        cv := specificHeatCapacityCv(state);
        gamma := cp / cv;
    end isentropicExponent;

    redeclare function extends isobaricExpansionCoefficient "Returns the approximate isobaric expansion coefficient beta=dV_dT/V"
        extends Modelica.Icons.Function;

      protected
        Real dd, d0, d1, d2, dv, vp1, vp2;

      algorithm
        d1 := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T - 0.1);
        d2 := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T + 0.1);
        d0 := (d1 + d2) / 2;
        beta := (d1 - d2) * d0 / (d1 * d2 * 0.2) "Coefficient along the saturation line, function only of T, not the real beta";
        vp1 := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T - 0.1);
        vp2 := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T + 0.1);
        dv := (vp2 - vp1) / 0.2;
        beta := beta + isothermalCompressibility(state) * dv "Beta at saturation, once corrected for pressure variation";
        if fluidConstants[1].lnuA > 0 then
          beta := beta * d0 * exp(-fluidConstants[1].lnuA * (state.d - d0)) / state.d + (1 - exp(-fluidConstants[1].lnuA * (state.d - d0))) / (fluidConstants[1].lnuA * state.d * state.T) "Beta at the given T and d. Postnikov et alt. 2016";
        end if;
    end isobaricExpansionCoefficient;

    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;

      protected
        Real rhoRed1 "reduced density - 1";

      algorithm
        if fluidConstants[1].lBulkModRCorr > 0 then
          kappa := fluidConstants[1].MW / (exp(fluidConstants[1].lBulkModRCoef[1] + fluidConstants[1].lBulkModRCoef[2] * state.d + fluidConstants[1].lBulkModRCoef[3] * state.d * state.d + fluidConstants[1].lBulkModRCoef[4] * state.d * state.d * state.d + fluidConstants[1].lBulkModRCoef[5] * state.d * state.d * state.d * state.d) * Modelica.Constants.R * state.T * state.d * 1000);
        elseif fluidConstants[1].lnuA > 0.0 then
          kappa := fluidConstants[1].MW / (exp(fluidConstants[1].lnuA * state.d + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T * state.d * 1000);
        elseif fluidConstants[1].Vc > 0.0 then
          rhoRed1 := state.d * 1000 * fluidConstants[1].Vc / fluidConstants[1].MW - 1;
          kappa := fluidConstants[1].MW / ((exp((-0.42704 * rhoRed1) + 2.089 * rhoRed1 * rhoRed1 + 0.42367 * rhoRed1 * rhoRed1 * rhoRed1) - 1) * Modelica.Constants.R * state.T * state.d * 1000);
        else
          kappa := fluidConstants[1].IsothComp;
        end if;
    end isothermalCompressibility;

    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;

      protected
        Real cp, kappa, beta;
        Real cv;

      algorithm
        kappa := isothermalCompressibility(state);
        beta := isobaricExpansionCoefficient(state);
        cp := specificHeatCapacityCp(state);
        cv := specificHeatCapacityCv(state);
        a := sqrt(max(0, cp / (cp * kappa * state.d - state.T * beta * beta)));
      annotation(
        Inline = true,
        smoothOrder = 2);
    end velocityOfSound;

    function jouleThomsonCoefficient
      input ThermodynamicState state;
      output Real J;
    algorithm
      J := (state.T * isobaricExpansionCoefficient(state) - 1) / (state.d * specificHeatCapacityCp(state));
    end jouleThomsonCoefficient;
  
    redeclare function extends dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;

      protected
        Real lpVisc;

      algorithm
        if fluidConstants[1].CAS == "7732-18-5" then
          eta := FreeFluids.MediaCommon.Functions.waterViscosity(state.T, state.d);
        elseif fluidConstants[1].lViscCorr > 0 then
          eta := if highPressure == false then PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, fluidConstants[1].MW, state.T) else liqViscPcorLucas(state.T, state.p, PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T), PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, fluidConstants[1].MW, state.T));
        else
          assert(false, "there is no correlation for liquid viscosity calculation");
        end if;
    end dynamicViscosity;

    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;

      algorithm
        lambda := if fluidConstants[1].lThCondCorr > 0 then PhysPropCorr(fluidConstants[1].lThCondCorr, fluidConstants[1].lThCondCoef, fluidConstants[1].MW, state.T) else liqThCondLatini(state.T);
        if highPressure == true then
          lambda := lambda * exp(state.d / PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) - 1);
        end if;
    end thermalConductivity;

    redeclare function extends molarMass
        "Return the molar mass of the medium"
    algorithm
        MM := fluidConstants[1].molarMass;
    end molarMass;
  
    function surfaceTension "Return surface tension"
      input ThermodynamicState state;
      output SurfaceTension sigma;
    algorithm
        sigma := if fluidConstants[1].lSurfTensCorr > 0 then PhysPropCorr(fluidConstants[1].lSurfTensCorr, fluidConstants[1].lSurfTensCoef, fluidConstants[1].MW, state.T) else liqSurfTensSastriRao(state.T);
    end surfaceTension;
  

    redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;

      algorithm
        state := ThermodynamicState(p = Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), d = Modelica.Media.Common.smoothStep(x, state_a.d, state_b.d, x_small), h = Modelica.Media.Common.smoothStep(x, state_a.h, state_b.h, x_small), T = Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
      annotation(
        Inline = true);
    end setSmoothState;

    function density_pT
      input AbsolutePressure p;
      input Temperature T;
      input FixedPhase phase = 0;
      output Density d;
    algorithm
      d := density_pT_state(p, T, setState_pTX(p, T));
      annotation(
        Inline = true);
        //Inline = true,
        //inverse(p = pressure_dT(d, T), T = temperature_pd(p, d)));
    end density_pT;

    function density_pT_state
      input AbsolutePressure p;
      input Temperature T;
      input ThermodynamicState state;
      output Density d;
    algorithm
      d := density(state);
      annotation(
        Inline = false,
        LateInline = true,
        inverse(p = pressure_dT_state(d = d, T = T, state = state)),
        derivative(noDerivative = state) = density_pT_der);
    end density_pT_state;

    function specificEnthalpy_pT
      input AbsolutePressure p;
      input Temperature T;
      input FixedPhase phase = 0;
      output SpecificEnthalpy h;
    algorithm
      h := specificEnthalpy_pT_state(p, T, setState_pTX(p, T));
      annotation(
        Inline = true,
        inverse(T = temperature_ph(p, h)));
    end specificEnthalpy_pT;

    function specificEnthalpy_pT_state
      input AbsolutePressure p;
      input Temperature T;
      input ThermodynamicState state;
      output SpecificEnthalpy h;
    algorithm
      h := specificEnthalpy(state);
      annotation(
        Inline = false,
        LateInline = true,
        inverse(T = temperature_ph_state(p = p, h = h, state = state)),
        derivative(noDerivative = state) = specificEnthalpy_pT_der);
    end specificEnthalpy_pT_state;

    function density_pT_der
      input AbsolutePressure p;
      input Temperature T;
      input ThermodynamicState state;
      input Real p_der "time derivative of pressure";
      input Real T_der "time derivative of temperature";
      output Real d_der "time derivative of density";
    algorithm
      d_der := p_der * density_derp_T(state) + T_der * density_derT_p(state);
      annotation(
        Inline = true);
    end density_pT_der;

    function specificEnthalpy_pT_der
      input AbsolutePressure p;
      input Temperature T;
      input ThermodynamicState state;
      input Real p_der "time derivative of pressure";
      input Real T_der "time derivative of temperature";
      output Real h_der "time derivative of specific Enthalpy";
    protected
      Real cp = specificHeatCapacityCp(state);
    algorithm
      h_der := (-p_der * ((state.T * isobaricExpansionCoefficient(state) - 1) / state.d)) + T_der * cp;
      annotation(
        Inline = true);
    end specificEnthalpy_pT_der;

    function density_ph
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input FixedPhase phase = 0;
      output Density d;
    algorithm
      d := density_ph_state(p, h, setState_phX(p, h));
      annotation(
        Inline = true);
        //Inline = true,
        //inverse(h = specificEnthalpy_pd(p, d)));
        
    end density_ph;

    function density_ph_state
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input ThermodynamicState state;
      output Density d;
    algorithm
      d := density(state);
      annotation(
        Inline = false,
        LateInline = true,
        derivative(noDerivative = state) = density_ph_der);
    end density_ph_state;

    function temperature_ph
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input FixedPhase phase = 0;
      output Temperature T;
    algorithm
      T := temperature_ph_state(p, h, setState_phX(p, h));
      annotation(
        Inline = true,
        inverse(h = specificEnthalpy_pT(p, T)));
    end temperature_ph;

    function temperature_ph_state
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input ThermodynamicState state;
      output Temperature T;
    algorithm
      T := temperature(state);
      annotation(
        Inline = false,
        LateInline = true,
        inverse(h = specificEnthalpy_pT_state(p = p, T = T, state = state)),
        derivative(noDerivative = state) = temperature_ph_der);
    end temperature_ph_state;

    function density_ph_der
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input ThermodynamicState state;
      input Real p_der "time derivative of pressure";
      input Real h_der "time derivative of specific enthalpy";
      output Real d_der "time derivative of density";
    algorithm
      d_der := p_der * density_derp_h(state) + h_der * density_derh_p(state);
      annotation(
        Inline = true);
    end density_ph_der;

    function temperature_ph_der
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      input ThermodynamicState state;
      input Real p_der "time derivative of pressure";
      input Real h_der "time derivative of specific enthalpy";
      output Real T_der "time derivative of density";
    algorithm
      T_der := p_der * jouleThomsonCoefficient(state) + h_der / specificHeatCapacityCp(state);
      annotation(
        Inline = true);
    end temperature_ph_der;

    function pressure_dT
      input Density d;
      input Temperature T;
      input FixedPhase phase = 0;
      output AbsolutePressure p;
    algorithm
      p := pressure_dT_state(d, T, setState_dTX(d, T));
      annotation(
        Inline = true);
        //Inline = true,
        //inverse(d = density_pT(p, T), T = temperature_pd(p, d)));
    end pressure_dT;

    function pressure_dT_state
      input Density d;
      input Temperature T;
      input ThermodynamicState state;
      output AbsolutePressure p;
    algorithm
      p := pressure(state);
      annotation(
        Inline = false,
        LateInline = true,
        inverse(d = density_pT_state(p = p, T = T, state = state)),
        derivative(noDerivative = state) = pressure_dT_der);
    end pressure_dT_state;

    function specificEnthalpy_dT
      input Density d;
      input Temperature T;
      input FixedPhase phase = 0;
      output SpecificEnthalpy h;
    algorithm
      h := specificEnthalpy_dT_state(d, T, setState_dTX(d, T));
      annotation(
        Inline = true);
    end specificEnthalpy_dT;

    function specificEnthalpy_dT_state
      input Density d;
      input Temperature T;
      input ThermodynamicState state;
      output SpecificEnthalpy h;
    algorithm
      h := specificEnthalpy(state);
      annotation(
        Inline = false,
        LateInline = true);
        //derivative(noDerivative = state) = specificEnthalpy_dT_der);
    end specificEnthalpy_dT_state;

    function pressure_dT_der
      input Density d;
      input Temperature T;
      input ThermodynamicState state;
      input Real d_der "time derivative of density";
      input Real T_der "time derivative of temperature";
      output Real p_der "time derivative of pressure";
    algorithm
      p_der := d_der * pressure_derd_T(state) + T_der * pressure_derT_d(state);
      annotation(
        Inline = true);
    end pressure_dT_der;

    function pressure_derd_T "pressure derivative w.r.t. density at constant T. OK"
      input ThermodynamicState state;
      output Real dpdT;
    algorithm
      dpdT := 1 / (state.d * isothermalCompressibility(state));
    end pressure_derd_T;

    function pressure_derT_d "pressure derivative w.r.t. T at constant density. OK"
      input ThermodynamicState state;
      output Real dpTd;
    algorithm
      dpTd := isobaricExpansionCoefficient(state) / isothermalCompressibility(state);
    end pressure_derT_d;

    function specifiEnthalpy_dT_der
      input Density d;
      input Temperature T;
      input ThermodynamicState state;
      input Real d_der "time derivative of density";
      input Real T_der "time derivative of temperature";
      output Real h_der "time derivative of specific enthalpy ";
    algorithm
      h_der := d_der * specificEnthalpy_derd_T(state) + T_der * specificEnthalpy_derT_d(state);
      annotation(
        Inline = true);
    end specifiEnthalpy_dT_der;

    function specificEnthalpy_derd_T "pressure derivative w.r.t. density at constant T. OK, but it seems clear that beta and kappa are bad calculated for the gas"
      input ThermodynamicState state;
      output Real dhdT;
    protected
      SaturationProperties sat;
    algorithm
      dhdT := (1 - isobaricExpansionCoefficient(state) * state.T) / (isothermalCompressibility(state) * state.d ^ 2);
    end specificEnthalpy_derd_T;

    function specificEnthalpy_derT_d "pressure derivative w.r.t. density at constant T. To be finished analyticaly for two phases"
      input ThermodynamicState state;
      output Real dhTd;
    protected
      Real beta;
      ThermodynamicState statePlus;
    algorithm
      beta := isobaricExpansionCoefficient(state);
      dhTd := specificHeatCapacityCp(state) + beta * (1 - beta * state.T) / state.d;
    end specificEnthalpy_derT_d;
  end LMedium;


  annotation(
    Documentation(info = "<html>
    <body>
    <p>The package TMedia is covering liquid, biphasic and gas phases, but this has increased the complexity, the computation time, and the convergence issues. In order to recover speed and compatibility, a reduced version, covering only the liquid phase has been implemented in this package.</p>
    </body>
    </html>"));
end LMedia;
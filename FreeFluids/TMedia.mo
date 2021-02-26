within FreeFluids;

package TMedia "TMedia.mo by Carlos Trujillo
  This file is part of the Free Fluids application
  Copyright (C) 2008-2020  Carlos Trujillo Gonzalez
    
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
  partial package TMedium
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(singleState = false, onePhase = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, reference_T = 298.15, redeclare record FluidConstants = FreeFluids.MediaCommon.DataRecord);
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
    //Auxiliary functions based in correlations
    //-----------------------------------------
  
    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;
  
      algorithm
        assert(T < fluidConstants[1].Tc, "The media can´t be used over Tc");
        p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
    end saturationPressure;
  
    function saturationPressureInv "Compute temperature from property value"
      extends FreeFluids.MediaCommon.Functions.CorrelationSolver(redeclare function f = PhysPropCorr);
    end saturationPressureInv;
  
    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;
  
      algorithm
        assert(p < fluidConstants[1].criticalPressure, "No saturated media possible over Pc", AssertionLevel.warning);
        T := if p >= fluidConstants[1].criticalPressure then fluidConstants[1].Tc - 0.01 else if fluidConstants[1].BtCorr > 0 then PhysPropCorr(fluidConstants[1].BtCorr, fluidConstants[1].BtCoef, fluidConstants[1].MW, p) else FreeFluids.MediaCommon.Functions.PhysPropCorrInv(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, p, fluidConstants[1].VpLimS - 0.01, fluidConstants[1].VpLimI+0.01);
    end saturationTemperature;
  
    redeclare function extends saturationTemperature_derp "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        dTp := (saturationTemperature(p) - saturationTemperature(0.999 * p)) / (0.001 * p);
    end saturationTemperature_derp;
  
    //BaseProperties model
    //--------------------
    /*redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), T(stateSelect =  StateSelect.default), d(stateSelect = StateSelect.default))
                constant String localInputChoice = inputChoice;
                Integer phase(min=0, max=2, start=1, fixed=false);
          
              algorithm
                sat := setSat_p(p);
                state := setState_phX(p, h);
                MM := fluidConstants[1].MW / 1000 "in kg/mol";
                R := Modelica.Constants.R / MM;
                T := state.T;
                d := state.d;
                u := h - p / d;
                phase:=state.phase;
                
                
            end BaseProperties;*/
  
    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates and localInputChoice == "ph" then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates and (localInputChoice == "ph" or localInputChoice == "pT") then StateSelect.prefer else StateSelect.default), T(stateSelect = if preferredMediumStates and (localInputChoice == "pT" or localInputChoice == "dT") then StateSelect.prefer else StateSelect.default), d(stateSelect = if preferredMediumStates and localInputChoice == "dT" then StateSelect.prefer else StateSelect.default))
        constant String localInputChoice = inputChoice;
        Integer phase(min = 0, max = 2, start = 1, fixed = false);
  
      algorithm
        MM := fluidConstants[1].MW / 1000 "in kg/mol";
        R := Modelica.Constants.R / MM;
        if localInputChoice == "ph" then
          state := setState_phX(p, h);
          T := state.T;
          d := state.d;
        elseif localInputChoice == "pT" then
          state := setState_pTX(p, T);
          d := state.d;
          h := state.h;
        elseif localInputChoice == "dT" then
          state := setState_dTX(d, T);
          p := state.p;
          h := state.h;
        else
          assert(false, "Invalid choice for BaseProperties inputChoice");
        end if;
        phase := state.phase;
        u := h - p / d;
        sat := setSat_p(p);
    end BaseProperties;
  
    //Thermodynamic state definition and constructors
    //-----------------------------------------------
  
    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        AbsolutePressure p "Pressure in Pa";
        Temperature T "Kelvin temperature";
        Density d(displayUnit = "kg/m3") "density in kg/m3";
        SpecificEnthalpy h "specific enthalpy";
        Fraction gf "gas fraction";
    end ThermodynamicState;
  
    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        AbsolutePressure Vp "saturated pressure at T";
        Temperature Tb "saturated temperature at p";
        Density Dls "saturated liquid density at T";
        Density DgsP "saturated gas density at p(Tb)";
        Density Dgi "ideal gas density at T and p";
        SpecificEnthalpy Href;
  
      algorithm
        assert(T < fluidConstants[1].Tc, "The media can´t be used over Tc.");
        state.p := p;
        state.T := T;
        Vp := saturationPressure(T);
        if p > Vp then
          state.phase := 1 "liquid";
          state.gf := 0;
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, T) "saturated liquid density";
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "liquid density pressure correction";
          state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy at T.";
          if highPressure == true then
            state.h := if T <= 0.45 * fluidConstants[1].Tc then state.h + state.p / state.d - Vp / Dls else if T < 0.85 * fluidConstants[1].Tc then state.h + (0.85 * fluidConstants[1].Tc - T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else state.h "pressure correction for liquid enthalpy";
          end if;
        elseif p < Vp then
          state.phase := 1 "gas";
          state.gf := 1.0;
          Tb := saturationTemperature(p);
          DgsP := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
          Dgi := p * fluidConstants[1].MW / (1000 * R * T);
          if p > 0.1 * Vp then
            state.d := Dgi + (DgsP - Dgi) / (0.9 * Vp) * (p - 0.1 * Vp) "ideal gas density at 0 pressure, saturated at Vp";
          else
            state.d := Dgi;
          end if;
          state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
          state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
          state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1+1.5e-7*p) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
        else
          assert(false, "A two phases state can't be constructed from p and T");
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
        Temperature Tb "boiling temperature at teste pressure";
        Density Dig "ideal gas density at T and p";
        Density Dsim "retrieved density at tested pressure";
        Density DsimN "retieved density at slightly higher pressure";
        Real Slope "density change with pressure";
        Integer i;
  
      algorithm
        assert(T < fluidConstants[1].Tc, "The media can´t be used over Tc");
        state.T := T;
        state.d := d;
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, T) "saturated liquid density at T";
        Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy";
        Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
        if d >= Dls then
          state.phase := 1 "liquid";
          state.gf := 0;
          state.p := if fluidConstants[1].lnuA > 0 then Vp + (exp((d - Dls) * fluidConstants[1].lnuA) - 1) * 1000 * Modelica.Constants.R * T * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) / (fluidConstants[1].lnuA * fluidConstants[1].MW) else (1 - Dls / d) / fluidConstants[1].IsothComp + Vp;
          state.h := if T <= 0.45 * fluidConstants[1].Tc then Hls + (state.p - Vp) / state.d else if T < 0.85 * fluidConstants[1].Tc then Hls + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p - Vp) / state.d else Hls;
        else
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, T);
          if d >= Dgs then
            state.p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
            state.h := if T <= 0.45 * fluidConstants[1].Tc then Hls + state.p / state.d - Vp / Dls else if T < 0.85 * fluidConstants[1].Tc then Hls + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else Hls "pressure correction for liquid enthalpy";
            if d > Dgs then
              state.phase := 2;
              state.gf := Dgs * (Dls - d) / (d * (Dls - Dgs));
              state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, T);
            elseif d == Dgs then
              state.phase := 1 "saturated gas";
              state.gf := 1;
              state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, T);
            end if;
          else
            state.phase := 1 "gas phase";
            state.gf := 1;
            Dig := 0.1 * Vp * fluidConstants[1].MW / (1000 * R * T);
            if d > Dig then
              state.p := 0.1 * Vp;
              for i in 1:5 loop
                Dig := state.p * fluidConstants[1].MW / (1000 * R * T);
                Tb := saturationTemperature(state.p);
                Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
                Dsim := Dig + (Dgs - Dig) / (0.9 * Vp) * (state.p - 0.1 * Vp);
                Dig := state.p * 1.001 * fluidConstants[1].MW / (1000 * R * T);
                Tb := saturationTemperature(state.p * 1.001);
                Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
                DsimN := Dig + (Dgs - Dig) / (0.9 * Vp) * (state.p * 1.001 - 0.1 * Vp);
                Slope := (DsimN - Dsim) / (0.001 * state.p);
                state.p := state.p + (d - Dsim) / Slope;
              end for;
            else
              state.p := d * 1000 * R * T / fluidConstants[1].MW "first aprox. of pressure as ideal gas";
            end if;
            Tb := saturationTemperature(state.p);
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
            state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
            state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1+1.5e-7*state.p) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
          end if;
        end if;
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
        SpecificEnthalpy H0b "ideal gas enthalpy at boiling p";
        SpecificEnthalpy H0t "ideal gas enthalpy at state.T";
  
      algorithm
        state.p := p;
        state.h := h;
        Tb := saturationTemperature(p);
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Hraw := h + Href;
        Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific enthalpy at Tb";
        if Hls >= Hraw then
          state.phase := 1 "liquid";
          state.gf := 0;
          state.T := if fluidConstants[1].lTfromHsatCorr > 0 then PhysPropCorr(fluidConstants[1].lTfromHsatCorr, fluidConstants[1].lTfromHsatCoef, fluidConstants[1].MW, Hraw) else SpecificEnthalpyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Hraw, fluidConstants[1].lCpLimS - 0.01);
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
        else
          Hgs := Hls + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb);
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
          if Hraw <= Hgs then
            state.phase := 2 "biphasic";
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, Tb);
            state.T := Tb;
            state.gf := (Hraw - Hls) / (Hgs - Hls);
            state.d := Dls * Dgs / (state.gf * Dls + (1 - state.gf) * Dgs);
          else
            state.phase := 1 "gas";
            state.gf := 1.0;
            H0b := SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb);
            H0t := (Hraw - Hgs) / (1 + 1.5e-7*p) + H0b;
            state.T := SpecificEnthalpyCorrInv(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, H0t, fluidConstants[1].Tc - 0.01);
            Dgi := p * fluidConstants[1].MW / (1000 * R * state.T);
            Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
            if p > 0.1 * Vp then
              state.d := Dgi + (Dgs - Dgi) * (p - 0.1 * Vp) / (0.9 * Vp) "ideal gas density at 0 pressure, saturated at Vp";
            else
              state.d := Dgi;
            end if;
          end if;
        end if;
  /*assert(Hgs >= Hraw, "The media can´t work with non-saturated gases");*/
    end setState_phX;
  
    redeclare function extends setState_psX "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        Temperature Tb;
        SpecificEntropy Sl;
        SpecificEntropy Sg;
        SpecificEntropy Sexc;
        SpecificEntropy Sig;
        Density Dls;
        Density Dgs;
        Density Dgi;
        SpecificEntropy Sraw "Input entropy without reference correction";
        SpecificEntropy Sref;
        SpecificEnthalpy Href;
        AbsolutePressure Vp;
  
      algorithm
        state.p := p;
        if p > fluidConstants[1].criticalPressure then
          Tb := fluidConstants[1].Tc;
        else
          Tb := saturationTemperature(p);
        end if;
        Sref := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        Sraw := s + Sref;
        Sl := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific entropy at T";    
        if Sl >= Sraw then
          state.phase := 1 "liquid";
          state.gf := 0;
          state.T := SpecificEntropyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Sraw, fluidConstants[1].Tc-0.01);
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
        else
          Sg := Sl + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) / Tb;
          if Sg >= Sraw then
            state.phase := 2 "biphasic";
            state.T := Tb;
            state.gf := (Sraw - Sl) / (Sg - Sl);
            Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, Tb);
            state.d := Dls * Dgs / (state.gf * Dls + (1 - state.gf) * Dgs);
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid enthalpy";
            state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) / Tb "unreferenced mix enthalpy";
          else
            state.phase := 1 "gas";
            state.gf := 1.0;
            Sexc := Sraw - Sg;
            Sig := SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb);
            state.T := SpecificEntropyCorrInv(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Sig + Sexc, 1500.0);
            Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
            Dgi := p * fluidConstants[1].MW / (1000 * R * state.T);
            Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
            if p > 0.1 * Vp then
              state.d := Dgi + (Dgs - Dgi) / (0.9 * Vp) * (p - 0.1 * Vp) "ideal gas density at 0 pressure, saturated at Vp";
            else
              state.d := Dgi;
            end if;
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
            state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
            state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1+1.5e-7*p) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
          end if;
        end if;
  //state.d := Dgi + (Dgs - Dgi) / Vp * p "ideal gas density at 0 pressure, saturated at Vp";
        state.h := state.h - Href "referenced enthalpy";
    end setState_psX;
  
    //Special thermodynamic states constructors
    //------------------------------------------
  
    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEnthalpy Href;
  
      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) - Href "referenced enthalpy";
        state.d := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T);
        state.gf := 0.0 "liquid";
        state.phase := 1;
    end setBubbleState;
  
    redeclare function extends setDewState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the dew point"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEnthalpy Href;
  
      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) "unreferenced saturated liquid enthalpy";
        state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat) - Href "referenced saturated gas enthalpy";
        state.d := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, sat.Tsat);
        state.gf := 1.0 "gas";
        state.phase := 1;
    end setDewState;
  
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
  
    redeclare function extends density_derp_T "Return density derivative w.r.t. pressure at const temperature"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.gf < 0.000001 or state.gf > 0.999999 then
          ddpT := state.d * isothermalCompressibility(state);
        else
          ddpT := 0;
        end if;
    end density_derp_T;
  
    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.gf < 0.000001 or state.gf > 0.999999 then
          ddTp := -state.d * isobaricExpansionCoefficient(state);
        else
          ddTp := 0;
        end if;
    end density_derT_p;
  
    redeclare function extends density_derp_h "Return density derivative w.r.t. pressure at const specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        Real dvt "(dV/dT)P";
        Real dvp "(dV/dT)P";
  
      algorithm
        if state.gf < 0.000001 or state.gf > 0.999999 then
          dvt := isobaricExpansionCoefficient(state) / state.d;
          dvp := -isothermalCompressibility(state) / state.d;
          ddph := -state.d * state.d * (dvp + (state.T * dvt * dvt - dvt / state.d) / specificHeatCapacityCp(state));
        else
          ddph := (density(setState_phX(1.02 * state.p, state.h)) - state.d) / (0.02 * state.p);
        end if;
    end density_derp_h;
  
    redeclare function extends density_derh_p "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.gf < 0.000001 or state.gf > 0.999999 then
          ddhp := -state.d * isobaricExpansionCoefficient(state) / specificHeatCapacityCp(state);
        else
          ddhp := -state.d * state.d * (saturationTemperature(state.p * 1.05) - state.T) / (state.p * 0.05 * state.T);
        end if;
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
        if state.gf <= 0 then//liquid
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid entropy at T";
          if highPressure == true then
            Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
            s := if state.T <= 0.45 * fluidConstants[1].Tc then s else s - (state.T - 0.45 * fluidConstants[1].Tc) * (state.p - Vp) / (0.4 * fluidConstants[1].Tc * state.d * state.T) "pressure correction for liquid entropy";
          end if;
        elseif state.gf < 1 then//biphasic
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid entropy at T";
          s := s + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, state.T) / state.T "mix liquid-vapor entropy";
        else//gas
          Tb := saturationTemperature(state.p);
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid entropy at p(Tb)";
          s := s + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) / Tb "saturated gas entropy at p(Tb)";
          s := s + SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T) - SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb) "gas entropy at T";
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
  
    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
  
      protected
        Real ds;
        Real Tb;
      algorithm
        if state.gf == 0.0 then
          cp := PhysPropCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T);
          if highPressure == true then
            ds := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T);
            cp := if fluidConstants[1].family == 17 then cp * exp(-5.0 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5)) else cp * exp(-2.8 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5));
          end if;
        elseif state.gf == 1.0 then
          cp := PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
          Tb:=saturationTemperature(state.p);
          //cp:=cp*(1+1.5e-7*state.p);
          cp:=cp+cp*3.5e-7*state.p^1.001*(1-(state.T-Tb)/(fluidConstants[1].Tc-Tb))^1.5;
        else
          cp := 0;
        end if;
    end specificHeatCapacityCp;
  
    redeclare function extends specificHeatCapacityCv "Return heat capacity at constant volume"
        extends Modelica.Icons.Function;
  
      protected
        Real beta, kappa;
  
      algorithm
        if state.gf == 0.0 then
          beta := isobaricExpansionCoefficient(state);
          kappa := isothermalCompressibility(state);
          cv := specificHeatCapacityCp(state) - state.T * beta * beta / (kappa * state.d);
        elseif state.gf == 1.0 then
          //cv := PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T) - R * 1000 / fluidConstants[1].MW;
          cv:=specificHeatCapacityCp(state)- R * 1000 / fluidConstants[1].MW;
        else
          cv := 0;
        end if;
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
  
    redeclare function extends isobaricExpansionCoefficient "Returns the approximate isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
  
      protected
        Real dd, d0, d1, d2, dv, vp1, vp2;
  
      algorithm
        if state.gf == 0 then
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
        elseif state.gf == 1.0 then
          beta := 1 / state.T "ideal gas expansion coefficient";
        else
          beta := 0;
        end if;
  //at saturation
    end isobaricExpansionCoefficient;
  
    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
  
      protected
        Real rhoRed1 "reduced density - 1";
  
      algorithm
        if state.gf == 0.0 then
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
        elseif state.gf == 1.0 then
          kappa := 1.0 / state.p "ideal gas compressibility";
        else
          kappa := 0.0;
        end if;
    end isothermalCompressibility;
  
    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
  
      protected
        Real cp, kappa, beta;
        Real cv;
  
      algorithm
        kappa := isothermalCompressibility(state);
        if state.gf == 0.0 then
          beta := isobaricExpansionCoefficient(state);
          cp := specificHeatCapacityCp(state);
          cv := specificHeatCapacityCv(state);
          a := sqrt(max(0, cp / (cp * kappa * state.d - state.T * beta * beta)));
        elseif state.gf == 1.0 then
          cp := PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
          cv:= cp - R * 1000 / fluidConstants[1].MW;
          a := sqrt(max(0, cp / (cv * kappa * state.d)));
        else
          a := 0;
        end if;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end velocityOfSound;
  
    redeclare function extends dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
  
      protected
        Real lpVisc;
  
      algorithm
        if state.gf == 0.0 then
          if fluidConstants[1].CAS == "7732-18-5" then
            eta := FreeFluids.MediaCommon.Functions.waterViscosity(state.T, state.d);
          elseif fluidConstants[1].lViscCorr > 0 then
            eta := if highPressure == false then PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, fluidConstants[1].MW, state.T) else FreeFluids.MediaCommon.Functions.liqViscPcorLucas(fluidConstants[1], state.T, state.p, PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T), PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, fluidConstants[1].MW, state.T));
          else
            assert(false, "there is no correlation for liquid viscosity calculation");
          end if;
        elseif state.phase == 1.0 then
          if fluidConstants[1].CAS == "7732-18-5" then
            eta := FreeFluids.MediaCommon.Functions.waterViscosity(state.T, state.d);
          else
            lpVisc := if fluidConstants[1].gViscCorr > 0 then PhysPropCorr(fluidConstants[1].gViscCorr, fluidConstants[1].gViscCoef, fluidConstants[1].MW, state.T) else FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(fluidConstants[1], state.T);
            eta := lpVisc "no pressure correction seems necessary for gas phase";
          end if;
        else
          eta := 0 "no viscosity calculation is done for two phases situation";
          assert(false, "no viscosity calculation is done for two phases situation", AssertionLevel.warning);
        end if;
  //eta:= if highPressure==false then lpVisc else  FreeFluids.MediaCommon.Functions.gasViscPcorLucas(fluidConstants[1], state.T, state.p, lpVisc);
    end dynamicViscosity;
  
    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.gf == 0.0 then
          lambda := if fluidConstants[1].lThCondCorr > 0 then PhysPropCorr(fluidConstants[1].lThCondCorr, fluidConstants[1].lThCondCoef, fluidConstants[1].MW, state.T) else FreeFluids.MediaCommon.Functions.liqThCondLatini(fluidConstants[1], state.T);
          if highPressure == true then
            lambda := lambda * exp(state.d / PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) - 1);
          end if;
        elseif state.gf == 1.0 then
          lambda := if fluidConstants[1].gThCondCorr > 0 then PhysPropCorr(fluidConstants[1].gThCondCorr, fluidConstants[1].gThCondCoef, fluidConstants[1].MW, state.T) else FreeFluids.MediaCommon.Functions.gasThCondLowPressureChung(fluidConstants[1], specificHeatCapacityCp(state), dynamicViscosity(state), state.T);
        else
          lambda := 0 "no thermal conductivity calculation is done in two phases situation";
          assert(false, "no thermal conductivity calculation is done for two phases situation", AssertionLevel.warning);
        end if;
    end thermalConductivity;
  
    redeclare function extends surfaceTension "Return surface tension"
        extends Modelica.Icons.Function;
  
      algorithm
        sigma := if fluidConstants[1].lSurfTensCorr > 0 then PhysPropCorr(fluidConstants[1].lSurfTensCorr, fluidConstants[1].lSurfTensCoef, fluidConstants[1].MW, sat.Tsat) else FreeFluids.MediaCommon.Functions.liqSurfTensSastriRao(fluidConstants[1], sat.Tsat);
    end surfaceTension;
  
    //Saturation properties
    //-------------------------
  
    redeclare function extends bubbleDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;
  
      algorithm
        dl := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, sat.Tsat);
    end bubbleDensity;
  
    redeclare function extends dBubbleDensity_dPressure "Return bubble point density derivative"
        extends Modelica.Icons.Function;
  
      algorithm
        ddldp := (PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, sat.Tsat) - PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, sat.Tsat - 0.01)) / (sat.psat - saturationPressure(sat.Tsat - 0.01));
    end dBubbleDensity_dPressure;
  
    redeclare function extends dewDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;
  
      algorithm
        dv := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, sat.Tsat);
    end dewDensity;
  
    redeclare function extends dDewDensity_dPressure "Return dew point density derivative"
        extends Modelica.Icons.Function;
  
      algorithm
        ddvdp := (PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, sat.Tsat) - PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, sat.Tsat - 0.01)) / (sat.psat - saturationPressure(sat.Tsat - 0.01));
    end dDewDensity_dPressure;
  
    redeclare function extends bubbleEnthalpy "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEnthalpy Href;
  
      algorithm
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        hl := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) - Href "referenced enthalpy";
    end bubbleEnthalpy;
  
    redeclare function extends dBubbleEnthalpy_dPressure "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;
  
      algorithm
        dhldp := (SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) - SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat - 0.01)) / (sat.psat - saturationPressure(sat.Tsat - 0.01));
    end dBubbleEnthalpy_dPressure;
  
    redeclare function extends dewEnthalpy "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEnthalpy Href;
  
      algorithm
        Href := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        hv := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat) - Href "referenced enthalpy";
    end dewEnthalpy;
  
    redeclare function extends dDewEnthalpy_dPressure "Return dew point specific enthalpy derivative"
        extends Modelica.Icons.Function;
  
      algorithm
  //dhvdp := (SpecificEnthalpyCorr(sat.Tsat) + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat) - SpecificEnthalpyCorr(sat.Tsat-0.01)-PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat-0.01))/(sat.psat-saturationPressure(sat.Tsat-0.01))"bad precision";
        dhvdp := ((-PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat)) + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat - 0.01)) / (sat.psat - saturationPressure(sat.Tsat - 0.01));
  //dhvdp:= (-PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat) + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, saturationTemperature(0.9*sat.psat)))/(0.1*sat.psat)"needs a high pressure range. As correlations are T functions better not to use.";
    end dDewEnthalpy_dPressure;
  
    redeclare function extends bubbleEntropy "Return bubble point specific entropy"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEntropy Sref;
  
      algorithm
        Sref := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        sl := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) - Sref "referenced entropy";
    end bubbleEntropy;
  
    redeclare function extends dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;
  
      protected
        SpecificEntropy Sref;
  
      algorithm
        Sref := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
        sv := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, sat.Tsat) + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat) / sat.Tsat - Sref "referenced entropy";
    end dewEntropy;
  
    redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
  
      algorithm
        state := ThermodynamicState(p = Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), s = Modelica.Media.Common.smoothStep(x, state_a.s, state_b.s, x_small), h = Modelica.Media.Common.smoothStep(x, state_a.h, state_b.h, x_small), T = Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small), ld = Modelica.Media.Common.smoothStep(x, state_a.ld, state_b.ld, x_small), gd = Modelica.Media.Common.smoothStep(x, state_a.gd, state_b.gd, x_small), gf = Modelica.Media.Common.smoothStep(x, state_a.gf, state_b.gf, x_small), phase = 0);
      annotation(
        Inline = true);
    end setSmoothState;
  end TMedium;

  //*****MEDIUMS OF TEMPERATURE ONLY DEPENDENT PROPERTIES*****
  //==========================================================

  package Fluids
    package Acetone
      extends FreeFluids.TMedia.TMedium(final mediumName = "Acetone", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Acetone}, refState = "IIR");
    end Acetone;

    package Ammonia
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ammonia", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ammonia}, reference_T = 273.15);
    end Ammonia;

    package Butane_n
      extends TMedium(final mediumName = "Butane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Butane_n}, refState = "IIR", reference_T = 273.15);
    end Butane_n;

    package Butanol_n
      extends FreeFluids.TMedia.TMedium(final mediumName = "n-Butanol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Butanol_n});
    end Butanol_n;

    package CO2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Carbon dioxide", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.CO2});
    end CO2;

    package DecanoicAcid
      extends FreeFluids.TMedia.TMedium(final mediumName = "DecanoicAcid", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.DecanoicAcid}, refState = "IIR", reference_T = 273.15);
    end DecanoicAcid;

    package Dichlorodifluormethane
      extends FreeFluids.TMedia.TMedium(final mediumName = "Dichlorodifluormethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Dichlorodifluormethane}, reference_T = 273.15, refState = "ASHRAE");
    end Dichlorodifluormethane;

    package EG
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ethylene glycol", final singleState = false, p_default = 3.0e5, T_default = 423.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.EG});
    end EG;

    package Ethane
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethane}, refState = "IIR", reference_T = 273.15);
    end Ethane;

    package Ethanol
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ethanol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethanol}, reference_T = 273.15);
    end Ethanol;

    package Heptane_n
      extends FreeFluids.TMedia.TMedium(final mediumName = "n-Heptane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Heptane_n}, refState = "IIR", reference_T = 273.15);
    end Heptane_n;

    package Hexane_n
      extends FreeFluids.TMedia.TMedium(final mediumName = "n-Hexane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Hexane_n}, refState = "IIR", reference_T = 273.15, p_default = 5.0e5, T_default = 293.15);
    end Hexane_n;

    package Isobutane
      extends FreeFluids.TMedia.TMedium(final mediumName = "Isobutane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Isobutane}, reference_T = 273.15);
    end Isobutane;

    package MarlothermSH
      extends FreeFluids.TMedia.TMedium(final mediumName = "Marlotherm SH", final singleState = false, p_default = 3.0e5, T_default = 473.15, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.MarlothermSH});
    end MarlothermSH;

    package Methane
      extends TMedium(final mediumName = "Methane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Methane}, refState = "IIR", reference_T = 273.15);
    end Methane;

    package N2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Nitrogen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.N2});
    end N2;

    package O2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Oxygen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.O2});
    end O2;

    package Pentane_n
      extends TMedium(final mediumName = "Pentane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Pentane_n}, refState = "IIR", reference_T = 273.15);
    end Pentane_n;

    package Propane
      extends TMedium(final mediumName = "Propane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Propane}, refState = "NBP");
    end Propane;

    package R134A
      extends TMedium(final mediumName = "R134A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R134A}, refState = "ASHRAE", reference_T = 273.15);
    end R134A;

    package R410A
      extends TMedium(final mediumName = "R410A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R410A}, refState = "ASHRAE", reference_T = 273.15);
    end R410A;

    package ShellS2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Shell S2", final singleState = false, p_default = 3.0e5, T_default = 473.15, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.ShellS2});
    end ShellS2;

    package Toluene
      extends FreeFluids.TMedia.TMedium(final mediumName = "Toluene", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Toluene}, reference_T = 273.15, p_default = 5.0e5, T_default = 393.15);
    end Toluene;

    package Water
      extends FreeFluids.TMedia.TMedium(final mediumName = "Water", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Water}, reference_T = 273.16, refState = "User", p_default = 5.0e5, T_default = 393.15);
    end Water;

    package Xylene_m
      extends FreeFluids.TMedia.TMedium(final mediumName = "Xylene_m", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Xylene_m}, refState = "NBP", reference_T = 273.15);
    end Xylene_m;
  end Fluids;

  package Tests
    partial model FluidTestingA
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.AbsolutePressure p = 1.0e5;
      parameter Medium.Temperature initialT = 273.19;
      parameter Medium.Temperature finalT = 372.15;
      parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
      Medium.Temperature T(start = initialT) "We will ramp the temperature";
      Medium.MolarMass MM;
      Medium.ThermodynamicState StateP "original state from p and T";
      Medium.ThermodynamicState StateD "state reproduced from d,T";
      Medium.ThermodynamicState StateH "state reproduced from p,h";
      Medium.ThermodynamicState StateS "state reproduced from p,s";
      Medium.ThermodynamicState StateHalfH "state half way between Bub and Dew";
      //Properties of StateP
      Real Tb;
      Real H;
      Real D;
      Real S;
      Real U "internal energy";
      Real G "Gibbs energy";
      Real A "Helmholtz energy";
      Real Cp;
      Real Cv;
      Real Gamma;
      Real SS "speed of sound";
      Real Beta "isobaric expansion coefficient";
      Real Kappa "isothermal compressibility";
      //Real DerD_p_T;
      //Real DerD_T_p;
      Real DerD_p_h;
      Real DerD_h_p;
      Real Mu;
      Real Th;
      //Properties of StateDlow
      Medium.ThermodynamicState StateDlow "state at fraction of original density";
      Real DlowGasFract "gas fraction of StateDlow";
      Real DlowD, DlowH;
      Real DlowS;
      Real DlowU;
      Real DlowA;
      Real DlowCv;
      Real DlowDerD_p_h;
      Real DlowDerD_h_p;
      //Saturation properties
      Medium.SaturationProperties sat "saturation point at given T";
      Medium.AbsolutePressure Vp;
      Real DerTb_p;
      Medium.ThermodynamicState StateBub "bubble state at sat";
      Medium.ThermodynamicState StateDew "dew state at sat";
      Real BubD;
      Real BubH "bubble properties";
      Real BubS;
      Real BubDerD_p;
      Real BubDerH_p;
      Real DewD;
      Real DewH "dew properties";
      Real DewS;
      Real DewDerD_p;
      Real DewDerH_p;
      Real Hv "vaporization enthalpy";
      Real DewMu;
      Real Sigma;
      Medium.Temperature Tsat "temperature recovered from saturation pressure";
      //BaseProperties
      Medium.BaseProperties BaseProp;
    algorithm
    //Construction of StateP and calculation of properties
      for i in 1:1 loop
        StateP := Medium.setState_pTX(p, T);
        MM := Medium.molarMass(StateP);
        Tb := Medium.saturationTemperature(p);
        H := Medium.specificEnthalpy(StateP);
        D := Medium.density(StateP);
        S := Medium.specificEntropy(StateP);
        U := Medium.specificInternalEnergy(StateP);
        G := Medium.specificGibbsEnergy(StateP);
        A := Medium.specificHelmholtzEnergy(StateP);
        Cp := Medium.specificHeatCapacityCp(StateP);
        Cv := Medium.specificHeatCapacityCv(StateP);
        Gamma := Cp / Cv;
        SS := Medium.velocityOfSound(StateP);
        Beta := Medium.isobaricExpansionCoefficient(StateP);
        Kappa := Medium.isothermalCompressibility(StateP);
        DerD_p_h := Medium.density_derp_h(StateP);
        DerD_h_p := Medium.density_derh_p(StateP);
        Mu := Medium.dynamicViscosity(StateP);
        Th := Medium.thermalConductivity(StateP);
        StateD := Medium.setState_dTX(D, T);
        StateH := Medium.setState_phX(p, H);
        StateDlow := Medium.setState_dTX(D * fract, T, fill(0, 0));
        DlowGasFract := Medium.vapourQuality(StateDlow);
        DlowD := Medium.density(StateDlow);
        DlowH := Medium.specificEnthalpy(StateDlow);
        DlowU := Medium.specificInternalEnergy(StateDlow);
        DlowS := Medium.specificEntropy(StateDlow);
        DlowA := Medium.specificHelmholtzEnergy(StateDlow);
        DlowCv := Medium.specificHeatCapacityCv(StateDlow);
        DlowDerD_p_h := Medium.density_derp_h(StateDlow);
        DlowDerD_h_p := Medium.density_derh_p(StateDlow);
        sat := Medium.setSat_T(T);
        Vp := Medium.saturationPressure_sat(sat);
        DerTb_p := Medium.saturationTemperature_derp_sat(sat);
        StateBub := Medium.setBubbleState(sat);
        StateDew := Medium.setDewState(sat);
        BubD := Medium.bubbleDensity(sat);
        BubDerD_p := Medium.dBubbleDensity_dPressure(sat);
        BubH := Medium.bubbleEnthalpy(sat);
        BubDerH_p := Medium.dBubbleEnthalpy_dPressure(sat);
        BubS := Medium.bubbleEntropy(sat);
        DewD := Medium.dewDensity(sat);
        DewDerD_p := Medium.dDewDensity_dPressure(sat);
        DewH := Medium.dewEnthalpy(sat);
        DewDerH_p := Medium.dDewEnthalpy_dPressure(sat);
        Hv := DewH - BubH;
        DewS := Medium.dewEntropy(sat);
        DewMu := Medium.dynamicViscosity(StateDew);
        Tsat := Medium.saturationTemperature(sat.psat);
        Sigma := Medium.surfaceTension(sat);
      //DerD_p_T := Medium.density_derp_T(StateP);
      //DerD_T_p := Medium.density_derT_p(StateP);
        StateS := Medium.setState_psX(p, S);
        StateHalfH := Medium.setState_phX(Vp, (BubH + DewH) / 2, fill(0, 0));    
      end for;
    equation
    //Construction of BaseProperties
      BaseProp.p = p;
      BaseProp.T = T;
      der(T) = finalT - initialT;
      annotation(
        Documentation(info = "<html>
        <body>
        <p>This test maintains a constant pressure, and ramps temperature between the selected values. At each temperature, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p>
        </body>
        </html>"));
    end FluidTestingA;

    model TestA1A
      extends FluidTestingA(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), p = 20e5, initialT = 0.1 + 273.15, fract = 0.1, finalT = 380 + 273.15);
    end TestA1A;

    model TestA1B
      extends TestA1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_pT);
    end TestA1B;

    model TestA1C
      extends TestA1A(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.WaterRef(refState = 4, reference_T = 273.16, reference_p = 101325.0, inputChoice = "pT", thermoModel = 3));
    end TestA1C;

    model TestA2A
      extends FluidTestingA(redeclare replaceable package Medium = TMedia.Fluids.R134A(highPressure = true, refState = "IIR", reference_T = 100.0, inputChoice = "pT"), p = 20.0e5, initialT = (-50.0) + 273.15, finalT = 110.0 + 273.15, fract = 0.5);
    end TestA2A;

    model TestA2B
      extends TestA2A(redeclare package Medium = Modelica.Media.R134a.R134a_ph);
    end TestA2B;

    model TestA2C
      extends TestA2A(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.R134A(thermoModel = 3, refState = 2, inputChoice = "pT"));
    end TestA2C;
  
    partial model FluidTestingB
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.Temperature T = 273.2;
      parameter Medium.AbsolutePressure initialP = 1.0e5;
      parameter Medium.AbsolutePressure finalP = 1.0e7;
      parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
      Medium.AbsolutePressure p(start = initialP) "We will ramp the temperature";
      Medium.ThermodynamicState StateP "original state from p and T";
      Medium.ThermodynamicState StateH "state reproduced from p,h";
      Medium.ThermodynamicState StateS "state reproduced from p,s";
      Medium.ThermodynamicState StateD "state reproduced from d,T";
      Medium.ThermodynamicState StateDlow "state at fraction of original density";
      Medium.ThermodynamicState StateBub "bubble state at sat";
      Medium.ThermodynamicState StateDew "dew state at sat";
      //Properties of StateP
      Real D, H, Mu, Th;
      Real S;
      Real Beta "isobaric expansion coefficient";
      Real Kappa "isothermal compressibility";
      Real Cp;
      Real Cv;
      Real Gamma;
      Real SS "speed of sound";
      //Properties of StateDlow
      Real DlowGasFract "gas fraction of StateDlow";
      Real DlowD, DlowH;
      Real DlowS;
      //Saturation properties
      Medium.SaturationProperties sat "saturation point at given T";
      Real BubD, BubH "bubble properties";
      Real BubS;
      Real DewD, DewH "bubble properties";
      Real DewS;
      Real Hv "vaporization enthalpy";
      Real Sigma;
      Medium.Temperature Tsat "temperature recovered from saturation pressure";
      //BaseProperties
      Medium.BaseProperties BaseProp;
    algorithm
//Construction of StateP and calculation of properties
      StateP := Medium.setState_pTX(p, T, fill(0, 0));
      H := Medium.specificEnthalpy(StateP);
      D := Medium.density(StateP);
      S := Medium.specificEntropy(StateP);
      Cp := Medium.specificHeatCapacityCp(StateP);
      Cv := Medium.specificHeatCapacityCv(StateP);
      Gamma := Cp / Cv;
      SS := Medium.velocityOfSound(StateP);
      Beta := Medium.isobaricExpansionCoefficient(StateP);
      Kappa := Medium.isothermalCompressibility(StateP);
      Mu := Medium.dynamicViscosity(StateP);
      Th := Medium.thermalConductivity(StateP);
//Reconstruction of states from properties
      StateD := Medium.setState_dTX(D, T, fill(0, 0));
      StateH := Medium.setState_phX(p, H, fill(0, 0));
      StateS := Medium.setState_psX(p, S, fill(0, 0));
//Construction of state Dlow and get properties
      StateDlow := Medium.setState_dTX(D * fract, T, fill(0, 0));
      DlowGasFract := Medium.vapourQuality(StateDlow);
      DlowD := Medium.density(StateDlow);
      DlowH := Medium.specificEnthalpy(StateDlow);
      DlowS := Medium.specificEntropy(StateDlow);
//Calculations at saturation
      sat := Medium.setSat_T(T);
      StateBub := Medium.setBubbleState(sat);
      StateDew := Medium.setDewState(sat);
      BubD := Medium.bubbleDensity(sat);
      BubH := Medium.bubbleEnthalpy(sat);
      DewD := Medium.dewDensity(sat);
      DewH := Medium.dewEnthalpy(sat);
      Hv := DewH - BubH;
      BubS := Medium.bubbleEntropy(sat);
      DewS := Medium.dewEntropy(sat);
      Tsat := Medium.saturationTemperature(sat.psat);
      Sigma := Medium.surfaceTension(sat);
    equation
//Construction of BaseProperties
      BaseProp.p = p;
      BaseProp.h = H;
      der(p) = finalP - initialP;
      annotation(
        Documentation(info = "<html>
        <body>
        <p>This test maintains a constant temperature, and ramps pressure between the selected values. At each pressure, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p>
        </body>
        </html>"));
    end FluidTestingB;

    model TestB1A
      extends FluidTestingB(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), T = 100 + 273.15, initialP = 100.0e5, fract = 1.0, finalP = 1.1e3);
    end TestB1A;

    model TestB1B
      extends TestB1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model");
    end TestB1B;

    model TestWaterFlow1DFV_A_FF "Test case for Water.Flow1DFV. Copied from ThermoPower. Added full path to elements. Changed Medium to FreeFluids.TMedia.Fluids.Water"
      extends Modelica.Icons.Example;
      replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = false) constrainedby Modelica.Media.Interfaces.PartialMedium;
      parameter Integer Nnodes = 20 "Number of nodes";
      parameter Modelica.SIunits.Length Lhex = 10 "Total length";
      parameter Modelica.SIunits.Diameter Dihex = 0.02 "Internal diameter";
      final parameter Modelica.SIunits.Radius rhex = Dihex / 2 "Internal radius";
      final parameter Modelica.SIunits.Length omegahex = Modelica.Constants.pi * Dihex "Internal perimeter";
      final parameter Modelica.SIunits.Area Ahex = Modelica.Constants.pi * rhex ^ 2 "Internal cross section";
      parameter Real Cfhex = 0.005 "Friction coefficient";
      parameter Modelica.SIunits.MassFlowRate whex = 0.31 "Nominal mass flow rate";
      parameter Modelica.SIunits.Pressure phex = 2e5 "Initial pressure";
      parameter Modelica.SIunits.SpecificEnthalpy hinhex = 1e5 "Initial inlet specific enthalpy";
      parameter Modelica.SIunits.SpecificEnthalpy houthex = 1e5 "Initial outlet specific enthalpy";
      parameter Modelica.SIunits.SpecificEnthalpy deltah = 41800 "Height of enthalpy step";
      parameter Modelica.SIunits.EnergyFlowRate W = 41800 * whex "Height of power step";
      Modelica.SIunits.Time tau "Transport time delay";
      ThermoPower.Water.SourceMassFlow fluidSource(p0 = phex, h = hinhex, w0 = whex, use_in_w0 = true, use_in_h = true, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-78, -10}, {-58, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure fluidSink(p0 = phex / 2, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{70, -10}, {90, 10}}, rotation = 0)));
      ThermoPower.Water.ValveLin valve(Kv = 3e-6, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{10, -10}, {30, 10}}, rotation = 0)));
      ThermoPower.Water.Flow1DFV hex(N = Nnodes, L = Lhex, omega = omegahex, Dhyd = Dihex, A = Ahex, wnom = whex, Cfnom = Cfhex, DynamicMomentum = false, hstartin = hinhex, hstartout = houthex, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, initOpt = ThermoPower.Choices.Init.Options.steadyState, pstart = phex, dpnom = 1000, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-20, -10}, {0, 10}}, rotation = 0)));
      ThermoPower.Water.SensT T_in(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-50, -6}, {-30, 14}}, rotation = 0)));
      ThermoPower.Thermal.HeatSource1DFV heatSource(Nw = Nnodes - 1) annotation(
        Placement(transformation(extent = {{-20, 22}, {0, 42}}, rotation = 0)));
      Modelica.Blocks.Sources.Step MassFlowRate(offset = whex, startTime = 50, height = -0.03) annotation(
        Placement(transformation(extent = {{-98, 28}, {-78, 48}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant Constant1(k = 1) annotation(
        Placement(transformation(extent = {{-10, 60}, {10, 80}}, rotation = 0)));
      Modelica.Blocks.Sources.Step InSpecEnthalpy(height = deltah, offset = hinhex, startTime = 1) annotation(
        Placement(transformation(extent = {{-90, 60}, {-70, 80}}, rotation = 0)));
      Modelica.Blocks.Sources.Step ExtPower(height = W, startTime = 30) annotation(
        Placement(transformation(extent = {{-40, 40}, {-20, 60}}, rotation = 0)));
      ThermoPower.Water.SensT T_out(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{40, -6}, {60, 14}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(transformation(extent = {{80, 80}, {100, 100}})));
    equation
      tau = sum(hex.rho) / Nnodes * Lhex * Ahex / whex;
      connect(hex.outfl, valve.inlet) annotation(
        Line(points = {{0, 0}, {6, 0}, {10, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(T_in.outlet, hex.infl) annotation(
        Line(points = {{-34, 0}, {-28, 0}, {-20, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(fluidSource.flange, T_in.inlet) annotation(
        Line(points = {{-58, 0}, {-46, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(heatSource.wall, hex.wall) annotation(
        Line(points = {{-10, 29}, {-10, 5}}, color = {255, 127, 0}));
      connect(T_out.outlet, fluidSink.flange) annotation(
        Line(points = {{56, 0}, {70, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(valve.outlet, T_out.inlet) annotation(
        Line(points = {{30, 0}, {44, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(MassFlowRate.y, fluidSource.in_w0) annotation(
        Line(points = {{-77, 38}, {-72, 38}, {-72, 6}}, color = {0, 0, 127}));
      connect(InSpecEnthalpy.y, fluidSource.in_h) annotation(
        Line(points = {{-69, 70}, {-64, 70}, {-64, 6}}, color = {0, 0, 127}));
      connect(ExtPower.y, heatSource.power) annotation(
        Line(points = {{-19, 50}, {-10, 50}, {-10, 36}}, color = {0, 0, 127}));
      connect(Constant1.y, valve.cmd) annotation(
        Line(points = {{11, 70}, {20, 70}, {20, 8}}, color = {0, 0, 127}));
      annotation(
        Diagram(graphics),
        experiment(StopTime = 80, Tolerance = 1e-006),
        Documentation(info = "<html>
    <p>The model is designed to test the component <code>Water.Flow1DFV</code> (fluid side of a heat exchanger, finite volumes).</p>
    <p>This model represent the fluid side of a heat exchanger with an applied external heat flow. The operating fluid is liquid water.</p>
    <p>During the simulation, the inlet specific enthalpy, heat flux and mass flow rate are changed. The outlet temperature can be predicted analytically assuming incompressible flow and constant cp.</p>
    <p><ul>
    <li>t=0 s, Step variation of the specific enthalpy of the fluid entering the heat exchanger. The outlet temperature should undergo a step increase of 10 degrees 10 s later. </li>
    <li>t=30 s, Step variation of the thermal flow entering the heat exchanger lateral surface. The outlet temperature should undergo a ramp increase of 10 degrees lasting 10 s </li>
    <li>t=50 s, Step reduction of the mass flow rate entering the heat exchanger. The outlet temperature should undergo a ramp change of one degree lasting 10s</li>
    </ul></p>
    <p>Simulation Interval = [0...80] sec </p>
    <p>Integration Algorithm = DASSL </p>
    <p>Algorithm Tolerance = 1e-6 </p>
    </html>", revisions = "<html>
    <p><ul>
    <li>12 Sep 2013 by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br/>Updated to new FV structure. Updated parameters.</li></li>
    <li><i>1 Oct 2003</i> by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br/>First release.</li>
    </ul></p>
    </html>"));
    end TestWaterFlow1DFV_A_FF;
  end Tests;
  annotation(
    Documentation(info = "<html>
    <body>
    <p>The medium is designed for liquid, liquid/vapor, or gas phases, at temperatures lower than 0.85 Tc, and pressure not higher than 200 bars, because at higher values the influence of pressure on properties becomes very difficult to correct for. It extends the Modelica PartialTwoPhaseMedium. It is somewhat similar to the Modelica TableBased medium, but uses specific correlations for each physical property, allows to work with gas phase, and adds a density dependent correlation for the reduced bulk modulus of the liquid, that improves a lot the calculation of liquid density at high pressure, isothermal compressibility, and isobaric expansion coefficient. Improving also the calculation of liquid heat capacity at constant volume (Cv) and the speed of sound. The medium properties are obtained using correlations that are mainly functions of T, but different pressure corrections are also used. It uses the substances data stored in the MediaCommon package.</p>
    <p>The use of pressure correction is controlled by the constant Boolean 'highPressure'. Its default value is false. If switched to true, pressure correction will be applied (in plus than to liquid density) to liquid specific enthalpy, specific entropy, heat capacity, viscosity and thermal conductivity. It is interesting to make highPressure=true if we need to work over 20 or 30 bars, but the price is a slower simulation.</p>
    <p>The values of enthalpy and entropy are calculated from a reference state. The reference state to use can be selected giving value to the constant string 'refState'. The values can be: 'ASHRAE', 'IIR', 'NBP' or 'User'. Any other value will eliminate any correction for the reference state. When using 'User', the raw values at reference_T will be used as zero for both enthalpy and entropy.</p>
    <p>The thermodynamic record contains: p,T,gas fraction, d and h. Care must be taken in limiting the use to the temperature limits of the correlations used, as only few checks are done by the media, in order not to interfere with the solver process.</p>
    <p>A constant string 'localInputChoice' has been added to the BaseProperties model in order to specify the independent variables to use in each instance of the BaseProperties model. The default value for this constant is the value given to the constant string 'inputChoice' at package level. The valid alternatives are: 'ph', 'pT', 'dT'.</p>
    <p>In the package Tests there is a comparison between the medium performance with water and the Modelica WaterIF97_ph medium model (TestA1A/B and TestB1A/B). And with R134A and the Modelica R134a_ph model. The Modelica R134a_ph model seems incorrect. There is also an example taken from ThermoPower(you will need to load the ThermoPower package).</p>
    <p>The global idea has been not to use the Modelica files for the storage of substances data, but to store the data in a database, from which we can recover and use them when needed. A database is provided with more than 400 substances that can be enlarged, and the FreeFluids GUI can retrieve the data from the data base, treat it as needed (for example creating EOS from saturated vapor pressure and/or densities, or creating correlations from the EOS), store the results in the database, and export the data in Modelica format when needed. Nevertheless, in order to make life easier for users, many common substances have been exported, and their packages included in the TMedia.Fluids package.</p>
    <p>As a resume: The medium is for fast calculation of liquid phase, condensation, evaporation, and gas phase below the critical point. In the liquid and saturated phases, the results are quite good. In the gas phase, the results are better than the ideal gas approach in density and enthalpy. The medium is compatible with OpenModelica 1.14 old and new frontends. The medium is also compatible, since the addition of derivative functions calculation, with the ThermoPower library.</p>
    <b><p>Liquid phase properties</p></b>                
    <p>The saturated density is calculated using a dedicated correlation. This density is corrected for pressure influence. If the coefficients for the reduced bulk modulus calculation, as function of density, are available, the correction is done using a very accurate Tait like equation. In other case, a substance specific isothermal compressibility factor (with a default value of 6.667 e-10) is used. The parameters for the reduced bulk modulus correlation (with liquid density as independent variable) are normally not available, but can be calculated from a good equation of state of the multiparameter or SAFT types. This can be done easily with the FreeFluids GUI: you make the calculation with the EOS, transfer the results (density and the natural logarithm of the reduced bulk modulus) to the correlations tab, make the regression of the coefficients, and store the result in the database. It is good to calculate the reduced bulk modulus at a pressure close to 50 bars but, if necessary in order to have liquid phase at the temperature of interest, it can be done at higher pressure. Check that all density data used correspond to the liquid state.</p>
    <p>A dedicated correlation is used also for the saturated heat capacity. The use of the saturated liquid Cp correlation, instead of ideal Cp correlation plus vaporization enthalpy, makes possible the use of the medium with substances for which we do not have Cp0 data, and improves the liquid phase thermal properties calculation. This correlation can be also a problem, as many times we only find it with a temperature limit of the normal boiling point. This can be solved using a Cp correlation constructed from a good EOS, using the FreeFluids GUI. It is important not to use data too close to the Tc for the regression (use data till 10 K below the Tc). The best equation for the regression of the liquid heat capacity is the PPDS15 equation. Do not use the ChemSep equations as they are not integrated by the medium to obtain enthalpy or entropy. If highPressure has been made equal to true, a density dependent correction is applied to the Cp.</p>
    <p>The liquid enthalpy is calculated from the liquid Cp correlation at saturation. If highPressure has been made equal to true, PV correction is applied. The correction interpolates linearly from full correction below 0.45Tc to none at 0.85Tc.</p>
    <p>When going back from liquid enthalpy, or entropy, to temperature, we have two ways: One of them is to fit a correlation that makes this calculation, that can again be made with FreeFluidsGui. The second, that will be used if we left the value of lTfromHsatCorr=0, will calculate in situ the temperature using a solving function, as is done in the IdealGasMedia package.</p>
    <p>For water, viscosity is calculated, independent of the phase, using the reference equation as function of temperature and density. For other substances, the saturated liquid viscosity is calculated using a dedicated temperature dependent correlation. Pressure correction, according to Lucas, is applied if highPressure is set to true.</p> 
    <p>Saturated thermal conductivity is calculated using the corresponding correlation, if available. Otherwise the Latini method is used. If highPressure is set to true, a pressure(in fact density) correction is used.</p>
    <p>Surface tension is calculated also by correlation, or using the Sastri-Rao method, if the first is not available.</p> 
    <b><p>Gas phase properties</p></b>
    <p>Saturated gas density is calculated using a dedicated correlation. At lower pressure, a linear extrapolation is performed between the saturated density and the ideal gas density. The extrapolation range is between the vapor pressure and 10% of its value.</p>
    <p>Enthalpy and entropy are calculated from the saturated values at the given pressure (obtained from liquid Cp and vaporization enthalpy correlations), by adding the increase from the saturated temperature to the given temperature, calculated according to the ideal gas Cp. A pressure correction is used for enthalpy, but not for entropy or others variables.</p> 
    <p>Cp is calculated according to the ideal gas heat capacity correlation with a pressure correction.</p>
    <p>Viscosity and thermal conductivity are calculated by correlations or, when not available, by the Chung method.</p>
    <b><p>Two phases properties</p></b>
    <p>Transport properties are not calculated for the two phases situation.</p>
    </body>
    </html>"));
end TMedia;

within FreeFluids;

package TMedia "TMedia.mo by Carlos Trujillo
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
  partial package TMedium
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(onePhase = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, reference_T = 298.15, redeclare record FluidConstants = FreeFluids.MediaCommon.DataRecord);
    extends FreeFluids.MediaCommon.Types;
    constant String refState = "IIR" "Enthalpy/entropy reference state. Alternatives: ASHRAE, IIR, NBP, otherwise the value at reference_T as 0";
    constant InputChoice inputChoice = InputChoice.ph "Allows to choose the input choise to use for the construction of a BaseProperties object";
    constant Boolean highPressure=false;
    constant SpecificEnthalpy UserRefH = SpecificEnthalpyCorr(reference_T);
    constant SpecificEnthalpy ASHRAErefH = SpecificEnthalpyCorr(233.15);
    constant SpecificEnthalpy NBPrefH = SpecificEnthalpyCorr(fluidConstants[1].Tb);
    constant SpecificEnthalpy IIRrefH = SpecificEnthalpyCorr(273.15) - 2.0e5;
    constant SpecificEntropy UserRefS = SpecificEntropyCorr(reference_T);
    constant SpecificEntropy ASHRAErefS = SpecificEntropyCorr(233.15);
    constant SpecificEntropy NBPrefS = SpecificEntropyCorr(fluidConstants[1].Tb);
    constant SpecificEntropy IIRrefS = SpecificEntropyCorr(273.15) - 1.0e3;
    //Auxiliary functions based in correlations
    //-----------------------------------------

    function PhysPropCorr "Calculates a physical property from a given correlation data and a given temperature, or pressure, using an external function."
      input Integer corr;
      input Real coef[:];
      input Real x;
      output Real y;
    
      external "C" FF_PhysPropCorr(corr, coef, fluidConstants[1].MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end PhysPropCorr;

    function SpecificEnthalpyCorr "Calculates specific enthalpy from a given Cp correlation at a given temperature from 0K"
      input Real x;
      output Real y;
    
      external "C" FF_SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
      annotation(
        Inline = true,
        smoothOrder = 2);
    end SpecificEnthalpyCorr;

    function SpecificEntropyCorr "Calculates specific entropy from a given Cp correlation at a given temperature from 0K"
      input Real x;
      output Real y;
    
      external "C" FF_SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end SpecificEntropyCorr;

    //Function SpecificEnthalpyCorrInv

    function SpecificEnthalpyCorrInv "Compute temperature from property value"
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal "Solve h(fluidConstants[1],T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Media.Common.OneNonLinearEquation;

        redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
        end f_nonlinear_Data;

        redeclare function extends f_nonlinear
          algorithm
            y := SpecificEnthalpyCorr(x);
        end f_nonlinear;

        redeclare function extends solve
        end solve;
      end Internal;

      Internal.f_nonlinear_Data fd;
    algorithm
      T := Internal.solve(y, 50.0, fluidConstants[1].Tc - 0.01, 1.0e5, {1}, fd);
    end SpecificEnthalpyCorrInv;

    //Function SpecificEntropyCorrInv

    function SpecificEntropyCorrInv "Compute temperature from property value"
      //input Correlation corrData;
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal "Solve h(fluidConstants[1],T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Media.Common.OneNonLinearEquation;

        redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
        end f_nonlinear_Data;

        redeclare function extends f_nonlinear
          algorithm
            y := SpecificEntropyCorr(x);
        end f_nonlinear;

        redeclare function extends solve
        end solve;
      end Internal;

      Internal.f_nonlinear_Data fd;
    algorithm
      T := Internal.solve(y, 50.0, fluidConstants[1].Tc - 0.01, 1.0e5, {1}, fd);
    end SpecificEntropyCorrInv;

    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;

      algorithm
        assert(T < fluidConstants[1].Tc, "The media can´t be used over Tc");
        p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, T);
    end saturationPressure;

    function saturationPressureInv "Compute temperature from property value"
      input Real y;
      output Real T;
    protected
      package Internal
        extends Modelica.Media.Common.OneNonLinearEquation;

        redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
        end f_nonlinear_Data;

        redeclare function extends f_nonlinear
          algorithm
            y := saturationPressure(x);
        end f_nonlinear;

        // Dummy definition has to be added for current Dymola

        redeclare function extends solve
        end solve;
      end Internal;

      Internal.f_nonlinear_Data fd;
    algorithm
      T := Internal.solve(y, 50.0, fluidConstants[1].Tc - 0.01, 1.0e5, {1}, fd);
    end saturationPressureInv;

    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;

      algorithm
        assert(p < fluidConstants[1].Pc, "No saturated media possible over Pc");
        T := if fluidConstants[1].BtCorr > 0 then PhysPropCorr(fluidConstants[1].BtCorr, fluidConstants[1].BtCoef, p) else saturationPressureInv(p);

    end saturationTemperature;

    //BaseProperties model
    //--------------------

    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default))
        parameter InputChoice localInputChoice = inputChoice;

      equation
        MM = fluidConstants[1].MW / 1000 "in kg/mol";
        R = Modelica.Constants.R / MM;
        u = h - p / d;
        sat = setSat_p(p);
        if localInputChoice == InputChoice.ph then
          state = setState_phX(p, h, fill(0, 0));
          state.T = T;
          state.d = d;
        elseif localInputChoice == InputChoice.pT then
          state = setState_pTX(p, T, fill(0, 0));
          state.d = d;
          state.h = h;
        elseif localInputChoice == InputChoice.dT then
          state = setState_dTX(d, T, fill(0, 0));
          state.p = p;
          state.h = h;
        else
          assert(false, "Invalid choice for BaseProperties inputChoice");
        end if;
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
        AbsolutePressure Vp;
        Density Dls "saturated liquid density at T";

      algorithm
        assert(T <= fluidConstants[1].Tc, "The media can´t be used over Tc.");
        state.p := p;
        state.T := T;
        Vp := saturationPressure(T);
        if p > Vp then
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp));
          state.h := if refState == "IIR" then SpecificEnthalpyCorr(T) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(T) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(T) - ASHRAErefH else SpecificEnthalpyCorr(T) - UserRefH "saturated liquid enthalpy";
          if highPressure==true then
            state.h:= if (T<=0.45*fluidConstants[1].Tc) then state.h+state.p/state.d-Vp/Dls else if  (T<0.85*fluidConstants[1].Tc) then state.h+(0.85*fluidConstants[1].Tc-T)/(0.4*fluidConstants[1].Tc)*(state.p/state.d-Vp/Dls) else state.h "pressure correction for liquid enthalpy";
          end if;
          state.phase := 1;
          state.gf := 0;
        elseif p < Vp then
          assert(false, "The medium can't work with non-saturated gas");
        else
          assert(false, "A two phases state can't be constructed from p and T");
        end if;
    end setState_pTX;

    redeclare function extends setState_dTX "Return ThermodynamicState record as function of T,d and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Density Dls;
        Density Dgs;
        SpecificEnthalpy Hls "enthalpy of the saturated liquid at given T";
        AbsolutePressure Vp;

      algorithm
        assert(T <= fluidConstants[1].Tc, "The media can´t be used over Tc");
        state.T := T;
        state.d := d;
        Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, T) "saturated liquid density at T";
        Hls := if refState == "IIR" then SpecificEnthalpyCorr(T) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(T) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(T) - ASHRAErefH else SpecificEnthalpyCorr(T) - UserRefH;
        Vp:=PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, T);
        if d >= Dls then
          state.p := if fluidConstants[1].lnuA > 0 then Vp + (exp((d - Dls) * fluidConstants[1].lnuA) - 1) * 1000 * Modelica.Constants.R * T * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) / (fluidConstants[1].lnuA * fluidConstants[1].MW) else (1 - Dls / d) / fluidConstants[1].IsothComp + Vp;
          state.h:= if (T<=0.45*fluidConstants[1].Tc) then Hls+(state.p-Vp)/state.d else if  (T<0.85*fluidConstants[1].Tc) then Hls+(0.85*fluidConstants[1].Tc-state.T)/(0.4*fluidConstants[1].Tc)*(state.p-Vp)/state.d else Hls;
          state.phase := 1;
          state.gf := 0;
        else
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, T);
          assert(d >= Dgs, "The media can´t work with non-saturated gases");
          state.p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, T);
          state.h:= if (T<=0.45*fluidConstants[1].Tc) then Hls+state.p/state.d-Vp/Dls else if  (T<0.85*fluidConstants[1].Tc) then Hls+(0.85*fluidConstants[1].Tc-state.T)/(0.4*fluidConstants[1].Tc)*(state.p/state.d-Vp/Dls) else Hls "pressure correction for liquid enthalpy";
          if d > Dgs then
            state.phase := 2;
            state.gf := Dgs * (Dls - d) / (d * (Dls - Dgs));
            state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, T);
          elseif d == Dgs then
            state.phase := 1;
            state.gf := 1;
            state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, T);
          end if;
        end if;
    end setState_dTX;

    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,H and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEnthalpy Hls "saturated liquid enthalpy";
        SpecificEnthalpy Hgs "saturated gas enthalpy";
        Density Dls "saturated liquid density";
        Density Dgs "saturated gas density";
        SpecificEnthalpy Hraw "Input enthalpy without reference correction";
        AbsolutePressure Vp;

      algorithm
        state.p := p;
        state.h := h;
        if p > fluidConstants[1].Pc then
          Tb := fluidConstants[1].Tc;
        else
          Tb := saturationTemperature(p);
        end if;
        Hraw := if refState == "IIR" then h + IIRrefH else if refState == "NBP" then h + NBPrefH else if refState == "ASHRAE" then h + ASHRAErefH else h + UserRefH;
        Hls := SpecificEnthalpyCorr(Tb) "unreferenced saturated liquid specific enthalpy at T";
//Hls := SpecificEnthalpyCorr(Tb) - UserRefH "liquid specific enthalpy at boiling T";
        if Hls >= Hraw then
          state.T := if fluidConstants[1].lTfromHsatCorr > 0 then PhysPropCorr(fluidConstants[1].lTfromHsatCorr, fluidConstants[1].lTfromHsatCoef, Hraw) else SpecificEnthalpyCorrInv(Hraw);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T) "saturated liquid density";
          Vp:=PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          if highPressure==true then
            Hraw:= if (state.T<=0.45*fluidConstants[1].Tc) then Hraw-state.p/state.d+Vp/Dls else if  (state.T<0.85*fluidConstants[1].Tc) then Hraw-(0.85*fluidConstants[1].Tc-state.T)/(0.4*fluidConstants[1].Tc)*(state.p/state.d-Vp/Dls) else Hraw "pressure correction for liquid enthalpy";
            //Hraw:=Hraw-(state.p-Vp)/state.d "pressure correction for liquid enthalpy";
            state.T := if fluidConstants[1].lTfromHsatCorr > 0 then PhysPropCorr(fluidConstants[1].lTfromHsatCorr, fluidConstants[1].lTfromHsatCoef, Hraw) else SpecificEnthalpyCorrInv(Hraw);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T) "saturated liquid density";
            Vp:=PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T);
            state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          end if;
          state.phase := 1;
          state.gf := 0;
        else
          Hgs := Hls + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, Tb);
          assert(Hgs >= Hraw, "The media can´t work with non-saturated gases");
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, Tb);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, Tb);
          state.T := Tb;
          state.phase := 2;
          state.gf := (Hraw - Hls) / (Hgs - Hls);
          state.d := Dls * Dgs / (state.gf * Dls + (1 - state.gf) * Dgs);
        end if;
    end setState_phX;

    redeclare function extends setState_psX "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEntropy Sl;
        SpecificEntropy Sg;
        Density Dls;
        Density Dgs;
        SpecificEntropy Sraw "Input entropy without reference correction";
        AbsolutePressure Vp;

      algorithm
        state.p := p;
        if p > fluidConstants[1].Pc then
          Tb := fluidConstants[1].Tc;
        else
          Tb := saturationTemperature(p);
        end if;
        Sraw := if refState == "IIR" then s + IIRrefS else if refState == "NBP" then s + NBPrefS else if refState == "ASHRAE" then s + ASHRAErefS else s + UserRefS;
        Sl := SpecificEntropyCorr(Tb) "unreferenced saturated liquid specific entropy at T";
        if Sl >= Sraw then
          state.T := SpecificEntropyCorrInv(Sraw);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T) "saturated liquid density";
          Vp:= PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          if highPressure==true then
            Sraw:=if state.T<=0.45*fluidConstants[1].Tc then Sraw else Sraw+(state.T-0.45*fluidConstants[1].Tc) * (state.p-Vp)/(0.4*fluidConstants[1].Tc*state.d*state.T) "pressure correction for liquid entropy";
            state.T := SpecificEntropyCorrInv(Sraw);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T) "saturated liquid density";
            Vp:= PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T);
            state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          end if;
          state.h := if refState == "IIR" then SpecificEnthalpyCorr(state.T) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(state.T) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(state.T) - ASHRAErefH else SpecificEnthalpyCorr(state.T) - UserRefH "saturated liquid enthalpy";
          if highPressure==true then         
            state.h:= if (state.T<=0.45*fluidConstants[1].Tc) then state.h+state.p/state.d-Vp/Dls else if  (state.T<0.85*fluidConstants[1].Tc) then state.h+(0.85*fluidConstants[1].Tc-state.T)/(0.4*fluidConstants[1].Tc)*(state.p/state.d-Vp/Dls) else state.h "pressure correction for liquid enthalpy";
          end if;
          state.phase := 1;
          state.gf := 0;
        else
          Sg := Sl + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, Tb) / Tb;
          assert(Sg >= Sraw, "The media can´t work with non-saturated gases");
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, Tb);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, Tb);
          state.T := Tb;
          state.phase := 2;
          state.gf := (Sraw - Sl) / (Sg - Sl);
          state.h := if refState == "IIR" then SpecificEnthalpyCorr(state.T) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(state.T) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(state.T) - ASHRAErefH else SpecificEnthalpyCorr(state.T) - UserRefH "saturated liquid enthalpy";
          state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, Tb) / Tb;
          state.d := Dls * Dgs / (state.gf * Dls + (1 - state.gf) * Dgs);
        end if;
    end setState_psX;

    //Special thermodynamic states constructors
    //------------------------------------------

    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;

      protected
        AbsolutePressure Vp;

      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.h := if refState == "User" then SpecificEnthalpyCorr(sat.Tsat) - UserRefH else if refState == "NBP" then SpecificEnthalpyCorr(sat.Tsat) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(sat.Tsat) - ASHRAErefH else if refState == "IIR" then SpecificEnthalpyCorr(sat.Tsat) - IIRrefH else 0.0;
        Vp := saturationPressure(state.T);
        state.d := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T);
        state.gf := 0.0 "liquid";
        state.phase := 1;
    end setBubbleState;

    redeclare function extends setDewState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the dew point"
        extends Modelica.Icons.Function;

      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.h := if refState == "User" then SpecificEnthalpyCorr(sat.Tsat) - UserRefH else if refState == "NBP" then SpecificEnthalpyCorr(sat.Tsat) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(sat.Tsat) - ASHRAErefH else if refState == "IIR" then SpecificEnthalpyCorr(sat.Tsat) - IIRrefH else 0.0;
        state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, sat.Tsat);
        state.d := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, sat.Tsat);
        state.gf := 1.0 "gas";
        state.phase := 1;
    end setDewState;

    //Properties calculation from thermodynamic state
    //-----------------------------------------------

    function molarMass "Return the molar mass of the medium"
      extends Modelica.Icons.Function;
      input ThermodynamicState state;
      output Real MM;
    algorithm
      MM := fluidConstants[1].MW * 1e-3;
    end molarMass;

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
        u := state.h - state.p / state.d;
    end specificInternalEnergy;

    redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;

      protected
        SpecificEntropy Sl "entropy of the liquid phase";
        AbsolutePressure Vp;

      algorithm
        Sl := if refState == "IIR" then SpecificEntropyCorr(state.T) - IIRrefS else if refState == "NBP" then SpecificEntropyCorr(state.T) - NBPrefS else if refState == "ASHRAE" then SpecificEntropyCorr(state.T) - ASHRAErefS else SpecificEntropyCorr(state.T) - UserRefS;
        if state.gf == 0.0 then
          s:=Sl;
          if highPressure==true then
            Vp:=PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T);
            s := if state.T<=0.45*fluidConstants[1].Tc then s else s-(state.T-0.45*fluidConstants[1].Tc) * (state.p-Vp)/(0.4*fluidConstants[1].Tc*state.d*state.T)"pressure correction for liquid entropy";
          end if;
        else
          s := Sl + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, state.T) / state.T;
        end if;
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
        Real vp1, vp2, dv, ds;

      algorithm
        if state.gf == 0.0 then
          cp := PhysPropCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, state.T);
          if highPressure==true then
            ds := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T);
            cp:= if fluidConstants[1].family==17 then cp*exp(-5.0*(state.d-ds)^0.5/(ds*(1-state.T/fluidConstants[1].Tc)^0.5)) else cp*exp(-2.8*(state.d-ds)^0.5/(ds*(1-state.T/fluidConstants[1].Tc)^0.5));
          end if;
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
          d1 := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T - 0.1);
          d2 := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T + 0.1);
          d0 := (d1 + d2) / 2;
          beta := (d1 - d2) * d0 / (d1 * d2 * 0.2) "Coefficient along the saturation line, function only of T, not the real beta";
          vp1 := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T - 0.1);
          vp2 := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T + 0.1);
          dv := (vp2 - vp1) / 0.2;
          beta := beta + isothermalCompressibility(state) * dv "Beta at saturation, once corrected for pressure variation";
          if fluidConstants[1].lnuA > 0 then
            beta := beta * d0 * exp(-fluidConstants[1].lnuA * (state.d - d0)) / state.d + (1 - exp(-fluidConstants[1].lnuA * (state.d - d0))) / (fluidConstants[1].lnuA * state.d * state.T) "Beta at the given T and d. Postnikov et alt. 2016";
          end if;
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
        else
          kappa := 0.0;
        end if;
    end isothermalCompressibility;

    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;

      protected
        Real cp, kappa, beta;
        //Real cv;

      algorithm
        beta := isobaricExpansionCoefficient(state);
        cp := specificHeatCapacityCp(state);
//cv := specificHeatCapacityCv(state);
        kappa := isothermalCompressibility(state);
//a := sqrt(max(0, cp / (cv * kappa * state.d)));
        a := sqrt(max(0, cp / (cp * kappa * state.d - state.T * beta * beta)));
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
          if fluidConstants[1].CAS=="7732-18-5" then
            eta:=FreeFluids.MediaCommon.Functions.waterViscosity(state.T,state.d);
          elseif fluidConstants[1].lViscCorr > 0 then
            eta := if highPressure==false  then PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, state.T) else FreeFluids.MediaCommon.Functions.liqViscPcorLucas(fluidConstants[1], state.T, state.p, PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, state.T), PhysPropCorr(fluidConstants[1].lViscCorr, fluidConstants[1].lViscCoef, state.T));
          else
            assert(false, "there is no correlation for liquid viscosity calculation");
          end if;
        elseif state.phase == 1.0 then
          if fluidConstants[1].CAS=="7732-18-5" then
            eta:=FreeFluids.MediaCommon.Functions.waterViscosity(state.T,state.d);
          else
            lpVisc:= if fluidConstants[1].gViscCorr > 0 then PhysPropCorr(fluidConstants[1].gViscCorr, fluidConstants[1].gViscCoef, state.T) else FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(fluidConstants[1], state.T);
            eta:=lpVisc "no pressure correction seems necessary for gas phase";
            //eta:= if highPressure==false then lpVisc else  FreeFluids.MediaCommon.Functions.gasViscPcorLucas(fluidConstants[1], state.T, state.p, lpVisc);
          end if;
        else
          eta := 0 "no viscosity calculation is done for two phases situation";
        end if;
    end dynamicViscosity;

    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf == 0.0 then
          lambda := if fluidConstants[1].lThCondCorr > 0 then PhysPropCorr(fluidConstants[1].lThCondCorr, fluidConstants[1].lThCondCoef, state.T) else FreeFluids.MediaCommon.Functions.liqThCondLatini(fluidConstants[1], state.T);
          if highPressure==true then 
            lambda:=lambda*exp(state.d/PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, state.T)-1);
 end if;
        elseif state.gf == 1.0 then
          lambda := if fluidConstants[1].gThCondCorr > 0 then PhysPropCorr(fluidConstants[1].gThCondCorr, fluidConstants[1].gThCondCoef, state.T) else FreeFluids.MediaCommon.Functions.gasThCondLowPressureChung(fluidConstants[1], specificHeatCapacityCp(state), dynamicViscosity(state), state.T);
        else
          lambda := 0 "no thermal conductivity calculation is done in two phases situation";
        end if;
    end thermalConductivity;

    redeclare function extends surfaceTension "Return surface tension"
        extends Modelica.Icons.Function;

      algorithm
        sigma := if fluidConstants[1].lSurfTensCorr>0 then PhysPropCorr(fluidConstants[1].lSurfTensCorr, fluidConstants[1].lSurfTensCoef, sat.Tsat) else FreeFluids.MediaCommon.Functions.liqSurfTensSastriRao(fluidConstants[1], sat.Tsat);
    end surfaceTension;

    //Saturation properties
    //-------------------------

    redeclare function extends bubbleDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;

      algorithm
        dl := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, sat.Tsat);
    end bubbleDensity;

    redeclare function extends dewDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;

      algorithm
        dv := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, sat.Tsat);
    end dewDensity;

    redeclare function extends bubbleEnthalpy "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        hl := if refState == "IIR" then SpecificEnthalpyCorr(sat.Tsat) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(sat.Tsat) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(sat.Tsat) - ASHRAErefH else SpecificEnthalpyCorr(sat.Tsat) - UserRefH;
    end bubbleEnthalpy;

    redeclare function extends dewEnthalpy "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        hv := if refState == "IIR" then SpecificEnthalpyCorr(sat.Tsat) - IIRrefH else if refState == "NBP" then SpecificEnthalpyCorr(sat.Tsat) - NBPrefH else if refState == "ASHRAE" then SpecificEnthalpyCorr(sat.Tsat) - ASHRAErefH else SpecificEnthalpyCorr(sat.Tsat) - UserRefH;
        hv := hv + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, sat.Tsat);
    end dewEnthalpy;

    redeclare function extends bubbleEntropy "Return bubble point specific entropy"
        extends Modelica.Icons.Function;

      algorithm
        sl := if refState == "IIR" then SpecificEntropyCorr(sat.Tsat) - IIRrefS else if refState == "NBP" then SpecificEntropyCorr(sat.Tsat) - NBPrefS else if refState == "ASHRAE" then SpecificEntropyCorr(sat.Tsat) - ASHRAErefS else SpecificEntropyCorr(sat.Tsat) - UserRefS;
    end bubbleEntropy;

    redeclare function extends dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;

      algorithm
        sv := if refState == "IIR" then SpecificEntropyCorr(sat.Tsat) - IIRrefS else if refState == "NBP" then SpecificEntropyCorr(sat.Tsat) - NBPrefS else if refState == "ASHRAE" then SpecificEntropyCorr(sat.Tsat) - ASHRAErefS else SpecificEntropyCorr(sat.Tsat) - UserRefS;
        sv := sv + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, sat.Tsat) / sat.Tsat;
    end dewEntropy;
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
  
    package Dichlorodifluormethane
      extends FreeFluids.TMedia.TMedium(final mediumName = "Dichlorodifluormethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Dichlorodifluormethane}, reference_T = 273.15, refState = "ASHRAE");
    end Dichlorodifluormethane;
  
    package EG
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ethylene glycol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.EG});
    end EG;
  
    package Ethane
      extends TMedium(final mediumName = "Ethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethane}, refState="IIR", reference_T=273.15);
    end Ethane;
  
    package Ethanol
      extends FreeFluids.TMedia.TMedium(final mediumName = "Ethanol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethanol}, reference_T = 273.15);
    end Ethanol;
  
    package Isobutane
      extends FreeFluids.TMedia.TMedium(final mediumName = "Isobutane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Isobutane}, reference_T = 273.15);
    end Isobutane;
  
    package MarlothermSH
      extends FreeFluids.TMedia.TMedium(final mediumName = "Marlotherm SH", final singleState = false, p_default = 1.0e5, T_default = 425.0, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.MarlothermSH});
    end MarlothermSH;
  
    package Methane
      extends TMedium(final mediumName = "Methane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Methane}, refState="IIR", reference_T=273.15);
    end Methane;
  
    package N2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Nitrogen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.N2});
    end N2;
  
    package O2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Oxygen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.O2});
    end O2;
  
    package Pentane_n
      extends TMedium(final mediumName = "Pentane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Pentane_n}, refState="IIR", reference_T=273.15);
    end Pentane_n;
  
    package Propane
      extends TMedium(final mediumName = "Propane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Propane}, refState = "NBP");
    end Propane;
  
    package R134A
      extends TMedium(final mediumName = "R134A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R134A}, refState="ASHRAE", reference_T=273.15);
    end R134A;
  
    package R410A
      extends TMedium(final mediumName = "R410A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R410A}, refState="ASHRAE", reference_T=273.15);
    end R410A;
  
    package ShellS2
      extends FreeFluids.TMedia.TMedium(final mediumName = "Shell S2", final singleState = false, p_default = 2.0e5, T_default = 425.0, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.ShellS2});
    end ShellS2;
  
    package Toluene
      extends FreeFluids.TMedia.TMedium(final mediumName = "Toluene", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Toluene}, reference_T = 273.15);
    end Toluene;
  
    package Water
      extends FreeFluids.TMedia.TMedium(final mediumName = "Water", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Water}, reference_T = 273.15, refState = "User");
    end Water;
  end Fluids;

  package Tests
    partial model FluidTest1
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.AbsolutePressure p = 1.0e5;
      parameter Medium.Temperature initialT = 273.19;
      parameter Medium.Temperature finalT = 372.15;
      parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
      Medium.Temperature T(start = initialT) "We will ramp the temperature";
      Medium.MolarMass MM;
      Medium.ThermodynamicState StateP "original state from p and T";
      Medium.ThermodynamicState StateH "state reproduced from p,h";
      Medium.ThermodynamicState StateS "state reproduced from p,s";
      Medium.ThermodynamicState StateD "state reproduced from d,T";
      Medium.ThermodynamicState StateDlow "state at fraction of original density";
      Medium.ThermodynamicState StateBub "bubble state at sat";
      Medium.ThermodynamicState StateDew "dew state at sat";
      Medium.ThermodynamicState StateHalfH "state half way between Bub and Dew";
      //Properties of StateP
      Real D, H;
      Real S;
      Real Beta "isobaric expansion coefficient";
      Real Kappa "isothermal compressibility";
      Real Cp;
      Real Cv;
      Real Gamma;
      Real SS "speed of sound";
      Real U "internal energy";
      Real G "Gibbs energy";
      Real A "Helmholtz energy";
      Real Mu;
      Real Th;
      //Properties of StateDlow
      Real DlowGasFract "gas fraction of StateDlow";
      Real DlowD, DlowH;
      Real DlowS;
      Real DlowU;
      Real DlowA;
      //Saturation properties
      Medium.SaturationProperties sat "saturation point at given T";
      Real Vp;
      Real BubD, BubH "bubble properties";
      Real BubS;
      Real DewD, DewH "dew properties";
      Real DewS;
      Real DewMu;
      Real Hv "vaporization enthalpy";
      Real Sigma;
      Medium.Temperature Tsat "temperature recovered from saturation pressure";
      //BaseProperties
      Medium.BaseProperties BaseProp;
    algorithm
//Construction of StateP and calculation of properties
    for i in 1:100 loop
      StateP := Medium.setState_pTX(p, T, fill(0, 0));
      MM := Medium.molarMass(StateP);
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
      DlowU := Medium.specificInternalEnergy(StateDlow);
      DlowS := Medium.specificEntropy(StateDlow);
      DlowA := Medium.specificHelmholtzEnergy(StateDlow);
//Calculations at saturation
      sat := Medium.setSat_T(T);
      Vp:=Medium.saturationPressure_sat(sat);
      StateBub := Medium.setBubbleState(sat);
      StateDew := Medium.setDewState(sat);
      BubD := Medium.bubbleDensity(sat);
      BubH := Medium.bubbleEnthalpy(sat);
      DewD := Medium.dewDensity(sat);
      DewH := Medium.dewEnthalpy(sat);
      Hv := DewH - BubH;
      BubS := Medium.bubbleEntropy(sat);
      DewS := Medium.dewEntropy(sat);
      DewMu:= Medium.dynamicViscosity(StateDew);
      Tsat := Medium.saturationTemperature(sat.psat);
      Sigma := Medium.surfaceTension(sat);
    //construction of state HalfH
      StateHalfH:=Medium.setState_phX(Vp, (BubH+DewH)/2, fill(0, 0));
      end for;
    equation
//Construction of BaseProperties
      BaseProp.p = p;
      BaseProp.h = H;
      der(T) = finalT - initialT;
    end FluidTest1;

    model Test1A
      extends FluidTest1(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState="User", reference_T=273.15, highPressure=true), p = 100.0e5, initialT = 0.1 + 273.15, fract = 0.3, finalT = 250.0 + 273.15);
    end Test1A;

    model Test1B
      extends Test1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model");
      //extends Test1A(redeclare package Medium = Modelica.Media.R134a.R134a_ph "Medium model");
    end Test1B;

    model Test2A
      extends FluidTest1(redeclare replaceable package Medium = TMedia.Fluids.O2(refState="NBP", highPressure=true), p = 50.0e5, initialT = 55.0, finalT = 120.0);
    end Test2A;

    partial model FluidTestingB
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.Temperature T = 273.2;
      parameter Medium.Temperature initialP = 1.0e5;
      parameter Medium.Temperature finalP = 1.0e7;
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
    end FluidTestingB;

    model TestB1A
      extends FluidTestingB(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState="User", reference_T=273.15, highPressure=true), T = 210 + 273.2, initialP = 20.0e5, fract = 1.0, finalP = 2.0e7);
    end TestB1A;

    model TestB1B
      extends TestB1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model");
    end TestB1B;
    
    package ModelicaTests
      model ModelicaR134ATest
        package Medium = Modelica.Media.R134a.R134a_ph;
        parameter Real T = 20 + 273.15;
        Medium.SaturationProperties sat;
        Medium.ThermodynamicState StateBub;
        Medium.ThermodynamicState StateDew;
        Medium.ThermodynamicState StateBi "a biphasic state with density halfway between bubble and dew";
        Medium.Density dBi "density half way between bubble and dew";
        Real gasFraction;
      algorithm
        sat := Medium.setSat_T(T);
        StateBub := Medium.setBubbleState(sat);
        StateDew := Medium.setDewState(sat);
        dBi := (StateBub.d + StateDew.d) / 2;
        StateBi := Medium.setState_dTX(dBi, T);
        gasFraction := Medium.vapourQuality(StateBi);
      end ModelicaR134ATest;
  end ModelicaTests;
  end Tests;
  annotation(
    Documentation(info = "<html>
    <body>
    <p>The medium is designed for liquid, and liquid/vapor phases, at a temperature lower than 0.85 Tc, and a pressure not higher than 200 bars, because in those situations properties are highly dependent on pressure. It extends the Modelica  PartialTwoPhaseMedium. The medium properties are obtained using correlations that are mainly functions of T, but different pressure corrections are also used. It is somewhat similar to the Modelica TableBased medium, but uses specific correlations for each physical property, allows to work with saturated gas phase, and adds a density dependent correlation for the reduced bulk modulus of the liquid, that improves a lot the calculation of liquid density at high pressure, isothermal compressiblity, and isobaric expansion coefficient. Improving also the calculation of heat capacity at constant volume (Cv) and the speed of sound.</p>
    <p>The obtained properties correspond to the liquid phase, or to the saturated gas phase. It can´t be used for non-saturated gas phase.</p>
    <p>If available, the reduced bulk modulus correlation for the liquid is used. In other case, a substance specific isothermal compressibility factor (with a default value of 6.667 e-10) is used. The parameters for the reduced bulk modulus correlation(with liquid density as independent variable) are normally not available, but can be calculated from a good equation of state of the multiparameter or SAFT types. This can be done easily with the FreeFluids GUI: you make the calculation with the EOS, transfer the results (density and the natural logarithm of the reduced bulk modulus) to the correlations tab, make the regresion of the coefficients, and store the result in the database. It is good to calculate the reduced bulk modulus at a pressure close to 50 bars but, if necessary in order to have liquid phase at the temperature of interest, it can be done at higher pressure. Check that all density data correspond to the liquid state.</p>
    <p>The liquid heat capacity correlation can be also a problem, as many times we only find it with a temperature limit of the normal boiling point. This can be solved using a Cp correlation constructed from a good EOS, using the FreeFluids GUI. It is important not to use data too close to the Tc for the regression (use data till 10 K below the Tc). The best equation for the regression of the liquid heat capacity is the PPDS15 equation. Do not use the Chemsep equations as they are not integrated by the medium to obtain enthalpy or entropy</p>
    <p>A The use of pressure correction is controlled by the constant Boolean 'highPressure'. Its default value is false. If switched to true, pressure correction will be applyied to liquid specific enthalpy, specific entropy, heat capcity, viscosity and thermal conductivity. No pressure correction is used for the gas phase. It is interesting to make highPressure=true if we need to work over 30 bars.</p>
    <p> The liquid enthalpy is calculated from the liquid Cp correlation at saturation. If highPressure has been made equal to true, a PV correction is applied. The correction comes from full correction below 0.45Tc to 0 at 0.85Tc.</p>
    <p>The values of enthalpy and entropy are calculated from a reference state. The reference state to use can be select giving value to the constant string 'refState'. The values can be:  'ASHRAE', 'IIR', 'NBP' or 'User'. When using 'User', the values at reference_T will be used as zero for both enthalpy and entropy.</p>
    <p>The use of the liquid Cp correlation is an alternative to the use of ideal Cp correlation plus vaporization enthalpy. This makes possible the use of the medium with substances for which we do not have Cp0 data.</p> 
    <p>When going back from liquid enthalpy, or entropy, to temperature, we have two ways: One of them is to fit a correlation that makes this calculation, that can again be made with FreeFluidsGui. The second, that will be used if we left the value of lTfromHsatCorr=0, will calculate insitu the temperature using a solving function, as is done in the IdealGasMedia package.</p>
    <p>For liquid transport properties, the calculation is done directly from T, but pressure corrections have been introduced for the liquid phase, as already stated.</p> 
    <p>The thermodynamic record contains: p,T,gas fraction, d and h. Care must be taken in limiting the use to the temperature limits of the correlations used, as only few checks are done by the media, in order not to interfer with the solver process</p>
    <p>A new enumeration type called InputChoice has been defined in order to specify the independent variables to use in each instance of the BaseProperties model. This is done assigning value to the parameter localInputChoice in the model, that has as default the value given to the constant InputChoice inputChoice at package level. </p>
    <p>In the package Tests there is a comparison between the medium performance with water and the WaterIF97_ph medium model.</p>
    <p>The global idea has been not to use the Modelica files for the storage of substances data, but to store the data in a database, from which we can recover and use them when needed. A database is provided with more than 400 substances that can be enlarged, and the FreeFluids GUI can retrieve the data from the data base, treat it as needed (for example creating EOS from saturated vapor pressure and/or densities, or creating correlations from the EOS), store the results in the database, and export the data in Modelica format when needed.</p>
    <p>As a resume: The medium is for fast calculation of liquid phase, condensation, and evaporation, but not for non-saturated gases.</p>
    </body>
    </html>"));
end TMedia;

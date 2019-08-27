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
  //*****MEDIUMS OF TEMPERATURE ONLY DEPENDENT PROPERTIES*****
  //==========================================================
  //***** PACKAGE TMedium*****
  //**************************

  partial package TMedium "Media for liquids or saturated gases, based in correlations"
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(onePhase = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
    constant MediaCommon.InputChoice inputChoice = MediaCommon.InputChoice.ph;
    //Data storage definitions
    //------------------------
    constant FreeFluids.MediaCommon.DataRecord data;
    //Auxiliary functions based in correlations
    //-----------------------------------------
    constant SpecificEnthalpy reference_h = SpecificEnthalpyCorr(reference_T);
    constant SpecificEntropy reference_s = SpecificEntropyCorr(reference_T);

    function PhysPropCorr "Calculates a physical property from a given correlation data and a given temperature, or pressure, using an external function."
      input Integer corr;
      input Real coef[:];
      input Real x;
      output Real y;
    
      external "C" FF_PhysPropCorr(corr, coef, data.MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end PhysPropCorr;

    function SpecificEnthalpyCorr "Calculates specific enthalpy from a given Cp correlation at a given temperature from 0K"
      input Real x;
      output Real y;
    
      external "C" FF_SpecificEnthalpyCorr(data.lCpCorr, data.lCpCoef, data.MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
      annotation(
        Inline = true,
        smoothOrder = 2);
    end SpecificEnthalpyCorr;

    function SpecificEntropyCorr "Calculates specific entropy from a given Cp correlation at a given temperature from 0K"
      input Real x;
      output Real y;
    
      external "C" FF_SpecificEntropyCorr(data.lCpCorr, data.lCpCoef, data.MW, x, y) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFphysprop.c\"");
    end SpecificEntropyCorr;

    //Function SpecificEnthalpyCorrInv

    function SpecificEnthalpyCorrInv "Compute temperature from property value"
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
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
      T := Internal.solve(y, 50.0, data.Tc, 1.0e5, {1}, fd);
    end SpecificEnthalpyCorrInv;

    //Function SpecificEntropyCorrInv

    function SpecificEntropyCorrInv "Compute temperature from property value"
      //input Correlation corrData;
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
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
      T := Internal.solve(y, 50.0, data.Tc, 1.0e5, {1}, fd);
    end SpecificEntropyCorrInv;

    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;

      algorithm
        assert(T <= data.Tc, "The media can´t be used over Tc");
        p := PhysPropCorr(data.VpCorr, data.VpCoef, T);
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
      T := Internal.solve(y, 50.0, data.Tc-0.01, 1.0e5, {1}, fd);
    end saturationPressureInv;

    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;

      algorithm
        if p<data.Pc then
          T := if data.BtCorr > 0 then PhysPropCorr(data.BtCorr, data.BtCoef, p) else saturationPressureInv(p);
        else
          T:=data.Tc;
        end if;
    end saturationTemperature;

    //BaseProperties model
    //--------------------

    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default))
        parameter MediaCommon.InputChoice localInputChoice = inputChoice;

      equation
        MM = data.MW / 1000 "in kg/mol";
        R = Modelica.Constants.R / MM;
        u = h - p / d;
        sat = setSat_p(p);
        if localInputChoice == MediaCommon.InputChoice.ph then
          state = setState_phX(p, h, fill(0, 0));
          state.T = T;
          state.d = d;
        elseif localInputChoice == MediaCommon.InputChoice.pT then
          state = setState_pTX(p, T, fill(0, 0));
          state.d = d;
          state.h = h;
        elseif localInputChoice == MediaCommon.InputChoice.dT then
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
        Density rho0;

      algorithm
        assert(T <= data.Tc, "The media can´t be used over Tc");
        state.p := p;
        state.T := T;
        Vp := saturationPressure(T);
        if p > Vp then
          rho0 := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T);
          state.d := if data.lnuA>0 then rho0 + log(data.lnuA * data.MW / (1000 * exp(data.lnuA * rho0 + data.lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / data.lnuA else rho0/ (1 - data.IsothComp*(state.p - Vp));
          state.h := SpecificEnthalpyCorr(T) - reference_h;
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
        Density Dl;
        Density Dg;
      algorithm
        assert(T <= data.Tc, "The media can´t be used over Tc");
        state.T := T;
        state.d := d;
        Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, T) "liquid density at boiling T";
        if d >= Dl then
          state.p:= if data.lnuA>0 then PhysPropCorr(data.VpCorr, data.VpCoef, T) + (exp((d-Dl)*data.lnuA)-1)*1000*Modelica.Constants.R*T*exp(data.lnuA*Dl+data.lnuB)/(data.lnuA*data.MW) else (1-Dl/d)/data.IsothComp+PhysPropCorr(data.VpCorr, data.VpCoef, T);
          state.h := SpecificEnthalpyCorr(T) - reference_h;
          state.phase := 1;
          state.gf := 0;
        else
          Dg := PhysPropCorr(data.gSatDensCorr, data.gSatDensCoef, T);
          assert(d >= Dg, "The media can´t work with non-saturated gases");
          state.p := PhysPropCorr(data.VpCorr, data.VpCoef, T);
          if d > Dg then
            state.phase := 2;
            state.gf := Dg * (Dl - d) / (d * (Dl - Dg));
            state.h := SpecificEnthalpyCorr(T) - reference_h + state.gf * PhysPropCorr(data.HvCorr, data.HvCoef, T);
          elseif d == Dg then
            state.phase := 1;
            state.gf := 1;
            state.h := SpecificEnthalpyCorr(T) - reference_h + PhysPropCorr(data.HvCorr, data.HvCoef, T);
          end if;
        end if;
    end setState_dTX;

    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,H and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEnthalpy Hl;
        SpecificEnthalpy Hg;
        Density Dl;
        Density Dg;

      algorithm
        state.p := p;
        state.h := h;
        if p > data.Pc then
          Tb := data.Tc;
        else
          Tb := saturationTemperature(p);
        end if;
        Hl := SpecificEnthalpyCorr(Tb) - reference_h "liquid specific enthalpy at boiling T";
        if Hl >= h then
          state.T := if data.lTfromHsatCorr > 0 then PhysPropCorr(data.lTfromHsatCorr, data.lTfromHsatCoef, h + reference_h) else SpecificEnthalpyCorrInv(h + reference_h);
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T) "saturated liquid density";
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T);
          state.d := if data.lnuA>0 then Dl + log(data.lnuA * data.MW / (1000 * exp(data.lnuA * Dl + data.lnuB) * Modelica.Constants.R * state.T) * (state.p - PhysPropCorr(data.VpCorr, data.VpCoef, state.T)) + 1) / data.lnuA else Dl/ (1 - data.IsothComp*(state.p - PhysPropCorr(data.VpCorr, data.VpCoef, state.T))) "saturated liquid density, pressure corrected";
          state.phase := 1;
          state.gf := 0;
        else
          Hg := Hl + PhysPropCorr(data.HvCorr, data.HvCoef, Tb);
          assert(Hg >= h, "The media can´t work with non-saturated gases");
          Dg := PhysPropCorr(data.gSatDensCorr, data.gSatDensCoef, Tb);
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, Tb);
          state.T := Tb;
          state.phase := 2;
          state.gf := (h - Hl) / (Hg - Hl);
          state.d := Dl * Dg / (state.gf * Dl + (1 - state.gf) * Dg);
        end if;
    end setState_phX;

    redeclare function extends setState_psX "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;

      protected
        Temperature Tb;
        SpecificEntropy Sl;
        SpecificEntropy Sg;
        Density Dl;
        Density Dg;

      algorithm
        state.p := p;
        if p > data.Pc then
          Tb := data.Tc;
        else
          Tb := saturationTemperature(p);
        end if;
        Sl := SpecificEntropyCorr(Tb) - reference_s "liquid specific entropy at boiling T";
        if Sl >= s then
          state.T := SpecificEntropyCorrInv(s + reference_s);
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T) "saturated liquid density";
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T);
          state.d := if data.lnuA>0 then Dl + log(data.lnuA * data.MW / (1000 * exp(data.lnuA * Dl + data.lnuB) * Modelica.Constants.R * state.T) * (state.p - PhysPropCorr(data.VpCorr, data.VpCoef, state.T)) + 1) / data.lnuA else Dl/ (1 - data.IsothComp*(state.p - PhysPropCorr(data.VpCorr, data.VpCoef, state.T))) "saturated liquid density, pressure corrected";
          state.h := SpecificEnthalpyCorr(state.T) - reference_h;
          state.phase := 1;
          state.gf := 0;
        else
          Sg := Sl + PhysPropCorr(data.HvCorr, data.HvCoef, Tb) / Tb;
          assert(Sg >= s, "The media can´t work with non-saturated gases");
          Dg := PhysPropCorr(data.gSatDensCorr, data.gSatDensCoef, Tb);
          Dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, Tb);
          state.T := Tb;
          state.phase := 2;
          state.gf := (s - Sl) / (Sg - Sl);
          state.h := SpecificEnthalpyCorr(state.T) - reference_h + state.gf * PhysPropCorr(data.HvCorr, data.HvCoef, Tb) / Tb;
          state.d := Dl * Dg / (state.gf * Dl + (1 - state.gf) * Dg);
        end if;
    end setState_psX;

    //Special thermodynamic states constructors
    //------------------------------------------

    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;
      protected
        AbsolutePressure Vp;
        Density rho0;
      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.h := SpecificEnthalpyCorr(sat.Tsat) - reference_h;
        Vp := saturationPressure(state.T);
        rho0 := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T);
        state.d := rho0 + log(data.lnuA * data.MW / (1000 * exp(data.lnuA * rho0 + data.lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / data.lnuA;
        state.gf := 0.0 "liquid";
        state.phase := 1;
    end setBubbleState;

    redeclare function extends setDewState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the dew point"
        extends Modelica.Icons.Function;

      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.h := SpecificEnthalpyCorr(sat.Tsat) - reference_h + PhysPropCorr(data.HvCorr, data.HvCoef, sat.Tsat);
        state.d := PhysPropCorr(data.gSatDensCorr, data.gSatDensCoef, sat.Tsat);
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
      MM := data.MW * 1e-3;
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

    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;

      protected
        Real vp1, vp2, dv;

      algorithm
        if state.gf == 0.0 then
          cp := PhysPropCorr(data.lCpCorr, data.lCpCoef, state.T);
        else
          cp := 0;
        end if;
//vp2:= PhysPropCorr(data.VpCorr, data.VpCoef, state.T);
//vp1:= PhysPropCorr(data.VpCorr, data.VpCoef, state.T-0.1);
//dv:=(vp2-vp1)/0.1;
//cp:=cp-isobaricExpansionCoefficient(state)*dv;
//cp:=cp-3e-7*(state.p-vp1)*state.T/data.MW;
//-0.00029936(P-Vp)*1e-6*T*1000/MW //J/kg·K
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

    redeclare function extends specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        h := state.h;
    end specificEnthalpy;

    redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf == 0.0 then
          s := SpecificEntropyCorr(state.T) - reference_s;
        else
          s := SpecificEntropyCorr(state.T) - reference_s + state.gf * PhysPropCorr(data.HvCorr, data.HvCoef, state.T) / state.T;
        end if;
    end specificEntropy;

    redeclare function extends isobaricExpansionCoefficient "Returns the approximate isobaric expansion coefficient beta at saturation pressure, neglecting the influence of the increment in Vp"
        extends Modelica.Icons.Function;

      protected
        Real dd, d0, d1, d2, dv, vp1, vp2;

      algorithm
        if state.gf == 0 then
          d1 := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T - 0.1);
          d2 := PhysPropCorr(data.lDensCorr, data.lDensCoef, state.T + 0.1);
          d0 := (d1 + d2) / 2;
          beta := (d1 - d2) * d0 / (d1 * d2 * 0.2);
          vp1 := PhysPropCorr(data.VpCorr, data.VpCoef, state.T - 0.1);
          vp2 := PhysPropCorr(data.VpCorr, data.VpCoef, state.T + 0.1);
          dv := (vp2 - vp1) / 0.2;
          beta := beta + isothermalCompressibility(state) * dv;
          if data.lnuA>0 then
            beta := beta * d0 * exp(-data.lnuA * (state.d - d0)) / state.d + (1 - exp(-data.lnuA * (state.d - d0))) / (data.lnuA * state.d * state.T);
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
          if data.lBulkModRCorr > 0 then
            kappa := data.MW / (exp(data.lBulkModRCoef[1] + data.lBulkModRCoef[2] * state.d +data.lBulkModRCoef[3] * state.d * state.d +data.lBulkModRCoef[4]* state.d * state.d * state.d +data.lBulkModRCoef[5]*state.d * state.d * state.d * state.d) * Modelica.Constants.R * state.T * state.d * 1000);
          elseif data.lnuA>0.0 then
            kappa := data.MW / (exp(data.lnuA * state.d + data.lnuB) * Modelica.Constants.R * state.T * state.d * 1000);
          elseif data.Vc>0.0 then
            rhoRed1:=state.d*1000*data.Vc/data.MW-1;
            kappa := data.MW / ((exp(-0.42704*rhoRed1+2.089*rhoRed1*rhoRed1 + 0.42367* rhoRed1* rhoRed1* rhoRed1)-1)* Modelica.Constants.R * state.T * state.d * 1000);
          else
            kappa := data.IsothComp;
          end if;
        else
          kappa := 0.0;
        end if;
    end isothermalCompressibility;

    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;

      protected
        Real cp, cv, kappa;

      algorithm
        cp := specificHeatCapacityCp(state);
        cv := specificHeatCapacityCv(state);
        kappa := isothermalCompressibility(state);
        a := sqrt(max(0, cp / (cv * kappa * state.d)));
      annotation(
        Inline = true,
        smoothOrder = 2);
    end velocityOfSound;

    redeclare function extends dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf == 0.0 then
          if data.lViscCorr>0 then
            eta := PhysPropCorr(data.lViscCorr, data.lViscCoef, state.T);
          else
            assert(false,"there is no correlation for liquid viscosity calculation");
          end if;
        elseif state.phase == 1.0 then
          eta := if data.gViscCorr > 0 then PhysPropCorr(data.gViscCorr, data.gViscCoef, state.T) else FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(data, state.T);
        else
          eta := 0;
        end if;
    end dynamicViscosity;

    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf == 0.0 then
          lambda := if data.lThCondCorr > 0 then PhysPropCorr(data.lThCondCorr, data.lThCondCoef, state.T) else FreeFluids.MediaCommon.Functions.liqThCondLatini(data, state.T);
        elseif state.gf == 1.0 then
          lambda := if data.gThCondCorr > 0 then PhysPropCorr(data.gThCondCorr, data.gThCondCoef, state.T) else FreeFluids.MediaCommon.Functions.gasThCondLowPressureChung(data, specificHeatCapacityCp(state), dynamicViscosity(state), state.T);
        else
          lambda := 0;
        end if;
    end thermalConductivity;

    redeclare function extends surfaceTension "Return surface tension"
        extends Modelica.Icons.Function;

      algorithm
        sigma := PhysPropCorr(data.lSurfTensCorr, data.lSurfTensCoef, sat.Tsat);
    end surfaceTension;

    //Saturation properties
    //-------------------------

    redeclare function extends bubbleDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;

      algorithm
        dl := PhysPropCorr(data.lDensCorr, data.lDensCoef, sat.Tsat);
    end bubbleDensity;

    redeclare function extends dewDensity "Return bubble point density from a sat record"
        extends Modelica.Icons.Function;

      algorithm
        dv := PhysPropCorr(data.gSatDensCorr, data.gSatDensCoef, sat.Tsat);
    end dewDensity;

    redeclare function extends bubbleEnthalpy "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        hl := SpecificEnthalpyCorr(sat.Tsat) - reference_h;
    end bubbleEnthalpy;

    redeclare function extends dewEnthalpy "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;

      algorithm
        hv := SpecificEnthalpyCorr(sat.Tsat) - reference_h + PhysPropCorr(data.HvCorr, data.HvCoef, sat.Tsat);
    end dewEnthalpy;

    redeclare function extends bubbleEntropy "Return bubble point specific entropy"
        extends Modelica.Icons.Function;

      algorithm
        sl := SpecificEntropyCorr(sat.Tsat) - reference_s;
    end bubbleEntropy;

    redeclare function extends dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;

      algorithm
        sv := SpecificEntropyCorr(sat.Tsat) - reference_s + PhysPropCorr(data.HvCorr, data.HvCoef, sat.Tsat) / sat.Tsat;
    end dewEntropy;
  end TMedium;

  //***** Substance Packages*****
  //*****************************

  package Acetone
    extends TMedium(final mediumName = "Acetone", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Acetone);
  end Acetone;

  package Ammonia
    extends TMedium(final mediumName = "Ammonia", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Ammonia,reference_T=273.15);
  end Ammonia;

  package n_Butanol
    extends TMedium(final mediumName = "n-Butanol", final singleState = false, data = FreeFluids.MediaCommon.MediaData.n_Butanol);
  end n_Butanol;

  package CO2
    extends TMedium(final mediumName = "Carbon dioxide", final singleState = false, data = FreeFluids.MediaCommon.MediaData.CO2);
  end CO2;

  package Dichlorodifluormethane
    extends TMedium(final mediumName = "Dichlorodifluormethane", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Dichlorodifluormethane,reference_T=273.15);
  end Dichlorodifluormethane;

  package EG
    extends TMedium(final mediumName = "Ethylene glycol", final singleState = false, data = FreeFluids.MediaCommon.MediaData.EG);
  end EG;

  package Ethanol
    extends TMedium(final mediumName = "Ethanol", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Ethanol,reference_T=273.15);
  end Ethanol;

  package Isobutane
    extends TMedium(final mediumName = "Isobutane", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Isobutane,reference_T=273.15);
  end Isobutane;

  package MarlothermSH
    extends TMedium(final mediumName = "Marlotherm SH", final singleState = false, p_default = 1.0e5, T_default = 425.0, reference_T=273.15, data = FreeFluids.MediaCommon.MediaData.MarlothermSH);
  end MarlothermSH;

  package N2
    extends TMedium(final mediumName = "Nitrogen", final singleState = false, data = FreeFluids.MediaCommon.MediaData.N2);
  end N2;

  package O2
    extends TMedium(final mediumName = "Oxygen", final singleState = false, data = FreeFluids.MediaCommon.MediaData.O2);
  end O2;

  package ShellS2
    extends TMedium(final mediumName = "Shell S2", final singleState = false, p_default = 2.0e5, T_default = 425.0, reference_T=273.15,data = FreeFluids.MediaCommon.MediaData.ShellS2);
  end ShellS2;

  package Toluene
    extends TMedium(final mediumName = "Toluene", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Toluene, reference_T=273.15);
  end Toluene;

  package Water
    extends TMedium(final mediumName = "Water saturated liquid", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Water, reference_T = 273.15);
  end Water;


  package Tests
    partial model FluidTesting
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.AbsolutePressure p = 1.0e5;
      parameter Medium.Temperature initialT = 273.19;
      parameter Medium.Temperature finalT = 372.15;
      parameter Real fract=0.1 "fraction of the density of stateP to be used in state Dlow";
      Medium.Temperature T(start = initialT) "We will ramp the temperature";
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
      StateDlow := Medium.setState_dTX(D*fract, T, fill(0, 0));
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
      BaseProp.p=p;
      BaseProp.h=H;
      der(T) = finalT - initialT;
    end FluidTesting;

    model Test1A
      extends FluidTesting(redeclare replaceable package Medium = TMedia.Water, p = 50.0e5, initialT = 0.1+273.15, finalT = 220+273.15);
    end Test1A;

    model Test1B
      extends Test1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model");
    end Test1B;

    model Test2A
      extends FluidTesting(redeclare replaceable package Medium = TMedia.O2, p = 50.0e5, initialT = 55.0, finalT = 155.0);
    end Test2A;
  end Tests;
  annotation(
    Documentation(info = "<html>
    <body>
    <p>The medium is designed for liquid, and liquid/vapor phases, at a temperature lower than 0.85 Tc, and a pressure not higher than 200 bars, because in those situations the liquid heat capacity coefficient is highly dependent on the pressure. It extends the Modelica  PartialTwoPhaseMedium. The medium properties are obtained using correlations that are mainly functions of T. It is somewhat similar to the Modelica TableBased medium, but uses specific correlations for each physical property, allows to work with saturated gas phase, and adds a density dependent correlation for the reduced bulk modulus of the liquid, that improves a lot the calculation of liquid density at high pressure, isothermal compressiblity, and isobaric expansion coefficient. Improving also the calculation of heat capacity at constant volume (Cv) and the speed of sound. It makes no use of the fluidConstants record, instead it defines a new record called data for storage of the substance data. This makes that the function molarMass doesn't work, as it takes the data from the fluidConstants record, and it is not replaceable.</p>
    <p>The obtained properties correspond to the liquid phase, or to the saturated gas phase. It can´t be used for non-saturated gas phase.</p>
    <p>If available, the reduced bulk modulus correlation for the liquid is used. In other, case a substance specific isothermal compressibility factor (with a default value of 6.667 e-10) is used. The parameters for the reduced bulk modulus correlation are normally not available, but can be calculated from a good equation of state of the multiparameter or SAFT types. This can be done easily with the FreeFluids GUI: you make the calculation with the EOS, transfer the results (density and the natural logarithm of the reduced bulk modulus) to the correlations tab, make the regresion of the coefficients, and store the result in the database. It is good to calculate the reduced bulk modulus at a pressure close to 50 bars but, if necessary in order to have liquid phase at the temperature of interest, it can be done at higher pressure. Check that all density data correspond to the liquid state.</p>
    <p>The liquid heat capacity correlation can be also a problem, as many times we only find it with a temperature limit of the normal boiling point. This can be solved using a Cp correlation constructed from a good EOS, using the FreeFluids GUI. It is important not to use data too close to the Tc for the regression (use data till 20 K below the Tc). The best equation for the regression of the liquid heat capacity is the PPDS15 equation. Do not use the Chemsep equations as they are not integrated by the medium to obtain enthalpy or entropy</p>
    <p> The liquid enthalpy is always calculated from the liquid Cp correlation as saturated. From the saturated specific enthalpy we can obtain an approximate value of T. It is inexact, but 40 bars of pressure difference from saturation point can generate in a liquid an error close to only 1 degree.</p>
    <p></p>The values of enthalpy and entropy are referenced to a liquid at reference_T.
    <p>The use of the liquid Cp correlation is an alternative to the use of ideal Cp correlation plus vaporization enthalpy. This makes possible the use of the medium with substances for which we do not have Cp0 data.</p> 
    <p>When going back from liquid enthalpy, or entropy, to temperature, we have two ways: One of them is to fit a correlation that makes this calculation, that can again be made with FreeFluidsGui. The second, that will be used if we left the value of lTfromHsatCorr=0, will calculate insitu the temperature using a solving function, as is done in the IdealGasMedia package.</p>
    <p>For liquid transport properties, the calculation is done directly from T, but pressure corrections could be introduced. For transport properties of gases the calculation is from T at low pressure. Here also, pressure correction could be introduced.</p> 
    <p>The thermodynamic record contains: p,T,gas fraction, d and h. Care must be taken in limiting the use to the temperature limits of the correlations used, as only few checks are done by the media, in order not to interfer with the solver process</p>
    <p>A new enumeration type called InputChoice has been defined in order to specify the independent variables to use in each instance of the BaseProperties model. This is done assigning value to the parameter localInputChoice in the model, that has as default the value given to the constant InputChoice inputChoice at package level. </p>
    <p>In the package Tests there is a comparison between the medium performance with water and the WaterIF97_ph medium model.</p>
    <p>The global idea has been not to use the Modelica files for the storage of substances data, but to store the data in a database, from which we can recover and use them when needed. A database is provided with more than 400 substances that can be enlarged, and the FreeFluids GUI can retrieve the data from the data base, treat it as needed (for example creating EOS from saturated vapor pressure and/or densities, or creating correlations from the EOS), store the results in the database, and export the data in Modelica format when needed.</p>
    <p>As a resume: The medium is for fast calculation of liquid phase, condensation, and evaporation, but not for non-saturated gases.</p>
    </body>
    </html>"));
end TMedia;

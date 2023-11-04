within FreeFluids;

package ExternalPure "ExternalPure.mo by Carlos Trujillo
  This file is part of the Free Fluids application
  Copyright (C) 2008-2023  Carlos Trujillo Gonzalez
    
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
  //***** PACKAGE ExternalMedium*****
  //*********************************

  package ExternalMedium
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(onePhase = false, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, reference_T = 298.15, reference_p = 101325.0, fluidConstants = {fluidK});
  
    redeclare record FluidConstants
      extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
      extends Modelica.Icons.Record;
      String description = " " "Description of the EOS and correlations used for the substance";
      Real criticalTemperature = 0.0 "critical temperature in K";
      Real criticalPressure = 0.0 "critical pressure in Pa";
      //criticalPressure must be added for compatibility with ThermoPower
    end FluidConstants;
  
  
  
    constant String refName = "Propane" "name of the reference fluid for ECS calculations. Defaults to propane";
    constant String resDir = Modelica.Utilities.Files.loadResource("modelica://FreeFluids/Resources") "resources directory";
    constant FluidConstants fluidK "just for minimum data out of the SubstanceData external object";
    constant String description = " " "Description of the EOS and correlations used for the substance";
    constant Integer refState = 1 "1=ASHRAE, 2=IIR, 3=NBP, >3=User. When user, the value at reference_T and reference_p as 0";
    constant String inputChoice = "ph" "Allows to choose the input variables to use for the construction of a BaseProperties object. alternatives: ph,pT,dT";
    constant Integer thermoModel = 1 "1 for cubic EOS; 2 for SAFT; 3 Saltzmann and Wagner. All transport properties will be computed by correlations";
    constant Integer ancillaries = 0 "if 0 no ancillary function will be used, if 1 they will be use for cubic EOS, if 2 for SAFT, and if 3 for SW. It indicates for which EOS are the functions fitted";
    
    constant Real T_min=150.0;
    constant Real TCri=fluidK.criticalTemperature;
    constant Real pCri=fluidK.criticalPressure;
    
    //Definition of ThermodynamicState, BaseProperties and SubstanceData
    //------------------------------------------------------------------
  
    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        //Real nMols "number of moles in a kg";
        AbsolutePressure p(displayUnit = "bar") "Pressure in Pa";
        Temperature T(displayUnit = "degC") "Kelvin temperature of medium";
        SpecificEnthalpy h "Overall specific enthalpy";
        SpecificEntropy s "Overall specific entropy";
        Density d(displayUnit = "kg/m3") "Overall density in kgr/m3";
        MoleFraction gf "gas molar fraction";
        //additional variables
        Density ld(displayUnit = "kg/m3") "Liquid phase density in kgr/m3";
        Density gd(displayUnit = "kg/m3") "Gas phase density in kgr/m3";
        SpecificEnthalpy lh "liquid phase specific enthalpy";
        SpecificEnthalpy gh "gas phase specific enthalpy";
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
  
    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates and localInputChoice == "ph" then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates and (localInputChoice == "ph" or localInputChoice == "pT") then StateSelect.prefer else StateSelect.default), T(stateSelect = if preferredMediumStates and (localInputChoice == "pT" or localInputChoice == "dT") then StateSelect.prefer else StateSelect.default), d(stateSelect = if preferredMediumStates and localInputChoice == "dT" then StateSelect.prefer else StateSelect.default))
        constant String localInputChoice = inputChoice;
        Integer phase(min = 0, max = 2, start = 1, fixed = false);
  
      algorithm
        MM := fluidK.molarMass "in kg/mol";
        R_s := Modelica.Constants.R/MM;
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
        u := h - p/d;
        sat := setSat_p(p);
    end BaseProperties;
  
    //Basic functions
    //---------------
  
    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;
  
  
        external "C" FF_saturationPressureM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, T, p) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end saturationPressure;
  
    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;
  
  
        external "C" FF_saturationTemperatureM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, p, T) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end saturationTemperature;
  
    redeclare function extends saturationTemperature_derp "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        dTp := if p < fluidConstants[1].criticalPressure then (saturationTemperature(p) - saturationTemperature(0.999*p))/(0.001*p) else 0;
    end saturationTemperature_derp;
  
    function eosCriticalConstants
      output Temperature Tc;
      output AbsolutePressure Pc;
    
      external "C" FF_eosCriticalConstantsM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, Tc, Pc) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end eosCriticalConstants;
  
    function pressureEOS_dT "Return pressure given temperature and density, by EOS"
      input Density d;
      input Temperature T;
      output AbsolutePressure p;
    
      external "C" FF_pressureEOS_dT(mediumName, resDir, thermoModel, refState, reference_T, reference_p, d, T, p) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end pressureEOS_dT;
  
    function densities_pT "Return liquid and gas densities at give temperature and pressure, by EOS, bellow the critical temperature"
      input AbsolutePressure p;
      input Temperature T;
      input String var;
      //input String var "b"=both phases, "g"= only gas, "l"=only liquid";
      output Density ld;
      output Density gd;
    
      external "C" FF_densities_pTM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, var, p, T, ld, gd) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end densities_pT;
  
    function dhsFrom_pT "Return liquid and gas densities at give temperature and pressure, by EOS, bellow the critical temperature"
      input AbsolutePressure p;
      input Temperature T;
      input String var;
      //"e"=bubble, "d"=dew, "b"=both phases, "l"= liquid phase, "g"=gas phase
      output Density ld;
      output SpecificEnthalpy lh;
      output SpecificEntropy ls;
      output Density gd;
      output SpecificEnthalpy gh;
      output SpecificEntropy gs;
    
      external "C" FF_dhsFrom_pTM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, var, p, T, ld, lh, ls, gd, gh, gs) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end dhsFrom_pT;
  
    function solveEOS "Calls the external function that calculates some basic thermodynamic properties of the phases, used later for the definition of the ThermodynamicState and the calculation of  all thermodynamic properties"
      input String opt;
      input Temperature x;
      input AbsolutePressure y;
      output Real T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT, gf;
    
      external "C" FF_solveEosM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, opt, x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT, gf) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS;
  
    //Establish general states
    //------------------------
  
    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi. The function from T and P is unable to compute any gas fraction different from 0 or 1"
        extends Modelica.Icons.Function;
  
      algorithm
        state.phase := 1;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("p", T, p);
        if state.ld < 1.0 then
  //state.gf := 1.0;
          state.d := state.gd;
          state.h := state.gh;
          state.s := state.gs;
        else
  //state.gf := 0.0;
          state.d := state.ld;
          state.h := state.lh;
          state.s := state.ls;
        end if;
    end setState_pTX;
  
    redeclare function extends setState_dTX "Return thermodynamic state as function of density, T and composition X or Xi"
        extends Modelica.Icons.Function;
  
      algorithm
        state.d := d;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("d", T, d);
  //state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else state.gd * (state.ld - d) / (d * (state.ld - state.gd));
        state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
        state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        state.phase := if state.gf == 0 then 1 else if state.gf == 1 then 1 else 2;
    end setState_dTX;
  
    //here is necessary to arrive to calculate the thermoprops, so better to do every thing inside C
  
    redeclare function extends setState_phX "Return thermodynamic state as function of pressure, enthalpy and composition X or Xi"
        extends Modelica.Icons.Function;
  
      algorithm
        state.h := h;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("h", p, h);
  //state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else (state.h - state.lh) / (state.gh - state.lh);
        state.d := if state.gf == 1.0 then state.gd else if state.gf == 0.0 then state.ld else state.ld*state.gd/(state.gf*state.ld + (1.0 - state.gf)*state.gd);
        state.s := state.ls*(1 - state.gf) + state.gs*state.gf;
        state.phase := if state.gf == 0 then 1 else if state.gf == 1 then 1 else 2;
    end setState_phX;
  
    redeclare function extends setState_psX "Return thermodynamic state as function of pressure, entropy and composition X or Xi"
        extends Modelica.Icons.Function;
  
      algorithm
        state.s := s;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("s", p, s);
  //state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else (state.s - state.ls) / (state.gs - state.ls);
        state.d := if state.gf > 0.999999 then state.gd else if state.gf < 0.000001 then state.ld else state.ld*state.gd/(state.gf*state.ld + (1.0 - state.gf)*state.gd);
        state.h := state.lh*(1 - state.gf) + state.gh*state.gf;
        state.phase := if state.gf == 0 then 1 else if state.gf == 1 then 1 else 2;
    end setState_psX;
  
    //Establish special states
    //------------------------
  
    redeclare function extends setDewState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the dew point"
        extends Modelica.Icons.Function;
  
      algorithm
        state.phase := 2;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("w", sat.Tsat, sat.psat);
        state.d := state.gd;
        state.h := state.gh;
        state.s := state.gs;
    end setDewState;
  
    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;
  
      algorithm
        state.phase := 2;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT, state.gf) := solveEOS("e", sat.Tsat, sat.psat);
        state.d := state.ld;
        state.h := state.lh;
        state.s := state.ls;
    end setBubbleState;
  
    redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
        import Modelica.Media.Common.smoothStep;
  
      algorithm
        state := ThermodynamicState(
        p = smoothStep(x, state_a.p, state_b.p, x_small),
        h = smoothStep(x, state_a.h, state_b.h, x_small), 
        s = smoothStep(x, state_a.s, state_b.s, x_small), 
        T = smoothStep(x, state_a.T, state_b.T, x_small), 
        d = smoothStep(x, state_a.d, state_b.d, x_small), 
        ld = smoothStep(x, state_a.ld, state_b.ld, x_small), 
        gd = smoothStep(x, state_a.gd, state_b.gd, x_small), 
        lh = smoothStep(x, state_a.lh, state_b.lh, x_small), 
        gh = smoothStep(x, state_a.gh, state_b.gh, x_small), 
        ls = smoothStep(x, state_a.ls, state_b.ls, x_small), 
        gs = smoothStep(x, state_a.gs, state_b.gs, x_small),
        lCv = smoothStep(x, state_a.lCv, state_b.lCv, x_small),
        gCv = smoothStep(x, state_a.gCv, state_b.gCv, x_small),
        lCp = smoothStep(x, state_a.lCp, state_b.lCp, x_small),
        gCp = smoothStep(x, state_a.gCp, state_b.gCp, x_small),
        lDvp = smoothStep(x, state_a.lDvp, state_b.lDvp, x_small),
        gDvp = smoothStep(x, state_a.gDvp, state_b.gDvp, x_small),
        lDvT = smoothStep(x, state_a.lDvT, state_b.lDvT, x_small),
        gDvT = smoothStep(x, state_a.gDvT, state_b.gDvT, x_small),
        gf = smoothStep(x, state_a.gf, state_b.gf, x_small), phase = 0);
      annotation(
        Inline = true);
    end setSmoothState;
  
    //Getting properties using the ThermodynamicState record
    //------------------------------------------------------------------------------------------------------------
  
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
        if state.gf == 0 then
          ddpT := state.ld*isothermalCompressibility(state);
        elseif state.gf == 1 then
          ddpT := state.gd*isothermalCompressibility(state);
        else
          ddpT := 0;
        end if;
    end density_derp_T;
  
    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.gf == 0 then
          ddTp := -state.ld*isobaricExpansionCoefficient(state);
        elseif state.gf == 1 then
          ddTp := -state.gd*isobaricExpansionCoefficient(state);
        else
          ddTp := 0;
        end if;
    end density_derT_p;
  
    redeclare function extends density_derp_h "Return density derivative w.r.t. pressure at const specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        Real lDvp, gDvp, DTp, lDhp, gDhp "total derivatives of v, T, and h along the saturation line, regarding p";
        Real Dgfp_h "partial derivative of gas fraction regarding p at constant h";
  
      algorithm
        if state.gf == 0 then
          ddph := -state.ld*state.ld*(state.lDvp + (state.T*state.lDvT*state.lDvT - state.lDvT/state.ld)/state.lCp);
        elseif state.gf == 1 then
          ddph := -state.gd*state.gd*(state.gDvp + (state.T*state.gDvT*state.gDvT - state.gDvT/state.gd)/state.gCp);
        else
          DTp := (1/state.gd - 1/state.ld)/(state.gs - state.ls);
          lDvp := state.lDvp + state.lDvT*DTp;
          gDvp := state.gDvp + state.gDvT*DTp;
          lDhp := 1/state.ld - state.T*state.lDvT + state.lCp*DTp;
          gDhp := 1/state.gd - state.T*state.gDvT + state.gCp*DTp;
          Dgfp_h := (state.gf*gDhp + (1 - state.gf)*lDhp)/(state.lh - state.gh);
          ddph := -state.d*state.d*(lDvp + Dgfp_h*(1/state.gd - 1/state.ld) + state.gf*(gDvp - lDvp));
        end if;
    end density_derp_h;
  
    redeclare function extends density_derh_p "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        if (state.gf == 0) then
          ddhp := -state.ld*state.ld*state.lDvT/state.lCp;
        elseif (state.gf == 1) then
          ddhp := -state.gd*state.gd*state.gDvT/state.gCp;
        else
          ddhp := -state.d*state.d*(1/state.gd - 1/state.ld)/(state.gh - state.lh);
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
        cp := if state.gf == 0 then state.lCp else if state.gf == 1 then state.gCp else Modelica.Constants.inf;
    end specificHeatCapacityCp;
  
    redeclare         function extends specificHeatCapacityCv "Return specific heat capacity at constant volume"
            extends Modelica.Icons.Function;
    
          protected
            Real lDsT, gDsT, dP_dT, ldV_dT, gdV_dT "total derivatives of s, p and v along the saturation line, regarding T";
            Real DgfT_v "partial derivative of gas fraction regarding T at constant v";
    
          algorithm
            if state.gf == 0 then
              cv := state.lCv;
            elseif state.gf == 1 then
              cv := state.gCv;
            else
  //dP_dT := (state.gs - state.ls) / (1 / state.gd - 1 / state.ld);
  //lDsT := state.lCp / state.T - state.lDvT * dP_dT;
  //gDsT := state.gCp / state.T - state.gDvT * dP_dT;
  //ldV_dT := state.lDvT + state.lDvp * dP_dT;
  //gdV_dT := state.gDvT + state.gDvp * dP_dT;
  //DgfT_v :=(state.gf * gdV_dT + (1 - state.gf) * ldV_dT) / (1 / state.ld - 1 / state.gd);
  //cv := state.T * lDsT + state.T * DgfT_v * (state.gs - state.ls) + state.gf * state.T * (gDsT - lDsT);
  //cv:=state.lCv+DgfT_v*((state.gh-state.p/state.gd)-(state.lh-state.p/state.ld))+state.gf*(state.gCv-state.lCv);
              cv := state.lCv*(1 - state.gf) + state.gCv*state.gf;
    
            end if;
        end specificHeatCapacityCv;
  
    redeclare function extends isentropicExponent "Return isentropic exponent"
        extends Modelica.Icons.Function;
  
      protected
        Real Vinv, cv, cp;
  
      algorithm
        gamma := if state.gf == 0 then state.lCp/state.lCv else if state.gf == 1 then state.gCp/state.gCv else 0.0;
    end isentropicExponent;
  
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
  
    function jouleThomsonCoefficient
      input ThermodynamicState state;
      output Real J;
    algorithm
      if state.phase == 1 then
        J := (state.T*isobaricExpansionCoefficient(state) - 1)/(state.d*specificHeatCapacityCp(state));
      else
        J := Modelica.Constants.inf;
      end if;
    end jouleThomsonCoefficient;
  
    function isothermalThrottlingCoefficient "OK. The ExternalMedia function is not found !!"
      input ThermodynamicState state;
      output Real I;
    algorithm
      if state.phase == 1 then
        I := -(state.T*isobaricExpansionCoefficient(state) - 1)/state.d;
      else
        I := Modelica.Constants.inf;
      end if;
    end isothermalThrottlingCoefficient;
  
    redeclare function extends dynamicViscosity "Return dynamic viscosity from a ThermodynamicState record"
        extends Modelica.Icons.Function;
        external "C" FF_ViscosityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, refName, state.T, state.d, state.p, state.gf, eta) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end dynamicViscosity;
  
    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;
        external "C" FF_thermalConductivityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, state.T, state.p, state.gf, lambda) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end thermalConductivity;
  
    function liquidDynamicViscosity "Return liquid dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    
      external "C" FF_dynamicViscosityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, state.T, state.p, 0.0, eta) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end liquidDynamicViscosity;
  
    function gasDynamicViscosity "Return gas dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    
      external "C" FF_dynamicViscosityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, state.T, state.p, 1.0, eta) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end gasDynamicViscosity;
  
    function liquidThermalConductivity "Return liquid thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    
      external "C" FF_thermalConductivityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, state.T, state.p, 0.0, lambda) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end liquidThermalConductivity;
  
    function gasThermalConductivity "Return gas thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    
      external "C" FF_thermalConductivityM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, state.T, state.p, 1.0, lambda) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end gasThermalConductivity;
  
    redeclare function extends surfaceTension "Return surface tension"
        extends Modelica.Icons.Function;
  
  
        external "C" FF_surfaceTensionM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, sat.Tsat, sigma) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end surfaceTension;
  
    //Functions at the Dew and Bubble points, from a SaturationProperties record linking SetSat functions are the originals
    //---------------------------------------------------------------------------------------------------------------------
  
    redeclare function extends bubbleDensity "Return bubble point density"
        extends Modelica.Icons.Function;
      protected
        Real dg;
  
      algorithm
        (dl, dg) := densities_pT(sat.psat, sat.Tsat, "e");
    end bubbleDensity;
  
    redeclare function extends dBubbleDensity_dPressure "Return bubble point density derivative"
        extends Modelica.Icons.Function;
  
      protected
        Real Tc, Pc;
        Real dl, dg, pl, dll, dgl;
        //pressure and liquid density at sat.Tsat-0.01
        Real dT_dp;
        //Along the saturation line
        ThermodynamicState stateB, stateD;
        //bubble and dewstate
  
      algorithm
        (Tc, Pc) := eosCriticalConstants();
        if (sat.Tsat >= Tc) then
          ddldp := 0;
        else
          if (0 == 0) then
            stateB := setBubbleState(sat);
            stateD := setDewState(sat);
            dT_dp := (1/stateD.d - 1/stateB.d)/(stateD.s - stateB.s);
            ddldp := -stateB.lDvp*stateB.ld^2 - stateB.lDvT*stateB.ld^2*dT_dp;
          else
  //numeric derivarive
            (dl, dg) := densities_pT(sat.psat, sat.Tsat, "l");
            pl := saturationPressure(sat.Tsat - 0.3);
            (dll, dgl) := densities_pT(pl, sat.Tsat - 0.3, "l");
            ddldp := if ((sat.psat - pl)/sat.psat > 0.000001) then (dl - dll)/(sat.psat - pl) else 0;
          end if;
          if (ddldp > 0) then
            ddldp := 0;
          end if;
        end if;
    end dBubbleDensity_dPressure;
  
    redeclare function extends dewDensity "Return dew point density"
        extends Modelica.Icons.Function;
  
      protected
        Real dl;
  
      algorithm
        (dl, dv) := densities_pT(sat.psat, sat.Tsat, "w");
    end dewDensity;
  
    redeclare function extends dDewDensity_dPressure "Return dew point density derivative"
        extends Modelica.Icons.Function;
  
      protected
        Real Tc, Pc;
        Real dl, dg, pl, dll, dgl;
        //pressure and liquid density at sat.Tsat-0.01
        Real dT_dp;
        //derivative of temperature regarding pressure along the saturation line
        ThermodynamicState stateB, stateD;
        //bubble and dewstate
  
      algorithm
        (Tc, Pc) := eosCriticalConstants();
        if (sat.Tsat >= Tc) then
          ddvdp := 0;
        else
          if (0 == 0) then
            stateB := setBubbleState(sat);
            stateD := setDewState(sat);
            dT_dp := (1/stateD.d - 1/stateB.d)/(stateD.s - stateB.s);
            ddvdp := -stateD.gDvp*stateD.gd^2 - stateD.gDvT*stateD.gd^2*dT_dp;
          else
  //numeric derivative
            (dl, dg) := densities_pT(sat.psat, sat.Tsat, "g");
            pl := saturationPressure(sat.Tsat - 1);
            (dll, dgl) := densities_pT(pl, sat.Tsat - 0.5, "g");
            ddvdp := if ((sat.psat - pl)/sat.psat > 0.001) then (dg - dgl)/(sat.psat - pl) else 0;
          end if;
        end if;
    end dDewDensity_dPressure;
  
    redeclare function extends bubbleEnthalpy "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        Real ld, ls, gd, gh, gs;
  
      algorithm
  //external "C" FF_bubbleEnthalpyM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, sat.psat, sat.Tsat, hl) annotation(
  //IncludeDirectory = "modelica://FreeFluids/Resources",
  //Include = "#include \"FFmodelicaMedium.c\"");
        (ld, hl, ls, gd, gh, gs) := dhsFrom_pT(sat.psat, sat.Tsat, "e");
    end bubbleEnthalpy;
  
    redeclare function extends dBubbleEnthalpy_dPressure "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;
  
      protected
        Real Tc, Pc;
        Real pl;
        Real dT_dp;
        //Along the saturation line
        ThermodynamicState stateB, stateD;
        //bubble and dewstate
  
      algorithm
        (Tc, Pc) := eosCriticalConstants();
        if (sat.Tsat >= Tc) then
          dhldp := 0;
        else
          if (0 == 0) then
  //symbolic calculation
            stateB := setBubbleState(sat);
            stateD := setDewState(sat);
            dT_dp := (1/stateD.d - 1/stateB.d)/(stateD.s - stateB.s);
            dhldp := -stateB.T*stateB.lDvT + 1/stateB.d + stateB.lCp*dT_dp;
          else
            pl := saturationPressure(sat.Tsat - 0.3);
            dhldp := (bubbleEnthalpy(sat) - bubbleEnthalpy(setSat_T(sat.Tsat - 0.3)))/(sat.psat - pl);
          end if;
          if (dhldp < 0) then
            dhldp := 0;
          end if;
        end if;
    end dBubbleEnthalpy_dPressure;
  
    redeclare function extends dewEnthalpy "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;
  
      protected
        Real ld, lh, ls, gd, gs;
  
      algorithm
        (ld, lh, ls, gd, hv, gs) := dhsFrom_pT(sat.psat, sat.Tsat, "w");
    end dewEnthalpy;
  
    redeclare function extends dDewEnthalpy_dPressure "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;
  
      protected
        Real Tc, Pc;
        Real pl;
        Real dT_dp;
        //Along the saturation line
        ThermodynamicState stateB, stateD;
        //bubble and dewstate
  
      algorithm
        (Tc, Pc) := eosCriticalConstants();
        if (sat.Tsat >= Tc) then
          dhvdp := 0;
        else
          if (0 == 0) then
            stateB := setBubbleState(sat);
            stateD := setDewState(sat);
            dT_dp := (1/stateD.d - 1/stateB.d)/(stateD.s - stateB.s);
            dhvdp := -stateD.T*stateD.gDvT + 1/stateD.d + stateD.gCp*dT_dp;
          else
  //numeric calculation
            pl := saturationPressure(sat.Tsat - 0.5);
            dhvdp := (dewEnthalpy(sat) - dewEnthalpy(setSat_T(sat.Tsat - 0.5)))/(sat.psat - pl);
          end if;
        end if;
  //algorithm
  //  dhvdp := if sat.Tsat < fluidConstants[1].criticalTemperature then (dewEnthalpy(setSat_T(sat.Tsat - 0.02)) - dewEnthalpy(sat)) / (saturationPressure(sat.Tsat - 0.02) - sat.psat) else 0;
    end dDewEnthalpy_dPressure;
  
    function vaporizationEnthalpy
      input SaturationProperties sat;
      output SpecificEnthalpy Hv;
    algorithm
      Hv := dewEnthalpy(sat) - bubbleEnthalpy(sat);
    end vaporizationEnthalpy;
  
    redeclare function extends bubbleEntropy "Return bubble point specific entropy"
        extends Modelica.Icons.Function;
  
      protected
        Real ld, lh, gd, gh, gs;
  
      algorithm
  //external "C" FF_bubbleEntropyM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, sat.psat, sat.Tsat, sl) annotation(
  //IncludeDirectory = "modelica://FreeFluids/Resources",
  //Include = "#include \"FFmodelicaMedium.c\"");
        (ld, lh, sl, gd, gh, gs) := dhsFrom_pT(sat.psat, sat.Tsat, "e");
    end bubbleEntropy;
  
    redeclare function extends dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;
  
      protected
        Real ld, lh, ls, gd, gh;
  
      algorithm
  //external "C" FF_dewEntropyM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, sat.psat, sat.Tsat, sv) annotation(
  //IncludeDirectory = "modelica://FreeFluids/Resources",
  //Include = "#include \"FFmodelicaMedium.c\"");
        (ld, lh, ls, gd, gh, sv) := dhsFrom_pT(sat.psat, sat.Tsat, "w");
    end dewEntropy;
  
  
  //Functions needed for compatibility with Buildings.Media.Refrigerants
    function dPressureVap_dSpecificVolume_Tv "Derivative of pressure with regards to specific volume, at constant T"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Real dpdv(final unit="Pa.kg/m3") "Derivative of pressure with regards to specific volume";
    protected
      ThermodynamicState state;
    algorithm
      state:=setState_dTX(1/v,T);
      dpdv:= if(state.gf==0) then 1/state.lDvp elseif(state.gf==1) then 1/state.gDvp else 0;
    end dPressureVap_dSpecificVolume_Tv;
  
    function dPressureVap_dTemperature_Tv "Derivative of pressure with regards to temperature, at constant v"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Real dpdT(final unit="Pa/K") "Derivative of pressure with regards to temperature";
    protected
      ThermodynamicState state;
    algorithm
      state:=setState_dTX(1/v,T);
      dpdT:= if(state.gf==0) then -state.lDvT/state.lDvp elseif(state.gf==1) then -state.gDvT/state.gDvp else 0;
    end dPressureVap_dTemperature_Tv;
    
    function enthalpySatLiq_T "Function that calculates the enthalpy of saturated liquid based on temperature"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      output Modelica.Units.SI.SpecificEnthalpy h;
    algorithm
      h:=bubbleEnthalpy(setSat_T(T));
    end enthalpySatLiq_T;
    
    function enthalpySatVap_T "Function that calculates the specific enthalpy of saturated vapor based on temperature"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      output Modelica.Units.SI.SpecificEnthalpy h;
    algorithm
      h:=dewEnthalpy(setSat_T(T));
    end enthalpySatVap_T;
    
    function isentropicExponentVap_Tv "Function that calculates the isentropic exponent of vapor based on temperature and specific volume"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Modelica.Units.SI.IsentropicExponent k;
    algorithm
      k:=isentropicExponent(setState_dTX(1/v,T));
    end isentropicExponentVap_Tv;
  
    function pressureSatVap_T "Function that calculates the pressure of saturated vapor based on temperature"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      output Modelica.Units.SI.AbsolutePressure p "Pressure of saturated refrigerant vapor";
        external "C" FF_saturationPressureM(mediumName, resDir, thermoModel, refState, reference_T, reference_p, ancillaries, T, p) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end pressureSatVap_T;
  
    function pressureVap_Tv "Function that calculates the pressure vapor based on temperature and specific volume"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Modelica.Units.SI.AbsolutePressure p "Pressure of refrigerant vapor";
    algorithm
      p:=pressure(setState_dTX(1/v,T));
    end pressureVap_Tv;
  
    function specificIsobaricHeatCapacityVap_Tv "Function that calculates the specific isobaric heat capacity of vapor based on temperature and specific volume"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Modelica.Units.SI.SpecificHeatCapacity cp "Specific isobaric heat capacity";
    algorithm
      cp:=specificHeatCapacityCp(setState_dTX(1/v,T));
    end specificIsobaricHeatCapacityVap_Tv;
  
    function specificIsochoricHeatCapacityVap_Tv "Function that calculates the specific isochoric heat capacity of vapor
    based on temperature and specific volume"
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      input Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
      output Modelica.Units.SI.SpecificHeatCapacity cv "Specific isochoric heat capacity";
    algorithm
      cv:=specificHeatCapacityCv(setState_dTX(1/v,T));
    end specificIsochoricHeatCapacityVap_Tv;
    
    function specificVolumeVap_pT
      "Function that calculates the specific volume of the vapor based on pressure and temperature"
      input Modelica.Units.SI.AbsolutePressure p "Pressure of refrigerant vapor";
      input Modelica.Units.SI.Temperature T "Temperature of refrigerant";
      output Modelica.Units.SI.SpecificVolume v "Specific volume of refrigerant";
    algorithm
      v:=1/density(setState_pTX(p,T));
    end specificVolumeVap_pT;  
  end ExternalMedium;






  annotation(
    Documentation(info = "<html><head></head><body>
    <p><b>Introduction</b></p><p>The medium is designed for pure or pseudo-pure substances in liquid, gas, or two phases. The thermodynamic calculations are performed using different equations of state(EOS), implemented in C language, that are called using external functions. The C code and the substances data are in the Resources folder.</p>
    <p>For each substance there is the possibility of choosing between three EOS types: multiparameter, PCSAFT, or cubic. The &nbsp;substance data is contained in a C structure, that is exported from the database using the FreeFluids GUI software.</p>
    <p>The quality of the results is very high when multiparameter (Seltzmann and Wagner, SW) EOS are available. If not, PCSAFT or different flavours of cubic EOS can be used.</p>
    <p>The transport properties are normally computed from temperature dependent correlations, with pressure correction. Nevertheless, for viscosity, the system will perform phase independent viscosity calculation (from temperature and density) for selected substances, using dedicated calculation or extended correspondent states with NIST correction polynomia, when available. In order to use the phase independent viscosity calculation the thermoModel parameter must be 3 (multiparameter EOS), as it is the only way to be sure that the supplied density is reliable. For thermal conductivity the same type of calculations are already programmed, but have been inactivated till finishing testing.</p>
    <p>The medium implements all the requirements of Modelica.Media.PartialTwoPhaseMedium. It is compatible with the new frontend of OpenModelica.&nbsp;</p><p>&nbsp;The main limitation of the medium, when using SW or PCSAFT EOS, is in the vecinity of the critical point, as no special technique, as for example splines, is used. The medium is much slower than the TMedia one, but is the price for the better precision and wider application.</p>
    <p><b>The C code</b></p><p>The C code is placed in the Resources folder.</p><p>The FFbasic.h file contains the basic definition of structures and enumerations used in the code.&nbsp;</p><p>The FFmodelicaMedium.c file contains the interface between Modelica and C. The FFeosPure.c and FFphysprop.c files contain the code that will perform the calculations in C.</p><p><b>Transport properties by dedicated calculation or ECS</b></p><p>The calculation of the substance viscosity is managed by the FF_Viscosity function in the FFphysprop.c file. If you have selected to work with the multiparameter EOS (thermoModel=3), that grants a good density calculation, it will check if there is a dedicated calculation from temperature and density (available only for few substances). If not, it will check if the correlation defined for gas viscosity calculation is the 112 (NIST coefficients for viscosity calculation using ECS) and that the data charged for the reference substance correspond to the needed one. If this is OK the ECS calculation will be applied for viscosity, otherwise temperature dependent correlations, with pressure correction, will be used.</p><p>If the EOS is of the cubic or PCSAFT type, correlations will be used if available. If not, the Lucas approximation will be used for gas. For liquid phase, if there is no correlation defined, the calculation will be done using ECS from the reference liquid defined.</p><p>The definition of the reference liquid is done at medium package level, giving value to the const String refName, which default value is 'Propane'.</p><p><b>Medium definition</b></p><p>We need to define the mediums to be used. The already made definitions are inside the Fluids subpackage. Each definition is very short, just with the minimum information necessary. The most important of it is the mediumName, as this name is the name of the file (with extension .sd) placed in the Resources/Fluids folder that will be used for charging the data to the C structure.</p><p></p><p>For each medium you need both a medium definition in Modelica, and a binary file with the data in the Resources/Fluids folder. For adding new media both can be made using the program FreeFluidsGui.exe that has been also placed in the Resources/Extra folder. It will access the data base, allow you to select the substance, with the EOS and correlations to be used, and later export it to a binary file, with the data as C structure, and to a text file. The data file (.sd extension) must be placed into the FreeFluids.ExternalPure.Fluids folder. The medium definition package in Modelica can be placed at your choice, but is recomended to be placed inside the ExternalPure/Fluids package.</p><p>It is important to understand that the Vp, liquid saturated density, and gas saturated density, to use are not the best one, but the ones that match the calculation of the EOS. For example, if you have a PCSAFT EOS you must first construct correlations that match this EOS, and export them along with the EOS as a C binary structure. In the package extension in Modelica you should make the 'ancillaries' Integer constant equal to 0,2 or 3 in order to: no use ancillaries, use with SAFT type EOS, or use with SW type EOS. Transport properties correlations should be the best availables.</p><p>The program allows also the exportation of the correlations in the format used by the TMedia package.</p><p><b>Medium configuration</b></p><p>When extending the medium, the following configuration must be done:</p>
    <p>The thermoModel Integer constant must be made equal to 1 for using the cubic EOS, 2 for the PCSAFT, or 3 for the SW.</p><p>The ancillaries Integer constant must be made equal to 3 for using ancillary functions (for saturated pressure and densities) with SW EOS, in order to obtain faster calculation. If made equal to 2 they will be used with the SAFT EOS. If made equal to 0, the calculations are made using only the EOS. It is slower, but can be necessary if the ancillary functions do not match the EOS with enough precision.&nbsp;</p>
    <p>The referenceState Integer constant  must be made equal to 1 for using ASHRAE reference state, to 2 for using IIR, to 3 for using NBP, to 4 for using user reference_T and reference_p. Any other number will produce an unreferenced calculation for enthalpy and entropy.</p>
    <p>If you are going to use the BaseProperties model, you must make the inputChoice constant String equal to 'ph', 'pT' or 'dT', as needed. When instantiating the model, you can still change the selection just for the object, specifying the value for the localInputChoice constant.</p>
    <p><b>How it compares with the ExternalMedia library</b></p>
    <p>ExternalPure is not a reference library as can be RefProp or CoolProp. It has a lot of multiparameter EOS, but some times are not the newest ones. The calculation of derivatives is in many times numerical and not symbolical. And only pure or pseudo pure substances are covered, although I think that this is also the case for ExternalMedia. Nevertheless the results are still of high quality and enough for the normal practical work.</p><p>It has also some advantages over ExternalMedia:</p><p>- The speed is higher if for any reason TTSE can't be used (if TTSE can be used ExternalMedia/CoolProp is faster)</p><p>- The calculation of transport properties is much more flexible.</p><p>- The addition of new substances by the user is much easier.</p><p><br></p>
    <p><br></p>
    
    </body></html>"));
end ExternalPure;

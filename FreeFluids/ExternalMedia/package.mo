within FreeFluids;

package ExternalMedia "ExternalMedia.mo by Carlos Trujillo
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
  //***** PACKAGE ExternalTPMedium*****
  //***********************************

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

    constant String resDir = Modelica.Utilities.Files.loadResource("modelica://FreeFluids/Resources") "resources directory";
    constant FluidConstants fluidK "just for minimum data out of the SubstanceData external object";
    constant String description = " " "Description of the EOS and correlations used for the substance";
    constant Integer refState = 1 "1=ASHRAE, 2=IIR, 3=NBP, >3=User. When user, the value at reference_T and reference_p as 0";
    constant String inputChoice = "ph" "Allows to choose the input variables to use for the construction of a BaseProperties object. alternatives: ph,pT,dT";
    constant Integer thermoModel = 1 "1 for cubic EOS; 2 for SAFT; 3 Saltzmann and Wagner. All transport properties will be computed by correlations";
    //constant String dataPath=Modelica.Utilities.System.getWorkDirectory();
    //constant String directory=Modelica.Utilities.Files.fullPath("modelica://Mediums") + mediumName;
    //Definition of ThermodynamicState, BaseProperties and SubstanceData
    //------------------------------------------------------------------

    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        //Real nMols "number of moles in a kg";
        AbsolutePressure p(displayUnit = "bar") "Pressure in Pa";
        Temperature T(displayUnit = "degC") "Kelvin temperature of medium";
        SpecificEnthalpy h "specific enthalpy";
        SpecificEntropy s;
        Density d(displayUnit = "kg/m3") "Liquid density in kgr/m3";
        MoleFraction gf "gas molar fraction";
        //additional variables, that should be records, but OpenModelica does not work well with nested records.
        Density ld(displayUnit = "kg/m3") "Liquid density in kgr/m3";
        Density gd(displayUnit = "kg/m3") "Gas density in kgr/m3";
        SpecificEnthalpy lh, gh "liquid and gas specific enthalpy";
        SpecificEntropy ls, gs;
        SpecificHeatCapacity lCv, lCp, gCv, gCp;
        Real lDvp, lDvT, gDvp, gDvT;
    end ThermodynamicState;

    /*redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), T(stateSelect =  StateSelect.default), d(stateSelect = StateSelect.default))
                    constant String localInputChoice = inputChoice;
                    Integer phase(min=0, max=2, start=1, fixed=false);
              
                  algorithm
                    sat := setSat_p(p);
                    state := setState_phX(p, h);
                    MM := fluidK.molarMass "in kg/mol";
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
        MM := fluidK.molarMass "in kg/mol";
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

    //Code for the creation and destruction of the external object with the substance data

    type SubstanceData "implemented as an ExternalObject"
      extends ExternalObject;

      function constructor
        input String name;
        input String resDir;
        input Integer thermoModel;
        input Real refT;
        input Real refP;
        output SubstanceData data;
      
        external "C" data = FF_createSubstanceData(name, resDir, thermoModel, refState, refT, refP) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
      end constructor;

      function destructor "Release storage"
        input SubstanceData data;
      
        external "C" FF_destroySubstanceData(data) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
      end destructor;
    end SubstanceData;

    constant SubstanceData data = SubstanceData(mediumName, resDir, thermoModel, reference_T, reference_p);
    //Basic functions
    //---------------

    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;


        external "C" FF_saturationPressure(data, T, p) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end saturationPressure;

    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;


        external "C" FF_saturationTemperature(data, p, T) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end saturationTemperature;

    redeclare function extends saturationTemperature_derp "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;

      algorithm
        dTp := if p < fluidConstants[1].criticalPressure then (saturationTemperature(p) - saturationTemperature(0.999 * p)) / (0.001 * p) else 0;
    end saturationTemperature_derp;

    function densitiesEOS_pT "Return liquid and gas densities at give temperature and pressure, by EOS"
      input AbsolutePressure p;
      input Temperature T;
      input Integer phase = 0 "0=both phases, 1= only gas, 2=only liquid";
      output Density ld;
      output Density gd;
    
      external "C" FF_densitiesEOS_pT(data, p, T, phase, ld, gd) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end densitiesEOS_pT;

    function solveEOS_Tp "Calculate densities, Arr and its T derivative, and ideal part, from T and p"
      input Temperature T;
      input AbsolutePressure p;
      input Integer ph;
      output Real gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEOS_Tp(data, T, p, ph, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_Tp;

    function solveEOS_Td "Calculate densities, Arr and its T derivative, and ideal part"
      input Temperature T;
      input Density d;
      output Real p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEOS_Td(data, T, d, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_Td;

    function solveEOS_ph "Calculates T, densities, Arr and its T derivative, and ideal part"
      input AbsolutePressure p;
      input SpecificEnthalpy h;
      output Real T, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEOS_ph(data, p, h, T, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_ph;

    function solveEOS_ps "Calculate T, densities, Arr and its T derivative, and ideal part"
      input AbsolutePressure p;
      input SpecificEntropy s;
      output Real T, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEOS_ps(data, p, s, T, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_ps;

    function thermoPropertiesEOS_Td "Return all thermo properties from T, densities and gas fraction"
      input Real T, ld, gd, gf;
      output Real h, u, s, cp, cv, dP_dT, dP_dV, ss, jt, it;
    
      external "C" FF_thermoPropertiesEOS_Td(data, T, ld, gd, gf, h, u, s, cp, cv, dP_dT, dP_dV, ss, jt, it) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end thermoPropertiesEOS_Td;

    function specificEnthalpyEOS_Td "Return specific enthalpy from T, densities and gas fraction"
      input Real T, ld, gd, gf;
      output Real h;
    
      external "C" FF_specificEnthalpy(data, T, ld, gd, gf, h) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end specificEnthalpyEOS_Td;

    function pressureEOS_dT "Return pressure given temperature and density, by EOS"
      input Density d;
      input Temperature T;
      output AbsolutePressure p;
    
      external "C" FF_pressureEOS_dT(data, d, T, p) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end pressureEOS_dT;

    //Establish general states
    //------------------------

    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi. The function from T and P is unable to compute any gas fraction different from 0 or 1"
        extends Modelica.Icons.Function;

      protected
        Integer ph = 5;

      algorithm
        state.p := p;
        state.T := T;
        state.phase := 1;
        (state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(T, p, ph);
        if state.ld < 1.0 then
          state.gf := 1.0;
          state.d := state.gd;
          state.h := state.gh;
          state.s := state.gs;
        else
          state.gf := 0.0;
          state.d := state.ld;
          state.h := state.lh;
          state.s := state.ls;
        end if;
    end setState_pTX;

    redeclare function extends setState_dTX "Return thermodynamic state as function of density, T and composition X or Xi"
        extends Modelica.Icons.Function;

      algorithm
        state.T := T;
        state.d := d;
        (state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Td(T, d);
        state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else state.gd * (state.ld - d) / (d * (state.ld - state.gd));
        state.h := state.lh * (1 - state.gf) + state.gh * state.gf;
        state.s := state.ls * (1 - state.gf) + state.gs * state.gf;
        state.phase := if state.gf < 0.000001 then 1 else if state.gf > 0.999999 then 1 else 2;
    end setState_dTX;

    //here is necessary to arrive to calculate the thermoprops, so better to do every thing inside C

    redeclare function extends setState_phX "Return thermodynamic state as function of pressure, enthalpy and composition X or Xi"
        extends Modelica.Icons.Function;

      algorithm
        state.p := p;
        state.h := h;
        (state.T, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_ph(p, h);
        state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else (state.h - state.lh) / (state.gh - state.lh);
        state.d := if state.gf > 0.999999 then state.gd else if state.gf < 0.000001 then state.ld else state.ld * state.gd / (state.gf * state.ld + (1.0 - state.gf) * state.gd);
        state.s := state.ls * (1 - state.gf) + state.gs * state.gf;
        state.phase := if state.gf < 0.000001 then 1 else if state.gf > 0.999999 then 1 else 2;
    end setState_phX;

    redeclare function extends setState_psX "Return thermodynamic state as function of pressure, entropy and composition X or Xi"
        extends Modelica.Icons.Function;

      algorithm
        state.p := p;
        state.s := s;
        (state.T, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_ps(p, s);
        state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else (state.s - state.ls) / (state.gs - state.ls);
        state.d := if state.gf > 0.999999 then state.gd else if state.gf < 0.000001 then state.ld else state.ld * state.gd / (state.gf * state.ld + (1.0 - state.gf) * state.gd);
        state.h := state.lh * (1 - state.gf) + state.gh * state.gf;
        state.phase := if state.gf < 0.000001 then 1 else if state.gf > 0.999999 then 1 else 2;
    end setState_psX;

    //Establish special states
    //------------------------

    redeclare function extends setDewState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the dew point"
        extends Modelica.Icons.Function;

      protected
        Integer ph = 1 "1 indicates gas phase";

      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.phase := 1;
        state.gf := 1.0 "gas";
        (state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(sat.Tsat, sat.psat, ph);
        state.d := state.gd;
        state.h := state.gh;
        state.s := state.gs;
    end setDewState;

    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;

      protected
        Integer ph = 2;

      algorithm
        state.p := sat.psat;
        state.T := sat.Tsat;
        state.phase := 1;
        state.gf := 0.0 "liquid";
        (state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(sat.Tsat, sat.psat, ph);
        state.d := state.ld;
        state.h := state.lh;
        state.s := state.ls;
    end setBubbleState;

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
        if state.gf < 0.000001 then
          ddpT := state.ld * isothermalCompressibility(state);
        elseif state.gf > 0.999999 then
          ddpT := state.gd * isothermalCompressibility(state);
        else
          ddpT := 0;
        end if;
    end density_derp_T;

    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf < 0.000001 then
          ddTp := -state.ld * isobaricExpansionCoefficient(state);
        elseif state.gf > 0.999999 then
          ddTp := -state.gd * isobaricExpansionCoefficient(state);
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
        if state.gf < 0.000001 then
          ddph := -state.ld * state.ld * (state.lDvp + (state.T * state.lDvT * state.lDvT - state.lDvT / state.ld) / state.lCp);
        elseif state.gf > 0.999999 then
          ddph := -state.gd * state.gd * (state.gDvp + (state.T * state.gDvT * state.gDvT - state.gDvT / state.gd) / state.gCp);
        else
          DTp := (1 / state.gd - 1 / state.ld) / (state.gs - state.ls);
          lDvp := state.lDvp + state.lDvT * DTp;
          gDvp := state.gDvp + state.gDvT * DTp;
          lDhp := 1 / state.ld - state.T * state.lDvT + state.lCp * DTp;
          gDhp := 1 / state.gd - state.T * state.gDvT + state.gCp * DTp;
          Dgfp_h := (state.gf * gDhp + (1 - state.gf) * lDhp) / (state.lh - state.gh);
          ddph := -state.d * state.d * (lDvp + Dgfp_h * (1 / state.gd - 1 / state.ld) + state.gf * (gDvp - lDvp));
        end if;
    end density_derp_h;

    redeclare function extends density_derh_p "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;

      algorithm
        if state.gf < 0.000001 then
          ddhp := -state.ld * state.ld * state.lDvT / state.lCp;
        elseif state.gf > 0.999999 then
          ddhp := -state.gd * state.gd * state.gDvT / state.gCp;
        else
          ddhp := -state.d * state.d * (1 / state.gd - 1 / state.ld) / (state.gh - state.lh);
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

      algorithm
        s := state.s;
    end specificEntropy;

    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
        extends Modelica.Icons.Function;

      algorithm
        g := state.h - state.T * state.s;
    end specificGibbsEnergy;

    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;

      algorithm
        f := state.h - state.p / state.d - state.T * state.s;
    end specificHelmholtzEnergy;

    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;

      algorithm
        cp := if state.gf < 0.000001 then state.lCp else if state.gf > 0.999999 then state.gCp else 0.0;
    end specificHeatCapacityCp;

    redeclare function extends specificHeatCapacityCv "Return specific heat capacity at constant volume"
        extends Modelica.Icons.Function;

      protected
        Real lDsT, gDsT, DpT, lDvT, gDvT "total derivatives of s, p and v along the saturation line, regarding T";
        Real DgfT_v "partial derivative of gas fraction regarding T at constant v";

      algorithm
        if state.gf < 0.000001 then
          cv := state.lCv;
        elseif state.gf > 0.999999 then
          cv := state.gCv;
        else
          DpT := (state.gs - state.ls) / (1 / state.gd - 1 / state.ld);
          lDsT := state.lCp / state.T - state.lDvT * DpT;
          gDsT := state.gCp / state.T - state.gDvT * DpT;
          lDvT := state.lDvT + state.lDvp * DpT;
          gDvT := state.gDvT + state.gDvp * DpT;
          DgfT_v := (state.gf * gDvT + (1 - state.gf) * lDvT) / (1 / state.ld - 1 / state.gd);
          cv := state.T * lDsT + state.T * DgfT_v * (state.gs - state.ls) + state.gf * state.T * (gDsT - lDsT);
        end if;
    end specificHeatCapacityCv;

    redeclare function extends isentropicExponent "Return isentropic exponent"
        extends Modelica.Icons.Function;

      protected
        Real Vinv, cv, cp;

      algorithm
        gamma := if state.gf < 0.000001 then state.lCp / state.lCv else if state.gf > 0.999999 then state.gCp / state.gCv else 0.0;
    end isentropicExponent;

    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;

      algorithm
        a := if state.gf == 0.0 then (-state.lCp / (state.lDvp * state.lCv)) ^ 0.5 / state.ld else if state.gf == 1.0 then (-state.gCp / (state.gDvp * state.gCv)) ^ 0.5 / state.gd else 0.0;
    end velocityOfSound;

    redeclare function extends isobaricExpansionCoefficient "Returns the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;

      protected
        Real Vinv;

      algorithm
        beta := if state.gf < 0.000001 then state.lDvT * state.ld else if state.gf > 0.999999 then state.gDvT * state.gd else 0.0;
    end isobaricExpansionCoefficient;

    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;

      protected
        Real Vinv;

      algorithm
        kappa := if state.gf < 0.000001 then -state.lDvp * state.ld else if state.gf > 0.999999 then -state.gDvp * state.gd else 0.0;
    end isothermalCompressibility;

    redeclare function extends dynamicViscosity "Return dynamic viscosity from a ThermodynamicState record"
        extends Modelica.Icons.Function;


        external "C" FF_dynamicViscosity(data, state.T, state.p, state.gf, eta) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end dynamicViscosity;

    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;


        external "C" FF_thermalConductivity(data, state.T, state.p, state.gf, lambda) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end thermalConductivity;

    function liquidDynamicViscosity "Return liquid dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    
      external "C" FF_dynamicViscosity(data, state.T, state.p, 0.0, eta) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end liquidDynamicViscosity;

    function gasDynamicViscosity "Return gas dynamic viscosity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output DynamicViscosity eta "Dynamic viscosity";
    
      external "C" FF_dynamicViscosity(data, state.T, state.p, 1.0, eta) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end gasDynamicViscosity;

    function liquidThermalConductivity "Return liquid thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    
      external "C" FF_thermalConductivity(data, state.T, state.p, 0.0, lambda) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end liquidThermalConductivity;

    function gasThermalConductivity "Return gas thermal conductivity from a ThermodynamicState record"
      input ThermodynamicState state "Thermodynamic state record";
      output ThermalConductivity lambda "Thermal conductivity";
    
      external "C" FF_thermalConductivity(data, state.T, state.p, 1.0, lambda) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end gasThermalConductivity;

    redeclare function extends surfaceTension "Return surface tension"
        extends Modelica.Icons.Function;


        external "C" FF_surfaceTension(data, sat.Tsat, sigma) annotation(
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
        (dl, dg) := densitiesEOS_pT(sat.psat, sat.Tsat, 2);
    end bubbleDensity;

    redeclare function extends dBubbleDensity_dPressure "Return bubble point density derivative"
        extends Modelica.Icons.Function;

      algorithm
        ddldp := if sat.Tsat < fluidConstants[1].criticalTemperature then (bubbleDensity(setSat_T(sat.Tsat - 0.01)) - bubbleDensity(sat)) / (saturationPressure(sat.Tsat - 0.01) - sat.psat) else 0;
    end dBubbleDensity_dPressure;

    redeclare function extends dewDensity "Return dew point density"
        extends Modelica.Icons.Function;

      protected
        Real dl;

      algorithm
        (dl, dv) := densitiesEOS_pT(sat.psat, sat.Tsat, 1);
    end dewDensity;

    redeclare function extends dDewDensity_dPressure "Return dew point density derivative"
        extends Modelica.Icons.Function;

      algorithm
        ddvdp := if sat.Tsat < fluidConstants[1].criticalTemperature then (dewDensity(setSat_T(sat.Tsat - 0.02)) - dewDensity(sat)) / (saturationPressure(sat.Tsat - 0.02) - sat.psat) else 0;
    end dDewDensity_dPressure;

    redeclare function extends bubbleEnthalpy "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;


        external "C" FF_bubbleEnthalpy(data, sat.psat, sat.Tsat, hl) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end bubbleEnthalpy;

    redeclare function extends dBubbleEnthalpy_dPressure "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;

      algorithm
        dhldp := if sat.Tsat < fluidConstants[1].criticalTemperature then (bubbleEnthalpy(setSat_T(sat.Tsat - 0.01)) - bubbleEnthalpy(sat)) / (saturationPressure(sat.Tsat - 0.01) - sat.psat) else 0;
    end dBubbleEnthalpy_dPressure;

    redeclare function extends dewEnthalpy "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;


        external "C" FF_dewEnthalpy(data, sat.psat, sat.Tsat, hv) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end dewEnthalpy;

    redeclare function extends dDewEnthalpy_dPressure "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;

      algorithm
        dhvdp := if sat.Tsat < fluidConstants[1].criticalTemperature then (dewEnthalpy(setSat_T(sat.Tsat - 0.02)) - dewEnthalpy(sat)) / (saturationPressure(sat.Tsat - 0.02) - sat.psat) else 0;
    end dDewEnthalpy_dPressure;

    function vaporizationEnthalpy
      input SaturationProperties sat;
      output SpecificEnthalpy Hv;
    algorithm
      Hv := dewEnthalpy(sat) - bubbleEnthalpy(sat);
    end vaporizationEnthalpy;

    redeclare function extends bubbleEntropy "Return bubble point specific entropy"
        extends Modelica.Icons.Function;


        external "C" FF_bubbleEntropy(data, sat.psat, sat.Tsat, sl) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end bubbleEntropy;

    redeclare function extends dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;


        external "C" FF_dewEntropy(data, sat.psat, sat.Tsat, sv) annotation(
          IncludeDirectory = "modelica://FreeFluids/Resources",
          Include = "#include \"FFmodelicaMedium.c\"");
    end dewEntropy;
  end ExternalMedium;

  package Fluids


    package Acetone
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Acetone", fluidK(casRegistryNumber = "67-64-1", description = "Multiparameter:Lemmon&Span 2006.", molarMass = 0.058079, criticalTemperature = 5.081000e+02, criticalPressure = 4.700000e+06), onePhase = false, thermoModel = 3, refState = 2);
    end Acetone;

    package Ammonia
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Ammonia",
 fluidK(casRegistryNumber = "7664-41-7", description = "Multiparameter: Tillner-Roth. PCSAFT: Mejbri 2005, Cubic: PRMC C.Trujillo 2019. Cp0: Willhoit", molarMass = 1.703026e-02, criticalTemperature = 4.054000e+02, criticalPressure = 1.133300e+07),
 final onePhase=false, thermoModel=3, refState=2);
    end Ammonia;
  
    package Butane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Butane_n", fluidK(casRegistryNumber = "106-97-8", description = "Multiparameter:Buecker&Wagner 2006. PCSAFT:Gross&Sadowski 2001. Cubic:PRMC, C.Trujillo 2019", molarMass = 0.05812, criticalTemperature = 4.251250e+02, criticalPressure = 3.796000e+06), onePhase = false, thermoModel = 3, refState = 2);
    end Butane_n;
  
    package Butanol_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Butanol_n",
   fluidK(casRegistryNumber = "71-36-3", description = "Multiparameter: none. PCSAFT: 2B C.Trujillo 2019. Cubic: PRMC. Cp0: DIPPR107", molarMass = 7.414000e-02, criticalTemperature = 5.630000e+02, criticalPressure = 4.414000e+06),
   final onePhase=false, thermoModel=2, refState=2);
    end Butanol_n;
  
    package CO2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "CO2", fluidK(casRegistryNumber = "124-38-9", description = "Multiparameter: GERG2004. PCSAFT and PRMC: C.Trujillo from GERG2004. Cp0:Jaeschke", molarMass = 4.400980e-02, criticalTemperature = 3.041282e+02, criticalPressure = 7.377730e+06), final onePhase = false, thermoModel = 1, refState = 2);
    end CO2;
  
  package Dichlorodifluormethane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Dichlorodifluormethane",
   fluidK(casRegistryNumber = "75-71-8", description = "Multiparameter:Marx et alt. PCSAFT:C.Trujillo 2020, Cubic:PRMC. Cp0:from CoolProp", molarMass = 1.209130e-01, criticalTemperature = 3.851200e+02, criticalPressure = 4.136100e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Dichlorodifluormethane;
  
    package EG
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "EG",
   fluidK(casRegistryNumber = "107-21-1", description = "Multiparameter: none. PPCSAFT_GV: 2B C.Trujillo 2019. Cubic: PRMC", molarMass = 6.207000e-02, criticalTemperature = 7.200000e+02, criticalPressure = 8.200000e+06),
   final onePhase=false, thermoModel=2, refState=2);
    end EG;
  
    package Ethane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Ethane",
   fluidK(casRegistryNumber = "74-84-0", description = "Multiparameter: Buecker-Wagner 2006. PCSAAFT: Gross-Sadowski 2001. Cubic:SRKMC Chemsep", molarMass = 3.006904e-02, criticalTemperature = 3.053220e+02, criticalPressure = 4.872200e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Ethane;
  
    package Heptane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Heptane_n", fluidK(casRegistryNumber = "142-82-5", description = "Multiparameter: Span and Wagner 2003. PCSAFT: Gross 2001. Cubic: PRMC", molarMass = 1.002020e-01, criticalTemperature = 5.401300e+02, criticalPressure = 2.736000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Heptane_n;
  
    package Hexane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Hexane_n", fluidK(casRegistryNumber = "110-54-3", description = "Multiparameter: Spand and Wagner 2003. PCSAFT: Gross and Sadowski. Cubic: PRMC", molarMass = 8.617536e-02, criticalTemperature = 5.078200e+02, criticalPressure = 3.034000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Hexane_n;
  
    package Isobutane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Isobutane",
   fluidK(casRegistryNumber = "75-28-5", description = "Multiparameter: Buecker-Wagner 2006. PCSAFT: Gross-Sadowski 2001. Cubic: PRMC. Cp0: Cooper", molarMass = 5.812220e-02, criticalTemperature = 4.078170e+02, criticalPressure = 3.629000e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Isobutane;
  
    package Methane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Methane",
   fluidK(casRegistryNumber = "74-82-8", description = "Multiparameter: Seltzmann-Wagner 1991. PCSAF: Gross-Sadowski 2001. Cubic: SRKMC from Chemsep. Cp0: Cooper", molarMass = 1.604280e-02, criticalTemperature = 1.905640e+02, criticalPressure = 4.599200e+06),
   final onePhase=false, thermoModel=1, refState=2);
    end Methane;
  
    package N2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "N2",
   fluidK(casRegistryNumber = "7727-37-9", description = "Multiparameter: Span 2001, PCSAFT: Gross 2001. PRMC: Barragan 2002. Cp0: DIPPR107", molarMass = 2.801348e-02, criticalTemperature = 1.261920e+02, criticalPressure = 3.395800e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end N2;
  
    package O2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "O2",
   fluidK(casRegistryNumber = "7782-44-7", description = "Multiparameter: Schmidt 1985. PCSAFT: Economou 2007. Cubic:SRKMC Chemsep. Cp0: Cooper", molarMass = 3.199880e-02, criticalTemperature = 1.545810e+02, criticalPressure = 5.043000e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end O2;
  
    package Pentane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Pentane_n",
   fluidK(casRegistryNumber = "109-66-0", description = "Multiparameter:Span-Wagner 2003. PCSAFT: Gross-Sadowski 2001. Cubic PRMC. Cp0: Wilhoit", molarMass = 7.215000e-02, criticalTemperature = 4.697000e+02, criticalPressure = 3.370000e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Pentane_n;
  
    package Propane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Propane",
   fluidK(casRegistryNumber = "74-98-6", description = "Multiparameter: Lemmon 2009. PCSAFT: Gross-Sadowski 2001. Cubic: SRKMC Chemsep. Cp0:Cooper", molarMass = 4.409562e-02, criticalTemperature = 3.698900e+02, criticalPressure = 4.251200e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Propane;
  
    package R134A
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "R134A", fluidK(casRegistryNumber = "811-97-2", description = "Multiparameter:Tillner 1994. PCSAFT: non assoc.C.Trujillo 2019. Cubic:PRMC, C.Trujillo 2019", molarMass = 1.020310e-01, criticalTemperature = 3.741800e+02, criticalPressure = 4.901200e+06), final onePhase = false, thermoModel = 3, refState = 1);
    end R134A;
  
    package R410A
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "R410A",
   fluidK(casRegistryNumber = "R410A.PPF", description = "Multiparameter: Lemmon 2003. PCSAFT: non assoc. C.T.2019. Cubic: none. Cp0: Cooper", molarMass = 7.258540e-02, criticalTemperature = 3.444940e+02, criticalPressure = 4.901200e+06),
   final onePhase=false, thermoModel=1, refState=2);
    end R410A;
  
    package WaterRef
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "WaterRef", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:IAPWS95. PCSAFT:C.Trujillo 2B max.633K. Cubic:PRMC,C.Trujillo 2018. Cp0:Jaeschke", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), final onePhase = false, thermoModel = 3, refState = 4);
    end WaterRef;
  
    package Water
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Water", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:GERG2004. PCSAFT: Diamantotis 4C. Cubic: SRKMC C.Trujillo. Cp0:Jaechske", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), onePhase = false, thermoModel = 1, refState = 4, reference_T = 273.15, reference_p = 101325);
    end Water;

  end Fluids;

  package Tests
    model Test1Cubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.WaterRef(thermoModel = 1, refState = 4, reference_T = 273.15, reference_p = 101325.0, inputChoice = "ph"), p = 20.0e5, initialT = 0.1 + 273.15, finalT = 470 + 273.15);
    end Test1Cubic;

    model Test1PCSAFT
      extends Test1Cubic(Medium(thermoModel = 2));
    end Test1PCSAFT;

    model Test1SW
      extends Test1Cubic(Medium(thermoModel = 3));
    end Test1SW;

    model Test1SW2
      extends Test1Cubic(redeclare replaceable package Medium = ExternalMedia.Fluids.Water(thermoModel = 3, refState = 4, reference_T = 273.15, reference_p = 101325.0, inputChoice = "ph"));
    end Test1SW2;

    model Test1TMedia
      extends Test1Cubic(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", highPressure = true));
      //extends Test1Cubic(redeclare package Medium = FreeFluids.TMedia.Fluids.R134A(refState = "ASHRAE"));
    end Test1TMedia;

    model Test1IF97
      extends Test1Cubic(redeclare package Medium = Modelica.Media.Water.StandardWater);
      //extends Test1A(redeclare package Medium = FreeFluids.TMedia.Acetone(refState="IIR"));
    end Test1IF97;

    model Test1AdCubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.R410A(thermoModel = 1, refState = 2, reference_T = 273.15, reference_p = 101325.0, inputChoice = "ph"), p = 35.0e5, initialT = (-50) + 273.15, finalT = 55 + 273.15);
    end Test1AdCubic;
    
    model Test1AdPCSAFT
      extends Test1AdCubic(Medium(thermoModel = 2));
    end Test1AdPCSAFT;
  
    model Test1AdSW
      extends Test1AdCubic(Medium(thermoModel = 3));
    end Test1AdSW;

    model TestRefrigerantEvaporator "The test case has been copied from ThermoPower Library. Test case with once-through evaporator using Dittus-Boelter 2-phase heat transfer model"
      replaceable package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1) constrainedby Modelica.Media.Interfaces.PartialMedium;
      parameter Integer Nnodes = 10 "number of nodes";
      parameter Modelica.SIunits.Length Dint = 5e-3 "Internal diameter of refrigerant tube";
      parameter Modelica.SIunits.Length th = 0.5e-3 "Thickness of refrigerant tube";
      parameter Modelica.SIunits.Length Dext = Dint + 2 * th "External diameter of refrigerant tube";
      parameter Modelica.SIunits.Length Dgas = 20e-3 "External diameter of gas tube";
      parameter Modelica.SIunits.Length L = 10 "Length of refrigerant tube";
      constant Real pi = Modelica.Constants.pi;
      Modelica.SIunits.Temperature Tgas[:] = gasFlow.T[end:(-1):1] "Alias variable for direct comparisons";
      Modelica.SIunits.Temperature Tfluid[:] = fluidFlow.T "Alias variable for direct comparisons";
      Modelica.SIunits.Temperature Twall[:] = tubeWall.Tvol "Alias variable for direct comparisons";
      Modelica.SIunits.SpecificEnthalpy hl[:] = fluidFlow.hl * ones(Nnodes);
      Modelica.SIunits.SpecificEnthalpy hv[:] = fluidFlow.hv * ones(Nnodes);
      ThermoPower.Water.SinkPressure fluidSink(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1), p0 = 3e6) annotation(
        Placement(transformation(extent = {{50, -82}, {70, -62}}, rotation = 0)));
      ThermoPower.Gas.SinkPressure gasSink(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{-54, 18}, {-74, 38}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow fluidSource(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1), h = -100e3, p0 = 3e5, use_in_h = false, use_in_w0 = true) annotation(
        Placement(transformation(extent = {{-70, -82}, {-50, -62}}, rotation = 0)));
      ThermoPower.Water.SensT fluid_T_in(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1)) annotation(
        Placement(transformation(extent = {{-46, -78}, {-26, -58}}, rotation = 0)));
      ThermoPower.Gas.SensT gas_T_in(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{34, 22}, {14, 42}}, rotation = 0)));
      ThermoPower.Gas.SourceMassFlow gasSource(redeclare package Medium = ThermoPower.Media.Air, T = 380.15, use_in_w0 = true) annotation(
        Placement(transformation(extent = {{64, 18}, {44, 38}}, rotation = 0)));
      ThermoPower.Water.SensT fluid_T_out(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1)) annotation(
        Placement(transformation(extent = {{20, -78}, {40, -58}}, rotation = 0)));
      ThermoPower.Gas.SensT gas_T_out(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{-24, 22}, {-44, 42}}, rotation = 0)));
      inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.steadyState) annotation(
        Placement(transformation(extent = {{80, 80}, {100, 100}})));
      ThermoPower.Gas.Flow1DFV gasFlow(redeclare package Medium = ThermoPower.Media.Air, redeclare model HeatTransfer = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = 120), A = (Dgas ^ 2 - Dext ^ 2) / 4 * pi, Dhyd = Dgas, FFtype = ThermoPower.Choices.Flow1D.FFtypes.NoFriction, HydraulicCapacitance = ThermoPower.Choices.Flow1D.HCtypes.Downstream, L = L, N = Nnodes, Nt = 1, Tstartin = 273.15, Tstartout = 283.15, noInitialPressure = true, omega = Dext * pi, wnom = 0.02) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-8, 28})));
      ThermoPower.Water.Flow1DFV2ph fluidFlow(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2(thermoModel = 3, refState = 1), redeclare model HeatTransfer = ThermoPower.Thermal.HeatTransferFV.FlowDependentHeatTransferCoefficient2ph(gamma_nom_liq = 800, gamma_nom_2ph = 8000, gamma_nom_vap = 800), A = Dint ^ 2 / 4 * pi, Cfnom = 0.005, Dhyd = Dint, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, HydraulicCapacitance = ThermoPower.Choices.Flow1D.HCtypes.Downstream, L = L, N = Nnodes, dpnom = 1000, hstartin = 250e3, hstartout = 250e3, noInitialPressure = true, omega = Dint * pi, pstart = 3e6, wnom = 0.005) annotation(
        Placement(visible = true, transformation(extent = {{-18, -82}, {2, -62}}, rotation = 0)));
      ThermoPower.Thermal.CounterCurrentFV counterCurrentFV(Nw = Nnodes - 1) annotation(
        Placement(transformation(extent = {{-18, -10}, {2, 10}})));
      ThermoPower.Thermal.MetalTubeFV tubeWall(Nw = Nnodes - 1, initOpt = ThermoPower.Choices.Init.Options.steadyState, L = L, rint = Dint / 2, rhomcm = 7800 * 400, lambda = 300, WallRes = false, rext = Dext / 2) annotation(
        Placement(transformation(extent = {{-18, -20}, {2, -40}})));
      Modelica.Blocks.Sources.Ramp ramp(height = 0, duration = 1, offset = 0.01) annotation(
        Placement(transformation(extent = {{-20, 62}, {0, 82}})));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1e6, height = 0.001, offset = 0.0025) annotation(
        Placement(visible = true, transformation(origin = {-106, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ramp1.y, fluidSource.in_w0) annotation(
        Line(points = {{-95, -22}, {-64, -22}, {-64, -66}}, color = {0, 0, 127}));
      connect(fluidFlow.outfl, fluid_T_out.inlet) annotation(
        Line(points = {{2, -72}, {24, -72}}, color = {0, 0, 255}));
      connect(tubeWall.int, fluidFlow.wall) annotation(
        Line(points = {{-8, -33}, {-8, -67}}, color = {255, 127, 0}));
      connect(fluid_T_in.outlet, fluidFlow.infl) annotation(
        Line(points = {{-30, -72}, {-18, -72}}, color = {0, 0, 255}));
      connect(fluidSource.flange, fluid_T_in.inlet) annotation(
        Line(points = {{-50, -72}, {-42, -72}}, thickness = 0.5, color = {0, 0, 255}));
      connect(fluid_T_out.outlet, fluidSink.flange) annotation(
        Line(points = {{36, -72}, {50, -72}}, thickness = 0.5, color = {0, 0, 255}));
      connect(counterCurrentFV.side1, gasFlow.wall) annotation(
        Line(points = {{-8, 3}, {-8, 23}}, color = {255, 127, 0}, smooth = Smooth.None));
      connect(counterCurrentFV.side2, tubeWall.ext) annotation(
        Line(points = {{-8, -3.1}, {-8, -26.9}}, color = {255, 127, 0}, smooth = Smooth.None));
      connect(gas_T_out.inlet, gasFlow.outfl) annotation(
        Line(points = {{-28, 28}, {-18, 28}}, color = {159, 159, 223}, smooth = Smooth.None));
      connect(gasFlow.infl, gas_T_in.outlet) annotation(
        Line(points = {{2, 28}, {18, 28}}, color = {159, 159, 223}, smooth = Smooth.None));
      connect(gas_T_in.inlet, gasSource.flange) annotation(
        Line(points = {{30, 28}, {44, 28}}, color = {159, 159, 223}, smooth = Smooth.None));
      connect(gasSink.flange, gas_T_out.outlet) annotation(
        Line(points = {{-54, 28}, {-40, 28}}, color = {159, 159, 223}, smooth = Smooth.None));
      connect(ramp.y, gasSource.in_w0) annotation(
        Line(points = {{1, 72}, {60, 72}, {60, 33}}, color = {0, 0, 127}, smooth = Smooth.None));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics),
        __Dymola_experimentSetupOutput,
        Documentation(info = "<html>
  <p>The model is designed to test the component <code>ThermoPower.Water.Flow1DFV</code> (fluid side of a heat exchanger, model uses finite volumes).</p><p>This model represent the two fluid sides of a heat exchanger made by two concentric tubes in counterflow configuration. The thickness of the wall separating the two tubes is negligible. The operating fluid is liquid ThermoPower.Water. The mass flow rate during the experiment and initial conditions are the same for the two sides. </p><p>During the simulation, the inlet specific enthalpy for hexA (&QUOT;hot side&QUOT;) is changed at time t = 50 s. The outlet temperature of the hot side starts changing after the fluid transport time delay, while the outlet temperature of the cold side starts changing immediately. </p>
  <p>Simulation Interval = [0...1200] sec </p><p>Integration Algorithm = DASSL </p><p>Algorithm Tolerance = 1e-6 </p>
  </html>", revisions = "<html>
  <ul>
  <li>18 Sep 2013 by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br/>Updated to new FV structure. Updated parameters.</li></li>
  <li><i>1 Oct 2003</i> by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
  First release.</li>
  </ul>
  </html>"),
        experiment(StopTime = 1e+006, Interval = 2000),
        uses(Modelica(version = "3.2.2")));
    end TestRefrigerantEvaporator;
  end Tests;
  annotation(
    Documentation(info = "<html>
    <body>
    <p>The medium is designed for pure or pseudo-pure substances in liquid, gas, or two phases. The thermodynamic calculations are performed using different equations of state(EOS), implemented in C language, that are called using external object and functions.</p>
    <p>For each substance there is the possibility of choising between three EOS types: multiparameter, PCSAFT, or cubic. The external object contains the substance data, and is a C structure, that is exported from the database using the FreeFluids GUI.</p>
    <p>The quality of the results is very high when multiparameter (Seltzmann and Wagner, SW) EOS are available. If not, PCSAFT or different flavours of cubic EOS can be used. The transport properties are computed from temperature dependent correlations, with pressure correction.</p>
    <p>The medium implements all the requirements of Modelica.Media.PartialTwoPhaseMedium. It is compatible with the old frontend of OpenModelica 14.1, but not with the new frontend, as it doesn`t support external functions yet. It is also compatible, at least partially, with the ThermoPower library. Compared with other high quality libraries, it is very easy to use and, although it is not as complete as for example CooProp, or RefProp, the quality of the thermodynamic results is practically the same. In plus there is the possibility of PCSAFT and cubic EOS for substances for which no multiparameter EOS is available.</p>
    <p>The main limitation of the medium, when using SW or PCSAFT EOS, is in the vecinity of the critical point, as no special technique, as for example splines, is used. Nevertheless you can go normally quite close to the critical point. For transport properties, a pending point is to use multiphase, temperature and density dependent functions, when available. The medium is much slower than the TMedia one, but is the price for the better precision and wider application.</p>
    <p>The ThermodynamicState record is quite big, in order to allow all the thermodynamic calculations with just one call to external functions.</p>
    <p>The BaseProperties model has been implemented using algorithms instead of equations. I do not know if this is correct, but seems the logical decision when you know the order in which calculations must be performed.</p>
    <p>When extending the medium, the following configuration must be done:</p>
    <p>The thermoModel Integer constant must be made equal to 1 for using the cubic EOS, 2 for the PCSAFT, or 3 for the SW.</p>
    <p>The referenceState Integer constant  must be made equal to 1 for using ASHRAE reference state, to 2 for using IIR, to 3 for using NBP, to 4 for using user reference_T and reference_p. Any other number will produce an unreferenced calculation for enthalpy and entropy.</p>
    <p>If you are going to use the BaseProperties model, you must make the inputChoice constant String equal to 'ph', 'pT' or 'dT', as needed. When instantiating the model, you can still change the selection just for the object, specifying the value for the localInputChoice constant.</p>
    <p>In the package Tests you can compare the performance of the cubic, PCSAFT, IAPWS95,  GERG2004 EOS, and TMedia for water, against the standard Modelica.Media.Water.StandardWater. There is also the ThermoPower TestRefrigerantEvaporator test, modified to use CO2 as liquid (you need to load the ThermoPower library).</p>
    </body>
    </html>"));
end ExternalMedia;

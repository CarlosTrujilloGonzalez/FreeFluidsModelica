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

    constant String refName = "Propane" "name of the reference fluid for ECS calculations. Defaults to propane";
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
        SpecificEnthalpy h "Overall specific enthalpy";
        SpecificEntropy s "Overall specific entropy";
        Density d(displayUnit = "kg/m3") "Overall density in kgr/m3";
        MoleFraction gf "gas molar fraction";
        //additional variables, that should be records, but OpenModelica does not work well with nested records.
        Density ld(displayUnit = "kg/m3") "Liquid phase density in kgr/m3";
        Density gd(displayUnit = "kg/m3") "Gas phase density in kgr/m3";
        SpecificEnthalpy lh, gh "liquid and gas phases specific enthalpy";
        SpecificEntropy ls, gs "liquid and gas phases specific entropy";
        SpecificHeatCapacity lCv, lCp, gCv, gCp "liquid and gas phases specific heats";
        Real lDvp, lDvT, gDvp, gDvT "liquid and gas phases derivatives of volume regarding pressure and temperature";
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

    constant SubstanceData data = SubstanceData(mediumName, resDir, thermoModel, reference_T, reference_p) "creates de substance object";
    constant SubstanceData refData = SubstanceData(refName, resDir, thermoModel, reference_T, reference_p) "creates de reference object for ECS calculation";
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

    function pressureEOS_dT "Return pressure given temperature and density, by EOS"
      input Density d;
      input Temperature T;
      output AbsolutePressure p;
    
      external "C" FF_pressureEOS_dT(data, d, T, p) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end pressureEOS_dT;

    function densitiesEOS_pT "Return liquid and gas densities at give temperature and pressure, by EOS, bellow the critical temperature"
      input AbsolutePressure p;
      input Temperature T;
      input Integer phase = 0 "0=both phases, 1= only gas, 2=only liquid";
      output Density ld;
      output Density gd;
    
      external "C" FF_densitiesEOS_pT(data, p, T, phase, ld, gd) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end densitiesEOS_pT;

    function solveEOS_Tp "Calls the external function that calculates some basic thermodynamic properties of the phases, used later for the definition of the ThermodynamicState and the calculation of  all thermodynamic properties"
      input Temperature x;
      input AbsolutePressure y;
      output Real T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEos("p", data, x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_Tp;

    function solveEOS_Td "Same as solveEOS_Tp, but from temperature and density"
      input Temperature x;
      input Density y;
      output Real T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEos("d", data, x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_Td;

    function solveEOS_ph "Same as solveEOS_Tp, but from pressure and specific enthalpy"
      input AbsolutePressure x;
      input SpecificEnthalpy y;
      output Real T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEos("h", data, x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_ph;

    function solveEOS_ps "Same as solveEOS_Tp, but from pressure and specific entropy"
      input AbsolutePressure x;
      input SpecificEntropy y;
      output Real T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT;
    
      external "C" FF_solveEos("s", data, x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT) annotation(
        IncludeDirectory = "modelica://FreeFluids/Resources",
        Include = "#include \"FFmodelicaMedium.c\"");
    end solveEOS_ps;

    //Establish general states
    //------------------------

    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi. The function from T and P is unable to compute any gas fraction different from 0 or 1"
        extends Modelica.Icons.Function;

      algorithm
//state.p := p;
//state.T := T;
        state.phase := 1;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(T, p);
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
//state.T := T;
        state.d := d;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Td(T, d);
        state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else state.gd * (state.ld - d) / (d * (state.ld - state.gd));
        state.h := state.lh * (1 - state.gf) + state.gh * state.gf;
        state.s := state.ls * (1 - state.gf) + state.gs * state.gf;
        state.phase := if state.gf < 0.000001 then 1 else if state.gf > 0.999999 then 1 else 2;
    end setState_dTX;

    //here is necessary to arrive to calculate the thermoprops, so better to do every thing inside C

    redeclare function extends setState_phX "Return thermodynamic state as function of pressure, enthalpy and composition X or Xi"
        extends Modelica.Icons.Function;

      algorithm
//state.p := p;
        state.h := h;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_ph(p, h);
        state.gf := if state.ld == 0 then 1.0 else if state.gd == 0 then 0.0 else (state.h - state.lh) / (state.gh - state.lh);
        state.d := if state.gf > 0.999999 then state.gd else if state.gf < 0.000001 then state.ld else state.ld * state.gd / (state.gf * state.ld + (1.0 - state.gf) * state.gd);
        state.s := state.ls * (1 - state.gf) + state.gs * state.gf;
        state.phase := if state.gf < 0.000001 then 1 else if state.gf > 0.999999 then 1 else 2;
    end setState_phX;

    redeclare function extends setState_psX "Return thermodynamic state as function of pressure, entropy and composition X or Xi"
        extends Modelica.Icons.Function;

      algorithm
//state.p := p;
        state.s := s;
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_ps(p, s);
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
//state.p := sat.psat;
//state.T := sat.Tsat;
        state.phase := 1;
        state.gf := 1.0 "gas";
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(sat.Tsat, sat.psat);
        state.d := state.gd;
        state.h := state.gh;
        state.s := state.gs;
    end setDewState;

    redeclare function extends setBubbleState "The input is a SaturationProperties record called sat. Returns the ThermodynamicState record at the bubble point"
        extends Modelica.Icons.Function;

      protected
        Integer ph = 2;

      algorithm
//state.p := sat.psat;
//state.T := sat.Tsat;
        state.phase := 1;
        state.gf := 0.0 "liquid";
        (state.T, state.p, state.gd, state.gh, state.gs, state.gCv, state.gCp, state.gDvp, state.gDvT, state.ld, state.lh, state.ls, state.lCv, state.lCp, state.lDvp, state.lDvT) := solveEOS_Tp(sat.Tsat, sat.psat);
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
        /*external "C" FF_dynamicViscosity(data, state.T, state.p, state.gf, eta) annotation(
                  IncludeDirectory = "modelica://FreeFluids/Resources",
                  Include = "#include \"FFmodelicaMedium.c\"");*/


        external "C" FF_Viscosity(data, refData, state.T, state.d, state.p, state.gf, eta) annotation(
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
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Acetone", fluidK(casRegistryNumber = "67-64-1", description = "Multiparameter:Lemmon&Span 2006. PCSAFT: Solms 2004. PRMC: Kleiman 2002. Cp0: Wilhoit", molarMass = 0.058079, criticalTemperature = 5.081000e+02, criticalPressure = 4.700000e+06), refName = "Propane", onePhase = false, thermoModel = 3, refState = 2);
    end Acetone;

    package Ammonia
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Ammonia", fluidK(casRegistryNumber = "7664-41-7", description = "Multiparameter: Tillner-Roth. PCSAFT: Mejbri 2005, Cubic: PRMC C.Trujillo 2019. Cp0: Willhoit", molarMass = 1.703026e-02, criticalTemperature = 4.054000e+02, criticalPressure = 1.133300e+07), final onePhase = false, thermoModel = 3, refState = 2);
    end Ammonia;

    package Butane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Butane_n", fluidK(casRegistryNumber = "106-97-8", description = "Multiparameter:Buecker&Wagner 2006. PCSAFT:Gross&Sadowski 2001. Cubic:PRMC, C.Trujillo 2019", molarMass = 0.05812, criticalTemperature = 4.251250e+02, criticalPressure = 3.796000e+06), onePhase = false, thermoModel = 3, refState = 2);
    end Butane_n;

    package Butanol_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Butanol_n", fluidK(casRegistryNumber = "71-36-3", description = "Multiparameter: none. PCSAFT: 2B C.Trujillo 2019. Cubic: PRMC. Cp0: DIPPR107", molarMass = 7.414000e-02, criticalTemperature = 5.630000e+02, criticalPressure = 4.414000e+06), final onePhase = false, thermoModel = 2, refState = 2);
    end Butanol_n;

    package CO2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "CO2", fluidK(casRegistryNumber = "124-38-9", description = "Multiparameter: Span and Wagner 1996. PCSAFT and PRMC: C.Trujillo from GERG2004. Cp0:Wilhoit", molarMass = 4.400980e-02, criticalTemperature = 3.041282e+02, criticalPressure = 7.377730e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end CO2;
    
    package CO2b
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "CO2b",
   fluidK(casRegistryNumber = "124-38-9", description = "SW: GERG2004. PCSAFT and Cubic: C.Trujillo from GERG2004. Cp0:Jaeschke", molarMass = 4.400980e-02, criticalTemperature = 3.041282e+02, criticalPressure = 7.377300e+06),
   final onePhase=false, thermoModel=1, refState=2);
    end CO2b;

    package Dichlorodifluormethane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Dichlorodifluormethane", fluidK(casRegistryNumber = "75-71-8", description = "Multiparameter:Marx et alt. PCSAFT:C.Trujillo 2020, Cubic:PRMC. Cp0:Cooper", molarMass = 1.209130e-01, criticalTemperature = 3.851200e+02, criticalPressure = 4.136100e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Dichlorodifluormethane;

    package EG
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "EG", fluidK(casRegistryNumber = "107-21-1", description = "Multiparameter: none. PPCSAFT_GV: 2B C.Trujillo 2019. Cubic: PRMC. Cp0: Wilhoit", molarMass = 6.207000e-02, criticalTemperature = 7.200000e+02, criticalPressure = 8.200000e+06), final onePhase = false, thermoModel = 2, refState = 2);
    end EG;

    package Eicosane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Eicosane_n", fluidK(casRegistryNumber = "112-95-8", description = "SW: none. PCSAFT: Gross. Cubic: PR. Cp0: Wilhoit ", molarMass = 2.825520e-01, criticalTemperature = 7.680000e+02, criticalPressure = 1.070000e+06), final onePhase = false, thermoModel = 2, refState = 2);
    end Eicosane_n;

    package Ethane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Ethane", fluidK(casRegistryNumber = "74-84-0", description = "Multiparameter: Buecker-Wagner 2006. PCSAAFT: Gross-Sadowski 2001. Cubic:SRKMC Chemsep. Cp0: Cooper", molarMass = 3.006904e-02, criticalTemperature = 3.053220e+02, criticalPressure = 4.872200e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Ethane;
    
    package Ethanol
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Ethanol",
   fluidK(casRegistryNumber = "64-17-5", description = "Multiparameter: Schroeder 2014. PPCSAFT-GV: C.Trujillo 2018. PRMC: C.Trujillo 2018. Cp0: Wilhoit", molarMass = 4.606844e-02, criticalTemperature = 5.147100e+02, criticalPressure = 6.268000e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Ethanol;

    package Heptane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Heptane_n", fluidK(casRegistryNumber = "142-82-5", description = "Multiparameter: Span and Wagner 2003. PCSAFT: Gross 2001. Cubic: PRMC. Cp0: Jaeschke", molarMass = 1.002020e-01, criticalTemperature = 5.401300e+02, criticalPressure = 2.736000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Heptane_n;

    package Hexane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Hexane_n", fluidK(casRegistryNumber = "110-54-3", description = "Multiparameter: Spand and Wagner 2003. PCSAFT: Gross and Sadowski. Cubic: PRMC. Cp0:Jaeschke", molarMass = 8.617536e-02, criticalTemperature = 5.078200e+02, criticalPressure = 3.034000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Hexane_n;

    package Isobutane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Isobutane", fluidK(casRegistryNumber = "75-28-5", description = "Multiparameter: Buecker-Wagner 2006. PCSAFT: Gross-Sadowski 2001. Cubic: PRMC. Cp0: Cooper", molarMass = 5.812220e-02, criticalTemperature = 4.078170e+02, criticalPressure = 3.629000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Isobutane;

    package Methane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Methane", fluidK(casRegistryNumber = "74-82-8", description = "Multiparameter: Seltzmann-Wagner 1991. PCSAF: Gross-Sadowski 2001. Cubic: SRKMC from Chemsep. Cp0: Cooper", molarMass = 1.604280e-02, criticalTemperature = 1.905640e+02, criticalPressure = 4.599200e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Methane;

    package N2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "N2", fluidK(casRegistryNumber = "7727-37-9", description = "Multiparameter: Span 2001, PCSAFT: Gross 2001. PRMC: Barragan 2002. Cp0: DIPPR107", molarMass = 2.801348e-02, criticalTemperature = 1.261920e+02, criticalPressure = 3.395800e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end N2;

    package O2
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "O2", fluidK(casRegistryNumber = "7782-44-7", description = "Multiparameter: Schmidt 1985. PCSAFT: Economou 2007. Cubic:SRKMC Chemsep. Cp0: Cooper", molarMass = 3.199880e-02, criticalTemperature = 1.545810e+02, criticalPressure = 5.043000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end O2;

    package Pentane_n
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Pentane_n", fluidK(casRegistryNumber = "109-66-0", description = "Multiparameter:Span-Wagner 2003. PCSAFT: Gross-Sadowski 2001. Cubic PRMC. Cp0: Wilhoit", molarMass = 7.215000e-02, criticalTemperature = 4.697000e+02, criticalPressure = 3.370000e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Pentane_n;

    package Propane
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Propane", fluidK(casRegistryNumber = "74-98-6", description = "Multiparameter: Lemmon 2009. PCSAFT: Gross-Sadowski 2001. Cubic: SRKMC Chemsep. Cp0:Cooper", molarMass = 4.409562e-02, criticalTemperature = 3.698900e+02, criticalPressure = 4.251200e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end Propane;

    package R134A
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "R134A", fluidK(casRegistryNumber = "811-97-2", description = "Multiparameter:Tillner 1994. PCSAFT: non assoc.C.Trujillo 2019. Cubic:PRMC, C.Trujillo 2019. Cp0: Wilhoit", molarMass = 1.020310e-01, criticalTemperature = 3.741800e+02, criticalPressure = 4.901200e+06), final onePhase = false, thermoModel = 3, refState = 1);
    end R134A;

    package R410A
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "R410A", fluidK(casRegistryNumber = "R410A.PPF", description = "Multiparameter: Lemmon 2003. PCSAFT: non assoc. C.T.2019. Cubic: none. Cp0: Cooper", molarMass = 7.258540e-02, criticalTemperature = 3.444940e+02, criticalPressure = 4.901200e+06), final onePhase = false, thermoModel = 3, refState = 2);
    end R410A;
    
    package Toluene
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Toluene",
   fluidK(casRegistryNumber = "108-88-3", description = "Multiparameter: Lemmon 2006. PCSAFT: Gross 2001. Cubic: SRKMC. Cp0: Wilhoit", molarMass = 9.214020e-02, criticalTemperature = 5.917500e+02, criticalPressure = 4.126000e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Toluene;

    package WaterRef
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "WaterRef", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:IAPWS95. PCSAFT:C.Trujillo 2B max.633K. Cubic:PRMC,C.Trujillo 2018. Cp0:Jaeschke", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), final onePhase = false, thermoModel = 3, refState = 4);
    end WaterRef;

    package Water
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Water", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:GERG2004. PCSAFT: Diamantotis 4C. Cubic: SRKMC C.Trujillo. Cp0:Jaechske", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), onePhase = false, thermoModel = 3, refState = 4, reference_T = 273.15, reference_p = 101325);
    end Water;
    
    package Xylene_m
      extends FreeFluids.ExternalMedia.ExternalMedium(final mediumName = "Xylene_m",
   fluidK(casRegistryNumber = "108-38-3", description = "Multiparameter: Zhou 2012. PCSAFT: Gross 2001. Cubic: PRSV Proust 1998. Cp0: Cooper", molarMass = 1.061670e-01, criticalTemperature = 6.168900e+02, criticalPressure = 3.534600e+06),
   final onePhase=false, thermoModel=3, refState=2);
    end Xylene_m;
  end Fluids;

  package Tests
    model Test1aCubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.WaterRef(thermoModel = 1, refState = 4, reference_T = 273.15, reference_p = 1.0e5, inputChoice = "pT"), p = 220.0e5, initialT = 600, finalT = 700);
    end Test1aCubic;

    model Test1aPCSAFT
      extends Test1aCubic(Medium(thermoModel = 2));
    end Test1aPCSAFT;

    model Test1aSW
      extends Test1aCubic(Medium(thermoModel = 3));
    end Test1aSW;

    model Test1aIF97
      extends Test1aCubic(redeclare package Medium = Modelica.Media.Water.StandardWater);
      end Test1aIF97;

    model Test1bCubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.R410A(thermoModel = 1, refState = 2, reference_T = 273.15, reference_p = 101325.0, inputChoice = "pT"), p = 10.0e5, initialT = (-50) + 273.15, finalT = 100 + 273.15);
    end Test1bCubic;

    model Test1bPCSAFT
      extends Test1bCubic(Medium(thermoModel = 2));
    end Test1bPCSAFT;

    model Test1bSW
      extends Test1bCubic(Medium(thermoModel = 3));
    end Test1bSW;
    
    model Test1bTMedia
      extends Test1bCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.R410A(refState = "IIR", highPressure = true, inputChoice="pT"));
      end Test1bTMedia;
      
    model Test1cCubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.Propane(thermoModel = 1, refState = 2, reference_T = 273.15, reference_p = 101325.0, inputChoice = "pT"), p = 5.0e5, initialT = (-100) + 273.15, finalT = 200 + 273.15);
    end Test1cCubic;
    
    model Test1cPCSAFT
      extends Test1cCubic(Medium(thermoModel = 2));
    end Test1cPCSAFT;
  
    model Test1cSW
      extends Test1cCubic(Medium(thermoModel = 3));
    end Test1cSW;
    
    model Test1cTMedia
      extends Test1cCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.Propane(refState = "IIR", highPressure = true, inputChoice="pT"));
      end Test1cTMedia;
      
    model Test1dCubic
      extends FreeFluids.TMedia.Tests.FluidTestingA(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.Ethanol(thermoModel = 1, refState = 2, reference_T = 273.15, reference_p = 101325.0, inputChoice = "pT"), p = 1.0e5, initialT = 273, finalT = 300 + 273.15);
    end Test1dCubic;
    
    model Test1dPCSAFT
      extends Test1dCubic(Medium(thermoModel = 2));
    end Test1dPCSAFT;
  
    model Test1dSW
      extends Test1dCubic(Medium(thermoModel = 3));
    end Test1dSW;
    
    model Test1dTMedia
      extends Test1dCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.Ethanol(refState = "IIR", highPressure = true, inputChoice="pT"));
      end Test1dTMedia;

    model TestRefrigerantEvaporator "The test case has been copied from ThermoPower Library. Test case with once-through evaporator using Dittus-Boelter 2-phase heat transfer model. The initial convergence takes some minuts, but should work with thermoModel=1(Cubic) or 2 (PCSAFT)"
      replaceable package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1) constrainedby Modelica.Media.Interfaces.PartialMedium;
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
      ThermoPower.Water.SinkPressure fluidSink(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1), p0 = 3e6) annotation(
        Placement(transformation(extent = {{50, -82}, {70, -62}}, rotation = 0)));
      ThermoPower.Gas.SinkPressure gasSink(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{-54, 18}, {-74, 38}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow fluidSource(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1), h = -100e3, p0 = 3e5, use_in_h = false, use_in_w0 = true) annotation(
        Placement(transformation(extent = {{-70, -82}, {-50, -62}}, rotation = 0)));
      ThermoPower.Water.SensT fluid_T_in(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1)) annotation(
        Placement(transformation(extent = {{-46, -78}, {-26, -58}}, rotation = 0)));
      ThermoPower.Gas.SensT gas_T_in(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{34, 22}, {14, 42}}, rotation = 0)));
      ThermoPower.Gas.SourceMassFlow gasSource(redeclare package Medium = ThermoPower.Media.Air, T = 380.15, use_in_w0 = true) annotation(
        Placement(transformation(extent = {{64, 18}, {44, 38}}, rotation = 0)));
      ThermoPower.Water.SensT fluid_T_out(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1)) annotation(
        Placement(transformation(extent = {{20, -78}, {40, -58}}, rotation = 0)));
      ThermoPower.Gas.SensT gas_T_out(redeclare package Medium = ThermoPower.Media.Air) annotation(
        Placement(transformation(extent = {{-24, 22}, {-44, 42}}, rotation = 0)));
      inner ThermoPower.System system(allowFlowReversal = false, initOpt = ThermoPower.Choices.Init.Options.steadyState) annotation(
        Placement(transformation(extent = {{80, 80}, {100, 100}})));
      ThermoPower.Gas.Flow1DFV gasFlow(redeclare package Medium = ThermoPower.Media.Air, redeclare model HeatTransfer = ThermoPower.Thermal.HeatTransferFV.ConstantHeatTransferCoefficient(gamma = 120), A = (Dgas ^ 2 - Dext ^ 2) / 4 * pi, Dhyd = Dgas, FFtype = ThermoPower.Choices.Flow1D.FFtypes.NoFriction, HydraulicCapacitance = ThermoPower.Choices.Flow1D.HCtypes.Downstream, L = L, N = Nnodes, Nt = 1, Tstartin = 273.15, Tstartout = 283.15, noInitialPressure = true, omega = Dext * pi, wnom = 0.02) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-8, 28})));
      ThermoPower.Water.Flow1DFV2ph fluidFlow(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.CO2b(thermoModel = 2, refState = 1), redeclare model HeatTransfer = ThermoPower.Thermal.HeatTransferFV.FlowDependentHeatTransferCoefficient2ph(gamma_nom_liq = 800, gamma_nom_2ph = 8000, gamma_nom_vap = 800), A = Dint ^ 2 / 4 * pi, Cfnom = 0.005, Dhyd = Dint, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, HydraulicCapacitance = ThermoPower.Choices.Flow1D.HCtypes.Downstream, L = L, N = Nnodes, dpnom = 1000, hstartin = 250e3, hstartout = 250e3, noInitialPressure = true, omega = Dint * pi, pstart = 3e6, wnom = 0.005) annotation(
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
    <p>The medium is designed for pure or pseudo-pure substances in liquid, gas, or two phases. The thermodynamic calculations are performed using different equations of state(EOS), implemented in C language, that are called using external object and functions. The C code and the substances data are in the Resources folder.</p>
    <p>For each substance there is the possibility of choosing between three EOS types: multiparameter, PCSAFT, or cubic. The external object contains the substance data, and is a C structure, that is exported from the database using the FreeFluids GUI software.</p>
    <p>The quality of the results is very high when multiparameter (Seltzmann and Wagner, SW) EOS are available. If not, PCSAFT or different flavours of cubic EOS can be used.</p>
    <p>The transport properties are normally computed from temperature dependent correlations, with pressure correction. Nevertheless, for viscosity, the system will perform phase independent viscosity calculation (from temperature and density) for selected substances, using dedicated calculation or extended correspondent states with NIST correction polynomia, when available. In order to use the phase independent viscosity calculation the thermoModel parameter must be 3 (multiparameter EOS), as it is the only way to be sure that the supplied density is reliable. For thermal conductivity the same type of calculations are already programmed, but have been inactivated till finishing testing.</p>
    <p>The medium implements all the requirements of Modelica.Media.PartialTwoPhaseMedium. It is compatible with the old frontend of OpenModelica, but not with the new frontend, as it does not support external functions yet. It is also compatible, at least partially, with the ThermoPower library. Compared with other high quality libraries, it is very easy to use and, although it is not as complete as for example CooProp, or RefProp, the quality of the thermodynamic results is practically the same. In plus there is the possibility of using PCSAFT and cubic EOS for substances for which no multiparameter EOS are not available.</p>
    <p>The main limitation of the medium, when using SW or PCSAFT EOS, is in the vecinity of the critical point, as no special technique, as for example splines, is used. Nevertheless you can go normally quite close to the critical point with SW, not so with PCSAFT as there is still some problems with its density solver. The medium is much slower than the TMedia one, but is the price for the better precision and wider application.</p>
    <p>The ThermodynamicState record is quite big, in order to allow all the thermodynamic calculations with just one call to external functions.</p>
    <p>The BaseProperties model has been implemented using algorithms instead of equations. I do not know if this is correct, but seems the logical decision when you know the order in which calculations must be performed.</p>
    <p>When extending the medium, the following configuration must be done:</p>
    <p>The thermoModel Integer constant must be made equal to 1 for using the cubic EOS, 2 for the PCSAFT, or 3 for the SW.</p>
    <p>The referenceState Integer constant  must be made equal to 1 for using ASHRAE reference state, to 2 for using IIR, to 3 for using NBP, to 4 for using user reference_T and reference_p. Any other number will produce an unreferenced calculation for enthalpy and entropy.</p>
    <p>If you are going to use the BaseProperties model, you must make the inputChoice constant String equal to 'ph', 'pT' or 'dT', as needed. When instantiating the model, you can still change the selection just for the object, specifying the value for the localInputChoice constant.</p>
    <p>With the Test1aXXX models of the Tests package, you can compare the performance of the cubic, PCSAFT, IAPWS95 and GERG2004 EOS, against the standard Modelica.Media.Water.StandardWater. Test1bXXX compare R410A, Test1cXXX compare propane, and Test1dXXX compare ethanol.</p>
    <p><b>The C code</b></p>
    <p>The C code is placed in the Resources folder.</p>
    <p>The FFbasic.h file contains the basic definition of structures and enumerations used in the code. </p>
    <p>The FFmodelicaMedium.c file contains the interface between Modelica and C. Basically the constructor and destructor of the structure that will contain the data of the selected substance. The reference to this structure is passed to the Modelica code, that uses this reference when calling the external functions written in C.</p>
    <p>The FFeosPure.c and FFphysprop.c files contain the code that will perform the calculation in C.</p>
    <p><b>The Medium data</b></p>
    <p>We need to define the mediums to be used. This definition is inside the Fluids subpackage. Each definition is very short, just with the minimum information necessary. The most important of it is the mediumName, as this name is the name of the file (with extension .sd) placed in the Resources folder that will be used for charging the data to the C structure.</p>
    <p></p>
    <p>For each medium you need both a medium definition in Modelica, and a binary file with the data in the Resources folder. Both can be made using the program FreeFluidsGui.exe that has been also placed in the Resources folder. It will access the data base, allow you to select the substance, with the EOS and correlations to be used, and later export it to a binary file, with the data as C structure, and to a text file, with the data to add to the FreeFluids.ExternalMedia.Fluids package.</p>
    <p>The program allows also the exportation of the correlations in the format used by the TMedia package.</p>
    <p><b>Transport properties by dedicated calculation or ECS</b></p>
    <p>The calculation of the substance viscosity is managed by the FF_Viscosity function in the FFphysprop.c file. If you have selected to work with the multiparameter EOS (thermoModel=3), that grants a good density calculation, it will check if there is a dedicated calculation from temperature and density (available only for few substances). If not, it will check if the correlation defined for gas viscosity calculation is the 112 (NIST coefficients for viscosity calculation using ECS) and that the data charged for the reference substance correspond to the needed one. If this is OK the ECS calculation will be applied for viscosity, otherwise temperature dependent correlations, with pressure correction, will be used.</p>
    <p>If the EOS is of the cubic or PCSAFT type, correlations will be used if available. If not, the Lucas approximation will be used for gas. For liquid phase, if there is no correlation defined, the calculation will be done using ECS from the reference liquid defined.</p>
    <p>The definition of the reference liquid is done at medium package level, giving value to the const String refName, which default value is 'Propane'.</p>
    </body>
    </html>"));
end ExternalMedia;

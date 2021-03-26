within FreeFluids;

package TMedia "TMedia.mo by Carlos Trujillo
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
  partial package TMedium
    extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
    singleState = false,
    onePhase = false,
    smoothModel=true,
    ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph,
    Temperature(min=fluidConstants[1].lCpLimI, max=if onePhase== true then fluidConstants[1].lCpLimS else fluidConstants[1].Cp0LimS, start=T_default, displayUnit="degC"),
    AbsolutePressure(min=Modelica.Constants.small, max=200.0e5, start=p_default, displayUnit="bar"),
    Density(displayUnit="kg/m3"),
    SpecificEntropy(min=-Modelica.Constants.inf, max=Modelica.Constants.inf),
    T_default=0.5*(fluidConstants[1].lCpLimS+fluidConstants[1].lCpLimI),
    p_default=101325,
    reference_T = 298.15,
    redeclare record FluidConstants = FreeFluids.MediaCommon.DataRecord);
    import PhysPropCorr = FreeFluids.MediaCommon.Functions.PhysPropCorr;
    import SpecificEnthalpyCorr = FreeFluids.MediaCommon.Functions.SpecificEnthalpyCorr;
    import SpecificEnthalpyCorrInv = FreeFluids.MediaCommon.Functions.SpecificEnthalpyCorrInv;
    import SpecificEntropyCorr = FreeFluids.MediaCommon.Functions.SpecificEntropyCorr;
    import SpecificEntropyCorrInv = FreeFluids.MediaCommon.Functions.SpecificEntropyCorrInv;
    constant String refState = "None" "Enthalpy/entropy reference state. Alternatives: ASHRAE, IIR, NBP, User(value at reference_T as 0), otherwise no reference";
    //constant FreeFluids.MediaCommon.Types.ReferenceState refState=FreeFluids.MediaCommon.Types.ReferenceState.None  "Enthalpy/entropy reference state. Alternatives: ASHRAE, IIR, NBP, User(value at reference_T as 0), otherwise no reference";
    constant String inputChoice = "ph" "Allows to choose the input choise to use for the construction of a BaseProperties object. Alternative are: pT, ph, dT";
    constant Boolean highPressure = false;
    /*constant SpecificEnthalpy reference_h= referenceEnthalpy();
    constant SpecificEntropy reference_s= referenceEntropy();
    constant SpecificEnthalpy critical_h = SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tc - 0.01);
    constant SpecificEnthalpy critical_h0 = SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, fluidConstants[1].Tc-0.01);*/
    
    //Auxiliary functions based in correlations
    //-----------------------------------------
  
    function referenceEnthalpy
      output SpecificEnthalpy reference_h;
    algorithm
      reference_h := if refState == "ASHRAE" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 2.0e5 else if refState == "NBP" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
    end referenceEnthalpy;
  
    function referenceEntropy
      output SpecificEntropy reference_s;
    algorithm
      reference_s := if refState == "ASHRAE" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 233.15) else if refState == "IIR" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, 273.15) - 1.0e3 else if refState == "NBP" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, fluidConstants[1].Tb) else if refState == "User" then SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, reference_T) else 0;
    end referenceEntropy;
      
    redeclare function extends saturationPressure "Return saturation pressure from T"
        extends Modelica.Icons.Function;
  
      algorithm
        //assert(T <= fluidConstants[1].Tc, "The media can´t be used over Tc");
        if T>=fluidConstants[1].Tc then
          p :=  fluidConstants[1].criticalPressure;
        elseif T > fluidConstants[1].VpLimS then
          assert(false, "T for saturation pressure: "+String(T)+" over Vp limit T") "Between upper limit of correlation and Tc we can't work";
        else
          p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
        end if;
    end saturationPressure;
    
    function saturationPressure_derT "Return derivative of saturation pressure w.r.t. temperature"
      input Temperature T;
      output Real dpT;
    algorithm
      if T<fluidConstants[1].Tc then
        dpT := (saturationPressure(T) - saturationPressure(T-0.01)) / 0.01;
      else
        dpT:=0;
      end if;
    end saturationPressure_derT;
  
    redeclare function extends saturationTemperature "Return saturation temperature from P"
        extends Modelica.Icons.Function;
  
      algorithm
        //assert(p < fluidConstants[1].criticalPressure, "No saturated media possible over Pc", AssertionLevel.warning);
        T := if p >= fluidConstants[1].criticalPressure then fluidConstants[1].Tc else FreeFluids.MediaCommon.Functions.PhysPropCorrInv(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, p, fluidConstants[1].VpLimS);
    end saturationTemperature;
  
    redeclare function extends saturationTemperature_derp "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
  
    algorithm
      if p<fluidConstants[1].criticalPressure then
        dTp := (saturationTemperature(p) - saturationTemperature(0.999 * p)) / (0.001 * p);
      else
        dTp:=0;
      end if;
    end saturationTemperature_derp;
  
  
    redeclare model extends BaseProperties(h(stateSelect = if preferredMediumStates and localInputChoice == "ph" then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates and (localInputChoice == "ph" or localInputChoice == "pT") then StateSelect.prefer else StateSelect.default), T(stateSelect = if preferredMediumStates and (localInputChoice == "pT" or localInputChoice == "dT") then StateSelect.prefer else StateSelect.default), d(stateSelect = if preferredMediumStates and localInputChoice == "dT" then StateSelect.prefer else StateSelect.default))
        parameter String localInputChoice = inputChoice;
  
      protected
        Real bubH, hv, dewH;
  
      equation
        MM = fluidConstants[1].MW / 1000 "in kg/mol";
        R = Modelica.Constants.R / MM;
        if localInputChoice == "ph" then
          T = temperature_ph(p=p, h=h);
          d = density_ph(p=p, h=h);
        elseif localInputChoice == "pT" then
          d = density_pT(p=p, T=T);
          h = specificEnthalpy_pT(p=p, T=T);
        elseif localInputChoice == "dT" then
          p = pressure_dT(d=d, T=T);
          h = specificEnthalpy_dT(d=d, T=T);
        else
          assert(false, "Invalid choice for BaseProperties inputChoice");
        end if;
        u = h - p / d;
        sat = setSat_p(p=p);
        bubH = bubbleEnthalpy(sat);
        hv = PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, sat.Tsat);
        dewH = bubH + hv;
        if p < fluidConstants[1].criticalPressure and h > bubH and h < dewH then
          state.phase = 2;
          state.gf = (h - bubH) / hv;
        else
          state.phase = 1;
          if h <= bubH then
            state.gf = 0.0;
          else
            state.gf = 1.0;
          end if;
        end if;
        state.p = p;
        state.T = T;
        state.h = h;
        state.d = d;
    end BaseProperties;
  
    //Thermodynamic state definition and constructors
    //-----------------------------------------------
  
    redeclare record extends ThermodynamicState(phase(start=0))
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
 
      algorithm
        state.p := p;
        state.T := T;
        Vp := saturationPressure(T);
        if p > Vp then
          state.phase := 1 "liquid";
          state.gf := 0;
          assert(T<=fluidConstants[1].lCpLimS, "Liquid at temperature over the limit of Cp correlation ("+String(fluidConstants[1].lCpLimS)+"). pTX1");
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
          assert((Tb <= fluidConstants[1].lCpLimS) and (Tb <= fluidConstants[1].HvLimS), "Gas at condensation temperature over correlations limits for enthalpy calculation. pTX2");
          DgsP := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
          Dgi := p * fluidConstants[1].MW / (1000 * R * T);
          state.d:=Dgi+(DgsP-Dgi)*(Tb/T)^9;
          state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
          state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
          state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1 + 1.2e-7 * p*(fluidConstants[1].Tc/T)^3) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
        else
          assert(false, "A two phases state can't be constructed from p and T. pTX3");
        end if;
        state.h := state.h - referenceEnthalpy() "referenced enthalpy";
    end setState_pTX;
  
    redeclare function extends setState_dTX "Return ThermodynamicState record as function of T,d and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        Density Dls;
        Density Dgs;
        SpecificEnthalpy Hls "enthalpy of the saturated liquid at given T";
        AbsolutePressure Vp "vapor pressure at given T";
        Temperature Tb "boiling temperature at tested pressure";
        Density Dig "ideal gas density at T and p";
        Density Dsim "retrieved density at tested pressure";
        Density DsimN "retieved density at slightly higher pressure";
        Real Error "density error";
        Real Slope "density change with pressure";
        Integer i;
  
      algorithm
        //assert(T < fluidConstants[1].Tc, "The media can´t be used over Tc");
        state.T := T;
        state.d := d;
        Vp := saturationPressure(T);
        if T<fluidConstants[1].Tc then
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, T) "saturated liquid density at T";
        else
          Dls := fluidConstants[1].molarMass/fluidConstants[1].Vc "saturated liquid density at T";
        end if;
        if d >= Dls then
          state.phase := 1 "liquid";
          state.gf := 0;
          assert(T<=fluidConstants[1].lCpLimS, "Liquid at temperature over the limit of Cp correlation dTX1");
          Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy"; 
          state.p := if fluidConstants[1].lnuA > 0 then Vp + (exp((d - Dls) * fluidConstants[1].lnuA) - 1) * 1000 * Modelica.Constants.R * T * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) / (fluidConstants[1].lnuA * fluidConstants[1].MW) else (1 - Dls / d) / fluidConstants[1].IsothComp + Vp;
          state.h := if T <= 0.45 * fluidConstants[1].Tc then Hls + (state.p - Vp) / state.d else if T < 0.85 * fluidConstants[1].Tc then Hls + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p - Vp) / state.d else Hls;
        else
          if T<fluidConstants[1].Tc then
            Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, T);
          else
            Dgs := Dls;
          end if;
          if d >= Dgs then
            assert(T<=fluidConstants[1].lCpLimS, "Liquid at temperature over the limit of Cp correlation  dTX2");
            Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, T) "unreferenced saturated liquid enthalpy";
            state.p := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, T);
            state.h := if T <= 0.45 * fluidConstants[1].Tc then Hls + state.p / state.d - Vp / Dls else if T < 0.85 * fluidConstants[1].Tc then Hls + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else Hls "pressure correction for liquid enthalpy";
            if d > Dgs then
              state.phase := 2;
              state.gf := Dgs * (Dls - d) / (d * (Dls - Dgs));
              state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid enthalpy";
              state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, T);
            elseif d == Dgs then
              state.phase := 1 "saturated gas";
              state.gf := 1;
              state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, T);
            end if;
          else
            state.phase := 1 "gas phase";
            state.gf := 1;
            Tb:=fluidConstants[1].lCpLimS/2;
            state.p:=saturationPressure(Tb) "pressure initialization at half of max. lCp)";
            Dig := state.p * fluidConstants[1].MW / (1000 * R * T);
            Dgs:=PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
            Dsim:=Dig+(Dgs-Dig)*(Tb/T)^9;
            Error:=Dsim-state.d;
            while abs(Error/state.d)>0.0001 loop
              Dig:= Dig*1.001;
              Tb:=saturationTemperature(state.p * 1.001);
              Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
              DsimN:=Dig+(Dgs-Dig)*(Tb/T)^9;
              Slope := (DsimN - Dsim) / (0.001 * state.p);
              state.p := state.p - Error/ Slope;
              Dig := state.p * fluidConstants[1].MW / (1000 * R * T);
              Tb:=saturationTemperature(state.p);
              Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
              Dsim:=Dig+(Dgs-Dig)*(Tb/T)^9;
              Error:=Dsim-state.d;
            end while;
            Tb := saturationTemperature(state.p);
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
            assert(Tb<=fluidConstants[1].HvLimS, "Liquid at temperature over the limit of Hv correlation dTX3");
            state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
            state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1 + 1.2e-7 * state.p*(fluidConstants[1].Tc/T)^3) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
          end if;
        end if;
      state.h := state.h - referenceEnthalpy() "referenced enthalpy";
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
        SpecificEnthalpy H0b "ideal gas enthalpy at boiling p";
        SpecificEnthalpy H0t "ideal gas enthalpy at state.T";
        Real error;
  
      algorithm
        state.p := p;
        state.h := h;
        Tb := saturationTemperature(p);
        Hraw := h + referenceEnthalpy();
        assert(Tb<=fluidConstants[1].lCpLimS, "Boiling temperature("+String(Tb)+") for pressure "+String(p)+" over the limit("+String(fluidConstants[1].lCpLimS)+") of Cp correlation  phX1", AssertionLevel.warning);
        Hls := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific enthalpy at Tb";
        if Hls >= Hraw then
          state.phase := 1 "liquid";
          state.gf := 0;
          state.T := SpecificEnthalpyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Hraw, fluidConstants[1].lCpLimS);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
          Vp := saturationPressure(state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          if highPressure == true then
            Hraw := if state.T <= 0.45 * fluidConstants[1].Tc then Hraw - state.p / state.d + Vp / Dls else if state.T < 0.85 * fluidConstants[1].Tc then Hraw - (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else Hraw "pressure correction for liquid enthalpy";
            state.T := SpecificEnthalpyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Hraw, fluidConstants[1].lCpLimS);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
            Vp := saturationPressure(state.T);
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
            state.T:=fluidConstants[1].Tc "begin search at Tc";
            H0t:=SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
            error:=(H0t-H0b)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^3)+Hgs-Hraw;
            
            while abs(error/Hraw)>0.00001 loop
              state.T := state.T-error/(PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^3)-3.6e-7*state.p*fluidConstants[1].Tc^3/state.T^4*(H0t-H0b));
              H0t:=SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
              error:=(H0t-H0b)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^3)+Hgs-Hraw;
            end while;
            Dgi := p * fluidConstants[1].MW / (1000 * R * state.T);
            state.d:=Dgi+(Dgs-Dgi)*(Tb/state.T)^9;
          end if;
        end if;
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
        SpecificEntropy S0b;
        SpecificEntropy S0t;
        SpecificEnthalpy Href;
        AbsolutePressure Vp;
        Real Error;
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
        assert(Tb<=fluidConstants[1].lCpLimS, "Boiling temperature over the limit of Cp correlation  psX1", AssertionLevel.warning);
        Sl := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid specific entropy at T";
        if Sl >= Sraw then
          state.phase := 1 "liquid";
          state.gf := 0;
          state.T := SpecificEntropyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Sraw, fluidConstants[1].lCpLimS - 0.01);
          Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
          Vp := saturationPressure(state.T);
          state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          if highPressure == true then
            Sraw := if state.T <= 0.45 * fluidConstants[1].Tc then Sraw else Sraw + (state.T - 0.45 * fluidConstants[1].Tc) * (state.p - Vp) / (0.4 * fluidConstants[1].Tc * state.d * state.T) "pressure correction for liquid entropy";
            state.T := SpecificEntropyCorrInv(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Sraw, fluidConstants[1].lCpLimS - 0.01);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) "saturated liquid density";
            Vp := saturationPressure(state.T);
            state.d := if fluidConstants[1].lnuA > 0 then Dls + log(fluidConstants[1].lnuA * fluidConstants[1].MW / (1000 * exp(fluidConstants[1].lnuA * Dls + fluidConstants[1].lnuB) * Modelica.Constants.R * state.T) * (state.p - Vp) + 1) / fluidConstants[1].lnuA else Dls / (1 - fluidConstants[1].IsothComp * (state.p - Vp)) "saturated liquid density, pressure corrected";
          end if;
          state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid enthalpy";
          if highPressure == true then
            state.h := if state.T <= 0.45 * fluidConstants[1].Tc then state.h + state.p / state.d - Vp / Dls else if state.T < 0.85 * fluidConstants[1].Tc then state.h + (0.85 * fluidConstants[1].Tc - state.T) / (0.4 * fluidConstants[1].Tc) * (state.p / state.d - Vp / Dls) else state.h "pressure correction for liquid enthalpy";
          end if;
        else
          Sg := Sl + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) / Tb;
          Dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, Tb);
          if Sg >= Sraw then
            state.phase := 2 "biphasic";
            state.T := Tb;
            state.gf := (Sraw - Sl) / (Sg - Sl);
            Dls := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, Tb);
            state.d := Dls * Dgs / (state.gf * Dls + (1 - state.gf) * Dgs);
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid enthalpy";
            state.h := state.h + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "unreferenced mix enthalpy";
          else
            state.phase := 1 "gas";
            state.gf := 1.0;
            S0b := SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb);
            state.T:=fluidConstants[1].Tc "begin search at Tc";
            S0t:=SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
            Error:=(S0t-S0b)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^4)+Sg-Sraw;
            while abs(Error/Sraw)>0.00001 loop
              state.T := state.T-Error/(PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^4)/state.T-4.8e-7*state.p*fluidConstants[1].Tc^4/state.T^5*(S0t-S0b));
              S0t:=SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
              Error:=(S0t-S0b)*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^4)+Sg-Sraw;
            end while;
            Dgi := p * fluidConstants[1].MW / (1000 * R * state.T);
            state.d:=Dgi+(Dgs-Dgi)*(Tb/state.T)^9;
            state.h := SpecificEnthalpyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreference saturated liquid enthalpy at p";
          state.h := state.h + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) "saturated gas enthalpy at p";
          state.h := state.h + (SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T) - SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1 + 1.2e-7 * p*(fluidConstants[1].Tc/state.T)^3) "Gas entahlpy at T and p. Adjustement from Tb to T using Cp0";
          end if;
        end if;
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
  
    redeclare function extends density_derp_T "Return density derivative w.r.t. pressure at const temperature. OK"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.phase == 1 then
          ddpT := state.d * isothermalCompressibility(state);
        else
          ddpT := Modelica.Constants.inf;
        end if;
    end density_derp_T;
  
    redeclare function extends density_derT_p "Return density derivative w.r.t. temperature at constant pressure. OK"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.phase == 1 then
          ddTp := -state.d * isobaricExpansionCoefficient(state);
        else
          ddTp := Modelica.Constants.inf;
        end if;
    end density_derT_p;
  
    redeclare function extends density_derp_h "Return density derivative w.r.t. pressure at const specific enthalpy. Numerical derivative in two phases differs from HelmholtzMedia, but is similar to ExternalMedia"
        extends Modelica.Icons.Function;
  
      protected
        Real dvt "(dV/dT)P";
        Real dvp "(dV/dT)P";
  
      algorithm
        if state.phase == 1 then
          dvt := isobaricExpansionCoefficient(state) / state.d;
          dvp := -isothermalCompressibility(state) / state.d;
          ddph := -state.d * state.d * (dvp + (state.T * dvt * dvt - dvt / state.d) / specificHeatCapacityCp(state));
        else
          ddph := (state.d- density(setState_phX(0.98 * state.p, state.h))) / (0.02 * state.p);
        end if;
    end density_derp_h;
  
    redeclare function extends density_derh_p "Return density derivative w.r.t. specific enthalpy at constant pressure. OK"
        extends Modelica.Icons.Function;
  
      algorithm
        if state.phase == 1 then
          ddhp := -state.d * isobaricExpansionCoefficient(state) / specificHeatCapacityCp(state);
        else
          ddhp := -state.d * state.d * (state.T-saturationTemperature(state.p * 0.95)) / (state.p * 0.05 * state.T);
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
        Temperature Tb;
  
      algorithm
        if state.gf <= 0 then
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid entropy at T";
          if highPressure == true then
            Vp := PhysPropCorr(fluidConstants[1].VpCorr, fluidConstants[1].VpCoef, fluidConstants[1].MW, state.T);
            s := if state.T <= 0.45 * fluidConstants[1].Tc then s else s - (state.T - 0.45 * fluidConstants[1].Tc) * (state.p - Vp) / (0.4 * fluidConstants[1].Tc * state.d * state.T) "pressure correction for liquid entropy";
          end if;
        elseif state.gf < 1 then
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T) "unreferenced saturated liquid entropy at T";
          s := s + state.gf * PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, state.T) / state.T "mix liquid-vapor entropy";
        else
          Tb := saturationTemperature(state.p);
          s := SpecificEntropyCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, Tb) "unreferenced saturated liquid entropy at p(Tb)";
          s := s + PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, Tb) / Tb "saturated gas entropy at p(Tb)";
          s:= s + (SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T) - SpecificEntropyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) * (1 + 1.2e-7 * state.p * (fluidConstants[1].Tc / state.T)^4);
        end if;
        s := s - referenceEntropy() "referenced entropy";
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
        if state.gf == 0.0 then
          cp := PhysPropCorr(fluidConstants[1].lCpCorr, fluidConstants[1].lCpCoef, fluidConstants[1].MW, state.T);
          if highPressure == true then
            ds := PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T);
            cp := if fluidConstants[1].family == 17 then cp * exp(-5.0 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5)) else cp * exp(-2.8 * (state.d - ds) ^ 0.5 / (ds * (1 - state.T / fluidConstants[1].Tc) ^ 0.5));
          end if;
        elseif state.gf == 1.0 then
          cp := PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
          Tb := saturationTemperature(state.p);
          //cp := cp * (1 + 3.0e-7 * state.p ^ 1.001 * (1 - (state.T - Tb) / (fluidConstants[1].Tc - Tb))^2);
          //cp:=cp*(1+1.5e-7*state.p) "alternative equivalent to enthalpy calculation";
          //cp:=cp*(1+4.0e-8*state.p*(fluidConstants[1].Tc/state.T)^6);
          cp:=cp*(1+1.2e-7*state.p*(fluidConstants[1].Tc/state.T)^3)-3*1.2e-7*state.p*fluidConstants[1].Tc^3/state.T^4*(SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T)-SpecificEnthalpyCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, Tb)) "derivative of enthalpy";
        else
          cp := Modelica.Constants.small;
        end if;
        
    end specificHeatCapacityCp;
  
    redeclare function extends specificHeatCapacityCv "Return heat capacity at constant volume"
        extends Modelica.Icons.Function;
  
      protected
        Real beta, kappa, cp0;
  
      algorithm
        if state.gf == 0.0 then
          beta := isobaricExpansionCoefficient(state);
          kappa := isothermalCompressibility(state);
          cv := specificHeatCapacityCp(state) - state.T * beta * beta / (kappa * state.d);
        elseif state.gf == 1.0 then
          cp0:=PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
          cv := specificHeatCapacityCp(state)*(cp0-R * 1000/ fluidConstants[1].MW)/cp0 "same relationship as ideal heat capacities";
        else
          cv := Modelica.Constants.small;
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
        if (state.phase==1) then
          gamma := cp / cv;
        else
          gamma:=Modelica.Constants.inf;
        end if;
    end isentropicExponent;
  
    redeclare function extends isobaricExpansionCoefficient "Returns the approximate isobaric expansion coefficient beta=dV_dT/V"
        extends Modelica.Icons.Function;
  
      protected
        Real dd, d0, d1, d2, dv, vp1, vp2;
        Real vp "vapor pressure at T";
        Real tb "boiling temperature at p";
        Real dgi;
        Real dgs;
  
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
          dgi := state.p * fluidConstants[1].MW / (1000 * R * (state.T));
          tb := saturationTemperature(state.p);
          dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, tb);
          beta:=-(-state.p * fluidConstants[1].MW / (1000 * R * (state.T)^2)*(1-(tb/state.T)^9)-9*(dgs-dgi)*tb^9/state.T^10)/state.d;
          /*
          vp := saturationPressure(state.T);
          if state.p > 0.1 * vp then
            vp := saturationPressure(state.T + 0.01);
            dgi := state.p * fluidConstants[1].MW / (1000 * R * (state.T + 0.01));
            tb := saturationTemperature(state.p);
            dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, tb);
            beta := -(dgi + (dgs - dgi) / (0.9 * vp) * (state.p - 0.1 * vp) - state.d) / (0.01 * state.d);
          else
            beta := 1 / state.T "ideal gas expansion coefficient";
          end if;
          */
          //beta := (1 / state.T)*(1+0.4*(tb/state.T)^6) "ideal gas expansion coefficient";
        else
          beta := Modelica.Constants.small;
        end if;
  //at saturation
    end isobaricExpansionCoefficient;
  
    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
  
      protected
        Real rhoRed1 "reduced density - 1";
        Real vp "vapor pressure at T";
        Real tb "boiling temperature at p";
        Real dgi;
        Real dgs;
        Real dgsMinus;
        Real dDgs;
        Real dTb;
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
          tb := saturationTemperature(state.p);
          dgi := state.p * fluidConstants[1].MW / (1000 * R * state.T);
          dgs := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, saturationTemperature(state.p));
          dgsMinus := PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, saturationTemperature(0.999 * state.p));
          dDgs:=(dgs-dgsMinus)/(0.001*state.p);
          dTb:=(tb-saturationTemperature(state.p*0.999))/(0.001*state.p);
          kappa:=(fluidConstants[1].MW / (1000 * R * state.T)+(dDgs-fluidConstants[1].MW / (1000 * R * state.T))*(tb/state.T)^9+9*(dgs-dgi)*tb^8*dTb/state.T^9)/state.d;
          //kappa := 1.0 / state.p "ideal gas compressibility";
        else
          kappa := Modelica.Constants.inf;
        end if;
    end isothermalCompressibility;
  
    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
  
      protected
        Real cp, kappa, beta;
        Real cv;
  
      algorithm
        if state.gf == 0.0 then
          kappa := isothermalCompressibility(state);
          beta := isobaricExpansionCoefficient(state);
          cp := specificHeatCapacityCp(state);
          cv := specificHeatCapacityCv(state);
          a := sqrt(max(0, cp / (cp * kappa * state.d - state.T * beta * beta)));
        elseif state.gf == 1.0 then
          cp := PhysPropCorr(fluidConstants[1].Cp0Corr, fluidConstants[1].Cp0Coef, fluidConstants[1].MW, state.T);
          cv := cp - R * 1000 / fluidConstants[1].MW;
          kappa := 1 / state.p "as ideal gas";
          a := sqrt(max(0, cp / (cv * kappa * state.d))) "velocity of sound as ideal gas" ;
        else
          a := 0;
        end if;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end velocityOfSound;
  
    function jouleThomsonCoefficient
      input ThermodynamicState state;
      output Real J;
    algorithm
      if state.phase==1 then
        J := (state.T * isobaricExpansionCoefficient(state) - 1) / (state.d * specificHeatCapacityCp(state));
      else
        J:=Modelica.Constants.inf;
      end if;
    end jouleThomsonCoefficient;
  
    function isothermalThrottlingCoefficient "OK. The ExternalMedia function is not found !!"
      input ThermodynamicState state;
      output Real I;
    algorithm
      if state.phase == 1 then
        I := -(state.T * isobaricExpansionCoefficient(state) - 1) / state.d;
      else
        I := Modelica.Constants.inf;
      end if;
    end isothermalThrottlingCoefficient;
  
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
            eta:=FreeFluids.MediaCommon.Functions.gasViscPcorLucas(fluidConstants[1], state.T,state.p,lpVisc);
            //eta := lpVisc "no pressure correction seems necessary for gas phase";
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
          lambda:=FreeFluids.MediaCommon.Functions.gasThCondTVcorChung(state.T, state.d,fluidConstants[1],lambda);     
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
        state := ThermodynamicState(p = Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), h = Modelica.Media.Common.smoothStep(x, state_a.h, state_b.h, x_small), T = Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small), d = Modelica.Media.Common.smoothStep(x, state_a.d, state_b.d, x_small), gf = Modelica.Media.Common.smoothStep(x, state_a.gf, state_b.gf, x_small), phase = 0);
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
        Inline = true,
        inverse(p = pressure_dT(d, T), T = temperature_pd(p, d)));
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
      if state.phase==1 then
        h_der := (-p_der * ((state.T * isobaricExpansionCoefficient(state) - 1) / state.d)) + T_der * cp;
      else
      h_der := Modelica.Constants.inf;
      end if;
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
        Inline = true,
        inverse(h = specificEnthalpy_pd(p, d)));
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
  T := temperature_ph_state(p=p, h=h, state=setState_phX(p=p, h=h, phase=phase));
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
        Inline = true,
        inverse(d = density_pT(p, T), T = temperature_pd(p, d)));
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
        LateInline = true,
        derivative(noDerivative = state) = specificEnthalpy_dT_der);
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
      if state.phase == 1 then
        dpdT := 1 / (state.d * isothermalCompressibility(state));
      else
        dpdT := Modelica.Constants.small;
      end if;
    end pressure_derd_T;
  
    function pressure_derT_d "pressure derivative w.r.t. T at constant density. OK"
      input ThermodynamicState state;
      output Real dpTd;
    algorithm
      if state.phase == 1 then
        dpTd := isobaricExpansionCoefficient(state) / isothermalCompressibility(state);
      else
        dpTd := -PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, state.T) / ((1 / PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) - 1 / PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, state.T)) * state.T);
      end if;
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
  
    function specificEnthalpy_derd_T "pressure derivative w.r.t. density at constant T."
      input ThermodynamicState state;
      output Real dhdT;
    protected
      SaturationProperties sat;
    algorithm
      if state.phase == 1 then
        dhdT := (1 - isobaricExpansionCoefficient(state) * state.T) / (isothermalCompressibility(state) * state.d ^ 2);
      elseif state.phase == 2 then
        sat := setSat_T(state.T);
        dhdT := PhysPropCorr(fluidConstants[1].HvCorr, fluidConstants[1].HvCoef, fluidConstants[1].MW, state.T) / ((1 / PhysPropCorr(fluidConstants[1].lDensCorr, fluidConstants[1].lDensCoef, fluidConstants[1].MW, state.T) - 1 / PhysPropCorr(fluidConstants[1].gSatDensCorr, fluidConstants[1].gSatDensCoef, fluidConstants[1].MW, state.T)) * state.d ^ 2);
      end if;
    end specificEnthalpy_derd_T;
  
    function specificEnthalpy_derT_d "pressure derivative w.r.t. density at constant T. To be finished analyticaly for two phases"
      input ThermodynamicState state;
      output Real dhTd;
    protected
      Real beta;
      ThermodynamicState statePlus;
    algorithm
      if state.phase == 1 then
        beta := isobaricExpansionCoefficient(state);
        dhTd := specificHeatCapacityCp(state) + beta * (1 - beta * state.T) / state.d;
      elseif state.phase == 2 then
        statePlus := setState_dTX(state.d, state.T + 0.02);
        dhTd := (statePlus.h - state.h) / 0.02;
      end if;
    end specificEnthalpy_derT_d;
  end TMedium;




  annotation(
    Documentation(info = "<html>
    <body>
    <p>There are several very reliable medium packages based en multiparameter equations of state; there are also very simple medium packages based in the ideal gas equation for gases, or in constant properties for liquids. But for many products we do not have multiparameter EOS, or we do not need their complexity, and we need need a good approximation of physical properties. The TMedia package is based in correlations an offers good accuracy for liquids till say 200 bars pressure and for gases till 20-30 bars.</p>
    <p>The medium is designed for liquid, liquid/vapor, or gas phases. For liquid or biphasic states, its application is limited to the higher temperature limits of the liquid heat capacity and vaporization enthalpy correlations used, but temperature should be lower than 0.85 Tc, and pressure not higher than 200 bars, because at higher values the influence of pressure on properties becomes very difficult to correct for. For gas state, it is limited to the maximum pressure which saturation temperature is below the maximum temperature limit of the liquid heat capacity correlation, but should be limited to a maximum of 20-30 bars. It can't work with supercritical states. It extends the Modelica PartialTwoPhaseMedium. It is somewhat similar to the Modelica TableBased medium, but uses specific correlations for each physical property, allows to work with gas phase, and adds a density dependent correlation for the reduced bulk modulus of the liquid, that improves a lot the calculation of liquid density at high pressure, isothermal compressibility, and isobaric expansion coefficient. Improving also the calculation of liquid heat capacity at constant volume (Cv) and the speed of sound. The medium properties are obtained using correlations that are mainly functions of T, but different pressure corrections are also used. It uses the substances data stored in the MediaCommon package.</p>
    <p>The use of pressure correction is controlled by the constant Boolean 'highPressure'. Its default value is false. If switched to true, pressure correction will be applied (in plus than to liquid density) to liquid specific enthalpy, specific entropy, heat capacity, viscosity and thermal conductivity. It is interesting to make highPressure=true if we need to work over 20 or 30 bars, but the price is a slower simulation.</p>
    <p>The values of enthalpy and entropy are calculated from a reference state. The reference state to use can be selected giving value to the constant string 'refState'. The values can be: 'ASHRAE', 'IIR', 'NBP' or 'User'. Any other value will eliminate any correction for the reference state. When using 'User', the raw values at reference_T will be used as zero for both enthalpy and entropy.</p>
    <p>The thermodynamic record contains: p,T,gas fraction, d and h. Care must be taken in limiting the use to the temperature limits of the correlations used, as only few checks are done by the media, in order not to interfere with the solver process.</p>
    <p>A constant string 'localInputChoice' has been added to the BaseProperties model in order to specify the independent variables to use in each instance of the BaseProperties model. The default value for this constant is the value given to the constant string 'inputChoice' at package level. The valid alternatives are: 'ph', 'pT', 'dT'.</p>
    <p>The global idea has been not to use the Modelica files for the storage of substances data, but to store the data in a database, from which we can recover and use them when needed. A database is provided with more than 400 substances that can be enlarged, and the FreeFluids GUI can retrieve the data from the data base, treat it as needed (for example creating EOS from saturated vapor pressure and/or densities, or creating correlations from the EOS), store the results in the database, and export the data in Modelica format when needed. Nevertheless, in order to make life easier for users, many common substances have been exported, and their packages included in the TMedia.Fluids package.</p>
    <p>As a resume: The medium is for fast calculation of liquid phase, condensation, evaporation, and gas phase below the critical point. In the liquid and saturated phases, the results are quite good. In the gas phase, the results are better than the ideal gas approach in density and enthalpy. The medium is compatible with OpenModelica 1.17 old and new frontends. The medium is also compatible, since the addition of derivative functions calculation, with the ThermoPower library and with Modelica.Fluid.</p>
    <b><p>Liquid phase properties</p></b>                
    <p>The saturated density is calculated using a dedicated correlation. This density is corrected for pressure influence. If the coefficients for the reduced bulk modulus calculation, as function of density, are available, the correction is done using a very accurate Tait like equation. In other case, a substance specific isothermal compressibility factor (with a default value of 6.667 e-10) is used. The parameters for the reduced bulk modulus correlation (with liquid density as independent variable) are normally not available, but can be calculated from a good equation of state of the multiparameter or SAFT types. This can be done easily with the FreeFluids GUI: you make the calculation with the EOS, transfer the results (density and the natural logarithm of the reduced bulk modulus) to the correlations tab, make the regression of the coefficients, and store the result in the database. It is good to calculate the reduced bulk modulus at a pressure close to 50 bars but, if necessary in order to have liquid phase at the temperature of interest, it can be done at higher pressure. Check that all density data used correspond to the liquid state.</p>
    <p>A dedicated correlation is used also for the saturated heat capacity. The use of the saturated liquid Cp correlation, instead of ideal Cp correlation plus vaporization enthalpy, makes possible the use of the medium with substances for which we do not have Cp0 data, and improves the liquid phase thermal properties calculation. This correlation can be also a problem, as many times we only find it with a temperature limit of the normal boiling point. This can be solved using a Cp correlation constructed from a good EOS, using the FreeFluids GUI. It is important not to use data too close to the Tc for the regression (use data till 10 K below the Tc). The best equation for the regression of the liquid heat capacity is the PPDS15 equation. Do not use the ChemSep equations as they are not integrated by the medium to obtain enthalpy or entropy. If highPressure has been made equal to true, a density dependent correction is applied to the Cp.</p>
    <p>The liquid enthalpy is calculated from the liquid Cp correlation at saturation. If highPressure has been made equal to true, PV correction is applied. The correction interpolates linearly from full correction below 0.45Tc to none at 0.85Tc.</p>
    <p>When going back from liquid enthalpy, or entropy, to temperature, we have two ways: One of them is to fit a correlation that makes this calculation, that can again be made with FreeFluidsGui. The second, that will be used if we left the value of lTfromHsatCorr=0, will calculate in situ the temperature using a solving function, as is done in the IdealGasMedia package.</p>
    <p>For water, viscosity is calculated, independent of the phase, using the reference equation as function of temperature and density. For other substances, the saturated liquid viscosity is calculated using a dedicated temperature dependent correlation. Pressure correction, according to Lucas, is applied if highPressure is set to true.</p> 
    <p>Saturated thermal conductivity is calculated using the corresponding correlation, if available. Otherwise the Latini method is used. If highPressure is set to true, a pressure(in fact density) correction is used.</p>
    <p>Surface tension is calculated also by correlation, or using the Sastri-Rao method, if the first is not available.</p> 
    <b><p>Gas phase properties</p></b>
    <p>Saturated gas density is calculated using a dedicated correlation. At temperature higher than saturation, for a given pressure, a temperature correction is introduced.</p>
    <p>Enthalpy and entropy are calculated from the saturated values at the given pressure (obtained from liquid Cp and vaporization enthalpy correlations), by adding the increase from the saturated temperature to the given temperature, calculated according to the ideal gas Cp, but with a pressure and temperature correction.</p> 
    <p>Viscosity and thermal conductivity are calculated by correlations or, when not available, by the Chung method. Pressure correction is later applied.</p>
    <b><p>Two phases properties</p></b>
    <p>Transport properties are not calculated for the two phases situation.</p>
    </body>
    </html>"));
end TMedia;

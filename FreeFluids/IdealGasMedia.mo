within FreeFluids;

package IdealGasMedia
  "IdealGasMedia.mo by Carlos Trujillo
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
  
  
  //*****MEDIUMS USING IDEAL GAS CALCULATIONS*****
  //==============================================
  
  //***** PACKAGE IdealGasMedium*****
  //*********************************
  partial package IdealGasMedium
    extends Modelica.Media.Interfaces.PartialPureSubstance(ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT);
  
      constant Boolean useTransportCorr=true "if true correlations will be used for transport properties
   calculation, if available";
  
    //Data storage definitions
    //------------------------
      constant FreeFluids.MediaCommon.DataRecord data;  
  
    //Auxiliary functions based in correlations
    //-----------------------------------------

    function SpecificEnthalpyCorr "Calculates specific enthalpy from a given Cp correlation at a given temperature, using an external function. No constant adjustment is used in the integration."
      input Real T;
      output Real h;
      external "C" FF_SpecificEnthalpyCorr(data.Cp0Corr, data.Cp0Coef, data.MW, T, h) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFphysprop.c\"");
    end SpecificEnthalpyCorr;
  
    //Function SpecificEnthalpyCorrInv
    function SpecificEnthalpyCorrInv "Compute temperature from property value"
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal
        extends Modelica.Media.Common.OneNonLinearEquation;

        redeclare record extends f_nonlinear_Data 
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
      T := Internal.solve(y, 50.0, 5000.0, 1.0e5, {1}, fd);
    end SpecificEnthalpyCorrInv;
    
    function SpecificEntropyCorr "Calculates specific entropy from a given Cp correlation at a given temperature, using an external function. No constant adjustment is used in the integration."
      input Real T;
      output Real s;
      external "C" FF_SpecificEntropyCorr(data.Cp0Corr, data.Cp0Coef, data.MW, T, s) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFphysprop.c\"");
    end SpecificEntropyCorr;
    
    //Function SpecificEntropyCorrInv
    function SpecificEntropyCorrInv "Compute temperature from property value"
      input Real y "Property value";
      output Temperature T "Temperature";
    protected
      package Internal
        extends Modelica.Media.Common.OneNonLinearEquation;
    
        redeclare record extends f_nonlinear_Data
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
      T := Internal.solve(y, 50.0, 5000.0, 1.0e5, {1}, fd);
    end SpecificEntropyCorrInv;
  
  //BaseProperties model
    //--------------------
  
    redeclare model extends BaseProperties(
     T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
     p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
  equation
  
    if (T<data.Tc) and (data.VpCorr>0) then
      assert(p<=FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, T),"The media can´t be used in liquid state", AssertionLevel.warning);
    end if;
      assert((state.T >= 200.0) and (state.T <= data.Cp0LimS), "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= Cp0 high limit of the medium model  \"" + mediumName + "\"");
  assert(p<=1.0e6,"Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
  
      MM = data.MW/1000 "in kg";
      R = Modelica.Constants.R/MM;
      h = SpecificEnthalpyCorr(T);
      u = h - R*state.T;
      d = p/(R*T);
  
      state.T = T;
      state.p = p;
      state.h = h;
      state.d = d;
  end BaseProperties;
  
  //Thermodynamic state definition and constructors
    //-----------------------------------------------
  
    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        AbsolutePressure p "Pressure in Pa";
        Temperature T "Kelvin temperature";
        Density d (displayUnit="kg/m3") "density in kg/m3";
        SpecificEnthalpy h "specific enthalpy";
    end ThermodynamicState;
  
    redeclare     function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi"
          extends Modelica.Icons.Function;
      
        protected
          AbsolutePressure vp "vapor pressure at given T";
        algorithm
          state.p := p;
          state.T := T;
          state.d:=p*data.MW/(1000*R*T);
          state.h:= SpecificEnthalpyCorr(T);
          
          vp := if state.T < data.Tc and data.VpCorr > 0 then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, state.T) else 1.0e10;
          assert(state.p<=vp,"The media can´t be used in liquid state", AssertionLevel.warning);
          assert((state.T >= 200.0) and (state.T <= data.Cp0LimS), "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= Cp0 high limit of the medium model  \"" + mediumName + "\"");
          assert(state.p<=1.0e6,"Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning); 
        end setState_pTX;
  
    redeclare function extends setState_dTX "Return ThermodynamicState record as function of T,d and composition X or Xi"
      extends Modelica.Icons.Function;
      protected
        AbsolutePressure vp "vapor pressure at given T";
      algorithm
        state.T := T;
        state.d := d;
        state.p:=d*1000*R*T/data.MW;
        state.h:= SpecificEnthalpyCorr(state.T);
        
        vp := if state.T < data.Tc and data.VpCorr > 0 then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, state.T) else 1.0e10;
        assert(state.p<=vp,"The media can´t be used in liquid state", AssertionLevel.warning);
        assert((state.T >= 200.0) and (state.T <= data.Cp0LimS), "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= Cp0 high limit of the medium model  \"" + mediumName + "\"");
        assert(state.p<=1.0e6,"Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning); 
        
    end setState_dTX;
    
    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,H and composition X or Xi"
      extends Modelica.Icons.Function;
      protected
        AbsolutePressure vp "vapor pressure at given T";
      algorithm
        state.p := p;
        state.h := h;
        state.T := SpecificEnthalpyCorrInv(h);
        state.d:=p*data.MW/(1000*R*state.T);
        
        vp := if state.T < data.Tc and data.VpCorr > 0 then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, state.T) else 1.0e10;
        assert(state.p<=vp,"The media can´t be used in liquid state", AssertionLevel.warning);
        assert((state.T >= 200.0) and (state.T <= data.Cp0LimS), "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= Cp0 high limit of the medium model  \"" + mediumName + "\""); 
        assert(state.p<=1.0e6,"Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning); 
    end setState_phX;
    
      redeclare function extends setState_psX "Return ThermodynamicState record as function of p,S and composition X or Xi"
      extends Modelica.Icons.Function;
        protected
          AbsolutePressure vp "vapor pressure at given T";
      algorithm
        state.p := p;
        state.T := SpecificEntropyCorrInv(s + R * Modelica.Math.log(state.p/reference_p)*1000/data.MW);
        if state.T<data.Tc then
          vp:=FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, state.T);
          assert(state.p<=vp,"The media can´t be used in liquid state"); 
        end if;
        state.d:=p*data.MW/(1000*R*state.T);
        state.h:= SpecificEnthalpyCorr(state.T);
        
        vp := if state.T < data.Tc and data.VpCorr > 0 then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.VpCorr, data.VpCoef, data.MW, state.T) else 1.0e10;
        assert(state.p<=vp,"The media can´t be used in liquid state", AssertionLevel.warning); 
        assert((state.T >= 200.0) and (state.T <= data.Cp0LimS), "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= Cp0 high limit of the medium model  \"" + mediumName + "\"");
        assert(state.p<=1.0e6,"Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning); 
    end setState_psX;
  
  //Properties calculation from thermodynamic state
    //-----------------------------------------------
  
    redeclare function extends molarMass
        "Return the molar mass of the medium"
      extends Modelica.Icons.Function;
    algorithm
      MM:=data.MW*1e-3;
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
        d :=state.p*data.MW/(1000*R*state.T);
    end density;
  
    redeclare function extends specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
  
      algorithm
        h := state.h;
    end specificEnthalpy;
  
    redeclare function extends specificInternalEnergy "Return specific internal energy"
      extends Modelica.Icons.Function;
    algorithm
      u := state.h - R*1000*state.T/data.MW;
      annotation(Inline=true,smoothOrder=2);
    end specificInternalEnergy;
    
    redeclare function extends specificEntropy "Return specific entropy"
      extends Modelica.Icons.Function;
    algorithm
      s := SpecificEntropyCorr(state.T) - R * Modelica.Math.log(state.p/reference_p)*1000/data.MW;
    end specificEntropy;
  
    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
      extends Modelica.Icons.Function;
    algorithm
      g := state.h - state.T*specificEntropy(state);
      annotation(Inline=true,smoothOrder=2);
    end specificGibbsEnergy;
  
    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
      extends Modelica.Icons.Function;
    algorithm
      f := state.h - R*state.T*1000/data.MW - state.T*specificEntropy(state);
      annotation(Inline=true,smoothOrder=2);
    end specificHelmholtzEnergy;
  
    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        cp := FreeFluids.MediaCommon.Functions.PhysPropCorr(data.Cp0Corr, data.Cp0Coef, data.MW, state.T);
    end specificHeatCapacityCp;
  
    redeclare function extends specificHeatCapacityCv "Compute specific heat capacity at constant volume from temperature and gas data"
        extends Modelica.Icons.Function;
    algorithm
      cv := FreeFluids.MediaCommon.Functions.PhysPropCorr(data.Cp0Corr, data.Cp0Coef, data.MW, state.T) - R*1000/data.MW;
      annotation(Inline=true,smoothOrder=2);
    end specificHeatCapacityCv;
  
    redeclare function extends isentropicExponent "Return isentropic exponent"
      extends Modelica.Icons.Function;
    protected
      Real cp;
    algorithm
      cp:= FreeFluids.MediaCommon.Functions.PhysPropCorr(data.Cp0Corr, data.Cp0Coef, data.MW, state.T);
      gamma := cp/(cp- R*1000/data.MW);
      annotation(Inline=true,smoothOrder=2);
    end isentropicExponent;
  
    redeclare function extends isobaricExpansionCoefficient
      "Returns overall the isobaric expansion coefficient beta"
      extends Modelica.Icons.Function;
    algorithm
      beta := 1/state.T;
      annotation(Inline=true,smoothOrder=2);
    end isobaricExpansionCoefficient;
  
    redeclare function extends isothermalCompressibility
      "Returns overall the isothermal compressibility factor"
      extends Modelica.Icons.Function;
    algorithm
      kappa := 1.0/state.p;
      annotation(Inline=true,smoothOrder=2);
    end isothermalCompressibility;
  
    redeclare function extends density_derp_T
      "Returns the partial derivative of density with respect to pressure at constant temperature"
      extends Modelica.Icons.Function;
    algorithm
      ddpT := data.MW/(state.T*R*1000);
      annotation(Inline=true,smoothOrder=2);
    end density_derp_T;
  
    redeclare function extends density_derT_p
      "Returns the partial derivative of density with respect to temperature at constant pressure"
      extends Modelica.Icons.Function;
    algorithm
      ddTp := -data.MW*state.p/(state.T*state.T*R*1000);
      annotation(Inline=true,smoothOrder=2);
    end density_derT_p;
    
    redeclare function extends isentropicEnthalpy "Return an approximation of isentropic enthalpy (gamma is considered constant)"
      extends Modelica.Icons.Function;
    protected
      IsentropicExponent gamma "Isentropic exponent";
      
    algorithm
      gamma := isentropicExponent(refState); 
      h_is := refState.h + gamma/(gamma - 1.0)*refState.p/refState.d*((p_downstream/refState.p)^((gamma - 1)/gamma) - 1.0);
    annotation(Inline=true,smoothOrder=2);
    end isentropicEnthalpy;
  
    redeclare function extends velocityOfSound "Return velocity of sound"
      extends Modelica.Icons.Function;
    protected
      Real cp;
    algorithm
      cp:= FreeFluids.MediaCommon.Functions.PhysPropCorr(data.Cp0Corr, data.Cp0Coef, data.MW, state.T);
      a := sqrt(max(0,R*state.T*cp*1000/data.MW/(cp-R*1000/data.MW)));
      annotation(Inline=true,smoothOrder=2);
    end velocityOfSound;
  
    redeclare function extends dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
    protected
      Real k= if data.MW>18.014 and data.MW<18.016 then 0.076 else 0.0;
      algorithm
        assert(useTransportCorr==false or data.gViscCorr==0 or state.T<=data.gViscLimS,"You are over the maximum temperature for the low pressure viscosity correlation", AssertionLevel.warning);
        eta:= if (data.gViscCorr>0 and useTransportCorr==true) then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.gViscCorr, data.gViscCoef, data.MW, state.T) else     FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(data, state.T);
    end dynamicViscosity;
  
    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;
  
      algorithm
        assert(useTransportCorr==false or data.gThCondCorr==0 or state.T<=data.gThCondLimS,"You are over the maximum temperature for the low pressure th.conductivity correlation", AssertionLevel.warning);
        lambda := if (data.gThCondCorr>0 and useTransportCorr==true) then FreeFluids.MediaCommon.Functions.PhysPropCorr(data.gThCondCorr, data.gThCondCoef, data.MW, state.T) else FreeFluids.MediaCommon.Functions.gasThCondLowPressureChung(data, specificHeatCapacityCp(state), dynamicViscosity(state),state.T);
    end thermalConductivity;
  
  redeclare function extends setSmoothState
      "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
        algorithm
          state := ThermodynamicState(p=Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small),
                                      d=Modelica.Media.Common.smoothStep(x, state_a.d, state_b.d, x_small),
                                      h=Modelica.Media.Common.smoothStep(x, state_a.h, state_b.h, x_small),
                                      T=Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
          annotation(Inline=true,smoothOrder=2);
        end setSmoothState;
  end IdealGasMedium;

  //***** PACKAGE Template*****

  package IdealGasMediumTemplate
    extends IdealGasMedium(final mediumName = "no name", final singleState = false, data=FreeFluids.MediaCommon.MediaData.MediaDataTemplate);
  end IdealGasMediumTemplate;

  //***** PACKAGE Acetone*****

  package Acetone
    extends IdealGasMedium(final mediumName = "Acetone", final singleState = false, data=FreeFluids.MediaCommon.MediaData.Acetone); 
  end Acetone;

  //***** PACKAGE Air*****

  package Air
    extends IdealGasMedium(final mediumName = "Air", final singleState = false, data=FreeFluids.MediaCommon.MediaData.Air); 
  end Air;

  package Ammonia
    extends IdealGasMedium(final mediumName = "Ammonia", final singleState = false, data = FreeFluids.MediaCommon.MediaData.Ammonia);
  end Ammonia;

  //***** PACKAGE n-Butanol*****
  
  package n_Butanol
    extends IdealGasMedium(final mediumName = "n-Butanol", final singleState = false, data = FreeFluids.MediaCommon.MediaData.n_Butanol);
  end n_Butanol;

  //***** PACKAGE Carbon dioxide*****

  package CO2
    extends IdealGasMedium(final mediumName = "Carbon dioxide", final singleState = false, data=FreeFluids.MediaCommon.MediaData.CO2); 
  end CO2;

  //***** PACKAGE Nitrogen*****

  package N2
    extends IdealGasMedium(final mediumName = "Nitrogen", final singleState = false, data=FreeFluids.MediaCommon.MediaData.N2); 
  end N2;

  //***** PACKAGE Nitrogen*****

  package O2
    extends IdealGasMedium(final mediumName = "Oxygen", final singleState = false, data=FreeFluids.MediaCommon.MediaData.O2); 
  end O2;
  
//***** PACKAGE Water*****

  package Water
    extends IdealGasMedium(final mediumName = "Water", final singleState = false, data=FreeFluids.MediaCommon.MediaData.Water); 
  end Water;

  package Tests
    partial model FluidTesting
      replaceable package Medium = Modelica.Media.Interfaces.PartialPureSubstance;
      parameter Medium.AbsolutePressure p = 1.0e5;
      parameter Medium.Temperature initialT=473.15;
      parameter Medium.Temperature finalT=273.15;
      parameter Medium.AbsolutePressure p2=10.0e5;
      Medium.Temperature T(start = initialT);
      Medium.ThermodynamicState StateP "state from p,T";
      Medium.ThermodynamicState StateH "state from p,h";
      Medium.ThermodynamicState StateS "state from p,s";
      Medium.ThermodynamicState StateD "state from d,T";
      Medium.SpecificEnthalpy H "StateP specific enthalpy";
      Medium.SpecificEntropy S "StateP specific entropy";
      Real Cp;
      Real Cv;
      Real G "isentropic exponent";
      SI.Velocity SS "StateP speed of sound";
      Medium.SpecificEnthalpy H2 "isentropic enthalpy at pressure p2";
      Medium.Density D(displayUnit="kg/m3") "StateP density";
      Medium.DynamicViscosity Mu "StateP dynamic viscosity";
      Real Th "StateP thermal conductivity";
      Real DD_T;
      Real DD_P;
      Medium.BaseProperties BaseProp;
    algorithm
      StateP := Medium.setState_pTX(p, T, fill(0, 0));
      H:=Medium.specificEnthalpy(StateP);
      S:=Medium.specificEntropy(StateP);
      Cp:=Medium.specificHeatCapacityCp(StateP);
      Cv:=Medium.specificHeatCapacityCv(StateP);
      G:=Medium.isentropicExponent(StateP);
      SS:=Medium.velocityOfSound(StateP);
      H2:=Medium.isentropicEnthalpy(p2,StateP);
      D := Medium.density(StateP);
      Mu := Medium.dynamicViscosity(StateP);
      Th := Medium.thermalConductivity(StateP);
      DD_T:=Medium.density_derT_p(StateP);
      DD_P:=Medium.density_derp_T(StateP);
      StateH := Medium.setState_phX(p, H, fill(0, 0));
      StateS:=Medium.setState_psX(p,S,fill(0,0));
      StateD:=Medium.setState_dTX(D,T,fill(0,0));
    equation
      BaseProp.p=p;
      BaseProp.T=T;
      der(T) = finalT-initialT;
      annotation(
        Documentation(info = "<html>
    <body>
    <p>It crates a thermodynamic state from a fixed pressure(defined by the parameter p), and a temperature T that ramps between the parameters initialT and finalT</p>
    <p>Some thermodynamic properties are calculated from this state.</p>
    <p>Later the state is reproduced, as StateH, from p and the calculated enthalpy. As StateS, from p and the calculated entropy. And as StateD from T and the calculated density. All states should have the same values for the state variables, if the reconstruction of the state has been acurate</p>
    <p>You can extend the model and run it, asigning the desired Medium and adjusting the parameters</p> 
    </body>
    </html>"));
    end FluidTesting;

    model Test1A
      extends FluidTesting(redeclare replaceable package Medium = FreeFluids.IdealGasMedia.N2(useTransportCorr=true), p=1.0e5, p2=1.0e6, initialT=1000.0, finalT=273.15);  
    end Test1A;
    model Test1B
      extends Test1A(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.N2);  
    end Test1B;
    
    model Test1C
      extends Test1A(redeclare package Medium = FreeFluids.IdealGasMedia.Air(useTransportCorr=true));  
    end Test1C;
    
    model Test1D
      extends Test1A(redeclare package Medium=Modelica.Media.Air.DryAirNasa);  
    end Test1D;

    model IdealGasH2OA "IdealGas H20 medium model"
      extends Modelica.Icons.Example;
      replaceable package Medium = FreeFluids.IdealGasMedia.Water "Medium model";
      Medium.ThermodynamicState state "Thermodynamic state record";
      Medium.ThermodynamicState state2;
      Medium.SpecificHeatCapacity cp = Medium.specificHeatCapacityCp(state);
      Medium.SpecificHeatCapacity cv = Medium.specificHeatCapacityCv(state);
      Medium.IsentropicExponent k = Medium.isentropicExponent(state);
      Medium.SpecificEntropy s = Medium.specificEntropy(state);
      //  Medium.SpecificEntropy s2=Medium.specificEntropy(state2);
      Medium.VelocityOfSound a = Medium.velocityOfSound(state);
      Real beta = Medium.isobaricExpansionCoefficient(state);
      Real gamma = Medium.isothermalCompressibility(state);
      Medium.SpecificEnthalpy h_is = Medium.isentropicEnthalpy(2.0e5, state);
      Medium.ThermodynamicState smoothState;
      Real m_flow_ext;
      Real der_p;
      Real der_T;
      Real P1, T1, P2, T2;
    equation
      P1 = 1.0e5;
      T1 = 280.0 + 1000.0 * time;
      P2 = 2.0e5;
      T2 = 700.0;
      state = Medium.setState_pTX(P1, T1, fill(0, 0));
      state2 = Medium.setState_pTX(P2, T2, fill(0, 0));
//  s2 = s;
// Smooth state
      m_flow_ext = time - 0.5;
      smoothState = Medium.setSmoothState(m_flow_ext, state, state2, 0.1);
      der_p = der(smoothState.p);
      der_T = der(smoothState.T);
      annotation(
        Documentation(info = "<html>
<p>An example for using ideal gas properties and how to compute
isentropic enthalpy changes.
The function that is implemented is approximate, but usually
 very good: the second medium record medium2
 is given to compare the approximation.
</p>
</html>"),
        experiment(StopTime = 1));
    end IdealGasH2OA;

    model IdealGasH2OB "IdealGas H20 medium model"
      extends IdealGasH2OA(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2O);  end IdealGasH2OB;
  
  end Tests;
    annotation(
    Documentation(info = "<html>
<body>
<p>The medium is designed for single substances in gas phase, at a pressure low enough and a temperature high enough, to allow for the ideal gas equation to be used. It extends the Modelica PartialPureSubstance medium. The definition is similar to that of the Modelica.Media.IdealGases.Common.SingleGasNasa, but uses several ecuations for the Cp0 correlation, as the Nasa Glenn coefficients are not available for many organic compounds. Look at the MediaData package information for details on how to use the database to create new substances. The use of a single equation for Cp0 limits somewhat the temperature range. The DIPPR107 equation in J/(kg·K) is recommended</p>
<p>It checks, if the vapour pressure correlation is supplied, that the media is in the gas state, warning if not.</p>
<p>Density is calculated using the ideal gas equation of state. Enthalpy and entropy are calculated from the ideal gas constant pressure heat capacity Cp0, using specific temperature correlations. No constant is added to the raw calculation.</p>
<p>For transport properties, correlations between temperature and the (low pressure) property are used if the constant useTransportCorr==true, and they are available. If not, the Chung method is used for their estimation.</p>
<p>The thermodynamic record contains: p,T,d and h.</p>
<p>As a resume: The medium is for fast calculation of gas phase at low pressure and not too low temperature.</p>
</body>
</html>"));
end IdealGasMedia;

within FreeFluids;

package IdealGasMixture "IdealGasMixture.mo by Carlos Trujillo
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
  //***** PACKAGE IdealGasMixture*****
  //*********************************

  partial package IdealGasMix
    extends Modelica.Media.Interfaces.PartialMixtureMedium(ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pTX, 
    substanceNames = data[:].name, 
    reducedX = false, 
    fixedX = true, 
    singleState = false, 
    reference_X = fill(1 / nX, nX),
      Temperature(min=200, max=6000, start=500, nominal=500), 
      SpecificEnthalpy(start=h_default, nominal=h_default),
      Density(start=10, nominal=10),
      AbsolutePressure(start=5e5, nominal=5e5));
    constant Boolean useTransportCorr=true;
    //Data storage definitions
    //------------------------
    constant FreeFluids.MediaCommon.DataRecord[:] data;
    constant Integer dippr107= min({if data[i].Cp0Corr==200 then 1 else 0 for i in 1:nS});
  
    //Auxiliary functions based in correlations
    //-----------------------------------------
    //Function SpecificEnthalpyMix
  
    function SpecificEnthalpyMix "Compute the ideal specific enthalpy of a mix from T"
      input Temperature T;
      input Real X[nS];
      output SpecificEnthalpy h;
    algorithm
  h:= X*FreeFluids.MediaCommon.Functions.SpecificEnthalpyCorr2(data, T);
    end SpecificEnthalpyMix;
  
      //Function SpecificEnthalpyCorrInv
  
      function SpecificEnthalpyMixInv "Compute temperature from property value"
        //input Correlation corrData;
        input Real y "Property value";
        input MassFraction[nS] X "composition as mass fractions";
        output Temperature T "Temperature that gives the property";
      protected
        package Internal "Solve y(data,T) for T with given y (use only indirectly via temperature_phX)"
          extends Modelica.Media.Common.OneNonLinearEquation;
  
          redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
          end f_nonlinear_Data;
  
          redeclare function extends f_nonlinear
            algorithm
              y := SpecificEnthalpyMix(x,X);
          end f_nonlinear;
  
          // Dummy definition has to be added for current Dymola??
  
          redeclare function extends solve
          end solve;
        end Internal;
  
        Internal.f_nonlinear_Data fd;
      algorithm
        T := Internal.solve(y, 50.0, min(data.Cp0LimS), 1.0e5, X, fd);
      end SpecificEnthalpyMixInv;
  
  
    //Function SpecificEntropyMix
  
    function SpecificEntropyMix "Compute the ideal specific entropy of a mix from T"
      input Temperature T;
      input Real X[nS];
      output SpecificEntropy s;
    algorithm
  s :=  X * FreeFluids.MediaCommon.Functions.SpecificEntropyCorr2(data, T);
    end SpecificEntropyMix;
  
      //Function SpecificEntropyCorrInv
  
      function SpecificEntropyMixInv "Compute temperature from property value"
        //input Correlation corrData;
        input Real y "Property value";
        input MassFraction[nS] X "composition as mass fractions";
        output Temperature T "Temperature that gives the property";
      protected
        package Internal "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
          extends Modelica.Media.Common.OneNonLinearEquation;
  
          redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
          end f_nonlinear_Data;
  
          redeclare function extends f_nonlinear
            algorithm
              y := SpecificEntropyMix(x,X);
          end f_nonlinear;
  
          // Dummy definition has to be added for current Dymola??
  
          redeclare function extends solve
          end solve;
        end Internal;
  
        Internal.f_nonlinear_Data fd;
      algorithm
        T := Internal.solve(y, 50.0, min(data.Cp0LimS), 1.0e5, X, fd);
      end SpecificEntropyMixInv;
  
  
    //BaseProperties model
    //--------------------
  
    redeclare model extends BaseProperties(T(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), final standardOrderComponents = true)
        Real nMols "number of moles in a kg";
  
      algorithm
        nMols := 1000 * sum(X ./ data.MW);
        h := SpecificEnthalpyMix(T, X);
      equation
        assert(T >= 200.0 and T <= 1500.0, "Temperature T = " + String(T) + " K is not in the allowed range 200 K <= T <= 1500 K of the medium model  \"" + mediumName + "\"",AssertionLevel.warning);
        assert(p <= 1.0e6, "Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
  //nMols={data[i].MW for i in 1:nS};
  //MM = 1.0/(1000*nMols) "in kg";
        MM = 1 / nMols;
        R = Modelica.Constants.R * nMols;
        u = h - R * state.T;
        d = p / (R * T);
        state.X = X;
        state.MM =MM;
        state.T = T;
        state.p = p;
        state.h = h;
        state.d = d;
    end BaseProperties;
  
    //Thermodynamic state definition and constructors
    //-----------------------------------------------
  
    redeclare record extends ThermodynamicState
        extends Modelica.Icons.Record;
        Density d(displayUnit = "kg/m3") "density in kg/m3";
        SpecificEnthalpy h "specific enthalpy";
        MolarMass MM;
    end ThermodynamicState;
  
    redeclare function extends setState_pTX "Return ThermodynamicState record as function of p,T and composition X or Xi"
        extends Modelica.Icons.Function;
      protected
        Real nMols "number of moles in a kg";
        SpecificHeatCapacity Rm "mix R in J/(kg·K)";
        AbsolutePressure vp "vapor pressure at given T";
  
      algorithm
        state.X := X;
        nMols := 1000 * sum(state.X ./ data.MW);
        Rm := Modelica.Constants.R * nMols;
        state.MM:=1/nMols;
        state.p := p;
        state.T := T;
        state.d := p / (Rm * T);
        state.h := SpecificEnthalpyMix(state.T, state.X);
        
        assert(state.T >= 200.0 and state.T <= 1500.0, "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= 1500 K of the medium model  \"" + mediumName + "\"");
        assert(state.p <= 1.0e6, "Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
    end setState_pTX;
  
    redeclare function extends setState_dTX "Return ThermodynamicState record as function of T,d and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        Real nMols "number of moles in a kg";
        SpecificHeatCapacity Rm "mix R in J/(kg·K)";
        AbsolutePressure vp "vapor pressure at given T";
  
      algorithm
        state.X := X;
        nMols := 1000 * sum(state.X ./ data.MW);
        Rm := Modelica.Constants.R * nMols;
        state.MM:=1/nMols;
        state.T := T;
        state.d := d;
        state.p := d * Rm * T;
        state.h := SpecificEnthalpyMix(state.T, state.X);
        assert(state.T >= 200.0 and state.T <= 1500.0, "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= 1500 K of the medium model  \"" + mediumName + "\"");
        assert(state.p <= 1.0e6, "Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
    end setState_dTX;
  
    redeclare function extends setState_phX "Return ThermodynamicState record as function of p,H and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        Real nMols "number of moles in a kg";
        SpecificHeatCapacity Rm "mix R in J/(kg·K)";
        AbsolutePressure vp "vapor pressure at given T";
  
      algorithm
        state.X := X;
        nMols := 1000 * sum(state.X ./ data.MW);
        Rm := Modelica.Constants.R * nMols;
        state.MM:=1/nMols;
        state.p := p;
        state.h := h;
        state.T := SpecificEnthalpyMixInv(state.h,state.X);
        state.d := state.p / (Rm * state.T);
        assert(state.T >= 200.0 and state.T <= 1500.0, "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= 1500 K of the medium model  \"" + mediumName + "\"");
        assert(state.p <= 1.0e6, "Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
    end setState_phX;
  
    redeclare function extends setState_psX "Return ThermodynamicState record as function of p,S and composition X or Xi"
        extends Modelica.Icons.Function;
  
      protected
        Real nMols "number of moles in a kg";
        SpecificHeatCapacity Rm "mix R in J/(kg·K)";
        AbsolutePressure vp "vapor pressure at given T";
  
      algorithm
        state.X := X;
        nMols := 1000 * sum(state.X ./ data.MW);
        Rm := Modelica.Constants.R * nMols;
        state.MM:=1/nMols;
        state.p := p;
        state.T := SpecificEntropyMixInv(s + R * Modelica.Math.log(state.p / reference_p) * nMols, state.X);
        state.d := state.p / (Rm * state.T);
        state.h := SpecificEnthalpyMix(state.T,state.X);
        assert(state.T >= 200.0 and state.T <= 1500.0, "Temperature T = " + String(state.T) + " K is not in the allowed range 200 K <= T <= 1500 K of the medium model  \"" + mediumName + "\"",AssertionLevel.warning);
        assert(state.p <= 1.0e6, "Ideal gas EOS is not adequate for pressure higher than 10 bars", AssertionLevel.warning);
    end setState_psX;
  
    //Properties calculation from thermodynamic state
    //-----------------------------------------------
  
    redeclare function extends molarMass "Return the molar mass of the medium"
        extends Modelica.Icons.Function;
  
      algorithm
        MM := state.MM;
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
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificInternalEnergy;
  
    redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
  
      algorithm
        s := SpecificEntropyMix(state.T, state.X) - R * Modelica.Math.log(state.p / reference_p) / state.MM;
    end specificEntropy;
  
    redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
  
      algorithm
        g := state.h - state.T * specificEntropy(state);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificGibbsEnergy;
  
    redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
  
      algorithm
        f := state.h - state.p / state.d - state.T * specificEntropy(state);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificHelmholtzEnergy;
  
    redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        cp:= state.X*FreeFluids.MediaCommon.Functions.Cp0Corr(data, state.T);
  
    end specificHeatCapacityCp;
  
    redeclare function extends specificHeatCapacityCv "Compute specific heat capacity at constant volume from temperature and gas data"
        extends Modelica.Icons.Function;
  
      algorithm
        cv := specificHeatCapacityCp(state) - R / state.MM;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end specificHeatCapacityCv;
  
    redeclare function extends isentropicExponent "Return isentropic exponent"
        extends Modelica.Icons.Function;
  
      protected
        Real cp;
  
      algorithm
        cp := specificHeatCapacityCp(state);
        gamma := cp / (cp - R / state.MM);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end isentropicExponent;
  
    redeclare function extends isobaricExpansionCoefficient "Returns overall the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
  
      algorithm
        beta := 1 / state.T;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end isobaricExpansionCoefficient;
  
    redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
  
      algorithm
        kappa := 1.0 / state.p;
      annotation(
        Inline = true,
        smoothOrder = 2);
    end isothermalCompressibility;
  
    redeclare function extends density_derp_T "Returns the partial derivative of density with respect to pressure at constant temperature"
        extends Modelica.Icons.Function;
  
      algorithm
        ddpT := state.MM / (state.T * R);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end density_derp_T;
  
    redeclare function extends density_derT_p "Returns the partial derivative of density with respect to temperature at constant pressure"
        extends Modelica.Icons.Function;
  
      algorithm
        ddTp := -state.MM * state.p / (state.T * state.T * R);
      annotation(
        Inline = true,
        smoothOrder = 2);
    end density_derT_p;
  
    redeclare function extends isentropicEnthalpy "Return an approximation of isentropic enthalpy (gamma is considered constant)"
        extends Modelica.Icons.Function;
  
      protected
        IsentropicExponent gamma "Isentropic exponent";
  
      algorithm
        gamma := isentropicExponent(refState);
        h_is := refState.h + gamma / (gamma - 1.0) * refState.p / refState.d * ((p_downstream / refState.p) ^ ((gamma - 1) / gamma) - 1.0);
      annotation(
        smoothOrder = 2);
    end isentropicEnthalpy;
  
    redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
  
      protected
        Real cp;
      algorithm
        cp := specificHeatCapacityCp(state);
        a := sqrt(max(0, R * state.T * cp / (state.MM * cp - R)));
      annotation(
        Inline = true,
        smoothOrder = 2);
    end velocityOfSound;
  
    redeclare function extends dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
    protected
        DynamicViscosity etaX[nS];
        Fraction MF[nS];
      algorithm
      for i in 1:nS loop
        etaX[i]:= if (data[i].gViscCorr>0 and useTransportCorr==true) then FreeFluids.MediaCommon.Functions.gasViscCorr(data[i], state.T) else FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(data[i],state.T);
        assert(useTransportCorr==false or data[i].gViscCorr==0 or state.T<=data[i].gViscLimS,"You are over the maximum temperature for the low pressure viscosity correlation", AssertionLevel.warning);        
      end for;
        MF:=massToMoleFractions(state.X, data[:].MW);
        eta:=FreeFluids.MediaCommon.Functions.gasMixViscosityWilke(state.T,MF,etaX,data[:].MW);
    end dynamicViscosity;
  
    redeclare function extends thermalConductivity "Return thermal conductivity"
        extends Modelica.Icons.Function;
    protected
        DynamicViscosity lambdaX[nS];
        DynamicViscosity etaX[nS];
        SpecificHeatCapacity Cp0X[nS];
        Fraction MF[nS];
    algorithm
      for i in 1:nS loop
        etaX[i]:= if (data[i].gViscCorr>0 and useTransportCorr==true) then FreeFluids.MediaCommon.Functions.gasViscCorr(data[i], state.T) else FreeFluids.MediaCommon.Functions.gasViscLowPressureChung(data[i],state.T);
        
        Cp0X[i]:=FreeFluids.MediaCommon.Functions.Cp0Corr(data[i], state.T);
        
        lambdaX[i]:= if (data[i].gThCondCorr>0 and useTransportCorr==true) then FreeFluids.MediaCommon.Functions.gasThCondCorr(data[i], state.T) else FreeFluids.MediaCommon.Functions.gasThCondLowPressureChung(data[i], Cp0X[i], etaX[i], state.T);
        assert(useTransportCorr==false or data[i].gThCondCorr==0 or state.T<=data[i].gThCondLimS,"You are over the maximum temperature for the low pressure thermal conductivity correlation", AssertionLevel.warning);        
      end for;
        
        MF:=massToMoleFractions(state.X, data[:].MW);
        lambda:=FreeFluids.MediaCommon.Functions.gasMixThCondMason(state.T,MF,lambdaX,data.MW,data.Tc,data.criticalPressure);
    end thermalConductivity;
    
    redeclare function extends setSmoothState
      "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
      extends Modelica.Icons.Function;
    algorithm
      state := ThermodynamicState(
      p=Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), 
      T=Media.Common.smoothStep(x, state_a.T, state_b.T, x_small),
      d=Media.Common.smoothStep(x, state_a.d, state_b.d, x_small),
      h=Media.Common.smoothStep(x, state_a.h, state_b.h, x_small),  
      X=Media.Common.smoothStep(x, state_a.X, state_b.X, x_small),
      X=Media.Common.smoothStep(x, state_a.MM, state_b.MM, x_small));
      annotation(Inline=true,smoothOrder=2);
    end setSmoothState;
  end IdealGasMix;

  //***** PACKAGE Template*****

  package IdealGasMediumTemplate
    extends IdealGasMix(final mediumName = "no name", final singleState = false, data = {FreeFluids.MediaData.MediaDataTemplate});
  end IdealGasMediumTemplate;

  //***** PACKAGE Air*****

  package Air
    extends IdealGasMix(final mediumName = "Air", final singleState = false, data = {FreeFluids.MediaCommon.MediaDataMZ.N2, FreeFluids.MediaCommon.MediaDataMZ.O2}, fixedX = true, reference_X = {0.768, 0.232});
  end Air;
  
  package FlueGas
    extends IdealGasMix(final mediumName = "Combustion gas", final singleState = false, data = {FreeFluids.MediaCommon.MediaDataMZ.Water, FreeFluids.MediaCommon.MediaDataMZ.N2, FreeFluids.MediaCommon.MediaDataMZ.O2,FreeFluids.MediaCommon.MediaDataAL.CO2}, fixedX = true, reference_X = {0.12, 0.71, 0.03, 0.14});
  end FlueGas;
  
  package Tests
    partial model FluidTesting
      replaceable package Medium = Modelica.Media.Interfaces.PartialMixtureMedium;
      parameter Medium.AbsolutePressure p = 1.0e5;
      parameter Medium.Temperature initialT = 473.15;
      parameter Medium.Temperature finalT = 273.15;
      parameter Medium.AbsolutePressure p2 = 1.0e6;
      Medium.Temperature T(start = initialT);
      Medium.ThermodynamicState StateP "state from p,T";
      Medium.ThermodynamicState StateH "state from p,h";
      Medium.ThermodynamicState StateS "state from p,s";
      Medium.ThermodynamicState StateD "state from d,T";
      Medium.SpecificEnthalpy H "StateP specific enthalpy";
      Medium.SpecificEntropy S "StateP specific entropy";
      Medium.SpecificInternalEnergy U;
      Medium.MolarMass MM;
      Real Cp;
      Real Cv;
      Medium.IsentropicExponent G;
      SI.Velocity SS "StateP speed of sound";
      Medium.SpecificEnthalpy H2 "isentropic enthalpy at pressure p2";
      Medium.Density D(displayUnit = "kg/m3") "StateP density";
      Real DD_T, DD_P;
      Medium.DynamicViscosity Mu "StateP dynamic viscosity";
      Real Th "StateP thermal conductivity";
      Medium.BaseProperties BaseProp;
    algorithm
      StateP := Medium.setState_pTX(p, T);
      H := Medium.specificEnthalpy(StateP);
      S := Medium.specificEntropy(StateP);
      U := Medium.specificInternalEnergy(StateP);
      MM := Medium.molarMass(StateP);
      Cp := Medium.specificHeatCapacityCp(StateP);
      Cv := Medium.specificHeatCapacityCv(StateP);
      G := Medium.isentropicExponent(StateP);
      SS := Medium.velocityOfSound(StateP);
      H2 := Medium.isentropicEnthalpy(p2, StateP);
      D := Medium.density(StateP);
      DD_T:=Medium.density_derT_p(StateP);
      DD_P:=Medium.density_derp_T(StateP);
      Mu := Medium.dynamicViscosity(StateP);
      Th := Medium.thermalConductivity(StateP);
      StateH := Medium.setState_phX(p, H);
      StateS := Medium.setState_psX(p, S);
      StateD := Medium.setState_dTX(D, T);
    equation
      BaseProp.p = p;
      BaseProp.T = T;
      der(T) = finalT - initialT;
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
      extends FluidTesting(redeclare replaceable package Medium = FreeFluids.IdealGasMixture.Air(useTransportCorr=true), p = 1.0e5, initialT = 1000.0, finalT = 700.0);
    end Test1A;

    model Test1B
      extends Test1A(redeclare package Medium = Modelica.Media.IdealGases.MixtureGases.CombustionAir(fixedX = true));
    end Test1B;

    model Test1C
      extends Test1A(redeclare package Medium = FreeFluids.IdealGasMixture.FlueGas(fixedX = true));
    end Test1C;

    model Test1D
      extends Test1A(redeclare package Medium = Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents(reference_X={0.71,0.0,0.0,0.03,0.12,0.14},fixedX = true));
    end Test1D;
    
  end Tests;

  annotation(
    Documentation(info = "<html>
<body>
<p>The medium is designed for a mixture of substances in the gas phase, at a pressure low enough to allow for the ideal gas equation to be used. It extends the Modelica PartialMixtureMedium. The definition is similar to that of the Modelica.Media.IdealGases.Common.MixtureGasNasa, but uses several equations for the Cp0 correlation, as the NASA Glenn coefficients are not available for many organic compounds. The use of a single equation for Cp0 limits somewhat the temperature range.</p>
<p> It shares the substances data with the other medium models. Look at the MediaData package information for details on how to use the database to create new substances.</p>
<p>It uses also correlations for gas viscosity and thermal conductivity.</p>
<p>Density is calculated using the ideal gas equation of state. Enthalpy and entropy are calculated from the ideal gas constant pressure heat capacity Cp0, using specific temperature correlations.</p>
<p>For transport properties, correlations between temperature and the (low pressure) property are used if available. If not, the Chung method is used for their estimation. Later mixing rules are applied.</p>
<p>The thermodynamic record contains: X,p,T,d and h.</p>
<p>As a resume: The medium is for fast calculation of gas phase at low pressure and not too low temperature.</p>
</body>
</html>"));
end IdealGasMixture;

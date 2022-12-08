within FreeFluids.TMedia;

package Tests "In Tests Axx are comparing TMedia with multiparameter EOS and ExternalMedia, using a temperature ramp.
In Tests Bxx a comparison of the gas phase only also with ideal gas."
  partial model FluidTestingA
    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
    parameter Medium.AbsolutePressure p = 1.0e5;
    parameter Medium.Temperature initialT = 273.19;
    parameter Medium.Temperature finalT = 372.15;
    parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
    Medium.Temperature T(start = initialT) "We will ramp the temperature";
    Medium.MolarMass MM;
    //creation of state from p and T and its properties
    Medium.ThermodynamicState StateP "original state from p and T";
    Real Tb;
    Real H;
    Real D;
    Real S;
    Real U "internal energy";
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
    //reconstruction of the state from other variables
    Medium.ThermodynamicState StateD "state reproduced from d,T";
    Medium.ThermodynamicState StateH "state reproduced from p,h";
    Medium.ThermodynamicState StateS "state reproduced from p,s";
    //creation of a lower density state, to force two phases, and its properties
    Medium.ThermodynamicState StateDlow "state at fraction of original density";
    Real DlowGasFract "gas fraction of StateDlow";
    Real DlowD, DlowH;
    Real DlowS;
    Real DlowU;
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
    //Real BubDerD_p;
    Real BubH "bubble properties";
    Real BubS;
    //Real BubDerH_p;
    Real DewD;
    Real DewH "dew properties";
    Real DewS;
    //Real DewDerD_p;
    //Real DewDerH_p;
    Real Hv "vaporization enthalpy";
    Real DewMu;
    Real Sigma;
    Medium.Temperature Tsat "temperature recovered from saturation pressure";
    //BaseProperties model
    Medium.BaseProperties BaseProp;
  algorithm
//Construction of StateP and calculation of properties
    for i in 1:10 loop
      StateP := Medium.setState_pTX(p, T);
      MM := Medium.molarMass(StateP);
      Tb := Medium.saturationTemperature(p);
      H := Medium.specificEnthalpy(StateP);
      D := Medium.density(StateP);
      S := Medium.specificEntropy(StateP);
      U := Medium.specificInternalEnergy(StateP);
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
      StateS := Medium.setState_psX(p, S);
      StateDlow := Medium.setState_dTX(D * fract, T, fill(0, 0));
      DlowGasFract := Medium.vapourQuality(StateDlow);
      DlowD := Medium.density(StateDlow);
      DlowH := Medium.specificEnthalpy(StateDlow);
      DlowU := Medium.specificInternalEnergy(StateDlow);
      DlowS := Medium.specificEntropy(StateDlow);
      DlowCv := Medium.specificHeatCapacityCv(StateDlow);
      DlowDerD_p_h := Medium.density_derp_h(StateDlow);
      DlowDerD_h_p := Medium.density_derh_p(StateDlow);
      sat := Medium.setSat_T(T);
      Vp := Medium.saturationPressure_sat(sat);
      DerTb_p := Medium.saturationTemperature_derp_sat(sat);
      StateBub := Medium.setBubbleState(sat);
      StateDew := Medium.setDewState(sat);
      BubD := Medium.bubbleDensity(sat);
      BubH := Medium.bubbleEnthalpy(sat);
      BubS := Medium.bubbleEntropy(sat);
      DewD := Medium.dewDensity(sat);
      DewH := Medium.dewEnthalpy(sat);
      Hv := DewH - BubH;
      DewS := Medium.dewEntropy(sat);
      DewMu := Medium.dynamicViscosity(StateDew);
      Tsat := Medium.saturationTemperature(sat.psat);
      Sigma := Medium.surfaceTension(sat);
    end for;
//DerD_p_T := Medium.density_derp_T(StateP);
//DerD_T_p := Medium.density_derT_p(StateP);
//BubDerD_p := Medium.dBubbleDensity_dPressure(sat);
//BubDerH_p := Medium.dBubbleEnthalpy_dPressure(sat);
//DewDerD_p := Medium.dDewDensity_dPressure(sat);
//DewDerH_p := Medium.dDewEnthalpy_dPressure(sat);
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
    extends FluidTestingA(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), p = 50e5, initialT = 0.1 + 273.15, fract = 0.1, finalT = 380 + 273.15);
  end TestA1A;

  model TestA1B
    extends TestA1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_pT);
  end TestA1B;

  model TestA2A
    extends FluidTestingA(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.R134A(highPressure = true, refState = "IIR", reference_T = 100.0, inputChoice = "pT"), p = 20.0e5, initialT = (-50.0) + 273.15, finalT = 110.0 + 273.15, fract = 0.5);
  end TestA2A;

  partial model FluidTestingB
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
    parameter Medium.Temperature T = 298.15;
    parameter Medium.AbsolutePressure initialP = 3e5;
    parameter Medium.AbsolutePressure finalP = 30e5;
    Medium.AbsolutePressure p(start = initialP) "We will ramp the temperature";
    Medium.ThermodynamicState StateP "original state from p and T";
    //Medium.ThermodynamicState StateH "state reproduced from p,h";
    //Medium.ThermodynamicState StateS "state reproduced from p,s";
    //Medium.ThermodynamicState StateD "state reproduced from d,T";
    //Properties of StateP
    Real D, H, Mu, Th;
    Real S;
    Real Beta "isobaric expansion coefficient";
    Real Kappa "isothermal compressibility";
    Real Cp;
    Real Cv;
    Real Gamma;
    Real SS "speed of sound";
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
//StateD := Medium.setState_dTX(D, T, fill(0, 0));
//StateH := Medium.setState_phX(p, H, fill(0, 0));
//StateS := Medium.setState_psX(p, S, fill(0, 0));
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
    extends FluidTestingB(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.CO2(refState = "IIR", reference_T = 273.15, highPressure = true, inputChoice = "ph"), T = 50 + 273.15, initialP = 30.0e5, finalP = 5e5);
  end TestB1A;

  model TestB1C
    extends TestB1A(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CO2);
  end TestB1C;

  model TestModel
    extends Modelica.Media.Examples.Utilities.PartialTestModel(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "pT"));
  equation

    annotation(
      Documentation(info = "<html><head></head><body>Run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst");
  end TestModel;

  model TestModel2
    extends Modelica.Media.Examples.Utilities.PartialTestModel2(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "pT"));
  equation

  annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");end TestModel2;

  model ThreeTanks
    extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "ph"));
  equation

    annotation(
      experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  Documentation(info = "<html><head></head><body>The ThreeTanks example from the Modelica standard library run with FreeFluids.TMedia.Fluids.Water.</body></html>"));
  end ThreeTanks;

  model HeatingSystem
    extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "ph"));
  equation
  
    annotation(
      experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-06, Interval = 12),
      Documentation(info = "<html><head></head><body>The Modelica example model using a FreeFluids.TMedia medium. Run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst");
  end HeatingSystem;

  model TestModelS
    extends Modelica.Media.Examples.Utilities.PartialTestModel(redeclare package Medium = Modelica.Media.Water.StandardWater);
  equation

  annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");end TestModelS;

  model TestModel2S
    extends Modelica.Media.Examples.Utilities.PartialTestModel2(redeclare package Medium = Modelica.Media.Water.StandardWater);
  equation

  end TestModel2S;

  model ThreeTanksS
    extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = Modelica.Media.Water.StandardWater);
  equation

    annotation(
      experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  Documentation);
  end ThreeTanksS;

  model HeatingSystemS
    extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = Modelica.Media.Water.StandardWater);
  equation

    annotation(
      experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-6, Interval = 12));
  end HeatingSystemS;

  model TestB2A
    extends FluidTestingB(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Propane(refState = "IIR", reference_T = 273.15, highPressure = true, inputChoice = "ph"), T = 100 + 273.15, initialP = 30.0e5, finalP = 0.1e5);
  end TestB2A;

  model TestB2C
    extends TestB2A(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C3H8);
  end TestB2C;
end Tests;

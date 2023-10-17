within FreeFluids.ExternalPure;

package Examples

  model TestingPressure
    package Medium = FreeFluids.ExternalPure.Fluids.Pentane_n(thermoModel = 2);
    parameter Real initD = 10, finalD = 690.8;
    parameter Real T1 = -50 + 273.15, T2 = 90 + 273.15;
    Real d(start = initD);
    Real V;
    Real p1;
    Real p2;
  algorithm
    for i in 1:1 loop
      p1 := Medium.pressureEOS_dT(d, T1);
      p2 := Medium.pressureEOS_dT(d, T2);
    end for;
  equation
    V=Medium.fluidK.molarMass/d;
    der(d) = finalD - initD;
  end TestingPressure;

  model TestingDensity
    package Medium = FreeFluids.ExternalPure.Fluids.Air(thermoModel = 3, refState = 3);
    parameter Medium.Temperature T(displayUnit = "K") = 640.15;
     parameter Medium.AbsolutePressure initialP = 5e5;
    parameter Medium.AbsolutePressure finalP = 5e5;
    Medium.AbsolutePressure P(start = initialP);
    Medium.Density ld(displayUnit = "kg/m3"), gd(displayUnit = "kg/m3");
  algorithm
    for i in 1:1 loop
      (ld, gd) := Medium.densities_pT(P, T, "b");
    end for;
  equation
    der(P) = finalP - initialP;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.005));
  end TestingDensity;

  model TestingVp
    package Medium1 = FreeFluids.ExternalPure.Fluids.Pentane_n(thermoModel = 1, refState = 2);
    package Medium2 = FreeFluids.ExternalPure.Fluids.Butene_1(thermoModel = 3, refState = 2);
    //package Medium3 = FreeFluids.ExternalPure.Fluids.Pentane_n(thermoModel = 3, refState = 2);
    package Medium3=ExternalMedia.Media.CoolPropMedium(mediumName = "1-Butene", substanceNames = {"1-Butene"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, inputChoice = ExternalMedia.Common.InputChoice.hs, SpecificEnthalpy(start = 2e5));
    parameter SI.Temperature initialT(displayUnit = "K") = -186+273.15;
    parameter SI.Temperature finalT(displayUnit = "K") = 150+273.15;
    SI.Temperature T(displayUnit = "K", start = initialT) "We will ramp the temperature";
    //SI.AbsolutePressure P1;
    SI.AbsolutePressure P2;
    SI.AbsolutePressure P3;
  algorithm
    //P1 := Medium1.saturationPressure(T);
    P2 := Medium2.saturationPressure(T);
    P3 := Medium3.saturationPressure(T);
  equation
    der(T) = finalT - initialT;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end TestingVp;

  model TestingTb
    package Medium1 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel=1, refState = 2);
    package Medium2 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 2, refState = 2);
    package Medium3 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 3, refState = 2);
    parameter SI.AbsolutePressure initialP(displayUnit = "bar") = 500;
    parameter SI.AbsolutePressure finalP(displayUnit = "bar") = 2e5;
    SI.AbsolutePressure P(displayUnit = "bar", start = initialP) "We will ramp the pressure";
    SI.Temperature T1(displayUnit = "K");
    SI.Temperature T2(displayUnit = "K");
    SI.Temperature T3(displayUnit = "K");
  algorithm
    T1 := Medium1.saturationTemperature(P);
  T2 := Medium2.saturationTemperature(P);
  T3 := Medium3.saturationTemperature(P);
  equation
    der(P) = finalP - initialP;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.005));
  end TestingTb;

  model TestingHfromTP
    package Medium1 = FreeFluids.TMedia.Fluids.MethylEthylKetone(refState = "IIR", highPressure = true, inputChoice = "pT");
    //package Medium2 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 2, refState = 2);
    //package Medium3 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 3, refState = 2);
    parameter SI.AbsolutePressure P1 = 1e5;
    parameter SI.AbsolutePressure P2 = 56.9e5;
    parameter SI.AbsolutePressure P3 = 200e5;
    parameter SI.Temperature initialT(displayUnit = "K") = 200.0;
    parameter SI.Temperature finalT(displayUnit = "K") = 800;
    SI.Temperature T(displayUnit = "K", start = initialT) "We will ramp the temperature";
    SI.SpecificEnthalpy H1;
    SI.SpecificEnthalpy H2;
    SI.SpecificEnthalpy H3;
  algorithm
    H1 := Medium1.specificEnthalpy(Medium1.setState_pTX(P1, T));
    H2 := Medium1.specificEnthalpy(Medium1.setState_pTX(P2, T));
    H3 := Medium1.specificEnthalpy(Medium1.setState_pTX(P3, T));
  equation
    der(T) = finalT - initialT;
  end TestingHfromTP;

  model TestingSfromTP
    package Medium1 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 2, refState = 2);
    //package Medium2 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 2, refState = 2);
    //package Medium3 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 3, refState = 2);
    parameter SI.AbsolutePressure P1 = 0.1e5;
    parameter SI.AbsolutePressure P2 = 56.9e5;
    parameter SI.AbsolutePressure P3 = 200e5;
    parameter SI.Temperature initialT(displayUnit = "K") = 200.0;
    parameter SI.Temperature finalT(displayUnit = "K") = 800;
    SI.Temperature T(displayUnit = "K", start = initialT) "We will ramp the temperature";
    SI.SpecificEntropy S1;
    SI.SpecificEntropy S2;
    SI.SpecificEntropy S3;
  algorithm
    S1 := Medium1.specificEntropy(Medium1.setState_pTX(P1, T));
    S2 := Medium1.specificEntropy(Medium1.setState_pTX(P2, T));
    S3 := Medium1.specificEntropy(Medium1.setState_pTX(P3, T));
  equation
    der(T) = finalT - initialT;
  end TestingSfromTP;

  model TestingHfromTV
    package Medium1 = FreeFluids.TMedia.Fluids.MethylEthylKetone(refState = "IIR", highPressure = true, inputChoice = "pT");
    //package Medium2 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 2, refState = 2);
    //package Medium3 = FreeFluids.ExternalPure.Fluids.Acetone(thermoModel = 3, refState = 2);
    parameter SI.Density D1 = 1;
    parameter SI.Density D2 = 350;
    parameter SI.Density D3 = 700;
    parameter SI.Temperature initialT(displayUnit = "K") = 200.0;
    parameter SI.Temperature finalT(displayUnit = "K") = 800;
    SI.Temperature T(displayUnit = "K", start = initialT) "We will ramp the temperature";
    SI.SpecificEnthalpy H1;
    SI.SpecificEnthalpy H2;
    SI.SpecificEnthalpy H3;
  algorithm
    H1 := Medium1.specificEnthalpy(Medium1.setState_dTX(D1, T));
    H2 := Medium1.specificEnthalpy(Medium1.setState_dTX(D2, T));
    H3 := Medium1.specificEnthalpy(Medium1.setState_dTX(D3, T));
  equation
    der(T) = finalT - initialT;
  end TestingHfromTV;

  model TestingBubbleDew
    package Medium = FreeFluids.ExternalPure.Fluids.MethylMethacrylate(thermoModel=2, ancillaries=2, inputChoice = "pT");
    parameter Medium.Temperature initialT = -30 + 273.15;
    parameter Medium.Temperature finalT = 320 + 273.15;
    Medium.Temperature T(start = initialT);
    //Medium.AbsolutePressure p;
    Real bd, dd;
    Real  bh, bs, dh, ds;
  algorithm
//p := Medium.saturationPressure(T);
    bd := Medium.bubbleDensity(Medium.setSat_T(T));
    bh := Medium.bubbleEnthalpy(Medium.setSat_T(T));
    bs := Medium.bubbleEntropy(Medium.setSat_T(T));
    dd := Medium.dewDensity(Medium.setSat_T(T));
    dh := Medium.dewEnthalpy(Medium.setSat_T(T));
    ds := Medium.dewEntropy(Medium.setSat_T(T));
  equation
    der(T) = finalT - initialT;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.1));
  end TestingBubbleDew;
  
  model TestingViscosity
    package Medium = FreeFluids.ExternalPure.Fluids.CarbonMonoxide(thermoModel=3);
    parameter Medium.AbsolutePressure P(displayUnit = "bar") = 5e5;
    parameter Medium.Temperature initialT = 100;
    parameter Medium.Temperature finalT = 300;
    Medium.Temperature T(start = initialT);
    Medium.ThermodynamicState StateP=Medium.setState_pT(P,T);
    Medium.DynamicViscosity Mu;
algorithm
    Mu:=Medium.dynamicViscosity(StateP);
equation
  der(T) = finalT - initialT;  
  end TestingViscosity;
end Examples;

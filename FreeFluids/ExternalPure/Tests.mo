within FreeFluids.ExternalPure;

package Tests
  partial model FluidTest
    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
    parameter Medium.AbsolutePressure Ps = 1.0e5;
    parameter Medium.Temperature initialT = 273.19;
    parameter Medium.Temperature finalT = 372.15;
    Medium.Temperature Ts(start = initialT) "We will ramp the temperature";
    //StateP construction and properties from it, in 1 phase
    Medium.ThermodynamicState StateP "original state from p and T";
    Medium.MolarMass MM;
    Medium.SpecificEnthalpy H;
    Medium.Density D(displayUnit = "kg/m3");
    Medium.SpecificEntropy S;
    Real U "internal energy";
    Real Cp;
    Real Cv;
    Real Gamma;
    Real SS "speed of sound";
    Real Beta "isobaric expansion coefficient";
    Real Kappa "isothermal compressibility";
    //Real J "JouleThomson coefficient";
    //Real I "Isothermal throttling coefficient";
    //Real Mu;
    //Real Th;
    //Medium.SurfaceTension Sigma;
    Real Density_derp_T;
    Real Density_derT_p;
    Real Density_derp_h;
    Real Density_derh_p;
    //Reconstruction of the state from other variables
    Medium.ThermodynamicState StateD "state reproduced from d,T";
    Medium.ThermodynamicState StateH "state reproduced from p,h";
    //Medium.ThermodynamicState StateS "state reproduced from p,s";
    //Saturation properties and states
    Medium.SaturationProperties sat "saturation point at given T";
    Medium.AbsolutePressure Vp;
    Medium.Temperature Tb;
    Medium.Temperature Tsat "temperature recovered from saturation pressure";
    Medium.ThermodynamicState StateBub "bubble state at sat";
    Medium.ThermodynamicState StateDew "dew state at sat";
    Real BetaBub "isobaric expansion coefficient";
    Real KappaBub "isothermal compressibility";
    Real BetaDew "isobaric expansion coefficient";
    Real KappaDew "isothermal compressibility";
    //properties from the sat record
    Real BubD;
    Real BubH "bubble properties";
    Real BubS;
    Real BubDerD_p;
    Real BubDerH_p;
    Real DewD;
    Real DewH "dew properties";
    Real DewS;
    Real DewDerD_p;
    Real DewDerH_p;
    //Real DewMu;
    Real Hv "vaporization enthalpy";
    //construction of a two phases state from p,h
    Medium.ThermodynamicState State2H "state half way between Bub and Dew";
    Real T2;
    Real D2;
    //density in the two phases state
    Real S2;
    //entropy in the two phases state
    Real Cv2;
    Real Cp2;
    Real Density_derp_h2;
    Real Density_derh_p2;
    //reconstruction of the state from other variables
    Medium.ThermodynamicState State2D;
    //BaseProperties model construction from different variables, also with two phases
    Medium.BaseProperties BasePropT;
    Medium.BaseProperties BasePropH;
    Medium.BaseProperties BasePropD;
    Medium.BaseProperties BasePropH2;
  algorithm
    for i in 1:1 loop
  //StateP construction and properties from it, in 1 phase
      StateP := Medium.setState_pTX(Ps, Ts);
      MM := Medium.molarMass(StateP);
      H := Medium.specificEnthalpy(StateP);
      D := Medium.density(StateP);
      S := Medium.specificEntropy(StateP);
      U := Medium.specificInternalEnergy(StateP);
      Cp := Medium.specificHeatCapacityCp(StateP);
      Cv := Medium.specificHeatCapacityCv(StateP);
      Gamma := Cp/Cv;
      SS := Medium.velocityOfSound(StateP);
      Beta := Medium.isobaricExpansionCoefficient(StateP);
      Kappa := Medium.isothermalCompressibility(StateP);
  //J:=Medium.jouleThomsonCoefficient(StateP);
  //I:=Medium.isothermalThrottlingCoefficient(StateP);
      //Mu := Medium.dynamicViscosity(StateP);
      //Th := Medium.thermalConductivity(StateP);
  //Sigma := Medium.surfaceTension(sat);
      Density_derp_T := Medium.density_derp_T(StateP);
      Density_derT_p := Medium.density_derT_p(StateP);
      Density_derp_h := Medium.density_derp_h(StateP);
      Density_derh_p := Medium.density_derh_p(StateP);
  //Reconstruction of the state from other variables
      StateD := Medium.setState_dTX(D, Ts);
      StateH := Medium.setState_phX(Ps, H);
  //StateS := Medium.setState_psX(Ps, S);
  //Saturation properties and states
      sat := Medium.setSat_T(Ts);
      Vp := Medium.saturationPressure_sat(sat);
      Tb := Medium.saturationTemperature(Ps);
      Tsat := Medium.saturationTemperature(sat.psat);
      StateBub := Medium.setBubbleState(sat);
      StateDew := Medium.setDewState(sat);
      BetaBub := Medium.isobaricExpansionCoefficient(StateBub);
      KappaBub := Medium.isothermalCompressibility(StateBub);
      BetaDew := Medium.isobaricExpansionCoefficient(StateDew);
      KappaDew := Medium.isothermalCompressibility(StateDew);
  //properties from the sat record
      BubD := Medium.bubbleDensity(sat);
      BubH := Medium.bubbleEnthalpy(sat);
      BubS := Medium.bubbleEntropy(sat);
      BubDerD_p := Medium.dBubbleDensity_dPressure(sat);
      BubDerH_p := Medium.dBubbleEnthalpy_dPressure(sat);
      DewD := Medium.dewDensity(sat);
      DewH := Medium.dewEnthalpy(sat);
      DewS := Medium.dewEntropy(sat);
      DewDerD_p := Medium.dDewDensity_dPressure(sat);
      DewDerH_p := Medium.dDewEnthalpy_dPressure(sat);
  //DewMu := Medium.dynamicViscosity(StateDew);
      Hv := DewH - BubH;
  //construction of a two phases state from p,h
      State2H := Medium.setState_phX(Vp, 0.3*BubH + 0.7*DewH, fill(0, 0));
      T2 := Medium.temperature(State2H);
      D2 := Medium.density(State2H);
      S2 := Medium.specificEntropy(State2H);
      Cv2 := Medium.specificHeatCapacityCv(State2H);
      Cp2 := Medium.specificHeatCapacityCp(State2H);
      Density_derp_h2 := Medium.density_derp_h(State2H);
      Density_derh_p2 := Medium.density_derh_p(State2H);
  //reconstruction of the state from other variables
      State2D := Medium.setState_dTX(D2, T2);
    end for;
  equation
  //Construction of BaseProperties
    BasePropT.p = Ps;
    BasePropT.T = Ts;
    BasePropH.p = Ps;
    BasePropH.h = H;
    BasePropD.d = D;
    BasePropD.T = Ts;
    BasePropH2.p = Ps;
    BasePropH2.h = 0.3*BubH + 0.7*DewH;
    der(Ts) = finalT - initialT;
    annotation(
      Documentation(info = "<html><head></head><body>
         <p>This test maintains a constant pressure, and ramps temperature between the selected values. At each temperature, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p><p>Also an intermediate state between bubble and dew is created and some of its properties calculate.</p><p>You can enable also the calculation of viscosity and thermal conductivity, but the CoolProp fluids will fail many times due to lack of implementation.</p>
         
         </body></html>"));
  end FluidTest;

  model TestA
    extends FluidTest(Ps = 1e5, initialT = 273.15+50, finalT = 273.15+50);
  end TestA;

  model TestACubic
    extends TestA(redeclare replaceable package Medium = FreeFluids.ExternalPure.Fluids.Methanol(thermoModel = 1, refState = 3, reference_T = 0, reference_p = 101325, inputChoice = "pT"), BasePropT(localInputChoice = "pT"), BasePropH(localInputChoice = "ph"), BasePropD(localInputChoice = "dT"), BasePropH2(localInputChoice = "ph"));
  end TestACubic;

  model TestAPCSAFT
    extends TestACubic(Medium(thermoModel = 2));
  end TestAPCSAFT;

  model TestASW
    extends TestACubic(Medium(thermoModel = 3, ancillaries = 3));
  end TestASW;

  model TestATMedia
    extends TestACubic(redeclare package Medium = FreeFluids.TMedia.Fluids.AceticAcid(refState = "IIR", highPressure = true, inputChoice = "ph"));
  end TestATMedia;

  model TestACoolProp
    extends TestA(redeclare package Medium = ExternalMedia.Media.CoolPropMedium(mediumName = "Pentane", substanceNames = {"Propylene"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, inputChoice = ExternalMedia.Common.InputChoice.hs, SpecificEnthalpy(start = 2e5)), BasePropT(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.pT), BasePropH(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.ph), BasePropD(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.dT), BasePropH2(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.ph));
    annotation(
      Documentation(info = "<html><head></head><body>Needs the system library ExternalMedia be loaded.</body></html>"));
  end TestACoolProp;

  model TestB
    extends FluidTest(Ps = 10.0e5, initialT = (20) + 273.15, finalT = 100 + 273.15);
  end TestB;

  model TestBCubic
    extends TestB(redeclare replaceable package Medium = FreeFluids.ExternalPure.Fluids.Octane_n(thermoModel = 1, refState = 3, reference_T = 273.15, reference_p = 101325.0, inputChoice = "pT"), BasePropT(localInputChoice = "pT"), BasePropH(localInputChoice = "ph"), BasePropD(localInputChoice = "dT"), BasePropH2(localInputChoice = "ph"));
  end TestBCubic;

  model TestBPCSAFT
    extends TestBCubic(Medium(thermoModel = 2));
  end TestBPCSAFT;

  model TestBSW
    extends TestBCubic(Medium(thermoModel = 3));
  end TestBSW;

  model TestBTMedia
    extends TestBCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.Octane_n(refState = "NBP", highPressure = true, inputChoice = "pT"));
    end TestBTMedia;

  model TestBCoolProp
  extends TestB(redeclare package Medium = ExternalMedia.Media.CoolPropMedium(mediumName = "Octane", substanceNames = {"Octane"}, ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, SpecificEnthalpy(start = 2e5)), BasePropT(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.pT), BasePropH(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.ph), BasePropD(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.dT), BasePropH2(basePropertiesInputChoice = ExternalMedia.Common.InputChoice.ph));
  //{"Octane|enable_TTSE=1"} enable_BICUBIC
    annotation(
      Documentation(info = "<html><head></head><body>Needs the system library ExternalMedia be loaded.</body></html>"));
  end TestBCoolProp;

  model TestC
    extends FluidTest(Ps = 5.0e5, initialT = -19+273.15, finalT = 400+273.15);
  end TestC;

  model TestCCubic
    extends TestC(redeclare replaceable package Medium = FreeFluids.ExternalPure.Fluids.Styrene(thermoModel = 1, refState = 3, reference_T = 273.15, reference_p = 101325.0, inputChoice = "ph"), BasePropT(localInputChoice = "pT"), BasePropH(localInputChoice = "ph"), BasePropD(localInputChoice = "dT"), BasePropH2(localInputChoice = "ph"));
  end TestCCubic;

  model TestCPCSAFT
    extends TestCCubic(Medium(thermoModel = 2));
  end TestCPCSAFT;

  model TestCTMedia
    extends TestCCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.Styrene(refState = "NBP", highPressure = true, inputChoice = "ph"));
  end TestCTMedia;

  model TestD
    extends FluidTest(Ps = 5.0e5, initialT = -30+273.15, finalT = 400+273.15);
  end TestD;

  model TestDCubic
    extends TestD(redeclare replaceable package Medium = FreeFluids.ExternalPure.Fluids.MethylMethacrylate(thermoModel = 1, refState = 2, reference_T = 273.15, reference_p = 101325.0, inputChoice = "pT"), BasePropT(localInputChoice = "pT"), BasePropH(localInputChoice = "ph"), BasePropD(localInputChoice = "dT"), BasePropH2(localInputChoice = "ph"));
  end TestDCubic;

  model TestDPCSAFT
    extends TestDCubic(Medium(thermoModel = 2, ancillaries = 2));
  end TestDPCSAFT;

  model TestDTMedia
    extends TestDCubic(redeclare package Medium = FreeFluids.TMedia.Fluids.MethylMethacrylate(refState = "IIR", highPressure = true, inputChoice = "pT"));
  end TestDTMedia;
end Tests;

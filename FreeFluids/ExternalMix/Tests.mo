within FreeFluids.ExternalMix;

package Tests
   partial model TestModelBase
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNames = "A,B", eosType = "PR", mixRule = "LCVM", activityModel = "UNIFACdort");
   parameter Real initialT = 373.15;
   parameter Real finalT = 373.15;
   parameter Real initialP = 1e5;
   parameter Real finalP = 1e5;
   Real Z[Medium.subsNum] "mole fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[Medium.subsNum];
   Medium.MolarMass MW "mix molecular weight in kg/mol)";
   Real Zmass[Medium.subsNum] "mass fraction of substances";
   Medium.ThermodynamicState stateBubbleP "bubble state from T and X";
   Medium.ThermodynamicState stateBubbleT;
   Medium.ThermodynamicState stateDewP;
   Medium.ThermodynamicState stateDewT;
   Medium.ThermodynamicState stateP "state at given p,T,Z";
   Medium.ThermodynamicState stateH "state at given p,h,Z. It shoud be a match or stateP";
   Real x[Medium.subsNum] "molar fractions of the liquid phase of stateP";
   Real y[Medium.subsNum] "molar fractions of the gas phase of stateP";
   Real gamma[Medium.subsNum];
   Medium.SpecificEnergy gE(displayUnit="J/mol");
   Medium.SpecificEnergy Hv;
   Medium.Density D(displayUnit = "kg/m3") "of stateP";
   Medium.SpecificEnthalpy H "of stateP";
   Medium.SpecificHeatCapacity Cp "of stateP";
   Real SS "speed of sound of stateP";
   Real Beta "isobaric expansion coefficient of stateP";
   Real Kappa "isothermal compressibility of stateP";
   Medium.DynamicViscosity Eta;
   
  algorithm
    MM := Medium.molarMasses();
    MW:=Z*MM;
    Zmass := Medium.moleToMassFractions(Z, MM);
    stateBubbleP := Medium.setBubbleState_TX(T, Zmass);
    stateBubbleT := Medium.setBubbleState_pX(p, Zmass);
    stateDewP := Medium.setDewState_TX(T, Zmass);
    stateDewT := Medium.setDewState_pX(p, Zmass);
    stateP := Medium.setState_pTX(p, T, Zmass);
    H:=Medium.specificEnthalpy(stateP);
    stateH:=Medium.setState_phX(p,H,Zmass);
    x:=Medium.massToMoleFractions(stateP.x,MM);
    y:=Medium.massToMoleFractions(stateP.y,MM);
    (gamma, gE) := Medium.activityCoefficients_TX(T, Zmass);
    Hv := Medium.specificVaporizationHeat(stateBubbleP);
    D:=Medium.density(stateP);
   
    Cp := Medium.specificHeatCapacityCp(stateP);
    SS := Medium.velocityOfSound(stateP);
    Beta:=Medium.isobaricExpansionCoefficient(stateP);
    Kappa:=Medium.isothermalCompressibility(stateP);
    Eta:=Medium.dynamicViscosity(stateP);
   
  equation
    der(T) = finalT - initialT;
    der(p) = finalP - initialP;
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
  Documentation(info = "<html><head></head><body>The model allows the testing of bubble and dew pressure and temperature, and the generation of thermodynamic states from p,T and p,h. You can ramp both pressure and temperature. You will need to fill the subsNames string, choose the thermodynamic models and give value for the molar fractions of the components.</body></html>"));
  end TestModelBase;
  
  partial model TestModel2
    extends TestModelBase(Z(start = {initialZ1, 1 - initialZ1}));
   parameter Real initialZ1 = 0.0001;
   parameter Real finalZ1 = 0.9999;
  equation
    der(Z[1]) = finalZ1 - initialZ1;
    Z[2] = 1.0 - Z[1];
  annotation(
      Documentation(info = "<html><head></head><body>For two components you can ramp the molar fraction of the first component, and it is not necessary to give the molar fraction of the second component.</body></html>"));
  end TestModel2;
    
  model EthanolWaterCubic
    extends TestModel2(Medium(subsNames = "Ethanol,WaterSRK", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 370, finalT = 370, initialP = 1.01325e5, finalP = 1.01325e5, initialZ1 = 0.001, finalZ1 = 0.999);
  
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolWaterCubic;

  model EthanolWaterPCSAFT
  extends EthanolWaterCubic(Medium(eosType = "PCSAFT"));
    annotation(
      Documentation(info = "<html><head></head><body>You can check (using WaterSRK) that the 4C schema for water predicts the azeotrope, but not the 2B schema</body></html>"));
  end EthanolWaterPCSAFT;

  model EthanolWaterHighPCubic
    extends TestModel2(Medium(mediumName = "Test", subsNames = "Ethanol,WaterRef", eosType = "PR", mixRule = "LCVM", activityModel = "NRTL"), initialT = 523, finalT = 523, initialP = 72.2e5, finalP = 72.2e5, initialZ1 = 0.0001, finalZ1 = 0.9999);
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolWaterHighPCubic;

  model EthanolWaterHighPPCSAFT
    extends EthanolWaterHighPCubic(Medium(eosType="PCSAFT"));
  end EthanolWaterHighPPCSAFT;

  model EthaneHeptaneCubic  extends TestModel2(Medium(subsNames = "Ethane,Heptane_n", eosType = "PR", mixRule = "VdW", activityModel = "UNIFACpsrk"), initialT = 403.15, finalT = 403.15, initialP = 29.2e5, finalP = 29.2e5, initialZ1 = 0.0, finalZ1 = 1.0);
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation(info = "<html><head></head><body>Here you can see the limitations of the equilibrium constants inizialization using the Wilson equation.</body></html>"));
  end EthaneHeptaneCubic;

  model EthaneHeptanePCSAFT
    extends EthaneHeptaneCubic(Medium(eosType="PCSAFT"));
  end EthaneHeptanePCSAFT;

  model EGWaterCubic
  extends TestModel2(Medium(subsNames = "EG,WaterRef", eosType = "PR", mixRule = "LCVM", activityModel = "NRTL"), initialT = 440, finalT = 440, initialP = 1.01325e5, finalP = 1.01325e5, initialZ1 = 0.001, finalZ1 = 0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
    
  end EGWaterCubic;

  model EGWaterPCSAFT
    extends EGWaterCubic(Medium(eosType="PCSAFT"));
  end EGWaterPCSAFT;

  model AcetoneChloroformCubic
    extends TestModel2(Medium(subsNames="Acetone,Chloroform",eosType="PR",mixRule="MHV2", activityModel = "UNIFACdort"),initialT=328.15,finalT=328.15,initialP=1.01325e5,finalP=1.01325e5);
  end AcetoneChloroformCubic;

  model AcetoneChloroformPCSAFT
    extends AcetoneChloroformCubic(Medium(eosType = "PCSAFT", mixRule = "IndAssoc"));
    algorithm
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end AcetoneChloroformPCSAFT;

  model EthaneButaneCubic
    extends TestModel2(Medium(subsNames = "Ethane,Butane_n", eosType = "PR", mixRule = "VdW", activityModel = "UNIFACpsrk"), initialT = 250, finalT = 250, initialP = 2e5, finalP = 2e5, initialZ1 = 0.0001, finalZ1 = 0.9999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthaneButaneCubic;

  model EthaneButanePCSAFT
    extends EthaneButaneCubic(Medium(eosType="PCSAFT"));
  end EthaneButanePCSAFT;
  
  model BenzeneHeptane_nCubic
  extends TestModel2(Medium(subsNames = "Benzene,Heptane_n", eosType = "PR", mixRule = "VdW", activityModel = "UNIFACpsrk"), initialT = 370, finalT = 370, initialP = 1.35e5, finalP = 1.35e5, initialZ1=0.0, finalZ1=1.0);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
      __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
  end BenzeneHeptane_nCubic;
  
  model AcetoneEthanolWaterCubic
  extends TestModelBase(Medium(subsNames = "Acetone,Ethanol,WaterRef", eosType = "PR", mixRule = "VdWnoInt", activityModel = "UNIFACpsrk"), initialT = 300, finalT = 350, initialP = 1.9e5, finalP = 1.5e5, Z={0.5,0.25,0.25});
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
      __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
  end AcetoneEthanolWaterCubic;
  
  model MethanolPentane_nCubic
  extends TestModel2(Medium(subsNames = "Methanol,Pentane_n", eosType = "PR", mixRule = "LCVM", activityModel = "NRTL"), initialT = 397.7, finalT = 397.7, initialP = 14e5, finalP = 14e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end MethanolPentane_nCubic;

  model MethanolPentane_nPCSAFT
  extends MethanolPentane_nCubic(Medium(eosType = "PCSAFT", mixRule = "None", activityModel = "None"));
  end MethanolPentane_nPCSAFT;

  model MethanolCyclohexaneCubic
  extends TestModel2(Medium(subsNames = "Methanol,cyclohexane", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 323.15, finalT = 323.15, initialP = 1e5, finalP = 1e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end MethanolCyclohexaneCubic;

  model MethanolCyclohexanePCSAFT
  extends MethanolCyclohexaneCubic(Medium(eosType = "PCSAFT", mixRule = "None", activityModel = "None"));
  annotation(
      Documentation(info = "<html><head></head><body>If you use Methanol2, you will see the the improvement that polar PCSAFT was doing in Methanol.</body></html>"));
  end MethanolCyclohexanePCSAFT;

  model TestLiquidStateCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNames = "EG,WaterRef", eosType = "PR", mixRule = "MHV2", activityModel = "UNIFACdort");
   parameter Real initialT = 273.15+10;   //93.33;
    //403
   parameter Real finalT = 273.15+10;//93.33;
   parameter Real initialP = 1e5;
   //29.2e5;
   parameter Real finalP = 1e5;
   //29.2e5;
   parameter Real initialZ1mass = 0.001;
   parameter Real finalZ1mass = 0.999;
   Real Z[2] "molar fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[2];
   Real Zmass[2](start = {initialZ1mass, 1 - initialZ1mass}) "mass fraction of substances";
   Medium.ThermodynamicState stateL;
   Medium.Density D(displayUnit="kg/m3");
   Real SS;
   Medium.SpecificHeatCapacity Cp;
   Medium.DynamicViscosity Eta;
  algorithm
   MM := Medium.molarMasses();
   Z := Medium.massToMoleFractions(Zmass, MM);
   stateL := Medium.setLiquidState_pTX(p, T, Zmass);
   D:=Medium.density(stateL);
   SS := Medium.velocityOfSound(stateL);
   Cp := Medium.specificHeatCapacityCp(stateL);
   Eta := Medium.dynamicViscosity(stateL);
  equation
   der(T) = finalT - initialT;
   der(p) = finalP - initialP;
   der(Zmass[1]) = finalZ1mass - initialZ1mass;
   Zmass[2] = 1.0 - Zmass[1];
  end TestLiquidStateCubic;
  
  model TestLiquiStatePCSAFT
    extends TestLiquidStateCubic(Medium(eosType="PCSAFT"));
  end TestLiquiStatePCSAFT;
end Tests;

within FreeFluids.ExternalMix;

package Tests
  partial model TestModelBase
  package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNames = "A,B", eosType = "PR", mixRule = "LCVM", activityModel = "UNIFACdort");
  parameter Real initialT = 373.15;
  parameter Real finalT = 373.15;
  parameter Real initialP = 1e5;
  parameter Real finalP = 1e5;
  Real Z[Medium.nX] "mole fraction of substances";
  Medium.Temperature T(start = initialT);
  Medium.AbsolutePressure p(start = initialP);
  Medium.MolarMass MM[Medium.nX];
  Medium.MolarMass MW "mix molecular weight in kg/mol)";
  Real Zmass[Medium.subsNum] "mass fraction of substances";
  Medium.ThermodynamicState stateBubbleP "bubble state from T and X";
  Medium.ThermodynamicState stateBubbleT;
  Medium.ThermodynamicState stateDewP;
  Medium.ThermodynamicState stateDewT;
  Medium.ThermodynamicState stateP "state at given p,T,Z";
  Medium.ThermodynamicState stateH "state at given p,h,Z. It shoud be a match of stateP";
  Medium.ThermodynamicState stateS "state at given p,s,Z. It shoud be a match of stateP";
  Medium.ThermodynamicState stateTheta "state at given p,gas fraction,Z. It shoud be a match of stateP";
  Real x[Medium.subsNum] "molar fractions of the liquid phase of stateP";
  Real y[Medium.subsNum] "molar fractions of the gas phase of stateP";
  Real gamma[Medium.subsNum];
  Medium.SpecificEnergy gE(displayUnit="J/mol");
  Medium.SpecificEnergy Hv;
  Medium.Density D(displayUnit = "kg/m3") "of stateP";
  Medium.SpecificEnthalpy H "of stateP";
  Medium.SpecificEntropy S "of stateP";
  Medium.SpecificHeatCapacity Cp "of stateP";
  Real SS "speed of sound of stateP";
  Real Beta "isobaric expansion coefficient of stateP";
  Real Kappa "isothermal compressibility of stateP";
  Medium.DynamicViscosity Eta;
  Medium.ThermalConductivity Lambda;
  Medium.SurfaceTension Sigma;
  
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
   S:=Medium.specificEntropy(stateP);
   stateH:=Medium.setState_phX(p,H,Zmass);
   stateS:=Medium.setState_psX(p,S,Zmass);
   stateTheta:=Medium.setState_pThetaX(p,stateP.gf,Zmass); 
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
   Lambda:=Medium.thermalConductivity(stateP);
   Sigma:=Medium.surfaceTension(stateP);
  
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
    extends TestModel2(Medium(subsNames = "Ethanol,WaterRef", eosType = "PR", mixRule = "VdW", activityModel = "UNIFACstd"), initialT = 273.15, finalT = 273.15, initialP = 1.0e5, finalP =1.0e5, initialZ1 = 0.001, finalZ1 = 0.999);
  
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation(info = "<html><head></head><body><br></body></html>"));
  end EthanolWaterCubic;

  model EthanolWaterPCSAFT
  extends EthanolWaterCubic(Medium(eosType = "PCSAFT",mixRule = "BL"));
    annotation(
      Documentation(info = "<html><head></head><body>You can check (using WaterSRK) that the 4C schema for water predicts the azeotrope, but not the 2B schema</body></html>"));
  end EthanolWaterPCSAFT;

  model EthanolWaterHighPCubic
    extends TestModel2(Medium(mediumName = "Test", subsNames = "Ethanol,WaterRef", eosType = "PR", mixRule = "LCVM", activityModel = "NRTL"), initialT = 523, finalT = 523, initialP = 50e5, finalP = 50e5, initialZ1 = 0.0001, finalZ1 = 0.9999);
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolWaterHighPCubic;

  model EthanolWaterHighPPCSAFT
    extends EthanolWaterHighPCubic(Medium(eosType="PCSAFT", mixRule = "BL"));
  end EthanolWaterHighPPCSAFT;

  model EthaneHeptaneCubic  extends TestModel2(Medium(subsNames = "Ethane,Heptane_n", eosType = "PR", mixRule = "VdW", activityModel = "NRTL"), initialT = 403.15, finalT = 403.15, initialP = 70e5, finalP = 70e5, initialZ1 = 0.0, finalZ1 = 1.0);
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation(info = "<html><head></head><body>Here you can see the limitations at high pressure of the equilibrium constants inizialization using the Wilson equation.</body></html>"));
  end EthaneHeptaneCubic;

  model EthaneHeptanePCSAFT
    extends EthaneHeptaneCubic(Medium(eosType="PCSAFT"));
  end EthaneHeptanePCSAFT;
  
  model AcetoneHeptaneCubic  extends TestModel2(Medium(subsNames = "Acetone,Heptane_n", eosType = "SRK", mixRule = "Reid", activityModel = "NRTL"), initialT = 400.15, finalT = 400.15, initialP = 20e5, finalP = 20e5, initialZ1 = 0.0, finalZ1 = 1.0);
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation(info = "<html><head></head><body>Here you can see the limitations of the equilibrium constants inizialization using the Wilson equation.</body></html>"));
  end AcetoneHeptaneCubic;

  model AcetoneHeptane2
    extends AcetoneHeptaneCubic(Medium(eosType="PR",mixRule="MHV2",activityModel="UNIFACpsrk"));
  end AcetoneHeptane2;

  model EGWaterCubic
  extends TestModel2(Medium(subsNames = "EthyleneGlycol,WaterRef", eosType = "PR", mixRule = "VdW", activityModel = "NRTL"), initialT = 440, finalT = 440, initialP = 1.01325e5, finalP = 1.01325e5, initialZ1 = 0.001, finalZ1 = 0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));  
  end EGWaterCubic;

  model EGWaterPCSAFT
    extends EGWaterCubic(Medium(eosType="PCSAFT", mixRule="BL"));
  end EGWaterPCSAFT;

  model ChloroformAcetoneCubic
    extends TestModel2(Medium(subsNames = "Chloroform,Acetone", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 318.15, finalT = 318.15, initialP = 1.01325e5, finalP = 1.01325e5, initialZ1 = 0.001, finalZ1 = 0.999);
  end ChloroformAcetoneCubic;

  model ChloroformAcetonePCSAFT
    extends ChloroformAcetoneCubic(Medium(eosType = "PCSAFT", mixRule = "IndAssoc"));
    algorithm
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end ChloroformAcetonePCSAFT;

  model EthaneButaneCubic
    extends TestModel2(Medium(subsNames = "Ethane,Butane_n", eosType = "PR", mixRule = "MHV2", activityModel = "UNIFACdort"), initialT = 250, finalT = 250, initialP = 2e5, finalP = 2e5, initialZ1 = 0.0001, finalZ1 = 0.9999);
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
    
  model MethanolPentane_nCubic
  extends TestModel2(Medium(subsNames = "Methanol,Pentane_n", eosType = "PR", mixRule = "LCVM", activityModel = "NRTL"), initialT = 397.7, finalT = 397.7, initialP = 14e5, finalP = 14e5, initialZ1=0.0001, finalZ1=0.9999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end MethanolPentane_nCubic;

  model MethanolPentane_nPCSAFT
  extends MethanolPentane_nCubic(Medium(eosType = "PCSAFT", mixRule = "None", activityModel = "None"));
  end MethanolPentane_nPCSAFT;

  model MethanolCyclohexaneCubic
  extends TestModel2(Medium(subsNames = "Methanol2,Cyclohexane", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 293.15, finalT = 293.15, initialP = 0.65e5, finalP = 0.65e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end MethanolCyclohexaneCubic;

  model MethanolCyclohexanePCSAFT
  extends MethanolCyclohexaneCubic(Medium(eosType = "PCSAFT", mixRule = "BL", activityModel = "None"));
  annotation(
      Documentation(info = "<html><head></head><body>If you use Methanol2, you will see the the improvement that polar PCSAFT was doing in Methanol.</body></html>"));
  end MethanolCyclohexanePCSAFT;
  
  model HeptaneMethanolCubic
  extends TestModel2(Medium(subsNames = "Heptane_n,Methanol", eosType = "PR", mixRule = "MHV2", activityModel = "UNIFACdort"), initialT = 313.15, finalT = 313.15, initialP = 0.65e5, finalP = 0.65e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end HeptaneMethanolCubic;

  model HeptaneMethanolPCSAFT
  extends HeptaneMethanolCubic(Medium(eosType = "PCSAFT", mixRule = "IndAssoc", activityModel = "None"));
  annotation(
      Documentation(info = "<html><head></head><body>If you use Methanol2, you will see the the improvement that polar PCSAFT was doing in Methanol.</body></html>"));
  end HeptaneMethanolPCSAFT;  
  
  
  model CyclohexaneAceticAcidCubic
  extends TestModel2(Medium(subsNames = "Cyclohexane,AceticAcid", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 323.15, finalT = 323.15, initialP = 1e5, finalP = 1e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end CyclohexaneAceticAcidCubic;
  
  model CyclohexaneAceticAcidPCSAFT
  extends CyclohexaneAceticAcidCubic(Medium(eosType = "PCSAFT", mixRule = "BL", activityModel = "None"));
  end CyclohexaneAceticAcidPCSAFT;
  
  model WaterAceticAcidCubic
  extends TestModel2(Medium(subsNames = "WaterRef,AceticAcid", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 400, finalT = 400, initialP = 2.734e5, finalP = 2.734e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end WaterAceticAcidCubic;
  
  model WaterAceticAcidPCSAFT
  extends WaterAceticAcidCubic(Medium(eosType = "PCSAFT", mixRule = "BL", activityModel = "None"));
  end WaterAceticAcidPCSAFT;

  model AcrylicAceticAcidCubic
  extends TestModel2(Medium(subsNames = "AcrylicAcid,AceticAcid", eosType = "PR", mixRule = "MHV2", activityModel = "UNIFACstd"), initialT = 400, finalT = 400, initialP = 1.0044e5, finalP = 1.0044e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01),
  Documentation);
  end AcrylicAceticAcidCubic;
  
  model AcrylicAceticAcidPCSAFT
  extends AcrylicAceticAcidCubic(Medium(eosType = "PCSAFT", mixRule = "None", activityModel = "None"));
  end AcrylicAceticAcidPCSAFT;

  model EthanolHeptaneCubic
  extends TestModel2(Medium(subsNames = "Ethanol,Heptane_n", eosType = "PR", mixRule = "LCVM", activityModel = "UNIFACpsrk"), initialT = 293.15, finalT = 353.15, initialP = 1.0e5, finalP = 20.0e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolHeptaneCubic;
  
  model WaterDioxaneCubic
  extends TestModel2(Medium(subsNames = "WaterRef,Dioxane_1_4", eosType = "PR", mixRule = "LCVM", activityModel = "UNIFACpsrk", viscMixRule="Teja"), initialT = 303.15, finalT = 303.15, initialP = 1.0e5, finalP = 1.0e5, initialZ1=0.001, finalZ1=0.999);
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end WaterDioxaneCubic;
  
  model WaterDioxane2
  extends WaterDioxaneCubic(Medium(viscMixRule="McAllister4"));
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end WaterDioxane2;

  model AcetoneEthanolWaterCubic
  extends TestModelBase(Medium(subsNames = "Acetone,Ethanol,WaterRef", eosType = "PR", mixRule = "MHV2", activityModel = "NRTL"), initialT = 300, finalT = 350, initialP = 1.9e5, finalP = 1.5e5, Z={0.5,0.25,0.25});
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end AcetoneEthanolWaterCubic;
  
  model n3i4n4i5n5
  extends TestModelBase(Medium(subsNames = "Propane,Isobutane,Butane_n,Isopentane,Pentane_n", eosType = "PR", mixRule = "VdW", activityModel = "UNIFACpsrk"), initialT = 300, finalT = 350, initialP = 8.5e5, finalP = 8.5e5, Z={0.05,0.15,0.25,0.2,0.35});
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end n3i4n4i5n5;
  
  

  model TestLiquidStateCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNames = "EthyleneGlycol,WaterRef", eosType = "PR", mixRule = "MHV2", activityModel = "UNIFACdort", viscMixRule="Andrade");
   parameter Real initialT = 273.15+10;           //93.33;
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
   Medium.ThermodynamicState stateH;
   Medium.Density D(displayUnit="kg/m3");
   Medium.SpecificEnthalpy H;
   Real SS;
   Medium.SpecificHeatCapacity Cp;
   Medium.DynamicViscosity Eta;
   Medium.ThermalConductivity Lambda;
  algorithm
   MM := Medium.molarMasses();
   Z := Medium.massToMoleFractions(Zmass, MM);
   stateL := Medium.setLiquidState_pTX(p, T, Zmass);
   D:=Medium.density(stateL);
   H:=Medium.specificEnthalpy(stateL);
   SS := Medium.velocityOfSound(stateL);
   Cp := Medium.specificHeatCapacityCp(stateL);
   Eta := Medium.dynamicViscosity(stateL);
   Lambda:=Medium.thermalConductivity(stateL);
   stateH:=Medium.setState_phX(p,H,Zmass);
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

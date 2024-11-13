within FreeFluids.ExternalMix;

package Tests
  model EthanolWaterCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNum = 2, subsNames = "Ethanol,WaterRef", eosType = "PR", cubicMixRule = "MHV2", activityModel = "UNIFACstd");
   parameter Real initialT = 353.15;//403
   parameter Real finalT = 453.15;//403
   parameter Real initialP = 10e5;//29.2e5;
   parameter Real finalP = 10e5;//29.2e5;
   parameter Real initialZ1 = 0.5;
   parameter Real finalZ1 = 0.5;
   Real Z[2](start = {initialZ1, 1 - initialZ1}) "mole fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[2];
   Real Zmass[2] "mass fraction of substances";
   Medium.ThermodynamicState stateP;
   Real gamma[2];
   Medium.SpecificEnergy gE;
   Medium.ThermodynamicState stateBubbleP "bubble state from T and X";
   Medium.ThermodynamicState stateBubbleT;
   Medium.ThermodynamicState stateDewP;
   Medium.ThermodynamicState stateDewT;
   Medium.SpecificEnergy Hv;
   Real SS;
   Real Cp;
  algorithm
    for i in 1:1 loop
      MM := Medium.molarMasses();
      Zmass := Medium.moleToMassFractions(Z, MM);
      stateP := Medium.setState_pTX(p, T, Z);
      (gamma, gE) := Medium.activityCoefficients_TX(T, Zmass);
      gE := gE/(Z*MM);
      stateBubbleP := Medium.setBubbleState_TX(T, Zmass);
      stateBubbleT := Medium.setBubbleState_pX(p, Zmass);
      stateDewP := Medium.setDewState_TX(T, Zmass);
      stateDewT := Medium.setDewState_pX(p, Zmass);
      Hv := Medium.specificVaporizationHeat(stateBubbleP);
      SS := Medium.velocityOfSound(stateP);
      Cp := Medium.specificHeatCapacityCp(stateP);
    end for;
  equation
    der(T) = finalT - initialT;
    der(p) = finalP - initialP;
    der(Z[1]) = finalZ1 - initialZ1;
    Z[2] = 1.0 - Z[1];
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolWaterCubic;

  model EthanolWaterPCSAFT
    extends EthanolWaterCubic(Medium(eosType="PCSAFT"));
  end EthanolWaterPCSAFT;

  model EthanolWaterHighPCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNum = 2, subsNames = "Ethanol,WaterRef", eosType = "PR", cubicMixRule = "LCVM", activityModel = "UNIFACpsrk");
   parameter Real initialT = 480;
   parameter Real finalT = 480;
   parameter Real initialP = 20e5;
   parameter Real finalP = 20e5;
   parameter Real initialZ1 = 0.0;
   parameter Real finalZ1 = 1.0;
   Real Z[2](start = {initialZ1, 1 - initialZ1}) "mole fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[2];
   Real Zmass[2] "mass fraction of substances";
   Medium.ThermodynamicState stateP;
   Real gamma[2];
   Medium.SpecificEnergy gE;
   Medium.ThermodynamicState stateBubbleP "bubble state from T and X";
   Medium.ThermodynamicState stateBubbleT;
   Medium.ThermodynamicState stateDewP;
   Medium.ThermodynamicState stateDewT;
   Medium.SpecificEnergy Hv;
   Real SS;
   Real Cp;
  algorithm
    for i in 1:1 loop
      MM := Medium.molarMasses();
      Zmass := Medium.moleToMassFractions(Z, MM);
      stateP := Medium.setState_pTX(p, T, Z);
      (gamma, gE) := Medium.activityCoefficients_TX(T, Zmass);
      gE := gE/(Z*MM);
      stateBubbleP := Medium.setBubbleState_TX(T, Zmass);
      stateBubbleT := Medium.setBubbleState_pX(p, Zmass);
      stateDewP := Medium.setDewState_TX(T, Zmass);
      stateDewT := Medium.setDewState_pX(p, Zmass);
      Hv := Medium.specificVaporizationHeat(stateBubbleP);
      SS := Medium.velocityOfSound(stateP);
      Cp := Medium.specificHeatCapacityCp(stateP);
    end for;
  equation
    der(T) = finalT - initialT;
    der(p) = finalP - initialP;
    der(Z[1]) = finalZ1 - initialZ1;
    Z[2] = 1.0 - Z[1];
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthanolWaterHighPCubic;

  model EthanolWaterHighPPCSAFT
    extends EthanolWaterHighPCubic(Medium(eosType="PCSAFT"));
  end EthanolWaterHighPPCSAFT;

  model EthaneHeptaneCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNum = 2, subsNames = "Ethane,Heptane_n", eosType = "PR", cubicMixRule = "HV", activityModel = "UNIFACpsrk");
   parameter Real initialT = 403;//403
   parameter Real finalT = 403;//403
   parameter Real initialP = 29.2e5;//29.2e5;
   parameter Real finalP = 29.2e5;//29.2e5;
   parameter Real initialZ1 = 0.0;
   parameter Real finalZ1 = 1.0;
   Real Z[2](start = {initialZ1, 1 - initialZ1}) "mole fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[2];
   Real Zmass[2] "mass fraction of substances";
   Medium.ThermodynamicState stateP;
   Real gamma[2];
   Medium.SpecificEnergy gE;
   Medium.ThermodynamicState stateBubbleP "bubble state from T and X";
   Medium.ThermodynamicState stateBubbleT;
   Medium.ThermodynamicState stateDewP;
   Medium.ThermodynamicState stateDewT;
   Medium.SpecificEnergy Hv;
   Real SS;
   Real Cp;
  algorithm
    for i in 1:1 loop
      MM := Medium.molarMasses();
      Zmass := Medium.moleToMassFractions(Z, MM);
      stateP := Medium.setState_pTX(p, T, Z);
      (gamma, gE) := Medium.activityCoefficients_TX(T, Zmass);
      gE := gE/(Z*MM);
      stateBubbleP := Medium.setBubbleState_TX(T, Zmass);
      stateBubbleT := Medium.setBubbleState_pX(p, Zmass);
      stateDewP := Medium.setDewState_TX(T, Zmass);
      stateDewT := Medium.setDewState_pX(p, Zmass);
      Hv := Medium.specificVaporizationHeat(stateBubbleP);
      SS := Medium.velocityOfSound(stateP);
      Cp := Medium.specificHeatCapacityCp(stateP);
    end for;
  equation
    der(T) = finalT - initialT;
    der(p) = finalP - initialP;
    der(Z[1]) = finalZ1 - initialZ1;
    Z[2] = 1.0 - Z[1];
    annotation(
     experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.01));
  end EthaneHeptaneCubic;

  model TestLiquidStateCubic
   package Medium = FreeFluids.ExternalMix.ExternalMixMedium(mediumName = "Test", subsNum = 2, subsNames = "EG,WaterRef", eosType = "PR", cubicMixRule = "LCVM", activityModel = "UNIFACpsrk");
   parameter Real initialT = 273.15;
   //403
   parameter Real finalT = 273.15;
   parameter Real initialP = 1e5;
   //29.2e5;
   parameter Real finalP = 1e5;
   //29.2e5;
   parameter Real initialZ1 = 0.01;
   parameter Real finalZ1 = 1.0;
   Real Z[2](start = {initialZ1, 1 - initialZ1}) "mass fraction of substances";
   Medium.Temperature T(start = initialT);
   Medium.AbsolutePressure p(start = initialP);
   Medium.MolarMass MM[2];
   Real Zmass[2] "mass fraction of substances";
   Medium.ThermodynamicState stateL;
   Real SS;
   Real Cp;
algorithm
   for i in 1:1 loop
      MM := Medium.molarMasses();
      Zmass := Medium.moleToMassFractions(Z, MM);
      stateL := Medium.setLiquidState_pTX(p, T, Z);
  
      SS := Medium.velocityOfSound(stateL);
      Cp := Medium.specificHeatCapacityCp(stateL);
   end for;
  equation
   der(T) = finalT - initialT;
   der(p) = finalP - initialP;
   der(Z[1]) = finalZ1 - initialZ1;
   Z[2] = 1.0 - Z[1];
  end TestLiquidStateCubic;
  
  model TestLiquiStatePCSAFT
    extends TestLiquidStateCubic(Medium(eosType="PCSAFT"));
  end TestLiquiStatePCSAFT;
end Tests;

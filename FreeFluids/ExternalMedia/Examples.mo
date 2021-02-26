within FreeFluids.ExternalMedia;

package Examples

  model dTState_RampTCorr "Makes some direct calculations from T and d, making a ramp on T"
    package Medium
      extends Mediums.ExternalTPMedium(mediumName = "waterG", thermoModel = 13);
      //extends Modelica.Media.Water.WaterIF97_pT;
    end Medium;

    //We define T and P and thermo record from them
    Medium.Temperature rawT;
    Medium.Density rawd;
    Medium.AbsolutePressure rawP;
    Medium.ThermodynamicState state;
    //now some properties of the thermo record
    Medium.Density rawDens;
    Medium.SpecificEnthalpy rawH;
    Medium.SpecificHeatCapacity rawCp;
    Medium.SpecificEntropy rawS;
    //Medium.DynamicViscosity mu;
    //Medium.ThermalConductivity k;
    //Medium.SurfaceTension sigma;
    //now a SaturationProperties record at the same T and some properties of the liquid and gas phases
    Medium.SaturationProperties sat;
    Medium.Density satLiqDens;
    Medium.Density satGasDens;
    Medium.SpecificEnthalpy satLiqH;
    Medium.SpecificEnthalpy satGasH;
    Medium.SpecificEnthalpy satHv;
  initial equation
    rawd = Medium.liquidDensityCorr_pT(4e7, 350);
  equation
    rawT = 350;
//we fix T
    der(rawd) = -24;
//kgr/m3/s
    state = Medium.setState_dT(rawd, rawT) "links T and P with a thermodynamic record in monophasic way (from P and T)";
//retrieve properties
    rawP = Medium.pressure(state);
    rawDens = Medium.density(state);
    rawH = Medium.specificEnthalpy(state);
    rawCp = Medium.specificHeatCapacityCp(state);
    rawS = Medium.specificEntropy(state);
//mu=Medium.dynamicViscosity(state);
//k=Medium.thermalConductivity(state);
//sigma=Medium.surfaceTension(state);
//link the saturated record with T
    sat = Medium.setSat_T(rawT);
//we check some properties at saturation to see the difference
    satLiqDens = Medium.bubbleDensity(sat);
    satGasDens = Medium.dewDensity(sat);
    satLiqH = Medium.bubbleEnthalpy(sat);
    satGasH = Medium.dewEnthalpy(sat);
    satHv = Medium.vaporizationEnthalpy(sat);
  end dTState_RampTCorr;

  model TemperatureRampTest "Temperature ramp with a thermodynamic state linked to P and T, and another to d and T. works OK"
    package Medium
      extends Mediums.ExternalTPMedium(mediumName = "acetone", thermoModel = 3);
      //extends Modelica.Media.Water.WaterIF97_pT;
    end Medium;

    //We make a ramp on T at the following speed
    parameter Real dT = 200 "K per second";
    //parameter Real P=101325;
    Medium.Temperature T;
    //We will calculate the vapor pressure inside de sat record and we will follow density and enthalpy of both saturated phases
    Medium.SaturationProperties sat "Will give us the saturation pressure at the given T";
    Medium.Density satLiqDens, satGasDens;
    Medium.Density satLiqDensCorr;
    Medium.SpecificEnthalpy satLiqH;
    Medium.SpecificEnthalpy satGasH;
    Medium.SpecificEntropy satLiqS;
    Medium.SpecificEntropy satGasS;
    Medium.SpecificEnthalpy satHdif;
    Medium.SpecificEnthalpy satHvCorr;
    //In order to get the liquid saturated viscosity we need to construct a liquid saturated ThermodynamicState record
    Medium.ThermodynamicState bubbleState " A thermodynamicState to be constructed for the saturated condition";
    //Now we test how to get data using the thermodynamicState record
    Medium.Temperature bubT;
    Medium.AbsolutePressure bubP;
    //Medium.DynamicViscosity bubLiqVisc "Liquid Viscosity in saturated conditions";
    //Medium.DynamicViscosity bubGasVisc "Gas Viscosity in saturated conditions";
    Medium.Density bubLiqDens "Will be the density in saturated conditions obtained from the ThermodynamicState insted than from the SaturationProperties";
    //I have added a function to get the data directly from the EOS
    //Medium.Density satLiqDensEos,satGasDensEos "Will be the densities of liquid and gas phases at saturated conditions by EOS";
    Medium.Density mixD "A given density, to which we will give a value that will obly a two phases situation";
    //Now we construct a thermodynamic state biphasic, from a density that oblies to biphasic, at least at the beginnings
    Medium.ThermodynamicState mixState "A thermodynamicState that corresponds to the given T and d mix. It is interesting to look at gas fraction";
    Medium.SpecificEnthalpy mixH;
    Medium.SpecificHeatCapacity mixCp;
    Medium.SpecificEntropy mixS;
  initial equation
    T = 273.15;
  equation
    der(T) = dT;
    sat = Medium.setSat_T(T);
    satLiqDens = Medium.bubbleDensity(sat);
    satGasDens = Medium.dewDensity(sat);
    satLiqDensCorr = Medium.liquidDensityCorr_pT(sat.psat, sat.Tsat);
    satLiqH = Medium.bubbleEnthalpy(sat);
    satGasH = Medium.dewEnthalpy(sat);
    satHdif = satGasH - satLiqH;
    satHvCorr = Medium.vaporizationEnthalpyCorr(sat.Tsat);
    satLiqS = Medium.bubbleEntropy(sat);
    satGasS = Medium.dewEntropy(sat);
    bubbleState = Medium.setBubbleState(sat);
    bubT = Medium.temperature(bubbleState);
    bubP = Medium.pressure(bubbleState);
//bubLiqVisc=Medium.liquidDynamicViscosity(bubbleState);
//bubGasVisc=Medium.gasDynamicViscosity(bubbleState);
    bubLiqDens = Medium.density(bubbleState);
//(satLiqDensEos,satGasDensEos)=Medium.densitiesEOS_pT(sat.psat,T);
    mixD = 600;
    mixState = Medium.setState_dT(mixD, T);
    mixH = Medium.specificEnthalpy(mixState);
    mixCp = Medium.specificHeatCapacityCp(mixState);
    mixS = Medium.specificEntropy(mixState);
  end TemperatureRampTest;

  model DensityRampTest "Density ramp with a thermodynamic state linked to d and T. works OK"
    package Medium
      extends Mediums.ExternalTPMedium(mediumName = "water", thermoModel = 13);
    end Medium;

    parameter Real dD = -4 "Joules per second";
    //Medium.AbsolutePressure p;
    Medium.Temperature T;
    Medium.Density mixD "A given density, to which we will give a value that will obly a two phases situation";
    //Now we construct a thermodynamic state biphasic, from a density that oblies to biphasic, at least at the beginnings
    Medium.ThermodynamicState mixState "A thermodynamicState that corresponds to the given T and d mix. It is interesting to look at gas fraction";
    Medium.SpecificEnthalpy mixH;
    Medium.SaturationProperties sat;
    Medium.ThermodynamicState bubbleState;
    Medium.AbsolutePressure mixP;
  initial equation
    mixD = 977;
  equation
    der(mixD) = dD;
    T = 348.15;
    mixState = Medium.setState_dT(mixD, T);
    mixP = Medium.pressure(mixState);
    mixH = Medium.specificEnthalpy(mixState);
    sat = Medium.setSat_T(T);
    bubbleState = Medium.setBubbleState(sat);
  end DensityRampTest;

  model DensityRampTest2 "Density ramp with a thermo record defined linked to T and P. Works OK"
    package Medium
      extends Mediums.ExternalTPMedium(mediumName = "water", thermoModel = 13);
    end Medium;

    parameter Real dD = -0.1 "kgr/m3 s";
    //Medium.AbsolutePressure p;
    Medium.Temperature T;
    Medium.Density mixD "A given density, to which we will give a value that will obly a two phases situation";
    //Now we construct a thermodynamic state biphasic, from a density that oblies to biphasic, at least at the beginnings
    Medium.ThermodynamicState mixState "A thermodynamicState that corresponds to the given T and d mix. It is interesting to look at gas fraction";
    Medium.SpecificEnthalpy mixH;
    Medium.SaturationProperties sat;
    Medium.ThermodynamicState bubbleState;
    Medium.AbsolutePressure mixP;
  initial equation
    mixD = 0.2;
  equation
    der(mixD) = dD;
    T = 348.15;
    mixState = Medium.setState_dT(mixP, T);
    mixD = Medium.density(mixState);
    mixH = Medium.specificEnthalpy(mixState);
    sat = Medium.setSat_T(T);
    bubbleState = Medium.setBubbleState(sat);
  end DensityRampTest2;

  model FluidTesting
    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
    //package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model";
    parameter Medium.AbsolutePressure p = 50.0e5;
    Medium.Temperature T(start = 273.17);
    Medium.ThermodynamicState state0;
    Medium.ThermodynamicState stateH "state from p,h";
    Medium.ThermodynamicState stateS "state from p,s";
    Medium.ThermodynamicState stateD "state from d,T";
    Medium.SaturationProperties sat;
    //Real MW,Tc,Zc;
    //Real MM;
    Real Rho, H, S, Mu, Th, Sigma;
    Real Beta "isobaric expansion coefficient";
    Real Kappa "isothermal compressibility";
    Real Cp;
    Real Cv;
    Real Gamma;
    Real SS "speed of sound";
  algorithm
//MW:=Medium.data.MW;
//Tc:=Medium.data.Tc;
//Zc:=Medium.data.Zc;
    state0 := Medium.setState_pTX(p, T, fill(0, 0));
    sat := Medium.setSat_T(T);
/*Pc:=Medium.data.Pc;
  Vp:=121325;*/
//MM:=Medium.molarMass(state0);
    H := Medium.specificEnthalpy(state0);
    S := Medium.specificEntropy(state0);
    Rho := Medium.density(state0);
    Beta := Medium.isobaricExpansionCoefficient(state0);
    Kappa := Medium.isothermalCompressibility(state0);
    Cp := Medium.specificHeatCapacityCp(state0);
    Cv := Medium.specificHeatCapacityCv(state0);
    Gamma := Cp / Cv;
    SS := Medium.velocityOfSound(state0);
    Sigma := Medium.surfaceTension(sat);
    Mu := Medium.dynamicViscosity(state0);
//CpL := Medium.cp(state0);
    Th := Medium.thermalConductivity(state0);
    stateH := Medium.setState_phX(p, H, fill(0, 0));
    stateS := Medium.setState_psX(p, S, fill(0, 0), 0);
    stateD := Medium.setState_dTX(Rho, T, fill(0, 0), 0);
  equation
    der(T) = 250.0;
//Vp=Medium.vp(Tb);
  end FluidTesting;

  model TestVfinder
    //extends FluidTesting(redeclare package Medium = FreeFluids.ExternalMedia.Water);
    package Medium = FreeFluids.ExternalMedia.Fluids.Eicosane(thermoModel = 2);
    parameter Real initD = 100, finalD = 600.0;
    parameter Real T=780;
    Real d(start = initD);
    Real p;
  algorithm
    p := Medium.pressureEOS_dT(d, T);
  equation
    der(d) = finalD - initD;
  end TestVfinder;

  model TestVfinderCubic
    //extends FluidTesting(redeclare package Medium = FreeFluids.ExternalMedia.Water);
    extends TestVfinder(Medium.thermoModel = 1);
  end TestVfinderCubic;

  model TestVfinderSW
    //extends FluidTesting(redeclare package Medium = FreeFluids.ExternalMedia.Water);
    extends TestVfinder(Medium.thermoModel = 3);
  end TestVfinderSW;

  model TestVfinder2
    //extends FluidTesting(redeclare package Medium = FreeFluids.ExternalMedia.Water);
    package Medium = FreeFluids.ExternalMedia.Acetone(thermoModel = 2);
    Real T;
    Real ld, gd;
    Real p;
  algorithm
    for i in 1:300 loop
      T := 230 + 273.15;
      p := 40e5;
      (ld, gd) := Medium.densitiesEOS_pT(p, T);
    end for;
  end TestVfinder2;
end Examples;

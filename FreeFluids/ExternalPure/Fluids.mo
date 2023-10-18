within FreeFluids.ExternalPure;
package Fluids
  package Acetone
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Acetone", fluidK(casRegistryNumber = "67-64-1", description = "Multiparameter:Lemmon&Span 2006. PCSAFT: Solms 2004. PRMC: Kleiman 2002. Cp0: Wilhoit", molarMass = 0.058079, criticalTemperature = 5.081000e+02, criticalPressure = 4.700000e+06), refName = "Propane", onePhase = false, thermoModel = 3, refState = 2);
  end Acetone;

  package Air
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Air",
 fluidK(casRegistryNumber = "132259-10-0", description = "SW: Lemmon 200, PR M.C. C.Trujillo 2023, Cp0: Perry", molarMass = 2.896546e-02, criticalTemperature = 1.326312e+02, criticalPressure = 3.785020e+06),
 final onePhase=false, thermoModel=1, refState=2, ancillaries=3);
  end Air;

  package Ammonia
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Ammonia", fluidK(casRegistryNumber = "7664-41-7", description = "Multiparameter: Tillner-Roth. PCSAFT: Mejbri 2005, Cubic: PRMC C.Trujillo 2019. Cp0: Willhoit", molarMass = 1.703026e-02, criticalTemperature = 4.054000e+02, criticalPressure = 1.133300e+07), final onePhase = false, thermoModel = 3, refState = 2);
  end Ammonia;

  package Argon
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Argon",
 fluidK(casRegistryNumber = "7440-37-1", description = "SW: Tegeler 1999, PCSAFT: Gross and Sadowski 2001, cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper", molarMass = 3.995000e-02, criticalTemperature = 1.506870e+02, criticalPressure = 4.863000e+06),
 final onePhase=false, thermoModel=1, refState=2, ancillaries=3);
  end Argon;

  package Benzene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Benzene",
   fluidK(casRegistryNumber = "71-43-2", description = "Multiparameter:Thol. PCSAFT: Gross 2001. Cubic:PRMC. Cp0: Wilhoit", molarMass = 7.811400e-02, criticalTemperature = 5.620200e+02, criticalPressure = 4.894000e+06),
   final onePhase=false, thermoModel=3, refState=2);
  end Benzene;
  
  package Butane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Butane_n", fluidK(casRegistryNumber = "106-97-8", description = "Multiparameter:Buecker&Wagner 2006. PCSAFT:Gross&Sadowski 2001. Cubic:PRMC, C.Trujillo 2019", molarMass = 0.05812, criticalTemperature = 4.251250e+02, criticalPressure = 3.796000e+06), onePhase = false, thermoModel = 3, refState = 2);
  end Butane_n;

  package Butanol_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Butanol_n", fluidK(casRegistryNumber = "71-36-3", description = "Multiparameter: none. PCSAFT: 2B C.Trujillo 2019. Cubic: PRMC. Cp0: DIPPR107", molarMass = 7.414000e-02, criticalTemperature = 5.630000e+02, criticalPressure = 4.414000e+06), final onePhase = false, thermoModel = 2, refState = 2);
  end Butanol_n;

  package Butene_1
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Butene_1",
 fluidK(casRegistryNumber = "106-98-9", description = "SW: Lemmon 2005, PCSAFT: Gross and Sadowski 2001, PR M.C. C.Trujillo 2023, Cp0: Cooper", molarMass = 5.610800e-02, criticalTemperature = 4.192900e+02, criticalPressure = 4.005100e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Butene_1;

  package Butene_cis_2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Butene_cis_2",
 fluidK(casRegistryNumber = "590-18-1", description = "SW: Lemmon 2005, Cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 5.610800e-02, criticalTemperature = 4.357500e+02, criticalPressure = 4.225500e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Butene_cis_2;

  package CarbonMonoxide
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "CarbonMonoxide",
 fluidK(casRegistryNumber = "630-08-0", description = "SW: Lemmon and Span 2006, PCSAFT: Gross and Sadowski 2001, PR M.C. C.Trujillo 2023, Cp0: Cooper", molarMass = 2.801800e-02, criticalTemperature = 1.328600e+02, criticalPressure = 3.494000e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="N2", ancillaries=3);
  end CarbonMonoxide;

  package CarbonylSulfide
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "CarbonylSulfide",
 fluidK(casRegistryNumber = "463-58-1", description = "SW: Lemmon and Span 2006, Cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 6.007510e-02, criticalTemperature = 3.787700e+02, criticalPressure = 6.370000e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end CarbonylSulfide;

  package CO2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "CO2", fluidK(casRegistryNumber = "124-38-9", description = "Multiparameter: Span and Wagner 1996. PCSAFT and PRMC: C.Trujillo from GERG2004. Cp0:Wilhoit", molarMass = 4.400980e-02, criticalTemperature = 3.041282e+02, criticalPressure = 7.377730e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end CO2;

  package CO2b
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "CO2b", fluidK(casRegistryNumber = "124-38-9", description = "SW: GERG2004. PCSAFT and Cubic: C.Trujillo from GERG2004. Cp0:Jaeschke", molarMass = 4.400980e-02, criticalTemperature = 3.041282e+02, criticalPressure = 7.377300e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end CO2b;

  package Cyclohexane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Cyclohexane",
 fluidK(casRegistryNumber = "110-82-7", description = "SW:Zhou 2018, SAFT: Gross and Sadowski 2001, Cubic: PR M.C. C.Trujillo 2018, Cp0: Cooper", molarMass = 8.418000e-02, criticalTemperature = 5.536000e+02, criticalPressure = 4.082400e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Cyclohexane;

  package Cyclopentane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Cyclopentane",
 fluidK(casRegistryNumber = "287-92-3", description = "SW: Gedanitz 2015, PCSAFT: Gross and Sadowski 2001, Cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper", molarMass = 7.013400e-02, criticalTemperature = 5.117200e+02, criticalPressure = 4.571200e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Cyclopentane;

  package D4
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "D4",
 fluidK(casRegistryNumber = "556-67-2", description = "SW: Thol 2015. Cubic: PR M.C. C.Trujillo 2018, Cp0:Cooper, Transport: Nist coefficients", molarMass = 2.966160e-01, criticalTemperature = 5.865000e+02, criticalPressure = 1.347215e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="N2", ancillaries=3);
  end D4;

  package D5
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "D5",
 fluidK(casRegistryNumber = "541-02-6", description = "SW: Thol 2019, Cubic: PR M.C. C.Trujillo 2023, Cp0: DIPPR 100, Transport: Nist coefficients", molarMass = 3.707700e-01, criticalTemperature = 6.183000e+02, criticalPressure = 1.077700e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="N2", ancillaries=3);
  end D5;

  package D6
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "D6",
 fluidK(casRegistryNumber = "540-97-6", description = "SW: Colonna 2008, Cubic: PR M.C. C.Trujillo 2023, Cp0: DIPPR 107, Transport: Nist coefficients", molarMass = 4.450200e-01, criticalTemperature = 6.457800e+02, criticalPressure = 9.610000e+05),
 final onePhase=false, thermoModel=3, refState=2, refName="N2", ancillaries=3);
  end D6;
  
  package Dichlorodifluormethane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Dichlorodifluormethane", fluidK(casRegistryNumber = "75-71-8", description = "Multiparameter:Marx et alt. PCSAFT:C.Trujillo 2020, Cubic:PRMC. Cp0:Cooper", molarMass = 1.209130e-01, criticalTemperature = 3.851200e+02, criticalPressure = 4.136100e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Dichlorodifluormethane;

  package Dichloroethane_1_2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Dichloroethane_1_2",
 fluidK(casRegistryNumber = "107-06-2", description = "SW: Thol 2015, SAFT non assoc C.Trujillo 2023, Cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 9.895900e-02, criticalTemperature = 5.616000e+02, criticalPressure = 5.227585e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end Dichloroethane_1_2;

  package DiethylEther
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "DiethylEther",
 fluidK(casRegistryNumber = "60-29-7", description = "SW: Thol 2014, PCSAFT Gross and Sadowski 2001, Cubic: PR M.C. C.Trujillo 2023, Cp0:Wilhoit, Transport: Nist coefficients", molarMass = 7.412000e-02, criticalTemperature = 4.667000e+02, criticalPressure = 3.649017e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end DiethylEther;

  package DimethylCarbonate
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "DimethylCarbonate",
 fluidK(casRegistryNumber = "616-38-6", description = "SW: Zhou 2011, Cubic: PR M.C. C.Trujillo 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 9.008400e-02, criticalTemperature = 5.570000e+02, criticalPressure = 4.908800e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end DimethylCarbonate;

  package DimethylEther
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "DimethylEther",
 fluidK(casRegistryNumber = "115-10-6", description = "SW: Wu 2011, PCSAFT: Gross and Sadowski 2001, Cubic: PR M.C. C.Trujillo 2018, Cp0: Cooper, Transport: Nist coefficients", molarMass = 4.607000e-02, criticalTemperature = 4.003780e+02, criticalPressure = 5.336800e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end DimethylEther;

  package EG
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "EG", fluidK(casRegistryNumber = "107-21-1", description = "Multiparameter: none. PPCSAFT_GV: 2B C.Trujillo 2019. Cubic: PRMC. Cp0: Wilhoit", molarMass = 6.207000e-02, criticalTemperature = 7.200000e+02, criticalPressure = 8.200000e+06), final onePhase = false, thermoModel = 2, refState = 2);
  end EG;

  package Eicosane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Eicosane_n", fluidK(casRegistryNumber = "112-95-8", description = "SW: none. PCSAFT: Gross. Cubic: PR. Cp0: Wilhoit ", molarMass = 2.825520e-01, criticalTemperature = 7.680000e+02, criticalPressure = 1.070000e+06), final onePhase = false, thermoModel = 2, refState = 2);
  end Eicosane_n;

  package Ethane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Ethane", fluidK(casRegistryNumber = "74-84-0", description = "Multiparameter: Buecker-Wagner 2006. PCSAFT: Gross-Sadowski 2001. Cubic:SRKMC Chemsep. Cp0: Cooper", molarMass = 3.006904e-02, criticalTemperature = 3.053220e+02, criticalPressure = 4.872200e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Ethane;

  package Ethanol
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Ethanol", fluidK(casRegistryNumber = "64-17-5", description = "Multiparameter: Schroeder 2014. PPCSAFT-GV: C.Trujillo 2018. PRMC: C.Trujillo 2018. Cp0: Wilhoit", molarMass = 4.606844e-02, criticalTemperature = 5.147100e+02, criticalPressure = 6.268000e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Ethanol;

  package Ethylbenzene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Ethylbenzene",
 fluidK(casRegistryNumber = "100-41-4", description = "SW: Zhou 2012, Cubic PR M.C. CTrujillo 2023, Cp0:Cooper, Visc.: Nist coefficients.", molarMass = 1.061600e-01, criticalTemperature = 6.171200e+02, criticalPressure = 3.622400e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end Ethylbenzene;

  package Ethylene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Ethylene",
 fluidK(casRegistryNumber = "74-85-1", description = "SW: Smukala 2000, PCSAFT non assoc.C.Trujillo 2023, Cubic PR M.C. 2023, Cp0: Cooper 2023", molarMass = 2.805400e-02, criticalTemperature = 2.823500e+02, criticalPressure = 5.041800e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Ethylene;

  package EthyleneOxide
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "EthyleneOxide",
 fluidK(casRegistryNumber = "75-21-8", description = "SW: thol 2015, PPCSAFT G.V. C.Trujillo 2023, Cubic: PR M.C. 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 4.405260e-02, criticalTemperature = 4.689200e+02, criticalPressure = 7.304686e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end EthyleneOxide;

  package Helium
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Helium",
 fluidK(casRegistryNumber = "7440-59-7", description = "SW: Ortiz-Vega 2010, Cubic: PR78, Cp0:Cooper", molarMass = 4.002602e-03, criticalTemperature = 5.195300e+00, criticalPressure = 2.276000e+05),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Helium;

  package Heptane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Heptane_n", fluidK(casRegistryNumber = "142-82-5", description = "Multiparameter: Span and Wagner 2003. PCSAFT: Gross 2001. Cubic: PRMC. Cp0: Jaeschke", molarMass = 1.002020e-01, criticalTemperature = 5.401300e+02, criticalPressure = 2.736000e+06), final onePhase = false, thermoModel = 3, refState = 2, ancillaries=3);
  end Heptane_n;

  package Hexane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Hexane_n", fluidK(casRegistryNumber = "110-54-3", description = "Multiparameter: Spand and Wagner 2003. PCSAFT: Gross and Sadowski. Cubic: PRMC. Cp0:Jaeschke", molarMass = 8.617536e-02, criticalTemperature = 5.078200e+02, criticalPressure = 3.034000e+06), final onePhase = false, thermoModel = 3, refState = 2, ancillaries=3);
  end Hexane_n;

  package Hydrogen
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Hydrogen",
 fluidK(casRegistryNumber = "1333-74-0", description = "SW: Leachman 2009, PCSAFT: IVC-SEP DTU, Cubic: PR, Cp0: Cooper", molarMass = 2.016000e-03, criticalTemperature = 3.314500e+01, criticalPressure = 1.296400e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Hydrogen;  

  package HydrogenChloride
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "HydrogenChloride",
 fluidK(casRegistryNumber = "7647-01-0", description = "SW: Thol 2018, PPCSAFT G.V.: C.Trujillo 2023, Cubic PR M.C. 2023, Cp0: Cooper, Transport: Nist coefficients", molarMass = 3.646000e-02, criticalTemperature = 3.245500e+02, criticalPressure = 8.288161e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end HydrogenChloride;

  package HydrogenSulfide
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "HydrogenSulfide",
 fluidK(casRegistryNumber = "7783-06-4", description = "SW: Lemmon and Span 2006, PCSAFT non assoc: IVC_SEP DTU, Cubic PR M.C. 2018, Cp0: Cooper, Transport: Nist Coefficients", molarMass = 3.408088e-02, criticalTemperature = 3.731000e+02, criticalPressure = 9.000000e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end HydrogenSulfide;

  package Isobutane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Isobutane", fluidK(casRegistryNumber = "75-28-5", description = "Multiparameter: Buecker-Wagner 2006. PCSAFT: Gross-Sadowski 2001. Cubic: PRMC. Cp0: Cooper", molarMass = 5.812220e-02, criticalTemperature = 4.078170e+02, criticalPressure = 3.629000e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Isobutane;
  
  package Isobutene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Isobutene",
 fluidK(casRegistryNumber = "115-11-7", description = "SW: Lemmon 2005, PCSAFT non assoc. C.Trujillo 2023, Cubic: PR M.C. 2023, Cp0: Cooper, Transport: Nist coefficients against propane", molarMass = 5.610800e-02, criticalTemperature = 4.180900e+02, criticalPressure = 4.009800e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end Isobutene;

  package Isohexane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Isohexane",
 fluidK(casRegistryNumber = "107-83-5", description = "SW: Lemmon 2006, PCSAFT non assoc. Gross and Sadowski 2001, Cubic: PR M.C. 2023, Cp0: Cooper, Transport Nist coefficients against propane", molarMass = 8.617660e-02, criticalTemperature = 4.977000e+02, criticalPressure = 3.040000e+06),
 final onePhase=false, thermoModel=3, refState=2, refName="Propane", ancillaries=3);
  end Isohexane;

  package Methane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Methane", fluidK(casRegistryNumber = "74-82-8", description = "Multiparameter: Seltzmann-Wagner 1991. PCSAF: Gross-Sadowski 2001. Cubic: SRKMC from Chemsep. Cp0: Cooper", molarMass = 1.604280e-02, criticalTemperature = 1.905640e+02, criticalPressure = 4.599200e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Methane;

  package Methanol
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Methanol",
 fluidK(casRegistryNumber = "67-56-1", description = "Multiparameter: Sun 2004. PPCSAFT-GV: C.Trujillo. Cubic: PRSV", molarMass = 3.205000e-02, criticalTemperature = 5.125000e+02, criticalPressure = 8.215850e+06),
 final onePhase=false, thermoModel=1, refState=2);
  end Methanol;
  
  package Methanol2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Methanol2",
 fluidK(casRegistryNumber = "67-56-1", description = "Multiparameter: Sun 2004. PCSAFT: von Solms. Cubic:PRSV", molarMass = 3.205000e-02, criticalTemperature = 5.125000e+02, criticalPressure = 8.215850e+06),
 final onePhase=false, thermoModel=1, refState=2);
  end Methanol2;

  package MethylEthylKetone
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "MethylEthylKetone",
 fluidK(casRegistryNumber = "78-93-3", description = "SW: none, PCSAFT 2B: C.Trujillo 2023, Cubic SRKMC: Chemsep, Cp0:Wilhoit, Ancillaries: for PCSAF 2B", molarMass = 7.212000e-02, criticalTemperature = 5.356000e+02, criticalPressure = 4.154325e+06),
 final onePhase=false, thermoModel=2, refState=2, ancillaries=2);
  end MethylEthylKetone;

  package MethylMethacrylate
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "MethylMethacrylate",
 fluidK(casRegistryNumber = "80-62-6", description = "SW: none, PCSAFT 2B, PR M.C. and Cp0 C.Trujillo 2023, Ancillaries: for PCSAF 2B", molarMass = 1.001300e-01, criticalTemperature = 5.660000e+02, criticalPressure = 3.680000e+06),
 final onePhase=false, thermoModel=1, refState=2, ancillaries=2);
  end MethylMethacrylate;

  package N2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "N2", fluidK(casRegistryNumber = "7727-37-9", description = "Multiparameter: Span 2001, PCSAFT: Gross 2001. PRMC: Barragan 2002. Cp0: DIPPR107", molarMass = 2.801348e-02, criticalTemperature = 1.261920e+02, criticalPressure = 3.395800e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end N2;

  package O2
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "O2", fluidK(casRegistryNumber = "7782-44-7", description = "Multiparameter: Schmidt 1985. PCSAFT: Economou 2007. Cubic:SRKMC Chemsep. Cp0: Cooper", molarMass = 3.199880e-02, criticalTemperature = 1.545810e+02, criticalPressure = 5.043000e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end O2;

  package Octane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Octane_n",
 fluidK(casRegistryNumber = "111-65-9", description = "SW: GERG 2004. PCSAFT: Gross and Sadowski 2003. Cubic: PR Twu 91. Cp0:Jaeschke", molarMass = 1.142285e-01, criticalTemperature = 5.693200e+02, criticalPressure = 2.497000e+06),
 final onePhase=false, thermoModel=1, refState=2, ancillaries=3);
  end Octane_n;

  package Pentane_n
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Pentane_n",
 fluidK(casRegistryNumber = "109-66-0", description = "SW: GERG2004. PCSAFT: GRoss and Sadowski 2001. Cubic P.R.M.C. Cp0:Jaeschke. Visc: NIST", molarMass = 7.215000e-02, criticalTemperature = 4.697000e+02, criticalPressure = 3.370000e+06),
 final onePhase=false, thermoModel=3, refState=2, ancillaries=3);
  end Pentane_n;

  package Propane
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Propane", fluidK(casRegistryNumber = "74-98-6", description = "Multiparameter: Lemmon 2009. PCSAFT: Gross-Sadowski 2001. Cubic: SRKMC Chemsep. Cp0:Cooper", molarMass = 4.409562e-02, criticalTemperature = 3.698900e+02, criticalPressure = 4.251200e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Propane;

  package R134A
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "R134A", fluidK(casRegistryNumber = "811-97-2", description = "Multiparameter:Tillner 1994. PCSAFT: non assoc.C.Trujillo 2019. Cubic:PRMC, C.Trujillo 2019. Cp0: Wilhoit", molarMass = 1.020310e-01, criticalTemperature = 3.741800e+02, criticalPressure = 4.901200e+06), final onePhase = false, thermoModel = 3, refState = 1);
  end R134A;

  package R410A
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "R410A", fluidK(casRegistryNumber = "R410A.PPF", description = "Multiparameter: Lemmon 2003. PCSAFT: non assoc. C.T.2019. Cubic: none. Cp0: Cooper", molarMass = 7.258540e-02, criticalTemperature = 3.444940e+02, criticalPressure = 4.901200e+06), final onePhase = false, thermoModel = 3, refState = 2, ancillaries=3);
  end R410A;

  package Styrene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Styrene",
 fluidK(casRegistryNumber = "100-42-5", description = "SW: none, PCSAFT: C.Trujillo 2023. Cubic: PRMC C.Trujillo 2023. Cp0: C.Trujillo 2020", molarMass = 1.041600e-01, criticalTemperature = 6.470000e+02, criticalPressure = 3.990000e+06),
 final onePhase=false, thermoModel=1, refState=2, ancillaries=2);
  end Styrene;

  package Toluene
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Toluene", fluidK(casRegistryNumber = "108-88-3", description = "Multiparameter: Lemmon 2006. PCSAFT: Gross 2001. Cubic: SRKMC. Cp0: Wilhoit", molarMass = 9.214020e-02, criticalTemperature = 5.917500e+02, criticalPressure = 4.126000e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Toluene;

  package WaterRef
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "WaterRef", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:IAPWS95. PCSAFT:C.Trujillo 2B max.633K. Cubic:PRMC,C.Trujillo 2018. Cp0:Jaeschke", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), final onePhase = false, thermoModel = 3, refState = 4);
  end WaterRef;

  package Water
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Water", fluidK(casRegistryNumber = "7732-18-5", description = "Multiparameter:GERG2004. PCSAFT: Diamantotis 4C. Cubic: SRKMC C.Trujillo. Cp0:Jaechske", molarMass = 1.801528e-02, criticalTemperature = 6.470960e+02, criticalPressure = 2.206400e+07), onePhase = false, thermoModel = 3, refState = 4, reference_T = 273.15, reference_p = 101325);
  end Water;

  package Xylene_m
    extends FreeFluids.ExternalPure.ExternalMedium(final mediumName = "Xylene_m", fluidK(casRegistryNumber = "108-38-3", description = "Multiparameter: Zhou 2012. PCSAFT: Gross 2001. Cubic: PRSV Proust 1998. Cp0: Cooper", molarMass = 1.061670e-01, criticalTemperature = 6.168900e+02, criticalPressure = 3.534600e+06), final onePhase = false, thermoModel = 3, refState = 2);
  end Xylene_m;
end Fluids;

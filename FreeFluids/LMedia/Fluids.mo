within FreeFluids.LMedia;
  package Fluids
    package Acetone
      extends FreeFluids.LMedia.LMedium(final mediumName = "Acetone", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Acetone}, refState = "IIR");
    end Acetone;

    package Ammonia
      extends FreeFluids.LMedia.LMedium(final mediumName = "Ammonia", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ammonia}, reference_T = 273.15);
    end Ammonia;

    package Butane_n
      extends LMedium(final mediumName = "Butane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Butane_n}, refState = "IIR", reference_T = 273.15);
    end Butane_n;

    package Butanol_n
      extends FreeFluids.LMedia.LMedium(final mediumName = "n-Butanol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Butanol_n});
    end Butanol_n;

    package CO2
      extends FreeFluids.LMedia.LMedium(final mediumName = "Carbon dioxide", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.CO2});
    end CO2;

  package D4
    extends FreeFluids.LMedia.LMedium(final mediumName = "D4", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.D4}, refState="IIR", reference_T=273.15);
  end D4;

  package D5
    extends FreeFluids.LMedia.LMedium(final mediumName = "D5", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.D5}, refState="IIR", reference_T=273.15);
  end D5;

    package DecanoicAcid
      extends FreeFluids.LMedia.LMedium(final mediumName = "DecanoicAcid", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.DecanoicAcid}, refState = "IIR", reference_T = 273.15);
    end DecanoicAcid;

    package Dichlorodifluormethane
      extends FreeFluids.LMedia.LMedium(final mediumName = "Dichlorodifluormethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Dichlorodifluormethane}, reference_T = 273.15, refState = "ASHRAE");
    end Dichlorodifluormethane;

    package Ethane
      extends FreeFluids.LMedia.LMedium(final mediumName = "Ethane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethane}, refState = "IIR", reference_T = 273.15);
    end Ethane;

    package Ethanol
      extends FreeFluids.LMedia.LMedium(final mediumName = "Ethanol", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Ethanol}, reference_T = 273.15);
    end Ethanol;

    package EthyleneGlycol
      extends FreeFluids.LMedia.LMedium(final mediumName = "Ethylene glycol", final singleState = false, p_default = 3.0e5, T_default = 423.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.EG});
    end EthyleneGlycol;

    package Heptane_n
      extends FreeFluids.LMedia.LMedium(final mediumName = "n-Heptane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Heptane_n}, refState = "IIR", reference_T = 273.15);
    end Heptane_n;

    package Hexane_n
      extends FreeFluids.LMedia.LMedium(final mediumName = "n-Hexane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Hexane_n}, refState = "IIR", reference_T = 273.15, p_default = 5.0e5, T_default = 293.15);
    end Hexane_n;

    package Isobutane
      extends FreeFluids.LMedia.LMedium(final mediumName = "Isobutane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataAL.Isobutane}, reference_T = 273.15);
    end Isobutane;

    package MarlothermSH
      extends FreeFluids.LMedia.LMedium(final mediumName = "Marlotherm SH", final singleState = false, p_default = 3.0e5, T_default = 473.15, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.MarlothermSH}, Temperature(min = 0 + 273.15, max = 300 + 273.15, nominal = 260 + 273.15));
    end MarlothermSH;

    package Methane
      extends LMedium(final mediumName = "Methane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Methane}, refState = "IIR", reference_T = 273.15);
    end Methane;

    package N2
      extends FreeFluids.LMedia.LMedium(final mediumName = "Nitrogen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.N2});
    end N2;

  package Nonane_n
    extends FreeFluids.TMedia.TMedium(final mediumName = "Nonane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Nonane_n}, refState="IIR", reference_T=273.15);
  end Nonane_n;
    package O2
      extends FreeFluids.LMedia.LMedium(final mediumName = "Oxygen", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.O2});
    end O2;

    package Pentane_n
      extends LMedium(final mediumName = "Pentane_n", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Pentane_n}, refState = "IIR", reference_T = 273.15);
    end Pentane_n;

  package Propanal
    extends FreeFluids.TMedia.TMedium(final mediumName = "Propanal", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Propanal}, refState="IIR", reference_T=273.15);
  end Propanal;

    package Propane
      extends LMedium(final mediumName = "Propane", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Propane}, refState = "NBP");
    end Propane;

    package R134A
      extends LMedium(final mediumName = "R134A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R134A}, refState = "ASHRAE", reference_T = 273.15);
    end R134A;

    package R410A
      extends LMedium(final mediumName = "R410A", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.R410A}, refState = "ASHRAE", reference_T = 273.15);
    end R410A;

    package ShellS2
      extends FreeFluids.LMedia.LMedium(final mediumName = "Shell S2", final singleState = false, p_default = 3.0e5, T_default = 473.15, reference_T = 273.15, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.ShellS2}, Temperature(min = 0 + 273.15, max = 300 + 273.15, nominal = 260 + 273.15));
    end ShellS2;

  package SulfurDioxide
    extends FreeFluids.TMedia.TMedium(final mediumName = "SulfurDioxide", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.SulfurDioxide}, refState="IIR", reference_T=273.15);
  end SulfurDioxide;
  
    package Toluene
      extends FreeFluids.LMedia.LMedium(final mediumName = "Toluene", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Toluene}, reference_T = 273.15, p_default = 5.0e5, T_default = 393.15);
    end Toluene;

    package Water
      extends FreeFluids.LMedia.LMedium(final mediumName = "Water", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Water}, reference_T = 273.16, refState = "User", p_default = 5.0e5, T_default = 393.15, Temperature(min = 0 + 273.15, max = 600 + 273.15, nominal = 50 + 273.15));
    end Water;

    package Xylene_m
      extends FreeFluids.LMedia.LMedium(final mediumName = "Xylene_m", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Xylene_m}, refState = "NBP", reference_T = 273.15);
    end Xylene_m;
    
  package Xylene_p
    extends FreeFluids.TMedia.TMedium(final mediumName = "Xylene_p", final singleState = false, fluidConstants = {FreeFluids.MediaCommon.MediaDataMZ.Xylene_p}, refState="IIR", reference_T=273.15);
  end Xylene_p;
  end Fluids;

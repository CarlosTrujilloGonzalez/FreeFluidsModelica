within FreeFluids;

package Columns "Columns.mo by Carlos Trujillo
  This file is part of the Free Fluids application
  Copyright (C) 2008-2022  Carlos Trujillo Gonzalez
    
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
  //***COLUMNS***
  //*************

  model DistillationBalance "Performs the total, and one component, mass balances, plus the heat balance at top and bottom."
    parameter Modelica.Units.SI.MolarMass mwL = 0 "light component molecular weight";
    parameter Modelica.Units.SI.MolarMass mwH = 0 "heavy component molecular weight";
    parameter Boolean isBinaryDistillation = true "if true only two components are considered" annotation(
      Dialog(tab = "Feed data"));
    parameter Boolean molarBasis=true "if true, use molar fraction and kmol/h, otherwise mass fraction and kg/h";
    parameter Real ffl(min=0) "fixed feed flow rate of light component in kmol/h or kg/h" annotation(
      Dialog(tab = "Feed data"));
    parameter Real ffh(min=0) "fixed feed flow rate of heavy component in kmol/h or kg/h" annotation(
      Dialog(tab = "Feed data"));
    parameter Modelica.Units.SI.TemperatureDifference fSuperHeating = 0 "feed superheating, subcooling if negative" annotation(
      Dialog(tab = "Feed data"));
    parameter Fraction fGfSat = 0 "gas fraction of the feed, if in saturated conditions" annotation(
      Dialog(tab = "Feed data"));
    parameter Modelica.Units.SI.MolarEnthalpy hvF(displayUnit = "kJ/mol") = 41000 "feed molar vaporization enthalpy" annotation(
      Dialog(tab = "Feed data"));
    parameter Modelica.Units.SI.MolarHeatCapacity cpFl = 0 "feed liquid heat capacity, if subcooled" annotation(
      Dialog(tab = "Feed data"));
    parameter Modelica.Units.SI.MolarHeatCapacity cpFv = 0 "feed vapor heat capacity, if superheated" annotation(
      Dialog(tab = "Feed data"));
    parameter FreeFluids.Types.PerformanceSpecification topFixed = FreeFluids.Types.PerformanceSpecification.FixedComposition "composition, flow, recovery or none are fixed" annotation(
      Dialog(tab = "Top data"));
    parameter Fraction dyl = 0.0 "fixed distillate light fraction" annotation(
      Dialog(tab = "Top data"));
    parameter Real dfl = 0.0 "fixed distillate light flow. kmol/h or kg/h" annotation(
      Dialog(tab = "Top data"));
    parameter Fraction dlr = 0.95 "fixed light recovery in distillate" annotation(
      Dialog(tab = "Top data"));
    parameter Modelica.Units.SI.TemperatureDifference dSubCooling = 0 "distillate subcooling, if total condenser used" annotation(
      Dialog(tab = "Top data"));
    parameter Modelica.Units.SI.MolarEnthalpy hvD(displayUnit = "kJ/mol") = 41000 "distillate molar vaporization enthalpy" annotation(
      Dialog(tab = "Top data"));
    parameter Modelica.Units.SI.MolarHeatCapacity cpDl = 0 "distillate liquid heat capacity" annotation(
      Dialog(tab = "Top data"));
    parameter FreeFluids.Types.PerformanceSpecification bottomFixed = FreeFluids.Types.PerformanceSpecification.FixedComposition "composition, flow, recovery or none are fixed" annotation(
      Dialog(tab = "Bottom data"));
    parameter Fraction bxh = 0.0 "fixed bottom heavy fraction" annotation(
      Dialog(tab = "Bottom data"));
    parameter Real bfh(min=0) = 0.0 "fixed bottom heavy flow. kmol/h or kg/h" annotation(
      Dialog(tab = "Bottom data"));
    parameter Fraction brh = 0.95 "fixed heavy recovery in bottom" annotation(
      Dialog(tab = "Bottom data"));
    parameter Modelica.Units.SI.MolarEnthalpy hvB(displayUnit = "kJ/mol") = 41000 "bottom molar vaporization enthalpy" annotation(
      Dialog(tab = "Bottom data"));
    parameter Boolean fixRatioToMinReflux=false "if true, the ratio reflux/minimum reflux is fixed" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Real minRratio=1.3 "ratio used, if previous is true" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Boolean fixRefluxRatio = false "if true, the ratio reflux/distillate is fixed" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Real rr(min = 0) = 0 "fixed internal reflux/distillate ratio" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Boolean fixCondenserPower = false annotation(
      Dialog(tab = "Operating conditions"));
    parameter Modelica.Units.SI.Power wc(min = 0, displayUnit = "kW") = 0 "fixed power removed at condenser" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Boolean fixReboilerPower = false annotation(
      Dialog(tab = "Operating conditions"));
    parameter Modelica.Units.SI.Power wb(min = 0, displayUnit = "kW") = 0 "fixed power added at reboiler" annotation(
      Dialog(tab = "Operating conditions"));
    parameter Boolean isPartialCondenser = false "if true the distillate leaves as gas" annotation(
      Dialog(tab = "Operating conditions"));
    Modelica.Units.SI.MolarFlowRate F(min = 0, displayUnit = "kmol/h") "feed molar flow rate";
    Modelica.Units.SI.MolarFlowRate D(min = 0, start = ffl/3.6, displayUnit = "kmol/h") "distillate molar flow rate";
    Modelica.Units.SI.MolarFlowRate B(min = 0, start = ffh/3.6, displayUnit = "kmol/h") "bottom purge molar flow rate";
    Modelica.Units.SI.MolarMass Fmw(displayUnit="g/mol") "feed molar mass";
    Modelica.Units.SI.MassFlowRate Fmass(displayUnit = "kg/h") "feed mass flow rate";
    Modelica.Units.SI.MassFlowRate FflMass(displayUnit = "kg/h") "feed light mass flow rate";
    Modelica.Units.SI.MassFlowRate FfhMass(displayUnit = "kg/h") "feed heavy mass flow rate";
    Modelica.Units.SI.MolarMass Dmw(displayUnit="g/mol") "distillate molar mass";
    Modelica.Units.SI.MassFlowRate Dmass(displayUnit = "kg/h") "distillate mass flow rate";
    Modelica.Units.SI.MassFlowRate DflMass(displayUnit = "kg/h") "distillate light mass flow rate";
    Modelica.Units.SI.MassFlowRate DfhMass(displayUnit = "kg/h") "distillate heavy mass flow rate";
    Modelica.Units.SI.MolarMass Bmw(displayUnit="g/mol") "bottom molar mass";
    Modelica.Units.SI.MassFlowRate Bmass(displayUnit = "kg/h") "bottom mass flow rate";
    Modelica.Units.SI.MassFlowRate BfhMass(displayUnit = "kg/h") "bottom heavy mass flow rate";
    Modelica.Units.SI.MassFlowRate BflMass(displayUnit = "kg/h") "bottom light mass flow rate";
    
    Fraction FZl "feed molar fraction of light component";
    Fraction FZlMass "feed mass fraction of light component";  
    Fraction FZh "feed molar fraction of heavy component";
    Fraction FZhMass "feed mass fraction of heavy component";
    Modelica.Units.SI.MolarFlowRate Ffl(min = 0, displayUnit = "kmol/h") "feed molar flow rate of light component";
    Modelica.Units.SI.MolarFlowRate Ffh(min = 0, displayUnit = "kmol/h") "feed molar flow rate of heavy component";
    Fraction DYl "distillate molar fraction of light comp.";
    Fraction DYlMass "distillate mass fraction of light comp.";
    Fraction DYh "distillate molar fraction of heavy comp.";
    Fraction DYhMass "distillate mass fraction of heavy comp.";
    Modelica.Units.SI.MolarFlowRate Dfl(min = 0, start=1, displayUnit = "kmol/h") "distillate molar flow rate of light component";
    Modelica.Units.SI.MolarFlowRate Dfh(min = 0, start=1, displayUnit = "kmol/h") "distillate molar flow rate of heavy component";
    Fraction BXl "bottom molar fraction of light comp.";
    Fraction BXlMass "bottom mass fraction of light comp.";
    Fraction BXh "bottom molar fraction of heavy comp.";
    Fraction BXhMass "bottom mass fraction of heavy comp.";
    Modelica.Units.SI.MolarFlowRate Bfl(min = 0, start=1, displayUnit = "kmol/h") "bottom molar flow rate of light component";
    Modelica.Units.SI.MolarFlowRate Bfh(min = 0, displayUnit = "kmol/h") "bottom molar flow rate of heavy component";
    Fraction Fgf "feed gas fraction";
    Real Q "feed enthalpy parameter";
    Real RR(min = 0) "internal reflux ratio, after heating of subcooled reflux";
    Real RRext(min = 0) "external reflux ratio, as observed at condenser";
    Real BR(min = 0) "boil ratio";
    Modelica.Units.SI.Power Wc(min = 0, displayUnit = "kW") "heat removed at condenser";
    Modelica.Units.SI.Power Wb(min = 0, displayUnit = "kW") "heat added at the reboiler";
    Modelica.Units.SI.MolarFlowRate Le(min = 0, displayUnit = "kmol/h") "liquid molar flow rate in the enrichment section";
    Modelica.Units.SI.MolarFlowRate Ve(min = 0, displayUnit = "kmol/h") "gas molar flow rate in the enrichment section";
    Modelica.Units.SI.MolarFlowRate Ls(min = 0, displayUnit = "kmol/h") "liquid molar flow rate in the stripping section";
    Modelica.Units.SI.MolarFlowRate Vs(min = 0, displayUnit = "kmol/h") "gas molar flow rate in the stripping section";
  equation
//global balance
    F = D + B;
    Ffl = Dfl + Bfl;
    if isBinaryDistillation == true then
      DYh = 1 - DYl;
      BXl = 1 - BXh;
    else
      Ffh = Dfh + Bfh;
    end if;
//feed equations
    if molarBasis == true then
      Ffl = ffl / 3.6 "in mol/s";
      Ffh = ffh / 3.6;
    else
      Ffl = ffl / (mwL * 3600);
      Ffh = ffh / (mwH * 3600);
    end if;
    FZl = Ffl/F;
    FZh=Ffh/F;
    
    if isBinaryDistillation == true then
      F = Ffl+Ffh;
      Fmw = FZl * mwL + FZh * mwH;
      Dmw = DYl * mwL + DYh * mwH;
      Bmw = BXl * mwL + (1 - BXl) * mwH;
    end if;
//molar-mass transformation
    Fmass = F * Fmw;
    FflMass = Ffl * mwL;
    FZlMass=FflMass/Fmass;
    FfhMass = Ffh * mwH;
    FZhMass=FfhMass/Fmass;
    Dmass = D * Dmw;
    DflMass = Dfl * mwL;
    DYlMass=DflMass/Dmass;
    DfhMass=Dfh*mwH;
    DYhMass=DfhMass/Dmass;
    Bmass = B * Bmw;   
    BfhMass = Bfh * mwH;
    BXhMass=BfhMass/Bmass;
    BflMass = Bfl * mwL;
    BXlMass=BflMass/Bmass;
  
    if fSuperHeating < 0 then
      Fgf = 0;
      Q = 1 - fSuperHeating * cpFl / hvF;
    elseif fSuperHeating > 0 then
      Fgf = 1;
      Q = 1 - fSuperHeating * cpFv / hvF;
    else
      Fgf = fGfSat;
      Q = 1 - fGfSat - 1e-7 "to avoid division by 0 when Q=1";
    end if;
    Ls = Le + Q * F;
//condenser equations
    if topFixed == FreeFluids.Types.PerformanceSpecification.FixedComposition then
      if molarBasis == true then
        DYl = dyl;
      else
        DYlMass = dyl;
      end if;
    elseif topFixed == FreeFluids.Types.PerformanceSpecification.FixedFlow then
      if molarBasis == true then
        Dfl = dfl / 3.6 "to mol/s";
      else
        Dfl = dfl / (mwL * 3600);
      end if;
    elseif topFixed == FreeFluids.Types.PerformanceSpecification.FixedRecovery then
      Dfl = dlr * Ffl;
    end if;
    Dfl = D * DYl;
    DYh = Dfh/D;  
    RR = Le / D;
    Ve = D * (RR + 1);
    if isPartialCondenser == true then
      Wc = Le * hvD;
      RRext = RR;
    else
      Wc = Ve * (hvD + cpDl * dSubCooling);
      RRext = Le / (D * (1 + cpDl * dSubCooling / hvD));
    end if;
    if fixCondenserPower == true then
      Wc = wc;
    end if;
//reflux ratio
    if fixRefluxRatio == true then
      RR = rr;
    end if;
//reboiler equations
    BR = Vs / B;
    Ls = B + Vs;
    Wb = Vs * hvB;
    if bottomFixed == FreeFluids.Types.PerformanceSpecification.FixedComposition then
      if molarBasis==true then 
        BXh = bxh;
      else
        BXhMass=bxh;
      end if;
    elseif bottomFixed == FreeFluids.Types.PerformanceSpecification.FixedFlow then
      if molarBasis==true then 
        Bfh = bfh/3.6;
      else
        Bfh=bfh/(mwH*3600);
      end if;
    elseif bottomFixed == FreeFluids.Types.PerformanceSpecification.FixedRecovery then
      Bfh = brh * Ffh;
    end if;
    Bfh = B * BXh;
    BXl = Bfl/B;
    if fixReboilerPower == true then
      Wb = wb;
    end if;
    annotation(
      defaultComponentName = "Column",
      Icon(graphics = {Rectangle(origin = {7, -10}, lineThickness = 0.5, extent = {{-47, 90}, {33, -70}}), Rectangle(fillColor = {0, 85, 255}, fillPattern = FillPattern.Backward, lineThickness = 0, extent = {{-40, 60}, {40, -60}}), Text(origin = {-1, 99}, extent = {{-59, 11}, {59, -11}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html><head></head><body>Performs the total, and one component, molar mass balance, plus the heat balance at top and bottom. It assumes that the molar vaporization heat is constant along the enrichment and the stripping sectors, but not necessarily equal between them.<div>Although the balance seems adequate only for binary mixtures, it can be used also with multicomponent ones, taking into account that the heavy component stands for all other components different from the light one.&nbsp;</div></body></html>"));
  end DistillationBalance;

  partial model DistillationBase
    extends DistillationBalance;
    parameter FreeFluids.Types.ReboilerType reboilerType = FreeFluids.Types.ReboilerType.Kettle annotation(
      Dialog(tab = "Operating conditions"));
    parameter Types.EquilibriumDataType dataType = Types.EquilibriumDataType.EquilibriumConstant annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real fkl = 0 "Feed equilibrium constant of light component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real fkh = 0 "Feed equilibrium constant of heavy component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real dkl = 0 "Distillate equilibrium constant of light component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real dkh = 0 "Distillate equilibrium constant of heavy component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real bkl = 0 "Bottom equilibrium constant of light component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real bkh = 0 "Bottom equilibrium constant of heavy component" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real frv = 0 "Feed relative volatility light/heavy" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real drv = 0 "Distillate relative volatility light/heavy" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real brv = 0 "Bottom relative volatility light/heavy" annotation(
      Dialog(tab = "Equilibrium data"));
    Real Rmin(min = 0, start = 1) "minimum reflux ratio";
    Fraction RYl(start = 0.5) "operation lines intersection Y value";
    Fraction RXl(start = 0.5) "operation lines intersection X value";
    Real Frv(min = 0) "relative volatility of light/heavy at feed conditions";
    Real Drv(min = 0) "top of column relative volatility l/h";
    Real Brv(min = 0) "bottom of column relative volatility l/h";
    Real RVavg "average relative volatility l/h";
    Real NstagesMin "number of states at total reflux";
    Real Nstages "number of stages at working reflux";
    Real Nplates "number of plates needed at working reflux";
    Real Fstage(min = 0) "optimum stage for fee, from top";
  equation
//equilibrium calculation
    if dataType == Types.EquilibriumDataType.EquilibriumConstant then
      Frv = fkl / fkh;
      Drv = dkl / dkh;
      Brv = bkl / bkh;
    elseif dataType == Types.EquilibriumDataType.RelativeVolatility then
      Frv = frv;
      Drv = drv;
      Brv = brv;
    end if;
    RVavg = (Drv * Brv) ^ 0.5 "average relative volatility";
//reflux ratio
    if fixRatioToMinReflux == true then
      RR = Rmin * minRratio;
    end if;
    if noEvent(Le == Ls) then
      RYl = FZl;
    elseif noEvent(Ve == Vs) then
      RXl = FZl;
    else
      RYl = (RXl * Q - FZl) / (Q - 1) "feed line equation";
    end if;
    RR / (RR + 1) = (DYl - RYl) / (DYl - RXl);
//Number of stages
    if topFixed == FreeFluids.Types.PerformanceSpecification.FixedComposition then
      NstagesMin = log(DYl * BXh / (DYh * BXl)) / log(RVavg) "Fenske method using concentrations";
    elseif topFixed <> FreeFluids.Types.PerformanceSpecification.FixedComposition then
      NstagesMin = log(Dfl / Dfh * (Bfh / Bfl)) / log(RVavg) "Fenske method using flows";
    end if;
    (Nstages - NstagesMin) / (Nstages + 1) = 1 - exp((1 + 54.4 * (RR - Rmin) / (RR + 1)) / (11 + 117.2 * (RR - Rmin) / (RR + 1)) * ((RR - Rmin) / (RR + 1) - 1) / ((RR - Rmin) / (RR + 1)) ^ 0.5) "Gilliland correlation";
//optimum feed stage
    Fstage / (Nstages - Fstage) = abs(B / D * (FZh / FZl) * (BXl / DYh) ^ 2) ^ 0.206 "Kirkbride equation";
  algorithm
    Nplates := Nstages;
    if reboilerType <> FreeFluids.Types.ReboilerType.Thermosyphon then
      Nplates := Nplates - 1;
    end if;
    if isPartialCondenser == true then
      Nplates := Nplates - 1;
    end if;
    annotation(
      Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})),
      Documentation(info = "<html><head></head><body>The model is partial.</body></html>"));
  end DistillationBase;

  model BinaryDistillation
    extends DistillationBase(final isBinaryDistillation = true);
    parameter Real a0 = 1 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real a1 = 0 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real a2 = 0 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real a3 = 0 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real a4 = 0 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real a5 = 0 "polynom.coeff. for relative volatility calc." annotation(
      Dialog(tab = "Equilibrium data"));
    Fraction FXl "liquid feed molar fraction of light component";
    Fraction FYl "gas feed molar fraction of light component";
    Real RminRV "relative volatility at minimum reflux ratio";
    Fraction RminYl(start = 0.5) "minimum reflux BXlgas light concentration";
    Fraction RminXl(start = 0.5) "minimum reflux liquid light concentration";
    Fraction DXl "light concentrations at top of column reflux";
    Integer Ne "number of enrichment stages needed";
    Integer Ns "number of stripping stages needed";
    Real Xt, Yt, RVt "variables for the iteration";
    Real Dm "mean slope of the equilibrium line at enrichment section";
    Real Dlambda "mean stripping factor at enrichment";
    Real Bm "mean slope of equilibrium line at stripping section";
    Real Blambda "mean stripping factor at stripping section";
  equation
    if dataType == Types.EquilibriumDataType.RVpolynomia then
      Frv = a0 + a1 * FXl + a2 * FXl ^ 2 + a3 * FXl ^ 3 + a4 * FXl ^ 4 + a5 * FXl ^ 5;
      Drv = a0 + a1 * DXl + a2 * DXl ^ 2 + a3 * DXl ^ 3 + a4 * DXl ^ 4 + a5 * DXl ^ 5;
      Brv = a0 + a1 * BXl + a2 * BXl ^ 2 + a3 * BXl ^ 3 + a4 * BXl ^ 4 + a5 * BXl ^ 5;
      RminRV = a0 + a1 * RminXl + a2 * RminXl ^ 2 + a3 * RminXl ^ 3 + a4 * RminXl ^ 4 + a5 * RminXl ^ 5;
    else
      if Frv > 0 then
        RminRV = Frv;
      else
        RminRV = RVavg;
      end if;
    end if;
  //feed equations
    FYl = Frv * FXl / (1 + (Frv - 1) * FXl);
    FZl = Fgf * FYl + (1 - Fgf) * FXl;
  //condenser equations
    DYl = Drv * DXl / (1 + (Drv - 1) * DXl);
  //minimum reflux
    if noEvent(Le == Ls) then
      RminYl = FZl;
    elseif noEvent(Ve == Vs) then
      RminXl = FZl;
    else
      RminYl = (RminXl * Q - FZl) / (Q - 1) "feed line equation";
    end if;
    RminYl = RminRV * RminXl / (1 + (RminRV - 1) * RminXl) "equilibrium equation";
    Rmin / (Rmin + 1) = (DYl - RminYl) / (DYl - RminXl);
  //stripping factor
    Dm = (Drv * Frv) ^ 0.5 / (1 + (DXl + FZl) / 2 * ((Drv * Frv) ^ 0.5 - 1)) ^ 2;
    Dlambda=Dm*Ve/Le;
    Bm=Brv/(1+BXl*(Brv-1))^2;
    Blambda=Bm*Vs/Ls;
  //Number of stages
  algorithm
    Ne := 0;
    Yt := RYl;
    while noEvent(Yt < DYl) loop
      Xt := (Yt * (1 + RR) - DYl) / RR;
      if dataType == Types.EquilibriumDataType.RVpolynomia then
        RVt := a0 + a1 * Xt + a2 * Xt ^ 2 + a3 * Xt ^ 3 + a4 * Xt ^ 4 + a5 * Xt ^ 5;
      else
        RVt := RVavg;
      end if;
      Yt := RVt * Xt / (1 + (RVt - 1) * Xt);
      Ne := Ne + 1;
    end while;
    Ns := 0;
    Yt := Brv * BXl / (1 + (Brv - 1) * BXl);
    while noEvent(Yt < RYl) loop
      Xt := (Yt * BR + BXl) / (BR + 1);
      if dataType == Types.EquilibriumDataType.RVpolynomia then
        RVt := a0 + a1 * Xt + a2 * Xt ^ 2 + a3 * Xt ^ 3 + a4 * Xt ^ 4 + a5 * Xt ^ 5;
      else
        RVt := RVavg;
      end if;
      Yt := RVt * Xt / (1 + (RVt - 1) * Xt);
      Ns := Ns + 1;
    end while;
  
    annotation(
      defaultComponentName = "Column",
      Documentation(info = "<html><head></head><body>Allows the calculation of the relative volatility using an approximation polynomia.<div>The possibililty of a tangent point to the equilibrium line is not considered in the minimum reflux calculation.<br><div><br><div><br></div></div></div></body></html>"));
  end BinaryDistillation;

  model MulticomponentDistillationFUG
    extends DistillationBase(final isBinaryDistillation = false, D(start = ffl), B(start = ffh));
    parameter Modelica.Units.SI.MolarMass mwI[nNonKey] = fill(0, nNonKey) "nonkey components molecular weight";
    parameter String fLight = "" "key light component name" annotation(
      Dialog(tab = "Feed data"));
    parameter String fHeavy = "" "key heavy component name" annotation(
      Dialog(tab = "Feed data"));
    parameter String fNonKeys = "" "feed nonkey components" annotation(
      Dialog(tab = "Feed data"));
    parameter Integer nNonKey = 1 "number of nonkey components" annotation(
      Dialog(tab = "Feed data"));
    parameter Real ffi[nNonKey] = fill(0, nNonKey) "fixed feed flow of nonkey components in kmol/h or kg/h " annotation(
      Dialog(tab = "Feed data"));
    parameter Real fki[nNonKey] = fill(0, nNonKey) "Feed equilibrium constant of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real dki[nNonKey] = fill(0, nNonKey) "Distillate equilibrium constant of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real bki[nNonKey] = fill(0, nNonKey) "Bottom equilibrium constant of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real frvi_h[nNonKey] = fill(0, nNonKey) "feed relative volatility of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real drvi_h[nNonKey] = fill(0, nNonKey) "distillate relative volatility of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    parameter Real brvi_h[nNonKey] = fill(0, nNonKey) "bottom relative volatility of nonkey components" annotation(
      Dialog(tab = "Equilibrium data"));
    //Modelica.Units.SI.MassFlowRate DfhMass(displayUnit = "kg/h") "distillate heavy mass flow rate";
    Fraction FZi[nNonKey] "molar fraction of nonkey components in feed";
    Modelica.Units.SI.MolarFlowRate Ffi[nNonKey](each min = 0, each displayUnit = "kmol/h") "molar fraction of nonkey components in feed";
    Fraction DYi[nNonKey] "molar fraction of nonkey components in distillate";
    Modelica.Units.SI.MolarFlowRate Dfi[nNonKey](each min = 0, each displayUnit = "kmol/h") "molar fraction of nonkey components in distillate";
    Fraction BXi[nNonKey] "molar fraction of nonkey components in bottom";
    Modelica.Units.SI.MolarFlowRate Bfi[nNonKey](each min = 0, each displayUnit = "kmol/h") "molar fraction of nonkey components in distillate";
    Real Frvi_h[nNonKey](min = fill(0, nNonKey)) "top of column relative volatility nonkey/heavy";
    Real Drvi_h[nNonKey](min = fill(0, nNonKey)) "top of column relative volatility nonkey/heavy";
    Real Brvi_h[nNonKey](min = fill(0, nNonKey)) "bottom of column relative volatility nonkey/heavy";
    Real RVavgi_h[nNonKey] "average relative volatility nonkey/heavy";
    Real Sigma1;
    Real Sigma1i[nNonKey];
    Real Sigma2;
    Real Sigma2i[nNonKey];
    Real Theta(min = 1.0, start = 1.001) "Underwood parameter";
    Real Rmin1 "minimum reflux for type 1 pinch point";
  equation
  
    assert(FZl + FZh + sum(FZi[i] for i in 1:nNonKey) > 0.999 and FZl + FZh + sum(FZi[i] for i in 1:nNonKey) < 1.001, "Feed fractions summatory is not 1.0");
//feed balance
    if molarBasis == true then
      for i in 1:nNonKey loop
        Ffi[i] = ffi[i] / 3.6;
      end for;
    else
      for i in 1:nNonKey loop
        Ffi[i] = ffi[i] / (mwI[i] * 3600);
      end for;
    end if;
    for i in 1:nNonKey loop
      Ffi[i] = F * FZi[i];
    end for;
    F=Ffl+Ffh+sum(Ffi[i] for i in 1:nNonKey);
    Fmw = FZl * mwL + FZh * mwH + sum(FZi[i] * mwI[i] for i in 1:nNonKey);
    Dmw = DYl * mwL + DYh * mwH + sum(DYi[i] * mwI[i] for i in 1:nNonKey);
    Bmw = BXl * mwL + BXh * mwH + sum(BXi[i] * mwI[i] for i in 1:nNonKey);
  
    1 = BXl + BXh + sum(BXi[i] for i in 1:nNonKey);
//nonkey components balances. Las concentraciones salen a partir BXl y BXh, pero faltan ellas.
    for i in 1:nNonKey loop
      F * FZi[i] = D * DYi[i] + B * BXi[i];
      if dataType == Types.EquilibriumDataType.EquilibriumConstant then
        if fki[i] > 0 then
          Frvi_h[i] = fki[i] / fkh;
        else
          Frvi_h[i] = 0.0;
        end if;
        Drvi_h[i] = dki[i] / dkh;
        Brvi_h[i] = bki[i] / bkh;
      elseif dataType == Types.EquilibriumDataType.RelativeVolatility then
        Frvi_h[i] = frvi_h[i];
        Drvi_h[i] = drvi_h[i];
        Brvi_h[i] = brvi_h[i];
      end if;
      RVavgi_h[i] = (Drvi_h[i] * Brvi_h[i]) ^ 0.5;
      DYi[i] / BXi[i] = DYh / BXh * RVavgi_h[i] ^ NstagesMin;
      Dfi[i] = D * DYi[i];
      Bfi[i] = B * BXi[i];
    end for;
//minimum reflux
    1 - Q = Sigma1;
    for i in 1:nNonKey loop
      if Frvi_h[i] > 0 then
        Sigma1i[i] = Frvi_h[i] * FZi[i] / (Frvi_h[i] - Theta);
        Sigma2i[i] = Frvi_h[i] * DYi[i] / (Frvi_h[i] - Theta);
      else
        Sigma1i[i] = RVavgi_h[i] * FZi[i] / (RVavgi_h[i] - Theta);
        Sigma2i[i] = RVavgi_h[i] * DYi[i] / (RVavgi_h[i] - Theta);
      end if;
    end for;
    if Frv > 0 then
      Sigma1 = Frv * FZl / (Frv - Theta) + FZh / (1 - Theta) + sum(Sigma1i[i] for i in 1:nNonKey);
      Sigma2 = Frv * DYl / (Frv - Theta) + DYh / (1 - Theta) + sum(Sigma2i[i] for i in 1:nNonKey);
      Rmin = Sigma2 - 1;
      Rmin1 = (DYl / FZl - Frv * DYh / FZh) / (Frv - 1);
    else
      Sigma1 = RVavg * FZl / (RVavg - Theta) + FZh / (1 - Theta) + sum(Sigma1i[i] for i in 1:nNonKey);
      Sigma2 = RVavg * DYl / (RVavg - Theta) + DYh / (1 - Theta) + sum(Sigma2i[i] for i in 1:nNonKey);
      Rmin = Sigma2 - 1;
      Rmin1 = (DYl / FZl - RVavg * DYh / FZh) / (RVavg - 1);
    end if;
    annotation(
      defaultComponentName = "Column",
      Documentation(info = "<html><head></head><body><div>It is a non standard implementation of the Fenske-Underwood-Gilliland shortcut method for multicomponent distillation.</div><div>The following simplifications are made:</div><div>- The distribution of components between the distillate and bottom, obtained by the Fenske equation for total reflux, is considered valid for the actual and minimum reflux.</div><div>- It is considered just one pinch point, close to the feed composition. As the distillate composition is taken from the Fenske's equation, only one Underwood set of equations is needed in order to calculate the minimum reflux. This may give poor calculation if there is an important distribution of any component between distillate and bottom.</div><div>The feed must be specified as flow, in kmol/h or kg/h. The top light, and the bottom heavy, can be specified as flow or concentration. Flow is recommended for best convergence.</div><div>The model is lacking one equation. You can configure the needed equation by selecting the required parameter. It is done at the graphical menu under the \"Operating conditions\" tab. You can select one of the options as true: fix the reflux to a constant value, or as a multiplier of the minimum reflux. Or fix the heat exchange at the condenser or the reboiler. Or you can add the extra equation (for example the value of the vpor flow at the stripping section) by hand. &nbsp;<br><div><br><div><br></div></div></div></body></html>"));
  end MulticomponentDistillationFUG;

  //***PACKED COLUMN MODELS***

  model PackedColumnBasic "24 variables,16 equations,8 lacking"
    replaceable package MediumL = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    replaceable package MediumG = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    parameter Integer physProp = 0 "origin of the physical properties data. 0= equation, 1=parameters, 2=Medium" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.Temperature T = 298.15 "absolute temperature of the medium, if used" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.AbsolutePressure P = 1e5 "absolute pressure of the medium, if used" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.Density rhoL(displayUnit = "kg/m3") = 0 "liquid density to use if physProp==1" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.Density rhoG(displayUnit = "kg/m3") = 0 "gas density to use if physProp==1" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.DynamicViscosity muL = 0 "liquid viscosity to use if physProp==1" annotation(
      Dialog(tab = "Physical prop."));
    parameter Real epsilon = 0.965 "void fraction of the package" annotation(
      Dialog(tab = "Physical data"));
    parameter Real av=0 "specific surface of the packed m2/m3" annotation(
      Dialog(tab = "Physical data"));
    Real Av;
    Modelica.Units.SI.Diameter Di "column internal diameter";
    Modelica.Units.SI.Area S "column section";
    Modelica.Units.SI.MassFlowRate Gl(start = 1, displayUnit = "kg/h") "Liquid massic flow rate";
    Modelica.Units.SI.MassFlowRate Gg(start = 0.1, displayUnit = "kg/h") "Gas massic flow rate";
    Modelica.Units.SI.VolumeFlowRate Ql(displayUnit = "m3/h") "Liquid volumetric flow rate";
    Modelica.Units.SI.VolumeFlowRate Qg(displayUnit = "m3/h") "Gas volumetric flow rate";
    Modelica.Units.SI.Velocity Vl "actual liquid velocity, referenced to the empty column";
    Modelica.Units.SI.Velocity Vg "actual gas velocity, referenced to the empty column";
    Modelica.Units.SI.Density RhoL(start = 1000, displayUnit = "kg/m3");
    Modelica.Units.SI.Density RhoG(start = 1, displayUnit = "kg/m3");
    Real Fg "Gas Load factor at operating conditions";
    Real FgFl1 "Gas Load factor at loading conditions";
    Modelica.Units.SI.DynamicViscosity MuL(start = 1e-3) "Dynamic viscosity of liquid phase";
    Modelica.Units.SI.ReynoldsNumber ReL "Reynolds number of liquid phase at operating conditions";
    Modelica.Units.SI.FroudeNumber FrL "Froude number of liquid at operating conditions";
    Modelica.Units.SI.MassFlowRate GlFl1(displayUnit = "kg/h") "Liquid massic flow rate at flooding";
    Modelica.Units.SI.MassFlowRate GgFl1(displayUnit = "kg/h") "Gas massic flow rate at flooding";
    Modelica.Units.SI.VolumeFlowRate QlFl1(displayUnit = "m3/h") "Liquid volumetric flow rate at flooding";
    Modelica.Units.SI.VolumeFlowRate QgFl1(displayUnit = "m3/h") "Gas volumetric flow rate at flooding";
    Modelica.Units.SI.Velocity VlFl1(start = 0.004) "Liquid flooding velocity";
    Modelica.Units.SI.Velocity VgFl1(start = 10) "Gas flooding velocity";
    Modelica.Units.SI.ReynoldsNumber ReLfl1 "Reynolds number of liquid phase at flooding conditions";
    Fraction Load(start = 0.5) "fraction of maximum gas Load, at operating conditions";
    Real DiffP_h "pressure drop per meter of packing in (Pa/m)";
    MediumL.ThermodynamicState StateL;
    MediumG.ThermodynamicState StateG;
  equation
    StateL = MediumL.setState_pTX(P, T, fill(0, 0));
    StateG = MediumG.setState_pTX(P, T, fill(0, 0));
    if physProp == 1 then
      RhoL = rhoL;
      RhoG = rhoG;
      MuL = muL;
    elseif physProp == 2 then
      RhoL = MediumL.density(StateL);
      RhoG = MediumG.density(StateG);
      MuL = MediumL.dynamicViscosity(StateL);
    end if;
    S = pi * Di ^ 2 / 4;
    Gl = Ql * RhoL;
    Gg = Qg * RhoG;
    Vl = Ql / S;
    Vg = Qg / S;
    Fg = Vg * RhoG ^ 0.5;
    ReL = Vl * RhoL / Av / MuL;
    FrL = Vl ^ 2 * Av / g_n;
    GlFl1 = QlFl1 * RhoL;
    GgFl1 = QgFl1 * RhoG;
    VlFl1 = QlFl1 / S;
    VgFl1 = QgFl1 / S;
    GlFl1 = Gl / Gg * GgFl1;
    FgFl1 = VgFl1 * RhoG ^ 0.5;
    ReLfl1 = VlFl1 * RhoL / Av / MuL;
    Load = Gg / GgFl1;
    annotation(
      defaultComponentName = "Column",
      Icon(graphics = {Rectangle(origin = {7, -10}, lineThickness = 0.5, extent = {{-47, 90}, {33, -70}}), Rectangle(fillColor = {0, 85, 255}, fillPattern = FillPattern.Backward, lineThickness = 0, extent = {{-40, 60}, {40, -60}}), Text(origin = {-1, 99}, extent = {{-59, 11}, {59, -11}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
  end PackedColumnBasic;

  model PCSherwoodLobo "Flood point calculation by a numerical aproximation made by N.H.Chen to the Sherwood-Lobo curve. Not too exact"
    extends PackedColumnBasic;
    parameter Real fPk = 0 "packing factor to use,if equal to 0 will be calculated as Av/epsilon^3" annotation(
      Dialog(tab = "Physical data"));
    Real Flv "flow parameter";
    Real X, Y;
    Real Fp "packing factor";
  algorithm
    DiffP_h := 0 "Sherwood-Lobo model has no pressure drop calculation";
  equation
    Av=av;
    Flv = Gl / Gg * (RhoG / RhoL) ^ 0.5;
    X = log(100 * Flv) / log(10);
    Y = 10 ^ ((-0.1922 * X ^ 2) - 0.2041 * X - 0.5584);
    if fPk > 0 then
      Fp = fPk;
    else
      Fp = Av / epsilon ^ 3;
    end if;
    VgFl1 = (Y / Fp / (MuL * 1000) ^ 0.2 * RhoL * RhoG * g_n) ^ 0.5 / RhoG;
  end PCSherwoodLobo;

  model PCLeva "Leva model has no flooding point calculation"
    //Pressure drop calculation according to Leva (modified version of generalized pressure drop correlation)
    //Pressure drop at flooding with water has been added, according to Kister and Gill
    extends PackedColumnBasic;
    parameter Real fPk = 0 "packing factor to use,if equal to 0 will be calculated as Av/epsilon^3" annotation(
      Dialog(tab = "Physical data"));
    Real Fp "packing factor used";
    Real Flv "flow parameter";
    Real DiffP_hFl1w "pressure drop at flooding at water viscosity";
  algorithm
    VgFl1 := 1e+099;
    if fPk > 0 then
      Fp := fPk;
    else
      Fp := Av / epsilon ^ 3;
    end if;
    DiffP_hFl1w := 12.718 * Fp ^ 0.7 "0.115*25.4*10*(Fp/3.28)^0.7";
  equation
    Av=av;
    Flv = Gl / Gg * (RhoG / RhoL) ^ 0.5;
    DiffP_h = 22.3 * Fp * MuL ^ 0.2 * 1000 / RhoL * (Vg * RhoG) ^ 2 * 10 ^ (0.035 * Vl * RhoL * 1000 / RhoL) / g_n / RhoG;
  end PCLeva;

  partial model PackedColumnExtended "24+20=44 variables,16+8=24 equations, 20 lacking"
    extends PackedColumnBasic;
    parameter Modelica.Units.SI.DynamicViscosity muG = 0 "Dynamic viscosity of gas phase if physProp==1" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.SurfaceTension sigmaL = 0 "Surface tension of the liquid if physProp==1" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.DiffusionCoefficient dL = 0 "diffusion coefficients in liquid phase if physProp==1 or 2" annotation(
      Dialog(tab = "Physical prop."));
    parameter Modelica.Units.SI.DiffusionCoefficient dG = 0 "diffusion coefficients in gas phase if physProp==1 or 2" annotation(
      Dialog(tab = "Physical prop."));
    parameter Real lambda=1 "stripping factor" annotation(
      Dialog(tab = "Operating conditions"));
    Modelica.Units.SI.DynamicViscosity MuG "Dynamic viscosity of gas phase";
    Modelica.Units.SI.SurfaceTension SigmaL "Surface tension of the liquid";
    Modelica.Units.SI.DiffusionCoefficient Dl "diffusion coefficients in liquid phase";
    Modelica.Units.SI.DiffusionCoefficient Dg "diffusion coefficients in gas phase";
    Modelica.Units.SI.ReynoldsNumber ReG "Reynolds number of gas phase at operating conditions";
    Modelica.Units.SI.WeberNumber WeL "Weber number of liquid phase at operating conditions";
    Modelica.Units.SI.SchmidtNumber ScL "Schmidt number of liquid phase at operating conditions";
    Modelica.Units.SI.SchmidtNumber ScG "Schmidt number of gas phase at operating conditions";
    Modelica.Units.SI.Velocity VlS1(start = 0.01) "Liquid velocity at loading point (fixed L/G ratio)";
    Modelica.Units.SI.Velocity VgS1(start = 2) "Gas velocity at loading point (fixed L/G ratio)";
    Modelica.Units.SI.Velocity VgS2(start = 2) "Gas velocity at loading point (fixed Vl)";
    Modelica.Units.SI.Velocity VgFl2(start = 5) "Gas velocity at flooding point (fixed Vl)";
    Real FgS1 "Gas loading factor at loading point(fixed L/G ratio)";
    Real FgS2 "Gas loading factor at loading point(fixed Vl)";
    Real FgFl2 "Gas loading factor at flooding point(fixed Vl)";
    Real Aph "Specific interface area between phases (m2/m3)";
    Real BetaL "Mass transfer coeff. in liquid phase m/s";
    Real BetaG "Mass transfer coeff. in gas phase m/s";
    Modelica.Units.SI.Height HTUl "Height of the transfer unit for liquid phase";
    Modelica.Units.SI.Height HTUg "Height of the transfer unit for gas phase";
    Modelica.Units.SI.Height HTUov "Overall height of the transfer unit";
    Modelica.Units.SI.Height HETP "Height equivalent to a theoretical plate";
  equation
    if physProp == 1 then
      MuG = muG;
      SigmaL = sigmaL;
      Dl = dL;
      Dg = dG;
    elseif physProp == 2 then
      MuG = MediumG.dynamicViscosity(StateG);
      SigmaL = MediumL.surfaceTension(MediumL.setSat_T(T));
      Dl = dL;
      Dg = dG;
    end if;
    ReG = Vg * RhoG / Av / MuG;
    WeL = Vl ^ 2 * RhoL / Av / SigmaL;
    ScL = MuL / RhoL / Dl;
    ScG = MuG / RhoG / Dg;
    VlS1 = VgS1 * Vl / Vg;
    FgS1 = VgS1 * RhoG ^ 0.5;
    FgS2 = VgS2 * RhoG ^ 0.5;
    FgFl2 = VgFl2 * RhoG ^ 0.5;
    HTUov = HTUg + HTUl * lambda;
    if lambda<1 then
      HETP = HTUov * log(lambda)/(lambda - 1);
    else
      HETP=HTUov;
    end if;
    
  end PackedColumnExtended;

  partial model PCOndaBravo "Onda or Bravo models has no hydraulic calculations"
    extends PackedColumnExtended;
    parameter String material = "metal" "alternatives are: metal, ceramic, plastic" annotation(
      Dialog(tab = "Physical data"));
    parameter Real size = 1.0 "packing size en inches" annotation(
      Dialog(tab = "Physical data"));
    Modelica.Units.SI.SurfaceTension SigmaP "surface tension of the packing";
    Real A;
  algorithm
    VgS1 := 1e+099;
    VgFl1 := 1e+099;
    VgS2 := 1e+099;
    VgFl2 := 1e+099;
    DiffP_h := 0;
    if material == "metal" then
      SigmaP := 0.075;
    elseif material == "ceramic" then
      SigmaP := 0.061;
    elseif material == "plastic" then
      SigmaP := 0.003;
    end if;
    if size <= 0.5 then
      A := 2.0;
    else
      A := 5.23;
    end if;
  equation
    Av=av;
    BetaL = 0.0051 * (ReL * Av / Aph) ^ (2 / 3) * ScL ^ (-0.5) * (Av * size * 0.0254) ^ 0.4 * (MuL * g_n / RhoL) ^ (1 / 3);
    BetaG = A * ReG ^ 0.7 * ScG ^ (1 / 3) * (Av * size * 0.0254) ^ (-2) * Av * Dg;
    HTUl = Vl / BetaL / Aph;
    HTUg = Vg / BetaG / Aph;
  end PCOndaBravo;

  model PCOnda
    extends PCOndaBravo;
  equation
    Aph = Av * (1 - exp(-1.45 * ReL ^ 0.1 * FrL ^ (-0.05) * WeL ^ 0.2 * (SigmaP / SigmaL) ^ 0.75));
  end PCOnda;

  model PCBravo
    extends PCOndaBravo;
    Real CaL "capillary number of the liquid";
    Modelica.Units.SI.Height H "packed height";
  equation
    CaL = Vl * MuL / SigmaL;
    Aph = Av * 19.78 * (CaL * 6 * ReG) ^ 0.392 * SigmaL ^ 0.5 / H ^ 0.4;
  end PCBravo;

  //model for unstructured and structured packing

  model PCBilletSchultes   "package characteristic parameters. Values are for 1,5 metallic Pall ring, replace them as needed"
    extends PackedColumnExtended(physProp=1);
    parameter Real cH = 0.644 "parameter for liquid holdup" annotation(
      Dialog(tab = "Physical data"));
    parameter Real cSt = 2.629 "parameter for loading point calculation" annotation(
      Dialog(tab = "Physical data"));
    parameter Real cFlt = 1.679 "parameter for flooding point calculation" annotation(
      Dialog(tab = "Physical data"));
    parameter Real cP0 = 0.967 "pressure drop coefficient" annotation(
      Dialog(tab = "Physical data"));
    parameter Real cL = 1.012 "parameter for liquid mass transfer" annotation(
      Dialog(tab = "Physical data"));
    parameter Real cG = 0.341 "parameter for gas mass transfer" annotation(
      Dialog(tab = "Physical data"));
    //liquid holdup, and loading and flooding velocities, related variables
    Real HlS "liquid holdup (m3/m3) at loading point";
    Real HlFl1(start = 0.4) "liquid holdup (m3/m3) at flooding point";
    Real DiffP_hFl1 "prerssure loss (Pa/m) at flooding point";
    //Modelica.Units.SI.Velocity vLs(start=0.01),vGs;//liquid and gas loading point velocities (initiation of liquid holdup)
    Real PsiS "resistance coefficient at loading point";
    Real PsiFl1 "resistance coefficient at flooding point";
    Real Ns, NFl1 "exponent for phiS and phiFl calculation";
    Real Cs, CFl1 "corrected parameters for loading and flooding points calculation";
    Real Hl "liquid holdup at operating conditions";
    //Pressure drop related variables
    Modelica.Units.SI.Diameter Dp "diameter of the particle of falling liquid";
    Real K "wall factor";
    Modelica.Units.SI.ReynoldsNumber ReGmod "modified gas Reynolds number at operating conditions";
    Real PsiL "resistance coefficient at operating conditions";
    //HTU related variables
    Modelica.Units.SI.Diameter Dh "hydraulic diameter of the package";
    Real AphS "specific interface area between phases in loading condition";
    Real AphFl1 "specific interface area between phases in flooding condition";
    Modelica.Units.SI.Velocity VmeanL "mean effective velocity of the liquid";
    //Real liquid holdup related variables
    Real Ah "specific surface of the liquid m2/m3 at operating conditions";
    Real AhS "specific surface of the liquid m2/m3 at loading conditions";
    Real AhFl1 "specific surface of the liquid m2/m3 at flooding conditions";
    Real HlSR "Real liquid holdup at loading point";
    Real HlFl1R "Real liquid holdup at flooding point";
    Real HlR "Real liquid holdup at operating point";
    Modelica.Units.SI.ReynoldsNumber ReLs "Reynolds number of the liquid at loading point";
    Modelica.Units.SI.FroudeNumber FrLs "Froude number of the liquid at loading point";
    Modelica.Units.SI.FroudeNumber FrLfl1 "Froude number of the liquid at flooding point";
  equation
    Av=av;
    VgS2 = 1e+99 "Billet and Schultes model does not need loading and flooding calculations at constant Vl, but could be done";
    VgFl2 = 1e+99;
//loading and flooding velocities, and liquid holdup, calculations
    if Gl / Gg * (RhoG / RhoL) ^ 0.499 <= 0.4 then
//if (Gl / Gg * (RhoG / RhoL) ^ 0.5) <= 0.4 then
      Ns = -0.326;
      Cs = cSt;
      NFl1 = -0.194;
      CFl1 = cFlt;
    else
      Ns = -0.723;
      Cs = 0.695 * cSt * (MuL / MuG) ^ 0.1588;
      NFl1 = -0.708;
      CFl1 = 0.6244 * cFlt * (MuL / MuG) ^ 0.1028;
    end if;
    PsiS = g_n / Cs ^ 2 * (Gl / Gg * (RhoG / RhoL) ^ 0.5 * (MuL / MuG) ^ 0.4) ^ (-2 * Ns);
    PsiFl1 = g_n / CFl1 ^ 2 * (Gl / Gg * (RhoG / RhoL) ^ 0.5 * (MuL / MuG) ^ 0.2) ^ (-2 * NFl1);
    VgS1 = (g_n / PsiS) ^ 0.5 * (epsilon / Av ^ (1 / 6) - Av ^ 0.5 * (12 * MuL * VlS1 / g_n / RhoL) ^ (1 / 3)) * (12 * MuL * VlS1 / g_n / RhoL) ^ (1 / 6) * (RhoL / RhoG) ^ 0.5;
//vLs=(RhoG/RhoL)*(Gl/Gg)*vGs;
    HlS = (12 * MuL * VlS1 * Av ^ 2 / g_n / RhoL) ^ (1 / 3);
    VgFl1 = 2 ^ 0.5 * (g_n / PsiFl1) ^ 0.5 * (epsilon - HlFl1) ^ (3 / 2) / epsilon ^ 0.5 * (HlFl1 / Av) ^ 0.5 * (RhoL / RhoG) ^ 0.5;
    HlFl1 ^ 3 * (3 * HlFl1 - epsilon) = 6 * Av ^ 2 * epsilon * MuL * Gl * RhoG * VgFl1 / (g_n * RhoL * Gg * RhoL);
    DiffP_hFl1 = PsiFl1 * Av * (VgFl1 * RhoG ^ 0.5) ^ 2 / ((epsilon - HlFl1R) ^ 3 * 2 * K);
//liquid holdup at operating conditions
    if Vg <= VgS1 then
      Hl = (12 * MuL * Vl * Av ^ 2 / g_n / RhoL) ^ (1 / 3);
    else
      Hl = HlS + (HlFl1 - HlS) * (Vg / VgFl1) ^ 13;
    end if;
//Pressure drop calculation
    Dp = 6 * (1 - epsilon) / Av;
    K = 1 / (1 + 2 * Dp / (3 * (1 - epsilon) * Di));
    ReGmod = Vg * Dp * K * RhoG / (1 - epsilon) / MuG;
    PsiL = cP0 * (64 / ReGmod + 1.8 / ReGmod ^ 0.08) * ((epsilon - Hl) / epsilon) ^ 1.5 * (Hl / HlS) ^ 0.3 * exp(13300 / Av ^ (3 / 2) * FrL ^ 0.5);
    DiffP_h = PsiL * Av * (Vg * RhoG ^ 0.5) ^ 2 / ((epsilon - Hl) ^ 3 * 2 * K);
//HTU calculation
    Dh = 4 * epsilon / Av;
    AphS / Av = 1.5 * (Av * Dh) ^ (-0.5) * (VlS1 * Dh * RhoL / MuL) ^ (-0.2) * (VlS1 ^ 2 * RhoL * Dh / SigmaL) ^ 0.75 * (VlS1 ^ 2 / g_n / Dh) ^ (-0.45);
    AphFl1 / Av = 10.5 * (SigmaL / 0.07199999999999999) ^ 0.5600000000000001 * (Av * Dh) ^ (-0.5) * (Vl * Dh * RhoL / MuL) ^ (-0.2) * (Vl ^ 2 * RhoL * Dh / SigmaL) ^ 0.75 * (Vl ^ 2 / g_n / Dh) ^ (-0.45);
    if Vg <= VgS1 then
      VmeanL = Vl / Hl;
      Aph / Av = 1.5 * (Av * Dh) ^ (-0.5) * (Vl * Dh * RhoL / MuL) ^ (-0.2) * (Vl ^ 2 * RhoL * Dh / SigmaL) ^ 0.75 * (Vl ^ 2 / g_n / Dh) ^ (-0.45);
    else
      VmeanL = (g_n * RhoG ^ 2 * Vg ^ 2 / 12 / MuL / Av ^ 2 / RhoL) ^ (1 / 3) * (Gl / Gg) ^ (2 / 3) * (1 - ((Vg - VgS1) / (VgFl1 - VgS1)) ^ 2);
      Aph = AphS + (AphFl1 - AphS) * (Vg / VgFl1) ^ 13;
    end if;
//BetaL=cL*12^(1/6)*VmeanL^0.5*(Dl/Dh)^0.5;
    BetaL = cL * (g_n * RhoL / MuL) ^ (1 / 6) * (Dl / Dh) ^ 0.5 * (Vl / Av) ^ (1 / 3);
    BetaG = cG * Dg * (Av / (epsilon - Hl) / Dh) ^ 0.5 * ReG ^ (3 / 4) * ScG ^ (1 / 3);
    HTUl = Vl / BetaL / Aph;
    HTUg = Vg / BetaG / Aph;
//Real liquid holdup at loading and flooding points. Not the same as theoretical, that are used for all the calculations
    ReLs = VlS1 * RhoL / Av / MuL;
//ReLfl1 = VlFl1 * RhoL / Av / MuL;
    FrLs = Av * VlS1 ^ 2 / g_n;
    FrLfl1 = Av * VlFl1 ^ 2 / g_n;
    if ReLs < 5 then
      AhS / Av = cH * ReLs ^ 0.15 * FrLs ^ 0.1;
    else
      AhS / Av = cH * 0.85 * ReLs ^ 0.25 * FrLs ^ 0.1;
    end if;
    if ReLfl1 < 5 then
      AhFl1 / Av = cH * ReLfl1 ^ 0.15 * FrLfl1 ^ 0.1;
    else
      AhFl1 / Av = cH * 0.85 * ReLfl1 ^ 0.25 * FrLfl1 ^ 0.1;
    end if;
    if ReL < 5 then
      Ah / Av = cH * ReL ^ 0.15 * FrL ^ 0.1;
    else
      Ah / Av = cH * 0.85 * ReL ^ 0.25 * FrL ^ 0.1;
    end if;
    HlSR = (12 * MuL * VlS1 * Av ^ 2 / g_n / RhoL) ^ (1 / 3) * (AhS / Av) ^ (2 / 3);
    HlFl1R = 2.2 * (12 * MuL * VlFl1 * Av ^ 2 / g_n / RhoL) ^ (1 / 3) * (AhFl1 / Av) ^ (2 / 3) * (MuL * 1000 / 0.001 / RhoL) ^ 0.05;
    if Vg < VgS1 then
      HlR = (12 * MuL * Vl * Av ^ 2 / g_n / RhoL) ^ (1 / 3) * (Ah / Av) ^ (2 / 3);
    else
      HlR = HlSR + (HlFl1R - HlSR) * (Vg / VgFl1) ^ 13;
    end if;
//HlR=(12*MuL*Vl*Av^2/g_n/RhoL)^(1/3)*(Ah/Av)^(2/3)+(HlFl1R-(12*MuL*Vl*Av^2/g_n/RhoL)^(1/3)*(Ah/Av)^(2/3))*(Vg/VgFl1)^13;
  end PCBilletSchultes;

  model PCDelft "Delft model 2004 for structured packing. Not very good for pressure drop (too high). Bad for liquid holdup after loading point (it is constant)"
    extends PackedColumnExtended;
    //packing characterization
    parameter Modelica.Units.SI.Angle alpha "angle of packing channels with horizontal" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.Units.SI.Height hPb = 0 "height of the packed bed" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.Units.SI.Height hPe = 0 "height of the packing element" annotation(
      Dialog(tab = "Physical data"));
    parameter Real omega = 0 "percentage of packing surface occupied by holes" annotation(
      Dialog(tab = "Physical data"));
    Modelica.Units.SI.Length Sw "width of the packing sheets between coarrugations";
    Modelica.Units.SI.Length h "height of the packing coarrugations";
    Modelica.Units.SI.Length B "distance between edges of coarrugation";
    Modelica.Units.SI.Angle AlphaL "efective flow angle of the liquid";
    Real phi "fraction of the gas flow channel occupied by the liquid";
    Modelica.Units.SI.Distance deltaL "mean thickness of liquid film";
    Real Hl "liquid holdup at operating conditions";
    Modelica.Units.SI.Velocity Vle, Vge "effective velocities of liquid and gas";
    Modelica.Units.SI.Diameter dhG "hydraulic diameter of the triangular gas flow channel";
    Modelica.Units.SI.ReynoldsNumber ReGe(start = 100) "Gas Reynods numbers effective";
    Modelica.Units.SI.ReynoldsNumber ReGr(start = 1000) "Gas Reynods numbers relative";
    Real XiGL "Gas-liquid friction factor";
    Real zetaGL, zetaGG, zetaDC "Overall gas-liquid, gas-gas, and direction change, pressure losses coefficients";
    Real xiWall, xiBulk "direction change factors for wall and bulk zones";
    Real Psi "fraction of gas flow channels ending at column walls";
    Real Fload "Correction factor for pressure loss";
    Real betaGlam "gas phase mass transfer coeff. under laminar conditions";
    Real betaGturb "gas phase mass transfer coeff. under turbulent conditions";
  algorithm
    AlphaL := atan(cos(pi / 2 - alpha) / sin(pi / 2 - alpha) / cos(atan(B / 2 / h)));
//The Delft model does not calculate the flooding point at constant G/L, but at fixed L
    VgFl1 := 1e+099;
    VgS1 := 1e+099;
  equation
    Sw = (B ^ 2 / 4 + h ^ 2) ^ 0.5 "relationship between dimensions of packing";
    Av = 4 * Sw / B / h "relationship between dimensions of packing";
    phi = 2 * Sw / (B + 2 * Sw) "relationship between dimensions of packing";
    deltaL = (3 * MuL * Vl / RhoL / g_n / Av / sin(AlphaL)) ^ (1 / 3);
    Hl = Av * deltaL;
    Vge = Vg / (epsilon - Hl) / sin(alpha);
    Vle = Vl / epsilon / Hl / sin(AlphaL);
    dhG = (B * h - 2 * deltaL * Sw) ^ 2 / B / h / ((((B * h - 2 * deltaL * Sw) / 2 / h) ^ 2 + ((B * h - 2 * deltaL * Sw) / B) ^ 2) ^ 0.5 + (B * h - 2 * deltaL * Sw) / 2 / h);
    ReGe = dhG * Vge * RhoG / MuG;
    ReGr = dhG * (Vge + Vle) * RhoG / MuG;
    FgS2 = (0.053 * epsilon ^ 2 * g_n * dhG * (RhoL - RhoG) / RhoG * (Vl * (RhoL / RhoG) ^ 0.5) ^ (-0.25) * sin(alpha) ^ 1.24) ^ 0.57 * RhoG ^ 0.5;
    VgFl2 = FgS2 / RhoG ^ 0.5 / 0.65;
    XiGL = (-2 * log(deltaL / dhG / 3.7 - 5.02 / ReGr * log(deltaL / dhG / 3.7 + 14.5 / ReGr))) ^ (-2);
    zetaGL = phi * XiGL * hPb / dhG / sin(alpha);
    zetaGG = (1 - phi) * 0.722 * cos(alpha) ^ 3.14 * hPb / dhG / sin(alpha);
    xiWall = (4092 * Vl ^ 0.31 + 4715 * cos(alpha) ^ 0.445) / ReGe + 34.19 * Vl ^ 0.44 * cos(alpha) ^ 0.779;
    xiBulk = 1.76 * cos(alpha) ^ 1.63;
    if Di >= hPe / tan(alpha) then
      Psi = 2 * hPe / pi / Di ^ 2 / tan(alpha) * (Di ^ 2 * hPe ^ 2 / tan(alpha) ^ 2) ^ 0.5 + 2 / pi * asin(hPe / Di / tan(alpha));
//ojo 2/pi y ojo * en vez de -
    else
      Psi = 2 * hPe / pi / Di ^ 2 / tan(alpha) * (Di ^ 2 * hPe ^ 2 / tan(alpha) ^ 2) ^ 0.5 + 1;
    end if;
    zetaDC = hPb / hPe * (xiBulk + Psi * xiWall);
    Fload = 3.8 * (Fg / FgS2) ^ (2 / sin(alpha)) * (Vl ^ 2 / epsilon ^ 2 / g_n / dhG) ^ 0.13;
    DiffP_h = Fload * (zetaGL + zetaGG + zetaDC) * RhoG * Vge ^ 2 / 2;
    Aph = (1 - omega) * Av * (1 - exp(-1.45 * ReL ^ 0.1 * FrL ^ (-0.05) * WeL ^ 0.2 * (0.075 / SigmaL) ^ 0.75));
    BetaL = 2 * (Dl * Vle / 0.9 / dhG / pi) ^ 0.5;
    betaGlam = 0.664 * ScG ^ (1 / 3) * (ReGr * dhG * sin(alpha) / hPe) ^ 0.5 * Dg / dhG;
    betaGturb = ReGr * ScG * XiGL * phi / 8 / (1 + 1.27 * (XiGL * phi / 8) ^ 0.5 * (ScG ^ (2 / 3) - 1)) * (1 + (dhG * sin(alpha) / hPe) ^ (2 / 3)) * Dg / dhG;
    BetaG = (betaGlam ^ 2 + betaGturb ^ 2) ^ 0.5;
    HTUl = Vl / BetaL / Aph;
    HTUg = Vg / BetaG / Aph;
  end PCDelft;

  package Examples
  extends Modelica.Icons.ExamplesPackage;
    model BinaryiC4nC4 "Packed towers. Billet 1995 pag 344"
      extends Modelica.Icons.Example;
      BinaryDistillation Column( bottomFixed = FreeFluids.Types.PerformanceSpecification.FixedComposition, brv = 1.35, bxh = 0.999, dataType = FreeFluids.Types.EquilibriumDataType.RelativeVolatility, drv = 1.35, dyl = 0.999, fGfSat = 0, ffh = 200, ffl = 300, fixRefluxRatio = true, frv = 1.35, isPartialCondenser = true, mwH = 0.058122, mwL = 0.058122, rr = 15, topFixed = FreeFluids.Types.PerformanceSpecification.FixedComposition, wb = 0) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    annotation(
        Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})),
        Documentation(info = "<html><head></head><body>Packed towers, R.Billet 1995. Example 14.2 pag.344</body></html>"));
    end BinaryiC4nC4;
    
    model FUGSeaderChap9
      extends Modelica.Icons.Example;
      MulticomponentDistillationFUG Column(bfh = 23, bottomFixed = FreeFluids.Types.PerformanceSpecification.FixedFlow, brv = 1.44, brvi_h = {1.6, 0.819, 0.5, 0.278, 0.167, 0.108}, dataType = FreeFluids.Types.EquilibriumDataType.RelativeVolatility, dfl = 442, drv = 2.08, drvi_h = {2.81, 0.737, 0.303, 0.123, 0.045, 0.02}, fGfSat = 0.1334, fHeavy = "iC5", fLight = "nC4", fNonKeys = "iC4,nC5,nC6,nC7,nC8,nC9", ffh = 36, ffi = {12, 15, 23, 39.1, 272.2, 31}, ffl = 448, fixRatioToMinReflux = true, frv = 1.93, frvi_h = {2.43, 0.765, 0.362, 0.164, 0.072, 0.0362}, minRratio = 1.3, mwH = 0.072149, mwI = {0.058122, 0.072149, 0.086175, 0.1002, 0.11423, 0.12825}, mwL = 0.058122, nNonKey = 6, topFixed = FreeFluids.Types.PerformanceSpecification.FixedFlow)  annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      annotation(
        Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})),
        Documentation(info = "<html><head></head><body>Example from Separation Process Principles 3rd edition, Seader et alt.&nbsp;</body></html>"));
    end FUGSeaderChap9;

    model UsagePCSherwoodLobo "Aalto university example. Water-Methanol, top of column"
      FreeFluids.Columns.PCSherwoodLobo Column(av=139.4,Gg(displayUnit = "kg/s"), Gl(displayUnit = "kg/s"), redeclare package MediumG = FreeFluids.TMedia.Fluids.Methanol, redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, T = 343.15, epsilon = 0.965, fPk = 155.125, muL = 4.06e-4, physProp = 2, rhoG = 0.198, rhoL = 978) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      //Column.Av := 139.4;
//Column.Di := 0.156;
      Column.Load := 0.6;
      Column.Gl := 0.35 * 0.0292;
      Column.Gg := 0.23 * 0.0292;
    end UsagePCSherwoodLobo;

    model UsagePCLeva "Aalto university example"
      FreeFluids.Columns.PCLeva Column(av = 139.4, Gg(displayUnit = "kg/s"), Gl(displayUnit = "kg/s"), redeclare package MediumG = FreeFluids.TMedia.Fluids.Water, redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, T = 343.15, epsilon = 0.965, muL = 4.06e-4, physProp = 1, rhoG = 0.198, rhoL = 978) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm

    end UsagePCLeva;

    model UsagePCOnda "Aalto university example"
      FreeFluids.Columns.PCOnda Column(av=256,redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, redeclare package MediumG = FreeFluids.TMedia.Fluids.Water, T = 343.15, epsilon = 0.89, material = "metal", size = 1.0, physProp = 1, rhoL = 978, rhoG = 0.198, muL = 4.06e-4, muG = 1.1e-5, sigmaL = 1.86e-2, dL = 5e-9, dG = 2e-5) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      //extends PCOnda(epsilon = 0.89, material = "metal", size = 1.0);
    algorithm
      //Column.Av := 256;
      Column.Di := 0.156;
      Column.Gl := 0.35 * 0.0292;
      Column.Gg := 0.23 * 0.0292;
    end UsagePCOnda;

    model UsagePCOnda99
      //pag 99 example bielacki rings metal 1"
      PCOnda Column(av=238, Gg(displayUnit = "kg/s"), Gl(displayUnit = "kg/s"), redeclare package MediumG = FreeFluids.TMedia.Fluids.Water, redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, VgFl1(start = 2), dG = 2e-5, dL = 5e-9, epsilon = 0.942, material = "metal", muG = 1.52e-5 * 1.17, muL = 1e-3, physProp = 1, rhoG = 1.17, rhoL = 998.2, sigmaL = 7.24e-2, size = 1.0) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      //Column.Av := 238;
      Column.Di := 0.15;
      Column.Vl := 11.1e-3;
      Column.Vg := 1;
    end UsagePCOnda99;

    model UsagePCBravo "Aalto university example"
      PCBravo Column(av = 256, epsilon = 0.89, material = "metal", size = 1.0, physProp = 1, rhoL = 978, rhoG = 0.198, muL = 4.06e-4, muG = 1.1e-5, sigmaL = 1.86e-2, dL = 5e-9, dG = 2e-5) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      //extends PCBravo(av = 256, epsilon = 0.89, material = "metal", size = 1.0);
    algorithm
      Column.Di := 0.156;
      Column.Gl := 0.35 * 0.0292;
      Column.Gg := 0.23 * 0.0292;
      Column.H := 2;
    end UsagePCBravo;

    model UsagePCBravo99
      //pag 99 example bielacki rings metal 1"
      PCBravo col(av = 238, Gg(displayUnit = "kg/s"), Gl(displayUnit = "kg/s"), VgFl1(start = 2), dG = 2e-5, dL = 5e-9, epsilon = 0.942, material = "metal", muG = 1.52e-5 * 1.17, muL = 1e-3, physProp = 1, rhoG = 1.17, rhoL = 998.2, sigmaL = 7.24e-2, size = 1.0) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      col.Di := 0.15;
      col.Vl := 11.1e-3;
      col.Vg := 1;
      col.H := 2;
    end UsagePCBravo99;

    model UsagePCBilletSchultes "Aalto university example. Water-Methanol, top of column"
      PCBilletSchultes Column(av = 256, redeclare package MediumG = FreeFluids.TMedia.Fluids.Water, redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, cFlt = 2, cH = 0.75, cP0 = 0.891, cSt = 2.5, epsilon = 0.89, physProp = 1, rhoL = 978, rhoG = 0.198, muL = 4.06e-4, muG = 1.1e-5, sigmaL = 1.86e-2, dL = 5e-9, dG = 2e-5) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      //extends PCBilletSchultes(av = 256, epsilon = 0.89, cSt = 2.5, cFlt = 2, cH = 0.75, cP0 = 0.891);
    algorithm
      Column.Di := 0.156;
      Column.Gl := 0.35 * 0.0292;
      Column.Gg := 0.23 * 0.0292;
    end UsagePCBilletSchultes;

    model UsagePCBilletSchultes2 "Billet example 14.4. Acetone in Water-Air. Pall metal 50 mm"
      PCBilletSchultes Column(redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, redeclare package MediumG = Modelica.Media.Air.DryAirNasa, av = 110, T = 300.15, cFlt = 1.58, cG = 0.41, cH = 0.784, cL = 1.192, cP0 = 0.763, cSt = 2.725, dG = 10.8e-6, dL = 1.18e-9, epsilon = 0.952, muG = 18.13e-6, muL = 0.857e-3, physProp = 1, rhoG = 1.162, rhoL = 997, sigmaL = 72e-3) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      //extends PCBilletSchultes(av = 110, epsilon = 0.952, cSt = 2.725, cFlt = 1.58, cH = 0.784, cP0 = 0.763, cL = 1.192, cG = 0.41);
    algorithm
//Column.Di := 0.6;
      Column.Vg := Column.VgS1;
      Column.Gl := 4617.7 / 3600;
//Column.Gg := 2301 / 3600;
      Column.Qg := 2000 / 3600;
    end UsagePCBilletSchultes2;

    model UsagePCBilletSchultes3 "Billet example 14.5. CO2 in Water_Air at 20ºC and 1 bar. Pall plastic 35 mm"
      PCBilletSchultes Column(redeclare package MediumL = FreeFluids.TMedia.Fluids.Water, redeclare package MediumG = Modelica.Media.Air.DryAirNasa, av = 145, cFlt = 1.742, cG = 0.38, cH = 0.718, cL = 0.856, cP0 = 0.927, cSt = 2.654, dG = 15.4e-6, dL = 1.82e-9, epsilon = 0.906) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      //extends PCBilletSchultes(av = 145, epsilon = 0.906, cSt = 2.654, cFlt = 1.742, cH = 0.718, cP0 = 0.927, cL = 0.856, cG = 0.38);
    algorithm
      Column.Di := 0.32;
      Column.Gl := 3989.45 / 3600;
      Column.Gg := 399.77 / 3600;
      Column.RhoL := 998;
      Column.RhoG := 1.188;
      Column.MuL := 0.998e-3;
      Column.MuG := 17.98e-6;
      Column.SigmaL := 72e-3;
      Column.Dl := 1.82e-9;
      Column.Dg := 15.4e-6;
    end UsagePCBilletSchultes3;

    model UsagePCBilletSchultes99
      //Mackowiak pag 99 example bielacki rings metal 1"
      PCBilletSchultes Column(av = 238, Gg(displayUnit = "kg/s"), Gl(displayUnit = "kg/s"), cFlt = 1.856, cG = 0.331, cH = 0.692, cL = 1.461, cP0 = 0.891, cSt = 2.521, epsilon = 0.942, physProp = 0) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      Column.Di := 0.15;
      Column.Vl := 11.1e-3;
      Column.Vg := 1;
      Column.RhoL := 998.2;
      Column.RhoG := 1.17;
      Column.MuL := 1e-3;
      Column.MuG := 1.52e-5 * 1.17;
      Column.SigmaL := 7.24e-2;
      Column.Dl := 5e-9;
      Column.Dg := 2e-5;
    end UsagePCBilletSchultes99;

    model UsagePCBilletSchultes231 "mackowiak pag 231 example Mellapack 250Y"
      extends PCBilletSchultes(av=250,epsilon = 0.975, cSt = 3.157, cFlt = 2.464, cH = 0.554, cP0 = 0.292, cL = 1.1, cG = 0.39, VgFl1(start = 2));
    algorithm
      //Av := 250;
      Di := 1.0;
      Ql := 15.7 / 3600;
      Fg := 1.6;
      RhoL := 998.2;
      RhoG := 1.22;
      MuL := 1e-3;
      MuG := 1.52e-5 * 1.22;
      SigmaL := 7.24e-2;
      Dl := 5e-9;
      Dg := 2e-5;
    end UsagePCBilletSchultes231;

    model PCBilletSchultes_iC4nC4 "Packet towers Billet 1995, pag.344"
      extends Modelica.Icons.Example;
      PCBilletSchultes Column(av=200,P = 360000,VgFl1(start = 2), cFlt = 2.339, cG = 0.39, cH = 0.547, cL = 0.971, cP0 = 0.355, cSt = 3.116, dG = 1.17e-6, dL = 8.44e-9, epsilon = 0.979, lambda = 0.8335, muG = 7.76e-6, muL = 0.13078e-3, physProp = 1, rhoG = 8.398, rhoL = 561.3, sigmaL = 8.88e-3) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    algorithm
      Column.Gg := 279088 / 3600;
      Column.Gl := 261645 / 3600;
    equation
      Column.Fg = Column.FgS1 "working at loading condition";
      annotation(
        Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})));
    end PCBilletSchultes_iC4nC4;

    model UsagePCDelft231 "Mackowiak pag 231 example Mellapack 250Y"
      PCDelft Column(epsilon = 0.975, alpha = pi / 4, hPe = 0.21, hPb = 2, omega = 0.1, VgFl1(start = 2)) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      Column.Sw := 0.017;
      Column.B := 0.0241;
      Column.Di := 1.0;
      Column.Ql := 15.7 / 3600;
      Column.Fg := 1.6;
      Column.RhoL := 998.2;
      Column.RhoG := 1.22;
      Column.MuL := 1e-3;
      Column.MuG := 1.52e-5 * 1.22;
      Column.SigmaL := 7.24e-2;
      Column.Dl := 5e-9;
      Column.Dg := 2e-5;
    end UsagePCDelft231;

    model UsagePCDelft230 "Thesis Erasmus pag 230 example Flexipac 350Y Water-Air"
      PCDelft col(epsilon = 0.985, alpha = pi / 4, hPe = 0.265, hPb = 0.8, omega = 0.1, VgFl1(start = 2)) annotation(
        Placement(visible = true, transformation(origin = {-1, 3}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    algorithm
      col.Sw := 0.0115;
      col.B := 0.0155;
      col.Di := 0.2;
      col.Vl := 35.6 / 3600;
      col.Fg := 2;
      col.RhoL := 763;
      col.RhoG := 1.15;
      col.MuL := 2.31e-3;
      col.MuG := 1.86e-5;
      col.SigmaL := 2.39e-2;
      col.Dl := 5e-9;
      col.Dg := 2e-5;
    end UsagePCDelft230;
  end Examples;
  annotation(
    Documentation(info = "<html><head></head><body>The package contains models for performing distillation and packed column calculations.<div>For distillation there is a basic model performing just the material and energy balances, without any treatment of the necessary number of contact stages. A binary distillation model, and a more complex multicomponent model.</div><div>For packed columns, the Sherwood-Lobo model will perform simple hydraulic calculations. In order to obtain the height of the transfer unit more complex models are needed. For unstructured packing there is the Onda model, and also the very complete Billet-Schultes model. For the last one, the packing characteristics to use can be found at \"Prediction of mass transfer columns with dumped and arranged packaging\", Trans IChemE, Vol.77 Part A, September 1999. For structured packing there is the Delft model.&nbsp;</div></body></html>"));
end Columns;

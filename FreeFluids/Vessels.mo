within FreeFluids;

package Vessels "Vessels.mo by Carlos Trujillo
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

  model VesselFlow
    /*LiquidV
    liquidL*/
    replaceable package MediumL = FreeFluids.LMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "liquid medium in the vessel";
    replaceable package MediumG = Modelica.Media.Air.DryAirNasa constrainedby Modelica.Media.Interfaces.PartialMedium "gas medium in the vessel";
    parameter Modelica.Units.SI.Height inletL=0 "inlet height over lowest internal point of the vessel" annotation(
      Dialog(tab = "Nozzles config"));
    parameter Modelica.Units.SI.Height outletL=0 "outlet height over lowest internal point of the vessel" annotation(
      Dialog(tab = "Nozzles config"));
    parameter Modelica.Units.SI.Height overflowL=0 "overflow height over lowest internal point of the vessel" annotation(
      Dialog(tab = "Nozzles config"));
    parameter Modelica.Units.SI.Height ventOutL=0 "overflow height over lowest internal point of the vessel" annotation(
      Dialog(tab = "Nozzles config"));
    parameter Modelica.Units.SI.ThermalConductance lgConductance=0 "thermal conductange at the liquid-gas interface" annotation(
      Dialog(tab = "Operations"));
    parameter SI.Volume vesselVolume = 1.0 "total volume of the vessel" annotation(
      Dialog(tab = "Initialization"));
    parameter SI.Volume initialLiquidVolume = 1.0 "initial volume of liquid in the vessel" annotation(
      Dialog(tab = "Initialization"));
    final parameter SI.Mass initialLiquidMass = initialLiquidVolume*MediumL.density(MediumL.setState_pTX(initialP,initialT)) "initial vessel liquid mass" annotation(
      Dialog(tab = "Initialization"));
    final parameter SI.Mass initialGasMass = (vesselVolume-initialLiquidVolume)*MediumG.density(MediumG.setState_pTX(initialP,initialT)) "initial vessel gas mass" annotation(
      Dialog(tab = "Initialization"));
    parameter SI.Temperature initialT=298 "initial temperature of the vessel liquid" annotation(
      Dialog(tab = "Initialization"));
    parameter SI.AbsolutePressure initialP = 101325 "initial vessel pressure in the gas phase" annotation(
      Dialog(tab = "Initialization"));
    parameter Modelica.Units.SI.Energy initialLiquidEnthalpy = initialLiquidMass*MediumL.specificEnthalpy(MediumL.setState_pTX(initialP,initialT)) annotation(
      Dialog(tab = "Initialization"));
    parameter Modelica.Units.SI.Energy initialGasEnthalpy = initialGasMass*MediumG.specificEnthalpy(MediumG.setState_pTX(initialP,initialT)) annotation(
      Dialog(tab = "Initialization"));
    parameter Boolean useVentIn=false annotation(
      Dialog(tab = "Nozzles config"));
    Modelica.Units.SI.Volume LiquidV "volume of liquid in the vessel";
    Modelica.Units.SI.Height LiquidL(min=0) "level of liquid in the vessel";
    Modelica.Units.SI.Volume ViTotal "total internal volume";
    FreeFluids.Interfaces.FluidPortA Inlet annotation(
      Placement(visible = true, transformation(origin = {-100, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, -72}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB Outlet annotation(
      Placement(visible = true, transformation(origin = {100, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {88, -72}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB Overflow annotation(
      Placement(visible = true, transformation(origin = {100, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {87, 81}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortA VentIn if useVentIn==true annotation(
      Placement(visible = true, transformation(origin = {-30, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-30, 92}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB VentOut annotation(
      Placement(visible = true, transformation(origin = {30, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {29, 93}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    SI.AbsolutePressure P(start=initialP) "vessel pressure at the gas phase";
    SI.Temperature Tl(start=initialT) "liquid temperature, considered uniform";
    SI.Temperature Tg(start=initialT) "gas temperature, considered uniform";
    MediumL.ThermodynamicState StateL "thermodynamic state of the liquid phase";
    MediumG.ThermodynamicState StateG "thermodynamic state of the gas medium";
    MediumL.SpecificEnthalpy Hl "specific enthalpy of the liquid in the tank";
    Modelica.Units.SI.Energy LiquidEnthalpy(start=initialLiquidEnthalpy) "total enthalpy of the liquid";
    MediumL.SpecificEnthalpy Hg "specific enthalpy of the gas in the tank";
    Modelica.Units.SI.Energy GasEnthalpy(start=initialGasEnthalpy) "total enthalpy of the gas";
    MediumL.Density RhoL "density of the liquid phase";
    MediumG.Density RhoG "density of the gas phase";
    Modelica.Units.SI.Volume GasV "volume of gas phase in the vessel";
    Modelica.Units.SI.Mass LiquidMass(start=initialLiquidMass) "liquid mass in the vessel";
    Modelica.Units.SI.Mass GasMass(start=initialGasMass) "gas mass in the vessel";
  
  equation
    Outlet.Elevation = Inlet.Elevation + outletL - inletL;
    Overflow.Elevation = Inlet.Elevation + overflowL - inletL;
    Outlet.H = Hl;
    Overflow.H = Hl;
    Outlet.X = Inlet.X;
    Overflow.X = Inlet.X;
    StateL = MediumL.setState_phX(P, Hl);
    Tl = MediumL.temperature(StateL);
    RhoL = MediumL.density(StateL);
    LiquidMass=LiquidV*RhoL;
    LiquidEnthalpy=LiquidMass*Hl;
    Inlet.P = if LiquidL > inletL then P + (LiquidL - inletL) * RhoL * g_n else P;
    Outlet.P = if LiquidL > outletL then P + (LiquidL - outletL) * RhoL * g_n else P;
    Overflow.P = if LiquidL > overflowL then P + (LiquidL - overflowL) * RhoL * g_n else P;
    GasV=ViTotal-LiquidV;
    StateG = MediumG.setState_phX(P, Hg);
    Tg = MediumG.temperature(StateG);
    RhoG = MediumG.density(StateG);
    GasMass=GasV*RhoG;
    GasEnthalpy=GasMass*Hg;
    if useVentIn==true then
      VentIn.P=P;
    end if;
    VentOut.P=P;
    VentOut.Elevation = Inlet.Elevation + ventOutL - inletL;
    VentOut.H=Hg;
    VentOut.X=MediumG.X_default;
    der(LiquidMass) = (Inlet.G + Outlet.G + Overflow.G);
    der(LiquidEnthalpy) = Inlet.G * Inlet.H +Outlet.G*Outlet.H +Overflow.G * Overflow.H + lgConductance*(Tg-Tl);
    if useVentIn==true then
      der(GasMass)=VentIn.G+VentOut.G;
      der(GasEnthalpy)=VentIn.G*VentIn.H+VentOut.G*VentOut.H+lgConductance*(Tl-Tg);
    else
      der(GasMass)=VentOut.G;
      der(GasEnthalpy)=VentOut.G*VentOut.H+lgConductance*(Tl-Tg);
    end if;
  annotation(
      Icon(graphics = {Text(origin = {-29, 76}, extent = {{-17, 6}, {17, -6}}, textString = "Vent in"), Text(origin = {30, 77}, extent = {{-18, 7}, {18, -7}}, textString = "Vent out"), Text(origin = {-59, -68}, lineColor = {236, 205, 0}, extent = {{-17, 6}, {17, -6}}, textString = "Inlet"), Text(origin = {52, -68}, lineColor = {236, 205, 0}, extent = {{-18, 6}, {18, -6}}, textString = "Outlet"), Text(origin = {80, 65}, extent = {{-20, 9}, {20, -9}}, textString = "Overflow")}),
      Documentation(info = "<html><head></head><body><div>It contains the equations for tracking the liquid and gas mass contained in the vessel and their thermodynamic properties, considering a complete mixing situation.</div>Declares inlet-outlet ports, for both gas and liquid.<div>The mass and enthalpy of the liquid and gas phases are calculated separately. For the mass the restriction that both phases must fill the vessel volume is introduced.</div><div>For enthalpy, the parameter lgConductance governs the enthalpy transfer between the two phases.</div><div>The model is initialized at initialLiquidMass, initialGasMass, initialP and initialT. In fact, initialT is just used in the calculation of initialLiquidEnthalpy, and initialGasEnthalpy.</div><div>InitialLiquidMass and initalGasMass are calculated from initialLiquidVolume and vesselVolume.</div><div>As the vessel volume can be unknown you can run the model with a vessel volume higher than the initialLiquidVolume, read the vessel volume from the result of the simulation, and run again the model with the real vessel volume.</div></body></html>"));
  end VesselFlow;
  
  model VesselSimple
    extends VesselFlow;
    parameter Modelica.Units.SI.Area section = 0 "uniform area of the vessel" annotation(
      Dialog(tab = "Shape"));
    parameter Modelica.Units.SI.Height height = 0 "height of the vessel" annotation(
      Dialog(tab = "Shape"));
  equation
    ViTotal=section*height;
    LiquidV=LiquidL*section;    
    annotation(
      Icon(graphics = {Rectangle(origin = {0, -30}, fillColor = {122, 208, 251}, fillPattern = FillPattern.VerticalCylinder, extent = {{-80, 50}, {80, -50}}), Rectangle(origin = {0, 10}, extent = {{-80, 90}, {80, -90}})}),
      Documentation(info = "<html><head></head><body>The VesselFlow model using just a fixed area and height for the vessel.</body></html>"));
  end VesselSimple;
  
  //Physical description of the vessel
  //----------------------------------
  model VesselPhysical
    parameter Boolean isVertical = true "if true, the shell is considered vertical, otherwise horizontal" annotation(
      Dialog(tab = "Shell data"));
    parameter Types.VesselForm vesselForm = Types.VesselForm.cylinder annotation(
      Dialog(tab = "Shell data"));
    parameter Modelica.Units.SI.Diameter di = 0 "internal diameter, if cylinder or sphere" annotation(
      Dialog(tab = "Shell data"));
    parameter Modelica.Units.SI.Area section = 0 "if prism" annotation(
      Dialog(tab = "Shell data"));
    parameter Modelica.Units.SI.Length perimeter = 0 "if prism" annotation(
      Dialog(tab = "Shell data"));
    parameter Modelica.Units.SI.Length shellLength = 0 "shell length" annotation(
      Dialog(tab = "Shell data"));
    parameter Types.HeadShape vesselBottom = Types.HeadShape.flat annotation(
      Dialog(tab = "Bottom data"));
    parameter Fraction bFd=if vesselBottom == Types.HeadShape.Klopper then 1.0 elseif vesselBottom == Types.HeadShape.Korbbogen then 0.8 else 0 "ratio of bottom dish radius to tank diameter for torispheral head" annotation(
      Dialog(tab = "Bottom data"));
    parameter Fraction bFk=if vesselBottom == Types.HeadShape.Klopper then 0.1 elseif vesselBottom == Types.HeadShape.Korbbogen then 0.154 else 0 "ratio of bottom knuckle radius to tank diameter for torispheral head" annotation(
      Dialog(tab = "Bottom data"));
    parameter Types.HeadShape vesselTop = Types.HeadShape.flat annotation(
      Dialog(tab = "Top data"));
    parameter Fraction tFd=if vesselTop == Types.HeadShape.Klopper then 1.0 elseif vesselTop == Types.HeadShape.Korbbogen then 0.8 else 0 "ratio of bottom dish radius to tank diameter for torispheral head" annotation(
      Dialog(tab = "Top data"));
    parameter Fraction tFk=if vesselTop == Types.HeadShape.Klopper then 0.1 elseif vesselTop == Types.HeadShape.Korbbogen then 0.154 else 0 "ratio of bottom knuckle radius to tank diameter for torispheral head" annotation(
      Dialog(tab = "Top data"));
    parameter SI.Distance topHeight = 0 "top height if ellipsoidal/conical/piramid. Measured inside" annotation(
      Dialog(tab = "Top data"));
    parameter SI.Distance bottomHeight = 0 "bottom height if ellipsoidal/conical/piramid" annotation(
      Dialog(tab = "Bottom data"));
    parameter SI.Distance shellThickness(displayUnit = "mm") = 0 "shell wall thickness" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Density shellRho(displayUnit = "kg/m3") = 8000 "typical for SS: 8000 kgr/m3, HDPE: 960" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.HeatCapacity shellCp = 500 "typical for SS: 500 J/(kgr·K), HDPE: 2260" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.ThermalConductivity shellK = 15 "shell wall thermal conductivity. typical value for SS=16.7, HDPE:0.33" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance bottomThickness(displayUnit = "mm") = shellThickness "bottom wall thickness" annotation(
      Dialog(tab = "Bottom data"));
    parameter SI.Density bottomRho(displayUnit = "kg/m3") = shellRho "typical for SS: 8000 kgr/m3, HDPE: 960" annotation(
      Dialog(tab = "Bottom data"));
    parameter SI.HeatCapacity bottomCp = shellCp "typical for SS: 500 J/(kgr·K), HDPE: 2260" annotation(
      Dialog(tab = "Bottom data"));
    parameter SI.ThermalConductivity bottompK = shellK "bottom wall thermal conductivity. typical value for SS=16.7, HDPE:0.33" annotation(
      Dialog(tab = "Bottom data"));
    parameter SI.Distance topThickness(displayUnit = "mm") = shellThickness "top wall thickness" annotation(
      Dialog(tab = "Top data"));
    parameter SI.Density topRho(displayUnit = "kg/m3") = shellRho "typical for SS: 8000 kgr/m3, HDPE: 960" annotation(
      Dialog(tab = "Top data"));
    parameter SI.HeatCapacity topCp = shellCp "typical for SS: 500 J/(kgr·K), HDPE: 2260" annotation(
      Dialog(tab = "Top data"));
    parameter SI.ThermalConductivity topK = shellK "top wall thermal conductivity. typical value for SS=16.7, HDPE:0.33" annotation(
      Dialog(tab = "Top data"));
    parameter SI.Distance insulationThickness(displayUnit="mm") = 0 annotation(
      Dialog(tab = "Insulation"));
    parameter SI.ThermalConductivity insulationK = 0.04 "Insulation thermal conductivity. Typical value for glass fiber=0.04 W/(K·m)" annotation(
      Dialog(tab = "Insulation"));
  
    Real Balpha1,Balpha2,Talpha1,Talpha2 "height of the dish and knuckle regions of torispheral heads" ;
    Modelica.Units.SI.Diameter Di "internal diameter, if cylinder or sphere";
    Modelica.Units.SI.Diameter Do "external diameter, if cylinder or sphere";
    Modelica.Units.SI.Area ShellSi "shell internal surface";
    Modelica.Units.SI.Volume ShellVi "shell internal volume";
    Modelica.Units.SI.Area ShellSo "shell external surface";
    Modelica.Units.SI.Volume ShellVo "shell external volume";
    Modelica.Units.SI.Mass ShellMass "shell mass";
    Modelica.Units.SI.Length BottomHi "bottom internal height";
    Modelica.Units.SI.Length BottomHo "bottom external height";
    Modelica.Units.SI.Area BottomSi "bottom internal surface";
    Modelica.Units.SI.Volume BottomVi "bottom internal volume";
    Modelica.Units.SI.Area BottomSo "bottom external surface";
    Modelica.Units.SI.Volume BottomVo "bottom external volume";
    Modelica.Units.SI.Mass BottomMass "bottom mass";
    Modelica.Units.SI.Length TopHi "top internal height";
    Modelica.Units.SI.Length TopHo "top external height";
    Modelica.Units.SI.Area TopSi "top internal surface";
    Modelica.Units.SI.Volume TopVi "top internal volume";
    Modelica.Units.SI.Area TopSo "top external surface";
    Modelica.Units.SI.Volume TopVo "top external volume";
    Modelica.Units.SI.Mass TopMass "top mass";
    Modelica.Units.SI.Area SiTotal "total internal surface";
    Modelica.Units.SI.Volume ViTotal "total internal volume";
    Modelica.Units.SI.Area SoTotal "total external surface";
    Modelica.Units.SI.Volume VoTotal "total external volume";
    Modelica.Units.SI.Mass MassTotal "total vessel mass, without containt";
  equation
    if vesselBottom == Types.HeadShape.Klopper or vesselBottom == Types.HeadShape.Korbbogen then
      Balpha1=bFd*(1-(1-((0.5-bFk)/(bFd-bFk))^2)^0.5);
      Balpha2=bFd-(bFd^2-2*bFd*bFk+bFk-0.25)^0.5;
    else
      Balpha1=0;
      Balpha2=0;
    end if;
    if vesselTop == Types.HeadShape.Klopper or vesselTop == Types.HeadShape.Korbbogen then
      Talpha1=tFd*(1-(1-((0.5-tFk)/(tFd-tFk))^2)^0.5);
      Talpha2=tFd-(tFd^2-2*tFd*tFk+tFk-0.25)^0.5;
    else
      Talpha1=0;
      Talpha2=0;
    end if;
    
    if vesselForm == Types.VesselForm.sphere then
      Di = di;
      Do = di + 2 * shellThickness;
      ShellSi = Modelica.Constants.pi * di ^ 2;
      ShellVi = Modelica.Constants.pi * di ^ 3 / 6;
      ShellSo = Modelica.Constants.pi * Do ^ 2;
      ShellVo = Modelica.Constants.pi * Do ^ 3 / 6;
      BottomHi = 0;
      BottomSi = 0;
      BottomVi = 0;
      BottomHo = 0;
      BottomSo = 0;
      BottomVo = 0;
      TopHi = 0;
      TopSi = 0;
      TopVi = 0;
      TopHo = 0;
      TopSo = 0;
      TopVo = 0;
    elseif vesselForm == Types.VesselForm.cylinder then
      Di = di;
      Do = di + 2 * shellThickness;
      ShellSi = Modelica.Constants.pi * di * shellLength;
      ShellVi = Modelica.Constants.pi * di ^ 2 * shellLength / 4;
      ShellSo = Modelica.Constants.pi * Do * shellLength;
      ShellVo = Modelica.Constants.pi * Do ^ 2 * shellLength / 4;
      if vesselTop == Types.HeadShape.open then
        TopHi = 0;
        TopSi = 0;
        TopVi = 0;
        TopHo = 0;
        TopSo = 0;
        TopVo = 0;
      elseif vesselTop == Types.HeadShape.flat then
        TopHi = 0;
        TopSi = Modelica.Constants.pi * di ^ 2 / 4;
        TopVi = 0;
        TopHo = topThickness;
        TopSo = Modelica.Constants.pi * Do ^ 2 / 4 "the intersect volume is considered as of the top";
        TopVo = Modelica.Constants.pi * Do ^ 2 / 4 * topThickness;
      elseif vesselTop == Types.HeadShape.conical then
        TopHi = topHeight;
        TopSi = Modelica.Constants.pi * di * (di ^ 2 / 4 + topHeight ^ 2) ^ 0.5;
        TopVi = Modelica.Constants.pi * di ^ 2 * topHeight / 12;
        TopHo = topHeight+topThickness*(di^2/4+topHeight^2)^0.5/(di/2);
        TopSo = Modelica.Constants.pi * (di + 2 * topThickness) * ((di + 2 * topThickness) ^ 2 / 4 + topHeight ^ 2) ^ 0.5;
        TopVo = Modelica.Constants.pi * (di + 2 * topThickness) ^ 2 * (topHeight + topThickness) / 12;
      //elseif vesselTop == Types.HeadShape.Klopper or vesselTop == Types.HeadShape.Korbbogenthen then
        
      elseif vesselTop == Types.HeadShape.Klopper then
        TopHi = Talpha2*di;//0.19377 * di;
        TopSi = 0.947 * di ^ 2;
        TopVi = 0.098966 * di ^ 3;
        TopHo = Talpha2*(di + 2 * topThickness);//0.19377 * (di + 2 * topThickness);
        TopSo = 0.947 * (di + 2 * topThickness) ^ 2;
        TopVo = 0.098966 * (di + 2 * topThickness) ^ 3;
      elseif vesselTop == Types.HeadShape.Korbbogen then
        TopHi = Talpha2*di;//0.2545 * di;
        TopSi = 0.986 * di ^ 2;
        TopVi = 0.1307 * di ^ 3;
        TopHo = Talpha2*(di + 2 * topThickness);//0.2545 * (di + 2 * topThickness);
        TopSo = 0.986 * (di + 2 * topThickness) ^ 2;
        TopVo = 0.1307 * (di + 2 * topThickness) ^ 3;
      elseif vesselTop == Types.HeadShape.semielliptical then
        TopHi = 0.25 * di;
        TopSi = 1.084 * di ^ 2;
        TopVi = Modelica.Constants.pi * di ^ 3 / 24;
        TopHo = 0.25 * (di + 2 * topThickness);
        TopSo = 1.084 * (di + 2 * topThickness) ^ 2;
        TopVo = Modelica.Constants.pi * (di + 2 * topThickness) ^ 3 / 24;
      elseif vesselTop == Types.HeadShape.hemispherical then
        TopHi = 0.5 * di;
        TopSi = Modelica.Constants.pi * di ^ 2 / 2;
        TopVi = Modelica.Constants.pi * di ^ 3 / 12;
        TopHo = 0.5 * (di + 2 * topThickness);
        TopSo = Modelica.Constants.pi * (di + 2 * topThickness) ^ 2 / 2;
        TopVo = Modelica.Constants.pi * (di + 2 * topThickness) ^ 3 / 12;
      end if;
      if vesselBottom == Types.HeadShape.flat then
        BottomHi = 0;
        BottomSi = Modelica.Constants.pi * di ^ 2 / 4;
        BottomVi = 0;
        BottomHo = bottomThickness;
        BottomSo = Modelica.Constants.pi * di ^ 2 / 4 "the intersect volume is considered as of the shell";
        BottomVo = Modelica.Constants.pi * di ^ 2 / 4 * bottomThickness;
      elseif vesselBottom == Types.HeadShape.conical then
        BottomHi = bottomHeight;
        BottomSi = Modelica.Constants.pi * di * (di ^ 2 / 4 + bottomHeight ^ 2) ^ 0.5;
        BottomVi = Modelica.Constants.pi * di ^ 2 * bottomHeight / 12;
        BottomHo = bottomHeight+bottomThickness*(di^2/4+bottomHeight^2)^0.5/(di/2);
        BottomSo = Modelica.Constants.pi * (di + 2 * bottomThickness) * ((di + 2 * bottomThickness) ^ 2 / 4 + bottomHeight ^ 2) ^ 0.5;
        BottomVo = Modelica.Constants.pi * (di + 2 * bottomThickness) ^ 2 * (bottomHeight + bottomThickness) / 12;
      elseif vesselBottom == Types.HeadShape.Klopper then
        BottomHi = Balpha2*di;//0.19377 * di;
        BottomSi = 0.947 * di ^ 2;
        BottomVi = 0.098966 * di ^ 3;
        BottomHo = Balpha2*(di + 2 * bottomThickness);//0.19377 * (di + 2 * bottomThickness);
        BottomSo = 0.947 * (di + 2 * bottomThickness) ^ 2;
        BottomVo = 0.098966 * (di + 2 * bottomThickness) ^ 3;
      elseif vesselBottom == Types.HeadShape.Korbbogen then
        BottomHi = Balpha2*di;//0.25447 * di;
        BottomSi = 0.986 * di ^ 2;
        BottomVi = 0.1307 * di ^ 3;
        BottomHo = Balpha2*(di + 2 * bottomThickness);//0.25447 * (di + 2 * bottomThickness);
        BottomSo = 0.986 * (di + 2 * bottomThickness) ^ 2;
        BottomVo = 0.1307 * (di + 2 * bottomThickness) ^ 3;
      elseif vesselBottom == Types.HeadShape.semielliptical then
        BottomHi = di / 4;
        BottomSi = 1.084 * di ^ 2;
        BottomVi = Modelica.Constants.pi * di ^ 3 / 24;
        BottomHo = di / 4 + bottomThickness;
        BottomSo = 1.084 * (di + 2 * bottomThickness) ^ 2;
        BottomVo = Modelica.Constants.pi * (di + 2 * bottomThickness) ^ 3 / 24;
      elseif vesselBottom == Types.HeadShape.hemispherical then
        BottomHi = 0.5 * di;
        BottomSi = Modelica.Constants.pi * di ^ 2 / 2;
        BottomVi = Modelica.Constants.pi * di ^ 3 / 12;
        BottomHo = 0.5 * di + bottomThickness;
        BottomSo = Modelica.Constants.pi * (di + 2 * bottomThickness) ^ 2 / 2;
        BottomVo = Modelica.Constants.pi * (di + 2 * bottomThickness) ^ 3 / 12;
      end if;
    else
      assert(false, "bad size specification for the vessel");
    end if;
    ShellMass = (ShellVo - ShellVi) * shellRho;
    BottomMass = (BottomVo - BottomVi) * bottomRho;
    TopMass = (TopVo - TopVi) * topRho;
    SiTotal = ShellSi + TopSi + BottomSi;
    ViTotal = ShellVi + TopVi + BottomVi;
    SoTotal = ShellSo + TopSo + BottomSo;
    VoTotal = ShellVo + TopVo + BottomVo;
    MassTotal = ShellMass + TopMass + BottomMass;

    annotation(
      Icon(graphics = {Rectangle(origin = {0, -30}, fillColor = {0, 170, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-80, 50}, {80, -50}}), Rectangle(origin = {0, 10}, extent = {{-80, 90}, {80, -90}})}),
      Documentation(info = "<html><head></head><body>Contains a detailed description of common vessel shapes. It allows the calculation of volumes and walls mass.</body></html>"));
  end VesselPhysical;
  
  model VesselLevel
    extends VesselPhysical;
    parameter Boolean fixLiquidLevel = true "if true, liquid height will be fixed to liquidL";
    parameter SI.Height liquidL = 0 "fixed liquid level over lowest vessel point";
    parameter Boolean fixLiquidVolume = false "if true, liquid volume will be fixed to liquidV. Do not fix height and volume simultaneously";
    parameter SI.Volume liquidV = 0 "fixed liquid volume";
    Modelica.Units.SI.Volume LiquidV "volume of liquid in the vessel";
    Modelica.Units.SI.Height LiquidL(min=0) "level of liquid in the vessel";
    Real Alpha(min=0) "adimensional level of liquid: LiquidL/Di";
  equation  
    if fixLiquidLevel == true then
      LiquidL = liquidL;
    elseif fixLiquidVolume == true then
      LiquidV = liquidV;
    end if;
  
    if vesselForm == Types.VesselForm.cylinder then
      Alpha=LiquidL/Di;
    else
      Alpha=0;
    end if;
    
    if vesselForm == Types.VesselForm.sphere then
      if LiquidL <= Di / 2 then
        LiquidV = Modelica.Constants.pi * LiquidL ^ 2 / 4 * (2 * Di - 4 * LiquidL / 3);
      elseif LiquidL <= Di then
        LiquidV = Modelica.Constants.pi * (Di ^ 3 / 6-(Di-LiquidL) ^ 2 / 4 * (2 * Di - 4 * (Di-LiquidL) / 3));
      else
        LiquidV=0;
      end if;
    elseif vesselForm == Types.VesselForm.cylinder then
      if isVertical == true then
        if vesselBottom == Types.HeadShape.flat then
          LiquidV = Modelica.Constants.pi * Di ^ 2 / 4 * LiquidL;
        elseif vesselBottom == Types.HeadShape.conical then
          if LiquidL <= bottomHeight then
            LiquidV = Modelica.Constants.pi / 4 * (Di * LiquidL / bottomHeight) ^ 2 * LiquidL / 3;
          elseif LiquidL <= bottomHeight + shellLength then
            LiquidV = Modelica.Constants.pi * Di ^ 2 / 4 * (LiquidL - 2 * bottomHeight / 3);
          else
            LiquidV=0;
          end if;
        elseif vesselBottom == Types.HeadShape.Klopper or vesselBottom == Types.HeadShape.Korbbogen then
          if LiquidL <= Balpha1 * Di then
            LiquidV = Modelica.Constants.pi *(bFd*LiquidL^2*Di-LiquidL^3/3) "pi*Di^3*(bFd*Alpha^2-Alpha^3/3)";
          elseif LiquidL <= Balpha2 * Di then
            LiquidV = Modelica.Constants.pi *Di^3*((bFd*Balpha1^2-Balpha1^3/3)+(((0.5-bFk)^2+bFk^2)*(Alpha-Balpha1)-((Alpha-Balpha2)^3-(Balpha1-Balpha2)^3)/3+(0.5-bFk)*((Alpha-Balpha2)*(bFk^2-(Alpha-Balpha2)^2)^0.5-(Balpha1-Balpha2)*(bFk^2-(Balpha1-Balpha2)^2)^0.5+bFk^2*asin((Alpha-Balpha2)/bFk)-bFk^2*asin((Balpha1-Balpha2)/bFk))));
          elseif LiquidL <= (Balpha2 * Di+shellLength) then
            LiquidV=Modelica.Constants.pi *Di^3*((bFd*Balpha1^2-Balpha1^3/3)+(((0.5-bFk)^2+bFk^2)*(Balpha2-Balpha1)+(Balpha1-Balpha2)^3/3+(0.5-bFk)*(-(Balpha1-Balpha2)*(bFk^2-(Balpha1-Balpha2)^2)^0.5-bFk^2*asin((Balpha1-Balpha2)/bFk)))+0.25*(Alpha-Balpha2));
          else
            LiquidV=0;
          end if;
        elseif vesselBottom == Types.HeadShape.semielliptical then
          if LiquidL <= bottomHeight then
            LiquidV = Modelica.Constants.pi / 4 * (Di * LiquidL / bottomHeight) ^ 2 * (bottomHeight - LiquidL / 3);
          elseif LiquidL <= bottomHeight + shellLength then
            LiquidV = Modelica.Constants.pi * Di ^ 2 / 4 * (LiquidL - bottomHeight / 3);
          else
            LiquidV = 0;
          end if;
        elseif vesselBottom == Types.HeadShape.hemispherical then
          if LiquidL <= Di / 2 then
            LiquidV = Modelica.Constants.pi * LiquidL ^ 2 / 4 * (2 * Di - 4 * LiquidL / 3);
          elseif LiquidL <= Di / 2 + shellLength then
            LiquidV = Modelica.Constants.pi / 4 * (Di ^ 3 / 12 - Di ^ 3 / 4 + LiquidL * Di ^ 2);
          else
            LiquidV = 0;
          end if;
        end if;
      else
        if vesselBottom == Types.HeadShape.flat then
        elseif vesselBottom == Types.HeadShape.conical then
  
        elseif vesselBottom == Types.HeadShape.Klopper then
  
        elseif vesselBottom == Types.HeadShape.Korbbogen then
  
        elseif vesselBottom == Types.HeadShape.semielliptical then
  
        elseif vesselBottom == Types.HeadShape.hemispherical then
        end if;
      end if;
    end if;
  annotation(
      Documentation(info = "<html><head></head><body>Adds the relationship between liquid filled level and volume to the VesselPhysical model.<div>For now the relationship between liquid level and liquid volume is coded only for spheric and vertical cylindrical vessels.<br><div>For vertical vessels the relationship covers only the bottom and the cylindrical shell.</div></div></body></html>"));end VesselLevel;

  model VesselDetailed
    extends VesselFlow;
    extends VesselLevel(final fixLiquidLevel = false, final fixLiquidVolume=false, final liquidL=0, final liquidV=0, LiquidV(start=initialLiquidMass));
  annotation(
      Documentation(info = "<html><head></head><body><div>Combines the tracking capabilities of the VesselFlow model with a detailed description of the vessel.</div></body></html>"));
  end VesselDetailed;

  model VesselCylVert
    parameter SI.Distance wallThickness(displayUnit = "mm") = 0 "wall thickness in m" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Distance thicknessInsul = 0 annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Density wallRho(displayUnit = "kg/m3") = 8000 "typical for SS: 8000 kgr/m3, HDPE: 960" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.HeatCapacity wallCp = 500 "typical for SS: 500 J/(kgr·K), HDPE: 2260" annotation(
      Dialog(tab = "Heat transfer"));
    parameter SI.ThermalConductivity wallK = 15 "Wall thermal conductivity. typical value for SS=16.7, HDPE:0.33" annotation(
      Dialog(tab = "Heat transfer"));
    parameter String bottomShape = "Klopper" "alternatives are Flat, Klopper, Korbbogen and Conical" annotation(
      Dialog(tab = "Physical data"));
    parameter String topShape = "Klopper" "alternatives are Flat, Korboggen and Conical" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Length coneH = 0 "cone internal heigth, if conical shape is used" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Diameter di = 0 "vessel internal diameter" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Length cylinderH = 0 "length of the cylindrical wall" annotation(
      Dialog(tab = "Physical data"));
    //parameter SI.Distance hLiquid = 0 "heigth of liquid from lower most point";
    //parameter SI.ThermalInsulance foulingF = 0 "Process side fouling factor";
    //parameter Integer numCoils = 1;
    parameter Integer nBaffles = 0 "number of baffles" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Distance baffleWidth = 0 annotation(
      Dialog(tab = "Physical data"));
    parameter Integer nConcCoils = 1 "number of concentric coils" annotation(
      Dialog(tab = "Heat transfer"));
    SI.Volume VtotalIn, VtotalOut "internal/external vessel volume";
    SI.Mass VesselMass;
    SI.Diameter Do "vessel external diameter";
    SI.Distance HbottomIn, HbottomOut "bottom head internal/external height";
    SI.Distance HtopIn, HtopOut "top head internal/external height";
    SI.Distance HtotalIn;
    SI.Area SbottomIn, SbottomOut "bottom head internal/external surface";
    SI.Area StopIn, StopOut "top head internal/external surface";
    SI.Volume VbottomIn, VbottomOut "bottom head volume";
    SI.Volume VtopIn, VtopOut "bottom head volume";
    SI.Area ScylinderIn, ScylinderOut "cylinder internal/external surface";
    SI.Volume VcylinderIn, VcylinderOut "cylinder internal/external volume";
  algorithm
    Do := di + 2 * wallThickness;
    if bottomShape == "Klopper" then
      HbottomIn := 0.19377 * di;
      SbottomIn := 0.947 * di ^ 2;
      VbottomIn := 0.098966 * di ^ 3;
      HbottomOut := 0.19377 * Do;
      SbottomOut := 0.947 * Do ^ 2;
      VbottomOut := 0.098966 * Do ^ 3;
    elseif bottomShape == "Korbbogen" then
      HbottomIn := 0.2544 * di;
      SbottomIn := 0.986 * di ^ 2;
      VbottomIn := 0.1307 * di ^ 3;
      HbottomOut := 0.2544 * Do;
      SbottomOut := 0.986 * Do ^ 2;
      VbottomOut := 0.1307 * Do ^ 3;
    elseif bottomShape == "Conical" then
      HbottomIn := coneH;
      SbottomIn := pi * di * (coneH ^ 2 + di ^ 2 / 4) ^ 0.5;
      VbottomIn := pi * di ^ 2 * coneH / 12;
      HbottomOut := coneH + wallThickness;
      SbottomOut := pi * Do * (HbottomOut ^ 2 + Do ^ 2 / 4) ^ 0.5;
      VbottomOut := pi * Do ^ 2 * HbottomOut / 12;
    else
      HbottomIn := 0;
      SbottomIn := pi * di ^ 2 / 4;
      VbottomIn := 0;
      HbottomOut := wallThickness;
      SbottomOut := pi * Do ^ 2 / 4;
      VbottomIn := SbottomOut * wallThickness;
    end if;
    if topShape == "Klopper" then
      HtopIn := 0.19377 * di;
      StopIn := 0.947 * di ^ 2;
      VtopIn := 0.098966 * di ^ 3;
      HtopOut := 0.19377 * Do;
      StopOut := 0.947 * Do ^ 2;
      VtopOut := 0.098966 * Do ^ 3;
    elseif topShape == "Korbbogen" then
      HtopIn := 0.2544 * di;
      StopIn := 0.986 * di ^ 2;
      VtopIn := 0.1307 * di ^ 3;
      HtopOut := 0.2544 * Do;
      StopOut := 0.986 * Do ^ 2;
      VtopOut := 0.1307 * Do ^ 3;
    elseif topShape == "Conical" then
      HtopIn := coneH;
      StopIn := pi * di * (coneH ^ 2 + di ^ 2 / 4) ^ 0.5;
      VtopIn := pi * di ^ 2 * coneH / 12;
      HtopOut := coneH + wallThickness;
      StopOut := pi * Do * (HtopOut ^ 2 + Do ^ 2 / 4) ^ 0.5;
      VtopOut := pi * Do ^ 2 * HtopOut / 12;
    else
      HtopIn := 0;
      StopIn := pi * di ^ 2 / 4;
      VtopIn := 0;
      HtopOut := wallThickness;
      StopOut := pi * Do ^ 2 / 4;
      VtopIn := StopOut * wallThickness;
    end if;
    HtotalIn := HbottomIn + cylinderH + HtopIn;
    ScylinderIn := pi * di * cylinderH;
    VcylinderIn := pi * di ^ 2 / 4 * cylinderH;
    ScylinderOut := pi * Do * cylinderH;
    VcylinderOut := pi * Do ^ 2 / 4 * cylinderH;
    VtotalIn := VbottomIn + VtopIn + VcylinderIn;
    VtotalOut := VbottomOut + VtopOut + VcylinderOut;
    VesselMass := (VtotalOut - VtotalIn) * wallRho;
    annotation(
      defaultComponentName = "vessel",
      Icon(coordinateSystem(extent = {{-150, -150}, {150, 150}}, initialScale = 0.2), graphics = {Rectangle(lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, extent = {{-150, 150}, {150, -150}}), Rectangle(fillColor = {85, 170, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-150, -150}, {150, -100}}), Line(points = {{-150, 150}, {-150, -150}, {150, -150}, {150, 150}})}),
      Documentation(info = "<html><head></head><body>Contains the physical description of a vertical cylindrical vessel. It is the base for the agitated tanks in the package. In the future it should be replaced by an extension of the model &nbsp;VesselDetailed.</body></html>"));
  end VesselCylVert;

  //***MIXERS***
  //************
  //General mixer

  partial model MixerBase
    parameter Modelica.Units.SI.Diameter d = 0 "impeller diameter";
    parameter Integer nImpellers = 0;
    parameter Modelica.Units.SI.Distance interImpellerDistance;
    parameter Modelica.Units.SI.Power motorPower = 0;
    parameter Boolean useFixedSpeed = true;
    parameter Modelica.Units.NonSI.AngularVelocity_rpm n(start = 60) "shaft rotational speed rpm";
    parameter Boolean doShaftCalc = false "perform or not shaft diameter calculation";
    parameter Boolean footBearing = false "bottom bearing for the shaft";
    parameter Real impellerHydraulicF = 4 "1.5 normal, 2.5 boiling, 3 gas sparge, 4 large solid volume, 6 startup with settled solids";
    parameter Modelica.Units.SI.ShearStress shaftDTS = 68.9e6 "shaft design tensile stress. AISI316L=60e6. AISI316=68.9e6";
    parameter Modelica.Units.SI.ShearStress shaftDSS = 41.4e6 "shaft design shear stress. AISI316L=35.9e6. AISI316=41.4e6";
    Modelica.Units.SI.Distance Htop "distance from top bearing to lower impeller";
    Modelica.Units.SI.Frequency N "rotational speed 1/sec";
    Modelica.Units.SI.ReynoldsNumber Re;
    Real Np(start = 1) "Power number";
    Modelica.Units.SI.Power W(start = 1e4) "nett absorbed power";
    Modelica.Units.SI.Power WimpellerCorr(start = 1e4) "corrected impeller power";
    Modelica.Units.SI.Force Fimpeller(start = 1000) "bending force at impeller";
    Modelica.Units.SI.Torque Torque(start = 3000) "shaftTorque";
    Modelica.Units.SI.Velocity V "peripheral speed";
    Modelica.Units.SI.Force ForceOnBottomBearing(start = 1000);
    Modelica.Units.SI.Torque BendingMomentMax(start = 1000) "maximum bending moment";
    Modelica.Units.SI.Diameter DiamShaftShear(start = 0.1) "minimum shaft diameter due to shear";
    Modelica.Units.SI.Diameter DiamShaftTension(start = 0.1) "minimum shaft diameter due to tension";
    Modelica.Units.SI.Diameter DiamShaftMinimum(start = 0.1) "minimum shaft diameter";
  equation
    if useFixedSpeed == true then
      N = n / 60;
    end if;
    V = N * pi * d;
    Torque = W / (2 * pi * N * 0.85) "0.85 is just a safety factor";
    WimpellerCorr = W / (nImpellers * 0.85);
//end if;
    if doShaftCalc == true then
      Fimpeller = 0.048 * WimpellerCorr * impellerHydraulicF / (N * d);
      DiamShaftShear = (16 * (BendingMomentMax ^ 2 + Torque ^ 2) ^ 0.5 / (pi * shaftDSS)) ^ (1 / 3);
      DiamShaftTension = (16 * (BendingMomentMax + (BendingMomentMax ^ 2 + Torque ^ 2) ^ 0.5) / (pi * shaftDTS)) ^ (1 / 3);
      if DiamShaftShear > DiamShaftTension then
        DiamShaftMinimum = DiamShaftShear;
      else
        DiamShaftMinimum = DiamShaftTension;
      end if;
    else
      Fimpeller = 0;
      DiamShaftShear = 0;
      DiamShaftTension = 0;
      DiamShaftMinimum = 0;
    end if;
    annotation(
      defaultComponentName = "Mixer",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(origin = {0, -2}, lineColor = {71, 71, 71}, fillColor = {191, 191, 191}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-90, -32}, {0, 2}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, 0}, {150, 50}}, textString = "%name"), Ellipse(origin = {90, -2}, lineColor = {71, 71, 71}, fillColor = {191, 191, 191}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-90, -32}, {0, 2}}, endAngle = 360)}),
      Documentation(info = "<html><head></head><body>Common description of a mixer.</body></html>"));
  end MixerBase;

  //Anchor mixer model
  //------------------

  model MixerAnchor "Anchor mixer"
    extends MixerBase(nImpellers = 1);
    parameter Integer nBlades = 2 "number of blades";
    parameter SI.Distance width "impeler frontal width(depending on impeller angle), not real";
    parameter SI.Height hTotal "impeller height";
    SI.Height Hwetted "wetted height of the anchor";
    Real Z "correction factor for Np for Young-Sei Lee correlation: Journal of Korean Institute of Chemical Engineering Vol.39 n.5 October 2001";
  algorithm
    Z := width / Hwetted + 0.684 * (nBlades * log(d / (d - 2 * width))) ^ 0.139;
  equation
    if doShaftCalc == true then
      if footBearing == false then
        ForceOnBottomBearing = 0;
        BendingMomentMax = Fimpeller * (Htop - hTotal / 2);
      else
        ForceOnBottomBearing = Fimpeller * (Htop - hTotal / 2) / Htop;
        BendingMomentMax = ForceOnBottomBearing * hTotal - Fimpeller * hTotal / 2;
      end if;
    else
      ForceOnBottomBearing = 0;
      BendingMomentMax = 0;
    end if;
  end MixerAnchor;

  //Kamei mixer model
  //-----------------

  model MixerKamei "From Furukawa et alt., International Journal of Chemical Engineering volume 2012, based on Kamei et al.(1995/6/7) model for many different impeller types"
    extends MixerBase;
    parameter Integer nBlades = 0 "number of blades. Not necessary for mig or intermig(2)";
    parameter SI.Distance width "average real width of the blades,not projected to the axis. Not necessary for mig or intermig";
    parameter SI.Height hBottom = 0 "Distance from lower impeller to vessel bottom";
    parameter SI.Angle angle(displayUnit = "deg") "average blade angle over horizontal. Not necessary for mig or intermig";
    parameter String reference = "hydrofoil" "alternatives are: hydrofoil, propeller, axialTurbine,radialTurbine,intermig,mig";
    //is used for Nusselt number calculation and for some algorithm selection
    Real NpMax "Power number in fully baffled condition";
    Real eta, beta, gamma, X, Ct, fInf, Ctr, m, Cl, ReG, f, Np0, NpB;
    //Np0=Power number in unbaffled conditions
    SI.VolumeFlowRate Q "flow induced by the impellers";
    SI.Velocity Vc "characteristic velocity of the flow produced by the impellers";
    Real ScaleOfAgitation "2 for 0.1 sg difference suspension, 10 for 1.0 difference";
    //Scale of agitation is char.velocity in ft/min divided by 6
    Real Nq "flow number";
  equation
    if footBearing == false then
      ForceOnBottomBearing = 0;
      BendingMomentMax = Fimpeller * (nImpellers * (Htop - (nImpellers - 1) / 2 * interImpellerDistance));
    else
      ForceOnBottomBearing = Fimpeller * (nImpellers * (Htop - (nImpellers - 1) / 2 * interImpellerDistance)) / (Htop + hBottom);
      BendingMomentMax = ForceOnBottomBearing * (hBottom + (nImpellers - 1) * interImpellerDistance) - Fimpeller * nImpellers * (nImpellers - 1) / 2 * interImpellerDistance;
    end if;
  annotation(
      Documentation(info = "<html><head></head><body>
<p><span>From Furukawa et alt., International Journal of Chemical Engineering volume 2012, based on Kamei et al.(1995/6/7) model for many different impeller types</span></p></body></html>"));
  end MixerKamei;

  //***TANKS***
  //***********

  partial model TankPM "tank is the content of a vessel"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    replaceable VesselCylVert Vessel;
    parameter SI.Distance hLiquid = 0 "heigth of liquid from lower most point";
    parameter SI.ThermalInsulance foulingF = 0 "Process side fouling factor" annotation(
      Dialog(tab = "Heat transfer"));
    parameter SI.AbsolutePressure fixedPressure = 30e5;
    parameter Integer numCoils = 0 "number of coils";
    parameter Integer numHalfCoils = 0 "number of halfcoils";
    parameter Boolean useHTWallCorrFactor = true annotation(
      Dialog(tab = "Heat transfer"));
    Medium.Temperature T "vessel liquid temperature";
    Medium.Density Rho "vessel liquid bulk density";
    Medium.DynamicViscosity Mu(min = 1.0e-6, start = 1.0e-3, max = 1.0e6) "vessel liquid bulk viscosity";
    SI.HeatCapacity Cp "vessel liquid Cp";
    Medium.ThermalConductivity K "vessel liquid thermal conductivity";
    SI.PrandtlNumber Pr;
    //Real HTWallCorrExp;
    SI.Volume Vliquid "without substracting volume of internal parts";
    SI.Volume Vuse "useful internal volume, bottom plus cylinder";
    Medium.ThermodynamicState StateV "thermodynamic state of vessel content";
    FreeFluids.Various.TemperatureOutput Tmass annotation(
      Placement(visible = true, transformation(origin = {-20, 40}, extent = {{-8, -8}, {8, 8}}, rotation = 90), iconTransformation(origin = {36, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //replaceable FreeFluids.Vessels.InitializationParamSP InitialValues;
  algorithm
    Vuse := Vessel.VbottomIn + Vessel.VcylinderIn;
    if hLiquid >= Vessel.HbottomIn then
      Vliquid := Vessel.VbottomIn + pi * Vessel.di ^ 2 / 4 * (hLiquid - Vessel.HbottomIn);
    else
      Vliquid := 0;
    end if;
  equation
    StateV = Medium.setState_pTX(fixedPressure, T);
    Rho = Medium.density(StateV);
    Mu = Medium.dynamicViscosity(StateV);
    Cp = Medium.specificHeatCapacityCp(StateV);
    K = Medium.thermalConductivity(StateV);
    Pr = Cp * Mu / K;
    Tmass = T;
//HTWallCorrExp:=0.1*(8.61e-2*Mu)^(-0.21);
    annotation(
      Icon(graphics = {Rectangle(origin = {-7, -1}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-55, 71}, {67, -55}}), Ellipse(origin = {-12, -55}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-50, 23}, {72, -23}}, endAngle = 360), Rectangle(origin = {-2, 22}, extent = {{-2, 58}, {6, -78}})}, coordinateSystem(initialScale = 0.1)));
  end TankPM;

  //***AGITATED TANK NAGATA MODELS***
  //---------------------------------
  //***AGITATED TANK ANCHOR MODELS***
  //---------------------------------

  partial model TankAgitatedAnchorPM "Anchor Vessel. Seems the best model for anchors"
    extends TankPM;
    replaceable MixerAnchor Mixer;
  algorithm
    if hLiquid > Mixer.hTotal then
      Mixer.Hwetted := Mixer.hTotal;
    else
      Mixer.Hwetted := hLiquid;
    end if;
  equation
    Mixer.Re = Mixer.d ^ 2 * Mixer.N * Rho / Mu;
    Mixer.W = Mixer.nImpellers * Mixer.Np * Rho * Mixer.d ^ 5 * Mixer.N ^ 3;
    Mixer.Htop = Vessel.HtotalIn;
//Power number calculation
//From ChemicalProcessing.com
    if Mixer.Re <= 15 then
      Mixer.Np = 113 / Mixer.Re * Mixer.Hwetted / Mixer.d * (Vessel.di / (Vessel.di - Mixer.d)) ^ 0.5 * (Mixer.width / Mixer.d) ^ 0.16 * (Mixer.nBlades / 2) ^ 0.67;
    elseif Mixer.Re < 10000 then
      Mixer.Np = 113 / Mixer.Re * Mixer.Hwetted / Mixer.d * (Vessel.di / (Vessel.di - Mixer.d)) ^ 0.5 * (Mixer.width / Mixer.d) ^ 0.16 * (Mixer.nBlades / 2) ^ 0.67 * exp(0.013 - 0.17 * log(Mixer.Re) + 0.061 * log(Mixer.Re) ^ 2);
    else
      Mixer.Np = 113 / Mixer.Re * Mixer.Hwetted / Mixer.d * (Vessel.di / (Vessel.di - Mixer.d)) ^ 0.5 * (Mixer.width / Mixer.d) ^ 0.16 * (Mixer.nBlades / 2) ^ 0.67 * 0.0065 * Mixer.Re ^ 0.94;
    end if;
    annotation(
      Icon(graphics = {Rectangle(origin = {-7, -53}, extent = {{-35, 1}, {47, -3}}), Rectangle(origin = {-40, -17}, extent = {{-2, -35}, {2, 59}}), Rectangle(origin = {38, -17}, extent = {{2, -35}, {-2, 59}})}));
  end TankAgitatedAnchorPM;

  model TankAgitatedAnchor "Anchor Vessel. Seems the best model for anchors"
    //Normally missing:N (rpm) from the mixer, and T from the process
    extends TankAgitatedAnchorPM;
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end TankAgitatedAnchor;

  model TankAgitAnchor1Coil "Anchor Vessel. Seems the best model for anchors"
    //Normally missing:N (rpm) from the mixer, and T from the process
    extends TankAgitatedAnchorPM(numCoils = 1, numHalfCoils = 0);
    replaceable VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorC1(start = 1.0);
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e3) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1 "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    NuCoil1 := homotopy(0.135 * HTWallCorrFactorC1 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.135 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33) "VDI atlas. Original leading constant 0.084 and 0.14 for viscosity corr. exponent";
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := homotopy((Coil1.Tb - Coil1.Ta) / log((T - Coil1.Ta) / (T - Coil1.Tb)), T - Coil1.Ta);
    end if;
  equation
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    if numCoils == 1 and numHalfCoils == 0 then
      W = -Coil1.W;
    end if;
  end TankAgitAnchor1Coil;

  model TankAgitAnchor2Coils
    extends TankAgitAnchor1Coil(Vessel(nConcCoils = 2), numCoils = 2, numHalfCoils = 0);
    replaceable VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorC2(start = 1.0);
    SI.DynamicViscosity MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil2 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil2(start = 1e2, min = 1, max = 1e5) "Coil2 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil2(min = 10, start = 1000) "Coil2 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil2(min = 10, start = 600) "Coil2 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil2 "Logarithmic mean temperature difference";
    Medium.ThermodynamicState StateC2;
  algorithm
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.13);
    end if;
    NuCoil2 := homotopy(0.135 * HTWallCorrFactorC2 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.135 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta - (T - Coil2.Tb) == 0 then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := homotopy((T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb)), T - Coil2.Ta);
    end if;
  equation
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    if numCoils == 2 and numHalfCoils == 0 then
      W = (-Coil1.W) - Coil2.W;
    end if;
  end TankAgitAnchor2Coils;

  model TankAgitAnchor2CoilsOld "Anchor Vessel. Seems the best model for anchors"
    //Normally missing:N (rpm) from the mixer, and T from the process
    extends TankAgitatedAnchorPM(final numCoils = 2, final numHalfCoils = 0);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorC1(start = 1.0), HTWallCorrFactorC2(start = 1.0);
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(start = 1e2, min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.NusseltNumber NuCoil2(start = 1e2, min = 1, max = 1e5) "Coil2 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Hcoil2(min = 10, start = 1000) "Coil2 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil2(min = 10, start = 600) "Coil2 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    NuCoil1 := homotopy(0.135 * HTWallCorrFactorC1 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.135 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33) "VDI atlas. Original leading constant 0.084 and 0.14 of viscosity corr. exponent";
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if noEvent((T - Coil1.Ta) * (T - Coil1.Tb) <= 0) then
      LMTDcoil1 := 0;
    elseif noEvent(T - Coil1.Ta == T - Coil1.Tb) then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := homotopy((Coil1.Tb - Coil1.Ta) / log((T - Coil1.Ta) / (T - Coil1.Tb)), T - Coil1.Ta);
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.13);
    end if;
    NuCoil2 := homotopy(0.135 * HTWallCorrFactorC2 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.135 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta - (T - Coil2.Tb) == 0 then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := homotopy((T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb)), T - Coil2.Ta);
    end if;
  equation
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
//StateC2 = Medium.setBubbleState(Medium.setSat_T(Coil2.Tsurf));
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    if numCoils == 2 and numHalfCoils == 0 then
      W = (-Coil1.W) - Coil2.W;
    end if;
  end TankAgitAnchor2CoilsOld;

  model TankAgitAnchor2Coils1Hc "Anchor Vessel. Seems the best model for anchors"
    extends TankAgitAnchor2Coils(numCoils = 2, numHalfCoils = 1);
    replaceable VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 12}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e6) "wall jacket process side viscosity at wall temperature";
    Real AuxSurfEffHc1(start = 1) "fraction of internal auxiliar surface used in heat transfer";
    SI.NusseltNumber NuHalfCoil1(start = 1e2, min = 1, max = 1e5) "jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 5, start = 1000) "jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 5, start = 600) "jacket global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    Medium.ThermodynamicState StateHc1;
  algorithm
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if Mixer.Re < 100 then
      NuHalfCoil1 := homotopy(0.69 * HTWallCorrFactorHc1 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3), 0.69 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3));
    else
      NuHalfCoil1 := homotopy(0.56 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3), 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3));
    end if;
//NuHalfCoil1 := 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mu / MuWallHc1) ^ 0.14 "McKetta: Heat Transfer Design Methods. Original 0.46";
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    AuxSurfEffHc1 := (HhalfCoil1 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc1 := (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) - 1) / (AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi) * (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) + 1));
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if noEvent((T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0) then
      LMTDhalfCoil1 := 0;
    elseif noEvent(T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0) then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := homotopy((T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb)), T - HalfCoil1.Ta);
    end if;
  equation
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT);
//HalfCoil1.W = HhalfCoil1 * (T-HalfCoil1.Tsurf) * (HalfCoil1.SactiveHT+AuxSurfEffHc1*HalfCoil1.SauxHT);
    end if;
//StateHc1 = Medium.setBubbleState(Medium.setSat_T(HalfCoil1.Tsurf));
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if numCoils == 2 and numHalfCoils == 1 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W;
    end if;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end TankAgitAnchor2Coils1Hc;

  model TankAgitAnchor2Coils1HcOld "Anchor Vessel. Seems the best model for anchors"
    //Normally missing:N (rpm) from the mixer, and T from the process
    extends TankAgitatedAnchorPM(final numCoils = 2, final numHalfCoils = 1);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0) annotation(
      Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 16}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 12}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0), HTWallCorrFactorC2(start = 1.0);
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(start = 1e2, min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.NusseltNumber NuCoil2(start = 1e2, min = 1, max = 1e5) "Coil2 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Hcoil2(min = 10, start = 1000) "Coil2 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff., referenced to internal surface";
    SI.CoefficientOfHeatTransfer Ucoil2(min = 10, start = 600) "Coil2 global heat transfer coeff., referenced to internal surface";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Coil logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e6) "wall jacket process side viscosity at wall temperature";
    Real AuxSurfEffHc1(start = 1) "fraction of internal auxiliar surface used in heat transfer";
    SI.NusseltNumber NuHalfCoil1(start = 1e2, min = 1, max = 1e5) "jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 5, start = 1000) "jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 5, start = 600) "jacket global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2, StateHc1;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
//HTWallCorrFactorC1 := (Mu / MuWall1) ^ 0.13;
    NuCoil1 := homotopy(0.13 * HTWallCorrFactorC1 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.13 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
//NuCoil1 := 0.135 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33 * (Mu / MuWall1) ^ 0.12 "VDI atlas. Original 0.084 and 0.14 of exponent";
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if noEvent((T - Coil1.Ta) * (T - Coil1.Tb) <= 0) then
      LMTDcoil1 := 0;
    elseif noEvent(T - Coil1.Ta == T - Coil1.Tb) then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := homotopy((Coil1.Tb - Coil1.Ta) / log((T - Coil1.Ta) / (T - Coil1.Tb)), T - Coil1.Ta);
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.13);
    end if;
    NuCoil2 := homotopy(0.13 * HTWallCorrFactorC2 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.13 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta - (T - Coil2.Tb) == 0 then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := homotopy((T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb)), T - Coil2.Ta);
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if Mixer.Re < 100 then
      NuHalfCoil1 := homotopy(0.69 * HTWallCorrFactorHc1 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3), 0.69 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3));
    else
      NuHalfCoil1 := homotopy(0.56 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3), 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3));
    end if;
//NuHalfCoil1 := 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mu / MuWallHc1) ^ 0.14 "McKetta: Heat Transfer Design Methods. Original 0.46";
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    AuxSurfEffHc1 := (HhalfCoil1 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc1 := (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) - 1) / (AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi) * (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) + 1));
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if noEvent((T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0) then
      LMTDhalfCoil1 := 0;
    elseif noEvent(T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0) then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := homotopy((T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb)), T - HalfCoil1.Ta);
    end if;
//HalfCoil1.W := UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
  equation
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SactiveHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SactiveHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT);
//HalfCoil1.W = HhalfCoil1 * (T-HalfCoil1.Tsurf) * (HalfCoil1.SactiveHT+AuxSurfEffHc1*HalfCoil1.SauxHT);
    end if;
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
//StateC2 = Medium.setBubbleState(Medium.setSat_T(Coil2.Tsurf));
//StateHc1 = Medium.setBubbleState(Medium.setSat_T(HalfCoil1.Tsurf));
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if numCoils == 2 and numHalfCoils == 1 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W;
    end if;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end TankAgitAnchor2Coils1HcOld;

  model TankAgitAnchor2Coils2Hc
    extends TankAgitAnchor2Coils1Hc(numCoils = 2, numHalfCoils = 2);
    replaceable VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2 constrainedby FreeFluids.Pipes.PipeThermalBase(final useThermalConnector = false) annotation(
       Placement(visible = true, transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 12}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false, isBottomJacket = true) annotation(
      Placement(visible = true, transformation(origin = {0, -38}, extent = {{-12, 12}, {12, -12}}, rotation = 0)));
    Real HTWallCorrFactorHc2(start = 1.0);
    SI.DynamicViscosity MuWallHc2(min = 1e-6, start = 1e-3, max = 1e6) "halfcoil process side viscosity at wall temperature";
    Real AuxSurfEffHc2(start = 1) "fraction of internal auxiliar surface used in heat transfer";
    SI.NusseltNumber NuHalfCoil2(start = 1e2, min = 1, max = 1e5) "halfcoil process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil2(min = 5, start = 1000) "halfcoil process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil2(min = 5, start = 600) "halfcoil global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateHc2;
  algorithm
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHc2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    else
      HTWallCorrFactorHc2 := max(0.4, (Mu / MuWallHc2) ^ 0.14);
    end if;
    if Mixer.Re < 100 then
      NuHalfCoil2 := homotopy(0.69 * HTWallCorrFactorHc2 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3), 0.69 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3));
    else
      NuHalfCoil2 := homotopy(0.56 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3), 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3));
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    AuxSurfEffHc2 := (HhalfCoil2 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc2 := (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) - 1) / (AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi) * (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) + 1));
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
    if noEvent((T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0) then
      LMTDhalfCoil2 := 0;
    elseif noEvent(T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0) then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := homotopy((T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb)), T - HalfCoil2.Ta);
    end if;
  equation
    if HalfCoil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT);
//HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    end if;
    StateHc2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    if numCoils == 2 and numHalfCoils == 2 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W - HalfCoil2.W;
    end if;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end TankAgitAnchor2Coils2Hc;

  model TankAgitAnchor2Coils2HcOld "Anchor Vessel. Seems the best model for anchors"
    //Normally missing:N (rpm) from the mixer, and T from the process
    extends TankAgitatedAnchorPM(final numCoils = 2, final numHalfCoils = 2);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerAnchor Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0) annotation(
      Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 16}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false, isBottomJacket = true) annotation(
      Placement(visible = true, transformation(origin = {0, -38}, extent = {{-12, 12}, {12, -12}}, rotation = 0)));
    FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-40, 12}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0), HTWallCorrFactorC2(start = 1.0);
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(start = 1e2, min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.NusseltNumber NuCoil2(start = 1e2, min = 1, max = 1e5) "Coil2 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Hcoil2(min = 10, start = 1000) "Coil2 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff., referenced to internal surface";
    SI.CoefficientOfHeatTransfer Ucoil2(min = 10, start = 600) "Coil2 global heat transfer coeff., referenced to internal surface";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Coil logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e6) "halfcoil process side viscosity at wall temperature";
    Real AuxSurfEffHc1(start = 1) "fraction of internal auxiliar surface used in heat transfer";
    SI.NusseltNumber NuHalfCoil1(start = 1e2, min = 1, max = 1e5) "halfcoil process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 5, start = 1000) "halfcoil process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 5, start = 600) "halfcoil global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    Real HTWallCorrFactorHc2(start = 1.0);
    SI.DynamicViscosity MuWallHc2(min = 1e-6, start = 1e-3, max = 1e6) "halfcoil process side viscosity at wall temperature";
    Real AuxSurfEffHc2(start = 1) "fraction of internal auxiliar surface used in heat transfer";
    SI.NusseltNumber NuHalfCoil2(start = 1e2, min = 1, max = 1e5) "halfcoil process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil2(min = 5, start = 1000) "halfcoil process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil2(min = 5, start = 600) "halfcoil global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2, StateHc1, StateHc2;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    NuCoil1 := homotopy(0.13 * HTWallCorrFactorC1 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.13 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
//NuCoil1 := 0.135 * (Coil1.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33 * (Mu / MuWall1) ^ 0.12 "VDI atlas. Original 0.084 and 0.14 of exponent";
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if noEvent((T - Coil1.Ta) * (T - Coil1.Tb) <= 0) then
      LMTDcoil1 := 0;
    elseif noEvent(T - Coil1.Ta == T - Coil1.Tb) then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := homotopy((Coil1.Tb - Coil1.Ta) / log((T - Coil1.Ta) / (T - Coil1.Tb)), T - Coil1.Ta);
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.14);
    end if;
    NuCoil2 := homotopy(0.13 * HTWallCorrFactorC2 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33, 0.13 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33);
//NuCoil2 := 0.135 * (Coil2.coilDiam / Vessel.di) ^ 0.74 * Mixer.Re ^ 0.55 * Pr ^ 0.33 * (Mu / MuWall2) ^ 0.13;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta - (T - Coil2.Tb) == 0 then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := homotopy((T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb)), T - Coil2.Ta);
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if Mixer.Re < 100 then
      NuHalfCoil1 := homotopy(0.69 * HTWallCorrFactorHc1 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3), 0.69 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3));
    else
      NuHalfCoil1 := homotopy(0.56 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3), 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3)) "McKetta: Heat Transfer Design Methods. Original 0.46";
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    AuxSurfEffHc1 := (HhalfCoil1 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc1 := (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) - 1) / (AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi) * (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) + 1));
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if noEvent((T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0) then
      LMTDhalfCoil1 := 0;
    elseif noEvent(T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0) then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := homotopy((T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb)), T - HalfCoil1.Ta);
    end if;
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHc2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    else
      HTWallCorrFactorHc2 := max(0.4, (Mu / MuWallHc2) ^ 0.14);
    end if;
    if Mixer.Re < 100 then
      NuHalfCoil2 := homotopy(0.69 * HTWallCorrFactorHc2 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3), 0.69 * Mixer.Re ^ 0.5 * Pr ^ (1 / 3));
    else
      NuHalfCoil2 := homotopy(0.56 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3), 0.56 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3));
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    AuxSurfEffHc2 := (HhalfCoil2 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc2 := (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) - 1) / (AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi) * (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) + 1));
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
    if noEvent((T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0) then
      LMTDhalfCoil2 := 0;
    elseif noEvent(T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0) then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := homotopy((T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb)), T - HalfCoil2.Ta);
    end if;
  equation
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SactiveHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SactiveHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT);
//HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if HalfCoil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT);
//HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    end if;
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    StateHc2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    if numCoils == 2 and numHalfCoils == 2 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W - HalfCoil2.W;
    end if;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end TankAgitAnchor2Coils2HcOld;

  partial model TankAgitatedKameiPM
    extends TankPM;
    replaceable MixerKamei Mixer;
  equation
//flow numbers must be reviewed
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" then
      Mixer.NpMax = 6.5 * (Mixer.nBlades ^ 0.7 * Mixer.width * sin(Mixer.angle) ^ 1.6 / Mixer.d) ^ 1.7;
      Mixer.Nq = 0.62 * (Mixer.d / Vessel.di) ^ (-0.7) * (Mixer.width / Vessel.di) ^ 0.82 * Mixer.nBlades ^ 0.6;
    elseif Mixer.reference == "radialTurbine" or Mixer.reference == "anchor" then
      if Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d <= 0.54 then
        Mixer.NpMax = 10 * (Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d) ^ 1.3;
      elseif Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d > 0.54 and Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d <= 1.6 then
        Mixer.NpMax = 8.3 * Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d;
      else
        Mixer.NpMax = 10 * (Mixer.nBlades ^ 0.7 * Mixer.width / Mixer.d) ^ 0.6;
      end if;
      Mixer.Nq = 1.3 * (Mixer.d / Vessel.di) ^ (-0.86) * (Mixer.width / Vessel.di) ^ 0.82 * Mixer.nBlades ^ 0.6;
    elseif Mixer.reference == "intermig" then
      Mixer.NpMax = 0.3;
      Mixer.Nq = 0.8 * (Mixer.d / Vessel.di) ^ (-0.7) * (0.2 * Mixer.d / Vessel.di) ^ 0.82 * Mixer.nBlades ^ 0.6;
    elseif Mixer.reference == "mig" then
      Mixer.NpMax = 0.2;
      Mixer.Nq = 0.8 * (Mixer.d / Vessel.di) ^ (-0.7) * (0.2 * Mixer.d / Vessel.di) ^ 0.82 * Mixer.nBlades ^ 0.6;
    else
      Mixer.NpMax = 8.3 * (2 * Mixer.angle / pi) ^ 0.9 * Mixer.nBlades ^ 0.7 * Mixer.width * sin(Mixer.angle) ^ 1.6 / Mixer.d "axial turbines";
      Mixer.Nq = 0.8 * (Mixer.d / Vessel.di) ^ (-0.7) * (Mixer.width * sin(Mixer.angle) / Vessel.di) ^ 0.82 * Mixer.nBlades ^ 0.6;
    end if;
//cambio explicado correcto? Va el ^1.7?
    Mixer.Re = Mixer.d ^ 2 * Mixer.N * Rho / Mu;
    if Mixer.reference == "intermig" or Mixer.reference == "mig" then
      Mixer.eta = 0;
      Mixer.beta = 0;
      Mixer.gamma = 0;
      Mixer.X = 0;
      Mixer.Ct = 0;
      Mixer.m = 0;
      Mixer.fInf = 0;
      Mixer.Ctr = 0;
      Mixer.Cl = 0;
      Mixer.ReG = 0;
      Mixer.f = 0;
      if Mixer.reference == "intermig" then
        if Mixer.Re < 500 then
          Mixer.Np0 = 0.5 * (0.39471382 + 100.84884 / Mixer.Re + 44.730788 / Mixer.Re ^ 1.5 - 45.978678 / Mixer.Re ^ 2);
        elseif Mixer.Re < 10000 then
          Mixer.Np0 = 0.5 * (0.25174729 - 3.8632725E-14 * Mixer.Re ^ 3 + 203533.57 * log(Mixer.Re) / Mixer.Re ^ 2 - 1177817 / Mixer.Re ^ 2);
        else
          Mixer.Np0 = 0.1;
        end if;
        if Mixer.Re < 10000 then
          Mixer.Np = 0.5 * (0.57957364 + 22.206363 * log(Mixer.Re) / Mixer.Re + 89.438236 / Mixer.Re ^ 1.5 + 56.420189 * log(Mixer.Re) / Mixer.Re ^ 2);
        else
          Mixer.Np = 0.5 * 0.6;
        end if;
      elseif Mixer.reference == "mig" then
        if Mixer.Re < 500 then
          Mixer.Np0 = (0.39471382 + 100.84884 / Mixer.Re + 44.730788 / Mixer.Re ^ 1.5 - 45.978678 / Mixer.Re ^ 2) / 3;
        elseif Mixer.Re < 10000 then
          Mixer.Np0 = (0.25174729 - 3.8632725E-14 * Mixer.Re ^ 3 + 203533.57 * log(Mixer.Re) / Mixer.Re ^ 2 - 1177817 / Mixer.Re ^ 2) / 3;
        else
          Mixer.Np0 = 0.067;
        end if;
        if Mixer.Re < 4000 then
          Mixer.Np = (0.56800494 - 8.1033728 * log(Mixer.Re) / Mixer.Re + 94.680892 / Mixer.Re - 15.271562 / Mixer.Re ^ 2) / 3;
        else
          Mixer.Np = 0.58 / 3;
        end if;
      end if;
      Mixer.NpB = 0;
    else
      Mixer.eta = 0.711 * (0.157 + (Mixer.nBlades * log(Vessel.di / Mixer.d)) ^ 0.611) / (Mixer.nBlades ^ 0.52 * (1 - (Mixer.d / Vessel.di) ^ 2)) "OK";
      Mixer.beta = 2 * log(Vessel.di / Mixer.d) / (Vessel.di / Mixer.d - Mixer.d / Vessel.di) "OK";
      Mixer.gamma = (Mixer.eta * log(Vessel.di / Mixer.d) / (Mixer.beta * Vessel.di / Mixer.d) ^ 5) ^ (1 / 3) "OK";
      Mixer.X = Mixer.gamma * Mixer.nBlades ^ 0.7 * Mixer.width * sin(Mixer.angle) ^ 1.6 / hLiquid "OK";
      if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" then
        Mixer.Ct = ((3 * Mixer.X ^ 1.5) ^ (-7.8) + 0.25 ^ (-7.8)) ^ (-1 / 7.8);
        Mixer.m = ((0.8 * Mixer.X ^ 0.373) ^ (-7.8) + 0.333 ^ (-7.8)) ^ (-1 / 7.8);
      else
        Mixer.Ct = ((1.96 * Mixer.X ^ 1.19) ^ (-7.8) + 0.25 ^ (-7.8)) ^ (-1 / 7.8);
        Mixer.m = ((0.71 * Mixer.X ^ 0.373) ^ (-7.8) + 0.333 ^ (-7.8)) ^ (-1 / 7.8);
      end if;
      Mixer.fInf = 0.0151 * Mixer.d / Vessel.di * Mixer.Ct ^ 0.308 "it was 0.00756 in original Kamei for axial turbine";
      Mixer.Ctr = 23.8 * (Mixer.d / Vessel.di) ^ (-3.24) * (Mixer.width * sin(Mixer.angle) / Vessel.di) ^ (-1.18) * Mixer.X ^ (-0.74) "OK";
      Mixer.Cl = 0.215 * Mixer.eta * Mixer.nBlades * Mixer.d / hLiquid * (1 - (Mixer.d / Vessel.di) ^ 2) + 1.83 * Mixer.width * sin(Mixer.angle) / hLiquid * (Mixer.nBlades / sin(Mixer.angle) / 2) ^ (1 / 3) "?";
      Mixer.ReG = pi * Mixer.eta * log(Vessel.di / Mixer.d) / (4 * Mixer.d / (Mixer.beta * Vessel.di)) * Mixer.Re "OK";
      Mixer.f = Mixer.Cl / Mixer.ReG + Mixer.Ct * ((Mixer.Ctr / Mixer.ReG + Mixer.ReG) ^ (-1) + (Mixer.fInf / Mixer.Ct) ^ (1 / Mixer.m)) ^ Mixer.m "friction factor OK";
      Mixer.Np0 = 1.2 * pi ^ 4 * Mixer.beta ^ 2 / (8 * Mixer.d ^ 3 / (Vessel.di ^ 2 * hLiquid)) * Mixer.f "OK";
      Mixer.NpB = (1 + (4.5 * Vessel.baffleWidth / Vessel.di * Vessel.nBaffles ^ 0.8 / ((2 * Mixer.angle / pi) ^ 0.72 * Mixer.NpMax ^ 0.2) + Mixer.Np0 / Mixer.NpMax) ^ (-3)) ^ (-1 / 3) * Mixer.NpMax "OK";
      if Mixer.NpB > Mixer.Np0 then
        Mixer.Np = Mixer.NpB;
      else
        Mixer.Np = Mixer.Np0;
      end if;
    end if;
    Mixer.W = Mixer.nImpellers * Mixer.Np * Rho * Mixer.d ^ 5 * Mixer.N ^ 3;
    Mixer.Htop = Vessel.HtotalIn - Mixer.hBottom;
    Mixer.Q = Mixer.nImpellers * Mixer.Nq * Mixer.N * Mixer.d ^ 3;
    Mixer.Vc = 4 * Mixer.Q / (pi * (Vessel.di ^ 2 * hLiquid) ^ (2 / 3));
    Mixer.ScaleOfAgitation = Mixer.Vc * 32.8;
    annotation(
      Icon(graphics = {Ellipse(origin = {-13, -57}, extent = {{13, -3}, {-15, 5}}, endAngle = 360), Ellipse(origin = {12, -55}, extent = {{-12, 3}, {16, -5}}, endAngle = 360)}),
      Documentation(info = "<html><head></head><body>From Furukawa et alt., International Journal of Chemical Engineering volume 2012, based on Kamei et al.(1995/6/7) model for many different impeller types</body></html>"));
  end TankAgitatedKameiPM;

  //***AGITATED TANK KAMEI MODELS***
  //--------------------------------

  model TankAgitatedKamei
    extends TankAgitatedKameiPM;
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end TankAgitatedKamei;

  model TankAgitKamei1Coil
    extends FreeFluids.Vessels.TankAgitatedKameiPM(numCoils = 1, numHalfCoils = 0);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if numCoils == 1 and numHalfCoils == 0 then
      W = -Coil1.W;
    end if;
  end TankAgitKamei1Coil;

  model TankAgitKamei1Coil1Hc
    extends FreeFluids.Vessels.TankAgitKamei1Coil(numCoils = 1, numHalfCoils = 1);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateHc1;
  algorithm
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
  equation
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if numCoils == 1 and numHalfCoils == 1 then
      W = (-Coil1.W) - HalfCoil1.W;
    end if;
  end TankAgitKamei1Coil1Hc;

  model TankAgitKamei1coils1hc
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 1, final numHalfCoils = 1);
    FreeFluids.Vessels.VesselCylVert Vessel(final nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1 "Logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateHc1;
  algorithm
//Nusselt halfcoil must be reviewed for mig and intermig
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
//StateHc1 = Medium.setBubbleState(Medium.setSat_T(HalfCoil1.Tsurf));
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if numCoils == 1 and numHalfCoils == 1 then
      W = (-Coil1.W) - HalfCoil1.W;
    end if;
  end TankAgitKamei1coils1hc;

  model TankAgitKamei1Coil2Hc
    extends TankAgitKamei1Coil1Hc(numCoils = 1, numHalfCoils = 2);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false, PLossFriction(displayUnit = "bar"), basePipeAngle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorHc2(start = 1.0);
    SI.DynamicViscosity MuWallHc2(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil2(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil2(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil2(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil2 "Logarithmic mean temperature difference";
    Medium.ThermodynamicState StateHc2;
  algorithm
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHc2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    else
      HTWallCorrFactorHc2 := max(0.4, (Mu / MuWallHc2) ^ 0.14);
    end if;
    if HalfCoil2.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil2 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil2 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil2 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        else
          NuHalfCoil2 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
    if (T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0 then
      LMTDhalfCoil2 := 0;
    elseif T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0 then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := (T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb));
    end if;
  equation
    StateHc2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    if HalfCoil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    end if;
    if numCoils == 1 and numHalfCoils == 2 then
      W = (-Coil1.W) - HalfCoil1.W - HalfCoil2.W;
    end if;
  end TankAgitKamei1Coil2Hc;

  model TankAgitKamei1coil2hc
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 1, final numHalfCoils = 2);
    FreeFluids.Vessels.VesselCylVert Vessel(final nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, PLossFriction(displayUnit = "bar"), fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false, isBottomJacket = true, PLossFriction(displayUnit = "bar"), basePipeAngle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false, PLossFriction(displayUnit = "bar"), basePipeAngle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1 "Logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0), HTWallCorrFactorHc2(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3), MuWallHc2(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5), NuHalfCoil2(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000), HhalfCoil2(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600), UhalfCoil2(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1, LMTDhalfCoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateHc1, StateHc2;
  algorithm
//Nusselt halfcoil must be reviewed for mig and intermig
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHc2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    else
      HTWallCorrFactorHc2 := max(0.4, (Mu / MuWallHc2) ^ 0.14);
    end if;
    if HalfCoil2.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil2 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil2 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil2 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        else
          NuHalfCoil2 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
    if (T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0 then
      LMTDhalfCoil2 := 0;
    elseif T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0 then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := (T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb));
    end if;
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    StateHc2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if HalfCoil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    end if;
    if numCoils == 1 and numHalfCoils == 2 then
      W = (-Coil1.W) - HalfCoil1.W - HalfCoil2.W;
    end if;
  end TankAgitKamei1coil2hc;

  model TankAgitKamei2Coils
    extends FreeFluids.Vessels.TankAgitKamei1Coil(numCoils = 2, numHalfCoils = 0);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 24}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Real HTWallCorrFactorC2(start = 1.0) "Mu/Mu Wall for coil 2";
    SI.DynamicViscosity MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil 2 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil2(min = 1, max = 1e5) "Coil 2 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil2(min = 10, start = 1000) "Coil 2 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil2(min = 10, start = 600) "Coil 2 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil2 "Logarithmic mean temperature difference";
    Medium.ThermodynamicState StateC2;
  algorithm
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    elseif noEvent(Mu / MuWall2 >= 1.0e4) then
      HTWallCorrFactorC2 := 10000.0;
    elseif noEvent(Mu / MuWall2 <= 0.000006) then
      HTWallCorrFactorC2 := 0.000006;
    else
      HTWallCorrFactorC2 := Mu / MuWall2;
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil2 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil2 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil2 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    end if;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta == T - Coil2.Tb then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := (T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb));
    end if;
  equation
    if Coil2.fullHTlength == false then
      if Coil2.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = 0.0;
      elseif Coil2.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = Coil2.SactiveHT;
      else
        Coil2.SusedHT = Coil2.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil2.heightInit) / Coil2.CoilHeigth;
      end if;
    end if;
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    if numCoils == 2 and numHalfCoils == 0 then
      W = (-Coil1.W) - Coil2.W;
    end if;
  end TankAgitKamei2Coils;

  model TankAgitKamei2coils
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 2, final numHalfCoils = 0);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 24}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    Real HTWallCorrFactorC2(start = 1.0) "Mu/Mu Wall for coil 2";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5), NuCoil2(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000), Hcoil2(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600), Ucoil2(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2;
  algorithm
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    elseif noEvent(Mu / MuWall1 >= 1.0e4) then
      HTWallCorrFactorC1 := 10000.0;
    elseif noEvent(Mu / MuWall1 <= 0.000006) then
      HTWallCorrFactorC1 := 0.000006;
    else
      HTWallCorrFactorC1 := Mu / MuWall1;
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    elseif noEvent(Mu / MuWall2 >= 1.0e4) then
      HTWallCorrFactorC2 := 10000.0;
    elseif noEvent(Mu / MuWall2 <= 0.000006) then
      HTWallCorrFactorC2 := 0.000006;
    else
      HTWallCorrFactorC2 := Mu / MuWall2;
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil2 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil2 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil2 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    end if;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta == T - Coil2.Tb then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := (T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb));
    end if;
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
    if Coil2.fullHTlength == false then
      if Coil2.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = 0.0;
      elseif Coil2.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = Coil2.SactiveHT;
      else
        Coil2.SusedHT = Coil2.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil2.heightInit) / Coil2.CoilHeigth;
      end if;
    end if;
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    if numCoils == 2 and numHalfCoils == 0 then
      W = (-Coil1.W) - Coil2.W;
    end if;
  end TankAgitKamei2coils;

  model TankAgitKamei2Coils1Hc
    extends FreeFluids.Vessels.TankAgitKamei2Coils(numCoils = 2, numHalfCoils = 1);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 2) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateHc1;
  algorithm
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.14);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
  equation
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if numCoils == 2 and numHalfCoils == 1 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W;
    end if;
  end TankAgitKamei2Coils1Hc;

  model TankAgitKamei2coils1hc
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 2, final numHalfCoils = 1);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 24}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    Real HTWallCorrFactorC2(start = 1.0) "Mu/Mu Wall for coil 2";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5), NuCoil2(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000), Hcoil2(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600), Ucoil2(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2, StateHc1;
  algorithm
//Nusselt halfcoil must be reviewed for mig and intermig
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    elseif noEvent(Mu / MuWall1 >= 1.0e4) then
      HTWallCorrFactorC1 := 10000.0;
    elseif noEvent(Mu / MuWall1 <= 0.000006) then
      HTWallCorrFactorC1 := 0.000006;
    else
      HTWallCorrFactorC1 := Mu / MuWall1;
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    elseif noEvent(Mu / MuWall2 >= 1.0e4) then
      HTWallCorrFactorC2 := 10000.0;
    elseif noEvent(Mu / MuWall2 <= 0.000006) then
      HTWallCorrFactorC2 := 0.000006;
    else
      HTWallCorrFactorC2 := Mu / MuWall2;
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil2 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil2 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil2 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    end if;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta == T - Coil2.Tb then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := (T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb));
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    elseif noEvent(Mu / MuWallHc1 >= 1e4) then
      HTWallCorrFactorHc1 := 10000 ^ 0.14;
    elseif noEvent(Mu / MuWallHc1 <= 0.000008) then
      HTWallCorrFactorHc1 := 0.000008 ^ 0.14;
    else
      HTWallCorrFactorHc1 := (Mu / MuWallHc1) ^ 0.14;
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
//HalfCoilW.W := UhalfCoilW * (T - HalfCoilW.Ta + T - HalfCoilW.Tb) / 2 * HalfCoilW.SactivedHT;
//HalfCoilW.W := UhalfCoilW * LMTDhalfCoilW * HalfCoilW.SactivedHT;
/*elseif Mixer.reference == "radialTurbine" then
        NuJacketC1 = 0.66 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / Vessel.hLiquid*Mixer.nImpellers) ^ 0.15 * (Mixer.nBlades * Mixer.width / (4 * 0.2    * Vessel.di)) ^ 0.2 * (Mu / MuWallJC1) ^ 0.14 "Handbook of Industrial Mixing";
        NuJacketC2 = 0.66 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / Vessel.hLiquid*Mixer.nImpellers) ^ 0.15 * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * (Mu / MuWallJC2) ^ 0.14 "Handbook of Industrial Mixing";
      elseif Mixer.reference == "axialTurbine" then
        NuJacketC1 = 0.45 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / Vessel.hLiquid*Mixer.nImpellers) ^ 0.15 * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * (Mu / MuWallJC1) ^ 0.14 "Handbook of Industrial Mixing";
        NuJacketC2 = 0.45 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / Vessel.hLiquid*Mixer.nImpellers) ^ 0.15 * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * (Mu / MuWallJC2) ^ 0.14 "Handbook of Industrial Mixing";
      end if;*/
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
    if Coil2.fullHTlength == false then
      if Coil2.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = 0.0;
      elseif Coil2.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = Coil2.SactiveHT;
      else
        Coil2.SusedHT = Coil2.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil2.heightInit) / Coil2.CoilHeigth;
      end if;
    end if;
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if numCoils == 2 and numHalfCoils == 1 then
      W = (-Coil1.W) - HalfCoil1.W - Coil2.W;
    end if;
  end TankAgitKamei2coils1hc;

  model TankAgitKamei2coil2hc
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 2, final numHalfCoils = 2);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer(angle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {0, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 24}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil1(final useThermalConnector = false, isBottomJacket = true, basePipeAngle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false, basePipeAngle(displayUnit = "deg")) annotation(
      Placement(visible = true, transformation(origin = {-34, 4}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    Real HTWallCorrFactorC2(start = 1.0) "Mu/Mu Wall for coil 2";
    SI.DynamicViscosity MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil1(min = 1, max = 1e5), NuCoil2(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil1(min = 10, start = 1000), Hcoil2(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil1(min = 10, start = 600), Ucoil2(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil1, LMTDcoil2 "Logarithmic mean temperature difference";
    Real HTWallCorrFactorHc1(start = 1.0), HTWallCorrFactorHc2(start = 1.0);
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3), MuWallHc2(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    SI.NusseltNumber NuHalfCoil1(min = 1, max = 1e5), NuHalfCoil2(min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000), HhalfCoil2(min = 10, start = 1000) "wall jacket process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600), UhalfCoil2(min = 10, start = 600) "wall jacket global heat transfer coeff.";
    SI.TemperatureDifference LMTDhalfCoil1, LMTDhalfCoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC1, StateC2, StateHc1, StateHc2;
  algorithm
//Nusselt halfcoil must be reviewed for mig and intermig
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil2 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil2 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil2 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    end if;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta == T - Coil2.Tb then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := (T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb));
    end if;
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHc1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.13);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc1 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc1 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc1 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
    if (T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0 then
      LMTDhalfCoil1 := 0;
    elseif T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0 then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHc2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    elseif noEvent(Mu / MuWallHc2 >= 1e4) then
      HTWallCorrFactorHc2 := 10000 ^ 0.14;
    elseif noEvent(Mu / MuWallHc2 <= 0.000008) then
      HTWallCorrFactorHc2 := 0.000008 ^ 0.14;
    else
      HTWallCorrFactorHc2 := (Mu / MuWallHc2) ^ 0.14;
    end if;
    if HalfCoil2.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.9 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 1.1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * HTWallCorrFactorHc2 "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil2 := 1.08 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 * HTWallCorrFactorHc2 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.31 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 0.5 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * HTWallCorrFactorHc2 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil2 := 0.74 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil2 := 0.54 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        else
          NuHalfCoil2 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
    if (T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0 then
      LMTDhalfCoil2 := 0;
    elseif T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0 then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := (T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb));
    end if;
  equation
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
    if Coil2.fullHTlength == false then
      if Coil2.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = 0.0;
      elseif Coil2.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = Coil2.SactiveHT;
      else
        Coil2.SusedHT = Coil2.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil2.heightInit) / Coil2.CoilHeigth;
      end if;
    end if;
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    StateHc1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    StateHc2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    if HalfCoil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    end if;
    if HalfCoil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    end if;
    if numCoils == 2 and numHalfCoils == 2 then
      W = (-Coil1.W) - Coil2.W - HalfCoil1.W - HalfCoil2.W;
    end if;
  end TankAgitKamei2coil2hc;

  model TankAgitKamei3Coils
    extends FreeFluids.Vessels.TankAgitKamei2Coils(numCoils = 3, numHalfCoils = 0);
    replaceable FreeFluids.Vessels.VesselCylVert Vessel(nConcCoils = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    replaceable FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 24}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil0(final thicknessInsul = 0, final useThermalConnector = false, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real HTWallCorrFactorC0(start = 1.0) "Mu/Mu Wall for coil 0";
    SI.DynamicViscosity MuWall0(min = 1e-6, start = 1e-3, max = 1e7) "Coil 0 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil0(min = 1, max = 1e5) "Coil 0 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil0(min = 10, start = 1000) "Coil 0 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil0(min = 10, start = 600) "Coil 0 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil0 "Logarithmic mean temperature difference";
    Medium.ThermodynamicState StateC0;
  algorithm
    MuWall0 := Medium.dynamicViscosity(StateC0);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC0 := 1.0;
    else
      HTWallCorrFactorC0 := max(0.4, (Mu / MuWall0) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil0 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil0 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil0 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.08;
    end if;
    Hcoil0 := NuCoil0 * K / Coil0.Do;
    Ucoil0 := 1 / ((foulingF + 1 / Hcoil0) * Coil0.Di / Coil0.Do + Coil0.Do * log(Coil0.Do / Coil0.Di) / (2 * Vessel.wallK) + 1 / Coil0.H + Coil0.foulingF);
    if (T - Coil0.Ta) * (T - Coil0.Tb) <= 0 then
      LMTDcoil0 := 0;
    elseif T - Coil0.Ta == T - Coil0.Tb then
      LMTDcoil0 := T - Coil0.Ta;
    else
      LMTDcoil0 := (T - Coil0.Ta - (T - Coil0.Tb)) / log((T - Coil0.Ta) / (T - Coil0.Tb));
    end if;
  equation
    if Coil0.fullHTlength == false then
      if Coil0.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil0.SusedHT = 0.0;
      elseif Coil0.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil0.SusedHT = Coil0.SactiveHT;
      else
        Coil0.SusedHT = Coil0.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil0.heightInit) / Coil0.CoilHeigth;
      end if;
    end if;
    StateC0 = Medium.setState_pTX(fixedPressure, Coil0.Tsurf);
    if Coil0.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil0.W = Ucoil0 * LMTDcoil0 * Coil0.SusedHT;
    end if;
    if numCoils == 3 then
      W = (-Coil0.W) - Coil1.W - Coil2.W;
    end if;
  end TankAgitKamei3Coils;

  model TankAgitKamei3coils
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 3, final numHalfCoils = 0);
    FreeFluids.Vessels.VesselCylVert Vessel annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil0(final thicknessInsul = 0, final useThermalConnector = false, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil1(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.CoilForcedConvection Coil2(final useThermalConnector = false, final thicknessInsul = 0, fullHTlength = false) annotation(
      Placement(visible = true, transformation(origin = {-8.88178e-16, 26}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Real HTWallCorrFactorC0(start = 1.0) "Mu/Mu Wall for coil 0";
    Real HTWallCorrFactorC1(start = 1.0) "Mu/Mu Wall for coil 1";
    Real HTWallCorrFactorC2(start = 1.0) "Mu/Mu Wall for coil 2";
    SI.DynamicViscosity MuWall0(min = 1e-6, start = 1e-3, max = 1e7), MuWall1(min = 1e-6, start = 1e-3, max = 1e7), MuWall2(min = 1e-6, start = 1e-3, max = 1e7) "Coil1 process side viscosity at wall temperature";
    SI.NusseltNumber NuCoil0(min = 1, max = 1e5), NuCoil1(min = 1, max = 1e5), NuCoil2(min = 1, max = 1e5) "Coil1 process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer Hcoil0(min = 10, start = 1000), Hcoil1(min = 10, start = 1000), Hcoil2(min = 10, start = 1000) "Coil1 process side heat transfer coeff.";
    SI.CoefficientOfHeatTransfer Ucoil0(min = 10, start = 600), Ucoil1(min = 10, start = 600), Ucoil2(min = 10, start = 600) "Coil1 global heat transfer coeff.";
    SI.TemperatureDifference LMTDcoil0, LMTDcoil1, LMTDcoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateC0, StateC1, StateC2;
  algorithm
    MuWall0 := Medium.dynamicViscosity(StateC0);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC0 := 1.0;
    else
      HTWallCorrFactorC0 := max(0.4, (Mu / MuWall0) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil0 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil0 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil0 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil0.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC0 ^ 0.08;
    end if;
    Hcoil0 := NuCoil0 * K / Coil0.Do;
    Ucoil0 := 1 / ((foulingF + 1 / Hcoil0) * Coil0.Di / Coil0.Do + Coil0.Do * log(Coil0.Do / Coil0.Di) / (2 * Vessel.wallK) + 1 / Coil0.H + Coil0.foulingF);
    if (T - Coil0.Ta) * (T - Coil0.Tb) <= 0 then
      LMTDcoil0 := 0;
    elseif T - Coil0.Ta == T - Coil0.Tb then
      LMTDcoil0 := T - Coil0.Ta;
    else
      LMTDcoil0 := (T - Coil0.Ta - (T - Coil0.Tb)) / log((T - Coil0.Ta) / (T - Coil0.Tb));
    end if;
    MuWall1 := Medium.dynamicViscosity(StateC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC1 := 1.0;
    else
      HTWallCorrFactorC1 := max(0.4, (Mu / MuWall1) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil1 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil1 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil1 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil1.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC1 ^ 0.08;
    end if;
    Hcoil1 := NuCoil1 * K / Coil1.Do;
    Ucoil1 := 1 / ((foulingF + 1 / Hcoil1) * Coil1.Di / Coil1.Do + Coil1.Do * log(Coil1.Do / Coil1.Di) / (2 * Vessel.wallK) + 1 / Coil1.H + Coil1.foulingF);
    if (T - Coil1.Ta) * (T - Coil1.Tb) <= 0 then
      LMTDcoil1 := 0;
    elseif T - Coil1.Ta == T - Coil1.Tb then
      LMTDcoil1 := T - Coil1.Ta;
    else
      LMTDcoil1 := (T - Coil1.Ta - (T - Coil1.Tb)) / log((T - Coil1.Ta) / (T - Coil1.Tb));
    end if;
    MuWall2 := Medium.dynamicViscosity(StateC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorC2 := 1.0;
    else
      HTWallCorrFactorC2 := max(0.4, (Mu / MuWall2) ^ 0.13);
    end if;
    if Mixer.reference == "hydrofoil" or Mixer.reference == "propeller" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
      NuCoil2 := 0.016 * Mixer.Re ^ 0.67 * Pr ^ 0.37 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.14;
    elseif Mixer.reference == "radialTurbine" then
      NuCoil2 := 0.03 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width / (0.2 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    elseif Mixer.reference == "axialTurbine" then
      NuCoil2 := 0.025 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (Mixer.width * sin(Mixer.angle) / (0.17 * Mixer.d)) ^ 0.2 * (3 * Mixer.d / Vessel.di) ^ 0.1 * (Coil2.Do / Vessel.di / 0.04) ^ 0.5 * (2 / Mixer.nBlades) ^ 0.2 * 0.82 ^ (Vessel.nConcCoils - 1) * HTWallCorrFactorC2 ^ 0.08;
    end if;
    Hcoil2 := NuCoil2 * K / Coil2.Do;
    Ucoil2 := 1 / ((foulingF + 1 / Hcoil2) * Coil2.Di / Coil2.Do + Coil2.Do * log(Coil2.Do / Coil2.Di) / (2 * Vessel.wallK) + 1 / Coil2.H + Coil2.foulingF);
    if (T - Coil2.Ta) * (T - Coil2.Tb) <= 0 then
      LMTDcoil2 := 0;
    elseif T - Coil2.Ta == T - Coil2.Tb then
      LMTDcoil2 := T - Coil2.Ta;
    else
      LMTDcoil2 := (T - Coil2.Ta - (T - Coil2.Tb)) / log((T - Coil2.Ta) / (T - Coil2.Tb));
    end if;
  equation
    if Coil0.fullHTlength == false then
      if Coil0.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil0.SusedHT = 0.0;
      elseif Coil0.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil0.SusedHT = Coil0.SactiveHT;
      else
        Coil0.SusedHT = Coil0.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil0.heightInit) / Coil0.CoilHeigth;
      end if;
    end if;
    if Coil1.fullHTlength == false then
      if Coil1.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = 0.0;
      elseif Coil1.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil1.SusedHT = Coil1.SactiveHT;
      else
        Coil1.SusedHT = Coil1.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil1.heightInit) / Coil1.CoilHeigth;
      end if;
    end if;
    if Coil2.fullHTlength == false then
      if Coil2.heightInit >= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = 0.0;
      elseif Coil2.CoilFinalHeight <= hLiquid - Vessel.HbottomIn then
        Coil2.SusedHT = Coil2.SactiveHT;
      else
        Coil2.SusedHT = Coil2.SactiveHT * (hLiquid - Vessel.HbottomIn - Coil2.heightInit) / Coil2.CoilHeigth;
      end if;
    end if;
//StateC0 = Medium.setBubbleState(Medium.setSat_T(Coil0.Tsurf));
//StateC1 = Medium.setBubbleState(Medium.setSat_T(Coil1.Tsurf));
//StateC2 = Medium.setBubbleState(Medium.setSat_T(Coil2.Tsurf));
    StateC0 = Medium.setState_pTX(fixedPressure, Coil0.Tsurf);
    StateC1 = Medium.setState_pTX(fixedPressure, Coil1.Tsurf);
    StateC2 = Medium.setState_pTX(fixedPressure, Coil2.Tsurf);
    if Coil0.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil0.W = Ucoil0 * LMTDcoil0 * Coil0.SusedHT;
    end if;
    if Coil1.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil1.W = Ucoil1 * LMTDcoil1 * Coil1.SusedHT;
    end if;
    if Coil2.thermalType == FreeFluids.Types.ThermalType.detailed then
      Coil2.W = Ucoil2 * LMTDcoil2 * Coil2.SusedHT;
    end if;
    if numCoils == 3 then
      W = (-Coil0.W) - Coil1.W - Coil2.W;
    end if;
  end TankAgitKamei3coils;

  model TankAgitKamei2hc
    extends FreeFluids.Vessels.TankAgitatedKameiPM(final numCoils = 0, final numHalfCoils = 2);
    FreeFluids.Vessels.VesselCylVert Vessel(final nConcCoils = 0) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    FreeFluids.Vessels.MixerKamei Mixer annotation(
      Placement(visible = true, transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilCondensing HalfCoil1(final useThermalConnector = false, isBottomJacket = true) constrainedby FreeFluids.Pipes.PipeThermalBase annotation(
       Placement(visible = true, transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil2(final useThermalConnector = false, isBottomJacket = false) constrainedby FreeFluids.Pipes.PipeThermalBase annotation(
       Placement(visible = true, transformation(origin = {-34, 12}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    SI.DynamicViscosity MuWallHc1(min = 1e-6, start = 1e-3, max = 1e3), MuWallHc2(min = 1e-6, start = 1e-3, max = 1e3) "wall jacket process side viscosity at wall temperature";
    Real HTWallCorrFactorHc1(start = 1.0), HTWallCorrFactorHc2(start = 1.0);
    Real AuxSurfEffHc1(start = 1), AuxSurfEffHc2(start = 1) "fraction of internal auxiliar surfaces used in heat transfer";
    SI.NusseltNumber NuHalfCoil1(start = 1e2, min = 1, max = 1e5), NuHalfCoil2(start = 1e2, min = 1, max = 1e5) "wall jacket process side Nusselt numbers";
    SI.CoefficientOfHeatTransfer HhalfCoil1(min = 10, start = 1000), HhalfCoil2(min = 10, start = 1000) "jacket process side heat transfer coeff.,referenced to vessel surface";
    SI.CoefficientOfHeatTransfer UhalfCoil1(min = 10, start = 600), UhalfCoil2(min = 10, start = 600) "jacket global heat transfer coeff.,referenced to directly heated surface";
    SI.TemperatureDifference LMTDhalfCoil1, LMTDhalfCoil2 "Logarithmic mean temperature difference";
    SI.Power W "heat transfer power possitive to the vessel";
    Medium.ThermodynamicState StateHC1, StateHC2;
  algorithm
    if HalfCoil1.isBottomJacket == true then
      HalfCoil1.HalfCoilDiam := (HalfCoil1.largerHalfCoilDiam + HalfCoil1.lowerHalfCoilDiam) / 2;
    else
      HalfCoil1.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc1 := Medium.dynamicViscosity(StateHC1);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc1 := 1.0;
    else
      HTWallCorrFactorHc1 := max(0.4, (Mu / MuWallHc1) ^ 0.13);
    end if;
    if HalfCoil1.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
        NuHalfCoil1 := 0.9 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 1.1 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil1 := 0.5 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil1 := 1.08 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil1 := 0.31 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil1 := 0.5 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "mig" then
        NuHalfCoil1 := 0.46 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "Ekato";
      elseif Mixer.reference == "intermig" then
        NuHalfCoil1 := 0.54 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "Ekato";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil1 := 0.74 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil1 := 0.54 * HTWallCorrFactorHc1 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 "DeltaT";
        else
          NuHalfCoil1 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc1 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil1 := NuHalfCoil1 * K / Vessel.di;
    AuxSurfEffHc1 := (HhalfCoil1 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc1 := (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) - 1) / (AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi) * (exp(2 * AuxSurfEffHc1 * (HalfCoil1.path - HalfCoil1.basePipeDi)) + 1));
    UhalfCoil1 := 1 / ((foulingF + 1 / HhalfCoil1) * HalfCoil1.SactiveHT / (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil1.H + HalfCoil1.foulingF);
//UhalfCoil1 := HhalfCoil1/(1+foulingF * HhalfCoil1);
    if noEvent((T - HalfCoil1.Ta) * (T - HalfCoil1.Tb) <= 0) then
      LMTDhalfCoil1 := 0;
    elseif noEvent(T - HalfCoil1.Ta - (T - HalfCoil1.Tb) == 0) then
      LMTDhalfCoil1 := T - HalfCoil1.Ta;
    else
      LMTDhalfCoil1 := (T - HalfCoil1.Ta - (T - HalfCoil1.Tb)) / log((T - HalfCoil1.Ta) / (T - HalfCoil1.Tb));
    end if;
    if HalfCoil2.isBottomJacket == true then
      HalfCoil2.HalfCoilDiam := (HalfCoil2.largerHalfCoilDiam + HalfCoil2.lowerHalfCoilDiam) / 2;
    else
      HalfCoil2.HalfCoilDiam := Vessel.Do;
    end if;
    MuWallHc2 := Medium.dynamicViscosity(StateHC2);
    if useHTWallCorrFactor == false then
      HTWallCorrFactorHc2 := 1.0;
    else
      HTWallCorrFactorHc2 := max(0.4, (Mu / MuWallHc2) ^ 0.13);
    end if;
    if HalfCoil2.isBottomJacket == true then
      if Mixer.reference == "hydrofoil" or Mixer.reference == "mig" or Mixer.reference == "intermig" then
        NuHalfCoil2 := 0.9 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 1.1 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "radialTurbine" then
        NuHalfCoil2 := 0.5 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.2 * Vessel.di)) ^ 0.2 "PDHengineer";
      elseif Mixer.reference == "axialTurbine" then
        NuHalfCoil2 := 1.08 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Mixer.nBlades * Mixer.width / (4 * 0.17 * Vessel.di)) ^ 0.2 "PDHengineer";
      end if;
    else
      if Mixer.reference == "hydrofoil" then
        NuHalfCoil2 := 0.31 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 "PDHengineer";
      elseif Mixer.reference == "propeller" then
        NuHalfCoil2 := 0.5 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 "Handbook of    Industrial Mixing. No he entrado la correcion para el angulo de la pala";
      elseif Mixer.reference == "mig" then
        NuHalfCoil2 := 0.46 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "Ekato";
      elseif Mixer.reference == "intermig" then
        NuHalfCoil2 := 0.54 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) "Ekato";
      elseif Mixer.reference == "radialTurbine" or Mixer.reference == "axialTurbine" then
        if Mixer.Re > 10000 then
          NuHalfCoil2 := 0.74 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 "DeltaT, liquid heigth correction from PDHengineer";
        elseif Mixer.Re < 400 then
          NuHalfCoil2 := 0.54 * HTWallCorrFactorHc2 * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 "DeltaT";
        else
          NuHalfCoil2 := (0.1679 + 0.0621 * log(Mixer.Re)) * Mixer.Re ^ (2 / 3) * Pr ^ (1 / 3) * (Vessel.di / hLiquid * Mixer.nImpellers) ^ 0.15 * (5 * Mixer.nBlades * Mixer.width * sin(Mixer.angle) / (6 * Mixer.d)) ^ 0.2 * HTWallCorrFactorHc2 "DeltaT";
        end if;
      end if;
    end if;
    HhalfCoil2 := NuHalfCoil2 * K / Vessel.di;
    AuxSurfEffHc2 := (HhalfCoil2 / (Vessel.wallK * Vessel.wallThickness)) ^ 0.5 "partial calculation";
    AuxSurfEffHc2 := (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) - 1) / (AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi) * (exp(2 * AuxSurfEffHc2 * (HalfCoil2.path - HalfCoil2.basePipeDi)) + 1));
    UhalfCoil2 := 1 / ((foulingF + 1 / HhalfCoil2) * HalfCoil2.SactiveHT / (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT) + Vessel.wallThickness / Vessel.wallK + 1 / HalfCoil2.H + HalfCoil2.foulingF);
//UhalfCoil2 := HhalfCoil2/(1+foulingF * HhalfCoil2);
    if noEvent((T - HalfCoil2.Ta) * (T - HalfCoil2.Tb) <= 0) then
      LMTDhalfCoil2 := 0;
    elseif noEvent(T - HalfCoil2.Ta - (T - HalfCoil2.Tb) == 0) then
      LMTDhalfCoil2 := T - HalfCoil2.Ta;
    else
      LMTDhalfCoil2 := (T - HalfCoil2.Ta - (T - HalfCoil2.Tb)) / log((T - HalfCoil2.Ta) / (T - HalfCoil2.Tb));
    end if;
  equation
//HalfCoil1.W = UhalfCoil1 * LMTDhalfCoil1 * HalfCoil1.SactiveHT;
    HalfCoil1.W = 1 / (1 / HhalfCoil1 + foulingF) * (T - HalfCoil1.Tsurf) * (HalfCoil1.SactiveHT + AuxSurfEffHc1 * HalfCoil1.SauxHT);
//HalfCoil2.W = UhalfCoil2 * LMTDhalfCoil2 * HalfCoil2.SactiveHT;
    HalfCoil2.W = 1 / (1 / HhalfCoil2 + foulingF) * (T - HalfCoil2.Tsurf) * (HalfCoil2.SactiveHT + AuxSurfEffHc2 * HalfCoil2.SauxHT);
//StateHC1 = Medium.setBubbleState(Medium.setSat_T(HalfCoil1.Tsurf));
//StateHC2 = Medium.setBubbleState(Medium.setSat_T(HalfCoil2.Tsurf));
    StateHC1 = Medium.setState_pTX(fixedPressure, HalfCoil1.Tsurf);
    StateHC2 = Medium.setState_pTX(fixedPressure, HalfCoil2.Tsurf);
    W = (-HalfCoil1.W) - HalfCoil2.W;
//it is necessary to improve the SauxHT efficiency factor, not leave it as 0.5
  end TankAgitKamei2hc;

  model InitializationParam
    replaceable package Medium1 = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    replaceable package Medium2 = FreeFluids.TMedia.Fluids.MarlothermSH constrainedby Modelica.Media.Interfaces.PartialMedium;
    parameter FreeFluids.Types.SourceOption TF1Option = FreeFluids.Types.SourceOption.useP_T "themo fluid 1 enthalpy calculation option" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.Temperature TF1initT(displayUnit = "degC") = 553.15 "initial temperature of the thermal fluid 1" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.AbsolutePressure TF1initP(displayUnit = "bar") = 9.0e5 "initial pressure of the thermal fluid 1" annotation(
      Dialog(tab = "Thermal fluid 1"));
    final parameter Modelica.Units.SI.SpecificEnthalpy TF1initH = if TF1Option == FreeFluids.Types.SourceOption.useP_T then Medium1.specificEnthalpy(Medium1.setState_pTX(TF1initP, TF1initT)) elseif TF1Option == FreeFluids.Types.SourceOption.useSatLiqT then Medium1.bubbleEnthalpy(Medium1.setSat_T(TF1initT)) else Medium1.specificEnthalpy(Medium1.setState_pTX(Medium1.p_default, Medium1.T_default)) "initial enthalpy of the thermal fluid 1. Leave it unchanged" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter FreeFluids.Types.SourceOption TF2Option = FreeFluids.Types.SourceOption.useP_T annotation(
      Dialog(tab = "Thermal fluid 2"));
    parameter Modelica.Units.SI.Temperature TF2initT(displayUnit = "degC") = 553.15 "initial temperature of the thermal fluid 2" annotation(
      Dialog(tab = "Thermal fluid 2"));
    parameter Modelica.Units.SI.AbsolutePressure TF2initP(displayUnit = "bar") = 9.0e5 "initial pressure of the thermal fluid 2" annotation(
      Dialog(tab = "Thermal fluid 2"));
    final parameter Modelica.Units.SI.SpecificEnthalpy TF2initH = if TF2Option == FreeFluids.Types.SourceOption.useP_T then Medium2.specificEnthalpy(Medium2.setState_pTX(TF2initP, TF2initT)) else Medium2.specificEnthalpy(Medium2.setState_pTX(Medium2.p_default, Medium2.T_default)) "initial enthalpy of the thermal fluid 2. Leave it unchanged" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.Temperature initT(displayUnit = "degC") = 473.15 "initial temperature of the vessel";
    parameter Modelica.Units.SI.ThermalInsulance foulingF(start = 0.0) "process side fouling factor";
    parameter Modelica.Units.SI.Height hLiquid(start = 1.0) "liquid height from bottom";
    parameter Boolean useHTWallCorr = true "use, or not, wall viscosity correction for heat transfer";
    parameter Modelica.Units.NonSI.AngularVelocity_rpm n(start = 60) "mixer shaft rotational speed rpm";
    parameter Modelica.Units.SI.Power distPower(displayUnit = "kW.h") = 0.0 "power needed for distillation";
    annotation(
      defaultComponentName = "InitialValues",
      Icon(graphics = {Rectangle(origin = {2, 59.4148}, extent = {{-162, 40.5852}, {158, -209.17}}), Text(origin = {-68, 134}, extent = {{-92, -44}, {226, -68}}, textString = "Vessel T: %initT"), Text(origin = {-14, 124}, lineColor = {0, 85, 255}, extent = {{-146, 22}, {174, -24}}, textString = "%name"), Text(origin = {-62, -26}, extent = {{-92, -42}, {224, -68}}, textString = "TF1.init.T: %TF1initT"), Text(origin = {-68, 56}, extent = {{-90, -44}, {226, -70}}, textString = "Liquid: %hLiquid"), Text(origin = {-68, 16}, extent = {{-88, -42}, {224, -70}}, textString = "Mixer: %n"), Text(origin = {-56, 94}, extent = {{-96, -44}, {218, -68}}, textString = "Fouling.F.: %foulingF"), Text(origin = {-62, -66}, extent = {{-92, -42}, {224, -68}}, textString = "distillation power: %distPower")}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html>
      <body>
      <p>This model allows the definition of the initial values for the vessel temperature, the process side fouling factor, the liquid height in the vessel, and the mixer speed.  You can indicate also if you want to use a correction factor for the heat transfer, based on the wall viscosity. And introduce a calculation for the power absorbed by distillation. Expected values for pressure and temperature of two fluids can be also introduced, these values will be used normally for the calculation of the corresponding fluid enthalpy, and this enthalpy can be used as the initial guess for fluid output at pipes or other equipment.</p>
      </body>
      </html>"));
  end InitializationParam;

  model InitializationParamSP
    replaceable package Medium1 = FreeFluids.LMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    replaceable package Medium2 = FreeFluids.LMedia.Fluids.MarlothermSH constrainedby Modelica.Media.Interfaces.PartialMedium;
    parameter FreeFluids.Types.SourceOption TF1Option = FreeFluids.Types.SourceOption.useP_T "themo fluid 1 enthalpy calculation option" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.Temperature TF1initT(displayUnit = "degC") = 553.15 "initial temperature of the thermal fluid 1" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.AbsolutePressure TF1initP(displayUnit = "bar") = 9.0e5 "initial pressure of the thermal fluid 1" annotation(
      Dialog(tab = "Thermal fluid 1"));
    final parameter Modelica.Units.SI.SpecificEnthalpy TF1initH = if TF1Option == FreeFluids.Types.SourceOption.useP_T then Medium1.specificEnthalpy(Medium1.setState_pTX(TF1initP, TF1initT)) else Medium1.specificEnthalpy(Medium1.setState_pTX(Medium1.p_default, Medium1.T_default)) "initial enthalpy of the thermal fluid 1. Leave it unchanged" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter FreeFluids.Types.SourceOption TF2Option = FreeFluids.Types.SourceOption.useP_T annotation(
      Dialog(tab = "Thermal fluid 2"));
    parameter Modelica.Units.SI.Temperature TF2initT(displayUnit = "degC") = 553.15 "initial temperature of the thermal fluid 2" annotation(
      Dialog(tab = "Thermal fluid 2"));
    parameter Modelica.Units.SI.AbsolutePressure TF2initP(displayUnit = "bar") = 9.0e5 "initial pressure of the thermal fluid 2" annotation(
      Dialog(tab = "Thermal fluid 2"));
    final parameter Modelica.Units.SI.SpecificEnthalpy TF2initH = if TF2Option == FreeFluids.Types.SourceOption.useP_T then Medium2.specificEnthalpy(Medium2.setState_pTX(TF2initP, TF2initT)) else Medium2.specificEnthalpy(Medium2.setState_pTX(Medium2.p_default, Medium2.T_default)) "initial enthalpy of the thermal fluid 2. Leave it unchanged" annotation(
      Dialog(tab = "Thermal fluid 1"));
    parameter Modelica.Units.SI.Temperature initT(displayUnit = "degC") = 473.15 "initial temperature of the vessel";
    parameter Modelica.Units.SI.ThermalInsulance foulingF(start = 0.0) "process side fouling factor";
    parameter Modelica.Units.SI.Height hLiquid(start = 1.0) "liquid height from bottom";
    parameter Boolean useHTWallCorr = true "use, or not, wall viscosity correction for heat transfer";
    parameter Modelica.Units.NonSI.AngularVelocity_rpm n(start = 60) "mixer shaft rotational speed rpm";
    parameter Modelica.Units.SI.Power distPower(displayUnit = "kW.h") = 0.0 "power needed for distillation";
    annotation(
      defaultComponentName = "InitialValues",
      Icon(graphics = {Rectangle(origin = {2, 59.4148}, extent = {{-162, 40.5852}, {158, -209.17}}), Text(origin = {-68, 134}, extent = {{-92, -44}, {226, -68}}, textString = "Vessel T: %initT"), Text(origin = {-14, 124}, lineColor = {0, 85, 255}, extent = {{-146, 22}, {174, -24}}, textString = "%name"), Text(origin = {-62, -26}, extent = {{-92, -42}, {224, -68}}, textString = "TF1.init.T: %TF1initT"), Text(origin = {-68, 56}, extent = {{-90, -44}, {226, -70}}, textString = "Liquid: %hLiquid"), Text(origin = {-68, 16}, extent = {{-88, -42}, {224, -70}}, textString = "Mixer: %n"), Text(origin = {-56, 94}, extent = {{-96, -44}, {218, -68}}, textString = "Fouling.F.: %foulingF"), Text(origin = {-62, -66}, extent = {{-92, -42}, {224, -68}}, textString = "distillation power: %distPower")}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html><head></head><body>The same as InitializationParam but for single phase mediums.</body></html>"));
  end InitializationParamSP;

  package Examples
    extends Modelica.Icons.ExamplesPackage;

    model VesselSimpleWithPump
      extends Modelica.Icons.Example;
      FreeFluids.Vessels.VesselSimple Vessel(redeclare package MediumL = FreeFluids.LMedia.Fluids.Water, height = 5, redeclare package MediumG = Modelica.Media.Air.DryAirNasa, initialLiquidVolume = 1, initialT = 298.15, inletL = 0.8, lgConductance = 2, overflowL = 1.5, section = 2, ventOutL = 2, vesselVolume = 10)  annotation(
        Placement(visible = true, transformation(origin = {44, -8}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, Elevation = 1, P(displayUnit = "bar") = 99999.99999999999, T = 323.15) annotation(
        Placement(visible = true, transformation(origin = {-56, -24}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      FreeFluids.Pumps.BumpPump Pump(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, forceSpeed = true, h0 = 29.1, h1 = 27, h2 = 23, n0 = 2900 / 60, numParallelUnits = 1, q1 = 0.0016667, q2 = 0.0025, r1 = 0.705, r2 = 0.73) annotation(
        Placement(visible = true, transformation(origin = {-16, -24}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    equation
  connect(Pump.PortB, Vessel.Inlet) annotation(
        Line(points = {{-10, -24}, {24, -24}}, color = {0, 127, 255}));
  connect(Source.PortB, Pump.PortA) annotation(
        Line(points = {{-50, -24}, {-22, -24}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-60, 20}, {80, -40}})),
        experiment(StartTime = 0, StopTime = 1700, Tolerance = 1e-06, Interval = 7.2),
  Documentation(info = "<html><head></head><body>Vessel preswsurization with a bump pump.</body></html>"));
    end VesselSimpleWithPump;
  
    model VesselLevelTest
  extends Modelica.Icons.Example;
      VesselLevel Vessel(bottomHeight = 0.25, di = 1, fixLiquidLevel = true, liquidL = 0.5, liquidV = 3, shellLength = 2, shellThickness = 0.001, topHeight = 0.25, vesselBottom = FreeFluids.Types.HeadShape.hemispherical, vesselTop = FreeFluids.Types.HeadShape.Korbbogen) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      annotation(
        Diagram(coordinateSystem(extent = {{-20, 20}, {20, -20}})));
    end VesselLevelTest;

    model VesselDetailedWithPVvalve
      extends Modelica.Icons.Example;
      FreeFluids.Vessels.VesselDetailed Vessel(bottomHeight = 0.25, di = 1.5, initialLiquidVolume = 1, initialP = 99999.99999999999, initialT = 298.15, inletL = 0.8, lgConductance = 5, outletL = 0.0, overflowL = 1.5, shellLength = 2, shellThickness = 0.001, topHeight = 0.25, ventOutL = 2, vesselBottom = FreeFluids.Types.HeadShape.semielliptical, vesselTop = FreeFluids.Types.HeadShape.semielliptical, vesselVolume = 4.41786) annotation(
        Placement(visible = true, transformation(origin = {36, -6}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, Elevation = 1, P(displayUnit = "bar"), T = 323.15, externalP = true) annotation(
        Placement(visible = true, transformation(origin = {-24, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  FreeFluids.Pipes.PipeFlow1Ph InletPipe(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.015, lTube = 10)  annotation(
        Placement(visible = true, transformation(origin = {-22, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Pipes.PipeFlow1Ph OutletPipe(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.022, lTube = 3) annotation(
        Placement(visible = true, transformation(origin = {94, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Interfaces.FlowSink Sink(P = 102000, fix = FreeFluids.Types.BoundaryOption.fixPressure)  annotation(
        Placement(visible = true, transformation(origin = {96, -6}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  FreeFluids.Pipes.PipeFlow1Ph VentPipe(redeclare package Medium = Modelica.Media.Air.DryAirNasa, di = 0.01, lTube = 4) annotation(
        Placement(visible = true, transformation(origin = {78, 54}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  FreeFluids.Interfaces.FlowSink Vent(P = 101325, fix = FreeFluids.Types.BoundaryOption.fixPressure)  annotation(
        Placement(visible = true, transformation(origin = {108, 54}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  FreeFluids.Valves.PressureVacuumValve PSV1(pressureFlow = 0.0005555555555555556, pressureSet = 103325, rho(displayUnit = "kg/m3"), vacuumFlow = 0.0005555555555555556, vacuumSet = 99324.99999999999)  annotation(
        Placement(visible = true, transformation(origin = {57.2319, 48.0414}, extent = {{-6.04142, 4.64725}, {6.04142, 16.7301}}, rotation = 90)));
  Modelica.Blocks.Sources.Ramp SourceRamp(duration = 3600, height = 2e5, offset = 1.1e5)  annotation(
        Placement(visible = true, transformation(origin = {-34, 42}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    equation
      connect(Source.PortB, InletPipe.PortA) annotation(
        Line(points = {{-18, 8}, {-18, 8.5}, {-10, 8.5}, {-10, -7}, {-38, -7}, {-38, -31.5}, {-32, -31.5}, {-32, -32}}, color = {0, 127, 255}));
      connect(InletPipe.PortB, Vessel.Inlet) annotation(
        Line(points = {{-12, -32}, {4, -32}}, color = {0, 127, 255}));
      connect(Vessel.Outlet, OutletPipe.PortA) annotation(
        Line(points = {{68, -32}, {84, -32}}, color = {0, 127, 255}));
      connect(OutletPipe.PortB, Sink.PortA) annotation(
        Line(points = {{104, -32}, {113, -32}, {113, -6}, {102, -6}}, color = {0, 127, 255}));
  connect(VentPipe.PortB, Vent.PortA) annotation(
        Line(points = {{86, 54}, {102, 54}}, color = {0, 127, 255}));
  connect(PSV1.PortB, VentPipe.PortA) annotation(
        Line(points = {{47, 54}, {70, 54}}, color = {0, 127, 255}));
  connect(SourceRamp.y, Source.Pext) annotation(
        Line(points = {{-28, 42}, {-20, 42}, {-20, 14}}, color = {0, 0, 127}));
  connect(Vessel.VentOut, PSV1.PortA) annotation(
        Line(points = {{46, 28}, {46, 42}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-40, 60}, {120, -40}})),
        experiment(StartTime = 0, StopTime = 3600, Tolerance = 1e-06, Interval = 7.2));
    end VesselDetailedWithPVvalve;
  
    model VesselDetailedClosed
      extends Modelica.Icons.Example;
      FreeFluids.Vessels.VesselDetailed Vessel(bottomHeight = 0.25, di = 1.5, initialLiquidVolume = 1, initialP = 99999.99999999999, initialT = 298.15, inletL = 0.8, outletL = 0.0, overflowL = 1.5, shellLength = 2, shellThickness = 0.001, topHeight = 0.25, ventOutL = 2, vesselBottom = FreeFluids.Types.HeadShape.semielliptical, vesselTop = FreeFluids.Types.HeadShape.semielliptical, vesselVolume = 4.41786) annotation(
        Placement(visible = true, transformation(origin = {84, 0}, extent = {{-34, -34}, {34, 34}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, Elevation = 1, P(displayUnit = "bar") = 150000, T = 323.15) annotation(
        Placement(visible = true, transformation(origin = {-20, -24}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  FreeFluids.Pipes.PipeFlow1Ph InletPipe(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.015, lTube = 10)  annotation(
        Placement(visible = true, transformation(origin = {14, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, InletPipe.PortA) annotation(
        Line(points = {{-14, -24}, {4, -24}}, color = {0, 127, 255}));
      connect(InletPipe.PortB, Vessel.Inlet) annotation(
        Line(points = {{24, -24}, {54, -24}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-40, 40}, {120, -40}})),
        experiment(StartTime = 0, StopTime = 3600, Tolerance = 1e-06, Interval = 7.2));
    end VesselDetailedClosed;
  
    model VesselDetailedPressureRelief
      extends Modelica.Icons.Example;
      FreeFluids.Vessels.VesselDetailed Vessel(bottomHeight = 0.25, di = 1.5, initialLiquidVolume = 1, initialP = 110000, initialT = 298.15, inletL = 0.8, outletL = 0.0, overflowL = 1.5, shellLength = 2, shellThickness = 0.001, topHeight = 0.25, ventOutL = 2, vesselBottom = FreeFluids.Types.HeadShape.semielliptical, vesselTop = FreeFluids.Types.HeadShape.semielliptical, vesselVolume = 4.41786) annotation(
        Placement(visible = true, transformation(origin = {36, -6}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, Elevation = 1, P(displayUnit = "bar") = 180000, T = 323.15) annotation(
        Placement(visible = true, transformation(origin = {-18, -6}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
    FreeFluids.Pipes.PipeFlow1Ph InletPipe(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.015, lTube = 10)  annotation(
        Placement(visible = true, transformation(origin = {-22, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph OutletPipe(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.015, lTube = 3) annotation(
        Placement(visible = true, transformation(origin = {94, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(P = 103000, fix = FreeFluids.Types.BoundaryOption.fixPressure)  annotation(
        Placement(visible = true, transformation(origin = {96, -6}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
    FreeFluids.Pipes.PipeFlow1Ph VentPipe(redeclare package Medium = Modelica.Media.Air.DryAirNasa, di = 0.01, lTube = 4) annotation(
        Placement(visible = true, transformation(origin = {76, 48}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Vent(P = 105000, fix = FreeFluids.Types.BoundaryOption.fixPressure)  annotation(
        Placement(visible = true, transformation(origin = {106, 48}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  FreeFluids.Valves.PressureReliefValve PSV(flowRating = 0.0002777777777777778, pSet = 110000)  annotation(
        Placement(visible = true, transformation(origin = {46, 42}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    equation
      connect(Source.PortB, InletPipe.PortA) annotation(
        Line(points = {{-24, -6}, {-24, -7}, {-38, -7}, {-38, -31.5}, {-32, -31.5}, {-32, -32}}, color = {0, 127, 255}));
      connect(InletPipe.PortB, Vessel.Inlet) annotation(
        Line(points = {{-12, -32}, {4, -32}}, color = {0, 127, 255}));
      connect(Vessel.Outlet, OutletPipe.PortA) annotation(
        Line(points = {{68, -32}, {84, -32}}, color = {0, 127, 255}));
      connect(OutletPipe.PortB, Sink.PortA) annotation(
        Line(points = {{104, -32}, {113, -32}, {113, -6}, {102, -6}}, color = {0, 127, 255}));
      connect(VentPipe.PortB, Vent.PortA) annotation(
        Line(points = {{84, 48}, {100, 48}}, color = {0, 127, 255}));
  connect(Vessel.VentOut, PSV.PortA) annotation(
        Line(points = {{46, 28}, {46, 36}}, color = {0, 127, 255}));
  connect(PSV.PortB, VentPipe.PortA) annotation(
        Line(points = {{46, 48}, {68, 48}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-40, 60}, {120, -40}})),
        experiment(StartTime = 0, StopTime = 3600, Tolerance = 1e-06, Interval = 7.2));
    end VesselDetailedPressureRelief;

    model AnchorVessel
      extends FreeFluids.Vessels.TankAgitAnchor2Coils1Hc(
      Coil1(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 0.0669, coilDiam = 1.7, num = 12, numTubes = 3, path = 0.142, rhoWall = 8000, thickness = 3.05e-3, numActiveTubes = 3, foulingF = 0.0004), 
      Coil2(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 0.0669, coilDiam = 2.0, num = 12, numTubes = 3, path = 0.142, rhoWall = 8000, thickness = 3.05e-3, numActiveTubes = 3, foulingF = 0.0004), 
      HalfCoil1(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, basePipeDi = 0.083, num = 12, numActiveTubes = 2, numTubes = 2, path = 0.115, thickness = 0.01, Ltube(start = 60.3186), foulingF = 0.0004, fixedW = -25000, thermalType = FreeFluids.Types.ThermalType.detailed), Vessel(cylinderH = 3.2, di = 2.5, nBaffles = 4, baffleWidth = 0.15, wallThickness = 0.01, topShape = "Korbbogen", nConcCoils = 2), Mixer(d = 1.5, width = 0.16, hTotal = 2.2, n = 35));
    end AnchorVessel;
    model AnchorKettle
      extends AnchorVessel(Coil1(Tsurf.start = InitialValues.TF1initT, PortB.H(start = InitialValues.TF1initH)), Coil2(Tsurf.start = InitialValues.TF1initT, PortB.H(start = InitialValues.TF1initH)), HalfCoil1(Tsurf.start = InitialValues.TF1initT, PortB.H(start = InitialValues.TF1initH)), hLiquid = InitialValues.hLiquid, useHTWallCorrFactor = InitialValues.useHTWallCorr, Mixer.n = InitialValues.n, T(start = InitialValues.initT));
      FreeFluids.Pipes.PipeFlow1Ph PipeToVessel(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 206.5e-3, equivL_Di = 2 * 16 + 340 + 2 * 16 + 2 * 340, kv = 150, lTube = 5, thermalType = FreeFluids.Types.ThermalType.isenthalpic, useTubeLength = true, PortB.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {-62, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, Elevation = 0, P = 231000, T = 553.15, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-137, 53}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, P = 200000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {67, -35}, extent = {{-9, -9}, {9, 9}}, rotation = 180)));
      FreeFluids.Pipes.Mixer3PH MixerOut(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, PortD.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {54, 0}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV1(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, fix = FreeFluids.Types.ValveFixOption.fixKv, fixedKv = 250, isLinear = false) annotation(
        Placement(visible = true, transformation(origin = {-114, -2}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible TV01(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, aperture = 0.7, fix = FreeFluids.Types.ValveFixOption.fixKv, fixedKv = 150, isLinear = false) annotation(
        Placement(visible = true, transformation(origin = {107, 19}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
      FreeFluids.Pipes.MixerPH MixerIn(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, PortC.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {104, -42}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      FreeFluids.Pumps.BumpPump Pump(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, Qtotal(start = 0.0463889), Qunit(start = 0.0463889), fixedSpeed = 0.85 * 2900 / 60, forceSpeed = true, h0 = 49.47, h1 = 46.82, h2 = 37.52, n0 = 2900 / 60, q1 = 0.0409639, q2 = 0.0689667, r1 = 0.729, r2 = 0.813, PortB.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {-138, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV2(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, fix = FreeFluids.Types.ValveFixOption.fixKv, fixedKv = 350, isLinear = false) annotation(
        Placement(visible = true, transformation(origin = {-128, -34}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
      FreeFluids.HeatExchangers.HEXgeneric1Ph Exchanger(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 0.02058, extG = 22.3439, extSratio = 20, extTin = 288.15, numPasses = 3, numRows = 6, numTubesTotal = 132, numVelocityHeads = 4.5, perimeter = 0.0735, section = 0.0002851, thickness = 0.00211, tubeLength = 4, u = 11.39) annotation(
        Placement(visible = true, transformation(origin = {-112, -66}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph PipeExchanger(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 0.1554, equivL_Di = 12 * 16 + 2 * 340 + 2 * 35, lTube = 24, thermalType = FreeFluids.Types.ThermalType.isenthalpic, thickness = 4.85e-3, useTubeLength = true, PortB.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {-96, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      FreeFluids.Pipes.MixerPH MixerExchanger(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, PortC.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {-90, -6}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Vessels.InitializationParamSP InitialValues(redeclare package Medium1 = FreeFluids.LMedia.Fluids.MarlothermSH, TF1initT = 523.15, foulingF = 0, hLiquid = 2.4, initT = 298.15, n = 27) annotation(
        Placement(visible = true, transformation(origin = {2, 78}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph PipeFromVessel(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 206.5e-3, lTube = 5, thermalType = FreeFluids.Types.ThermalType.isenthalpic, useTubeLength = true, PortB.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(PipeExchanger.PortB, MixerExchanger.PortB) annotation(
        Line(points = {{-96, -26}, {-96, -9}, {-97, -9}}, color = {0, 127, 255}));
      connect(Exchanger.PortB, PipeExchanger.PortA) annotation(
        Line(points = {{-104, -66}, {-96, -66}, {-96, -46}}, color = {0, 127, 255}));
      connect(PipeToVessel.PortA, MixerExchanger.PortC) annotation(
        Line(points = {{-72, -6}, {-83, -6}}, color = {0, 127, 255}));
      connect(FV1.PortB, MixerExchanger.PortA) annotation(
        Line(points = {{-109, -2}, {-98, -2}, {-98, -3}, {-97, -3}}, color = {0, 127, 255}));
      connect(Pump.PortB, FV1.PortA) annotation(
        Line(points = {{-128, -2}, {-119, -2}}, color = {0, 127, 255}));
      connect(Exchanger.PortA, FV2.PortB) annotation(
        Line(points = {{-120, -66}, {-120, -66.75}, {-128, -66.75}, {-128, -39}}, color = {0, 127, 255}));
      connect(Pump.PortB, FV2.PortA) annotation(
        Line(points = {{-128, -2}, {-128, -29}}, color = {0, 127, 255}));
      connect(Pump.PortA, MixerIn.PortC) annotation(
        Line(points = {{-148, -2}, {-148.5, -2}, {-148.5, -88}, {103.25, -88}, {103.25, -49}, {104, -49}}, color = {0, 127, 255}));
      connect(Source.PortB, TV01.PortA) annotation(
        Line(points = {{-128, 53}, {106.75, 53}, {106.75, 25}, {107, 25}}, color = {0, 127, 255}));
      connect(HalfCoil1.PortB, MixerOut.PortA) annotation(
        Line(points = {{-34, 24}, {-33, 24}, {-33, 36}, {46, 36}, {46, 5}, {47, 5}}, color = {0, 127, 255}));
      connect(Sink.PortA, MixerIn.PortB) annotation(
        Line(points = {{76, -35}, {101, -35}}, color = {0, 127, 255}));
      connect(TV01.PortB, MixerIn.PortA) annotation(
        Line(points = {{107, 13}, {107, -35}}, color = {0, 127, 255}));
      connect(PipeToVessel.PortB, Coil1.PortA) annotation(
        Line(points = {{-52, -6}, {-10, -6}}, color = {0, 127, 255}));
      connect(PipeToVessel.PortB, Coil2.PortA) annotation(
        Line(points = {{-52, -6}, {-12, -6}, {-12, 16}}, color = {0, 127, 255}));
      connect(PipeToVessel.PortB, HalfCoil1.PortA) annotation(
        Line(points = {{-52, -6}, {-34, -6}, {-34, 0}}, color = {0, 127, 255}));
      connect(Coil2.PortB, MixerOut.PortB) annotation(
        Line(points = {{10, 16}, {38, 16}, {38, 0}, {46, 0}}, color = {0, 127, 255}));
      connect(Coil1.PortB, MixerOut.PortC) annotation(
        Line(points = {{10, -6}, {48, -6}, {48, -4}, {46, -4}}, color = {0, 127, 255}));
      connect(MixerOut.PortD, PipeFromVessel.PortA) annotation(
        Line(points = {{62, 0}, {74, 0}, {74, 0}, {74, 0}}, color = {0, 127, 255}));
      connect(PipeFromVessel.PortB, MixerIn.PortB) annotation(
        Line(points = {{94, 0}, {102, 0}, {102, -35}, {101, -35}}, color = {0, 127, 255}));
      der(T) = W / (Vliquid * Rho * Cp + (Vessel.VesselMass + 2000.0) * Vessel.wallCp);
      annotation(
        Diagram(coordinateSystem(extent = {{-160, 100}, {140, -100}})));
    end AnchorKettle;

    model AnchorKettleHeating
      extends AnchorKettle(redeclare package Medium = Modelica.Media.Water.StandardWater, FV2.aperture = 0.0001, InitialValues.TF1initT = 403.15, Exchanger.thermalType = FreeFluids.Types.ThermalType.isenthalpic, Source.T = 453.15, Source.P.displayUnit = "Pa", TV01.aperture = 0.6);
      annotation(
        experiment(StartTime = 0, StopTime = 1800, Tolerance = 1e-06, Interval = 3.6));
    end AnchorKettleHeating;

    model AnchorKettleCooling
      extends AnchorKettle(redeclare package Medium = Modelica.Media.Water.StandardWater, FV2.aperture = 1, InitialValues.TF1initT = 368.15, Exchanger.thermalType = FreeFluids.Types.ThermalType.detailed, Source.T = 363.15, Source.P.displayUnit = "Pa", TV01.aperture = 0, FV1.aperture = 0, InitialValues.initT = 373.15);
      annotation(
        experiment(StartTime = 0, StopTime = 1800, Tolerance = 1e-06, Interval = 3.6));
    end AnchorKettleCooling;

    model TurbineVessel
      extends FreeFluids.Vessels.TankAgitKamei2Coils(Vessel(baffleWidth = 0.2, cylinderH = 2, di = 2.5, nBaffles = 4, wallThickness = 0.008, bottomShape = "Klopper", topShape = "Korbbogen", nConcCoils = 2), Mixer(angle(displayUnit = "rad") = 0.5235987755982988, d = 1.6, doShaftCalc = true, hBottom = 0.6, interImpellerDistance = 0.8, nBlades = 3, nImpellers = 2, width = 0.45, n = 60), Coil1(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, coilDiam = 1.7, di = 0.06693, num = 3, numTubes = 1, path = 0.146, thickness = 0.00305, SusedHT(start = 6.3), PLossFriction.displayUnit = "Pa", heightInit = 0), Coil2(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, coilDiam = 2.0, di = 0.06693, num = 3, numTubes = 1, path = 0.146, thickness = 0.00305, SusedHT(start = 8.9), PLossFriction.displayUnit = "Pa", heightInit = 0.4, thermalType = FreeFluids.Types.ThermalType.detailed));
    equation
    
    end TurbineVessel;

    model TurbineKettle
      extends TurbineVessel(Coil1(Tsurf.start = InitialValues.initT, PortB.H(start = InitialValues.TF1initH)), Coil2(Tsurf.start = InitialValues.initT, PortB.H(start = InitialValues.TF1initH)), hLiquid = InitialValues.hLiquid, useHTWallCorrFactor = InitialValues.useHTWallCorr, Mixer.n = InitialValues.n, T(start = InitialValues.initT));
      FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 0, redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, P(displayUnit = "bar") = 899999.9999999998, T(displayUnit = "degC") = 553.15) annotation(
        Placement(visible = true, transformation(origin = {74, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, P = 880000, PortA(P(start = Source.P - 0.5e5)), fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {84, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Vessels.InitializationParamSP InitialValues(redeclare package Medium1 = FreeFluids.LMedia.Fluids.ShellS2, n = 55) annotation(
        Placement(visible = true, transformation(origin = {12, 62}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      FreeFluids.Pipes.MixerPH MixerOut(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, PortC.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {66, 8}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Pipes.MixerPH MixerIn(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, PortC.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {99, 11}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      FreeFluids.Pumps.BumpPump P10(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, forceSpeed = true, h0 = 50, h1 = 47, h2 = 34, n0 = 2900 / 60, q1(displayUnit = "m3/h") = 0.01944444444444444, q2(displayUnit = "m3/h") = 0.02777777777777778, r1 = 0.76, r2 = 0.73) annotation(
        Placement(visible = true, transformation(origin = {-145, 3}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      FreeFluids.Pipes.PdiffSource DeltaP(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, dP(displayUnit = "bar") = -235000, refG(displayUnit = "kg/h") = 20.83333333333333, useFixedDiffP = false) annotation(
        Placement(visible = true, transformation(origin = {-58, 0}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible TV01(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, fixedKv = 160, isLinear = false, useFixedAperture = false) annotation(
        Placement(visible = true, transformation(origin = {92, 44}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant Tmax(k = 270) annotation(
        Placement(visible = true, transformation(origin = {-152, 94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Math.UnitConversions.From_degC From_degC annotation(
        Placement(visible = true, transformation(origin = {-132, 94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Math.Min Min annotation(
        Placement(visible = true, transformation(origin = {-105, 89}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant DeltaT(k = 60) annotation(
        Placement(visible = true, transformation(origin = {-152, 72}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Modelica.Blocks.Math.Add Add annotation(
        Placement(visible = true, transformation(origin = {-132, 68}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV01(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, fixedKv = 230, useFixedAperture = false) annotation(
        Placement(visible = true, transformation(origin = {-117, 3}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV02(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, aperture = 0.001, fixedKv = 230, useFixedAperture = false) annotation(
        Placement(visible = true, transformation(origin = {-136, -30}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      FreeFluids.HeatExchangers.HEXgeneric1Ph E01(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, di = 0.02058, extG = 6.944444444444445, extSratio = 20, extTin = 303.15, fixExternalFlow = true, isCircular = true, numPasses = 3, numRows = 6, numTubesTotal = 120, tubeLength = 4, u = 11.39, useDiameter = true) annotation(
        Placement(visible = true, transformation(origin = {-112, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.MixerPH MixerE01(redeclare package Medium = FreeFluids.LMedia.Fluids.ShellS2, PortC.H(start = InitialValues.TF1initH)) annotation(
        Placement(visible = true, transformation(origin = {-86, 0}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Instruments.SplitRange SplitRange annotation(
        Placement(visible = true, transformation(origin = {-32, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Instruments.PI TIC10(T = 400, initType = Modelica.Blocks.Types.Init.NoInit, k = 10) annotation(
        Placement(visible = true, transformation(origin = {-62, 84}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    equation
      der(T) = (W - InitialValues.distPower) / (Vliquid * Rho * Cp + (Vessel.VesselMass + 1920.0) * Vessel.wallCp);
      connect(Coil2.PortB, MixerOut.PortA) annotation(
        Line(points = {{12, 26}, {60, 26}, {60, 12}}, color = {0, 127, 255}));
      connect(MixerOut.PortC, Sink.PortA) annotation(
        Line(points = {{74, 8}, {73, 8}, {73, -18}, {74, -18}}));
      connect(DeltaP.PortB, Coil2.PortA) annotation(
        Line(points = {{-50, 0}, {-12, 0}, {-12, 26}}, color = {0, 127, 255}));
      connect(MixerOut.PortC, MixerIn.PortB) annotation(
        Line(points = {{73, 8}, {93, 8}}, color = {0, 127, 255}));
      connect(MixerIn.PortC, P10.PortA) annotation(
        Line(points = {{105, 11}, {113, 11}, {113, -84}, {-154, -84}, {-154, 3}}, color = {0, 127, 255}));
      connect(Source.PortB, TV01.PortA) annotation(
        Line(points = {{84, 60}, {92, 60}, {92, 51}}, color = {0, 127, 255}));
      connect(TV01.PortB, MixerIn.PortA) annotation(
        Line(points = {{92, 37}, {93, 37}, {93, 14}}, color = {0, 127, 255}));
      connect(Tmax.y, From_degC.u) annotation(
        Line(points = {{-145, 94}, {-139, 94}}, color = {0, 0, 127}));
      connect(From_degC.y, Min.u1) annotation(
        Line(points = {{-125, 94}, {-119, 94}, {-119, 93}, {-113, 93}}, color = {0, 0, 127}));
      connect(DeltaT.y, Add.u1) annotation(
        Line(points = {{-145, 72}, {-139, 72}}, color = {0, 0, 127}));
      connect(Tmass, Add.u2) annotation(
        Line(points = {{-20, 40}, {-139, 40}, {-139, 64}}, color = {0, 0, 127}));
      connect(Add.y, Min.u2) annotation(
        Line(points = {{-125, 68}, {-120, 68}, {-120, 85}, {-113, 85}}, color = {0, 0, 127}));
      connect(P10.PortB, FV01.PortA) annotation(
        Line(points = {{-136, 3}, {-123, 3}}));
      connect(P10.PortB, FV02.PortA) annotation(
        Line(points = {{-136, 3}, {-136, -23}}, color = {0, 127, 255}));
      connect(FV02.PortB, E01.PortA) annotation(
        Line(points = {{-136, -37}, {-136, -50}, {-122, -50}}, color = {0, 127, 255}));
      connect(FV01.PortB, MixerE01.PortA) annotation(
        Line(points = {{-111, 3}, {-93, 3}}, color = {0, 127, 255}));
      connect(E01.PortB, MixerE01.PortB) annotation(
        Line(points = {{-102, -50}, {-93, -50}, {-93, -3}}, color = {0, 127, 255}));
      connect(SplitRange.y1D, TV01.Opening) annotation(
        Line(points = {{-21, 91}, {106, 91}, {106, 44}, {98, 44}}, color = {0, 0, 127}));
      connect(SplitRange.y2D, FV02.Opening) annotation(
        Line(points = {{-21, 81}, {-15.5, 81}, {-15.5, 52}, {-130, 52}, {-130, -30}}, color = {0, 0, 127}));
      connect(SplitRange.y2I, FV01.Opening) annotation(
        Line(points = {{-21, 77}, {-19.5, 77}, {-19.5, 58}, {-117, 58}, {-117, 8}}, color = {0, 0, 127}));
      connect(Min.y, TIC10.u1) annotation(
        Line(points = {{-97, 89}, {-72, 89}}, color = {0, 0, 127}));
      connect(TIC10.y, SplitRange.u1) annotation(
        Line(points = {{-53, 84}, {-44, 84}}, color = {0, 0, 127}));
      connect(MixerE01.PortC, DeltaP.PortA) annotation(
        Line(points = {{-78, 0}, {-66, 0}}, color = {0, 127, 255}));
      connect(MixerE01.Tout, TIC10.u2) annotation(
        Line(points = {{-80, 4}, {-80, 79.5}, {-72, 79.5}, {-72, 79}}, color = {0, 0, 127}));
      connect(Coil1.PortB, MixerOut.PortB) annotation(
        Line(points = {{10, -14}, {60, -14}, {60, 4}}, color = {0, 127, 255}));
      connect(DeltaP.PortB, Coil1.PortA) annotation(
        Line(points = {{-50, 0}, {-10, 0}, {-10, -14}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-160, 100}, {120, -80}})));
    end TurbineKettle;

    model TurbineKettleHeating
      extends TurbineKettle(redeclare package Medium = FreeFluids.LMedia.Fluids.EG, InitialValues.TF1initT = 328.15, InitialValues.initT = 323.15, Source.T = 523.15, InitialValues.hLiquid = 1.8, InitialValues.foulingF = 0.0, Source.T.displayUnit = "degC", Source.P = 899999.9999999999, DeltaP.refG.displayUnit = "kg/h", DeltaT.k = 60, E01.thermalType = FreeFluids.Types.ThermalType.isenthalpic, InitialValues.distPower.displayUnit = "W", Tmax.k = 200, P10.q1 = 0.01777777777777778, P10.q2 = 0.02666666666666667);
    equation

      annotation(
        experiment(StartTime = 0, StopTime = 12000, Tolerance = 1e-06, Interval = 20),
        Documentation(info = "<html><head></head><body>This is an example of a controlled heating, where the thermal oil temperature is controlled to some degrees over the kettle temperature. The maximum oil temperature is also limited.<div>As this is a complex example, convergence is difficult, and you must play with the guessed temperature for oil outlet in the InitialValues template.</div><div>The heat exchanger has been configured to isenthalpic in order not to make any calculation regarding exchanged heat.</div></body></html>"));
    end TurbineKettleHeating;

    model TurbineKettleCooling
      extends TurbineKettle(redeclare package Medium = FreeFluids.LMedia.Fluids.EG, InitialValues.TF1initT = 453.15, InitialValues.initT = 473.15, Source.T = 523.15, InitialValues.hLiquid = 1.8, InitialValues.foulingF = 0.0, Source.T.displayUnit = "degC", Source.P = 899999.9999999999, DeltaP.refG.displayUnit = "kg/h", DeltaT.k = -60, E01.thermalType = FreeFluids.Types.ThermalType.detailed, Tmax.k = 200);
    equation

      annotation(
        experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 20));
    end TurbineKettleCooling;
    
  end Examples;
  annotation(
    Icon(coordinateSystem(initialScale = 0.1)),
    Documentation(info = "<html><head></head><body>The package contains models vessels, with special focus in heat transfer and power consumption in agitated vessels.</body></html>"));
end Vessels;

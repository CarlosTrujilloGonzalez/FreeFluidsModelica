within FreeFluids;

package Valves "Valves.mo by Carlos Trujillo
  This file is part of the Free Fluids application
  Copyright (C) 2008-2021  Carlos Trujillo Gonzalez
    
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
  import FreeFluids.Types.*;

  partial model ValveBase "General valve model"
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter FreeFluids.Types.ValveFixOption fix = FreeFluids.Types.ValveFixOption.fixKv "select the value to fix: Kv, pressure drop, or flow" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.Area fixedKv = 1.0 "Kv value if it is fixed" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.PressureDifference fixedDP(displayUnit = "bar") = -0.5e5 "negative if PortB.P<PortA.P" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.MassFlowRate fixedFlow(displayUnit = "kg/h")=1.0 "fixed mass flow to maintain at Port A. Possitive if flow is in" annotation(
      Dialog(tab = "Flow"));
    Modelica.Units.SI.Area Kv(min = 0.0001, max = 2000, start = 1);
    Modelica.Units.SI.Area Cv(min = 0.0001, max = 2000, start = 1);
    Modelica.Units.SI.VolumeFlowRate Q(displayUnit = "m3/h", start = 1);
    Medium.ThermodynamicState State;
    Medium.Density Rho(displayUnit = "kg/m3");
    Medium.Temperature T(displayUnit = "degC") "Temperature";
  equation
    Cv = 1.156 * Kv;
    if fix == FreeFluids.Types.ValveFixOption.fixDP then
      Pdiff = fixedDP;
    elseif fix == FreeFluids.Types.ValveFixOption.fixFlow then
      PortA.G = fixedFlow;
    else
      Kv = fixedKv;
    end if;
    Rho = Medium.density(State);
    Q = PortA.G / Rho;
    T = Medium.temperature(State);
  annotation(defaultComponentName = "Valve",
      Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-4498\">It is an extension of the Interfaces.TwoFluidPorts model. It adds the following parameters to configure a two ports valve:</div>

<ul class=\"itemize\" id=\"magicparlabel-4499\"><li class=\"itemize_item\">parameter FreeFluids.Types.ValveFixOption fix = FreeFluids.Types.ValveFixOption.fixKv, in order to select if we are going to specify a fized Kv, a fixed pressure drop or a fixed flow along the valve.</li>
<li class=\"itemize_item\">And the following parameters: fixedKv, fixedDP and fixedFlow, in order to allow the user specification of the needed value.</li>
</ul>
<div class=\"standard\" id=\"magicparlabel-4501\">It adds also the following variables:</div>

<div class=\"standard\" id=\"magicparlabel-4502\">Kv for the valve Kv</div>

<div class=\"standard\" id=\"magicparlabel-4503\">Cv for the valve Cv</div>

<div class=\"standard\" id=\"magicparlabel-4504\">Medium.ThermodynamicState State, in order to specify the state at which physical properties will be evaluated.</div>

<div class=\"standard\" id=\"magicparlabel-4505\">Q for the volumetric flow at State conditions.</div>

<div class=\"standard\" id=\"magicparlabel-4506\">T for fluid temperature at State conditions.</div>

<div class=\"standard\" id=\"magicparlabel-4507\">Rho for fluid density at State conditions.</div></body></html>"));
  end ValveBase;

  partial model ValvePartial "General control valve model"
    extends ValveBase;
    parameter Boolean isCompressibleFlow = false "if compressible, it will be treated always as adiabatic" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean isLinear = true "If true, linear characteristic is applied, otherwise isoporcentual" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean useFixedAperture = true "if true the aperture parameter is used. Otherwise aperture is supplied by the connector" annotation(
      Dialog(tab = "Flow"));
    parameter Real aperture = 1.0 "maximum just 1 valve should be closed (aperture=0) in a line" annotation(
      Dialog(tab = "Flow"));
    SI.Area KvFlow(start = 1.0) "Kv at actual aperture";
    Modelica.Blocks.Interfaces.RealInput Opening(min = 0.0, max = 1.0) if useFixedAperture == false annotation(
      Placement(visible = true, transformation(origin = {0, 94}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 88}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  protected
    Modelica.Blocks.Interfaces.RealInput Aperture;
    //Real Aperture;
  equation
    connect(Opening, Aperture);
    if useFixedAperture == true then
      Aperture = aperture;
    end if;
//For isoporcentual, it is equivalent to use rangeability=50(industry standard) and KvFlow=Kv*rangeability^(aperture-1). This would produce 2% of Kv at aperture=0. So, below 10% aperture a linear function is used. If we go down with rangeability we can go down with the limit at which we apply the linear function.
//Av=28e-6*Kv
    KvFlow = if isLinear == true then Kv * Aperture elseif Aperture < 0.1 then Kv * Aperture * 0.296 else Kv * exp(Aperture * 3.91 - 3.91) "exp(3.91)=50";
    if Aperture > 0 then
      if isCompressibleFlow == true and PortB.P < 0.5 * PortA.P then
        0.5 * PortA.P = sign(PortA.G) * abs(1.296e9 * (abs(PortA.G) / KvFlow) ^ 2 / Rho) "This allows the calculation of G or Kv given PortA.P and PortB.P, but not the calculation of PortB.P";
      elseif isCompressibleFlow == true and PortA.P < 0.5 * PortB.P then
        0.5 * PortB.P = sign(PortB.G) * abs(1.296e9 * (abs(PortA.G) / KvFlow) ^ 2 / Rho) "if the flow is reversed";
      else
        Pdiff + (PortB.Elevation - PortA.Elevation) * Rho * g_n = -sign(PortA.G) * abs(1.296e9 * (abs(PortA.G) / KvFlow) ^ 2 / Rho) "1.296e9=3600^2*100";
      end if;
    else
      PortA.G = 0;
    end if;
    annotation(
      Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)),
      defaultComponentName = "FV",
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Text(lineColor = {0, 0, 255}, extent = {{-154, -52}, {150, -120}}, textString = "%name"), Polygon(lineColor = {0, 57, 172}, fillColor = {85, 170, 255}, fillPattern = FillPattern.HorizontalCylinder, points = {{-90, 10}, {-70, 10}, {-70, 60}, {0, 0}, {70, 60}, {70, 10}, {90, 10}, {90, -10}, {70, -10}, {70, -60}, {0, 0}, {-70, -60}, {-70, -10}, {-90, -10}, {-90, 10}}), Line(origin = {0.328629, 1.31454}, points = {{0, 86}, {0, 0}}, color = {0, 0, 127}), Text(origin = {40, 121}, extent = {{-38, 25}, {110, -15}}, textString = "Opening")}),
  Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-3210\">It is the base model for control valves, and extends the ValveBase model. It adds the following parameters for configuration:</div>

<ul class=\"itemize\" id=\"magicparlabel-3211\"><li class=\"itemize_item\">parameter Boolean isCompressibleFlow = false. In order to specify the type of flow to consider. If it is made true the flow will be always considered as adiabatic.</li>
<li class=\"itemize_item\">parameter Boolean isLinear = true. To specify the valve characteristics. If made false, the valve will be considered as isopercentual.</li>
<li class=\"itemize_item\">parameter Boolean useFixedAperture = true. If true, the aperture of the valve will be fixed to the value entered manually at the aperture parameter. Otherwise the aperture will be taken from the RealInput connector Opening, making possible the calculation of the aperture outside of the model.</li>
<li class=\"itemize_item\">parameter Real aperture = 1.0. The fraction of aperture of the valve used if useFixedAperture = true.</li>
</ul>
<div class=\"standard\" id=\"magicparlabel-3215\">The following variables are also defined:</div>

<ul class=\"itemize\" id=\"magicparlabel-3216\"><li class=\"itemize_item\">SI.Area KvFlow. The acting Kv of the valve taking into account its aperture, different from the nominal Kv.</li>
<li class=\"itemize_item\">The conditional connector Modelica.Blocks.Interfaces.RealInput Opening, conditioned to useFixedAperture=false. In order to receive the valve aperture by connexion.</li>
<li class=\"itemize_item\">The protected connector Modelica.Blocks.Interfaces.RealInput Aperture.</li>
</ul>
<div class=\"standard\" id=\"magicparlabel-3219\">The implementation of the Aperture connector is done in the standard way: The connectors Opening and Aperture are connected, but if useFixedAperture is made true, the Opening connector disappears, and also the connection. And the equation Aperture=aperture (reading the manual entered value) is activated.</div>

<div class=\"standard\" id=\"magicparlabel-3220\">There is also an equation relating the Kv, KvFlow and Aperture, taking into account the selected characteristic of the valve.</div>

<div class=\"standard\" id=\"magicparlabel-3221\">And an equation relating flow with KvFlow and differential pressure, if Aperture&gt;0. Otherwise the flow is made equal to 0.</div></body></html>"));
  end ValvePartial;

  model ValveIncompressible "Control valve model for incompresible liquid flow"
    extends ValvePartial(final isCompressibleFlow = false);
    parameter Boolean isIsenthalpicFlow = true "if false, adiabatic flow applies" annotation(
      Dialog(tab = "Flow"));
  equation
    State = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    if calcEnthalpyDifference == true then
      if isIsenthalpicFlow == true then
        PortA.H = PortB.H "this means also adiabatic, incompressible flow, if no height change";
      else
        0 = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n "adiabatic incompressible";
      end if;
    end if;
  annotation(defaultComponentName = "Valve",
      Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-3276\">Extends the ValvePartial model, fixing isCompressibleFlow=false.</div>

<div class=\"standard\" id=\"magicparlabel-3277\">It adds the parameter Boolean isIsenthalpicFlow = true, in orther to specify the type of flow to consider. If made equal to false, adiabatic flow is considered.</div>

<div class=\"standard\" id=\"magicparlabel-3278\">It defines the State (for physical properties calculation) from the PortA variables.</div>

<div class=\"standard\" id=\"magicparlabel-3279\">An equation is added relating the two ports enthalpy, if calcEnthalpyDifference=true. If the flow is isenthalpic, they are the same. if it is not, only the differential height between ports is taken into account.</div></body></html>"));
  end ValveIncompressible;

  model ValveCompressible "Control valve model for compresible liquid flow"
    extends ValvePartial(final isCompressibleFlow = true);
    parameter SI.Length di = 0.0 "if >0, kinetic energy is taken into account"annotation(
      Dialog(tab = "Physical data"));
    Medium.ThermodynamicState StateA "thermodynamic state at PortA";
    Medium.ThermodynamicState StateB "thermodynamic state at PortB";
    Medium.Density RhoA(displayUnit = "kg/m3") "density at PortA";
    Medium.Density RhoB(displayUnit = "kg/m3") "density at PortB";
    Modelica.Units.SI.Velocity Va(start = 1);
    Modelica.Units.SI.Velocity Vb(start = -1);
    Medium.Temperature Ta(displayUnit = "degC") "Temperature at PortA";
    Medium.Temperature Tb(displayUnit = "degC") "Temperature at PortB";
  equation
    if PortA.P >= PortB.P then
      State = Medium.setState_phX(max(PortB.P, PortA.P / 2),PortA.H, PortA.X) "for compressible flow max. discharge is 0.5*PortA.P";
    else
      State = Medium.setState_phX(max(PortA.P, PortB.P / 2), PortB.H, PortA.X);
    end if;
    StateA = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    StateB = Medium.setState_phX(PortB.P, PortB.H, PortB.X);
    Ta = Medium.temperature(StateA);
    Tb = Medium.temperature(StateB);
    RhoA = Medium.density(StateA);
    RhoB = Medium.density(StateB);
    if di>0 then
      Va = PortA.G * 4 / (RhoA * pi * di ^ 2);
      Vb = PortB.G * 4 / (RhoB * pi * di ^ 2);
      if calcEnthalpyDifference == true then
        0 = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n + 0.5 * (abs(Vb) ^ 2 - abs(Va) ^ 2) "energy conservation for adiabatic compressible";
      end if;
    else
      Va=0;
      Vb=0;
      if calcEnthalpyDifference == true then
        0 = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n "energy conservation for adiabatic incompressible";
      end if;
    end if;
  annotation(defaultComponentName = "Valve",
      Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-3310\">Extends the ValvePartial model, fixing isCompressibleFlow=true.</div>

<div class=\"standard\" id=\"magicparlabel-3311\">It adds two Medium.ThermodynamicState StateA and StateB. In order to retrieve the physical properties at both ports. It defines also variables for temperature, density and velocity at both ports: Ta, RhoA,Va, Tb, RhoB, Vb.</div>

<div class=\"standard\" id=\"magicparlabel-3312\">The State variable (used for retrieving physical properties for the pressure drop calculation) is defined at the highest enthalpy of the two ports and at the lowest pressure of the ports (but not lowest than half of the highest pressure).</div>

<div class=\"standard\" id=\"magicparlabel-3313\">An equation is added relating the two ports enthalpy calcEnthalpyDifference=true. It will take into account the kinetic energy only if the valve diameter is higher than 0.</div></body></html>"));
  end ValveCompressible;

  model CheckValve
    extends ValveBase(useElevDifference = true, final elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
  equation
    State = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    if noEvent(PortA.P > PortB.P) then
      Pdiff + Hdiff * Rho * g_n = -sign(PortA.G) * abs(1.296e9 * (abs(PortA.G) / Kv) ^ 2 / Rho) "1.296e9=3600^2*100";
    else
      PortA.G = 0;
    end if;
    if calcEnthalpyDifference == true then
      PortA.H = PortB.H "this means adiabatic, incompressible flow, as there is no height change";
    end if;
    annotation(
      defaultComponentName = "CV",
      Icon(coordinateSystem(initialScale = 0.04), graphics = {Line(origin = {60, 0}, points = {{0, 64}, {0, -64}}, thickness = 5), Polygon(fillColor = {1, 115, 255}, fillPattern = FillPattern.Solid, points = {{-60, 70}, {60, 0}, {-60, -70}, {-60, 70}})}));
  end CheckValve;

  model SafetyValve "Single phase flow calculation using Medium's physical properties, and adiabatic, isentropic, flow"
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Modelica.Units.SI.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5;
    parameter Boolean balancedValve = false "If balanced, Kb must be applied";
    parameter Boolean useFixedArea = true "if false, the flow will be fixed";
    parameter Modelica.Units.SI.Area fixedArea = 0.0 "area to use if useFixedArea=true";
    parameter Modelica.Units.SI.MassFlowRate fixedFlow(displayUnit = "kg/h") "flow to use if useFixedArea=false";
    parameter Real kd = 0.975 "discharge coefficient. Default gas:0.975, rup.disk: 0.62, liquid=0.65, biphasic=0.85";
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk";
    Modelica.Units.SI.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Modelica.Units.SI.Area Aeff(displayUnit = "cm2", start = 0.00001) "Effective orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.Units.SI.Density RhoA(displayUnit = "kg/m3");
    Modelica.Units.SI.DynamicViscosity MuA;
    Modelica.Units.SI.AbsolutePressure Pc(displayUnit = "bar", start = 1.0e5) "orifice outlet pressure for critical flow";
    Medium.ThermodynamicState StateC "orifice thermodynamic state at critical flow";
    Modelica.Units.SI.SpecificEnthalpy Hc "orifice specific enthalpy at critical flow";
    Modelica.Units.SI.Density RhoC(displayUnit = "kg/m3") "orifice density at critical flow";
    Modelica.Units.SI.Velocity Vc(start = 340) "velocity of sound at orifice discharge";
    Modelica.Units.SI.MassFlowRate Gc(displayUnit = "kg/h") "critical flow";
    Medium.ThermodynamicState StateB "thermodynamic state at orifice at PortB pressure";
    Modelica.Units.SI.Density RhoB(displayUnit = "kg/m3") "orifice density at PortB pressure";
    Modelica.Units.SI.SpecificEnthalpy Hb "orifice specific enthalpy at P";
    Modelica.Units.SI.Velocity Vb(start = 340) "velocity at orifice discharge";
    Modelica.Units.SI.MassFlowRate Gb(displayUnit = "kg/h");
    Modelica.Units.SI.ReynoldsNumber Re "Reynolds number at orifice";
    Real Kb "correction coefficient for back pressure";
    Real Kv "viscosity correction factor for liquids";
  algorithm
    StateA := Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    RhoA := Medium.density(StateA);
    MuA := Medium.dynamicViscosity(StateA);
    if RhoA > 500 then
      if balancedValve == true and (PortB.P - 101325) / (pSet - 101325) > 0.169 then
        Kb := 1.0 - 0.988 * ((PortB.P - 101325) / (pSet - 101325) - 0.169);
      else
        Kb := 1.0;
      end if;
      if PortA.P >= pSet then
        Re := abs(PortA.G * (4 / (A * pi)) ^ 0.5 / MuA);
        Kv := 1.0 / (0.9935 + 2.878 / Re ^ 0.5 + 342.75 / Re ^ 1.5);
      else
        Re := 0.0;
        Kv := 1.0;
      end if;
    else
      if balancedValve == true and (PortB.P - 101325) / (pSet - 101325) > 0.3 then
        Kb := 1.0 - 0.31 / 0.2 * ((PortB.P - 101325) / (pSet - 101325) - 0.3);
      else
        Kb := 1.0;
      end if;
      Re := 0.0;
      Kv := 1.0;
    end if;
  equation
    Aeff = A * kd * Kb * Kv * kc;
    StateC = Medium.setState_psX(Pc, Medium.specificEntropy(StateA), PortA.X);
    RhoC = Medium.density(StateC);
    Hc = Medium.specificEnthalpy(StateC);
    Vc = Medium.velocityOfSound(StateC);
    if RhoA < 500 then
      Hc = PortA.H - 0.5 * Vc ^ 2;
      Gc = Aeff * Vc * RhoC;
    else
      Pc = 1.0;
      Gc = 0;
    end if;
    StateB = Medium.setState_psX(PortB.P, Medium.specificEntropy(StateA), PortA.X);
    RhoB = Medium.density(StateB);
    Hb = Medium.specificEnthalpy(StateB);
    Hb = PortA.H - 0.5 * Vb ^ 2;
    Gb = Aeff * Vb * RhoB;
    if PortA.P < pSet then
      PortA.G = 0.0;
      A = 0.001;
    else
      if useFixedArea == true then
        A = fixedArea;
      else
        PortA.G = fixedFlow;
      end if;
      if PortB.P >= Pc then
        PortA.G = Gb;
      else
        PortA.G = Gc;
      end if;
    end if;
    if calcEnthalpyDifference == true then
      PortB.H = PortA.H;
    end if;
    annotation(
      Diagram(coordinateSystem(initialScale = 0.1)),
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-1980\">Allows the calculation of the flow or the area of a safety relief valve or disk. It uses isentropic calculation between the inlet, supposed with velocity equal to 0, and the orifice. The gas flow can be chocked or not, and all the calculations are done with the physical properties calculated by the Medium, except the possibility of manual entering of the isentropic coefficient. Correction coefficients are applied for backpressure and for viscosity, according to the API 520 methodology. It is capable for biphasic flow calculation, but the model SafetyValveFlash, that conforms to the API recommenations, is probably more adequate.</div>

<div class=\"standard\" id=\"magicparlabel-1981\">Two thermodynamic states are declared at the orifice outlet, both isentropics with the inlet. One of them at the discharge pressure of the valve, the other adjusted to obtain sonic speed at the orifice. If the discharge pressure is below the obtained critical pressure, the flow is taken from the critical calculation, if not from the discharge pressure calculation.</div></body></html>"));
  end SafetyValve;

  model SafetyValveStd "Safety valve calculation according to ISO 4126-1 and API 520"
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Integer selectStd = 2 "1=ISO 4126-1, 2=API 520" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5 "set pressure" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedArea = true "if false, the flow will be fixed" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.Area fixedArea = 0.0 "area to use if useFixedArea=true" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.MassFlowRate fixedFlow(displayUnit = "kg/h") = 0.0 "flow to use if useFixedArea=false" annotation(
      Dialog(tab = "Basic"));
    parameter Real kd = 0.975 "discharge coefficient(API). Default gas:0.975, rup.disk: 0.62, liquid=0.65" annotation(
      Dialog(tab = "Basic"));
    parameter Real kdr = 0.8775 "derated discharge coefficient (ISO). Default gas:0.8775, liquid: 0.585" annotation(
      Dialog(tab = "Basic"));
    parameter Real kc = 1.0 "Only for API: 0.9 if it is a safety valve with an upstream bursting disk" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean balancedValve = false "If balanced, Kb must be applied for API critical flow" annotation(
      Dialog(tab = "Basic"));
    parameter Real kb = 1.0 "Only for API: backpressure correction factor from manufacturer" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedDensity = false "if true, density is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.Density rhoA(displayUnit = "kg/m3") = 1e3 "manually fixed density at inlet" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Real z = 1.0 "if useFixedDensity=true and rhoA=0, density is calculated from Z and MW" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Real mw = 18.015 "molar mass. Used if useFixedDensity=true and rhoA=0" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.IsentropicExponent gamma = 1.0 "isentropic coefficient" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Boolean useFixedViscosity = false "if true, viscosity is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.DynamicViscosity muA = 1e-3 "fixed viscosity at inlet" annotation(
      Dialog(tab = "User phys. prop."));
    Modelica.Units.SI.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.Units.SI.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.Units.SI.Density RhoA(displayUnit = "kg/m3");
    Modelica.Units.SI.DynamicViscosity MuA;
    Modelica.Units.SI.IsentropicExponent Gamma "isentropic exponent";
    Modelica.Units.SI.AbsolutePressure Pc(displayUnit = "bar") "discharge pressure for critical flow";
    Real Overpressure;
    Real C "isentropic exponent dependent coefficient";
    Real Kb "correction coefficient for gas back pressure";
    Real F2 "correction coefficent for subcritical gas flow (API)";
    Real Kw "correction coefficient for liquid back pressure";
    Real Kv "correction coefficient for liquid viscosity";
    Modelica.Units.SI.ReynoldsNumber Re "Reynolds number at orifice";
  algorithm
    StateA := Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    Ta := Medium.temperature(StateA);
    if useFixedDensity == true then
      if rhoA == 0.0 then
        RhoA := mw * PortA.P / (1000 * z * R * Ta);
      else
        RhoA := rhoA;
      end if;
    else
      RhoA := Medium.density(StateA);
    end if;
    Gamma := gamma;
    if useFixedViscosity == true then
      MuA := muA;
    else
      MuA := Medium.dynamicViscosity(StateA);
    end if;
    if RhoA > 500 then
      Re := abs(PortA.G * (4 / (A * pi)) ^ 0.5 / MuA);
      Kv := 1.0 / (0.9935 + 2.878 / Re ^ 0.5 + 342.75 / Re ^ 1.5);
    else
      Re := 0;
      Kv := 1.0;
    end if;
    Pc := PortA.P * (2 / (Gamma + 1)) ^ (Gamma / (Gamma - 1));
    C := (Gamma * (2 / (Gamma + 1)) ^ ((Gamma + 1) / (Gamma - 1))) ^ 0.5;
  equation
    Overpressure = (PortA.P - 101325) / (pSet - 101325) - 1.0;
    if selectStd == 1 then
      if PortB.P <= Pc or RhoA > 500 then
        Kb = 1.0;
      else
        Kb = (2 * Gamma * ((-(PortB.P / PortA.P) ^ (2 / Gamma)) + (PortB.P / PortA.P) ^ ((Gamma + 1) / Gamma)) / ((1 - Gamma) * Gamma * (2 / (Gamma + 1)) ^ ((Gamma + 1) / (Gamma - 1)))) ^ 0.5;
      end if;
      F2 = 1.0;
      Kw = 1.0;
    else
      if balancedValve == true and RhoA < 500 then
        if PortB.P > Pc then
          Kb = kb;
        elseif (PortB.P - 101325) / (pSet - 101325) > 0.3 then
          Kb = 1.0 - 0.31 / 0.2 * ((PortB.P - 101325) / (pSet - 101325) - 0.3);
        else
          Kb = 1.0;
        end if;
      else
        Kb = 1.0;
      end if;
      if balancedValve == true and RhoA >= 500 then
        if (PortB.P - 101325) / (pSet - 101325) > 0.169 then
          Kw = 1.0 - 0.988 * ((PortB.P - 101325) / (pSet - 101325) - 0.169);
        else
          Kw = 1.0;
        end if;
      else
        Kw = 1.0;
      end if;
      if PortB.P > Pc then
        F2 = (Gamma * (PortB.P / PortA.P) ^ (2 / Gamma) * (1 - (PortB.P / PortA.P) ^ ((Gamma - 1) / Gamma)) / ((Gamma - 1) * (1 - PortB.P / PortA.P))) ^ 0.5;
      else
        F2 = 1.0;
      end if;
    end if;
    if useFixedArea == true then
      A = fixedArea;
    else
      PortA.G = fixedFlow;
    end if;
    if PortA.P > pSet then
      if selectStd == 1 then
        if RhoA < 500 then
          A = PortA.G / (C * kdr * Kb * (PortA.P * RhoA) ^ 0.5) "gas all situations";
        else
          A = 0.7068 * PortA.G / (kdr * Kv * ((PortA.P - PortB.P) * RhoA) ^ 0.5) "liquid";
        end if;
      else
        if RhoA < 500 then
          if PortB.P <= Pc or balancedValve == true then
            A = PortA.G * (1 / (PortA.P * RhoA)) ^ 0.5 / (C * kd * Kb * kc) "gas critical flow";
          else
            A = 0.7068 * PortA.G / (kd * F2 * kc * ((PortA.P - PortB.P) * RhoA) ^ 0.5) "gas subcritical flow. 2.0378=17.9*3600*1e-6*(1000^2)^0.5/1000^0.5 KPa and 1000*Rho";
          end if;
        else
          A = 0.7068 * PortA.G * (1 / ((PortA.P - PortB.P) * RhoA)) ^ 0.5 / (kd * Kw * Kv * kc) "liquid. 0.7068=11.78*60000";
        end if;
      end if;
    else
      PortA.G = 0.0;
    end if;
    if calcEnthalpyDifference == true then
      PortB.H = PortA.H;
    end if;
    annotation(
      Diagram(coordinateSystem(initialScale = 0.1)),
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3), Text(origin = {0, -126}, lineColor = {0, 0, 255}, extent = {{-120, 26}, {120, -26}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-1991\">API 520 and ISO 4126-1 calculation methods are implemented, with the API one as default. It is limited to monophasic flow. Although a Medium is declared, and physical properties can be taken from it (except the isentropic coefficient for gases, that must be always entered by hand), there is the possibility of manual entering of all the physical properties. If you enter all physical properties manually, you can use data for a different fluid than that of the Medium. The Medium will be used just for generating P and H at the PortA, with H being passed to the PortB if the calculation has been activated. As H is not used in the calculation of the flow, if the pressures are OK the calculation should be fine.</div></body></html>"));
  end SafetyValveStd;

  model SafetyValveFlash "Safety valve calculation for flashing liquids, according to API standard 520 Annex C.2.1"
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Modelica.Units.SI.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5;
    parameter Boolean balancedValve = false "If balanced, Kb must be applied";
    parameter Boolean useFixedArea = true "if false, the flow will be fixed";
    parameter Modelica.Units.SI.Area fixedArea = 0.0 "area to use if useFixedArea=true";
    parameter Modelica.Units.SI.MassFlowRate fixedFlow(displayUnit = "kg/h") "flow to use if useFixedArea=false";
    parameter Real kd = 0.975 "discharge coefficient. Default gas:0.975, rup.disk: 0.62, liquid=0.65, biphasic=0.85";
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk";
    Modelica.Units.SI.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Modelica.Units.SI.Area Aeff(displayUnit = "cm2", start = 0.00001) "Effective orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.Units.SI.Density RhoA(displayUnit = "kg/m3") "inlet density";
    Modelica.Units.SI.DynamicViscosity MuA "inlet viscosity";
    Modelica.Units.SI.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.Units.SI.SpecificEnthalpy Ha;
    Real GFa "gas fraction at inlet";
    Real GFc[25] "orifice gas fractions at testing pressures";
    Modelica.Units.SI.AbsolutePressure Pc[25] "orifice outlet pressures to find for critical flow";
    Modelica.Units.SI.AbsolutePressure Pmax(displayUnit = "bar") "orifice outlet pressure for critical flow";
    Medium.ThermodynamicState StateC "orifice thermodynamic state";
    Modelica.Units.SI.SpecificEnthalpy Hc[25] "orifice specific enthalpy at testing pressures";
    Modelica.Units.SI.Density RhoC[25] "orifice density at testing pressures";
    Modelica.Units.SI.Velocity Vc[25] "orifice velocity at testing pressures";
    Real Gc[25] "orifice mass velocity at testing pressures";
    Real Gmax(displayUnit = "kg/h") "critical mass velocity";
    Modelica.Units.SI.ReynoldsNumber Re "Reynolds number at orifice";
    Real Kb "correction coefficient for back pressure";
    Real Kv "viscosity correction factor for liquids";
  algorithm
    StateA := Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    RhoA := Medium.density(StateA);
    MuA := Medium.dynamicViscosity(StateA);
    Ta := Medium.temperature(StateA);
    GFa := Medium.vapourQuality(StateA);
    Ha := Medium.specificEnthalpy(Medium.setState_psX(PortA.P, Medium.specificEntropy(StateA), PortA.X));
    if RhoA > 500 then
      if balancedValve == true and (PortB.P - 101325) / (pSet - 101325) > 0.169 then
        Kb := 1.0 - 0.988 * ((PortB.P - 101325) / (pSet - 101325) - 0.169);
      else
        Kb := 1.0;
      end if;
      if PortA.P >= pSet then
        Re := abs(PortA.G * (4 / (A * pi)) ^ 0.5 / MuA);
        Kv := 1.0 / (0.9935 + 2.878 / Re ^ 0.5 + 342.75 / Re ^ 1.5);
      else
        Re := 0.0;
        Kv := 1.0;
      end if;
    else
      if balancedValve == true and (PortB.P - 101325) / (pSet - 101325) > 0.3 then
        Kb := 1.0 - 0.31 / 0.2 * ((PortB.P - 101325) / (pSet - 101325) - 0.3);
      else
        Kb := 1.0;
      end if;
      Re := 0;
      Kv := 1.0;
    end if;
    for i in 1:25 loop
      Pc[i] := PortA.P - 0.04 * i * (PortA.P - PortB.P);
      StateC := Medium.setState_psX(Pc[i], Medium.specificEntropy(StateA), PortA.X);
      RhoC[i] := Medium.density(StateC);
      Hc[i] := Medium.specificEnthalpy(StateC);
      Vc[i] := abs(2 * (Ha - Hc[i])) ^ 0.5;
      Gc[i] := Vc[i] * RhoC[i];
      GFc[i] := Medium.vapourQuality(StateC);
    end for;
    Gmax := 0.0;
    for i in 1:25 loop
      if Gmax < Gc[i] then
        Gmax := Gc[i];
        Pmax := Pc[i];
      end if;
    end for;
  equation
    if PortA.P < pSet then
      A = 0.0;
      PortA.G = 0.0;
    else
      PortA.G = Aeff * Gmax;
      if useFixedArea == true then
        A = fixedArea;
      else
        PortA.G = fixedFlow;
      end if;
    end if;
    Aeff = A * kd * Kv * kc;
    if calcEnthalpyDifference == true then
      PortB.H = PortA.H;
    end if;
    annotation(
      Diagram(coordinateSystem(initialScale = 0.1)),
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-3135\">It implements the calculation methodology described at API 520 Annex C section C.2.1, based on isentropic flow. In order to get the flow, an array of 25 pressures, between inlet and outlet pressures, are tested, and the higher flow is taken. Instead of implementing the integration of <math xmlns=\"http://www.w3.org/1998/Math/MathML\">
 <mrow>
  <mfrac>
   <mrow>
    <mrow><mi>d</mi><mi>P</mi>
    </mrow>
   </mrow>
   <mrow><mi> ρ </mi>
   </mrow>
  </mfrac>
 </mrow></math> to obtain the velocity, the relationship <math xmlns=\"http://www.w3.org/1998/Math/MathML\">
 <mrow>
  <mrow>
   <mfrac>
    <mrow>
     <mrow><mi>d</mi><mo>(</mo>
      <msup>
       <mrow><mi>v</mi>
       </mrow>
       <mrow>
        <mrow><mn>2</mn><mo>)</mo>
        </mrow>
       </mrow>
      </msup>
     </mrow>
    </mrow>
    <mrow><mn>2</mn>
    </mrow>
   </mfrac><mo>=</mo><mo>-</mo><mi>d</mi><mi>H</mi>
  </mrow>
 </mrow></math> has been used.&nbsp;</div></body></html>"));
  end SafetyValveFlash;

  model SafetyValveOmega "Safety valve calculation for flashing liquids, according to API standard 520 Annex C.2.2-3"
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Boolean liquidInlet = true "if false, a biphasic inlet calculation is used" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5 "set pressure" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedArea = true "if false, the flow will be fixed" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.Area fixedArea = 0.0 "area to use if useFixedArea=true" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.Units.SI.MassFlowRate fixedFlow(displayUnit = "kg/h") = 0.0 "flow to use if useFixedArea=false" annotation(
      Dialog(tab = "Basic"));
    parameter Real kd = 0.85 "discharge coefficient(API). Default two phases: 0.85, gas:0.975, rup.disk: 0.62, liquid=0.65" annotation(
      Dialog(tab = "Basic"));
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean balancedValve = false "If balanced, Kb must be applied for API critical flow" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedDensities = false "if true, density is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.AbsolutePressure pSatur(displayUnit = "bar") = 0.0 "manually fixed saturation pressure of the inlet. Only if liquidInlet=true" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.Density rhoA(displayUnit = "kg/m3") = 1e3 "manually fixed density at inlet" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.Density rhoA9(displayUnit = "kg/m3") = 1e3 "manually fixed density at 0.9 of inlet pressure/saturation pressure" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Boolean useFixedViscosity = false "if true, viscosity is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "User phys. prop."));
    parameter Modelica.Units.SI.DynamicViscosity muA = 1e-3 "fixed viscosity at inlet" annotation(
      Dialog(tab = "User phys. prop."));
    Modelica.Units.SI.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.Units.SI.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.Units.SI.AbsolutePressure Ps(displayUnit = "bar") "saturation pressure of the inlet";
    Real Nst " transition saturation pressure ratio";
    Real Ns "saturation pressure ratio: Ps/P inlet";
    Modelica.Units.SI.Density RhoA(displayUnit = "kg/m3") "inlet density";
    Modelica.Units.SI.Density RhoA9(displayUnit = "kg/m3") "inlet density at 0.9 of inlet pressure/saturation pressure";
    Modelica.Units.SI.DynamicViscosity MuA "inlet viscosity";
    Real Omega;
    Real Nc "critical ratio of pressures";
    Modelica.Units.SI.AbsolutePressure Pc(displayUnit = "bar") "discharge pressure for critical flow";
    Real Overpressure;
    Real Kb "correction coefficient for gas back pressure";
    Real Kw "correction coefficient for liquid back pressure";
    Real Kv "correction coefficient for liquid viscosity";
    Modelica.Units.SI.ReynoldsNumber Re "Reynolds number at orifice";
  algorithm
    StateA := Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    Ta := Medium.temperature(StateA);
    if useFixedDensities == true then
      Ps := pSatur;
      RhoA := rhoA;
      RhoA9 := rhoA9;
    else
      RhoA := Medium.density(StateA);
      if liquidInlet == true then
        Ps := Medium.saturationPressure(Medium.temperature(StateA));
        RhoA9 := Medium.density(Medium.setState_psX(0.9 * Ps, Medium.specificEntropy(StateA), PortA.X));
      else
        Ps := 0.0;
        RhoA9 := Medium.density(Medium.setState_psX(0.9 * PortA.P, Medium.specificEntropy(StateA), PortA.X));
      end if;
    end if;
    Omega := 9 * (RhoA / RhoA9 - 1);
    if liquidInlet == true then
      Nst := 2 * Omega / (1 + 2 * Omega);
      Ns := Ps / PortA.P;
      if Ns > Nst then
        Nc := Ns * 2 * Omega * (1 - (1 - (2 * Omega - 1) / (Ns * 2 * Omega)) ^ 0.5) / (2 * Omega - 1);
      else
        Nc := Ns;
      end if;
    else
      Nst := 0.0;
      Ns := 0.0;
      Nc := (1 + (1.0446 - 0.0093431 * Omega ^ 0.5) * Omega ^ (-0.56261)) ^ ((-0.70356) + 0.014685 * log(Omega));
    end if;
    Pc := Nc * PortA.P;
    if useFixedViscosity == true then
      MuA := muA;
    else
      MuA := Medium.dynamicViscosity(StateA);
    end if;
    if RhoA > 500 then
      Re := abs(PortA.G * (4 / (A * pi)) ^ 0.5 / MuA);
      Kv := 1.0 / (0.9935 + 2.878 / Re ^ 0.5 + 342.75 / Re ^ 1.5);
    else
      Re := 0;
      Kv := 1.0;
    end if;
  equation
    Overpressure = (PortA.P - 101325) / (pSet - 101325) - 1.0;
    if balancedValve == true and RhoA < 500 then
      if (PortB.P - 101325) / (pSet - 101325) > 0.3 then
        Kb = 1.0 - 0.31 / 0.2 * ((PortB.P - 101325) / (pSet - 101325) - 0.3);
      else
        Kb = 1.0;
      end if;
    else
      Kb = 1.0;
    end if;
    if balancedValve == true and RhoA >= 500 then
      if (PortB.P - 101325) / (pSet - 101325) > 0.169 then
        Kw = 1.0 - 0.988 * ((PortB.P - 101325) / (pSet - 101325) - 0.169);
      else
        Kw = 1.0;
      end if;
    else
      Kw = 1.0;
    end if;
    if PortA.P >= pSet then
      if useFixedArea == true then
        A = fixedArea;
      else
        PortA.G = fixedFlow;
      end if;
      if liquidInlet == true then
        if Ps < Nst * PortA.P then
//high subcooling region
          if Ps < PortB.P then
            PortA.G / (A * kd * Kw * Kv * kc) = 1.414 * (RhoA * (PortA.P - PortB.P)) ^ 0.5 "subcritical all liquid";
          else
            PortA.G / (A * kd * Kw * Kv * kc) = 1.414 * (RhoA * (PortA.P - Ps)) ^ 0.5 "critical";
          end if;
        else
//low subcooling region
          if Pc < PortB.P then
            PortA.G / (A * kd * Kw * Kv * kc) = (2 * (1 - Ns + Omega * Ns * log(Ns * PortA.P / PortB.P) - (Omega - 1) * (Ns - PortB.P / PortA.P)) * (PortA.P * RhoA)) ^ 0.5 / (Omega * (Ns * PortA.P / PortB.P - 1) + 1) "subcritical";
          else
            PortA.G / (A * kd * Kw * Kv * kc) = (2 * (1 - Ns + Omega * Ns * log(Ns / Nc) - (Omega - 1) * (Ns - Nc)) * (PortA.P * RhoA)) ^ 0.5 / (Omega * (Ns / Nc - 1) + 1) "critical";
          end if;
        end if;
      else
        if PortB.P < Pc then
          PortA.G / (A * kd * Kb * Kv * kc) = Nc * (PortA.P * RhoA / Omega) ^ 0.5 "critical flow";
        else
          PortA.G / (A * kd * Kb * Kv * kc) = (-2 * (Omega * log(PortA.P / PortB.P) + (Omega - 1) * (1 - PortA.P / PortB.P)) * PortA.P * RhoA) ^ 0.5 / (Omega * (PortA.P / PortB.P - 1) + 1) "subcritical flow";
        end if;
      end if;
    else
      A = fixedArea;
      PortA.G = 0.0;
    end if;
    if calcEnthalpyDifference == true then
      PortB.H = PortA.H;
    end if;
    annotation(
      Diagram(coordinateSystem(initialScale = 0.1)),
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->









<div class=\"standard\" id=\"magicparlabel-3150\">It is a model for safety valves working with flashing liquids, according to API 520 Annex C sections C.2.2 and C.2.3. The parameter Boolean liquidInlet allows you so select the calculation method to use, if made true it will use the C.2.3 methodology, that works with an inlet 100% in liquid state. Otherwise the C.2.2 method will be used. The parameters allow you to choose between fixing the discharge area or the flow, configure the basic characteristics of the valve, and enter the physical properties or let them to be taken from the Medium.&nbsp;</div></body></html>"));
  end SafetyValveOmega;

  
  model PressureVacuumValve
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium=Modelica.Media.Air.DryAirNasa, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Modelica.Units.SI.AbsolutePressure pressureSet(displayUnit="bar")=102325;
    parameter Fraction overPressure=0.1 "as fraction of set pressure";
    parameter Modelica.Units.SI.AbsolutePressure pressureHysteresis(displayUnit="mbar")=500 "hysteresis for revert overpressure valve opening";
    parameter Modelica.Units.SI.VolumeFlowRate pressureFlow(displayUnit="m3/h") "air flow at the given overpressure";
    parameter Modelica.Units.SI.AbsolutePressure vacuumSet(displayUnit="bar")=100325;
    parameter Fraction underPressure=0.1 "as fraction of set vacuum";
    parameter Modelica.Units.SI.AbsolutePressure vacuumHysteresis(displayUnit="mbar")=500 "hysteresis for revert underpressure valve opening";
    parameter Modelica.Units.SI.VolumeFlowRate vacuumFlow(displayUnit="m3/h") "air flow at the given underpressure";
  
    parameter Modelica.Units.SI.Density rho(displayUnit="kg/m3")=1.2 "density used for rating";
    Boolean p;
    Boolean v;
  equation
    Hdiff=0;
    when {PortA.P>pressureSet, PortA.P<(pressureSet-pressureHysteresis)} then
      p=PortA.P>pressureSet;
    end when;
    when {(PortA.P<vacuumSet), (PortA.P>(vacuumSet+vacuumHysteresis))} then
      v= PortA.P<vacuumSet;
    end when;   
    if p then
      PortA.G= rho*pressureFlow*((PortA.P-PortB.P)/(pressureSet*(1+overPressure)-101325))^0.5;
    elseif v then
      PortA.G=-rho*vacuumFlow*((PortA.P-PortB.P)*Medium.density(Medium.setState_phX(PortB.P,PortB.H))/(vacuumSet*(1-underPressure)-101325))^0.5;
    else
      PortA.G=0;
    end if;
    annotation(
      defaultComponentName = "PSV",
      Icon(coordinateSystem(initialScale = 0.04), graphics = {Line(origin = {22.5477, 62.5318}, points = {{0, 34}, {0, -40}}, thickness = 5), Polygon(origin = {20, 60}, fillColor = {1, 115, 255}, fillPattern = FillPattern.Solid, points = {{-60, 40}, {0, 0}, {-60, -40}, {-60, 40}}), Line(origin = {-43.6944, -57.5479}, points = {{0, 34}, {0, -40}}, thickness = 5), Polygon(origin = {-40, -60}, rotation = 180, fillColor = {1, 115, 255}, fillPattern = FillPattern.Solid, points = {{-60, 40}, {0, 0}, {-60, -40}, {-60, 40}}), Line(origin = {-66, 33}, points = {{-26, -27}, {26, 27}, {26, 27}}, thickness = 0.5), Line(origin = {-69, -34}, points = {{-23, 26}, {23, -26}}, thickness = 0.5), Line(origin = {59, 34}, points = {{-33, 26}, {33, -26}}, thickness = 0.5), Line(origin = {57, -34}, points = {{37, 26}, {-37, -26}}, thickness = 0.5), Text(origin = {0, -134}, lineColor = {45, 12, 201}, extent = {{-100, 28}, {100, -28}}, textString = "%name")}),
  Documentation(info = "<html><head></head><body>The model is for a pressure/vacuum relief valve for tanks. The model performs no calculation from the media defined, and the supplied density is only used to convert from volumetric flow to mass flow. This can produce innexact results if the pressure (and so the density)is too diferent from the one used for the rating. This could be the case if the pressure in the tank is too high. For vacuum no problem should be expected, as normally we will be taking the gas from the atmosphere.<div>Elevation and specific enthalpy are passed from Port A to Port B, also when the flow is from PortB to PortA. So normally Port A should be connected to the tank, mainly for getting the correct elevation of the port, because specific enthalpy is not usable when the flow is from PortB to A.</div><div>If a more flexible model is needed, we can use separate valves, attached to separate ports of the tank. In this case the pressure correction can be applied and the specific enthalpy used.&nbsp;</div></body></html>"));
  end PressureVacuumValve;

  
  model PressureReliefValve
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium=Modelica.Media.Air.DryAirNasa, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Modelica.Units.SI.AbsolutePressure pSet(displayUnit="bar")=102325;
    parameter Fraction overPressure=0.1 "design overpressure as fraction of set pressure";
    parameter Boolean balancedValve=false "if true, a balanced or pilot governed valve is used"; 
    parameter Modelica.Units.SI.AbsolutePressure pHysteresis(displayUnit="mbar")=500 "hysteresis for revert valve opening";
    parameter Modelica.Units.SI.VolumeFlowRate flowRating(displayUnit="m3/h") "rating flow at the given overpressure";
    parameter Modelica.Units.SI.AbsolutePressure pRating(displayUnit="bar")=101325 "discharge pressure used for rating";
    parameter Modelica.Units.SI.Density rhoRating(displayUnit="kg/m3")=1.2 "density used for rating";
    Modelica.Units.SI.Density RhoA "density at PortA";
    Modelica.Units.SI.Density RhoB "density at PortB";  
    Modelica.Units.SI.VolumeFlowRate Qa(displayUnit="m3/h") "volume flow rate at Port A";
    Modelica.Units.SI.VolumeFlowRate Qb(displayUnit="m3/h") "volume flow rate at Port A";
    Boolean p;
  
  equation
    Hdiff=0;
    RhoA=Medium.density(Medium.setState_phX(PortA.P,PortA.H,PortA.X));
    RhoB=Medium.density(Medium.setState_phX(PortB.P,PortB.H,PortB.X));
    if balancedValve==true then
      when {PortA.P>pSet, PortA.P<(pSet-pHysteresis)} then
        p=PortA.P>pSet;
      end when;
    else
      when {PortA.P>(pSet+PortB.P-pRating), PortA.P<(pSet+PortB.P-pRating-pHysteresis)} then
        p=PortA.P>(pSet+PortB.P-pRating);
      end when;  
    end if; 
    if p then
      PortA.G= rhoRating*flowRating*((PortA.P-PortB.P)*RhoB/((pSet*(1+overPressure)-pRating)*rhoRating))^0.5;
    else
      PortA.G=0;
    end if;
    Qa=PortA.G/RhoA;
    Qb=PortA.G/RhoB;
    annotation(
      defaultComponentName = "PSV",
      Icon(coordinateSystem(initialScale = 0.04), graphics = {Line(origin = {17.3838, 3.43344}, points = {{0, 34}, {0, -40}}, thickness = 5), Polygon(origin = {14, 0}, fillColor = {1, 115, 255}, fillPattern = FillPattern.Solid, points = {{-60, 40}, {0, 0}, {-60, -40}, {-60, 40}}), Line(origin = {53.5492, -26.5328}, points = {{-33, 26}, {37, 26}}, thickness = 0.5), Line(origin = {-74.1557, -88.541}, points = {{-15, 88}, {27, 88}}, thickness = 0.5), Text(origin = {-2, -73}, lineColor = {12, 4, 199}, extent = {{-98, 29}, {98, -29}}, textString = "%name")}),
  Documentation(info = "<html><head></head><body>Similar to a safety valve but without the possibility of choked flow limitation.<div>The flow is calculated from a rated flow in a prescribed conditions.</div><div>The physical properties are taken from the Medium.</div></body></html>"));
  end PressureReliefValve;

  package Examples
  extends Modelica.Icons.ExamplesPackage;
    package Water1 = FreeFluids.TMedia.Fluids.Water(refState = "User", highPressure = true) "alias for TMedia water";
    package WaterS = Modelica.Media.Water.StandardWater;
    package Air2 = Modelica.Media.Air.DryAirNasa;
    package N2 = Modelica.Media.IdealGases.SingleGases.N2;
    package R134a1 = FreeFluids.TMedia.Fluids.R134A(refState = "User", reference_T = 100, highPressure = false);
    package MarlothermSH = FreeFluids.TMedia.Fluids.MarlothermSH;

    model ValveWaterTest1 "Very simple model using external connectors for flow"
      extends Modelica.Icons.Example;
      FreeFluids.Valves.ValveIncompressible FV(redeclare package Medium = Water1, Q(displayUnit = "m3/s"), aperture = 1, fix = FreeFluids.Types.ValveFixOption.fixFlow, fixedFlow(displayUnit = "kg/h") = 3.88889, isLinear = false, useFixedAperture = true) annotation(
        Placement(visible = true, transformation(origin = {-18, -1.77636e-15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, P = 160000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = Water1,Elevation = 1, P = 300000, T(displayUnit = "degC") = 298.15, externalG = false, externalP = true, externalT = true, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const1(k = 25 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {-94, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp RampP(duration = 1, height = 3e5, offset = 2e5) annotation(
        Placement(visible = true, transformation(origin = {-94, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, FV.PortA) annotation(
        Line(points = {{-40, 0}, {-25, 0}}, color = {0, 127, 255}));
      connect(const1.y, Source.Text) annotation(
        Line(points = {{-83, 54}, {-50, 54}, {-50, 11}}, color = {0, 0, 127}));
      connect(RampP.y, Source.Pext) annotation(
        Line(points = {{-82, 88}, {-44, 88}, {-44, 12}}, color = {0, 0, 127}));
      connect(FV.PortB, Sink.PortA) annotation(
        Line(points = {{-10, 0}, {26, 0}}, color = {0, 127, 255}));
      annotation(
        Documentation(info = "<html><head></head><body>
    <p>The model has the valve Kv as unknown. The inlet and outlet pressures, are fixed, as is the flow (in tab \"Flow\" of the valve).</p>
    
    </body></html>"),
        Diagram(coordinateSystem(extent = {{-120, 100}, {60, -20}})));
    end ValveWaterTest1;

    model ValveAirTest1
      extends Modelica.Icons.Example;
      FreeFluids.Valves.ValveCompressible FV(redeclare package Medium = Air2, Q(displayUnit = "m3/h"), di = 0.02, fix = FreeFluids.Types.ValveFixOption.fixKv, fixedFlow(displayUnit = "kg/h") = 0.006944444444444444, fixedKv = 2.3404) annotation(
        Placement(visible = true, transformation(origin = {-32, -1.77636e-15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink sink(redeclare package Medium = Air2, P = 499999.9999999999, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source( redeclare package Medium = Air2,Elevation = 1, P(displayUnit = "bar") , T(displayUnit = "degC") = 298.15, externalP = true) annotation(
        Placement(visible = true, transformation(origin = {-72, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Ramp(duration = 1, height = 12e5, offset = 2e5) annotation(
        Placement(visible = true, transformation(origin = {-96, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Instruments.Reader Reader1 annotation(
        Placement(visible = true, transformation(origin = {-2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
    connect(Source.PortB, FV.PortA) annotation(
        Line(points = {{-62, 0}, {-39, 0}}, color = {0, 127, 255}));
      connect(Ramp.y, Source.Pext) annotation(
        Line(points = {{-84, 52}, {-66, 52}, {-66, 12}, {-66, 12}}, color = {0, 0, 127}));
  connect(FV.PortB, Reader1.PortA) annotation(
        Line(points = {{-24, 0}, {-12, 0}}, color = {0, 127, 255}));
  connect(Reader1.PortB, sink.PortA) annotation(
        Line(points = {{8, 0}, {26, 0}}, color = {0, 127, 255}));
    annotation(
        Diagram(coordinateSystem(extent = {{-120, 80}, {60, -20}})));end ValveAirTest1;

    model SafetyValveLeser7_5_10_3 "example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by direct isentropic flow calculation"
      extends Modelica.Icons.Example;
      SafetyValve SaftValv1(A(displayUnit = "m2"), Aeff(displayUnit = "m2"), redeclare package Medium = WaterS, Pc(displayUnit = "Pa"), fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 19.44444444444444, kd = 0.84, pSet = 11040000) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 12245000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>
<p>Example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by direct isentropic flow calculation</p></body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveLeser7_5_10_3;

    model SafetyValveStdLeser7_5_10_3 "example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by standard calculation"
      extends Modelica.Icons.Example;
      FreeFluids.Valves.SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow = 19.3846, gamma = 0.966, kd = 0.84, kdr = 0.84, pSet = 11040000, rhoA = 72.02, selectStd = 1, useFixedArea = false, useFixedDensity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 12245000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveStdLeser7_5_10_3;

    model SafetyValveStdTest7_5_10_6 "example 7.5.10.6 from Leser handbook, viscous glycerine flow, solved by standard calculation"
      extends Modelica.Icons.Example;
      SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, fixedFlow = 6.3, kd = 0.65, kdr = 0.45, muA = 1.46, pSet = 1101000, rhoA = 1260, selectStd = 1, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 1201000) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveStdTest7_5_10_6;

    model SafetyValveStdAPI1 "example from API standard 520 5.6.3.2(hydrocarbon), critical flow"
      extends Modelica.Icons.Example;
      SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 6.741666666666666, gamma = 1.11, kd = 0.975, kdr = 0.45, muA = 0.001, mw = 51, pSet = 617999.9999999999, rhoA = 0, selectStd = 2, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true, z = 0.9) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 670000, T(displayUnit = "K") = 348) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example from API standard 520 5.6.3.2(hydrocarbon), critical flow.</body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveStdAPI1;

    model SafetyValveStdAPI2 "example from API standard 520 5.6.4.2(hydrocarbon), subcritical flow"
      SafetyValveStd SaftValv1( redeclare package Medium = WaterS,A(displayUnit = "m2"), balancedValve = false, fixedFlow(displayUnit = "kg/s") = 6.741666666666666, gamma = 1.11, kd = 0.975, kdr = 0.45, muA = 0.001, mw = 51, pSet = 617999.9999999999, rhoA = 0, selectStd = 2, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true, z = 0.9) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = WaterS, P = 670000, T(displayUnit = "K") = 348) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 532000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example from API standard 520 5.6.4.2 (hydrocarbon).</body></html>"));end SafetyValveStdAPI2;

    model SafetyValveStdAPI5 "example from API standard 520 5.8.2 (crude oil), liquid flow"
      FreeFluids.Valves.SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = true, fixedFlow = 102.21, kd = 0.65, kdr = 0.45, muA = 0.44, pSet = 1825000, rhoA = 900, selectStd = 2, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 102.21, redeclare package Medium = WaterS, P = 1997000) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 445999.9999999999, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example from API standard 520 5.8.2 (crude oil), liquid flow.</body></html>"));end SafetyValveStdAPI5;

    model SafetyValveTest1 "Example of flashing water using direct isentropic flow calculation."
      extends Modelica.Icons.Example;
      SafetyValve SaftValv1(redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 2.451666666666667, kd = 0.85, pSet = 399999.9999999999, useFixedArea = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(D = 50, redeclare package Medium = WaterS, T = 422.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example of flashing water using direct isentropic flow calculation.</body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveTest1;
  
    model SafetyValveFlashTest1 "example of flashing water using API Annex C 2.2.1 methodology"
      extends Modelica.Icons.Example;
      SafetyValveFlash SaftValv1(redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 2.451666666666667, kd = 0.85, pSet = 399999.9999999999, useFixedArea = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(D = 50, redeclare package Medium = WaterS, T = 422.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example of flashing water using API 520 standard Annex C 2.2.1 methodology. And using Modelica.Media.Water.StandardWater as the medium.</body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveFlashTest1;

    
    model SafetyValveFlashTest1TMedia "example of flashing water using API Annex C 2.2.1 methodology"
      extends Modelica.Icons.Example;
      SafetyValveFlash SaftValv1(redeclare package Medium = Water1, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 2.451666666666667, kd = 0.85, pSet = 399999.9999999999, useFixedArea = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(D = 50, redeclare package Medium = Water1, T = 422.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = Water1, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example of flashing water using API standard 520 AnneX C 2.2.1 methodology. And using FreeFluids.TMedia.Fluids.Water as the medium.</body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveFlashTest1TMedia;
  
    model SafetyValveOmegaTest1 "example of flashing water using API Annex C 2.2.2 methodology"
      extends Modelica.Icons.Example;
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedArea = 1298e-6, fixedFlow = 2.479166666666667, muA = 0.001, pSet = 399999.9999999999, useFixedArea = true, useFixedDensities = false, useFixedViscosity = false) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(D = 50, redeclare package Medium = WaterS, T(displayUnit = "degC") = 422.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>
<p>Example of flashing water using API 520 standard Annex C 2.2.2 methodology.</p></body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveOmegaTest1;

    model SafetyValveFlashTest2 "example of flashing water using API Annex C 2.2.3 methodology"
      SafetyValveFlash SaftValv1(redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 13.27444444444444, kd = 0.85, pSet = 990000, useFixedArea = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = WaterS, P = 999999.9999999999, T = 452.65, sourceOption = FreeFluids.Types.SourceOption.useP_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveFlashTest2;

    model SafetyValveOmegaTest2 "example of flashing water,low subcooling, critical, using API Annex C 2.2.3 methodology"
      extends Modelica.Icons.Example;
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedArea = 1298e-6, fixedFlow = 13.27444444444444, liquidInlet = true, muA = 0.001, pSet = 990000, useFixedArea = true, useFixedDensities = false, useFixedViscosity = false) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = WaterS, P = 999999.9999999999, T(displayUnit = "degC") = 452.65, sourceOption = FreeFluids.Types.SourceOption.useP_T) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example of flashing water, low subcooling, critical flow, using API Annex C.2.2.3 methodology.</body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveOmegaTest2;

    model SafetyValveOmegaC222 "Example C.2.2.2(crude oil) from API RP520 Annex C, biphasic inlet"
      extends Modelica.Icons.Example;
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 60.15555555555556, liquidInlet = false, muA = 0.001, pSet = 514700, rhoA = 1 / 0.01945, rhoA9 = 1 / 0.02265, useFixedArea = false, useFixedDensities = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 556400, T(displayUnit = "K") = 366.5) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 204700, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>
<pre style=\"margin-top: 0px; margin-bottom: 0px;\"><span style=\"font-family: 'Bitstream Vera Sans Mono'; font-size: 12px; white-space: normal;\">Example C.2.2.2 (crude oil) from API RP520 Annex C, biphasic inlet.</span><!--EndFragment--></pre></body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveOmegaC222;

    model SafetyValveOmegaC223 "Example C.2.2.3(Propane) from API standard 520 Annex C, liquid inlet"
      extends Modelica.Icons.Example;
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 378.5 / 60 * 511.3 * 1e-3, kd = 0.65, liquidInlet = true, muA = 0.001, pSatur = 741899.9999999999, pSet = 1893600, rhoA = 511.3, rhoA9 = 262.7, useFixedArea = false, useFixedDensities = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 2073300, T(displayUnit = "K") = 288.7) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 170300, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body><pre style=\"margin-top: 0px; margin-bottom: 0px;\"><span style=\"font-family: 'Bitstream Vera Sans Mono'; font-size: 12px; white-space: normal;\">Example C.2.2.3 (Propane) from API RP520 Annex C, liquid inlet.</span></pre></body></html>"),
        Diagram(coordinateSystem(extent = {{-80, 20}, {80, -20}})));end SafetyValveOmegaC223;
    
    model SVPlusPipeSteam "example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by direct isentropic flow calculation"
      extends Modelica.Icons.Example;
      FreeFluids.Valves.SafetyValveStd SaftValv1( redeclare package Medium = WaterS,fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 19.44444444444444, gamma = 0.966, kd = 0.84, pSet = 11040000) annotation(
        Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source( redeclare package Medium = WaterS, P = 12245000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
        Placement(visible = true, transformation(origin = {-96, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Pipes.PipeFlowChoked Pipe1(redeclare package Medium = WaterS, di = 0.05, fixedW(displayUnit = "kW") = 1000, lTube = 10, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
        Placement(visible = true, transformation(origin = {2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-86, 0}, {-62, 0}}, color = {0, 127, 255}));
  connect(SaftValv1.PortB, Pipe1.PortA) annotation(
        Line(points = {{-42, 0}, {-8, 0}}, color = {0, 127, 255}));
  connect(Pipe1.PortB, Sink.PortA) annotation(
        Line(points = {{12, 0}, {56, 0}}, color = {0, 127, 255}));
    annotation(
        Documentation(info = "<html><head></head><body>Example 7.5.10.3 from Leser handbook, critical flow, saturated steam, with discharge pipe.</body></html>"),
        Diagram(coordinateSystem(extent = {{-120, 20}, {80, -20}})));
    end SVPlusPipeSteam;
    
    model PresureReliefAir
      extends Modelica.Icons.Example;
      FreeFluids.Valves.PressureReliefValve PSV(redeclare package Medium = Modelica.Media.Air.DryAirNasa, flowRating = 0.04166666666666666, overPressure = 0, pRating = 99999.99999999999, pSet = 100600, rhoRating = 1.2) annotation(
        Placement(visible = true, transformation(origin = {-8, -8.88178e-16}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Modelica.Media.Air.DryAirNasa, T (displayUnit = "degC") = 293.15, externalP = true) annotation(
        Placement(visible = true, transformation(origin = {-64, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Interfaces.FlowSink Sink(redeclare package Medium = Modelica.Media.Air.DryAirNasa, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp InRamp(duration = 0.9, height = 0.3e5, offset = 9e4, startTime = 0.1)  annotation(
        Placement(visible = true, transformation(origin = {-96, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
      connect(Source.PortB, PSV.PortA) annotation(
        Line(points = {{-54, 0}, {-20, 0}}, color = {0, 127, 255}));
      connect(PSV.PortB, Sink.PortA) annotation(
        Line(points = {{4, 0}, {36, 0}}, color = {0, 127, 255}));
  connect(InRamp.y, Source.Pext) annotation(
        Line(points = {{-84, 46}, {-58, 46}, {-58, 12}}, color = {0, 0, 127}));
    annotation(
        Diagram(coordinateSystem(extent = {{-120, 60}, {60, -20}})));end PresureReliefAir;
  end Examples;

end Valves;

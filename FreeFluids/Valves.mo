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
    parameter Modelica.SIunits.Area fixedKv = 1.0 "Kv value if it is fixed" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.PressureDifference fixedDP(displayUnit = "bar") = -0.5e5 "negative if PortB.P<PortA.P" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.MassFlowRate fixedFlow(displayUnit = "kg/h") "fixed mass flow to maintain at Port A. Possitive if flow is in" annotation(
      Dialog(tab = "Flow"));
    Modelica.SIunits.Area Kv(min = 0.0001, max = 2000, start = 1);
    Modelica.SIunits.Area Cv(min = 0.0001, max = 2000, start = 1);
    Modelica.SIunits.VolumeFlowRate Q(displayUnit = "m3/h", start = 1);
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
    annotation();
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
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Text(lineColor = {0, 0, 255}, extent = {{-154, -52}, {150, -120}}, textString = "%name"), Polygon(lineColor = {0, 57, 172}, fillColor = {85, 170, 255}, fillPattern = FillPattern.HorizontalCylinder, points = {{-90, 10}, {-70, 10}, {-70, 60}, {0, 0}, {70, 60}, {70, 10}, {90, 10}, {90, -10}, {70, -10}, {70, -60}, {0, 0}, {-70, -60}, {-70, -10}, {-90, -10}, {-90, 10}}), Line(origin = {0.328629, 1.31454}, points = {{0, 86}, {0, 0}}, color = {0, 0, 127}), Text(origin = {35, 121}, extent = {{-31, 31}, {89, -19}}, textString = "Opening")}));
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
  end ValveIncompressible;

  model ValveCompressible "Control valve model for compresible liquid flow"
    extends ValvePartial(final isCompressibleFlow = true);
    parameter SI.Length di = 0.0 annotation(
      Dialog(tab = "Physical data"));
    Medium.ThermodynamicState StateA "thermodynamic state at PortA";
    Medium.ThermodynamicState StateB "thermodynamic state at PortB";
    Medium.Density RhoA(displayUnit = "kg/m3") "density at PortA";
    Medium.Density RhoB(displayUnit = "kg/m3") "density at PortB";
    Modelica.SIunits.Velocity Va(start = 1);
    Modelica.SIunits.Velocity Vb(start = -1);
    Medium.Temperature Ta(displayUnit = "degC") "Temperature at PortA";
    Medium.Temperature Tb(displayUnit = "degC") "Temperature at PortB";
  equation
    if PortA.P >= PortB.P then
      State = Medium.setState_phX(max(PortB.P, PortA.P / 2), max(PortA.H, PortB.H), PortA.X) "for compressible flow max. discharge is 0.5*PortA.P";
    else
      State = Medium.setState_phX(max(PortA.P, PortB.P / 2), max(PortA.H, PortB.H), PortA.X);
    end if;
    StateA = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    StateB = Medium.setState_phX(PortB.P, PortB.H, PortB.X);
    Ta = Medium.temperature(StateA);
    Tb = Medium.temperature(StateB);
    RhoA = Medium.density(StateA);
    RhoB = Medium.density(StateB);
    Va = PortA.G * 4 / (RhoA * pi * di ^ 2);
    Vb = PortB.G * 4 / (RhoB * pi * di ^ 2);
    if calcEnthalpyDifference == true then
      0 = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n + 0.5 * (abs(Vb) ^ 2 - abs(Va) ^ 2) "energy conservation for adiabatic compressible";
    end if;
  end ValveCompressible;

  model CheckValve
    extends ValveBase(useElevDifference = true, final elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
  equation
    State = Medium.setState_phX(PortA.P, PortA.H);
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
    parameter Modelica.SIunits.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5;
    parameter Boolean balancedValve = false "If balanced, Kb must be applied";
    parameter Boolean useFixedArea = true "if false, the flow will be fixed";
    parameter Modelica.SIunits.Area fixedArea = 0.0 "area to use if useFixedArea=true";
    parameter Modelica.SIunits.MassFlowRate fixedFlow(displayUnit = "kg/h") "flow to use if useFixedArea=false";
    parameter Real kd = 0.975 "discharge coefficient. Default gas:0.975, rup.disk: 0.62, liquid=0.65, biphasic=0.85";
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk";
    Modelica.SIunits.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Modelica.SIunits.Area Aeff(displayUnit = "cm2", start = 0.00001) "Effective orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.SIunits.Density RhoA(displayUnit = "kg/m3");
    Modelica.SIunits.DynamicViscosity MuA;
    Modelica.SIunits.AbsolutePressure Pc(displayUnit = "bar", start = 1.0e5) "orifice outlet pressure for critical flow";
    Medium.ThermodynamicState StateC "orifice thermodynamic state at critical flow";
    Modelica.SIunits.SpecificEnthalpy Hc "orifice specific enthalpy at critical flow";
    Modelica.SIunits.Density RhoC(displayUnit = "kg/m3") "orifice density at critical flow";
    Modelica.SIunits.Velocity Vc(start = 340) "velocity of sound at orifice discharge";
    Modelica.SIunits.MassFlowRate Gc(displayUnit = "kg/h") "critical flow";
    Medium.ThermodynamicState StateB "thermodynamic state at orifice at PortB pressure";
    Modelica.SIunits.Density RhoB(displayUnit = "kg/m3") "orifice density at PortB pressure";
    Modelica.SIunits.SpecificEnthalpy Hb "orifice specific enthalpy at P";
    Modelica.SIunits.Velocity Vb(start = 340) "velocity at orifice discharge";
    Modelica.SIunits.MassFlowRate Gb(displayUnit = "kg/h");
    Modelica.SIunits.ReynoldsNumber Re "Reynolds number at orifice";
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
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)));
  end SafetyValve;

  model SafetyValveStd "Safety valve calculation according to ISO 4126-1 and API 520"
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Integer selectStd = 2 "1=ISO 4126-1, 2=API 520" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5 "set pressure" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedArea = true "if false, the flow will be fixed" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.Area fixedArea = 0.0 "area to use if useFixedArea=true" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.MassFlowRate fixedFlow(displayUnit = "kg/h") = 0.0 "flow to use if useFixedArea=false" annotation(
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
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.Density rhoA(displayUnit = "kg/m3") = 1e3 "manually fixed density at inlet" annotation(
      Dialog(tab = "Physical properties"));
    parameter Real z = 1.0 "if useFixedDensity=true and rhoA=0, density is calculated from Z and MW" annotation(
      Dialog(tab = "Physical properties"));
    parameter Real mw = 18.015 "molar mass. Used if useFixedDensity=true and rhoA=0" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.IsentropicExponent gamma = 1.0 "isentropic coefficient" annotation(
      Dialog(tab = "Physical properties"));
    parameter Boolean useFixedViscosity = false "if true, viscosity is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.DynamicViscosity muA = 1e-3 "fixed viscosity at inlet" annotation(
      Dialog(tab = "Physical properties"));
    Modelica.SIunits.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.SIunits.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.SIunits.Density RhoA(displayUnit = "kg/m3");
    Modelica.SIunits.DynamicViscosity MuA;
    Modelica.SIunits.IsentropicExponent Gamma "isentropic exponent";
    Modelica.SIunits.AbsolutePressure Pc(displayUnit = "bar") "discharge pressure for critical flow";
    Real Overpressure;
    Real C "isentropic exponent dependent coefficient";
    Real Kb "correction coefficient for gas back pressure";
    Real F2 "correction coefficent for subcritical gas flow (API)";
    Real Kw "correction coefficient for liquid back pressure";
    Real Kv "correction coefficient for liquid viscosity";
    Modelica.SIunits.ReynoldsNumber Re "Reynolds number at orifice";
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
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3), Text(origin = {0, -126}, lineColor = {0, 0, 255}, extent = {{-120, 26}, {120, -26}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
  end SafetyValveStd;

  model SafetyValveFlash "Safety valve calculation for flashing liquids, according to API standard 520 Annex C.2.1"
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Modelica.SIunits.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5;
    parameter Boolean balancedValve = false "If balanced, Kb must be applied";
    parameter Boolean useFixedArea = true "if false, the flow will be fixed";
    parameter Modelica.SIunits.Area fixedArea = 0.0 "area to use if useFixedArea=true";
    parameter Modelica.SIunits.MassFlowRate fixedFlow(displayUnit = "kg/h") "flow to use if useFixedArea=false";
    parameter Real kd = 0.975 "discharge coefficient. Default gas:0.975, rup.disk: 0.62, liquid=0.65, biphasic=0.85";
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk";
    Modelica.SIunits.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Modelica.SIunits.Area Aeff(displayUnit = "cm2", start = 0.00001) "Effective orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.SIunits.Density RhoA(displayUnit = "kg/m3") "inlet density";
    Modelica.SIunits.DynamicViscosity MuA "inlet viscosity";
    Modelica.SIunits.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.SIunits.SpecificEnthalpy Ha;
    Real GFa "gas fraction at inlet";
    Real GFc[25] "orifice gas fractions at testing pressures";
    Modelica.SIunits.AbsolutePressure Pc[25] "orifice outlet pressures to find for critical flow";
    Modelica.SIunits.AbsolutePressure Pmax(displayUnit = "bar") "orifice outlet pressure for critical flow";
    Medium.ThermodynamicState StateC "orifice thermodynamic state";
    Modelica.SIunits.SpecificEnthalpy Hc[25] "orifice specific enthalpy at testing pressures";
    Modelica.SIunits.Density RhoC[25] "orifice density at testing pressures";
    Modelica.SIunits.Velocity Vc[25] "orifice velocity at testing pressures";
    Real Gc[25] "orifice mass velocity at testing pressures";
    Real Gmax(displayUnit = "kg/h") "critical mass velocity";
    Modelica.SIunits.ReynoldsNumber Re "Reynolds number at orifice";
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
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)));
  end SafetyValveFlash;

  model SafetyValveOmega "Safety valve calculation for flashing liquids, according to API standard 520 Annex C.2.2-3"
    extends FreeFluids.Interfaces.TwoFluidPorts(redeclare replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Boolean liquidInlet = true "if false, a biphasic inlet calculation is used" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.AbsolutePressure pSet(displayUnit = "bar") = 1.01325e5 "set pressure" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedArea = true "if false, the flow will be fixed" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.Area fixedArea = 0.0 "area to use if useFixedArea=true" annotation(
      Dialog(tab = "Basic"));
    parameter Modelica.SIunits.MassFlowRate fixedFlow(displayUnit = "kg/h") = 0.0 "flow to use if useFixedArea=false" annotation(
      Dialog(tab = "Basic"));
    parameter Real kd = 0.85 "discharge coefficient(API). Default two phases: 0.85, gas:0.975, rup.disk: 0.62, liquid=0.65" annotation(
      Dialog(tab = "Basic"));
    parameter Real kc = 1.0 "0.9 if it is a safety valve with an upstream bursting disk" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean balancedValve = false "If balanced, Kb must be applied for API critical flow" annotation(
      Dialog(tab = "Basic"));
    parameter Boolean useFixedDensities = false "if true, density is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.AbsolutePressure pSatur(displayUnit = "bar") = 0.0 "manually fixed saturation pressure of the inlet. Only if liquidInlet=true" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.Density rhoA(displayUnit = "kg/m3") = 1e3 "manually fixed density at inlet" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.Density rhoA9(displayUnit = "kg/m3") = 1e3 "manually fixed density at 0.9 of inlet pressure/saturation pressure" annotation(
      Dialog(tab = "Physical properties"));
    parameter Boolean useFixedViscosity = false "if true, viscosity is taken from manually fixed value. If false, from Medium calculation" annotation(
      Dialog(tab = "Physical properties"));
    parameter Modelica.SIunits.DynamicViscosity muA = 1e-3 "fixed viscosity at inlet" annotation(
      Dialog(tab = "Physical properties"));
    Modelica.SIunits.Area A(displayUnit = "cm2", start = 0.00001) "Discharge orifice area";
    Medium.ThermodynamicState StateA "thermodynamic state at inlet";
    Modelica.SIunits.Temperature Ta(displayUnit = "degC") "inlet temperature";
    Modelica.SIunits.AbsolutePressure Ps(displayUnit = "bar") "saturation pressure of the inlet";
    Real Nst " transition saturation pressure ratio";
    Real Ns "saturation pressure ratio: Ps/P inlet";
    Modelica.SIunits.Density RhoA(displayUnit = "kg/m3") "inlet density";
    Modelica.SIunits.Density RhoA9(displayUnit = "kg/m3") "inlet density at 0.9 of inlet pressure/saturation pressure";
    Modelica.SIunits.DynamicViscosity MuA "inlet viscosity";
    Real Omega;
    Real Nc "critical ratio of pressures";
    Modelica.SIunits.AbsolutePressure Pc(displayUnit = "bar") "discharge pressure for critical flow";
    Real Overpressure;
    Real Kb "correction coefficient for gas back pressure";
    Real Kw "correction coefficient for liquid back pressure";
    Real Kv "correction coefficient for liquid viscosity";
    Modelica.SIunits.ReynoldsNumber Re "Reynolds number at orifice";
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
        RhoA9 := Medium.density(Medium.setState_psX(0.9 * Ps, Medium.specificEntropy(StateA)));
      else
        Ps := 0.0;
        RhoA9 := Medium.density(Medium.setState_psX(0.9 * PortA.P, Medium.specificEntropy(StateA)));
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
      Icon(graphics = {Polygon(origin = {0, -40}, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Polygon(origin = {40, 0}, rotation = 90, fillColor = {1, 111, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {-40, -40}, {40, -40}, {0, 40}}), Line(origin = {0.32, 39.96}, points = {{-0.317267, -39.9637}, {-20.3173, -19.9637}, {19.6827, 0.0362864}, {-20.3173, 20.0363}, {19.6827, 40.0363}}, thickness = 2), Line(origin = {-50, -55}, points = {{-50, 45}, {-50, -45}, {50, -45}, {50, -25}}, thickness = 3)}, coordinateSystem(initialScale = 0.1)));
  end SafetyValveOmega;

model Dampener
  //replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
  replaceable FreeFluids.Interfaces.FluidPortA PortA annotation(
    Placement(visible = true, transformation(extent = {{-10, -110}, {10, -90}}, rotation = 0), iconTransformation(extent = {{-10, -110}, {10, -90}}, rotation = 0)));
  parameter Modelica.SIunits.Volume v0(displayUnit="l")=1e-3;
  parameter Modelica.SIunits.AbsolutePressure p0(displayUnit = "bar") = 1e5;
  parameter Real kv=1.0, rho=1000.0;
  Modelica.SIunits.AbsolutePressure P(displayUnit = "bar");
  Modelica.SIunits.Volume V(displayUnit="l", start=0);
equation
  P*(v0-V)=p0*v0;
  sign(PortA.G) * abs(1.296e9 * (abs(PortA.G) / kv) ^ 2 /rho)=PortA.P-P;
  der(V)=PortA.G*1e-3;
  annotation(
    defaultComponentName = "Volume",
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {0, 48, 144}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-90, -90}, {90, 90}}, endAngle = 360), Text(origin = {-33, 3}, extent = {{-45, 29}, {115, -33}}, textString = "Volume"), Text(origin = {54, -128}, lineColor = {0, 0, 255}, extent = {{-154, 40}, {44, -20}}, textString = "%name")}),
experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));

end Dampener;

  package Examples
    package Water1 = FreeFluids.TMedia.Fluids.Water(refState = "User", highPressure = false) "alias for TMedia water";
    package WaterExt = FreeFluids.ExternalMedia.Fluids.WaterRef(thermoModel = 1, inputChoice = "ph") "alias for ExternalMedia water";
    package WaterS = Modelica.Media.Water.StandardWater;
    package Air1 = FreeFluids.IdealGasMedia.Air;
    package Air2 = Modelica.Media.Air.DryAirNasa;
    package N2 = Modelica.Media.IdealGases.SingleGases.N2;
    package R134a1 = FreeFluids.TMedia.Fluids.R134A(refState = "User", reference_T = 100, highPressure = false);
    package MarlothermSH = FreeFluids.TMedia.Fluids.MarlothermSH;

    model ValveWaterTest1 "Very simple model using external connectors for flow"
      FreeFluids.Valves.ValveIncompressible FV(redeclare package Medium = Water1, Q(displayUnit = "m3/s"), aperture = 1, fix = FreeFluids.Types.ValveFixOption.fixFlow, fixedFlow(displayUnit = "kg/s") = 3.88889, isLinear = false, useFixedAperture = true) annotation(
        Placement(visible = true, transformation(origin = {-6, -1.77636e-15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, G = 3.88889, P = 160000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(Elevation = 1, G = 3.88889, redeclare package Medium = Water1, P = 300000, T(displayUnit = "degC") = 298.15, externalG = false, externalP = true, externalT = true, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const1(k = 25 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {-94, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp RampP(duration = 1, height = 3e5, offset = 2e5) annotation(
        Placement(visible = true, transformation(origin = {-94, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(FV.PortB, Sink.PortA) annotation(
        Line(points = {{1, 0}, {26, 0}}, color = {0, 127, 255}));
      connect(Source.PortB, FV.PortA) annotation(
        Line(points = {{-40, 0}, {-13, 0}}, color = {0, 127, 255}));
      connect(const1.y, Source.Text) annotation(
        Line(points = {{-83, 54}, {-50, 54}, {-50, 11}}, color = {0, 0, 127}));
      connect(RampP.y, Source.Pext) annotation(
        Line(points = {{-82, 88}, {-44, 88}, {-44, 12}}, color = {0, 0, 127}));
      annotation(
        Documentation(info = "<html>
    <body>
    <p>The model has the valve Kv as unknown, as flow, and the inlet and outlet pressures, are fixed. Observe that useFixedKv has been set to false in the flow tab of the valve. You can fix it to a value, but it will be necessary not to force the flow or one of the pressures.</p>
    </body>
    </html>"));
    end ValveWaterTest1;

    model ValveAirTest1
      FreeFluids.Valves.ValveCompressible FV(redeclare package Medium = Air2, Q(displayUnit = "m3/h"), di = 0.02, fix = FreeFluids.Types.ValveFixOption.fixKv, fixedFlow(displayUnit = "kg/h") = 0.006944444444444444, fixedKv = 2.3404) annotation(
        Placement(visible = true, transformation(origin = {-6, -1.77636e-15}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink sink(redeclare package Medium = Air2, P = 499999.9999999999, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 1, G = 0.0555556, redeclare package Medium = Air2, P(displayUnit = "bar") = 200000, T(displayUnit = "degC") = 298.15, externalP = true) annotation(
        Placement(visible = true, transformation(origin = {-72, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp Ramp(duration = 1, height = 10e5, offset = 2e5) annotation(
        Placement(visible = true, transformation(origin = {-96, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(FV.PortB, sink.PortA) annotation(
        Line(points = {{1, 0}, {26, 0}}, color = {0, 127, 255}));
      connect(Source.PortB, FV.PortA) annotation(
        Line(points = {{-62, 0}, {-13, 0}}, color = {0, 127, 255}));
      connect(Ramp.y, Source.Pext) annotation(
        Line(points = {{-84, 52}, {-66, 52}, {-66, 12}, {-66, 12}}, color = {0, 0, 127}));
    end ValveAirTest1;

    model SafetyValveLeser7_5_10_3 "example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by direct isentropic flow calculation"
      SafetyValve SaftValv1(A(displayUnit = "m2"), Aeff(displayUnit = "m2"), redeclare package Medium = WaterS, Pc(displayUnit = "Pa"), fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 19.44444444444444, kd = 0.84, pSet = 11040000) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 0.2777777777777778, redeclare package Medium = WaterS, P = 12245000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveLeser7_5_10_3;

    model SafetyValveStdLeser7_5_10_3 "example 7.5.10.3 from Leser handbook, critical flow saturated steam, solved by standard calculation"
      FreeFluids.Valves.SafetyValveStd SaftValv1(A(displayUnit = "m2"), kd = 0.84, kdr = 0.84, redeclare package Medium = WaterS, pSet = 11040000, fixedArea = 1298e-6, gamma = 0.966, rhoA = 72.02, selectStd = 1, useFixedDensity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 6.3, redeclare package Medium = WaterS, P = 12245000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveStdLeser7_5_10_3;

    model SafetyValveStdTest7_5_10_6 "example 7.5.10.6 from Leser handbook, viscous glycerine flow, solved by standard calculation"
      SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, fixedFlow = 6.3, kd = 0.65, kdr = 0.45, muA = 1.46, pSet = 1101000, rhoA = 1260, selectStd = 1, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 6.3, redeclare package Medium = WaterS, P = 1201000) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveStdTest7_5_10_6;

    model SafetyValveStdAPI1 "example from API standard 520 5.6.3.2(hydrocarbon), critical flow"
      SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 6.741666666666666, gamma = 1.11, kd = 0.975, kdr = 0.45, muA = 0.001, mw = 51, pSet = 617999.9999999999, rhoA = 0, selectStd = 2, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true, z = 0.9) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 6.741666666666666, redeclare package Medium = WaterS, P = 670000, T(displayUnit = "K") = 348) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveStdAPI1;

    model SafetyValveStdAPI2 "example from API standard 520 5.6.4.2(hydrocarbon), subcritical flow"
      SafetyValveStd SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow(displayUnit = "kg/s") = 6.741666666666666, gamma = 1.11, kd = 0.975, kdr = 0.45, muA = 0.001, mw = 51, pSet = 617999.9999999999, rhoA = 0, selectStd = 2, useFixedArea = false, useFixedDensity = true, useFixedViscosity = true, z = 0.9) annotation(
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
    end SafetyValveStdAPI2;

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
    end SafetyValveStdAPI5;

    model SafetyValveFlashTest1 "example of flashing water using API Annex C 2.2.1 methodology"
      SafetyValveFlash SaftValv1(A(displayUnit = "m2"), Aeff(displayUnit = "m2"), redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 2.451666666666667, kd = 0.85, pSet = 399999.9999999999, useFixedArea = true) annotation(
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
    end SafetyValveFlashTest1;

    model SafetyValveOmegaTest1 "example of flashing water using API Annex C 2.2.2 methodology"
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
    end SafetyValveOmegaTest1;

    model SafetyValveFlashTest2 "example of flashing water using API Annex C 2.2.3 methodology"
      SafetyValveFlash SaftValv1(A(displayUnit = "m2"), Aeff(displayUnit = "m2"), redeclare package Medium = WaterS, fixedArea = 1298e-6, fixedFlow(displayUnit = "kg/h") = 13.27444444444444, kd = 0.85, pSet = 990000, useFixedArea = true) annotation(
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
    end SafetyValveOmegaTest2;

    model SafetyValveOmegaC222 "Example C.2.2.2(crude oil) from API RP520 Annex C, biphasic inlet"
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 60.15555555555556, liquidInlet = false, muA = 0.001, pSet = 514700, rhoA = 1 / 0.01945, rhoA9 = 1 / 0.02265, useFixedArea = false, useFixedDensities = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 6.741666666666666, redeclare package Medium = WaterS, P = 556400, T(displayUnit = "K") = 366.5) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 204700, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveOmegaC222;

    model SafetyValveOmegaC232 "Example C.2.2.3(Propane) from API standard 520 Annex C, liquid inlet"
      SafetyValveOmega SaftValv1(A(displayUnit = "m2"), redeclare package Medium = WaterS, balancedValve = false, fixedFlow = 378.5 / 60 * 511.3 * 1e-3, kd = 0.65, liquidInlet = true, muA = 0.001, pSatur = 741899.9999999999, pSet = 1893600, rhoA = 511.3, rhoA9 = 262.7, useFixedArea = false, useFixedDensities = true, useFixedViscosity = true) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSource Source(G = 6.741666666666666, redeclare package Medium = WaterS, P = 2073300, T(displayUnit = "K") = 288.7) annotation(
        Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.FlowSink Sink(redeclare package Medium = WaterS, P = 170300, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Source.PortB, SaftValv1.PortA) annotation(
        Line(points = {{-56, 0}, {-10, 0}}, color = {0, 127, 255}));
      connect(SaftValv1.PortB, Sink.PortA) annotation(
        Line(points = {{10, 0}, {56, 0}}, color = {0, 127, 255}));
    end SafetyValveOmegaC232;
  end Examples;
end Valves;

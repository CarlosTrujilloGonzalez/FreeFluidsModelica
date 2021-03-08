within FreeFluids;

package Pumps "Pumps.mo by Carlos Trujillo
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
  //***PUMPS***
  //***********

  partial model PumpBase
    extends Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true);
    parameter Boolean forceSpeed = false "if true, the internal, or external, speed will be applied" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean useExternalSpeed = false "If true, the connector speed will be used" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.Frequency fixedSpeed = 2900.0 / 60.0 "internally fixed speed for the pump" annotation(
      Dialog(tab = "Flow"));
    parameter Integer numParallelUnits = 1 annotation(
      Dialog(tab = "Flow"));
    parameter Boolean directFlow = true "if false, reverse flow" annotation(
      Dialog(tab = "Flow"));
    Modelica.SIunits.Frequency N "Rotational speed in 1/s";
    Real Efficiency(start = 0.6) "efficiency";
    Modelica.SIunits.Power Wabs "total absorbed power";
    Modelica.SIunits.SpecificEnergy SE;
    Modelica.SIunits.VolumeFlowRate Qunit(start = 0.01, min = 0, displayUnit = "m3/h");
    Modelica.SIunits.VolumeFlowRate Qtotal(min = 0, displayUnit = "m3/h");
    Medium.ThermodynamicState State;
    Modelica.SIunits.Density Rho(start = 1000.0, displayUnit = "kg/m3");
    Medium.Temperature T(displayUnit = "degC") "Temperature at PortA";
    Modelica.SIunits.Height Dh "differential pressure as liquid height";
    Modelica.Blocks.Interfaces.RealInput Speed if useExternalSpeed == true annotation(
      Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 110}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput Ta annotation(
      Placement(visible = true, transformation(origin = {-100, 80}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {-120, 44}, extent = {{-16, -16}, {16, 16}}, rotation = 90)));
  protected
    Modelica.Blocks.Interfaces.RealInput SpeedIn;
  equation
    connect(Speed, SpeedIn);
    if useExternalSpeed == false then
      SpeedIn = fixedSpeed;
    end if;
    if forceSpeed then
      N = SpeedIn;
    end if;
    Wabs = PortA.G * Dh * g_n / Efficiency;
    SE = Wabs / abs(PortA.G);
    if calcEnthalpyDifference == true then
      Hdiff = homotopy(SE, 0);
    end if;
    Qtotal = Qunit * numParallelUnits;
    State = Medium.setState_phX(PortA.P, PortA.H, fill(0, 0));
    Rho = Medium.density(State);
    T = Medium.temperature(State);
    Ta = T;
    PortA.G = Qtotal * Rho;
    Dh = Pdiff / Rho / g_n;
    annotation(
      defaultComponentName = "P",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {170, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Polygon(lineColor = {0, 54, 162}, fillColor = {255, 255, 255}, fillPattern = FillPattern.HorizontalCylinder, points = {{-60, 68}, {90, 10}, {90, -10}, {-60, -68}, {-60, 68}}), Text(origin = {51, 106}, extent = {{-29, 30}, {43, -12}}, textString = "Speed"), Text(origin = {-159.083, 33.1429}, extent = {{-16.9167, 22.8571}, {25.0833, -9.14286}}, textString = "Ta")}));
  end PumpBase;

  partial model BumpPumpBase "General bump pump model.We have as parameters three points of the curve, includig that at 0 flow"
    extends PumpBase(final directFlow = true, fixedSpeed = n0, Qunit.start = q1, Qtotal.start = numParallelUnits * q1, Dh.start = h1, Efficiency.start = r1, PortB.G(start = -numParallelUnits * q1 * 1e-3));
    parameter Modelica.SIunits.Frequency n0 "rotational speed of the pump at which the next values are taken (1/s)" annotation(
      Dialog(tab = "Pump curve"));
    parameter Modelica.SIunits.Height h0 "pressure height at q=0" annotation(
      Dialog(tab = "Pump curve"));
    parameter Modelica.SIunits.Height h1 "pressure height at q1" annotation(
      Dialog(tab = "Pump curve"));
    parameter Modelica.SIunits.VolumeFlowRate q1(displayUnit = "m3/h") annotation(
      Dialog(tab = "Pump curve"));
    parameter Real r1 "efficiency at q1 point" annotation(
      Dialog(tab = "Pump curve"));
    parameter Modelica.SIunits.Height h2 "pressure height at q2" annotation(
      Dialog(tab = "Pump curve"));
    parameter Modelica.SIunits.VolumeFlowRate q2(displayUnit = "m3/h") annotation(
      Dialog(tab = "Pump curve"));
    parameter Real r2 "efficiency at q2 point" annotation(
      Dialog(tab = "Pump curve"));
    Real[5] coef "pump characteristic coefficients";
    Modelica.SIunits.VolumeFlowRate QunitNV(min = 0, start = q1, displayUnit = "m3/h") "pump flow in non viscous conditions";
    Modelica.SIunits.Height DhNV "differential pressure as liquid height in non viscous conditions";
    Real EfficiencyNV(start = 0.6) "efficiency in non viscous conditions";
  algorithm
    coef[1] := h0 / n0 ^ 2;
    coef[3] := (q1 * (h2 - h0) - q2 * (h1 - h0)) / (q2 ^ 2 * q1 - q1 ^ 2 * q2);
    coef[2] := (h1 - h0 - coef[3] * q1 ^ 2) / n0 / q1;
    coef[5] := (r2 / n0 ^ 0.08 - r1 * q2 / q1 / n0 ^ 0.08) / ((q2 / n0) ^ 1.5 - q2 / q1 * (q1 / n0) ^ 1.5);
    coef[4] := (r1 / n0 ^ 0.08 - coef[5] * (q1 / n0) ^ 1.5) * n0 / q1;
  equation
    DhNV = coef[1] * N ^ 2 + coef[2] * N * QunitNV + coef[3] * QunitNV ^ 2;
    if N > 0 then
      EfficiencyNV = N ^ 0.08 * (coef[4] * QunitNV / N + coef[5] * (QunitNV / N) ^ 1.5);
    else
      EfficiencyNV = 1;
    end if;
  end BumpPumpBase;

  model BumpPump
    extends BumpPumpBase(Qunit.start = q1);
  equation
    Qunit = QunitNV;
    Dh = DhNV;
    Efficiency = EfficiencyNV;
    annotation(
      Documentation(info = "<html>
        <body>
        <p>The model is ready to be connected after filling the parameters.</p>
        <p>If useHeightDifference is set to false, the height of both ports became undefined by the model, so they must be supplied by the connexions. If there are parallell paths, remember to pass only the height of one of the parallel paths, or of none and connect a height refrence.</p>
        <p>If equalH is set to false, the absorbed power will be added to the PortB enthalpy.</p>
        <p>If useFixedSpeed is set to false, you must supply value to the Speed connector.</p>
        </body>
        </html>"));
  end BumpPump;

  model BumpPumpViscous "Bump pump for viscous flow"
    extends BumpPumpBase(QunitNV.start = q1);
    Modelica.SIunits.DynamicViscosity Mu(start = 1.0);
    Modelica.SIunits.KinematicViscosity Nu(start = 0.001);
    Modelica.SIunits.VolumeFlowRate Qopt(start = q1, min = 0, displayUnit = "m3/h") "pump optimal flow";
    Possitive bhi(start = 6.0) "positive auxiliary variable";
    Fraction Fq(min = 0.01, max = 1.0, start = 0.9) "viscous correction factor for flow";
    Fraction Fh(min = 0.01, max = 1.0, start = 0.8) "viscous correction factor for pressure";
    Fraction Fr(min = 0.01, max = 1.0, start = 0.7) "viscous correction factor for efficiency";
  algorithm
    Mu := Medium.dynamicViscosity(State);
  equation
    Nu = Mu / Rho;
    0 = N ^ 0.08 * (coef[4] / N + coef[5] * 1.5 / N * (Qopt / N) ^ 0.5) "en el punto optimo la derivada del rendimiento es 0";
    bhi = 520 * Nu ^ 0.5 / QunitNV ^ 0.25 / (g_n * DhNV) ^ 0.125 "rectificado ligeramente sobre original que era 480";
    Fq = exp(-0.11 * abs(log(bhi) / log(10)) ^ 5.5);
    Qunit = Fq * QunitNV;
    Fh = (0.2 + 0.8 * Fq) * (1 - 0.014 * (bhi - 1) * (QunitNV / Qopt - 1)) "rectificado ligeramente sobre el original que era 0.25+0.75*Fq";
    Dh = Fh * DhNV;
    Fr = exp(-0.053 * exp(0.04 * (bhi - 0.5) ^ 0.5) * (bhi - 0.5) ^ 1.08) "rectificado ligeramente sobre el original que era -0.05";
    Efficiency = Fr * EfficiencyNV;
  end BumpPumpViscous;

  model PositivePump "Possitive desplacement pump with volumetric flow proportional to speed"
    extends PumpBase(directFlow = true, fixedSpeed = n0);
    parameter Modelica.SIunits.Frequency n0 "nominal rotational speed (1/s)" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.VolumeFlowRate q0(displayUnit = "m3/h") "individual pump flow at nominal speed" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.VolumeFlowRate qLeak(displayUnit = "m3/h", min = 0) = 0 "individual leak flow at given differential pressure" annotation(
      Dialog(tab = "Flow"));
    parameter Medium.AbsolutePressure pLeak(displayUnit = "bar", min = 1e-6) = 1e-3 "differential pressure at which the leak is referenced" annotation(
      Dialog(tab = "Flow"));
    parameter Fraction r = 0.7 "fixed pump efficiency" annotation(
      Dialog(tab = "Flow"));
    Modelica.SIunits.VolumeFlowRate Qleak(displayUnit = "m3/h", min = 0);
  equation
    Efficiency = r;
    Qleak = qLeak * Pdiff / pLeak;
    Qunit = q0 * N / n0 - Qleak;
    annotation(
      Documentation(info = "<html>
        <body>
        <p>The model is ready to be connected after filling the parameters.</p>
        <p>If useHeightDifference is set to false, the height of both ports became undefined by the model, so they must be supplied by the connexions. If there are parallell paths, remember to pass only the height of one of the parallel paths, or of none and connect a height refrence.</p>
        <p>If equalH is set to false, the absorbed power will be added to the PortB enthalpy.</p>
        <p>If useFixedSpeed is set to false, you must supply value to the Speed connector.</p>
        </body>
        </html>"));
  end PositivePump;

  model PistonPump "pump governed by the position of a mechanical flange"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    FreeFluids.Interfaces.FluidPortA PortA annotation(
      Placement(visible = true, transformation(origin = {-100, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB PortB annotation(
      Placement(visible = true, transformation(origin = {100, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter Modelica.SIunits.Area section(displayUnit = "mm2");
    parameter Modelica.SIunits.Volume v0(displayUnit = "l") "chamber volume at position 0 of piston. Must be higher than pulse volume";
    
    Modelica.Mechanics.Translational.Interfaces.Flange_b Flange annotation(
      Placement(visible = true, transformation(origin = {0, 94}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {1.77636e-15, 88}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.SIunits.Volume Vch (displayUnit="l") "chamber volume";
    Modelica.SIunits.Volume Vpumped (displayUnit="l", start=0) "pumped volume";
    Medium.ThermodynamicState State;
    Modelica.SIunits.Density rho;
    Modelica.SIunits.Mass Mpumped;
  
  equation
    PortB.H=PortA.H+(PortB.P-PortA.P)/rho;
    State=Medium.setState_phX(PortA.P,PortA.H);
    rho=Medium.density(State);
    Vch = v0 - section * Flange.s;
    PortA.G = max(rho * der(Vch), 0);
    PortB.G = min(rho * der(Vch), 0);
    PortB.Elevation=PortA.Elevation;
    PortB.X=PortA.X;
    Flange.f = PortB.P * section;
    der(Mpumped)=-PortB.G;
    Vpumped=Mpumped/rho;
    annotation(
      defaultComponentName = "Source",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(origin = {-84, 68}, extent = {{-56, 42}, {60, -22}}, textString = "Pump"), Rectangle(origin = {0, -49}, fillColor = {43, 111, 205}, fillPattern = FillPattern.Solid, extent = {{-32, 41}, {32, -41}}), Rectangle(origin = {-7.10543e-15, -6}, fillColor = {170, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-32, -2}, {32, 2}}), Rectangle(origin = {-0.169112, 20.2121}, extent = {{-31.8309, -24.2121}, {31.8309, 24.2121}}), Rectangle(origin = {0, 36}, fillColor = {170, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-2, -40}, {2, 40}}), Line(origin = {-61, -80}, points = {{29, 0}, {-29, 0}}, thickness = 1), Line(origin = {61, -80}, points = {{-29, 0}, {29, 0}}, thickness = 1), Polygon(origin = {65, -79}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, points = {{-13, 9}, {-13, -11}, {9, -1}, {-13, 9}})}),
      Documentation(info = "<html>
    <body>
    <p>The pump provides the function of a possitive displacement pump with a moving piston or membrane. the pump movement and speed is governed by the position of the mechanical interface. Reverse flow inhibition is incorporated inside the pump</p>
    </body>
    </html>"));
  end PistonPump;

  package Examples
    model BumpPumpTest
      FreeFluids.Pumps.BumpPump Pump(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, calcEnthalpyDifference = true, forceSpeed = true, h0 = 29.1, h1 = 27, h2 = 23, n0 = 2900 / 60, numParallelUnits = 2, q1(displayUnit = "m3/s") = 0.0166667, q2(displayUnit = "m3/s") = 0.025, r1 = 0.705, r2 = 0.73, useElevDifference = true) annotation(
        Placement(visible = true, transformation(origin = {-42, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph Pipe(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Rho(displayUnit = ""), RhoA(displayUnit = ""), Ta(displayUnit = ""), Tb(displayUnit = ""), calcEnthalpyDifference = false, di = 0.064, lTube = 100.0, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {6, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Elevation = 0, G = 15, P = 1e5, T = 523.15, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-106, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Rho(displayUnit = ""), RhoA(displayUnit = ""), Ta(displayUnit = ""), Tb(displayUnit = ""), calcEnthalpyDifference = false, di = 0.064, lTube = 200.0, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {8, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible Valve(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, fixedKv = 80.0) annotation(
        Placement(visible = true, transformation(origin = {-76, 8}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, G = 15, P = 80000, fix = FreeFluids.Types.BoundaryOption.fixFlow) annotation(
        Placement(visible = true, transformation(origin = {70, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Pipe1.PortB, Pump.PortA) annotation(
        Line(points = {{18, -24}, {18, -50}, {-52, -50}, {-52, 8}}, color = {0, 127, 255}));
      connect(Pipe1.PortB, Sink.PortA) annotation(
        Line(points = {{18, -24}, {60, -24}}, color = {0, 127, 255}));
      connect(Pipe.PortB, Sink.PortA) annotation(
        Line(points = {{16, 34}, {60, 34}, {60, -24}}, color = {0, 127, 255}));
      connect(Valve.PortB, Pump.PortA) annotation(
        Line(points = {{-70, 8}, {-52, 8}, {-52, 8}, {-52, 8}}, color = {0, 127, 255}));
      connect(Source.PortB, Valve.PortA) annotation(
        Line(points = {{-96, 8}, {-81, 8}}, color = {0, 127, 255}));
      connect(Pump.PortB, Pipe1.PortA) annotation(
        Line(points = {{-32, 8}, {-32, -25}, {-2, -25}, {-2, -24}}, color = {0, 127, 255}));
      connect(Pump.PortB, Pipe.PortA) annotation(
        Line(points = {{-32, 8}, {-32, 34}, {-4, 34}}, color = {0, 127, 255}));
    end BumpPumpTest;

    model BumpPumpViscousTest
      FreeFluids.Pumps.BumpPumpViscous Pump(redeclare package Medium = FreeFluids.TMedia.Fluids.EG, forceSpeed = true, h0 = 50, h1 = 30, h2 = 5, n0 = 1450 / 60, q1 = 3.5 / 3600, q2 = 7.5 / 3600, r1 = 0.259, r2 = 0.231) annotation(
        Placement(visible = true, transformation(origin = {-6, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 0, redeclare package Medium = FreeFluids.TMedia.Fluids.EG, P = 1e5, T = 273.15, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-80, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.TMedia.Fluids.EG, P = 247300, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {62, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Pump.Dhv=13.8;
      connect(Pump.PortB, Sink.PortA) annotation(
        Line(points = {{4, 4}, {52, 4}, {52, 4}, {52, 4}}, color = {0, 127, 255}));
      connect(Source.PortB, Pump.PortA) annotation(
        Line(points = {{-70, 4}, {-16, 4}}, color = {0, 127, 255}));
    end BumpPumpViscousTest;

    model PumpPipesTest
      FreeFluids.Pumps.BumpPump Pump(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, forceSpeed = true, h0 = 29.1, h1 = 27, h2 = 23, n0 = 2900 / 60, numParallelUnits = 1, q1 = 0.0166667, q2 = 0.025, r1 = 0.705, r2 = 0.73) annotation(
        Placement(visible = true, transformation(origin = {-78, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, di = 0.06, equivL_Di = 2 * 340 + 5 * 15, lTube = 25, thermalType = FreeFluids.Types.ThermalType.isenthalpic, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {12, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP Source( Elevation = 0, G = 15,redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, P = 150000, T = 393.15) annotation(
        Placement(visible = true, transformation(origin = {-127, 61}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      FreeFluids.Pipes.MixerPH MixerIn(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH) annotation(
        Placement(visible = true, transformation(origin = {-106, 12}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, P = 120000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {127, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV10(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, fixedKv = 80) annotation(
        Placement(visible = true, transformation(origin = {-99, 61}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      FreeFluids.Pipes.CoilForcedConvection Coil1(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, PLossFriction(displayUnit = "Pa"), di = 0.06, fixedW = 50e3, numActiveTubes = 2, numTubes = 2, thermalType = FreeFluids.Types.ThermalType.fixedPower, useThermalConnector = false) annotation(
        Placement(visible = true, transformation(origin = {46, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.CoilForcedConvection Coil2(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, PLossFriction(displayUnit = "Pa"), coilDiam = 2, di = 0.06, fixedW = 30e3, num = 7, numActiveTubes = 2, numTubes = 2, thermalType = FreeFluids.Types.ThermalType.fixedPower, useThermalConnector = false) annotation(
        Placement(visible = true, transformation(origin = {46, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.CoilForcedConvection Coil3(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, di = 0.05, fixedW = -80e3, numActiveTubes = 2, numTubes = 2, thermalType = FreeFluids.Types.ThermalType.fixedPower, useThermalConnector = false) annotation(
        Placement(visible = true, transformation(origin = {46, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.Mixer3PH MixerOut(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH) annotation(
        Placement(visible = true, transformation(origin = {76, 12}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV18(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Q(displayUnit = "m3/s"), fixedKv = 80) annotation(
        Placement(visible = true, transformation(origin = {-44, 12}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible FV19(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Q(displayUnit = "m3/s"), aperture = 0.001, fixedKv = 80) annotation(
        Placement(visible = true, transformation(origin = {-68, -20}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
      FreeFluids.Pipes.PipeFlow1Ph PipeCAA(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, di = 0.06, equivL_Di = 2 * 340 + 5 * 15, lTube = 25, thermalType = FreeFluids.Types.ThermalType.isenthalpic, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {-68, -52}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FreeFluids.Pipes.MixerPH MixerCAA(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH) annotation(
        Placement(visible = true, transformation(origin = {-20, 12}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    equation
      connect(Source.PortB, FV10.PortA) annotation(
        Line(points = {{-116, 61}, {-107, 61}}, color = {0, 127, 255}));
      connect(FV10.PortB, MixerIn.PortA) annotation(
        Line(points = {{-91, 61}, {-78, 61}, {-78, 33}, {-112.5, 33}, {-112.5, 15}, {-113, 15}}, color = {0, 127, 255}));
      connect(MixerIn.PortC, Pump.PortA) annotation(
        Line(points = {{-99, 12}, {-88, 12}}, color = {0, 127, 255}));
      connect(Coil1.PortB, MixerOut.PortA) annotation(
        Line(points = {{56, 34}, {68, 34}, {68, 18}, {69, 18}, {69, 17}}, color = {0, 127, 255}));
      connect(Coil2.PortB, MixerOut.PortB) annotation(
        Line(points = {{56, 12}, {69, 12}}, color = {0, 127, 255}));
      connect(Coil3.PortB, MixerOut.PortC) annotation(
        Line(points = {{56, -16}, {68, -16}, {68, 7}, {69, 7}}, color = {0, 127, 255}));
      connect(MixerOut.PortD, Sink.PortA) annotation(
        Line(points = {{83, 12}, {116, 12}, {116, 11}}, color = {0, 127, 255}));
      connect(Pipe1.PortB, Coil2.PortA) annotation(
        Line(points = {{22, 12}, {36, 12}}, color = {0, 127, 255}));
      connect(Pipe1.PortB, Coil1.PortA) annotation(
        Line(points = {{22, 12}, {36, 12}, {36, 34}}, color = {0, 127, 255}));
      connect(Pipe1.PortB, Coil3.PortA) annotation(
        Line(points = {{22, 12}, {36, 12}, {36, -16}}, color = {0, 127, 255}));
      connect(MixerOut.PortD, MixerIn.PortB) annotation(
        Line(points = {{83, 12}, {83.75, 12}, {83.75, -78}, {-111.75, -78}, {-111.75, 8}, {-112, 8}}, color = {0, 127, 255}));
      connect(Pump.PortB, FV18.PortA) annotation(
        Line(points = {{-68, 12}, {-51, 12}}, color = {0, 127, 255}));
      connect(Pump.PortB, FV19.PortA) annotation(
        Line(points = {{-68, 12}, {-68, -13}}, color = {0, 127, 255}));
      connect(FV19.PortB, PipeCAA.PortA) annotation(
        Line(points = {{-68, -27}, {-68, -42}}, color = {0, 127, 255}));
      connect(PipeCAA.PortB, MixerCAA.PortB) annotation(
        Line(points = {{-68, -62}, {-26, -62}, {-26, 8}, {-26, 8}}, color = {0, 127, 255}));
      connect(MixerCAA.PortC, Pipe1.PortA) annotation(
        Line(points = {{-12, 12}, {2, 12}, {2, 12}, {2, 12}}, color = {0, 127, 255}));
      connect(FV18.PortB, MixerCAA.PortA) annotation(
        Line(points = {{-37, 12}, {-31.5, 12}, {-31.5, 16}, {-26, 16}}, color = {0, 127, 255}));
    end PumpPipesTest;

    model PositivePumpTest
      FreeFluids.Pumps.PositivePump pump(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Qunit(displayUnit = "m3/s"), forceSpeed = true, n0 = 300 / 60, pLeak(displayUnit = ""), q0(displayUnit = "m3/s") = 0.00833333, r = 0.7) annotation(
        Placement(visible = true, transformation(origin = {-42, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph pipe(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Rho(displayUnit = ""), RhoA(displayUnit = ""), Ta(displayUnit = ""), Tb(displayUnit = ""), calcEnthalpyDifference = false, di = 0.064, lTube = 100.0, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {6, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSourceSP source(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Elevation = 0, G = 15, P = 110000, T = 513.15, isGsource = false) annotation(
        Placement(visible = true, transformation(origin = {-106, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Pipes.PipeFlow1Ph pipe1(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, Rho(displayUnit = ""), RhoA(displayUnit = ""), Ta(displayUnit = ""), Tb(displayUnit = ""), calcEnthalpyDifference = false, di = 0.064, lTube = 200.0, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
        Placement(visible = true, transformation(origin = {8, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FreeFluids.Valves.ValveIncompressible valve(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, fixedKv = 80.0) annotation(
        Placement(visible = true, transformation(origin = {-76, 8}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
      FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.TMedia.Fluids.MarlothermSH, P = 100000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
        Placement(visible = true, transformation(origin = {70, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Sink.PortA, pump.PortA) annotation(
        Line(points = {{60, 2}, {60, -42}, {-52, -42}, {-52, 8}}, color = {0, 127, 255}));
      connect(pipe1.PortB, Sink.PortA) annotation(
        Line(points = {{18, -24}, {60, -24}, {60, 2}, {60, 2}}, color = {0, 127, 255}));
      connect(pipe.PortB, Sink.PortA) annotation(
        Line(points = {{16, 34}, {60, 34}, {60, 2}, {60, 2}}, color = {0, 127, 255}));
      connect(valve.PortB, pump.PortA) annotation(
        Line(points = {{-70, 8}, {-52, 8}, {-52, 8}, {-52, 8}}, color = {0, 127, 255}));
      connect(source.PortB, valve.PortA) annotation(
        Line(points = {{-96, 8}, {-81, 8}}, color = {0, 127, 255}));
      connect(pump.PortB, pipe1.PortA) annotation(
        Line(points = {{-32, 8}, {-32, -25}, {-2, -25}, {-2, -24}}, color = {0, 127, 255}));
      connect(pump.PortB, pipe.PortA) annotation(
        Line(points = {{-32, 8}, {-32, 34}, {-4, 34}}, color = {0, 127, 255}));
    end PositivePumpTest;
    
    model PistonPumpTest
  replaceable package medium = FreeFluids.TMedia.Fluids.Water;
  constant Modelica.SIunits.Frequency freq = 46;
  parameter Modelica.SIunits.Length amplitude = 0.005567;
  Modelica.Blocks.Sources.Sine Pulse(amplitude = amplitude, freqHz = freq, offset = amplitude, phase(displayUnit = "rad")) annotation(
    Placement(visible = true, transformation(extent = {{-113.9, -2.89999}, {-90.2, 20.9}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Position Pump_Position(exact = true, useSupport = false, s(start = 0)) annotation(
    Placement(visible = true, transformation(extent = {{-53.8, -1.2}, {-27.9, 19.2}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sensors.PositionSensor Pump_PosMeasure annotation(
    Placement(visible = true, transformation(extent = {{10, 2.10001}, {34.2, 16.4}}, rotation = 0)));
  FreeFluids.Pumps.PistonPump Pump(redeclare package Medium = medium, Vpumped, section = 0.0001131, v0 = 2e-06) annotation(
    Placement(visible = true, transformation(origin = {-12, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Pipes.PipeFlow1Ph Pipe(redeclare package Medium = medium, PLossFriction, PortA(G(start = -0.001)), di = 0.005000000000000001, lTube = 0.24, roughness = 2.500000000000001e-05) annotation(
    Placement(visible = true, transformation(origin = {32, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Interfaces.FlowSourceSP Bondary1(redeclare package Medium = medium, P(displayUnit = "bar") = 99999.99999999999, T(displayUnit = "K") = 298) annotation(
    Placement(visible = true, transformation(origin = {-70, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = medium, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
    Placement(visible = true, transformation(origin = {90, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
//Connection(s)
  connect(Pulse.y, Pump_Position.s_ref) annotation(
    Line(points = {{-89, 9}, {-56, 9}}, color = {0, 0, 127}));
  connect(Pump_Position.flange, Pump_PosMeasure.flange) annotation(
    Line(points = {{-28, 9}, {10, 9}}, color = {0, 127, 0}));
  connect(Pump_Position.flange, Pump.Flange) annotation(
    Line(points = {{-28, 10}, {-12, 10}, {-12, -23}}, color = {0, 127, 0}));
  connect(Pipe.PortB, Sink.PortA) annotation(
    Line(points = {{42, -40}, {80, -40}}, color = {0, 127, 255}));
  connect(Bondary1.PortB, Pump.PortA) annotation(
    Line(points = {{-60, -40}, {-22, -40}}, color = {0, 127, 255}));
  connect(Pump.PortB, Pipe.PortA) annotation(
    Line(points = {{-2, -40}, {22, -40}}, color = {0, 127, 255}));
  annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.001));

end PistonPumpTest;
  end Examples;
end Pumps;

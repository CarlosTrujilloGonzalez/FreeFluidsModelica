within FreeFluids;

package Pipes "Pipes.mo by Carlos Trujillo
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

  model MixerPH "Inlets: ports A and B. Outlet: port C. Uses mix enthalpy. Adequate for liquids and gases"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    replaceable FreeFluids.Interfaces.FluidPortA PortA(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(transformation(extent = {{-110, 35}, {-90, 55}})));
    replaceable FreeFluids.Interfaces.FluidPortA PortB(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(transformation(extent = {{-110, -55}, {-90, -35}})));
    replaceable FreeFluids.Interfaces.FluidPortB PortC(redeclare package Medium = Medium, G(max = 0, start = -2), H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default)))) annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
  equation
    PortC.Elevation = PortA.Elevation;
    PortC.X = PortA.X;
    PortB.P = PortA.P;
    PortC.P = PortA.P;
    0 = PortA.G + PortB.G + PortC.G;
    0 = PortA.G * PortA.H + PortB.G * PortB.H + PortC.G * PortC.H;
    annotation(
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Text(origin = {30, 16}, extent = {{-150, 100}, {150, 30}}, textString = "%name"), Line(origin = {-46, -20}, points = {{-46, -20}, {46, 20}}, color = {0, 85, 255}, thickness = 1), Line(origin = {45, 0}, points = {{-45, 0}, {45, 0}, {45, 0}}, color = {0, 85, 255}, thickness = 1), Line(origin = {-45, 21}, points = {{-45, 21}, {45, -21}}, color = {0, 85, 255}, thickness = 1), Polygon(origin = {39.62, 0}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, points = {{-19.6213, 10}, {20.3787, 0}, {-19.6213, -10}, {-19.6213, 10}}), Text(origin = {-78, -66}, extent = {{-20, 20}, {20, -20}}, textString = "B"), Text(origin = {-80, 22}, extent = {{-20, 20}, {20, -20}}, textString = "A"), Text(origin = {84, -24}, extent = {{-20, 20}, {20, -20}}, textString = "C")}));
  end MixerPH;

  model Mixer3PH
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    replaceable FreeFluids.Interfaces.FluidPortA PortA(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, 57}, {-90, 77}}, rotation = 0), iconTransformation(extent = {{-110, 57}, {-90, 77}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortA PortB(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, -11}, {-90, 9}}, rotation = 0), iconTransformation(extent = {{-110, -11}, {-90, 9}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortA PortC(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, -79}, {-90, -59}}, rotation = 0), iconTransformation(extent = {{-110, -77}, {-90, -57}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortB PortD(redeclare package Medium = Medium, G(max = 0, start = -3), H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default)))) annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
  equation
    PortD.Elevation = PortA.Elevation;
    PortD.X = PortA.X;
    PortB.P = PortA.P;
    PortC.P = PortA.P;
    PortD.P = PortA.P;
    0 = PortA.G + PortB.G + PortC.G + PortD.G;
    0 = PortA.G * PortA.H + PortB.G * PortB.H + PortC.G * PortC.H + PortD.G * PortD.H;
    annotation(
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Text(origin = {44, 16}, extent = {{-150, 100}, {150, 30}}, textString = "%name"), Line(origin = {-46, -20}, points = {{-44, 20}, {46, 20}}, color = {0, 85, 255}, thickness = 1), Line(origin = {45, 0}, points = {{-45, 0}, {45, 0}, {45, 0}}, color = {0, 85, 255}, thickness = 1), Line(origin = {-45.258, 41.6388}, points = {{-45, 21}, {45, -41}}, color = {0, 85, 255}, thickness = 1), Polygon(origin = {39.62, 0}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, points = {{-19.6213, 10}, {20.3787, 0}, {-19.6213, -10}, {-19.6213, 10}}), Text(origin = {-80, 16}, extent = {{-20, 20}, {20, -20}}, textString = "B"), Text(origin = {-82, 88}, extent = {{-20, 20}, {20, -20}}, textString = "A"), Text(origin = {84, -24}, extent = {{-20, 20}, {20, -20}}, textString = "D"), Line(origin = {-41.9337, -53.8575}, points = {{-50, -8}, {42, 54}}, color = {0, 85, 255}, thickness = 1), Text(origin = {-80, -82}, extent = {{-20, 20}, {20, -20}}, textString = "C")}));
  end Mixer3PH;

  model Mixer4PH
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    replaceable FreeFluids.Interfaces.FluidPortA PortA(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, 71}, {-90, 91}}, rotation = 0), iconTransformation(extent = {{-110, 71}, {-90, 91}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortA PortB(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, 17}, {-90, 37}}, rotation = 0), iconTransformation(extent = {{-110, 17}, {-90, 37}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortA PortC(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, -39}, {-90, -19}}, rotation = 0), iconTransformation(extent = {{-110, -35}, {-90, -15}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortA PortD(redeclare package Medium = Medium, G(min = 0, start = 1)) annotation(
      Placement(visible = true, transformation(extent = {{-110, -91}, {-90, -71}}, rotation = 0), iconTransformation(extent = {{-110, -89}, {-90, -69}}, rotation = 0)));
    replaceable FreeFluids.Interfaces.FluidPortB PortE(redeclare package Medium = Medium, G(max = 0, start = -4), H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default)))) annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
  equation
    PortE.Elevation = PortA.Elevation;
    PortE.X = PortA.X;
    PortB.P = PortA.P;
    PortC.P = PortA.P;
    PortD.P = PortA.P;
    PortE.P = PortA.P;
    0 = PortA.G + PortB.G + PortC.G + PortD.G + PortE.G;
    0 = PortA.G * PortA.H + PortB.G * PortB.H + PortC.G * PortC.H + PortD.G * PortD.H + PortE.G * PortE.H;
    annotation(
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Text(origin = {38, 14}, extent = {{-150, 100}, {150, 30}}, textString = "%name"), Line(origin = {-46, -20}, points = {{-44, 46}, {46, 20}}, color = {0, 85, 255}, thickness = 1), Line(origin = {45, 0}, points = {{-45, 0}, {45, 0}, {45, 0}}, color = {0, 85, 255}, thickness = 1), Line(origin = {-45.258, 41.6388}, points = {{-47, 33}, {45, -41}}, color = {0, 85, 255}, thickness = 1), Polygon(origin = {39.62, 0}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, points = {{-19.6213, 10}, {20.3787, 0}, {-19.6213, -10}, {-19.6213, 10}}), Text(origin = {-78, 42}, extent = {{-20, 20}, {20, -20}}, textString = "B"), Text(origin = {-78, 88}, extent = {{-20, 20}, {20, -20}}, textString = "A"), Text(origin = {84, -24}, extent = {{-20, 20}, {20, -20}}, textString = "E"), Line(origin = {-42.1917, -54.1155}, points = {{-48, 30}, {42, 54}}, color = {0, 85, 255}, thickness = 1), Text(origin = {-76, -34}, extent = {{-20, 20}, {20, -20}}, textString = "C"), Text(origin = {-76, -84}, extent = {{-20, 20}, {20, -20}}, textString = "D"), Line(origin = {-15.4227, -98.2924}, points = {{-76, 24}, {16, 98}}, color = {0, 85, 255}, thickness = 1)}));
  end Mixer4PH;

  model Junction "Bidirectional mixer"
    replaceable FreeFluids.Interfaces.FluidPortA PortA annotation(
      Placement(transformation(extent = {{-110, -7}, {-90, 13}})));
    replaceable FreeFluids.Interfaces.FluidPortB PortB annotation(
      Placement(transformation(extent = {{90, 35}, {110, 55}})));
    replaceable FreeFluids.Interfaces.FluidPortB PortC annotation(
      Placement(transformation(extent = {{90, -55}, {110, -35}})));
  protected
    SI.MassFlowRate G;
    SI.Energy H;
  equation
    PortA.P = PortB.P;
    PortA.P = PortC.P;
    PortA.Elevation = PortB.Elevation;
    PortA.Elevation = PortC.Elevation;
    PortA.X = PortB.X;
    PortA.X = PortC.X;
    0 = PortA.G + PortB.G + PortC.G;
    if PortA.G > 0 then
      if PortB.G > 0 then
        PortC.H = -(PortA.G * PortA.H + PortB.G * PortB.H) / PortC.G;
      elseif PortC.G > 0 then
        PortB.H = -(PortA.G * PortA.H + PortC.G * PortC.H) / PortB.G;
      else
        PortB.H = PortA.H;
        PortC.H = PortA.H;
      end if;
    elseif PortA.G < 0 then
      if PortB.G > 0 then
        if PortC.G > 0 then
          PortA.H = -(PortB.G * PortB.H + PortC.G * PortC.H) / PortA.G;
        else
          PortA.H = PortB.H;
          PortC.H = PortB.H;
        end if;
      else
        PortA.H = PortC.H;
        PortB.H = PortC.H;
      end if;
    end if;
    annotation(
      defaultComponentName = "junct",
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Line(points = {{-100, 0}, {0, 0}}, thickness = 0.5), Line(points = {{0, 0}, {100, 50}}, thickness = 0.5), Line(points = {{0, 0}, {100, -50}}, thickness = 0.5), Text(extent = {{-173, 104}, {175, 62}}, lineColor = {0, 0, 0}, textString = "%name")}));
  end Junction;

  model AbruptAdaptor
    extends Interfaces.TwoFluidPorts(final useElevDifference = true, final elevDifference = 0.0, final calcEnthalpyDifference = true);
    parameter SI.Density rhoStart(displayUnit = "kg/m3") = 5.0 "start value for RhoA and RhoB calculation" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean isCircular = true "the pipe is circular in shape" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Distance diA = Modelica.Constants.inf "internal diameter of PortA if path is circular" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Distance diB = 0.0 "internal diameter of PortB if path is circular" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Area sectionA = 0 "pipe section of PortA if not circular" annotation(
      Dialog(tab = "Physical data"));
    parameter SI.Area sectionB = 0 "pipe section of PortB if not circular" annotation(
      Dialog(tab = "Physical data"));
    SI.Area SectionA;
    SI.Area SectionB;
    Medium.ThermodynamicState StateA, StateB;
    SI.Density RhoA(displayUnit = "kg/m3", start = rhoStart);
    SI.Density RhoB(displayUnit = "kg/m3", start = rhoStart);
    SI.Velocity Va(start = 1) "velocity at PortA";
    SI.Velocity Vb(start = -1) "velocity at PortB";
  algorithm
    if isCircular == true then
      SectionA := 0.25 * pi * diA * diA;
      SectionB := 0.25 * pi * diB * diB;
    else
      SectionA := sectionA;
      SectionB := sectionB;
    end if;
  equation
    StateA = Medium.setState_phX(PortA.P, PortA.H);
    StateB = Medium.setState_phX(PortB.P, PortB.H);
    RhoA = abs(Medium.density(StateA));
    RhoB = abs(Medium.density(StateB));
    Va = PortA.G / RhoA / SectionA;
    Vb = -PortA.G / RhoB / SectionB;
    PortB.P = homotopy(PortA.P + PortA.G * (Va + Vb) / min(SectionA, SectionB), PortA.P) "Mommentum conservation";
    if calcEnthalpyDifference == true then
      0 = PortB.H - PortA.H + 0.5 * (abs(Vb) ^ 2 - abs(Va) ^ 2) "energy conservation";
    end if;
    annotation(
      defaultComponentName = "Red",
      Icon(coordinateSystem(initialScale = 0.05), graphics = {Text(origin = {0, 8}, extent = {{-150, 100}, {150, 30}}, textString = "%name"), Line(points = {{80, 20}, {0, 20}, {0, 40}, {-80, 40}, {-80, -40}, {0, -40}, {0, -20}, {80, -20}, {80, 20}}, color = {7, 119, 255}, thickness = 3)}),
      Documentation(info = "<html>
      <body>
      <p>The pressure assigned to a flow port is that measured at the pipe end, taking into account the velocity of the flow, as this is the preassure needed for the calculation of physical properties like density. When making connections, we assume equal pressure at all connecting ports, without taking into accound that, if the streams have different velocities, there will be a change in pressure. This makes the balance inexact. If we want to improve the calculation, taking into account the velocity impact, it is necessary to reference pressure and enthalpy to the same velocity. This is the function of the abrupt adapter. It is a frictionless adapter.</p>
      <p>If we want to use the adapter, it is necessary to put one adapter at each end of the pipe, except for the connection of two only pipes of the same diameter. In this case the equal velocity is granted at the connection point.</p>
      </body>
      </html>"));
  end AbruptAdaptor;

  model PdiffSource "Introduces a pressure loss between PortA and PortB. Is negative if PortB.P<PortA.P. If applicable, ports are isenthalpic."
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Boolean useFixedDiffP = true "if true fixed dP will be used. Otherwise a  function of mass flow rate" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.PressureDifference dP(displayUnit = "bar") "pressure loss, constant or at reference flow rate. Negative if PortB.P < PortA.P" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.MassFlowRate refG(displayUnit = "kg/h", start = 1.0) "reference mass flow rate. If useFixedDiffP is false" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean isLaminarFlow = false "flow regime to apply. Only if useFixedDiffP is false" annotation(
      Dialog(tab = "Flow"));
  equation
    if useFixedDiffP == true then
      Pdiff = dP;
    else
      if isLaminarFlow == true then
        Pdiff = dP * PortA.G / refG;
      else
        Pdiff = dP * PortA.G ^ 1.77 / refG ^ 1.77 "approximation with correction for the Reynolds in turbulent flow";
      end if;
    end if;
    if calcEnthalpyDifference == true then
      PortA.H = PortB.H "this can lead to lower temperature in PortB than in PortA. alternative is to add the increase in P*V";
    end if;
    annotation(
      defaultComponentName = "Pdiff",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {170, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(extent = {{-68, 54}, {78, -48}}, textString = ">dP<")}));
  end PdiffSource;

  model PipePhysical "Physical description of a pipe of flexible shape, without fluid characteristics."
    parameter Boolean useTubeLength = true "use the supplied length" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance lTube = 0 "supplied individual tube length" annotation(
      Dialog(tab = "Physical data"));
    parameter Boolean fixNumTubes = true "if true fixes the number of tubes to numTubes" annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numTubes = 1 "number of parallel identical paths" annotation(
      Dialog(tab = "Physical data"));
    parameter Boolean isCircular = true "the pipe is circular in shape" annotation(
      Dialog(tab = "Physical data"));
    parameter Boolean useDiameter = true "the supplied diameter will be used for calculation of section and perimeter" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance di(displayUnit = "mm") = 0.0 "supplied internal diameter if path is circular" annotation(
      Dialog(tab = "Physical data"));
    parameter Boolean useSectionAndPerimeter = false "the supplied section and perimeter will be used" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Area section = 0 "supplied pipe section if not circular" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance perimeter = 0 "supplied perimeter if not circular" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance thickness(displayUnit = "mm") = 1e-3 "pipe thickness. Vessel wall thickness for halfcoils" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance roughness(displayUnit = "mm") = 1.5e-005 "pipe roughness. SS:1.5e-5, Steel new:4.6e-5, Steel old:2.0e-4, Concrete:1.5e-3" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance thicknessInsul(displayUnit = "mm") = 0 "insulation thickness" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Density rhoWall(displayUnit = "kg/m3") = 8000 annotation(
      Dialog(tab = "Physical data"));
    Modelica.SIunits.Diameter Di "pipe internal diameter only if circular, otherwise four times Section/Perimeter";
    Modelica.SIunits.Area PathSection "Internal section of the tube. But not the flow section if double pipe";
    Modelica.SIunits.Distance PathPerimeter "Internal perimeter of the tube. But not the flow perimeter if double pipe";
    Modelica.SIunits.Distance Ltotal "total pipe length, as necessary for obtainning its total weight";
    Modelica.SIunits.Distance Ltube(start = 1.0) "pipe path length, as necessary for obtainning friction loss";
    Integer NumTubes "number of parallel identical paths";
    Modelica.SIunits.Area SiTube "each tube internal surface";
    Modelica.SIunits.Area SiTotal "total pipe internal surface";
    Modelica.SIunits.Volume ViTube "each tube internal volume";
    Modelica.SIunits.Volume ViTotal "total pipe internal volume";
    Modelica.SIunits.Diameter Do "pipe external diameter";
    Modelica.SIunits.Diameter Dinsul "pipe diameter with insulation";
    Modelica.SIunits.Area SoTube "each tube external surface";
    Modelica.SIunits.Area SoTotal "total pipe external surface";
    Modelica.SIunits.Volume VoTube "each tube external volume";
    Modelica.SIunits.Volume VoTotal "total pipe external volume";
    Modelica.SIunits.Area SinsulTube "each tube insulation surface";
    Modelica.SIunits.Area SinsulTotal "total pipe insulation surface";
    Modelica.SIunits.Mass MassPipeTotal "total mass of the pipe";
  equation
    if useTubeLength == true then
      Ltube = lTube;
    end if;
    if fixNumTubes == true then
      NumTubes = numTubes;
    end if;
    if isCircular == true then
      PathPerimeter = pi * Di;
      if useDiameter == true then
        Di = di;
      end if;
    elseif useSectionAndPerimeter == true then
      PathSection = section;
      PathPerimeter = perimeter;
    end if;
    MassPipeTotal = (VoTotal - ViTotal) * rhoWall;
    Di = 4 * PathSection / PathPerimeter "Di is four times hydraulic radius, if it is not a double pipe";
    Ltube = Ltotal / NumTubes;
    SiTube = PathPerimeter * Ltube;
    SiTotal = SiTube * NumTubes;
    ViTube = PathSection * Ltube;
    ViTotal = ViTube * NumTubes;
    if isCircular == true then
      Do = Di + 2 * thickness;
      SoTube = pi * Do * Ltube;
      VoTube = pi * Do ^ 2 / 4 * Ltube;
      Dinsul = Do + 2 * thicknessInsul;
      SinsulTube = pi * Dinsul * Ltube;
    else
//If pipe is not circular we treat it as rectangular
      Do = 4 * (PathSection + PathPerimeter * thickness + 4 * thickness * thickness) / (PathPerimeter + 8 * thickness);
      SoTube = (PathPerimeter + 8 * thickness) * Ltube;
      VoTube = (PathSection + PathPerimeter * thickness + 4 * thickness * thickness) * Ltube;
      Dinsul = 4 * (PathSection + PathPerimeter * (thickness + thicknessInsul) + 4 * (thickness + thicknessInsul) ^ 2) / (PathPerimeter + 8 * (thickness + thicknessInsul));
      SinsulTube = (SoTube + 8 * thicknessInsul) * Ltube;
    end if;
    SoTotal = SoTube * NumTubes;
    VoTotal = VoTube * NumTubes;
    SinsulTotal = SinsulTube * NumTubes;
    annotation(
      defaultComponentName = "Pipe",
      Documentation(info = "<html>
  <body>
  <p>The model contains the physical definition of the pipe. Depending on the parameter values, between 0 and 3 equations are missing, related to the sectional shape and length of the pipe.</p>
  <p>If you make the parameter isCircular=true, the Perimeter will be calculated from Di and, as there is inside the model an equation relating Di, Perimeter, and Section, just 1 equation related with the sectional shape will be missing. If you make the parameter useDiameter=true, it will be taken from the parameter di, and the sectional description will be completed. Otherwise Di or Section must be supplied.</p>
  <p>As an alternative to isCircular=true, you can  make the parameter useSectionAndPerimeter=true. In this case Section and Perimeter will be taken from their values in parameters. If none of the previous alternatives has been used, it is necessary to supply equations to solve two of the following variables: Di, Section, Perimeter.</p>
  <p>Regarding the pipe length, is you make the parameter useTubeLength=true. The tube length will be taken from the parameter lTube. Otherwise an equation to solve for the pipe length must be supplied, </p>
  <p>The pipe allows for the use of multiple parallel paths, specified by the parameter numTubes.</p>
  </body>
  </html>"));
  end PipePhysical;

  //***PIPES BASE MODELS***
  //***********************

  partial model Pipe "Pipe with two ports fluid ports and an internal fluid"
    extends PipePhysical;
    extends Interfaces.TwoFluidPorts(calcEnthalpyDifference = true, useElevDifference = true, elevDifference = 0.0, passComposition = true);
    parameter Integer numActiveTubes(start = numTubes, max = numTubes) "number of parallel identical paths used with flow. Default=numTubes" annotation(
      Dialog(tab = "Flow"));
    parameter Real pipeComplexity = 0 "Approx. estimation of equivalent length,0=not used,0.25=supply lines,0.5=long runs, 1=normal,2=valves,4=complex valves" annotation(
      Dialog(tab = "Flow"));
    parameter Real equivL_Di = 0 "Equivalent length of accesories as multiples of Di" annotation(
      Dialog(tab = "Flow"));
    parameter Real numVelocityHeads = 0 "number of velocity heads for pressure loss" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.SIunits.Area kv = 0 "This allows for special head loses to be incorporated inside the pipe calculation via Kv" annotation(
      Dialog(tab = "Flow"));
    parameter Fraction aperture = 1 "Fraction of Kv used" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean fullBore = true "true if the full pipe section is used for flow" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean isCompressibleFlow = false "indicates if kinetic energy and momentum are taken into account or not" annotation(
      Dialog(tab = "Flow"));
    Integer NumActiveTubes "number of parallel identical paths used with flow.";
    Fraction Slope(start = 0.5) "pipe slope as abs(PortB.Elevation-PortA.Elevation)/Ltube";
    Modelica.SIunits.Area SiActive(start = 1.0) "active pipe internal surface";
    Modelica.SIunits.Volume ViActive(start = 1.0) "active pipe internal volume";
    Modelica.SIunits.Area SoActive(start = 1.0) "active pipe external surface";
    Modelica.SIunits.Volume VoActive(start = 1.0) "active pipe external volume";
    Modelica.SIunits.Area SinsulActive(start = 1.0) "active pipe external surface";
    Modelica.SIunits.Distance LeTube "equivalent length of the pipe for pressure loss";
    Modelica.SIunits.Area PathSectionActive "section of each tube used for flow calculation";
    Modelica.SIunits.Distance PathPerimeterActive "perimeter of each tube used for flow calculation";
    Modelica.SIunits.Diameter Dh "four time the hydraulic radius";
    Modelica.SIunits.Mass MassProductTotal "Mass of product inside all the tubes";
    Modelica.SIunits.Pressure PLossFriction(start = 0.01, displayUnit = "bar") "friction head loss";
  equation
    if fixNumTubes == true then
      NumActiveTubes = numActiveTubes;
    end if;
    Slope = abs(PortB.Elevation - PortA.Elevation) / Ltube;
    SiActive = SiTube * NumActiveTubes;
    ViActive = ViTube * NumActiveTubes;
    SoActive = SoTube * NumActiveTubes;
    VoActive = VoTube * NumActiveTubes;
    SinsulActive = SinsulTube * NumActiveTubes;
    if pipeComplexity > 0 then
      LeTube = Ltube * (1 + (0.347 * (Di / 0.0254) ^ 0.5 + 0.216) * pipeComplexity);
    else
      LeTube = Ltube + equivL_Di * Di;
    end if;
    Dh = 4 * PathSectionActive / PathPerimeterActive "This is different from Di, and used for double pipe";
    if fullBore == true then
      PathSectionActive = PathSection "In a double pipe the secction is reduced";
      PathPerimeterActive = PathPerimeter "In a double pipe the perimeter is enlarged";
    end if;
  end Pipe;

  partial model PipeFlowBase "Base model for a single medium pipe."
    extends Pipe;
    parameter FreeFluids.Types.ThermalType thermalType = FreeFluids.Types.ThermalType.detailed "Alternatives: detailed, isenthalpic, adiabatic, isothermal, fixed W, fixed deltaT. If detailed, W calculation must be supplied" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.HeatFlowRate fixedW = 0 "Heat exchange if thermaltype = fixedPower. Positive if heat enters the tube" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.HeatFlowRate fixedDeltaT = 0 "fixed T diff. between ports if thermalType = fixedDeltaT. Positive if Tb > Ta" annotation(
      Dialog(tab = "Heat transfer"));
    Modelica.SIunits.Velocity Va(start = 1);
    Modelica.SIunits.Velocity Vb(start = -1);
    Modelica.SIunits.Density RhoA(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "density at Port A";
    Modelica.SIunits.Density RhoB(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "density at PortB";
    Medium.ThermodynamicState StateA "state at PortA";
    Medium.ThermodynamicState StateB "state at PortB";
    Modelica.SIunits.VolumeFlowRate Qa(displayUnit = "m3/h") "volume flow rate at PortA";
    Modelica.SIunits.VolumeFlowRate Qb(displayUnit = "m3/h") "volume flow rate at PortB";
    Modelica.SIunits.Temperature Ta(displayUnit = "degC", start = Medium.T_default);
    Modelica.SIunits.Temperature Tb(displayUnit = "degC", start = Medium.T_default);
    Modelica.SIunits.Power W "heat transfer between fluid and wall. Positive if fluid inputs heat";
  equation
    StateA = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    StateB = Medium.setState_phX(PortB.P, PortB.H, PortB.X);
    RhoA = abs(Medium.density(StateA));
    RhoB = abs(Medium.density(StateB));
    Va = PortA.G / RhoA / PathSectionActive / NumActiveTubes;
    Vb = PortB.G / RhoB / PathSectionActive / NumActiveTubes;
    Qa = PortA.G / RhoA;
    Qb = PortA.G / RhoB;
    Ta = Medium.temperature(StateA);
    Tb = Medium.temperature(StateB);
    if isCompressibleFlow == true then
      Pdiff = (-sign(PortA.G) * PLossFriction) + (PortA.Elevation - PortB.Elevation) * g_n * (RhoA + RhoB) / 2 + PortA.G * (Va + Vb) / PathSectionActive / NumActiveTubes "Momentum conservation";
      W / PortA.G = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n + 0.5 * (abs(Vb) ^ 2 - abs(Va) ^ 2) "energy conservation";
    else
      Pdiff = (-sign(PortA.G) * PLossFriction) + (PortA.Elevation - PortB.Elevation) * g_n * (RhoA + RhoB) / 2 "momentum change is not taken into account";
      W / PortA.G = PortB.H - PortA.H + (PortB.Elevation - PortA.Elevation) * g_n "kinetic energy is not taken into account";
    end if;
    if calcEnthalpyDifference == true then
      if thermalType == ThermalType.isenthalpic then
        PortA.H = PortB.H;
      elseif thermalType == ThermalType.isothermal then
        Ta = Tb;
      elseif thermalType == ThermalType.adiabatic then
        W = 0;
      elseif thermalType == ThermalType.fixedPower then
        W = fixedW;
      elseif thermalType == ThermalType.fixedDeltaT then
        Tb = Ta + fixedDeltaT;
      end if;
    end if;
    assert(isCompressibleFlow == false or PLossFriction < 0.4 * max(PortA.P, PortB.P), "Too high pressure drop for the pipe model in gas application", AssertionLevel.warning);
//assert(isCompressibleFlow == false or F * Ltube / Dh > 2.0, "Too short pipe for the pipe model in gas application", AssertionLevel.warning);
    annotation(
      defaultComponentName = "Pipe",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 63, 191}, fillColor = {0, 170, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-90, 20}, {90, -20}}), Text(origin = {0, 10}, lineColor = {0, 0, 255}, extent = {{-150, 100}, {146, 42}}, textString = "%name")}),
      Documentation(info = "<html>
  <body>
  <p>The solution for the following variables are missing: MassProductTotal, StateAvg, Re, F, PLossFriction. Plus three more for flow calculation, normally: Ltube, PortA.P, PortB.P. Two of them are supplied by the connectors, so normally only one of them must be specified internally. </p>
  </body>
  </html>"));
  end PipeFlowBase;

  model DoublePipeFlowBase
    parameter Boolean useTubeLengthC = true "use the supplied length" annotation(
      Dialog(tab = "Common data", group = "Physical data"));
    parameter Modelica.SIunits.Distance lTubeC = 0 "supplied individual tube length" annotation(
      Dialog(tab = "Common data", group = "Physical data"));
    parameter Integer numTubesC = 1 "number of parallel identical paths" annotation(
      Dialog(tab = "Common data", group = "Physical data"));
    parameter Integer numActiveTubesC(max = numTubesC) = 1 "number of parallel identical paths used with flow" annotation(
      Dialog(tab = "Common data", group = "Flow"));
    parameter Real pipeComplexityC = 0 "Approx. estimation of equivalent length,0=not used,0.25=supply lines,0.5=long runs, 1=normal,2=valves,4=complex valves" annotation(
      Dialog(tab = "Common data", group = "Flow"));
    parameter Real equivL_DiC = 0 "Equivalent length of accesories as multiples of Di" annotation(
      Dialog(tab = "Common data", group = "Flow"));
    parameter Boolean isCircularE = true "the pipe is circular in shape" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Boolean useDiameterE = true "the supplied diameter will be used for calculation of section and perimeter" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance diE(displayUnit = "mm") = 0.0 "supplied internal diameter if path is circular" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Boolean useSectionAndPerimeterE = false "the supplied section and perimeter will be used" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Area sectionE = 0 "supplied pipe section if not circular" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance perimeterE = 0 "supplied wetted perimeter if not circular" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance roughnessE(displayUnit = "mm") = 1.5e-005 "pipe roughness. SS:1.5e-5, Steel new:4.6e-5, Steel old:2.0e-4, Concrete:1.5e-3" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance thicknessE(displayUnit = "mm") = 1e-3 "pipe thickness" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance thicknessInsulE(displayUnit = "mm") = 0 "insulation thickness" annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Modelica.SIunits.Density rhoWallE(displayUnit = "kg/m3") = 8000 annotation(
      Dialog(tab = "External pipe", group = "Physical data"));
    parameter Boolean isCompressibleFlowE = false annotation(
      Dialog(tab = "External pipe", group = "Flow"));
    parameter FreeFluids.Types.ThermalType thermalTypeE = FreeFluids.Types.ThermalType.detailed "Alternatives: detailed, isenthalpic, adiabatic, isothermal. If detailed, W must be supplied" annotation(
      Dialog(tab = "External pipe", group = "Heat transfer"));
    parameter Boolean isCircularI = true "the pipe is circular in shape" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Boolean useDiameterI = true "the supplied diameter will be used for calculation of section and perimeter" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance diI(displayUnit = "mm") = 0.0 "supplied internal diameter if path is circular" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Boolean useSectionAndPerimeterI = false "the supplied section and perimeter will be used" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Area sectionI = 0 "supplied pipe section if not circular" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance perimeterI = 0 "supplied wetted perimeter if not circular" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance roughnessI(displayUnit = "mm") = 1.5e-005 "pipe roughness. SS:1.5e-5, Steel new:4.6e-5, Steel old:2.0e-4, Concrete:1.5e-3" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Distance thicknessI(displayUnit = "mm") = 1e-3 "pipe thickness" annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Modelica.SIunits.Density rhoWallI(displayUnit = "kg/m3") = 8000 annotation(
      Dialog(tab = "Internal pipe", group = "Physical data"));
    parameter Boolean isCompressibleFlowI = false annotation(
      Dialog(tab = "Internal pipe", group = "Flow"));
    parameter FreeFluids.Types.ThermalType thermalTypeI = FreeFluids.Types.ThermalType.detailed "Alternatives: detailed, isenthalpic, adiabatic, isothermal. If detailed, W must be supplied" annotation(
      Dialog(tab = "Internal pipe", group = "Heat transfer"));
    replaceable FreeFluids.Pipes.PipeFlowBase Ext(final useTubeLength = false, numTubes = numTubesC, isCircular = isCircularE, useDiameter = useDiameterE, di = diE, useSectionAndPerimeter = useSectionAndPerimeterE, section = sectionE, perimeter = perimeterE, roughness = roughnessE, thickness = thicknessE, thicknessInsul = thicknessInsulE, numActiveTubes = numActiveTubesC, pipeComplexity = pipeComplexityC, equivL_Di = equivL_DiC, final kv = 0, final fullBore = false, isCompressibleFlow = isCompressibleFlowE, thermalType = thermalTypeE) annotation(
      Placement(visible = false, transformation(origin = {0, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    replaceable FreeFluids.Pipes.PipeFlowBase Int(final useTubeLength = false, numTubes = numTubesC, isCircular = isCircularI, useDiameter = useDiameterI, di = diI, useSectionAndPerimeter = useSectionAndPerimeterI, section = sectionI, perimeter = perimeterI, roughness = roughnessI, thickness = thicknessI, final thicknessInsul = 0, numActiveTubes = numActiveTubesC, pipeComplexity = pipeComplexityC, equivL_Di = equivL_DiC, final kv = 0, final fullBore = true, isCompressibleFlow = isCompressibleFlowI, thermalType = thermalTypeI) annotation(
      Placement(visible = false, transformation(origin = {-2, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortA PortAi annotation(
      Placement(visible = true, transformation(origin = {-80, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB PortBi annotation(
      Placement(visible = true, transformation(origin = {78, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortA PortAe annotation(
      Placement(visible = true, transformation(origin = {80, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB PortBe annotation(
      Placement(visible = true, transformation(origin = {-80, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    if useTubeLengthC == true then
      Ext.Ltube = lTubeC;
    end if;
    Int.Ltube = Ext.Ltube;
    Ext.PathSectionActive = if Int.isCircular == true then Ext.PathSection - pi * Int.Do ^ 2 / 4 else Ext.PathSection - Int.PathSection - Int.PathPerimeter * Int.thickness + 4 * Int.thickness * Int.thickness "In a double pipe the secction is reduced";
    Ext.PathPerimeterActive = if Int.isCircular == true then Ext.PathPerimeter + pi * Int.Do else Ext.PathPerimeter + Int.PathPerimeter + 8 * Int.Do "In a double pipe the perimeter is enlarged";
    connect(PortBe, Ext.PortB) annotation(
      Line(points = {{-80, 20}, {-10, 20}, {-10, 20}, {-10, 20}}));
    connect(PortAe, Ext.PortA) annotation(
      Line(points = {{80, 20}, {10, 20}, {10, 20}, {10, 20}}));
    connect(PortBi, Int.PortB) annotation(
      Line(points = {{78, -20}, {8, -20}, {8, -20}, {8, -20}}));
    connect(PortAi, Int.PortA) annotation(
      Line(points = {{-80, -20}, {-12, -20}, {-12, -20}, {-12, -20}}));
    annotation(
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(origin = {-20, 14}, lineColor = {0, 85, 255}, extent = {{-60, 6}, {100, -34}}), Rectangle(origin = {0, 35}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 5}, {80, -15}}), Rectangle(origin = {0, -35}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 15}, {80, -5}}), Text(origin = {0, 10}, lineColor = {0, 0, 255}, extent = {{-150, 100}, {146, 42}}, textString = "%name")}));
  end DoublePipeFlowBase;

  model PipeFlow1Ph "Single pipe(not double) for single phase, full bore flow. With no detailed heat transfer model"
    extends PipeFlowBase(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, isCompressibleFlow = false, thermalType = FreeFluids.Types.ThermalType.adiabatic);
    Medium.ThermodynamicState StateAvg "average state for physical properties calculation";
    Modelica.SIunits.Density Rho(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "average density used in calculations";
    Medium.DynamicViscosity Mu(min = 1e-6, start = 1e-3, max = 1e6) "average dynamic viscosity used in calculations";
    Modelica.SIunits.ReynoldsNumber Re(min = 0.01, start = 20000) "Reynolds number at average conditions";
    Real F(start = 0.01) "Darcy's friction factor at average conditions";
    Modelica.SIunits.VolumeFlowRate Q(displayUnit = "m3/h") "volume flow rate in each active tube at average conditions";
    Modelica.SIunits.Velocity V(start = 1) "velocity at average conditions. Normally between 0.9 and 3.0 m/s for liquids";
  equation
    StateAvg = Medium.setState_phX((PortA.P + PortB.P) / 2, (PortA.H + PortB.H) / 2, PortA.X);
    Rho = abs(Medium.density(StateAvg));
    Mu = Medium.dynamicViscosity(StateAvg);
    Q = PortA.G / Rho;
    V = abs(Q / (PathSectionActive * NumActiveTubes));
    MassProductTotal = ViTotal * Rho;
    if noEvent(PortA.G < (-1e-7) or PortA.G > 1e-7) then
      Re = Dh * abs(PortA.G) / (PathSectionActive * NumActiveTubes * Mu) + 0.01 "in order to avoid division by 0 calculating F";
      F = 8 * ((8 / Re) ^ 12 + ((37530 / Re) ^ 16 + (-2.457 * log((7 / Re) ^ 0.9 + 0.27 * roughness / Dh)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
      if kv > 0 then
        PLossFriction = homotopy(0.5 * V ^ 2 * Rho * (F * LeTube / Dh + numVelocityHeads) + 1296000000.0 * (PortA.G / (kv * aperture)) ^ 2 / Rho, 0.5 * V ^ 2 * Rho * (F * LeTube / Dh + numVelocityHeads)) "1.296e9=3600^2*100";
      else
        PLossFriction = 0.5 * V ^ 2 * Rho * (F * LeTube / Dh + numVelocityHeads);
      end if;
    else
      Re = 0.01;
      F = 6400;
      PLossFriction = abs(PortA.G) * 1e3;
    end if;
    annotation(
      Icon(coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html>
  <body>
  <p>The model is mainly for adiabatic (almost isothermal) liquid flow. It can be used also for gas flow, provided that the pressure loss is low (say less than 30% of inlet pressure) and  f*L/D is higher than 2.</p>
  </body>
  </html>"));
  end PipeFlow1Ph;

  model DoublePipeFlow1Ph
    extends DoublePipeFlowBase(redeclare FreeFluids.Pipes.PipeFlow1Ph Ext(final thermalType = FreeFluids.Types.ThermalType.adiabatic), redeclare FreeFluids.Pipes.PipeFlow1Ph Int(final thermalType = FreeFluids.Types.ThermalType.adiabatic));
  end DoublePipeFlow1Ph;

  model PipeFlow2Ph "Single pipe(not double) for two phases flow. With no detailed heat transfer model. Muller-Steinhagen and Heck method"
    extends PipeFlowBase(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, isCompressibleFlow = true, thermalType = ThermalType.adiabatic);
    Modelica.SIunits.Density RhoL(displayUnit = "kg/m3") "density of liquid phase";
    Modelica.SIunits.DynamicViscosity MuL "dynamic viscosity of liquid phase";
    Modelica.SIunits.Velocity Vl(start = 1) "velocity in each active tube if all was liquid";
    Modelica.SIunits.ReynoldsNumber ReL(min = 0.1, start = 20000) "Reynolds number if all was liquid";
    Real Fl(start = 0.01) "Darcy's friction factor if all liquid";
    Modelica.SIunits.Pressure PLossFrictionL(start = 0, displayUnit = "bar") "friction head loss if all was liquid";
    Modelica.SIunits.Density RhoG(displayUnit = "kg/m3") "density of gas phase";
    Modelica.SIunits.DynamicViscosity MuG "dynamic viscosity of gas phase";
    Modelica.SIunits.Velocity Vg(start = 1) "velocity in each active tube if all was gas";
    Modelica.SIunits.ReynoldsNumber ReG(min = 0.1, start = 20000) "Reynolds number if all was gas";
    Real Fg(start = 0.01) "Darcy's friction factor if all gas";
    Modelica.SIunits.Pressure PLossFrictionG(start = 0, displayUnit = "bar") "friction head loss if all was gas";
    Fraction X "gas fraction at average density";
    Medium.ThermodynamicState StateL;
    Medium.ThermodynamicState StateG;
  algorithm
    X := RhoG * (RhoL - (RhoA + RhoB) / 2) / ((RhoA + RhoB) / 2 * (RhoL - RhoG));
    if X < 0 then
      X := 0;
    end if;
  equation
    MassProductTotal = ViTotal * (RhoA + RhoB) / 2;
    StateL = Medium.setBubbleState(Medium.setSat_p(PortA.P));
    RhoL = Medium.density(StateL);
    MuL = Medium.dynamicViscosity(StateL);
    Vl = PortA.G / RhoL / PathSectionActive / NumActiveTubes;
    StateG = Medium.setDewState(Medium.setSat_p(PortA.P));
    RhoG = Medium.density(StateG);
    MuG = Medium.dynamicViscosity(StateG);
    Vg = PortA.G / RhoG / PathSectionActive / NumActiveTubes;
    ReL = Dh * (abs(PortA.G) / PathSectionActive / NumActiveTubes) / MuL;
    Fl = 8 * ((8 / ReL) ^ 12 + ((37530 / ReL) ^ 16 + (-2.457 * log((7 / ReL) ^ 0.9 + 0.27 * roughness / Dh)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
    ReG = Dh * (abs(PortA.G) / PathSectionActive / NumActiveTubes) / MuG;
    Fg = 8 * ((8 / ReL) ^ 12 + ((37530 / ReL) ^ 16 + (-2.457 * log((7 / ReG) ^ 0.9 + 0.27 * roughness / Dh)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
    if kv > 0 then
      PLossFrictionL = 0.5 * Vl ^ 2 * RhoL * Fl / Dh * LeTube + abs(1296000000.0 * PortA.G ^ 2 / RhoL) / (kv * aperture) ^ 2 "1.296e9=3600^2*100";
      PLossFrictionG = 0.5 * Vg ^ 2 * RhoG * Fg / Dh * LeTube + abs(1296000000.0 * PortA.G ^ 2 / RhoG) / (kv * aperture) ^ 2 "1.296e9=3600^2*100";
    else
      PLossFrictionL = 0.5 * Vl ^ 2 * RhoL * Fl / Dh * LeTube;
      PLossFrictionG = 0.5 * Vg ^ 2 * RhoG * Fg / Dh * LeTube;
    end if;
    PLossFriction = (PLossFrictionL + 2 * (PLossFrictionG - PLossFrictionL) * X) * (1 - X) ^ 0.33333 + PLossFrictionG * X ^ 3 "Muller-Steinhagen and Heck";
  end PipeFlow2Ph;

  partial model PipeThermalBase "Base model for pipes with thermal flow"
    extends PipeFlowBase(final pipeComplexity = 0, final equivL_Di = 0, final kv = 0, final aperture = 1, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true, isCompressibleFlow = false, thermalType = ThermalType.detailed, PortB.H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))));
    parameter Modelica.SIunits.ThermalConductivity kWall = 16 "Wall thermal conductivity. typical value for SS=16.7" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.ThermalConductivity kInsul = 0.04 "Insulation thermal conductivity. Typical value for glass fiber=0.04 W/(K·m)" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.ThermalInsulance foulingF = 0.0002 "Tube side fouling factor.Typical: 0.00018 for thermal oil, or treated cooling water" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean useWallsResistance = false "if true, will calculate Twall-Tsurf according to wall conduction, otherwise Twall=Tsurf" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Real emissionCoef = 0.26 "emission coefficient of surface typical values:  Aluminium:0.1, SS:0.15, dirty metal: 0.45, non-metallic:0.94" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fullHTperimeter = true "indicates if all active perimeter is used in heat transfer. It is false for halfpipes" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fullHTlength = true "indicates if all tube length is used in heat transfer. For partial immersion or condensation" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean useThermalConnector = true "if true, the thermal connector will be used for heat transfer. Otherwise a specific calculation must be coded" annotation(
      Dialog(tab = "Heat transfer"));
    Medium.Temperature Twall(start = Medium.T_default) "average wall temperature (tube inside)";
    Medium.Temperature Tsurf(start = Medium.T_default) "average superficial temperature (tube outside)";
    Modelica.SIunits.TemperatureDifference LMTD "Logarithmic mean temperature difference, referenced to Tsurf";
    Modelica.SIunits.Area SactiveHT "internal surface potentially active for heat transfer, due to partial perimeter usage";
    Modelica.SIunits.Area SusedHT(start = 1.0) "internal surface used for heat transfer, due to partial perimeter or partial length usage. Probably to be deprecated";
    Modelica.SIunits.CoefficientOfHeatTransfer H(min = 1, start = 1000) "average heat transfer coefficient of unique or main mechanism";
    FreeFluids.Interfaces.HeatPortB PortH annotation(
      Placement(visible = true, transformation(origin = {1.9984e-15, -48}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.9984e-15, -42}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  equation
    if fullHTperimeter == true then
      SactiveHT = SiActive;
    end if;
    if fullHTlength == true then
      SusedHT = SactiveHT;
    end if;
    if noEvent((Tsurf - Ta) * (Tsurf - Tb) <= 0) then
      LMTD = 0;
    elseif noEvent(not (Ta > Tb or Tb > Ta)) then
      LMTD = Tsurf - Ta "Ta==Tb";
    else
      LMTD = (Tb - Ta) / log((Tsurf - Ta) / (Tsurf - Tb));
    end if;
    if useWallsResistance == true then
      W = SusedHT * 2 / (Di * (log(Do / Di) / kWall + log(Dinsul / Do) / kInsul)) * (Tsurf - Twall) "The definition here of LMTD and the calculation of wall transfer makes them usable only for one transfer mode. Not for example for condensation and subcooling";
    else
      Tsurf = Twall;
    end if;
    PortH.T = Tsurf;
    PortH.S = SusedHT * Dinsul / Di "we pass the external surface to the connector";
    PortH.L = Ltube;
    PortH.D = Dinsul;
    PortH.Slope = Slope;
    if useThermalConnector == true then
      PortH.W = W;
    end if;
    annotation(
      defaultComponentName = "Pipe",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 0, 255}, fillColor = {255, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-90, 20}, {90, -20}})}));
  end PipeThermalBase;

  model PipeForcedConvection "Model for pipes with thermal,single phase, flow"
    extends PipeThermalBase(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, isCompressibleFlow = false, thermalType = ThermalType.detailed, fullHTperimeter = true, fullHTlength = true);
    parameter Boolean usePLossWallCorr = true annotation(
      Dialog(tab = "Flow"));
    parameter Boolean useHTWallCorrFactor = true annotation(
      Dialog(tab = "Heat transfer"));
    Medium.ThermodynamicState StateAvg "average state for physical properties calculation";
    Modelica.SIunits.Density Rho(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "average density used in calculations";
    Medium.DynamicViscosity Mu(min = 1e-6, start = 1e-3, max = 1e6) "average dynamic viscosity used in calculations";
    Modelica.SIunits.ReynoldsNumber Re(min = 0.01, start = 20000) "Reynolds number at average conditions";
    Real F(start = 0.01) "Darcy's friction factor at average conditions";
    Modelica.SIunits.VolumeFlowRate Q(displayUnit = "m3/h") "volume flow rate in each active tube at average conditions";
    Modelica.SIunits.Velocity V(start = 1) "velocity at average conditions. Normally between 0.9 and 3.0 m/s for liquids";
    Modelica.SIunits.HeatCapacity Cp(start = 2000.0);
    Modelica.SIunits.ThermalConductivity K(start = 0.1);
    Modelica.SIunits.PrandtlNumber Pr;
    Real PLossWallCorr(start = 1.0);
    Medium.ThermodynamicState StateW "Thermodynamic state at wall";
    Medium.Temperature T "absolute temperature used for physical properties calc.";
    Medium.DynamicViscosity MuWall(min = 1.0e-6, start = 1.0e-3, max = 1.0e6) "average wall viscosity";
    Real Fsmooth "smooth pipe friction factor";
    Real HTWallCorrFactor(start = 1.0);
  algorithm
    StateAvg := Medium.setState_phX((PortA.P + PortB.P) / 2, (PortA.H + PortB.H) / 2, PortA.X);
    Rho := abs(Medium.density(StateAvg));
    Mu := Medium.dynamicViscosity(StateAvg);
    Q := PortA.G / Rho;
    V := abs(Q / (PathSectionActive * NumActiveTubes));
    MassProductTotal := ViTotal * Rho;
    T := Medium.temperature(StateAvg);
    StateW := Medium.setState_pTX(PortA.P, Twall, PortA.X);
    MuWall := Medium.dynamicViscosity(StateW);
    Cp := Medium.specificHeatCapacityCp(StateAvg);
    K := Medium.thermalConductivity(StateAvg);
    Pr := Cp * Mu / K;
    if usePLossWallCorr == true then
      PLossWallCorr := max((MuWall / Mu) ^ 0.14, 0.25);
    else
      PLossWallCorr := 1;
    end if;
    Re := Dh * abs(PortA.G) / PathSectionActive / NumActiveTubes / Mu + 0.01;
    F := 8 * ((8 / Re) ^ 12 + ((37530 / Re) ^ 16 + (-2.457 * log((7 / Re) ^ 0.9 + 0.27 * roughness / Di)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
//F=0.0791/(Re)^0.25"alternative turbulent";
    PLossFriction := homotopy(0.5 * V ^ 2 * Rho * PLossWallCorr * (F * LeTube / Dh + numVelocityHeads), 0.5 * V ^ 2 * Rho * (F * LeTube / Dh + numVelocityHeads));
    if Rho > 300 and useHTWallCorrFactor == true then
      HTWallCorrFactor := (Mu / MuWall) ^ 0.11 "liquids wall temperature correction factor";
    else
      HTWallCorrFactor := 1;
    end if;
//HTWallCorrFactor := max((Mu / MuWall) ^ 0.11, 0.25) "liquids wall temperature correction factor";
  equation
    if noEvent(Re > 10000) then
      Fsmooth = (0.78173 * log(Re) - 1.5) ^ (-2);
      H = OpenModelica.Internal.realAbs(K / Di * Fsmooth / 8 * Re * Pr / (1 + 12.7 * (Fsmooth / 8) ^ 0.5 * (Pr ^ 0.667 - 1)) * (1 + (Di / Ltube) ^ 0.667) * HTWallCorrFactor) "Pethukov/Gnielinsky equation for smooth tubes";
//H = abs(K / Di * 0.023 * Re ^ 0.8 * Pr ^ 0.333 * HTWallCorrFactor) "Sieder-Tate equation for turbulent flow in smooth pipes";
    elseif noEvent(Re < 2100) then
      Fsmooth = 0;
      H = abs(K / Di * (3.66 ^ 3 + 0.7 ^ 3 + (1.65 * (Re * Pr * Di / Ltube) ^ 0.333 - 0.7) ^ 3 + ((2 / (1 + 22 * Pr)) ^ (1 / 6) * (Re * Pr * Di / Ltube) ^ 0.5) ^ 3) ^ 0.33) "Gnielinsky-Martin correlation: VDI mean";
//H = abs(K / Di * (3.657 + 0.0668 * Re * Pr * Di / Ltube / (1 + 0.04 * (Re * Pr * Di / Ltube) ^ 0.667)) *HTWallCorrFactor) "Hausen correlation";
//H = abs(K / Di * 1.86 * (Re * Pr * Di / Ltube) ^ 0.33 * HTWallCorrFactor) "Sieder-Tate equation for laminar flow";
    else
//interpolation between turbulent and laminar flow
      Fsmooth = (0.78173 * log(Re) - 1.5) ^ (-2);
      H = abs(K / Di * Fsmooth / 8 * Re * Pr / (1 + 12.7 * (Fsmooth / 8) ^ 0.5 * (Pr ^ 0.667 - 1)) * (1 + (Di / Ltube) ^ 0.667) * HTWallCorrFactor) * (Re - 2100) / 7900 + abs(K / Di * (3.66 ^ 3 + 0.7 ^ 3 + (1.65 * (Re * Pr * Di / Ltube) ^ 0.333 - 0.7) ^ 3 + ((2 / (1 + 22 * Pr)) ^ (1 / 6) * (Re * Pr * Di / Ltube) ^ 0.5) ^ 3) ^ 0.33) * HTWallCorrFactor * (10000 - Re) / 7900;
//H = abs(K / Di * Fsmooth / 8 * Re * Pr / (1 + 12.7 * (Fsmooth / 8) ^ 0.5 * (Pr ^ 0.667 - 1)) * (1 + (Di / Ltube) ^ 0.667) * HTWallCorrFactor) * (Re - 2100) / 7900 + abs(K / Di * (3.65 + 0.0668 * Re * Pr * Di / Ltube / (1 + 0.04 * (Re * Pr * Di / Ltube) ^ 0.667)) * (Mu / MuWall) ^ 0.14) * (10000 - Re) / 7900;
//H = abs(K / Di * 0.023 * 10000 ^ 0.8 * Pr ^ 0.333 * PLossWallCorr) * (Re - 2100) / 7900 + abs(K / Di * 1.86 * (2100 * Pr * Di / Ltube) ^ 0.33 * HTWallCorrFactor)* (10000 - Re) / 7900;
    end if;
    W = SactiveHT * 1 / (1 / H + foulingF + Di * (log(Do / Di) / kWall + log(Dinsul / Do) / kInsul) / 2) * LMTD "Twall is only an approximate average, in order to evaluate wall viscosity correction factor";
    annotation(
      defaultComponentName = "Pipe");
  end PipeForcedConvection;

  model PipeFallingFilm "heat exchanged between a flowing fluid and the internal wall. It uses all the active pipe surface"
    extends PipeThermalBase(useElevDifference = true, calcEnthalpyDifference = true, fullBore = false, final isCompressibleFlow = false, fullHTperimeter = true, fullHTlength = true);
    Medium.ThermodynamicState StateAvg "average state for physical properties calculation";
    Modelica.SIunits.Density Rho(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "average density used in calculations";
    Medium.DynamicViscosity Mu(min = 1e-6, start = 1e-3, max = 1e6) "average dynamic viscosity used in calculations";
    Modelica.SIunits.ReynoldsNumber Re(min = 0.01, start = 20000) "Reynolds number at average conditions";
    Modelica.SIunits.HeatCapacity Cp(start = 2000.0);
    Modelica.SIunits.ThermalConductivity K(start = 0.1);
    Modelica.SIunits.PrandtlNumber Pr;
    Modelica.SIunits.DynamicViscosity MuWall(min = 1e-6, start = 1e-3, max = 1e3) "average wall viscosity";
    Modelica.SIunits.ReynoldsNumber ReVDI(min = 0.1, start = 20000) "Uses mass flow rate per unit length = normal Reynolds/4";
    Modelica.SIunits.Distance FilmThickness "falling film thickness";
    Modelica.SIunits.ReynoldsNumber ReCrit "Reynolds for transition zone";
    Real Nup[4] "Nusselt(special) number for laminar, development, transition and turbulent";
    Real HTWallCorrFactor;
    Medium.ThermodynamicState StateW "Thermodynamic state at wall";
  algorithm
    Rho := abs(Medium.density(StateAvg));
    Mu := Medium.dynamicViscosity(StateAvg);
    Re := Dh * abs(PortA.G) / PathSectionActive / NumActiveTubes / Mu + 0.01;
    Cp := Medium.specificHeatCapacityCp(StateAvg);
    K := Medium.thermalConductivity(StateAvg);
    Pr := Cp * Mu / K;
    ReCrit := 2460 * Pr ^ (-0.65);
    PLossFriction := (PortA.Elevation - PortB.Elevation) * Rho * g_n "The friction loss is the same as the gravitational energy loss";
  equation
    StateAvg = Medium.setState_phX((PortA.P + PortB.P) / 2, (PortA.H + PortB.H) / 2, PortA.X);
    PathSectionActive = PathPerimeterActive * FilmThickness "approximate flow section";
    PathPerimeterActive = PathPerimeter;
    MassProductTotal = SiActive * FilmThickness * Rho;
//StateAvg = Medium.setState_phX(PortA.P, (PortA.H + PortB.H) / 2, PortA.X);
    StateW = Medium.setState_pTX(PortA.P, Twall, PortA.X);
    MuWall = Medium.dynamicViscosity(StateW);
    ReVDI = Re / 4;
    FilmThickness = 0.451 * ((Mu / Rho) ^ 2 / 9.81) ^ 0.3333 * ReVDI ^ 0.538 "Karapantsios et al. 1989";
//Re = Di * abs(PortA.G) / (NumActiveTubes * PathSection * Mu);
    HTWallCorrFactor = (Mu / MuWall) ^ 0.25;
    if Slope > 0.25 then
      if FilmThickness < Di / 2 then
        Nup[1] = 1.3 * ReVDI ^ (-0.3333) "laminar Nusselt";
        Nup[2] = 0.912 * (Pr * ((Mu / Rho) ^ 2 / 9.81 * ReVDI) ^ 0.3333 / Ltube) ^ 0.3333 "development Nusselt";
        Nup[3] = 0.0425 * ReVDI ^ 0.2 * Pr ^ 0.344 "transition Nusselt";
        Nup[4] = 0.0136 * ReVDI ^ 0.4 * Pr ^ 0.344 "turbulent Nusselt";
//H=Nu * K * (Rho ^ 2 * 9.81 / Mu ^ 2) ^ 0.3333 "VDI atlas, using maximum Nusselt. Problems with solver";
        if Re < ReCrit then
          H = 0.78 * K * (Cp * Mu / (K * Ltube)) ^ 0.3333 * (Mu ^ 2 / Rho ^ 2 / 9.81) ^ (-2 / 9) * Re ^ (1 / 9) * HTWallCorrFactor "Mueller for laminar";
        elseif Re < 10000 then
          H = 0.032 * K * Re ^ 0.2 * Pr ^ 0.34 * (Mu ^ 2 / Rho ^ 2 / 9.81) ^ (-0.3333) * HTWallCorrFactor "Mueller for transition";
        else
          H = Nup[4] * K * (Rho ^ 2 * 9.81 / Mu ^ 2) ^ 0.3333 * HTWallCorrFactor "VDI atlas for turbulent";
        end if;
      else
        Nup = {0, 0, 0, 0};
        H = 0;
      end if;
    else
      Nup = {0, 0, 0, 0};
      H = 0;
    end if;
    W = SactiveHT * 1 / (1 / H + foulingF + Di * (log(Do / Di) / kWall + log(Dinsul / Do) / kInsul) / 2) * LMTD "Twall is only approximate";
  end PipeFallingFilm;

  partial model PipeCondensingBase "Condensation of pure vapors inside pipes. Twall, Tsurf and H are those of condensation step"
    //Jung:2344; Boyko:4.66; Shah:2344; VDI:P08S; Kutateladze:P08S
    extends PipeFlowBase(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium, final pipeComplexity = 0, final equivL_Di = 0, final kv = 0, final aperture = 1, useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true, isCompressibleFlow = true, thermalType = ThermalType.detailed, PortB.H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))));
    parameter Modelica.SIunits.ThermalConductivity kWall = 16 "Wall thermal conductivity. typical value for SS=16.7" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.ThermalConductivity kInsul = 0.04 "Insulation thermal conductivity. Typical value for glass fiber=0.04 W/(K·m)" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.SIunits.ThermalInsulance foulingF = 0.0002 "Tube side fouling factor.Typical: 0.00018 for thermal oil, or treated cooling water" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean useWallsResistance = false "if true, will calculate Twall-Tsurf according to wall conduction, otherwise Twall=Tsurf" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Real emissionCoef = 0.26 "emission coefficient of surface typical values:  Aluminium:0.1, SS:0.15, dirty metal: 0.45, non-metallic:0.94" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fullHTperimeter = true "indicates if all active perimeter is used in heat transfer. It is false for halfpipes" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fullHTlength = true "indicates if all tube length is used in heat transfer. For partial immersion or condensation" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean useThermalConnector = true "if true, the thermal connector will be used for heat transfer. Otherwise a specific calculation must be coded" annotation(
      Dialog(tab = "Heat transfer"));
    Medium.Temperature Twall(start = Medium.T_default) "average wall temperature (tube inside) along condensation";
    Medium.Temperature Tsurf(start = Medium.T_default) "average superficial temperature (tube outside) along condensation";
    Modelica.SIunits.Area SactiveHT "internal surface potentially active for heat transfer, due to partial perimeter usage";
    Modelica.SIunits.Area SusedHT(start = 1.0) "internal surface used for heat transfer, due to partial perimeter or partial length usage.";
    Modelica.SIunits.CoefficientOfHeatTransfer H(min = 1, start = 1000) "average heat transfer coefficient along condensation";
    FreeFluids.Interfaces.HeatPortB PortH annotation(
      Placement(visible = true, transformation(origin = {1.9984e-15, -48}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.9984e-15, -42}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    parameter CondensationOption condensationOption = FreeFluids.Types.CondensationOption.totalCondensation annotation(
      Dialog(tab = "Heat transfer"));
    Medium.ThermodynamicState StateFilm "average state for condensation film physical properties calculation";
    Modelica.SIunits.Density RhoFilm(displayUnit = "kg/m3", start = if isCompressibleFlow == true then 5.0 else 1000.0) "liquid density at film temperature";
    Medium.Density RhoG(displayUnit = "kg/m3") "saturated gas density";
    Medium.DynamicViscosity MuG "dynamic viscosity of saturated gas";
    Modelica.SIunits.SpecificEnergy Hvc "gas cooling plus condensation heat(J/kg). Inlet can be sligthly superheated";
    Medium.MassFraction Xin "vapor quality at inlet as fraction in mass";
    Medium.MassFraction Xout(start = 0.5) "vapor quality at outlet as fraction in mass";
    Medium.Temperature Tfilm "liquid film average absolute temperature";
    Modelica.SIunits.SpecificHeatCapacity CpFilm(start = 2000.0) "liquid specific heat capacity at film temperature";
    Modelica.SIunits.ThermalConductivity Kfilm(start = 0.1) "liquid thermal conductivity at film temperature";
    Modelica.SIunits.PrandtlNumber PrFilm "liquid Prandt number at film temperature";
    Medium.DynamicViscosity MuFilm(min = 1e-6, start = 1e-3, max = 1e6) "liquid viscosity at film temperature";
    Medium.DynamicViscosity MuWallC(min = 1e-6, start = 1e-3, max = 1e3) "average wall viscosity along condensation";
    Modelica.SIunits.ReynoldsNumber ReC(min = 0.01, start = 20000) "Reynolds number of the liquid phase at the end of condensation";
    Modelica.SIunits.ReynoldsNumber ReVDIc(min = 0.1, start = 200) "At end of condensation. Uses mass flow per unit perimeter length =condensate outlet Reynolds/4";
    Real Xm "mean vapor quality";
    Real NuL, NuT, Nu "Nusselt(special) number laminar, turbulent, combined, for gravity governed";
    Modelica.SIunits.CoefficientOfHeatTransfer Hcg "gravity governed heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer Hcs "gas shear governed heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer Uc(min = 1, start = 1000) "overqall heat transfer coefficient";
    Modelica.SIunits.Power Wc(start = -1) "heat transfer between fluid and wall along condensation. Positive if fluid inputs heat";
    Medium.ThermodynamicState StateGas "saturated gas state";
    Medium.ThermodynamicState StateWallC "Thermodynamic state at wall along condensation";
  equation
    if fullHTperimeter == true then
      SactiveHT = SiActive;
    end if;
    if fullHTlength == true then
      SusedHT = SactiveHT;
    end if;
    if useWallsResistance == true then
      Wc = SusedHT * 2 / (Di * (log(Do / Di) / kWall + log(Dinsul / Do) / kInsul)) * (Tsurf - Twall);
    else
      Tsurf = Twall;
    end if;
    PortH.T = Tsurf;
    PortH.S = SusedHT * Dinsul / Di "we pass the external surface to the connector";
    PortH.L = Ltube;
    PortH.D = Dinsul;
    PortH.Slope = Slope;
    if useThermalConnector == true then
      PortH.W = W;
    end if;
    Xin = Medium.vapourQuality(StateA);
//PLossFriction = (PortA.Elevation - PortB.Elevation) * RhoFilm * g_n "The friction loss is the same as the gravitational energy loss";
    Tfilm = (Ta + Twall) / 2;
    StateFilm = Medium.setState_pTX(PortA.P, Tfilm, PortA.X) "state film is at inlet pressure and film temperature";
    StateWallC = Medium.setState_pTX(PortA.P, Twall, PortA.X);
    RhoFilm = abs(Medium.density(StateFilm));
    MuFilm = Medium.dynamicViscosity(StateFilm);
    CpFilm = Medium.specificHeatCapacityCp(StateFilm);
    Kfilm = Medium.thermalConductivity(StateFilm);
    PrFilm = CpFilm * MuFilm / Kfilm;
    MuWallC = Medium.dynamicViscosity(StateWallC);
    StateGas = Medium.setDewState(Medium.setSat_p(PortA.P));
    RhoG = abs(Medium.density(StateGas));
    MuG = Medium.dynamicViscosity(StateGas);
    Hvc = Medium.specificEnthalpy(StateA) - Medium.specificEnthalpy(StateFilm) "Liquid outlet is at film temperature";
    ReC = Di * abs(PortA.G) * min(1, 1 - Xout) / (NumActiveTubes * PathSection * MuFilm) "Reynolds of the liquid phase at outlet";
    Xm = (Xin + max(0, Xout)) / 2;
    ReVDIc = abs(PortA.G) * min(1, 1 - Xout) / (NumActiveTubes * PathPerimeter * MuFilm) "should be equivalent to ReC / 4";
    if Slope > 0.25 then
      NuL = 0.925 * ((1 - RhoG / RhoFilm) / ReVDIc) ^ 0.333;
      NuT = 0.02 * ReVDIc ^ (7 / 24) * PrFilm ^ 0.333 / (1 + 20.52 * ReVDIc ^ (-3 / 8) * PrFilm ^ (-1 / 6));
      Nu = ((ReVDIc ^ 0.04 * NuL) ^ 1.2 + NuT ^ 1.2) ^ (1 / 1.2) * (MuFilm / MuWallC) ^ 0.25 "VDI mean Nusselt for vertical tubes";
      Hcg = Nu * Kfilm * (RhoFilm ^ 2 * 9.81 / MuFilm ^ 2) ^ 0.333;
//Alternative for laminar HgL=1.47*Kfilm*(RhoFilm*(RhoFilm-RhoG)*9.81/(MuFilm^2*ReC))^(1/3);
    else
      NuL = 0;
      NuT = 0;
      Nu = 0;
      Hcg = 0.761 * Kfilm * (RhoFilm * (RhoFilm - RhoG) * 9.81 * Ltube * NumActiveTubes / (PortA.G * (1 - Xout) * MuFilm)) ^ 0.333 "horizontal tubes. Take into account partial length use of the tube is missing";
    end if;
//For Xout=0, Hs=abs(Kfilm / MuFilm * 0.065 * (RhoFilm * PrFilm * 0.078*(MuG/Di/Vm)^0.25 * Vm ^2/2/RhoG)^0.5), with Vm = abs(0.58 * PortA.G / NumActiveTubes / PathSection)
    Hcs = 0.021 * Kfilm / Di * ReC ^ 0.8 * PrFilm ^ 0.43 * (1 + Xm * (RhoFilm / RhoG - 1)) ^ 0.5 "Boyko equation valid between ReC=1500-15000";
//Dobson=0.023*Kfilm/Di*ReLM^0.8*PrFilm^0.4*(1+2.22/(((1-Xm)/Xm)^0.9*(RhoG/RhoFilm)^0.5*(MuFilm/MuG)^0.1))^0.89;
//Shah=0.023*Kfilm/Di*ReC*0.5^0.8*PrFilm^0.4*(1-Xm)^0.8*(1.8/((1/Xm-1)^0.8*(RhoG/RhoFilm)^0.5)^0.8)"ReC>350";
    if Hcg > Hcs then
      H = Hcg "H is that of condensation";
    else
      H = Hcs "H is that of condensation";
    end if;
    1 / Uc = 1 / H + foulingF + Di * 0.5 * (log(Do / Di) / kWall + log(Dinsul / Do) / kInsul);
    Wc = SusedHT * Uc * (Tsurf - Ta) "main relation between Tsurf and W";
    Wc = -abs(PortA.G) * (Xin - Xout) * Hvc;
//If we force the flow as per having Xout>0. The gas fraction in PortB will be sligthly lower than Xout. Due to condensation for heating the liquid from Tfilm to T boil.
    annotation(
      defaultComponentName = "pipe",
      Documentation(info = "<html>
  <body>
  <p>As the pressure drop has been set to 0, it is not possible the calculation of the massic flow from it. But it is possible its calculation outlet gas fraction.</p>
  
  </body>
  </html>"));
  end PipeCondensingBase;

  model PipeCondensing
    extends PipeCondensingBase(fullHTlength = true, condensationOption = FreeFluids.Types.CondensationOption.totalCondensation);
  equation
    if condensationOption == FreeFluids.Types.CondensationOption.partialCondensation then
      Pdiff = 0 "This makes the friction pressure loss equal to velocity and height loss";
      assert(Xout >= 0.0, "Partial condensation has been selected, but total happened: change model to PipeCondSubcool", AssertionLevel.warning);
    elseif condensationOption == FreeFluids.Types.CondensationOption.totalCondensation then
      Xout = 0;
    else
      assert(false, "you must use total or partial condensation");
    end if;
    MassProductTotal = ViTotal * SusedHT / SactiveHT * RhoFilm * RhoG / (Xm * RhoFilm + (1 - Xm) * RhoG);
    W = Wc "all exchanged heat is due to condensation";
  end PipeCondensing;

  partial model CoilPhysical "Gives pipes additional information regarding the coil shape"
    parameter Modelica.SIunits.Diameter coilDiam = 1.5 "Coil diameter" annotation(
      Dialog(tab = "Physical data"));
    parameter Real num = 10 "total number of spires" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance path = 0.15 "distance between centers of coil spires in meters(pitch)" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance heightInit = 0 "initial heigh of coil over reference (normally tangent) line" annotation(
      Dialog(tab = "Physical data"));
    Modelica.SIunits.Distance CoilHeigth;
    Modelica.SIunits.Distance CoilFinalHeight;
    //duplicate definitions
    Modelica.SIunits.Distance Ltotal "total pipe length, as necessary for obtainning its total weight";
  algorithm
    Ltotal := pi * coilDiam * num;
    CoilHeigth := path * num;
    CoilFinalHeight := heightInit + CoilHeigth;
  end CoilPhysical;

  model CoilForcedConvection "Power exchanged between a flowing fluid and the internal wall. Needs an external heat flow model. It uses all the active coil surface"
    //Normally lacking:PortA elevation,G,P,H and some way to calculate W
    extends CoilPhysical;
    extends PipeForcedConvection(final useTubeLength = false, final lTube = 0, final isCircular = true, final useDiameter = true, final useSectionAndPerimeter = false, final section = 0, final perimeter = 0);
    annotation(
      defaultComponentName = "coil");
  end CoilForcedConvection;

  //***COIL MODELS***

  partial model HalfCoilPhysical
    parameter Modelica.SIunits.Diameter basePipeDi(displayUnit = "mm") = 0 "Internal diameter of the base pipe used" annotation(
      Dialog(tab = "Physical data"));
    //parameter Integer basePipeAngle = 180 "Alternative is 120" annotation(
    //  Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Angle basePipeAngle(displayUnit = "deg") = Modelica.Constants.pi "sector of the full pipe used" annotation(
      Dialog(tab = "Physical data"));
    parameter Real num = 0 "total number of spires" annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Distance path = 0 "distance between centers of coil spires in meters" annotation(
      Dialog(tab = "Physical data"));
    parameter Boolean isBottomJacket = false annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Diameter lowerHalfCoilDiam = 0 annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.SIunits.Diameter largerHalfCoilDiam = 0 annotation(
      Dialog(tab = "Physical data"));
    Modelica.SIunits.Diameter HalfCoilDiam "Half coil diameter";
    Modelica.SIunits.Distance HalfCoilHeigth;
    Modelica.SIunits.Area SauxHT "free surface between active spires";
    //duplicate definitions
    Modelica.SIunits.Area SiActive(start = 1.0) "active pipe internal surface";
    Modelica.SIunits.Area PathSection "Internal section of the tube. But not the flow section if double pipe";
    Modelica.SIunits.Distance PathPerimeter "Internal perimeter of the tube. But not the flow perimeter if double pipe";
    Modelica.SIunits.Distance Ltotal "total pipe length, as necessary for obtainning its total weight";
    Modelica.SIunits.Area SactiveHT "internal surface potentially active for heat transfer, due to partial perimeter usage";
  algorithm
    HalfCoilHeigth := path * num;
/*
  if basePipeAngle == 180 then
    PathSection := 0.3927 * basePipeDi ^ 2;
    PathPerimeter := basePipeDi * 2.5708;
    SactiveHT := SiActive * 0.389;
    SauxHT := SactiveHT * (path - basePipeDi) / basePipeDi;
  elseif basePipeAngle == 120 then
    PathSection := 0.1535 * basePipeDi ^ 2;
    PathPerimeter := basePipeDi * 1.9132;
    SactiveHT := SiActive * 0.4527;
    SauxHT := SactiveHT * (path - 0.866 * basePipeDi) / (0.866 * basePipeDi);
  else
    PathSection := 0;
    PathPerimeter := 0;
    SactiveHT := 0;
    SauxHT := 0;
  end if;*/
    PathSection := basePipeDi ^ 2 * (basePipeAngle / 2 - sin(basePipeAngle / 2) * cos(basePipeAngle / 2)) / 4;
    PathPerimeter := basePipeDi * (basePipeAngle / 2 + sin(basePipeAngle / 2));
    SactiveHT := SiActive * basePipeDi * sin(basePipeAngle / 2) / PathPerimeter;
    SauxHT := SactiveHT * (path - basePipeDi * sin(basePipeAngle / 2)) / basePipeDi * sin(basePipeAngle / 2);
  equation
    Ltotal = pi * HalfCoilDiam * num;
  end HalfCoilPhysical;

  model HalfCoilForcedConvection "same as PipeThermalStd, but using as active heat transfer surface only de tank side"
    extends HalfCoilPhysical;
    extends PipeForcedConvection(final useTubeLength = false, final lTube = 0, final isCircular = false, final useDiameter = false, final fullHTperimeter = false, final useSectionAndPerimeter = false, final di = 0, final section = 0, final perimeter = 0, final thicknessInsul = 0);
    annotation(
      defaultComponentName = "HalfCoil",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Line(origin = {0, -19.72}, points = {{-90, 19.7236}, {-80, 19.7236}, {-60, -20.2764}, {-40, 19.7236}, {-20, -20.2764}, {0, 19.7236}, {20, -20.2764}, {40, 19.7236}, {60, -20.2764}, {80, 19.7236}, {90, 19.7236}, {90, 19.7236}}, color = {255, 85, 255}, thickness = 1)}));
  end HalfCoilForcedConvection;

  model HalfCoilCondensing
    extends HalfCoilPhysical;
    extends PipeCondensing(final useTubeLength = false, final lTube = 0, final isCircular = false, final useDiameter = false, final fullHTperimeter = false, final useSectionAndPerimeter = false, final di = 0, final section = 0, final perimeter = 0, final thicknessInsul = 0, condensationOption = FreeFluids.Types.CondensationOption.totalCondensation);
    annotation(
      defaultComponentName = "HalfCoil",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Line(origin = {0, -19.72}, points = {{-90, 19.7236}, {-80, 19.7236}, {-60, -20.2764}, {-40, 19.7236}, {-20, -20.2764}, {0, 19.7236}, {20, -20.2764}, {40, 19.7236}, {60, -20.2764}, {80, 19.7236}, {90, 19.7236}, {90, 19.7236}}, color = {255, 85, 255}, thickness = 1)}));
  end HalfCoilCondensing;

  model PipeFluidBoundary "Models the external heat exchange against a fluid"
    FreeFluids.Interfaces.HeatPortA PortH annotation(
      Placement(visible = true, transformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable package Medium = Modelica.Media.Air.DryAirNasa constrainedby Modelica.Media.Interfaces.PartialMedium;
    parameter Modelica.SIunits.ThermalInsulance foulingF = 0.0002 "Fouling factor of the surface";
    parameter Modelica.SIunits.Temperature tMedia(displayUnit = "degC") = 298.15 "Absolute temperature of sorrounding media";
    parameter Modelica.SIunits.AbsolutePressure pMedia(displayUnit = "bar") = 1e5;
    parameter Modelica.SIunits.Velocity vMedia = 3.5 "sourranding media velocity";
    parameter Real emissionCoef = 0.26 "emission coefficient typical values:  Aluminium:0.1, SS:0.15, dirty metal: 0.45, non-metallic:0.94";
    Medium.ThermodynamicState State "thermodynamic state to get access to properties of the medium";
    Modelica.SIunits.Density Rho(displayUnit = "kg/m3") "medium density";
    Modelica.SIunits.DynamicViscosity Mu "medium viscosity";
    Modelica.SIunits.SpecificHeatCapacityAtConstantPressure Cp "medium heat capacity";
    Modelica.SIunits.ThermalConductivity Lambda;
    Medium.IsobaricExpansionCoefficient Beta;
    //Modelica.SIunits.Length Lcn "characteristic length for natural convection";
    Modelica.SIunits.PrandtlNumber Pr "Prandtl number";
    Modelica.SIunits.ReynoldsNumber Re "Reynolds number";
    Modelica.SIunits.RayleighNumber Ra "Rayleigh number";
    Modelica.SIunits.NusseltNumber NuN "Nusselt number for natural convection";
    Modelica.SIunits.NusseltNumber NuF "Nusselt number for forced convection";
    Modelica.SIunits.CoefficientOfHeatTransfer Hn(start = 5) "natural convection heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer Hf(start = 100) "forced convection heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer Hc(start = 100.0) "combined convection heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer Hr(start = 10.0) "radiation heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer H(start = 100.0) "total heat transfer coefficient";
    Modelica.SIunits.CoefficientOfHeatTransfer U(min = 1, start = 100) "overall heat transfer coeff. at mean conditions";
  equation
    State = Medium.setState_pTX(pMedia, (PortH.T + tMedia) / 2);
    Rho = Medium.density(State);
    Mu = Medium.dynamicViscosity(State);
    Cp = Medium.specificHeatCapacityCp(State);
    Lambda = Medium.thermalConductivity(State);
    Beta = Medium.isobaricExpansionCoefficient(State);
    Pr = Cp * Mu / Lambda;
//Lcn=(PortH.L*PortH.D/(PortH.L*cos(asin(PortH.Slope))/PortH.D+PortH.D*PortH.Slope/PortH.L))^0.5 "characteristic length for inclined pipes";
    Re = PortH.D * vMedia * Rho / Mu;
    Ra = g_n * Beta * Rho * Rho * Cp * PortH.D * PortH.D * PortH.D * abs(PortH.T - tMedia) / (Mu * Lambda);
//Ra=g_n*Beta*Rho*Rho*Cp*Lcn*Lcn*Lcn*abs(PortH.T-tMedia)/(Mu*Lambda) "Rayleigh number in inclined pipes";
    NuN = (0.6 + 0.387 * Ra ^ 0.166667 / (1 + (0.559 / Pr) ^ 0.5625) ^ 0.296296) ^ 2 "Nusselt number for natural convection in horizontal. Churchill-Chu 1975";
//NuN=(0.54+0.39*(Ra/(1+(0.559/Pr)^0.5625)^1.777778)^0.1685)^2 "Nusselt number for natural convection in inclined pipes";
    NuF = 0.3 + 0.62 * Re ^ 0.5 * Pr ^ 0.33333 * (1 + (Re / 282000) ^ 0.625) ^ 0.8 / (1 + (0.4 / Pr) ^ 0.666667) ^ 0.25;
    Hn = NuN * Lambda / PortH.D;
//Hn=NuN*Lambda/Lcn "natural convection heat transfer coefficient in inclined pipes";
    Hf = NuF * Lambda / PortH.D;
    Hc = (Hn ^ 4 + Hf ^ 4) ^ 0.25;
    if Rho < 300 then
      Hr = abs((PortH.T ^ 4 - tMedia ^ 4) / (PortH.T - tMedia) * emissionCoef * 5.67e-8);
    else
      Hr = 0.0;
    end if;
    H = Hc + Hr;
/*if vAir * PortH.D < 8.55e-3 then
    Hc = 8.1e-3 / PortH.D + 3.14 * (vAir / PortH.D) ^ 0.5;
  else
    Hc = 8.9 * vAir ^ 0.9 / PortH.D ^ 0.1;
  end if;*/
    1 / U = 1 / H + foulingF;
    PortH.W = PortH.S * U * (PortH.T - tMedia);
/*PortH.S = 0;
  PortH.L = 100;
  PortH.D = 0.0909;
  PortH.Slope = 0;
  PortH.T = 129.4 + 273.15;*/
    annotation(
      Icon(graphics = {Line(origin = {37.5412, 33.3067}, points = {{-56, 17}, {44, -11}}), Line(origin = {29.8566, 17.2874}, points = {{-58, 11}, {52, -1}}), Line(origin = {21.7386, 9.49529}, points = {{-54, 0}, {60, 0}, {58, 0}}), Line(origin = {25.2782, -23.999}, points = {{-44, -7}, {56, 21}}), Line(origin = {21.3346, 1.11031}, points = {{-50, -12}, {60, 2}}), Ellipse(origin = {-55, 55}, rotation = 70, extent = {{19, 17}, {-5, -7}}, endAngle = 180), Ellipse(origin = {-63, 21}, rotation = 80, extent = {{19, 17}, {-5, -7}}, endAngle = 180), Ellipse(origin = {-59, -15}, rotation = 100, extent = {{19, 17}, {-5, -7}}, endAngle = 180), Ellipse(origin = {-51, -47}, rotation = 110, extent = {{19, 17}, {-5, -7}}, endAngle = 180), Text(origin = {11, -81}, extent = {{-57, 37}, {37, -23}}, textString = "Fluid")}, coordinateSystem(initialScale = 0.1)));
  end PipeFluidBoundary;

end Pipes;

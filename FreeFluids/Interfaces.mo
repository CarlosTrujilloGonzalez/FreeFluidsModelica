within FreeFluids;

package Interfaces "Interfaces.mo by Carlos Trujillo
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


  connector FluidPort "for the connection of fluid flow, without reference to the mechanical characteristics of the connection"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    Medium.AbsolutePressure P(displayUnit = "bar") "Pressure";
    flow Modelica.Units.SI.MassFlowRate G(displayUnit = "kg/h") "Mass flow";
    replaceable Modelica.Units.SI.Height Elevation(start = 0) "Port Height";
    replaceable Medium.SpecificEnthalpy H "Specific enthalpy";
    replaceable Medium.MassFraction X[Medium.nX] "mass fractions composition";
    annotation(
      Documentation(info = "<html><head></head><body>
  <div class=\"standard\" id=\"magicparlabel-8207\">The FluidPort connector is the base of all the Interfaces package. Its composition is:</div><ul class=\"itemize\" id=\"magicparlabel-8208\"><li class=\"itemize_item\">replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;</li>
<li class=\"itemize_item\">SI.AbsolutePressure P(displayUnit = \"bar\") \"Pressure\";</li>
<li class=\"itemize_item\">flow SI.MassFlowRate G(displayUnit = \"kg/h\") \"Mass flow\";</li>
<li class=\"itemize_item\">replaceable SI.Height Elevation(start = 0) \"Port Height\";</li>
<li class=\"itemize_item\">replaceable SI.SpecificEnthalpy H \"Specific enthalpy\";</li>
<li class=\"itemize_item\">replaceable Medium.MassFraction X[Medium.nX] \"mass fractions composition\";</li>
</ul><div class=\"standard\" id=\"magicparlabel-8214\">An explanation is needed regarding why and how each element is defined.</div><div class=\"standard\" id=\"magicparlabel-8215\">The package Medium is included in the connector mainly because it seems the standard in Modelica.Fluid. It will be useful if the propagation of the Medium package is done, in the future, by the connectors. It is useful also for initialization.</div><div class=\"standard\" id=\"magicparlabel-8216\">P and G are according to the Modelica implementation of a connector: one potential variable and one flow variable. This aids to grant that the number of equations and the number of variables are the same if the relationship between elements is done only by connections.</div><div class=\"standard\" id=\"magicparlabel-8217\">The variables Elevation and H don't conform to the Modelica.Fluid way, that uses stream variables. The Elevation could be avoid, as the only information needed inside the elements with connectors is the elevation difference between ports, not the absolute elevation. But its absence would avoid to enter information, based on elevation data at the ends, directly.</div><div class=\"standard\" id=\"magicparlabel-8218\">The enthalpy is absolutely necessary if we want to pass energy information between connectors. The solution adopted by Modelica 3.0 was the creation of stream variables, but I think that this complicates too much the situation. The alternative is to declare them as causal (input or output). This seems to prevent from allowing reverse flow, but it is not totally true, as can be seen in several examples. They are declared replaceable in order to allow its declaration as causal in derived connectors.</div><p><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->






















</p><div class=\"standard\" id=\"magicparlabel-8219\">From the FluidPort connector, a FluidPortA as normal input, and FluidPortB as normal output are derived. The difference is the appearance, the default value for G (positive for the normal input, and negative for the normal output), and the causality of the extra variables. Initialization is also provided for the H and P of FluidPortB, based in the default values of the Medium.</div>
  
  </body></html>"));
  end FluidPort;

  connector FluidPortA "Is the normal fluid inlet port"
    extends FluidPort(redeclare input SI.Height Elevation(start = 0), redeclare input Medium.SpecificEnthalpy H, redeclare input Medium.MassFraction X[Medium.nX](start = Medium.X_default), G(start = 1.0));
    annotation(
      defaultComponentName = "PortA",
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid)}));
  end FluidPortA;

  connector FluidPortB "Is the normal fluid outlet port"
    extends FluidPort(redeclare output SI.Height Elevation(start = 0), redeclare output SI.SpecificEnthalpy H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))), redeclare output Medium.MassFraction X[Medium.nX](start = Medium.X_default), G(start = -1.0), P(start = Medium.p_default));
    annotation(
      defaultComponentName = "PortB",
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}),
      Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
  end FluidPortB;

  partial model TwoFluidPorts "To be used by equipment with two fluid ports."
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    FreeFluids.Interfaces.FluidPortA PortA(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
    FreeFluids.Interfaces.FluidPortB PortB(redeclare package Medium = Medium, P(start = Medium.p_default), H(start = Medium.specificEnthalpy(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))), X(start = Medium.X_default)) annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
    parameter Boolean useElevDifference = true "introduces, or not, an equation for the ports elevation. Use only in one connector" annotation(
      Dialog(tab = "Ports links"));
    parameter FreeFluids.Types.ElevationOption elevCalcMethod = FreeFluids.Types.ElevationOption.differential "if useElevDifference==true, selects the equation to use: differential, or absolute for PortB" annotation(
      Dialog(tab = "Ports links"));
    parameter Modelica.Units.SI.Length elevDifference = 0.0 "used only if previously selected. Positive if PortB higher than PortA" annotation(
      Dialog(tab = "Ports links"));
    parameter Modelica.Units.SI.Length portBelevation = 0.0 "used only if previously selected." annotation(
      Dialog(tab = "Ports links"));
    parameter Boolean calcEnthalpyDifference = true "if true, will calculate the enthalpy difference between ports. Use only in one connector" annotation(
      Dialog(tab = "Ports links"));
    parameter Boolean passComposition = true "if true, the equation PortA.X=PortB.X will be activated" annotation(
      Dialog(tab = "Ports links"));
    Modelica.Units.SI.SpecificEnthalpy Hdiff "specific enthalpy difference: Port B - Port A";
    Modelica.Units.SI.PressureDifference Pdiff(displayUnit = "bar") "pressure difference: Port B - Port A";
  equation
    0 = PortA.G + PortB.G;
    if useElevDifference == true then
      if elevCalcMethod == FreeFluids.Types.ElevationOption.absolute then
        PortB.Elevation = portBelevation;
      else
        PortB.Elevation = PortA.Elevation + elevDifference;
      end if;
    end if;
    if passComposition == true then
      PortA.X = PortB.X;
    end if;
    Hdiff = PortB.H - PortA.H;
    Pdiff = PortB.P - PortA.P;
    annotation(
      Documentation(info = "<html><head></head><body>
  <div class=\"standard\" id=\"magicparlabel-8270\">This partial model contains a FluidPortA PortA, and a FluidPortB PortB connectors, and some equations linking the values of both connectors variables. The more fundamental equation is:</div><div class=\"standard\" id=\"magicparlabel-8271\">0=PortA.G+PortB.G. This grants that all massic flow entering one port goes out by the other. G is considered positive if the flow enters the port.</div><div class=\"standard\" id=\"magicparlabel-8272\">Some parameters will configure more possible links between the ports variables:</div><ul class=\"itemize\" id=\"magicparlabel-8273\"><li class=\"itemize_item\">parameter Boolean useElevDifference=true will activate the use of an extra equation for ports elevation.</li>
<li class=\"itemize_item\">parameter FreeFluids.Types.ElevationOption elevCalcMethod. If useElevDifference=true, selects the equation to use for the elevation calculation. If the selection is “differential” the equation PortB.elevation=PortA.Elevation+elevDifference is applied. If the selection is “absolute” the equation PortB.Elevation= portBelevation is applied. Being elevDifference and portBelevation parameters with default value of 0.</li>
<li class=\"itemize_item\">parameter Boolean calcEnthalpyDifference=true will control if a calculation of the difference of enthalpy between ports is applied or not. But the calculation to apply remains undefined.</li>
<li class=\"itemize_item\">parameter Boolean passComposition = true will control if the composition of both ports is made equal or not.</li>
</ul><div class=\"standard\" id=\"magicparlabel-8277\">Finally two variables, with the equations to solve them are added: Hdiff = PortB.H-PortA.H, and Pdiff=PortB.P-PortA.P.</div><div class=\"standard\" id=\"magicparlabel-8278\">With this implementation, the model is imbalanced in two ways: For the two pairs of potential/flow variables, there are three equations(PortA.G=0, PortB.G=0, PortA.G=-PortB.G), so an equation linking these variables is missing. Normally it will be the calculation of the pressure drop as function of flow. The second point of imbalance is that there is no equation for the calculation of the output PortB.H. The extending models must provide the missing equations.</div><div class=\"standard\" id=\"magicparlabel-8279\">Although not correct according to the Modelica standard, the model is prepared for the connection of more than one PortB to a PortA, with the idea of making simpler and faster some connections. When doing so we must take into account that, in a connection, point elevation, enthalpy, and composition, must be supplied only by one connector. If not, we will have an over specification for the value. The situation is different for elevation than for enthalpy and composition. In a connection point all the ports must have the same elevation, so the only thing that we have to do is to inhibit the duplicate propagation of the elevation, making useHeightDifference=true only in one of the elements that can supply the elevation value to the connection. For enthalpy and composition, it must be possible to connect elements with different enthalpy or composition output, and this can be done using mixers, that for simplification are coded as no reverse flow. But, if the enthalpy variation inside the elements is only due to elevation changes, we can grant that all connectors have the same enthalpy at the connection point, and solve the problem in the same way that we did with elevation, allowing reverse flow. The same is applicable for composition, if composition is constant.</div><div class=\"standard\" id=\"magicparlabel-8280\">A second point must be considered: In a connection of more than two connectors, regardless the physical size of the connections could be the same, there is normally a change in velocity, so in kinetic energy and enthalpy. This means that the equal enthalpy solution is only an approximate solution valid when the velocity is low (liquids and gases at not high velocity).</div><p><!--?xml version=\"1.0\" encoding=\"UTF-8\"?-->
























</p><div class=\"standard\" id=\"magicparlabel-8281\">As a resume, we will develop elements allowing for reverse flow, based in a fixed enthalpy and composition at each connection point, and elements without reverse flow capability, valid for situations with different enthalpy, or composition, of the connectors.</div>
  
  </body></html>"));
  end TwoFluidPorts;

  model FlowSource "Single or two phases flow source. Mainly necessary to generate pressure, enthalpy, and composition values at a connection"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium;
    replaceable FreeFluids.Interfaces.FluidPortB PortB(redeclare package Medium=Medium) annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
    parameter FreeFluids.Types.SourceOption sourceOption = FreeFluids.Types.SourceOption.useP_T "how to generate port pressure and enthalpy";
    parameter Boolean isElevSource = true "use specified gravitational height";
    parameter Boolean isGsource = false "use specified mass flow";
    parameter Modelica.Units.SI.Height Elevation = 0 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.MassFlowRate G(displayUnit = "kg/h") = 0 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.AbsolutePressure P(displayUnit = "bar") = 1e5 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.Temperature T(displayUnit = "degC") = 298.15 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.SpecificEnthalpy H = if sourceOption == FreeFluids.Types.SourceOption.useP_T then Medium.specificEnthalpy(Medium.setState_pTX(P, T, X)) elseif sourceOption == FreeFluids.Types.SourceOption.useD_T then Medium.specificEnthalpy(Medium.setState_dTX(D, T, X)) elseif sourceOption == FreeFluids.Types.SourceOption.useSatGasP then Medium.dewEnthalpy(Medium.setSat_p(P))  elseif sourceOption == FreeFluids.Types.SourceOption.useSatGasT then Medium.dewEnthalpy(Medium.setSat_T(T)) else Medium.h_default annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.Density D(displayUnit = "kg/m3") = 100 annotation(
      Dialog(tab = "Local values"));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default annotation(
      Dialog(tab = "Local values"));
    /*conditional connectors*/
    parameter Boolean externalG = false annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalT = false annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalP = false annotation(
      Dialog(tab = "External values"));
    Modelica.Blocks.Interfaces.RealInput Gext(unit = "kg/s", displayUnit = "kg/h") if externalG == true annotation(
      Placement(visible = true, transformation(origin = {-60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {-68, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealInput Text(unit = "K", displayUnit = "degC") if externalT == true annotation(
      Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealInput Pext(unit = "Pa", displayUnit = "bar") if externalP == true annotation(
      Placement(visible = true, transformation(origin = {60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {66, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  protected
    Modelica.Blocks.Interfaces.RealInput GextIn "to connect with the conditional input block";
    Modelica.Blocks.Interfaces.RealInput TextIn;
    Modelica.Blocks.Interfaces.RealInput PextIn;
  equation
    connect(Gext, GextIn);
    connect(Text, TextIn);
    connect(Pext, PextIn);
    if externalG == false then
      GextIn = G;
    end if;
    if externalT == false then
      TextIn = T;
    end if;
    if externalP == false then
      PextIn = P;
    end if;
    if isGsource == true then
      PortB.G = -GextIn;
    end if;
    if isElevSource == true then
      PortB.Elevation = Elevation;
    end if;
    if sourceOption == FreeFluids.Types.SourceOption.useP_H then
      PortB.P = PextIn;
      PortB.H = H;
    elseif sourceOption == FreeFluids.Types.SourceOption.useP_T then
      PortB.P = PextIn;
      PortB.H = Medium.specificEnthalpy(Medium.setState_pTX(PextIn, TextIn, X));
    elseif sourceOption == FreeFluids.Types.SourceOption.useD_T then
      PortB.P = Medium.pressure(Medium.setState_dTX(D, TextIn, X));
      PortB.H = Medium.specificEnthalpy(Medium.setState_dTX(D, TextIn, X));
    elseif sourceOption == FreeFluids.Types.SourceOption.useSatLiqT then
      PortB.P = Medium.saturationPressure(TextIn);
      PortB.H = Medium.bubbleEnthalpy(Medium.setSat_T(TextIn));
    elseif sourceOption == FreeFluids.Types.SourceOption.useSatGasT then
      PortB.P = Medium.saturationPressure(TextIn);
      PortB.H = Medium.dewEnthalpy(Medium.setSat_T(TextIn));
    elseif sourceOption == FreeFluids.Types.SourceOption.useSatLiqP then
      PortB.P = PextIn;
      PortB.H = Medium.bubbleEnthalpy(Medium.setSat_p(PextIn));
    elseif sourceOption == FreeFluids.Types.SourceOption.useSatGasP then
      PortB.P = PextIn;
      PortB.H = Medium.dewEnthalpy(Medium.setSat_p(PextIn));
    end if;
    PortB.X = X;
    annotation(
      defaultComponentName = "Source",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {0, 48, 144}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(origin = {-4, 14}, extent = {{-54, 46}, {58, -40}}, textString = "Source"), Text(origin = {0, -40}, extent = {{-54, 46}, {20, -8}}, textString = "2 Ph"), Text(origin = {-100, 110}, extent = {{-20, 24}, {16, -22}}, textString = "G"), Text(origin = {98, 114}, extent = {{-16, 22}, {22, -26}}, textString = "P"), Text(origin = {-34, 114}, extent = {{-16, 22}, {22, -28}}, textString = "T")}),
      Documentation(info = "<html>
  <body>
  <p>The main function of the FlowSource is to provide pressure, specific enthalpy and composition to a connection. In plus, elevation and flow can be added also to the connection. In the 'General' tab you can select if you want to use, or not, the connection as a flow and/or elevation source.</p>
  <p>In the sourceOption parameter, you can select between several calculation methods for pressure and specific enthalpy. Later, in the 'External values' tab you can select, for some variables only, if you want to use the values supplied through the Real input connectors, or the values entered manually at the 'Local values' tab. And in the 'Local values' tab you enter the values to use.</p>
  
  </body>
  </html>"));
  end FlowSource;

  model FlowSourceSP "Single phase flow source. Mainly necessary to generate pressure, enthalpy, and composition values at a connection"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium;
    replaceable FreeFluids.Interfaces.FluidPortB PortB(redeclare package Medium=Medium, P.displayUnit="bar") annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
    parameter FreeFluids.Types.SourceOptionS sourceOption = FreeFluids.Types.SourceOptionS.useP_T "how to generate port pressure and enthalpy";
    parameter Boolean isElevSource = true "use specified gravitational height";
    parameter Boolean isGsource = false "use specified mass flow";
    parameter Modelica.Units.SI.Height Elevation = 0.0 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.MassFlowRate G(displayUnit = "kg/h") = 0.0 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.AbsolutePressure P(displayUnit = "bar") = if sourceOption == FreeFluids.Types.SourceOptionS.useD_T then Medium.pressure(Medium.setState_dTX(D, T, X)) else Medium.p_default annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.Temperature T(displayUnit = "degC") annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.SpecificEnthalpy H = if sourceOption == FreeFluids.Types.SourceOptionS.useP_T then Medium.specificEnthalpy(Medium.setState_pTX(P, T, X)) elseif sourceOption == FreeFluids.Types.SourceOptionS.useD_T then Medium.specificEnthalpy(Medium.setState_dTX(D, T, X)) else Medium.h_default annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.Density D(displayUnit = "kg/m3") = 100 annotation(
      Dialog(tab = "Local values"));
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default annotation(
      Dialog(tab = "Local values"));
    /*conditional connectors*/
    parameter Boolean externalG = false annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalT = false annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalP = false annotation(
      Dialog(tab = "External values"));
    Modelica.Blocks.Interfaces.RealInput Gext(unit = "kg/s", displayUnit = "kg/h") if externalG == true annotation(
      Placement(visible = true, transformation(origin = {-60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {-60, 114}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealInput Text(unit = "K", displayUnit = "degC") if externalT == true annotation(
      Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 114}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealInput Pext(unit = "Pa", displayUnit = "bar") if externalP == true annotation(
      Placement(visible = true, transformation(origin = {60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {60, 114}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  protected
    Modelica.Blocks.Interfaces.RealInput GextIn "to connect with the conditiona input block";
    Modelica.Blocks.Interfaces.RealInput TextIn;
    Modelica.Blocks.Interfaces.RealInput PextIn;
  equation
    connect(Gext, GextIn);
    connect(Text, TextIn);
    connect(Pext, PextIn);
    if externalG == false then
      GextIn = G;
    end if;
    if externalT == false then
      TextIn = T;
    end if;
    if externalP == false then
      PextIn = P;
    end if;
    if isElevSource == true then
      PortB.Elevation = Elevation;
    end if;
    if isGsource == true then
      PortB.G = -GextIn;
    end if;
    if sourceOption == FreeFluids.Types.SourceOptionS.useP_H then
      PortB.P = PextIn;
      PortB.H = H;
    elseif sourceOption == FreeFluids.Types.SourceOptionS.useP_T then
      PortB.P = PextIn;
      PortB.H = Medium.specificEnthalpy(Medium.setState_pTX(PextIn, TextIn, X));
    elseif sourceOption == FreeFluids.Types.SourceOptionS.useD_T then
      PortB.P = Medium.pressure(Medium.setState_dTX(D, TextIn, X));
      PortB.H = Medium.specificEnthalpy(Medium.setState_dTX(D, TextIn, X));
    end if;
    PortB.X = X;
    annotation(
      defaultComponentName = "Source",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {0, 51, 154}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(origin = {-2, 8}, extent = {{-56, 42}, {60, -22}}, textString = "Source"), Text(origin = {-8, -38}, extent = {{-56, 42}, {36, -6}}, textString = "1 Ph"), Text(origin = {-88, 103}, extent = {{-34, 31}, {20, -15}}, textString = "G"), Text(origin = {-20, 103}, extent = {{-34, 31}, {20, -15}}, textString = "T"), Text(origin = {94, 103}, extent = {{-34, 31}, {20, -15}}, textString = "P")}),
      Documentation(info = "<html>
  <body>
  <p>It is simplified model of the Flow source model, to be used with single phase media.</p>
  </body>
  </html>"));
  end FlowSourceSP;

  model FreeReference
    replaceable FreeFluids.Interfaces.FluidPortB PortB annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}})));
  equation
    PortB.G = 0;
  end FreeReference;

  model PressureReference "Introduces a pressure reference for the connected ports"
    extends FreeReference;
    parameter Modelica.Units.SI.AbsolutePressure refP(displayUnit = "bar") = 1e5;
  equation
    PortB.P = refP;
    annotation(
      defaultComponentName = "Pref",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {170, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(extent = {{-64, 54}, {62, -46}}, textString = "P.ref")}));
  end PressureReference;

  model ElevationReference "Introduces a elevation reference for the connected ports"
    extends FreeReference;
    parameter Modelica.Units.SI.Height refElev;
  equation
    PortB.Elevation = refElev;
    annotation(
      defaultComponentName = "ElevRef",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {170, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(origin = {-16, -2}, extent = {{-54, 46}, {94, -40}}, textString = "El.Ref.")}));
  end ElevationReference;

  model PressureElevationReference
    extends FreeReference;
    parameter Modelica.Units.SI.AbsolutePressure refP(displayUnit = "bar") = 1e5;
    parameter Modelica.Units.SI.Height refElev = 0.0;
  equation
    PortB.P = refP;
    PortB.Elevation = refElev;
    annotation(
      defaultComponentName = "PElevRef",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {170, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-90, 90}, {90, -90}}, endAngle = 360), Text(lineColor = {0, 0, 255}, extent = {{-150, -90}, {150, -150}}, textString = "%name"), Text(origin = {-2, 18}, extent = {{-82, 56}, {96, -10}}, textString = "P,El."), Text(origin = {-4, -50}, extent = {{-82, 56}, {96, -10}}, textString = "ref.")}));
  end PressureElevationReference;

  model FlowSink "Just to specify the pressure, or the flow, at one boundary"
    replaceable package Medium = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    replaceable FreeFluids.Interfaces.FluidPortA PortA(redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
    parameter FreeFluids.Types.BoundaryOption fix = FreeFluids.Types.BoundaryOption.fixFlow "indicates if pressure or flow must be fixed";
    parameter Modelica.Units.SI.AbsolutePressure P(displayUnit = "bar") = 1e5;
    parameter Modelica.Units.SI.MassFlowRate G(displayUnit = "kg/h") = 0 "possitive if flow enters the sink";
    parameter Boolean externalG = false annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalP = false annotation(
      Dialog(tab = "External values"));
    Modelica.Blocks.Interfaces.RealInput Gext(unit = "kg/s", displayUnit = "kg/h") if externalG == true annotation(
      Placement(visible = true, transformation(origin = {-60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {-68, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealInput Pext(unit = "Pa", displayUnit = "bar") if externalP == true annotation(
      Placement(visible = true, transformation(origin = {60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {66, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  protected
    Modelica.Blocks.Interfaces.RealInput GextIn "to connect with the conditional input block";
    Modelica.Blocks.Interfaces.RealInput PextIn;
  equation
    connect(Gext, GextIn);
    connect(Pext, PextIn);
    if externalG==false then
        GextIn=G;
    end if;
    if externalP==false then
        PextIn=P;
    end if;
    if fix == FreeFluids.Types.BoundaryOption.fixPressure then
      PortA.P=PextIn;
    elseif fix == FreeFluids.Types.BoundaryOption.fixFlow then
      PortA.G=GextIn;
    end if;
    annotation(
      defaultComponentName = "Sink",
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Ellipse(lineColor = {0, 48, 144}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-90, -90}, {90, 90}}, endAngle = 360), Text(origin = {-33, 7}, extent = {{-45, 29}, {115, -33}}, textString = "Sink"), Text(origin = {54, -128}, lineColor = {0, 0, 255}, extent = {{-154, 40}, {44, -20}}, textString = "%name"), Text(origin = {-111, 115}, extent = {{-33, 23}, {33, -23}}, textString = "G"), Text(origin = {107, 115}, extent = {{-33, 23}, {33, -23}}, textString = "P")}));
  end FlowSink;

  connector HeatPort "for the connection of heat interchange. With additional information regarding a pipe"
    Modelica.Units.SI.Temperature T(displayUnit = "degC", start = 298.15) "Port temperature";
    flow Modelica.Units.SI.HeatFlowRate W(start = 0.0) "Heat flow. Possitive if entering the element";
    replaceable Modelica.Units.SI.Area S(start = 1) "surface of heat transfer";
    replaceable Modelica.Units.SI.Length L(start = 1) "individual tube length";
    replaceable Modelica.Units.SI.Diameter D(start = 0.05) "external pipe diameter";
    replaceable Fraction Slope "pipe slope H/L";
    annotation(
      defaultComponentName = "HeatPort",
      Documentation(info = "<html>
  <body>
  <p></p>
  </body>
  </html>"));
  end HeatPort;

  connector HeatPortA "will receive the causal variables"
    extends HeatPort(redeclare input Modelica.Units.SI.Area S, redeclare input Modelica.Units.SI.Length L, redeclare input Modelica.Units.SI.Diameter D, redeclare input Fraction Slope);
    annotation(
      Icon(graphics = {Rectangle(lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}, coordinateSystem(initialScale = 0.1)));
  end HeatPortA;

  connector HeatPortB "will supply the causal variables"
    extends HeatPort(redeclare output Modelica.Units.SI.Area S, redeclare output Modelica.Units.SI.Length L, redeclare output Modelica.Units.SI.Diameter D, redeclare output Fraction Slope);
    annotation(
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(origin = {-2, 4}, lineColor = {255, 0, 0}, lineThickness = 10, extent = {{-78, 78}, {82, -84}})}));
  end HeatPortB;

  model ThermalSource
    parameter Boolean isTsource = true "if true, temperature is forced, otherwise power is forced";
    parameter Modelica.Units.SI.Temperature T(displayUnit = "degC") = 298.15 annotation(
      Dialog(tab = "Local values"));
    parameter Modelica.Units.SI.HeatFlowRate W = 0 "if positive will inject heat in the connection" annotation(
      Dialog(tab = "Local values"));
    FreeFluids.Interfaces.HeatPortA PortH annotation(
      Placement(visible = true, transformation(origin = {1.77636e-15, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    /*conditional connectors*/
    parameter Boolean externalT = false "if true will accept the real connector value for T" annotation(
      Dialog(tab = "External values"));
    parameter Boolean externalW = false "if true will accept the real connector value for W" annotation(
      Dialog(tab = "External values"));
    Modelica.Blocks.Interfaces.RealInput Text if externalT == true annotation(
      Placement(visible = true, transformation(origin = {-106, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-106, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput Wext if externalW == true annotation(
      Placement(visible = true, transformation(origin = {-106, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-106, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  protected
    Modelica.Blocks.Interfaces.RealInput TextIn;
    Modelica.Blocks.Interfaces.RealInput WextIn;
  equation
    connect(Text, TextIn);
    connect(Wext, WextIn);
    if externalT == false then
      TextIn = T;
    end if;
    if externalW == false then
      WextIn = W;
    end if;
    if isTsource == true then
      PortH.T = TextIn;
    else
      PortH.W = -WextIn;
    end if;
    annotation(
      defaultComponentName = "ThSource",
      Icon(coordinateSystem(initialScale = 0.07), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {-116, 83}, extent = {{-24, 27}, {20, -15}}, textString = "T"), Text(origin = {-122, -89}, extent = {{-20, 23}, {22, -17}}, textString = "W"), Text(origin = {26, -149}, lineColor = {0, 0, 255}, extent = {{-130, 83}, {76, -41}}, textString = "%name"), Text(origin = {-54, 37}, extent = {{-24, 27}, {136, -61}}, textString = "Thermal"), Text(origin = {-40, -35}, extent = {{-24, 27}, {110, -41}}, textString = "Source")}));
  end ThermalSource;

  model ThermalBridge
    FreeFluids.Interfaces.HeatPortA PortA1 annotation(
      Placement(visible = true, transformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.HeatPortA PortA2 annotation(
      Placement(visible = true, transformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    PortA1.T = PortA2.T;
    0 = PortA1.W + PortA2.W;
    annotation(
      defaultComponentName = "ThBridge",
      Icon(graphics = {Text(origin = {-36, 41}, extent = {{-58, 31}, {130, -37}}, textString = "Thermal bridge"), Text(origin = {-11, -44}, lineColor = {0, 0, 255}, extent = {{-69, 60}, {89, -34}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
  end ThermalBridge;
end Interfaces;
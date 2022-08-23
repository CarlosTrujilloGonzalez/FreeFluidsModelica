within FreeFluids;

package Instruments
  "This file is part of the Free Fluids application
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
	
  model Reader "Allows to read Temperature, pressure,density, enthalpy and flow at PortA. PortB can be optionally connected"
    replaceable package Medium = FreeFluids.TMedia.Fluids.MarlothermSH constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
    Interfaces.FluidPortA PortA(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Interfaces.FluidPortB PortB(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput T(displayUnit="K") annotation(
      Placement(visible = true, transformation(origin = {-72, 44}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {-71, 45}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealOutput P(displayUnit="Pa") annotation(
      Placement(visible = true, transformation(origin = {-38, 60}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {-39, 75}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealOutput D(displayUnit="kg/m3") annotation(
      Placement(visible = true, transformation(origin = {38, 60}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {41, 75}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealOutput H(displayUnit="J/kg") annotation(
      Placement(visible = true, transformation(origin = {72, 44}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {71, 45}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealOutput G(displayUnit="kg/s") annotation(
      Placement(visible = true, transformation(origin = {-2.66454e-15, 84}, extent = {{-16, -16}, {16, 16}}, rotation = 90), iconTransformation(origin = {1, 87}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
    Medium.ThermodynamicState State;
  equation
    State = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    T = Medium.temperature(State);
    P = PortA.P;
    D = Medium.density(State);
    H = PortA.H;
    G = PortA.G;
    0 = PortA.G + PortB.G;
    PortA.P = PortB.P;
    PortA.H = PortB.H;
    PortA.X = PortB.X;
    PortA.Elevation = PortB.Elevation;
    annotation(
      defaultComponentName = "Reader1",
      Icon(graphics = {Line(points = {{-90, 0}, {90, 0}}, thickness = 1), Text(origin = {3, -22}, extent = {{-55, 20}, {55, -20}}, textString = "%name"), Text(origin = {-69, 22}, extent = {{-17, 8}, {17, -8}}, textString = "T"), Text(origin = {-39, 52}, extent = {{-17, 8}, {17, -8}}, textString = "P"), Text(origin = {43, 52}, extent = {{-17, 8}, {17, -8}}, textString = "D"), Text(origin = {71, 22}, extent = {{-17, 8}, {17, -8}}, textString = "H"), Text(origin = {1, 64}, extent = {{-17, 8}, {17, -8}}, textString = "G")}));
  end Reader;
  
  model ReaderExtended
    extends Reader;
    Modelica.Blocks.Interfaces.RealOutput Mu(displayUnit="Pa·s") annotation(
      Placement(visible = true, transformation(origin = {2.22045e-16, -80}, extent = {{-16, -16}, {16, 16}}, rotation = -90), iconTransformation(origin = {-4.44089e-16, -78}, extent = {{-12, -12}, {12, 12}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput Cp(displayUnit="J/(kg·K)") annotation(
      Placement(visible = true, transformation(origin = {60, -60}, extent = {{-16, -16}, {16, 16}}, rotation = -90), iconTransformation(origin = {60, -60}, extent = {{-12, -12}, {12, 12}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput K (displayUnit="W/(m·K)") annotation(
      Placement(visible = true, transformation(origin = {-60, -60}, extent = {{-16, -16}, {16, 16}}, rotation = -90), iconTransformation(origin = {-60, -60}, extent = {{-12, -12}, {12, 12}}, rotation = -90)));
  equation
    Mu = Medium.dynamicViscosity(State);
    Cp = Medium.specificHeatCapacityCp(State);
    K = Medium.thermalConductivity(State);
  annotation(
      defaultComponentName = "Reader1",
      Icon(graphics = {Text(origin = {-79, -60}, extent = {{-17, 8}, {17, -8}}, textString = "K"), Text(origin = {32, -80}, extent = {{-22, 8}, {22, -8}}, textString = "Mu"), Text(origin = {92, -60}, extent = {{-22, 8}, {22, -8}}, textString = "Cp")}));
  end ReaderExtended;
  
  model Reader2Ph
    extends Reader(redeclare replaceable package Medium=FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    Modelica.Blocks.Interfaces.RealOutput VQ annotation(
      Placement(visible = true, transformation(origin = {2.22045e-16, -84}, extent = {{-16, -16}, {16, 16}}, rotation = -90), iconTransformation(origin = {1, -87}, extent = {{-13, -13}, {13, 13}}, rotation = -90)));
  equation
    VQ=Medium.vapourQuality(State);
  annotation(
      Icon(graphics = {Text(origin = {1, -59}, extent = {{-71, 9}, {71, -9}}, textString = "Vap.Qual.")}));
  end Reader2Ph;

block PI "Proportional-Integral controller"
  extends Modelica.Blocks.Interfaces.SI2SO(y.start=0.5);
  import Modelica.Blocks.Types.Init;
  parameter Real k(unit="1")=1 "Gain. Output is 0..1";
  parameter Modelica.Units.SI.Time T(start=1,min=Modelica.Constants.small)=1
    "Integral time (seconds)";
  parameter Modelica.Blocks.Types.Init initType=Modelica.Blocks.Types.Init.NoInit
    "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)"
                                                                            annotation(Evaluate=true,
      Dialog(group="Initialization"));
  parameter Real x_start=0 "Initial or guess value of integral action"
    annotation (Dialog(group="Initialization"));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(enable=initType == Init.SteadyState or initType == Init.InitialOutput, group=
          "Initialization"));
  parameter Real range=300.0 "Process variable range";
  parameter Boolean directAction=false;
  parameter Boolean allowIntegralSaturation=true "if not allowed, simulation can be slower";
  Real error "PV-SP"; 
  Real x(start=x_start) "Integral action";

initial equation
  if initType == Init.SteadyState then
    der(x) = 0;
  elseif initType == Init.InitialState then
    x = x_start;
  elseif initType == Init.InitialOutput then
    y = y_start;
  end if;
equation
  error=if directAction==true then (u2-u1)/range else (u1-u2)/range;
  if allowIntegralSaturation==true then
    der(x) = error/T;
  else
    if y<1 and y>0 then
      der(x) = error/T;
    else
      der(x)=0;
    end if;
  end if;
  y = min(max(k*(x + error),0-1e-6),1+1e-6);
  annotation (defaultComponentName="PI",
    Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u and
the output y as <em>PI</em> system:
</p>
<pre>
               1
 y = k * (1 + ---) * u
              T*s
         T*s + 1
   = k * ------- * u
           T*s
</pre>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general model class <strong>TransferFunction</strong>
instead and model a PI SISO system with parameters<br>
b = {k*T, k}, a = {T, 0}.
</p>
<pre>
Example:

 parameter: k = 0.3,  T = 0.4

 results in:
             0.4 s + 1
    y = 0.3 ----------- * u
               0.4 s
</pre>

<p>
It might be difficult to initialize the PI component in steady state
due to the integrator part.
This is discussed in the description of package
<a href=\"modelica://Modelica.Blocks.Continuous#info\">Continuous</a>.
</p>

</html>"), Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
        Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
        Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}),
        Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
        Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{90, -80}, {68, -72}, {68, -88}, {90, -80}}),
        Line(points = {{-80, -80}, {-80, -20}, {60, 80}}, color = {0, 0, 127}),
        Text(lineColor = {192, 192, 192}, extent = {{0, 6}, {60, -56}}, textString = "PI"),
        Text(
          extent={{-150,-150},{150,-110}},
          textString="T=%T"), Text(origin = {-129, 114}, extent = {{-25, 18}, {25, -18}}, textString = "SP"), Text(origin = {-131, -106}, extent = {{-25, 18}, {25, -18}}, textString = "PV")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}}), graphics={
        Rectangle( lineColor={0,0,255},extent={{-60,60},{60,-60}}),
        Text(
          extent={{-68,24},{-24,-18}},
          textString="k"),
        Text(
          extent={{-32,48},{60,0}},
          textString="T s + 1"),
        Text(
          extent={{-30,-8},{52,-40}},
          textString="T s"),
        Line(points={{-24,0},{54,0}}),
        Line(points={{-100,0},{-60,0}}, color={0,0,255}),
        Line(points={{62,0},{100,0}}, color={0,0,255})}));

end PI;
  
  model SplitRange
  "1 Single Input / 2 Single Output continuous control block"
    extends Modelica.Blocks.Icons.Block;
  
    Modelica.Blocks.Interfaces.RealInput u1 "input from controller" annotation (Placement(
          transformation(extent={{-140,-20},{-100,20}})));
  
    Modelica.Blocks.Interfaces.RealOutput y1D "Upper half output for valve, direct action" annotation (Placement(
          transformation(extent={{100,60},{120,80}})));
      Modelica.Blocks.Interfaces.RealOutput y1I "Upper output for valve, reverse action" annotation (Placement(
          transformation(extent={{100,20},{120,40}})));
  
    Modelica.Blocks.Interfaces.RealOutput y2D "Lower half output for valve, direct action" annotation (Placement(
          transformation(extent={{100,-40},{120,-20}})));
    Modelica.Blocks.Interfaces.RealOutput y2I "lower half output for valve, reverse action" annotation (Placement(
          transformation(extent={{100,-80},{120,-60}})));
  
  equation
    y1D=min(max(0,u1*2-1),1);
    y1I=1-y1D;
    y2D=min(max(0,-u1*2+1),1);
    y2I=1-y2D;
    annotation (Documentation(info = "<html>
  <p>
  Block has two continuous Real input signals u1 and u2 and one
  continuous Real output signal y.
  </p>
  </html>"),
      Icon(graphics = {Line(origin = {-4, 1}, points = {{0, 47}, {0, -47}}), Text(origin = {-18, 46}, extent = {{-6, 4}, {6, -4}}, textString = "1"), Text(origin = {-18, -42}, extent = {{-6, 4}, {6, -4}}, textString = "0"), Text(origin = {-37, 4}, extent = {{-33, 4}, {33, -4}}, textString = "0.5"), Line(origin = {41, 36}, points = {{-41, -12}, {59, 14}}), Line(origin = {41, -31}, points = {{-41, 11}, {59, -19}}), Line(origin = {-4.24817, 2.64373}, points = {{-6, 0}, {6, 0}}), Text(origin = {-14, 87}, extent = {{-100, 9}, {100, -9}}, textString = "split range"), Text(origin = {86, 69}, extent = {{-8, 7}, {8, -7}}, textString = "D"), Text(origin = {88, -25}, extent = {{-8, 7}, {8, -7}}, textString = "D"), Text(origin = {86, 29}, extent = {{-8, 7}, {8, -7}}, textString = "I"), Text(origin = {86, -69}, extent = {{-8, 7}, {8, -7}}, textString = "I")}));
  end SplitRange;
end Instruments;

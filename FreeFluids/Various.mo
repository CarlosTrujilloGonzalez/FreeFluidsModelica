within FreeFluids;

package Various "Various.mo by Carlos Trujillo
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
  //import Media.*;
  import FreeFluids.Interfaces.*;
  import FreeFluids.Pipes.*;
  //***MISCELANEOUS CALCULATIONS***
  //*******************************
 
  import SI = Modelica.Units.SI;
  connector TemperatureOutput = output Modelica.Units.SI.Temperature "'output Real' as connector" annotation(
    defaultComponentName = "y",
    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100.0, -100.0}, {100.0, 100.0}}), graphics = {Polygon(lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{-100.0, 100.0}, {100.0, 0.0}, {-100.0, -100.0}})}),
    Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100.0, -100.0}, {100.0, 100.0}}), graphics = {Polygon(lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{-100.0, 50.0}, {0.0, 0.0}, {-100.0, -50.0}}), Text(textColor = {0, 0, 127}, extent = {{30.0, 60.0}, {30.0, 110.0}}, textString = "%name")}),
    Documentation(info = "<html>
  <p>
  Connector with one output signal of type Real.
  </p>
  </html>"));

end Various;
within;
package FreeFluids "FreeFluids.mo by Carlos Trujillo
      This file is part of the Free Fluids application
      Copyright (C) 2008-2024  Carlos Trujillo Gonzalez
        
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
  import SI = Modelica.Units.SI;
  import NonSI = Modelica.Conversions.NonSIunits;
  import pi = Modelica.Constants.pi;
  import g_n = Modelica.Constants.g_n;
  import R = Modelica.Constants.R;
  type Fraction = Real(min = 0.0, max = 1.0, nominal = 0.5, start=0.5);
  type MassFlux = Real(quantity = "MassFlux", final unit = "kg/(s.m2)");
  type Possitive = Real(min = 0.0);

  annotation(
    version="4.2.2",
    uses(Modelica(version = "4.0.0")),
  Documentation);
end FreeFluids;

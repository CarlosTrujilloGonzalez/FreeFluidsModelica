within FreeFluids;

package Types "Types.mo by Carlos Trujillo
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

  type Fraction = Real(min = 0.0, max = 1.0, nominal = 0.5, start=0.5);
  type MassFlux = Real(quantity = "MassFlux", final unit = "kg/(s.m2)");
  type Possitive = Real(min = 0.0);
  type ElevationOption=enumeration(differential "adds equation relation ports elevation", absolute "asigns the elevation value to PortB");
  type BoundaryOption=enumeration(fixPressure, fixFlow, fixNone);
  type ValveFixOption=enumeration(fixKv, fixDP, fixFlow);
  type ThermalType = enumeration(detailed "a detailed heat transfer model is used", isenthalpic , adiabatic , isothermal, fixedPower "a constant heat exchange is specified", fixedDeltaT);
  type CondensationOption = enumeration(partialCondensation "", totalCondensation , subcooling );

  type SourceOption = enumeration(useP_T, useSatLiqT, useSatLiqP, useSatGasT, useSatGasP, useP_H, useD_T );
  type SourceOptionS = enumeration(useP_T, useP_H, useD_T );
  type ExchangerOption = enumeration(fixedPower, fixedConductance, fixedTransferCoef);
  type TemaShell = enumeration(E "one pass",F "two pass longitudinal",G "split flow",H "double split flow",J "divided flow",K ,X "cross flow");
  type ExchangerType = enumeration(undefined, shellAndTubes, crossflow);
  type VesselForm = enumeration(cylinder, prism, sphere); 
  type HeadShape = enumeration(open, flat, conical, Klopper, Korbbogen,semielliptical, hemispherical);
  type EquilibriumDataType = enumeration(EquilibriumConstant, RelativeVolatility, RVpolynomia);
  type PerformanceSpecification = enumeration(FixedComposition, FixedFlow, FixedRecovery, None);
  type ReboilerType=enumeration(Kettle, Thermosyphon, ThermosyphonFromTray);
end Types;

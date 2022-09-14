within FreeFluids;

package HeatExchangers "HeatExchangers.mo by Carlos Trujillo
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
  //***EXCHANGERS***
  //****************
  //***NOT DETAILED EXCHANGER MODELS***
  //-----------------------------------

  model HEXsimple "Generic base for heat exchanger. Calculates pressure drop. Defines main variables"
    extends FreeFluids.Interfaces.TwoFluidPorts(useElevDifference = true, elevDifference = 0.0, calcEnthalpyDifference = true, passComposition = true);
    parameter Boolean useFixedDiffP = true "if true fixed dP will be used. Otherwise a  function of mass flow rate" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.PressureDifference dP(displayUnit = "bar") "pressure loss, constant or at reference flow rate. Negative if PortB.P < PortA.P" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.MassFlowRate refG(displayUnit = "kg/h", start = 1.0) "reference mass flow rate. If useFixedDiffP is false" annotation(
      Dialog(tab = "Flow"));
    parameter Boolean isLaminarFlow = false "flow regime to apply. Only if useFixedDiffP is false" annotation(
      Dialog(tab = "Flow"));
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer u annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fixSurface = true "if true, please fix the surface. If false, fix the power" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Area s = 0 "fixed heat exchange surface" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Power w "fixed heat transfer power. Positive if fluid inputs heat" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.SpecificHeatCapacity extCp = 1003 "fixed specific heat capacity of external media. Air=1003. Water=4194" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTin(displayUnit = "degC") = 298.15 "fixed external media inlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fixExternalFlow = false "if false the external media outlet temp. must be given. If true the external flow." annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.MassFlowRate extG(displayUnit = "kg/h") "fixed external media massic flow" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTout(displayUnit = "degC") = 318.15 "fixed external media outlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean counterCurrentFlow = true "if false, concurrent flow is considered" annotation(
      Dialog(tab = "Heat transfer"));
    Modelica.Units.SI.TemperatureDifference LMTD;
    Modelica.Units.SI.Area S;
    Modelica.Units.SI.MassFlowRate Gext(displayUnit = "kg/h") "external media mass flow rate";
    Modelica.Units.SI.Temperature TextOut(displayUnit = "degC");
    Modelica.Units.SI.ThermalConductance ThConductance "W/K";
    Modelica.Units.SI.Power W "heat transfer. Positive if fluid inputs heat";
    Medium.ThermodynamicState StateA "state at PortA";
    Medium.ThermodynamicState StateB "state at PortB";
    Medium.Temperature Ta;
    Medium.Temperature Tb;
  equation
    StateA = Medium.setState_phX(PortA.P, PortA.H, PortA.X);
    StateB = Medium.setState_phX(PortB.P, PortB.H, PortB.X);
    Ta = Medium.temperature(StateA);
    Tb = Medium.temperature(StateB);
//Flow
    if useFixedDiffP == true then
      Pdiff = dP;
    else
      if isLaminarFlow == true then
        Pdiff = dP * PortA.G / refG;
      else
        Pdiff = dP * PortA.G ^ 1.77 / refG ^ 1.77 "approximation with correction for the Reynolds in turbulent flow";
      end if;
    end if;
//Heat transfer
    if fixExternalFlow == true then
      Gext = extG;
    else
      TextOut = extTout;
    end if;
    if counterCurrentFlow == true then
      if noEvent((extTin - Tb) * (TextOut - Ta) <= 0) then
        LMTD = 0;
      elseif noEvent(extTin - Tb == TextOut - Ta) then
        LMTD = extTin - Tb;
      else
        LMTD = (extTin - Tb - TextOut + Ta) / log((extTin - Tb) / (TextOut - Ta));
      end if;
    else
      if noEvent((extTin - Ta) * (TextOut - Tb) <= 0) then
        LMTD = 0;
      elseif noEvent(extTin - Ta == TextOut - Tb) then
        LMTD = extTin - Ta;
      else
        LMTD = (extTin - Ta - TextOut + Tb) / log((extTin - Ta) / (TextOut - Tb));
      end if;
    end if;
    if fixSurface == true then
      S = s;
    else
      W = w;
    end if;
    ThConductance = S * u;
    W = ThConductance * LMTD;
    W = (extTin - TextOut) * extCp * Gext;
    W = PortA.G * (PortB.H - PortA.H);
    annotation(
      Dialog(tab = "Heat transfer"),
      defaultComponentName = "hex",
      Icon(graphics = {Line(origin = {-0.19, 0}, points = {{-89.8066, 0}, {-59.8066, -20}, {-19.8066, 20}, {20.1934, -20}, {60.1934, 20}, {90.1934, 0}}, color = {255, 0, 0}, thickness = 2), Line(origin = {-20, -18}, points = {{-50, -50}, {50, 50}}), Polygon(origin = {35.1208, 37.1208}, points = {{4.87921, -15.1208}, {-15.1208, 4.87921}, {14.8792, 14.8792}, {4.87921, -15.1208}}), Text(origin = {-52, 100}, extent = {{114, -42}, {-2, -2}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html><head></head><body>A very simple, and undefined, heat exchanger.<div>For the pressure drop you can consider it constant or variable with the 1.77 power of the flow.</div><div>Regarding heat transfer, you just have to supply the global heat transfer coefficient and the surface (or the power) of the exchanger, plus the inlet and outlet temperatures for the medium at the other side of the exchanger. The heat exchanged is calculated using the plain LMTD.</div></body></html>"));
  end HEXsimple;

  model HEXgeneric1Ph "Exchanger model without heat transfer coefficient calculation, for forced convection in tubes"
    extends FreeFluids.Pipes.PipeFlow1Ph(final useTubeLength = true, final lTube = tubeLength * numPasses, final fixNumTubes = true, final numTubes = div(numTubesTotal, numPasses), final pipeComplexity = 0, final kv = 0, final aperture = 1, final fullBore = true, final thicknessInsul = 0, final isCompressibleFlow = false, final twoPhaseFlow = false, thermalType = FreeFluids.Types.ThermalType.detailed);
    parameter FreeFluids.Types.ExchangerType exchangerType = FreeFluids.Types.ExchangerType.undefined "if undefined, raw LMTD will be used, otherwise mixed LMTD/NTU method" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean counterCurrent = true "for undefined type of exchanger select if flow is counter-current, otherwise will be considered co-current" annotation(
      Dialog(tab = "Heat transfer"));
    final parameter FreeFluids.Types.TemaShell shellType = FreeFluids.Types.TemaShell.E annotation(
      Dialog(tab = "Physical data"));
    parameter Modelica.Units.SI.Length tubeLength annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numTubesTotal annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numPasses = 1 annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numRows = 1 "necessary only for crossflow" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.SpecificHeatCapacity extCp = 1003 "fixed specific heat capacity of external media. Air=1003. Water=4194" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTin(displayUnit = "degC") = 298.15 "fixed external media inlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fixExternalFlow = true "if false the external media outlet temp. will be fixed" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.MassFlowRate extG(displayUnit = "kg/h") "fixed external media massic flow" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTout(displayUnit = "degC") = 318.15 "fixed external media outlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer u(start = 100) "fixed overall heat transfer coefficient referenced to external surface" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Real extSratio = (di + 2 * thickness) / di "ratio of external surface to internal bare tube surface" annotation(
      Dialog(tab = "Heat transfer"));
    Modelica.Units.SI.Area Sext;
    Modelica.Units.SI.MassFlowRate Gext "external media mass flow rate";
    Modelica.Units.SI.Temperature TextOut(displayUnit = "degC");
    Modelica.Units.SI.ThermalConductance ThConductance;
    Modelica.Units.SI.TemperatureDifference LMTD;
    Real Pntu(min = 0.0);
    Real Rntu(min = 0.0);
    Real NTU(min = 0.0, max = 1.0);
    Real Flmtd(min = 0.5, max = 1.0);
  equation
    if fixExternalFlow == true then
      Gext = extG;
    else
      TextOut = extTout;
    end if;
    Sext = SiActive * extSratio;
    ThConductance = u * Sext;
    W = (extTin - TextOut) * extCp * Gext;
    if exchangerType == FreeFluids.Types.ExchangerType.undefined then
      Pntu = 0;
      Rntu = 0;
      NTU = 0;
      Flmtd = 0;
      if counterCurrent == true then
        if noEvent((extTin - Tb) * (TextOut - Ta) <= 0) then
          LMTD = 0;
        elseif noEvent(extTin - Tb == TextOut - Ta) then
          LMTD = extTin - Tb;
        else
          LMTD = (extTin - Tb - TextOut + Ta) / log((extTin - Tb) / (TextOut - Ta));
        end if;
      else
        if noEvent((extTin - Ta) * (TextOut - Tb) <= 0) then
          LMTD = 0;
        elseif noEvent(extTin - Ta == TextOut - Tb) then
          LMTD = extTin - Ta;
        else
          LMTD = (extTin - Ta - TextOut + Tb) / log((extTin - Ta) / (TextOut - Tb));
        end if;
      end if;
      if thermalType == FreeFluids.Types.ThermalType.detailed then
        W = ThConductance * LMTD;
      end if;
    else
      Pntu = W / (PortA.G * Medium.specificHeatCapacityCp(StateAvg) * (extTin - Ta));
      Rntu = (TextOut - extTin) / (Ta - Tb);
      NTU = ThConductance / (PortA.G * Medium.specificHeatCapacityCp(StateAvg));
      if exchangerType == FreeFluids.Types.ExchangerType.shellAndTubes then
        Flmtd = FreeFluids.HeatExchangers.Functions.ShellLMTDfactor(shellType, 1, numPasses, Rntu, NTU);
      else
        Flmtd = FreeFluids.HeatExchangers.Functions.CrossLMTDfactor(numPasses, numRows, Rntu, NTU);
      end if;
      if noEvent((extTin - Tb) * (TextOut - Ta) <= 0) then
        LMTD = 0;
      elseif noEvent(extTin - Tb == TextOut - Ta) then
        LMTD = extTin - Tb;
      else
        LMTD = (extTin - Tb - TextOut + Ta) / log((extTin - Tb) / (TextOut - Ta));
      end if;
      if thermalType == FreeFluids.Types.ThermalType.detailed then
        W = ThConductance * Flmtd * LMTD;
      end if;
    end if;
    annotation(
      defaultComponentName = "HEX",
      Icon(graphics = {Line(origin = {-20, -18}, points = {{-60, -60}, {54, 54}}, color = {251, 2, 20}, thickness = 0.5), Polygon(origin = {39.12, 41.12}, fillColor = {250, 2, 17}, fillPattern = FillPattern.Solid, points = {{4.87921, -15.1208}, {-15.1208, 4.87921}, {14.8792, 14.8792}, {4.87921, -15.1208}})}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html><head></head><body>The model is for a tube heat exchanger that can have several passes and extended surface, in an undefined external flowing media. We need to specify the global U, as it is impossible its calculation. We define also the internal media. For the external media only Cp, inlet temperature, and outlet temperature (or mass flow rate) are needed.<div>As an alternative it is also possible to use a fixed exchanged power, or a fixed temperature for the tubes external surface. And, for better convergence when the exchanger is not in use, switch the pipes to isenthalpic or adiabatic.</div><div>If the exchanger type is undefined, a global LMTD will be used for the exchanged heat calculation, otherwise a mixed NTU/LTMD method will be used.</div></body></html>"));
  end HEXgeneric1Ph;

  model HEXgeneric2Ph "Exchanger model without heat transfer coefficient calculation, for forced convection in tubes"
    extends FreeFluids.Pipes.PipeFlow2Ph(final useTubeLength = true, final lTube = tubeLength * numPasses, final fixNumTubes = true, final numTubes = div(numTubesTotal, numPasses), final pipeComplexity = 0, final kv = 0, final aperture = 1, final fullBore = true, final thicknessInsul = 0, final isCompressibleFlow = true, final twoPhaseFlow = true, final rhoL = 0, final muL = 0, final rhoG = 0, final muG = 0, final x = 0, thermalType = FreeFluids.Types.ThermalType.detailed);
    parameter Modelica.Units.SI.Length tubeLength annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numTubesTotal annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numPasses = 1 annotation(
      Dialog(tab = "Physical data"));
    parameter Integer numRows = 1 "necessary only for crossflow" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.SpecificHeatCapacity extCp = 1003 "fixed specific heat capacity of external media. Air=1003. Water=4194" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTin(displayUnit = "degC") = 298.15 "fixed external media inlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean fixExternalFlow = true "if false the external media outlet temp. will be fixed" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.MassFlowRate extG(displayUnit = "kg/h") "fixed external media massic flow" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Temperature extTout(displayUnit = "degC") = 318.15 "fixed external media outlet temperature" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer u(start = 100) "fixed overall heat transfer coefficient referenced to external surface" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Real extSratio = (di + 2 * thickness) / di "ratio of external surface to internal bare tube surface" annotation(
      Dialog(tab = "Heat transfer"));
    Modelica.Units.SI.Area Sext;
    Modelica.Units.SI.MassFlowRate Gext "external media mass flow rate";
    Modelica.Units.SI.Temperature TextOut(displayUnit = "degC");
    Modelica.Units.SI.ThermalConductance ThConductance;
    Modelica.Units.SI.TemperatureDifference LMTD;
  equation
    if fixExternalFlow == true then
      Gext = extG;
    else
      TextOut = extTout;
    end if;
    Sext = SiActive * extSratio;
    ThConductance = u * Sext;
    W = (extTin - TextOut) * extCp * Gext;
    if noEvent((extTin - Tb) * (TextOut - Ta) <= 0) then
      LMTD = 0;
    elseif noEvent(extTin - Tb == TextOut - Ta) then
      LMTD = extTin - Tb;
    else
      LMTD = (extTin - Tb - TextOut + Ta) / log((extTin - Tb) / (TextOut - Ta));
    end if;
    if thermalType == FreeFluids.Types.ThermalType.detailed then
      W = ThConductance * LMTD;
    end if;
    annotation(
      defaultComponentName = "HEX",
      Icon(graphics = {Line(origin = {-20, -18}, points = {{-60, -60}, {54, 54}}, color = {249, 2, 2}, thickness = 0.5), Polygon(origin = {39.12, 41.12}, fillColor = {251, 0, 0}, fillPattern = FillPattern.Solid, points = {{4.87921, -15.1208}, {-15.1208, 4.87921}, {14.8792, 14.8792}, {4.87921, -15.1208}})}, coordinateSystem(initialScale = 0.1)),
      Documentation(info = "<html><head></head><body>The model is for a tube heat exchanger, that can have several passes and extended surface, in an undefined external flowing media. We need to specify the global U, as it is impossible its calculation. We define also the inner media, and the external media Cp, inlet temperature, and outlet temperature or mass flow rate.<div>As an alternative it is also possible to use a fixed exchanged power, or a fixed temperature for the tubes external surface.</div><div>As two phases are coexisting along the exchanger a straight LMTD method is used.</div><div><br></div></body></html>"));
  end HEXgeneric2Ph;

  //***DETAILED EXCHANGER MODELS***
  //-------------------------------
  //Shell with sensible heat transfer. tubes with condensing and subcooling
  //-----------------------------------------------------------------------

  package Functions
    extends Modelica.Icons.FunctionsPackage;

    function ShellLMTDfactor "As per VDI Atlas"
      input FreeFluids.Types.TemaShell shellType;
      input Integer numPassShell;
      input Integer numPassTubes;
      input Real R;
      input Real NTU;
      output Real F;
    protected
      Real NTUm, nTm;
    algorithm
      NTUm := NTU / numPassShell;
      nTm := numPassTubes / numPassShell;
      if shellType == FreeFluids.Types.TemaShell.E then
        if nTm == 1 then
          F := 1 / (1 + 0.671 * R ^ 1.055 * NTUm ^ 2.11) ^ 0.534;
        elseif nTm == 2 then
          F := 1 / (1 + 0.317 * R ^ 1.45 * NTUm ^ 2.09) ^ 0.543;
        elseif nTm == 3 then
          F := 1 / (1 + 0.431 * R ^ 1.0485 * NTUm ^ 2.33) ^ 0.371;
        elseif nTm == 4 then
          F := 1 / (1 + 0.274 * R ^ 1.05664 * NTUm ^ 2.08) ^ 0.624;
        end if;
      end if;
    end ShellLMTDfactor;

    function CrossLMTDfactor "As per VDI Atlas"
      input Integer numPass;
      input Integer numRows;
      input Real R;
      input Real NTU;
      output Real F;
    algorithm
      if numPass == 1 then
        if numRows == 1 then
          F := 1 / (1 + 0.234 * R ^ 1.27588 * NTU ^ 1.91) ^ 0.597;
        elseif numRows == 2 then
          F := 1 / (1 + 0.158 * R ^ 0.94401 * NTU ^ 1.53) ^ 0.705;
        elseif numRows == 3 then
          F := 1 / (1 + 0.15 * R ^ 0.82248 * NTU ^ 1.38) ^ 0.722;
        elseif numRows == 4 then
          F := 1 / (1 + 0.167 * R ^ 0.78122 * NTU ^ 1.34) ^ 0.648;
        elseif numRows == 5 then
          F := 1 / (1 + 0.195 * R ^ 0.76815 * NTU ^ 1.35) ^ 0.560;
        elseif numRows == 6 then
          F := 1 / (1 + 0.226 * R ^ 0.75465 * NTU ^ 1.37) ^ 0.486;
        end if;
      elseif numPass == 2 then
        F := 1 / (1 + 0.149 * R ^ 0.87472 * NTU ^ 1.76) ^ 0.264;
      elseif numPass == 3 then
        F := 1 / (1 + 0.0711 * R ^ 0.7807 * NTU ^ 1.85) ^ 0.253;
      elseif numPass == 4 then
        F := 1 / (1 + 0.0419 * R ^ 0.75411 * NTU ^ 1.89) ^ 0.246;
      else
        F := 1 / (1 + 0.251 * R ^ 1.03 * NTU ^ 2.06) ^ 0.677;
      end if;
//F:=1.0;
    end CrossLMTDfactor;
  end Functions;


end HeatExchangers;

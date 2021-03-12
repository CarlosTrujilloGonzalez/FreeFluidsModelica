within FreeFluids.TMedia;
package Tests
        partial model FluidTestingA
          replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          parameter Medium.AbsolutePressure p = 1.0e5;
          parameter Medium.Temperature initialT = 273.19;
          parameter Medium.Temperature finalT = 372.15;
          parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
          Medium.Temperature T(start = initialT) "We will ramp the temperature";
          Medium.MolarMass MM;
          Medium.ThermodynamicState StateP "original state from p and T";
          Medium.ThermodynamicState StateD "state reproduced from d,T";
          Medium.ThermodynamicState StateH "state reproduced from p,h";
          Medium.ThermodynamicState StateS "state reproduced from p,s";
          Medium.ThermodynamicState StateHalfH "state half way between Bub and Dew";
          //Properties of StateP
          Real Tb;
          Real H;
          Real D;
          Real S;
          Real U "internal energy";
          Real G "Gibbs energy";
          Real A "Helmholtz energy";
          Real Cp;
          Real Cv;
          Real Gamma;
          Real SS "speed of sound";
          Real Beta "isobaric expansion coefficient";
          Real Kappa "isothermal compressibility";
          //Real DerD_p_T;
          //Real DerD_T_p;
          Real DerD_p_h;
          Real DerD_h_p;
          Real Mu;
          Real Th;
          //Properties of StateDlow
          Medium.ThermodynamicState StateDlow "state at fraction of original density";
          Real DlowGasFract "gas fraction of StateDlow";
          Real DlowD, DlowH;
          Real DlowS;
          Real DlowU;
          Real DlowA;
          Real DlowCv;
          Real DlowDerD_p_h;
          Real DlowDerD_h_p;
          //Saturation properties
          Medium.SaturationProperties sat "saturation point at given T";
          Medium.AbsolutePressure Vp;
          Real DerTb_p;
          Medium.ThermodynamicState StateBub "bubble state at sat";
          Medium.ThermodynamicState StateDew "dew state at sat";
          Real BubD;
          Real BubH "bubble properties";
          Real BubS;
          Real BubDerD_p;
          Real BubDerH_p;
          Real DewD;
          Real DewH "dew properties";
          Real DewS;
          Real DewDerD_p;
          Real DewDerH_p;
          Real Hv "vaporization enthalpy";
          Real DewMu;
          Real Sigma;
          Medium.Temperature Tsat "temperature recovered from saturation pressure";
          //BaseProperties
          Medium.BaseProperties BaseProp;
        algorithm
        //Construction of StateP and calculation of properties
          for i in 1:1 loop
            StateP := Medium.setState_pTX(p, T);
            MM := Medium.molarMass(StateP);
            Tb := Medium.saturationTemperature(p);
            H := Medium.specificEnthalpy(StateP);
            D := Medium.density(StateP);
            S := Medium.specificEntropy(StateP);
            U := Medium.specificInternalEnergy(StateP);
            G := Medium.specificGibbsEnergy(StateP);
            A := Medium.specificHelmholtzEnergy(StateP);
            Cp := Medium.specificHeatCapacityCp(StateP);
            Cv := Medium.specificHeatCapacityCv(StateP);
            Gamma := Cp / Cv;
            SS := Medium.velocityOfSound(StateP);
            Beta := Medium.isobaricExpansionCoefficient(StateP);
            Kappa := Medium.isothermalCompressibility(StateP);
            DerD_p_h := Medium.density_derp_h(StateP);
            DerD_h_p := Medium.density_derh_p(StateP);
            Mu := Medium.dynamicViscosity(StateP);
            Th := Medium.thermalConductivity(StateP);
            StateD := Medium.setState_dTX(D, T);
            StateH := Medium.setState_phX(p, H);
            StateDlow := Medium.setState_dTX(D * fract, T, fill(0, 0));
            DlowGasFract := Medium.vapourQuality(StateDlow);
            DlowD := Medium.density(StateDlow);
            DlowH := Medium.specificEnthalpy(StateDlow);
            DlowU := Medium.specificInternalEnergy(StateDlow);
            DlowS := Medium.specificEntropy(StateDlow);
            DlowA := Medium.specificHelmholtzEnergy(StateDlow);
            DlowCv := Medium.specificHeatCapacityCv(StateDlow);
            DlowDerD_p_h := Medium.density_derp_h(StateDlow);
            DlowDerD_h_p := Medium.density_derh_p(StateDlow);
            sat := Medium.setSat_T(T);
            Vp := Medium.saturationPressure_sat(sat);
            DerTb_p := Medium.saturationTemperature_derp_sat(sat);
            StateBub := Medium.setBubbleState(sat);
            StateDew := Medium.setDewState(sat);
            BubD := Medium.bubbleDensity(sat);
            BubDerD_p := Medium.dBubbleDensity_dPressure(sat);
            BubH := Medium.bubbleEnthalpy(sat);
            BubDerH_p := Medium.dBubbleEnthalpy_dPressure(sat);
            BubS := Medium.bubbleEntropy(sat);
            DewD := Medium.dewDensity(sat);
            DewDerD_p := Medium.dDewDensity_dPressure(sat);
            DewH := Medium.dewEnthalpy(sat);
            DewDerH_p := Medium.dDewEnthalpy_dPressure(sat);
            Hv := DewH - BubH;
            DewS := Medium.dewEntropy(sat);
            DewMu := Medium.dynamicViscosity(StateDew);
            Tsat := Medium.saturationTemperature(sat.psat);
            Sigma := Medium.surfaceTension(sat);
          //DerD_p_T := Medium.density_derp_T(StateP);
          //DerD_T_p := Medium.density_derT_p(StateP);
            StateS := Medium.setState_psX(p, S);
            StateHalfH := Medium.setState_phX(Vp, (BubH + DewH) / 2, fill(0, 0));    
          end for;
        equation
        //Construction of BaseProperties
          BaseProp.p = p;
          BaseProp.T = T;
          der(T) = finalT - initialT;
          annotation(
            Documentation(info = "<html>
            <body>
            <p>This test maintains a constant pressure, and ramps temperature between the selected values. At each temperature, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p>
            </body>
            </html>"));
        end FluidTestingA;

        model TestA1A
          extends FluidTestingA(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), p = 20e5, initialT = 0.1 + 273.15, fract = 0.1, finalT = 380 + 273.15);
        end TestA1A;

        model TestA1B
          extends TestA1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_pT);
        end TestA1B;

        model TestA1C
          extends TestA1A(redeclare package Medium = FreeFluids.ExternalMedia.Fluids.WaterRef(refState = 4, reference_T = 273.16, reference_p = 101325.0, inputChoice = "pT", thermoModel = 3));
        end TestA1C;

    model TestA2A
      extends FluidTestingA(redeclare replaceable package Medium = TMedia.Fluids.R134A(highPressure = true, refState = "IIR", reference_T = 100.0, inputChoice = "pT"), p = 20.0e5, initialT = (-50.0) + 273.15, finalT = 110.0 + 273.15, fract = 0.5);
    end TestA2A;

        model TestA2B
          extends TestA2A(redeclare package Medium = Modelica.Media.R134a.R134a_ph);
        end TestA2B;

    model TestA2C
      extends TestA2A(redeclare replaceable package Medium = FreeFluids.ExternalMedia.Fluids.R134A(thermoModel = 3, refState = 2, inputChoice = "pT"));
    end TestA2C;

    partial model FluidTestingB
      replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
      parameter Medium.Temperature T = 273.2;
      parameter Medium.AbsolutePressure initialP = 1.0e5;
      parameter Medium.AbsolutePressure finalP = 1.0e7;
      parameter Real fract = 0.1 "fraction of the density of stateP to be used in state Dlow";
      Medium.AbsolutePressure p(start = initialP) "We will ramp the temperature";
      Medium.ThermodynamicState StateP "original state from p and T";
      Medium.ThermodynamicState StateH "state reproduced from p,h";
      Medium.ThermodynamicState StateS "state reproduced from p,s";
      Medium.ThermodynamicState StateD "state reproduced from d,T";
      Medium.ThermodynamicState StateDlow "state at fraction of original density";
      Medium.ThermodynamicState StateBub "bubble state at sat";
      Medium.ThermodynamicState StateDew "dew state at sat";
      //Properties of StateP
      Real D, H, Mu, Th;
      Real S;
      Real Beta "isobaric expansion coefficient";
      Real Kappa "isothermal compressibility";
      Real Cp;
      Real Cv;
      Real Gamma;
      Real SS "speed of sound";
      //Properties of StateDlow
      Real DlowGasFract "gas fraction of StateDlow";
      Real DlowD, DlowH;
      Real DlowS;
      //Saturation properties
      Medium.SaturationProperties sat "saturation point at given T";
      Real BubD, BubH "bubble properties";
      Real BubS;
      Real DewD, DewH "bubble properties";
      Real DewS;
      Real Hv "vaporization enthalpy";
      Real Sigma;
      Medium.Temperature Tsat "temperature recovered from saturation pressure";
      //BaseProperties
      Medium.BaseProperties BaseProp;
    algorithm
//Construction of StateP and calculation of properties
      StateP := Medium.setState_pTX(p, T, fill(0, 0));
      H := Medium.specificEnthalpy(StateP);
      D := Medium.density(StateP);
      S := Medium.specificEntropy(StateP);
      Cp := Medium.specificHeatCapacityCp(StateP);
      Cv := Medium.specificHeatCapacityCv(StateP);
      Gamma := Cp / Cv;
      SS := Medium.velocityOfSound(StateP);
      Beta := Medium.isobaricExpansionCoefficient(StateP);
      Kappa := Medium.isothermalCompressibility(StateP);
      Mu := Medium.dynamicViscosity(StateP);
      Th := Medium.thermalConductivity(StateP);
//Reconstruction of states from properties
      StateD := Medium.setState_dTX(D, T, fill(0, 0));
      StateH := Medium.setState_phX(p, H, fill(0, 0));
      StateS := Medium.setState_psX(p, S, fill(0, 0));
//Construction of state Dlow and get properties
      StateDlow := Medium.setState_dTX(D * fract, T, fill(0, 0));
      DlowGasFract := Medium.vapourQuality(StateDlow);
      DlowD := Medium.density(StateDlow);
      DlowH := Medium.specificEnthalpy(StateDlow);
      DlowS := Medium.specificEntropy(StateDlow);
//Calculations at saturation
      sat := Medium.setSat_T(T);
      StateBub := Medium.setBubbleState(sat);
      StateDew := Medium.setDewState(sat);
      BubD := Medium.bubbleDensity(sat);
      BubH := Medium.bubbleEnthalpy(sat);
      DewD := Medium.dewDensity(sat);
      DewH := Medium.dewEnthalpy(sat);
      Hv := DewH - BubH;
      BubS := Medium.bubbleEntropy(sat);
      DewS := Medium.dewEntropy(sat);
      Tsat := Medium.saturationTemperature(sat.psat);
      Sigma := Medium.surfaceTension(sat);
    equation
//Construction of BaseProperties
      BaseProp.p = p;
      BaseProp.h = H;
      der(p) = finalP - initialP;
      annotation(
        Documentation(info = "<html>
        <body>
        <p>This test maintains a constant temperature, and ramps pressure between the selected values. At each pressure, a thermodynamic state is created from p and T and, from it, different properties are calculated. The state is recovered also from d_T, p_H, and p_S, in order to check if the state is correctly reconstructed. </p>
        </body>
        </html>"));
    end FluidTestingB;

    model TestB1A
      extends FluidTestingB(redeclare replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = true, inputChoice = "pT"), T = 100 + 273.15, initialP = 100.0e5, fract = 1.0, finalP = 1.1e3);
    end TestB1A;

    model TestB1B
      extends TestB1A(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph "Medium model");
    end TestB1B;

    model TestWaterFlow1DFV_A_FF "Test case for Water.Flow1DFV. Copied from ThermoPower. Added full path to elements. Changed Medium to FreeFluids.TMedia.Fluids.Water"
      extends Modelica.Icons.Example;
      replaceable package Medium = FreeFluids.TMedia.Fluids.Water(refState = "User", reference_T = 273.15, highPressure = false) constrainedby Modelica.Media.Interfaces.PartialMedium;
      parameter Integer Nnodes = 20 "Number of nodes";
      parameter Modelica.SIunits.Length Lhex = 10 "Total length";
      parameter Modelica.SIunits.Diameter Dihex = 0.02 "Internal diameter";
      final parameter Modelica.SIunits.Radius rhex = Dihex / 2 "Internal radius";
      final parameter Modelica.SIunits.Length omegahex = Modelica.Constants.pi * Dihex "Internal perimeter";
      final parameter Modelica.SIunits.Area Ahex = Modelica.Constants.pi * rhex ^ 2 "Internal cross section";
      parameter Real Cfhex = 0.005 "Friction coefficient";
      parameter Modelica.SIunits.MassFlowRate whex = 0.31 "Nominal mass flow rate";
      parameter Modelica.SIunits.Pressure phex = 2e5 "Initial pressure";
      parameter Modelica.SIunits.SpecificEnthalpy hinhex = 1e5 "Initial inlet specific enthalpy";
      parameter Modelica.SIunits.SpecificEnthalpy houthex = 1e5 "Initial outlet specific enthalpy";
      parameter Modelica.SIunits.SpecificEnthalpy deltah = 41800 "Height of enthalpy step";
      parameter Modelica.SIunits.EnergyFlowRate W = 41800 * whex "Height of power step";
      Modelica.SIunits.Time tau "Transport time delay";
      ThermoPower.Water.SourceMassFlow fluidSource(p0 = phex, h = hinhex, w0 = whex, use_in_w0 = true, use_in_h = true, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-78, -10}, {-58, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure fluidSink(p0 = phex / 2, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{70, -10}, {90, 10}}, rotation = 0)));
      ThermoPower.Water.ValveLin valve(Kv = 3e-6, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{10, -10}, {30, 10}}, rotation = 0)));
      ThermoPower.Water.Flow1DFV hex(N = Nnodes, L = Lhex, omega = omegahex, Dhyd = Dihex, A = Ahex, wnom = whex, Cfnom = Cfhex, DynamicMomentum = false, hstartin = hinhex, hstartout = houthex, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, initOpt = ThermoPower.Choices.Init.Options.steadyState, pstart = phex, dpnom = 1000, redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-20, -10}, {0, 10}}, rotation = 0)));
      ThermoPower.Water.SensT T_in(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{-50, -6}, {-30, 14}}, rotation = 0)));
      ThermoPower.Thermal.HeatSource1DFV heatSource(Nw = Nnodes - 1) annotation(
        Placement(transformation(extent = {{-20, 22}, {0, 42}}, rotation = 0)));
      Modelica.Blocks.Sources.Step MassFlowRate(offset = whex, startTime = 50, height = -0.03) annotation(
        Placement(transformation(extent = {{-98, 28}, {-78, 48}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant Constant1(k = 1) annotation(
        Placement(transformation(extent = {{-10, 60}, {10, 80}}, rotation = 0)));
      Modelica.Blocks.Sources.Step InSpecEnthalpy(height = deltah, offset = hinhex, startTime = 1) annotation(
        Placement(transformation(extent = {{-90, 60}, {-70, 80}}, rotation = 0)));
      Modelica.Blocks.Sources.Step ExtPower(height = W, startTime = 30) annotation(
        Placement(transformation(extent = {{-40, 40}, {-20, 60}}, rotation = 0)));
      ThermoPower.Water.SensT T_out(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{40, -6}, {60, 14}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(transformation(extent = {{80, 80}, {100, 100}})));
    equation
      tau = sum(hex.rho) / Nnodes * Lhex * Ahex / whex;
      connect(hex.outfl, valve.inlet) annotation(
        Line(points = {{0, 0}, {6, 0}, {10, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(T_in.outlet, hex.infl) annotation(
        Line(points = {{-34, 0}, {-28, 0}, {-20, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(fluidSource.flange, T_in.inlet) annotation(
        Line(points = {{-58, 0}, {-46, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(heatSource.wall, hex.wall) annotation(
        Line(points = {{-10, 29}, {-10, 5}}, color = {255, 127, 0}));
      connect(T_out.outlet, fluidSink.flange) annotation(
        Line(points = {{56, 0}, {70, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(valve.outlet, T_out.inlet) annotation(
        Line(points = {{30, 0}, {44, 0}}, thickness = 0.5, color = {0, 0, 255}));
      connect(MassFlowRate.y, fluidSource.in_w0) annotation(
        Line(points = {{-77, 38}, {-72, 38}, {-72, 6}}, color = {0, 0, 127}));
      connect(InSpecEnthalpy.y, fluidSource.in_h) annotation(
        Line(points = {{-69, 70}, {-64, 70}, {-64, 6}}, color = {0, 0, 127}));
      connect(ExtPower.y, heatSource.power) annotation(
        Line(points = {{-19, 50}, {-10, 50}, {-10, 36}}, color = {0, 0, 127}));
      connect(Constant1.y, valve.cmd) annotation(
        Line(points = {{11, 70}, {20, 70}, {20, 8}}, color = {0, 0, 127}));
      annotation(
        Diagram(graphics),
        experiment(StopTime = 80, Tolerance = 1e-06, StartTime = 0, Interval = 0.16),
        Documentation(info = "<html>
    <p>The model is designed to test the component <code>Water.Flow1DFV</code> (fluid side of a heat exchanger, finite volumes).</p>
    <p>This model represent the fluid side of a heat exchanger with an applied external heat flow. The operating fluid is liquid water.</p>
    <p>During the simulation, the inlet specific enthalpy, heat flux and mass flow rate are changed. The outlet temperature can be predicted analytically assuming incompressible flow and constant cp.</p>
    <p><ul>
    <li>t=0 s, Step variation of the specific enthalpy of the fluid entering the heat exchanger. The outlet temperature should undergo a step increase of 10 degrees 10 s later. </li>
    <li>t=30 s, Step variation of the thermal flow entering the heat exchanger lateral surface. The outlet temperature should undergo a ramp increase of 10 degrees lasting 10 s </li>
    <li>t=50 s, Step reduction of the mass flow rate entering the heat exchanger. The outlet temperature should undergo a ramp change of one degree lasting 10s</li>
    </ul></p>
    <p>Simulation Interval = [0...80] sec </p>
    <p>Integration Algorithm = DASSL </p>
    <p>Algorithm Tolerance = 1e-6 </p>
    </html>", revisions = "<html>
    <p><ul>
    <li>12 Sep 2013 by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br/>Updated to new FV structure. Updated parameters.</li></li>
    <li><i>1 Oct 2003</i> by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br/>First release.</li>
    </ul></p>
    </html>"),
        uses(thermoPower(version = "3.1")));
    end TestWaterFlow1DFV_A_FF;

    model TestModel
      extends Modelica.Media.Examples.Tests.Components.PartialTestModel(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "pT"));
    equation

    end TestModel;

    model TestModel2
      extends Modelica.Media.Examples.Tests.Components.PartialTestModel2(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "pT"));
    equation

    end TestModel2;

    model ThreeTanks
      extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "ph"));
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4));
    end ThreeTanks;

    model HeatingSystem
      extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = FreeFluids.TMedia.Fluids.Water(inputChoice = "ph"));
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-6, Interval = 12));
    end HeatingSystem;
    
    model TestModelS
        extends Modelica.Media.Examples.Tests.Components.PartialTestModel(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
    end TestModelS;
    
    model TestModel2S
      extends Modelica.Media.Examples.Tests.Components.PartialTestModel2(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
    end TestModel2S;
    
    model ThreeTanksS
      extends Modelica.Fluid.Examples.Tanks.ThreeTanks(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.4));
    end ThreeTanksS;
    
    model HeatingSystemS
      extends Modelica.Fluid.Examples.HeatingSystem(redeclare package Medium = Modelica.Media.Water.StandardWater);
    equation
    
      annotation(
        experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-6, Interval = 12));
    end HeatingSystemS;
  end Tests;
within FreeFluids.Pipes;

package Examples
  extends Modelica.Icons.ExamplesPackage;

  package SimpleMedium
    extends Modelica.Media.Interfaces.PartialSimpleMedium(mediumName = "Simple", cp_const = 2850, cv_const = 2850, d_const = 930, eta_const = 0.08, lambda_const = 0.15, a_const = 1000, T_min = 273.15, T_max = 573.15, T0 = 273.15, MM_const = 3 "eta= dynamic visc. lambda= thermal cond. a= sound velocity, MM= molecular mas in kg");
  end SimpleMedium;

  package Glycol47 "1,2-Propylene glycol, 47% mixture with water"
    extends Modelica.Media.Incompressible.TableBased(mediumName = "Glycol-Water 47%", T_min = Modelica.Units.Conversions.from_degC(-30), T_max = Modelica.Units.Conversions.from_degC(100), TinK = false, T0 = 273.15, tableDensity = [-30, 1066; -20, 1062; -10, 1058; 0, 1054; 20, 1044; 40, 1030; 60, 1015; 80, 999; 100, 984], tableHeatCapacity = [-30, 3450; -20, 3490; -10, 3520; 0, 3560; 20, 3620; 40, 3690; 60, 3760; 80, 3820; 100, 3890], tableConductivity = [-30, 0.397; -20, 0.396; -10, 0.395; 0, 0.395; 20, 0.394; 40, 0.393; 60, 0.392; 80, 0.391; 100, 0.390], tableViscosity = [-30, 0.160; -20, 0.0743; -10, 0.0317; 0, 0.0190; 20, 0.00626; 40, 0.00299; 60, 0.00162; 80, 0.00110; 100, 0.00081], tableVaporPressure = [0, 500; 20, 1.9e3; 40, 5.3e3; 60, 16e3; 80, 37e3; 100, 80e3]);
    annotation(
      Documentation(info = "<html><head></head><body>Copied from the Modelica standard library in order to show how to define a simple table based medium.</body></html>"));
  end Glycol47;

  package Water1 = FreeFluids.TMedia.Fluids.Water(refState = "User", highPressure = false) "alias for TMedia water";
  package WaterS = Modelica.Media.Water.StandardWater;
  package Air2 = Modelica.Media.Air.DryAirNasa;
  package R134a1 = FreeFluids.TMedia.Fluids.R134A(refState = "User", reference_T = 100, highPressure = false);
  package R134aS = Modelica.Media.R134a.R134a_ph;
  package MarlothermSH = FreeFluids.TMedia.Fluids.MarlothermSH;

  model PhysicalPipeTest
    extends Modelica.Icons.Example;
    FreeFluids.Pipes.PipePhysical pipe(di = 0.05, isCircular = true, lTube = 6, numTubes = 3, thicknessInsul = 0.05, useDiameter = true, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end PhysicalPipeTest;

  model WaterFlow1PhIsothermal
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 1, G = 2, redeclare package Medium = Water1(highPressure = true), P(displayUnit = "Pa") = 2e+07, T(displayUnit = "K") = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-82, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1(highPressure = true), G = 2) annotation(
      Placement(visible = true, transformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Water1(highPressure = true), PLossFriction(displayUnit = "Pa"), di = 0.05, elevDifference = 1000, kv = 15, lTube = 1500, thermalType = FreeFluids.Types.ThermalType.isothermal, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Water1(highPressure = true), PLossFriction(displayUnit = "Pa"), di = 0.03, lTube = 100, thermalType = FreeFluids.Types.ThermalType.isothermal, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Pipe2.PortB, Sink.PortA) annotation(
      Line(points = {{40, 0}, {74, 0}, {74, 0}, {74, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Pipe2.PortA) annotation(
      Line(points = {{-14, 0}, {20, 0}, {20, 0}, {20, 0}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe1.PortA) annotation(
      Line(points = {{-72, 0}, {-34, 0}}, color = {0, 127, 255}));
  end WaterFlow1PhIsothermal;

  model WaterFlow1PhAbruptAdaptor
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Water1(highPressure = true), Elevation = 1, G = 2, P = 2e+07, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-82, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1(highPressure = true), G = 2) annotation(
      Placement(visible = true, transformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Water1(highPressure = true), PLossFriction(displayUnit = "Pa"), Ta(displayUnit = ""), Tb(displayUnit = ""), di = 0.05, elevDifference = 1000, kv = 15, lTube = 1500, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Water1(highPressure = true), Ta(displayUnit = ""), Tb(displayUnit = ""), di = 0.03, lTube = 100, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedM(redeclare package Medium = Water1(highPressure = true), diA = 0.05, diB = 0.03, rhoStart(displayUnit = "kg/m3") = 1000) annotation(
      Placement(visible = true, transformation(origin = {2, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedI(redeclare package Medium = Water1(highPressure = true), diB = 0.05, rhoStart(displayUnit = "kg/m3") = 1000) annotation(
      Placement(visible = true, transformation(origin = {-54, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedF(redeclare package Medium = Water1(highPressure = true), diA = 0.03, diB = Modelica.Constants.inf, rhoStart(displayUnit = "kg/m3") = 1000) annotation(
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  equation
    connect(RedF.PortB, Sink.PortA) annotation(
      Line(points = {{64, 0}, {74, 0}, {74, 0}, {74, 0}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, RedF.PortA) annotation(
      Line(points = {{40, 0}, {52, 0}, {52, 0}, {54, 0}}, color = {0, 127, 255}));
    connect(RedI.PortB, Pipe1.PortA) annotation(
      Line(points = {{-48, 0}, {-34, 0}, {-34, 0}, {-34, 0}}, color = {0, 127, 255}));
    connect(Source.PortB, RedI.PortA) annotation(
      Line(points = {{-72, 0}, {-60, 0}, {-60, 0}, {-58, 0}}, color = {0, 127, 255}));
    connect(RedM.PortB, Pipe2.PortA) annotation(
      Line(points = {{8, 0}, {20, 0}, {20, 0}, {20, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, RedM.PortA) annotation(
      Line(points = {{-14, 0}, {-2, 0}, {-2, 0}, {-4, 0}}, color = {0, 127, 255}));
  end WaterFlow1PhAbruptAdaptor;

  model WaterFlow1PhParallel
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP flowSource1(redeclare package Medium = WaterS, D = 1015, Elevation = 1, G = 2, H = 2.5e6, P = 500000, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = WaterS, G = 2) annotation(
      Placement(visible = true, transformation(origin = {34, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph pipe(redeclare package Medium = WaterS, di = 0.05, fixedDeltaT = 12, kv = 15, lTube = 20, thermalType = FreeFluids.Types.ThermalType.adiabatic, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph pipe1(redeclare package Medium = WaterS, di = 0.02, lTube = 30, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-6, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph pipe2(redeclare package Medium = WaterS, calcEnthalpyDifference = false, di = 0.02, lTube = 35, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-6, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.PortB, pipe2.PortA) annotation(
      Line(points = {{-32, 0}, {-16, 0}, {-16, -14}, {-16, -14}, {-16, -14}, {-16, -14}}, color = {0, 127, 255}));
    connect(pipe2.PortB, Sink.PortA) annotation(
      Line(points = {{4, -14}, {14, -14}, {14, -14}, {24, -14}, {24, 0}, {24, 0}, {24, 0}, {24, 0}}, color = {0, 127, 255}));
    connect(pipe.PortB, pipe1.PortA) annotation(
      Line(points = {{-32, 0}, {-24, 0}, {-24, 0}, {-16, 0}, {-16, 12}, {-16, 12}}, color = {0, 127, 255}));
    connect(pipe1.PortB, Sink.PortA) annotation(
      Line(points = {{4, 12}, {24, 12}, {24, 0}, {24, 0}}, color = {0, 127, 255}));
    connect(flowSource1.PortB, pipe.PortA) annotation(
      Line(points = {{-68, 0}, {-60, 0}, {-60, 0}, {-52, 0}, {-52, 0}, {-52, 0}, {-52, 0}, {-52, 0}}, color = {0, 127, 255}));
  end WaterFlow1PhParallel;

  model AirFlowAdiabaticParallel
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Air2, Elevation = 1, G = 2.36111, P(displayUnit = "Pa") = 700000, T(displayUnit = "K") = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-92, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(G = 2.36111, redeclare package Medium = Air2) annotation(
      Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), PortB(H(start = Source.H), P(start = Source.P)), di = 0.06, isCompressibleFlow = true, lTube = 15, roughness = 4.6e-5, thermalType = ThermalType.adiabatic, useElevDifference = true, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), di = 0.06, isCompressibleFlow = true, lTube = 20, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {0, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), calcEnthalpyDifference = false, di = 0.06, isCompressibleFlow = true, lTube = 22, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {0, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Pipe.PortB, Pipe1.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}, {-10, 14}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Pipe2.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}, {-10, -18}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Sink.PortA) annotation(
      Line(points = {{10, -18}, {40, -18}, {40, 0}, {40, 0}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-82, 0}, {-52, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Sink.PortA) annotation(
      Line(points = {{10, 14}, {40, 14}, {40, 0}, {40, 0}}, color = {0, 127, 255}));
//Pipe.Q = 1.988 / 60;
  end AirFlowAdiabaticParallel;

  model AirFlowAbruptAdaptor
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(G = 0.25, redeclare package Medium = Air2) annotation(
      Placement(visible = true, transformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), Ta(displayUnit = ""), Tb(displayUnit = ""), di = 0.05, elevDifference = 10, isCompressibleFlow = true, kv = 15, lTube = 100, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), Ta(displayUnit = ""), Tb(displayUnit = ""), di = 0.03, isCompressibleFlow = true, lTube = 50, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedM(redeclare package Medium = Air2, diA = 0.05, diB = 0.03, rhoStart(displayUnit = "kg/m3") = 7) annotation(
      Placement(visible = true, transformation(origin = {2, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedI(redeclare package Medium = Air2, diB = 0.05, rhoStart(displayUnit = "kg/m3") = 7) annotation(
      Placement(visible = true, transformation(origin = {-54, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    FreeFluids.Pipes.AbruptAdaptor RedF(redeclare package Medium = Air2, diA = 0.03, diB = Modelica.Constants.inf, rhoStart(displayUnit = "kg/m3") = 7) annotation(
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 1, G = 0.25, redeclare package Medium = Air2, P(displayUnit = "Pa") = 700000, T(displayUnit = "K") = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-86, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, RedI.PortA) annotation(
      Line(points = {{-76, 0}, {-60, 0}, {-60, 0}, {-58, 0}}, color = {0, 127, 255}));
    connect(RedF.PortB, Sink.PortA) annotation(
      Line(points = {{64, 0}, {74, 0}, {74, 0}, {74, 0}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, RedF.PortA) annotation(
      Line(points = {{40, 0}, {52, 0}, {52, 0}, {54, 0}}, color = {0, 127, 255}));
    connect(RedI.PortB, Pipe1.PortA) annotation(
      Line(points = {{-48, 0}, {-34, 0}, {-34, 0}, {-34, 0}}, color = {0, 127, 255}));
    connect(RedM.PortB, Pipe2.PortA) annotation(
      Line(points = {{8, 0}, {20, 0}, {20, 0}, {20, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, RedM.PortA) annotation(
      Line(points = {{-14, 0}, {-2, 0}, {-2, 0}, {-4, 0}}, color = {0, 127, 255}));
  end AirFlowAbruptAdaptor;

  model PipeFlow1PhAir3Test
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Air2, Elevation = 1, G = 1.38889, P(displayUnit = "bar") = 899999.9999999999, T(displayUnit = "K") = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Air2, G = 1.38889, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {36, -6.66134e-16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), di = 0.0703, isCompressibleFlow = true, kv = 0, lTube = 20, roughness = 4.6e-5, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), di = 0.06, isCompressibleFlow = true, lTube = 30, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-6, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Air2, PLossFriction(displayUnit = "Pa"), calcEnthalpyDifference = false, di = 0.06, isCompressibleFlow = true, lTube = 35, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-6, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Pipe.PortB, Pipe2.PortA) annotation(
      Line(points = {{-32, 0}, {-16, 0}, {-16, -14}, {-16, -14}, {-16, -14}, {-16, -14}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Sink.PortA) annotation(
      Line(points = {{4, -14}, {26, -14}, {26, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Pipe1.PortA) annotation(
      Line(points = {{-32, 0}, {-24, 0}, {-24, 0}, {-16, 0}, {-16, 12}, {-16, 12}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Sink.PortA) annotation(
      Line(points = {{4, 12}, {26, 12}, {26, 0}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-68, 0}, {-60, 0}, {-60, 0}, {-52, 0}, {-52, 0}, {-52, 0}, {-52, 0}, {-52, 0}}, color = {0, 127, 255}));
//Pipe.Q = 1.988 / 60;
  end PipeFlow1PhAir3Test;

  model WaterSPhReverseFlow
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Water1, D = 1015, Elevation = 1, G = 2, P = 500000, T = 298.15, externalP = true, externalT = true, isGsource = false) annotation(
      Placement(visible = true, transformation(origin = {-88, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, P = 400000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {62, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe(redeclare package Medium = Water1, di = 0.05, kv = 15, lTube = 20, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-58, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Water1, di = 0.04, lTube = 30, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {26, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Water1, di = 0.04, lTube = 35, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-16, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Valves.ValveIncompressible FV1(redeclare package Medium = Water1, aperture = 0, fixedKv = 10, useFixedAperture = false) annotation(
      Placement(visible = true, transformation(origin = {-14, 24}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    FreeFluids.Valves.CheckValve CV1(redeclare package Medium = Water1, Q(displayUnit = "m3/s"), calcEnthalpyDifference = false, fixedKv = 15, passComposition = false, useElevDifference = false) annotation(
      Placement(visible = true, transformation(origin = {38, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 0.1) annotation(
      Placement(visible = true, transformation(origin = {-42, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp2(duration = 1, height = -2e5, offset = 5e5) annotation(
      Placement(visible = true, transformation(origin = {-106, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const1(k = 25 + 273.15) annotation(
      Placement(visible = true, transformation(origin = {-124, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Instruments.Reader Reader1(redeclare package Medium = Water1) annotation(
      Placement(visible = true, transformation(origin = {14, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Instruments.ReaderExtended Reade2(redeclare package Medium = Water1) annotation(
      Placement(visible = true, transformation(origin = {11, 49}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    FreeFluids.Pipes.Dampener Volume(kv = 1.0) annotation(
      Placement(visible = true, transformation(origin = {-6, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  equation
    connect(const1.y, Source.Text) annotation(
      Line(points = {{-113, 30}, {-87, 30}, {-87, 10}, {-88, 10}, {-88, 10}, {-89, 10}}, color = {0, 0, 127}));
    connect(ramp2.y, Source.Pext) annotation(
      Line(points = {{-95, 58}, {-81, 58}, {-81, 10}, {-82, 10}, {-82, 10}, {-83, 10}}, color = {0, 0, 127}));
    connect(ramp1.y, FV1.Opening) annotation(
      Line(points = {{-31, 58}, {-14, 58}, {-14, 30}}, color = {0, 0, 127}));
    connect(CV1.PortB, Sink.PortA) annotation(
      Line(points = {{44, -26}, {52, -26}, {52, -2}}, color = {0, 127, 255}));
    connect(Pipe.PortB, FV1.PortA) annotation(
      Line(points = {{-48, -2}, {-21, -2}, {-21, 24}}, color = {0, 127, 255}));
    connect(FV1.PortB, Pipe1.PortA) annotation(
      Line(points = {{-7, 24}, {16, 24}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Pipe2.PortA) annotation(
      Line(points = {{-48, -2}, {-26, -2}, {-26, -26}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Sink.PortA) annotation(
      Line(points = {{36, 24}, {46, 24}, {46, 24}, {52, 24}, {52, 11}, {52, 11}, {52, -2}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-78, -2}, {-68, -2}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Reader1.PortA) annotation(
      Line(points = {{-6, -26}, {4, -26}}, color = {0, 127, 255}));
    connect(Reader1.PortB, CV1.PortA) annotation(
      Line(points = {{24, -26}, {32, -26}}, color = {0, 127, 255}));
    connect(FV1.PortB, Reade2.PortA) annotation(
      Line(points = {{-6, 24}, {2, 24}, {2, 49}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Volume.PortA) annotation(
      Line(points = {{-6, -26}, {-6, -56}}, color = {0, 127, 255}));
  end WaterSPhReverseFlow;

  model AirFlowChoked
    extends Modelica.Icons.Example;
    FreeFluids.Pipes.PipeFlowChoked Pipe(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 100000, lTube = 10, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, P(displayUnit = "Pa") = 2000000, T = 323.15, sourceOption = FreeFluids.Types.SourceOptionS.useP_T) annotation(
      Placement(visible = true, transformation(origin = {-106, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(P = 49999.99999999999, externalP = true, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp Ramp(duration = 500, height = -19e5, offset = 19.5e5) annotation(
      Placement(visible = true, transformation(origin = {18, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink2(P = 49999.99999999999, externalP = true, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {64, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2_1(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 25000, isCompressibleFlow = true, lTube = 2, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {-86, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2_2(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 25000, isCompressibleFlow = true, lTube = 2, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {-54, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2_3(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 25000, isCompressibleFlow = true, lTube = 2, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {-24, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlowChoked Pipe2_5(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 25000, lTube = 2, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {36, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FreeFluids.Pipes.PipeFlow1Ph Pipe2_4(redeclare package Medium = Modelica.Media.Air.ReferenceAir.Air_ph, di = 0.05000000000000001, fixedW(displayUnit = "kW") = 25000, isCompressibleFlow = true, lTube = 2, roughness(displayUnit = "mm") = 1.8e-06, thermalType = FreeFluids.Types.ThermalType.adiabatic) annotation(
      Placement(visible = true, transformation(origin = {6, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
  connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-96, 0}, {-22, 0}}, color = {0, 127, 255}));
  connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{-2, 0}, {40, 0}}));
  connect(Ramp.y, Sink.Pext) annotation(
      Line(points = {{29, 42}, {57, 42}, {57, 11}}, color = {0, 0, 127}));
    connect(Ramp.y, Sink2.Pext) annotation(
      Line(points = {{29, 42}, {71, 42}, {71, -25}}, color = {0, 0, 127}));
  connect(Pipe2_1.PortB, Pipe2_2.PortA) annotation(
      Line(points = {{-76, -36}, {-64, -36}}, color = {0, 127, 255}));
  connect(Pipe2_2.PortB, Pipe2_3.PortA) annotation(
      Line(points = {{-44, -36}, {-34, -36}}, color = {0, 127, 255}));
  connect(Pipe2_5.PortB, Sink2.PortA) annotation(
      Line(points = {{46, -36}, {54, -36}}, color = {0, 127, 255}));
  connect(Source.PortB, Pipe2_1.PortA) annotation(
      Line(points = {{-96, 0}, {-96, -36}}, color = {0, 127, 255}));
  connect(Pipe2_3.PortB, Pipe2_4.PortA) annotation(
      Line(points = {{-14, -36}, {-4, -36}}, color = {0, 127, 255}));
  connect(Pipe2_4.PortB, Pipe2_5.PortA) annotation(
      Line(points = {{16, -36}, {26, -36}}, color = {0, 127, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-120, 60}, {80, -60}})),
      experiment(StartTime = 0, StopTime = 500, Tolerance = 1e-6, Interval = 1),
  Documentation(info = "<html><head></head><body>An example of compressible flow that becomes choked. In Pipe2_x the length of Pipe has been splitted in five, in order to illustrate the influence of using average properties in the calculation of the pressure drop.</body></html>"));
  end AirFlowChoked;

  model WaterFlow2PhParallel
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Water1, D = 20, Elevation = 1, G = 0.1, H = 2.567e6, P = 200000, T = 388.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(G = 0.1, redeclare package Medium = Water1, externalG = true) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow2Ph Pipe(redeclare package Medium = Water1, PLossFrictionG(displayUnit = ""), PLossFrictionL(displayUnit = ""), di = 0.045, kv = 15, lTube = 20, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow2Ph Pipe1(redeclare package Medium = Water1, PLossFriction(displayUnit = "Pa"), PLossFrictionG(displayUnit = ""), PLossFrictionL(displayUnit = ""), di = 0.06, lTube = 30, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {14, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow2Ph Pipe2(redeclare package Medium = Water1, calcEnthalpyDifference = false, di = 0.06, lTube = 35, passComposition = false, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {14, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp FlowRamp(duration = 1, height = 0.1, offset = 0.1) annotation(
      Placement(visible = true, transformation(origin = {26, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Instruments.Reader2Ph Reader1 annotation(
      Placement(visible = true, transformation(origin = {-18, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Pipe2.PortB, Sink.PortA) annotation(
      Line(points = {{24, -14}, {50, -14}, {50, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Sink.PortA) annotation(
      Line(points = {{24, 14}, {50, 14}, {50, 0}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-68, 0}, {-58, 0}}, color = {0, 127, 255}));
    connect(FlowRamp.y, Sink.Gext) annotation(
      Line(points = {{37, 66}, {53, 66}, {53, 11}}, color = {0, 0, 127}));
    connect(Pipe.PortB, Reader1.PortA) annotation(
      Line(points = {{-38, 0}, {-28, 0}}, color = {0, 127, 255}));
    connect(Reader1.PortB, Pipe1.PortA) annotation(
      Line(points = {{-8, 0}, {4, 0}, {4, 14}}, color = {0, 127, 255}));
    connect(Reader1.PortB, Pipe2.PortA) annotation(
      Line(points = {{-8, 0}, {4, 0}, {4, -14}}, color = {0, 127, 255}));
  end WaterFlow2PhParallel;

  model WaterFlow2PhConstantPower
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Water1, D = 10, Elevation = 1, G = 0.0416667, H = 2.567e6, P = 200000, T = 423.15, sourceOption = FreeFluids.Types.SourceOption.useD_T) annotation(
      Placement(visible = true, transformation(origin = {-70, -1.77636e-15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, G = 0.0416667) annotation(
      Placement(visible = true, transformation(origin = {34, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow2Ph Pipe(redeclare package Medium = Water1, PLossFriction(displayUnit = "Pa"), PLossFrictionG(displayUnit = ""), PLossFrictionL(displayUnit = ""), di = 0.01, elevCalcMethod = FreeFluids.Types.ElevationOption.absolute, fixedW = 1e3, isCompressibleFlow = true, lTube = 2, portBelevation = -3.0, thermalType = FreeFluids.Types.ThermalType.fixedPower, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-60, 0}, {-30, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{-10, 0}, {24, 0}}, color = {0, 127, 255}));
  end WaterFlow2PhConstantPower;

  model SteamCondensing2Ph
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = WaterS, G = 0.04722222222222222, P = 180000, fix = FreeFluids.Types.BoundaryOption.fixFlow) annotation(
      Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow2Ph Pipe(redeclare package Medium = WaterS, PLossFrictionG(displayUnit = ""), PLossFrictionL(displayUnit = ""), Tb(start = 60 + 273), di = 0.05, elevDifference = -2.0, fixedW(displayUnit = "kW") = -99400, isCompressibleFlow = true, lTube = 10, thermalType = FreeFluids.Types.ThermalType.fixedPower, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = WaterS, Elevation = 1, G = 0.0277778, P = 2e5, T = 373.15, isGsource = false, sourceOption = FreeFluids.Types.SourceOption.useSatGasT) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {30, 0}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body>I using FreeFluids.TMedia you will need to use the old frontend. This is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 1));
  end SteamCondensing2Ph;

  model ForcedConvectionConstantPower
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.TMedia.Fluids.EG, D = 55, Elevation = 1, G(displayUnit = "kg/h") = 0.555556, H = 2.567e6, P = 1e+06, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-66, 4.44089e-16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.TMedia.Fluids.EG, G = 0.555556) annotation(
      Placement(visible = true, transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeForcedConvectionNew Pipe(redeclare package Medium = FreeFluids.TMedia.Fluids.EG, PLossFriction(displayUnit = "Pa"), di = 0.04, lTube = 60, thermalType = FreeFluids.Types.ThermalType.detailed, twistRatio = 20, useElevDifference = false, useThermalConnector = true, useTubeLength = true, useTwistedTapeInserts = true) annotation(
      Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ElevationReference ElevRef(refElev = 50) annotation(
      Placement(visible = true, transformation(origin = {-10, 44}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 283.15, W = 500, isTsource = false) annotation(
      Placement(visible = true, transformation(origin = {-20, -26}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  equation
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-56, 0}, {-30, 0}}, color = {0, 127, 255}));
    connect(ElevRef.PortB, Sink.PortA) annotation(
      Line(points = {{-2, 44}, {26, 44}, {26, 0}}, color = {0, 127, 255}));
    connect(ThSource.PortH, Pipe.PortH) annotation(
      Line(points = {{-20, -20}, {-20, -20}, {-20, -4}, {-20, -4}}, color = {255, 0, 0}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{-10, 0}, {26, 0}}, color = {0, 127, 255}));
  end ForcedConvectionConstantPower;

  model PipeWaterFallingFilmTest
    extends Modelica.Icons.Example;
    FreeFluids.Pipes.PipeFallingFilm Pipe(redeclare package Medium = Water1, PLossFriction(displayUnit = "Pa"), di = 23e-3, elevDifference = -1.0, foulingF = 0.000136, kWall = 17, lTube = 1, useThermalConnector = true, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {-2, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Water1, Elevation = 2, G = 0.0163333, P = 500000, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-54, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, G = 0.0163333) annotation(
      Placement(visible = true, transformation(origin = {48, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 403.15, isTsource = true) annotation(
      Placement(visible = true, transformation(origin = {-2, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  algorithm
//Pipe.Tb := 110 + 273.15;
  equation
    connect(ThSource.PortH, Pipe.PortH) annotation(
      Line(points = {{-2, -16}, {-2, -16}, {-2, 8}, {-2, 8}}, color = {255, 0, 0}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{8, 12}, {38, 12}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-44, 12}, {-12, 12}, {-12, 12}, {-12, 12}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body>Example from VDI Heat Atlas M3 Example 1. The heat transfer coefficient of the condensing vapor has been simulated by the fouling factor.</body></html>"));
  end PipeWaterFallingFilmTest;

  model SteamCondensingTotal
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(G = 0.04722222222222222, redeclare package Medium = Water1, P = 101000, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeCondensing Pipe(redeclare package Medium = Water1, PLossFriction(displayUnit = "Pa"), Tb(start = 60 + 273), TwallC(start = 20 + 273), condensationOption = FreeFluids.Types.CondensationOption.totalCondensation, di = 0.05, elevDifference = -5.0, lTube = 10, useElevDifference = true, useThermalConnector = true, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Water1, Elevation = 1, G = 0.0277778, P = 2e5, T = 373.15, isGsource = false, sourceOption = FreeFluids.Types.SourceOption.useSatGasT) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 350, isTsource = true) annotation(
      Placement(visible = true, transformation(origin = {0, -36}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  equation
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {30, 0}, {30, 0}, {30, 0}}, color = {0, 127, 255}));
    connect(ThSource.PortH, Pipe.PortH) annotation(
      Line(points = {{0, -30}, {0, -4}}, color = {255, 0, 0}));
    annotation(
      Documentation(info = "<html><head></head><body>Run with the old frontend. This is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst");
  end SteamCondensingTotal;

  model SteamCondensingPartial
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, G = 0.04722222222222222, P = 180000, fix = FreeFluids.Types.BoundaryOption.fixFlow) annotation(
      Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeCondensing Pipe(redeclare package Medium = Water1, Tb(start = 60 + 273), TwallC(start = 20 + 273), condensationOption = FreeFluids.Types.CondensationOption.partialCondensation, di = 0.05, elevDifference = -2.0, lTube = 10, useThermalConnector = true, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Water1, Elevation = 1, G = 0.0277778, P = 2e5, T = 373.15, isGsource = false, sourceOption = FreeFluids.Types.SourceOption.useSatGasT) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 350, isTsource = true) annotation(
      Placement(visible = true, transformation(origin = {0, -36}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  equation
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {30, 0}}, color = {0, 127, 255}));
    connect(ThSource.PortH, Pipe.PortH) annotation(
      Line(points = {{0, -30}, {0, -4}}, color = {255, 0, 0}));
    annotation(
      Documentation(info = "<html><head></head><body>Run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 1),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian,nonewInst");
  end SteamCondensingPartial;

  model SteamCondensingSubcooling
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = WaterS, G = 0.03333333333333333, P = 101100, fix = FreeFluids.Types.BoundaryOption.fixFlow) annotation(
      Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeCondSubcool Pipe(redeclare package Medium = WaterS, Tb(start = 60 + 273), TwallC(start = 20 + 273), di = 0.05, elevDifference = -10, lTube = 10, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = WaterS, Elevation = 1, P = 2e5, T = 373.15, isGsource = false, sourceOption = FreeFluids.Types.SourceOption.useSatGasT) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    Pipe.TsurfC = 76.85 + 273.15;
    Pipe.Ws = -3000;
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-32, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {30, 0}, {30, 0}, {30, 0}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
      Documentation(info = "<html><head></head><body><span style=\"font-family: 'Bitstream Vera Sans Mono'; font-size: 12px;\">If used with FreeFluids.TMedia, run with the old frontend. This is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</span></body></html>"));
  end SteamCondensingSubcooling;

  model ThermalPipesContact
    extends Modelica.Icons.Example;
    FreeFluids.Pipes.PipeForcedConvection Pipe1(redeclare package Medium = Air2, di = 0.019, lTube = 10, thickness = 0.001, useThermalConnector = true, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {0, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeForcedConvection Pipe2(redeclare package Medium = R134a1, PortB(P(start = 10e5)), di = 0.005, lTube = 10, thickness = 0.0005, useThermalConnector = true, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source1(redeclare package Medium = Air2, G = 0.222222, P(displayUnit = "bar") = 99999.99999999999, T(displayUnit = "degC") = 353.15, externalG = true, isGsource = true) annotation(
      Placement(visible = true, transformation(origin = {-58, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source2(redeclare package Medium = R134a1, G = 0.0555556, H = 250e3, P = 1500000, T(displayUnit = "degC") = 253.15, externalG = true, isGsource = true, sourceOption = FreeFluids.Types.SourceOption.useP_T) annotation(
      Placement(visible = true, transformation(origin = {-58, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink1(G = 0.222222, redeclare package Medium = Air2, fix = FreeFluids.Types.BoundaryOption.fixNone) annotation(
      Placement(visible = true, transformation(origin = {64, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink2(redeclare package Medium = R134a1, fix = FreeFluids.Types.BoundaryOption.fixNone) annotation(
      Placement(visible = true, transformation(origin = {64, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Interfaces.ThermalBridge thBridge annotation(
      Placement(visible = true, transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.001, offset = 0.0025) annotation(
      Placement(visible = true, transformation(extent = {{-110, -14}, {-90, 6}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = 0, offset = 0.01) annotation(
      Placement(visible = true, transformation(extent = {{-110, 58}, {-90, 78}}, rotation = 0)));
  equation
    connect(Pipe1.PortB, Sink1.PortA) annotation(
      Line(points = {{10, 26}, {54, 26}, {54, 26}, {54, 26}}));
    connect(Source2.PortB, Pipe2.PortA) annotation(
      Line(points = {{-48, -40}, {-10, -40}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Sink2.PortA) annotation(
      Line(points = {{10, -40}, {54, -40}}, color = {0, 127, 255}));
    connect(Source1.PortB, Pipe1.PortA) annotation(
      Line(points = {{-48, 26}, {-10, 26}}, color = {0, 127, 255}));
    connect(Pipe2.PortH, thBridge.PortA2) annotation(
      Line(points = {{0, -36}, {0, -36}, {0, -16}, {0, -16}}));
    connect(Pipe1.PortH, thBridge.PortA1) annotation(
      Line(points = {{0, 22}, {0, 22}, {0, 4}, {0, 4}}));
    connect(ramp1.y, Source2.Gext) annotation(
      Line(points = {{-89, -4}, {-66, -4}, {-66, -28}, {-64, -28}}, color = {0, 0, 127}));
    connect(ramp.y, Source1.Gext) annotation(
      Line(points = {{-88, 68}, {-64, 68}, {-64, 38}, {-64, 38}}, color = {0, 0, 127}));
  end ThermalPipesContact;

  model MarlothermSHInAir
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(G = 27.7778, redeclare package Medium = MarlothermSH) annotation(
      Placement(visible = true, transformation(origin = {36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 0, G = 27.7778, redeclare package Medium = MarlothermSH, P(displayUnit = "bar") = 499999.9999999999, T = 280 + 273) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeForcedConvection Pipe(redeclare package Medium = MarlothermSH, PLossFriction(displayUnit = "Pa"), PortB(H(start = Source.H), P(start = Source.P)), di = 0.324, kInsul = 0.063, lTube = 150 * 2, thickness = 3e-3, thicknessInsul = 0.1, useHTWallCorrFactor = false, useTubeLength = true, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {-14, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFluidBoundary externalHFMair1(emissionCoef = 0.35, foulingF = 0.0, tMedia = 293.15, vMedia = 4) annotation(
      Placement(visible = true, transformation(origin = {-14, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(externalHFMair1.PortH, Pipe.PortH) annotation(
      Line(points = {{-14, -27}, {-14, -14}}, color = {255, 0, 0}));
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{-4, -10}, {26, -10}, {26, -10}, {26, -10}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-64, -10}, {-24, -10}, {-24, -10}, {-24, -10}}, color = {0, 127, 255}));
  end MarlothermSHInAir;

  model WaterInAir
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink sink1(redeclare package Medium = Water1, G = 87750 / 3600) annotation(
      Placement(visible = true, transformation(origin = {36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP source1(Elevation = 0, G = 87750 / 3600, redeclare package Medium = Water1, P(displayUnit = "Pa") = 5e5, T(displayUnit = "K") = 343.15) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeForcedConvection pipe(redeclare package Medium = Water1, PLossFriction(displayUnit = "Pa"), di = 0.0889, kInsul = 0.052, lTube = 70, thickness = 3e-3, thicknessInsul = 0.03, useWallsResistance = true) annotation(
      Placement(visible = true, transformation(origin = {-14, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    PipeFluidBoundary Air(emissionCoef = 0.46, foulingF = 0, tMedia = 268, vMedia = 5) annotation(
      Placement(visible = true, transformation(origin = {-14, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.PortB, sink1.PortA) annotation(
      Line(points = {{-4, -10}, {26, -10}, {26, -10}, {26, -10}}, color = {0, 127, 255}));
    connect(source1.PortB, pipe.PortA) annotation(
      Line(points = {{-64, -10}, {-24, -10}, {-24, -10}, {-24, -10}}, color = {0, 127, 255}));
    connect(Air.PortH, pipe.PortH) annotation(
      Line(points = {{-14, -30}, {-14, -14}}, color = {255, 0, 0}));
  end WaterInAir;

  model SteamCondensInAir
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSink Sink(G = 0.01388888888888889, redeclare package Medium = WaterS) annotation(
      Placement(visible = true, transformation(origin = {38, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(Elevation = 0, redeclare package Medium = WaterS, P = 200000, T = 120 + 273, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeCondensing Pipe(redeclare package Medium = WaterS, PLossFriction(displayUnit = "Pa"), condensationOption = FreeFluids.Types.CondensationOption.partialCondensation, di = 0.02, emissionCoef = 0.46, kInsul = 0.052, lTube = 100, thickness = 0.001, thicknessInsul = 0, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-14, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFluidBoundary Boundary(redeclare package Medium = Air2, pMedia = 100000, tMedia = 323.15, vMedia = 0) annotation(
      Placement(visible = true, transformation(origin = {-14, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Pipe.PortB, Sink.PortA) annotation(
      Line(points = {{-4, -10}, {28, -10}}, color = {0, 127, 255}));
    connect(Source.PortB, Pipe.PortA) annotation(
      Line(points = {{-64, -10}, {-24, -10}, {-24, -10}, {-24, -10}}, color = {0, 127, 255}));
    connect(Boundary.PortH, Pipe.PortH) annotation(
      Line(points = {{-14, -30}, {-14, -14}}, color = {255, 0, 0}));
    annotation(
      Documentation(info = "<html><head></head><body>If used with FreeFluids.TMedia, run with the old frontend. This is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"));
  end SteamCondensInAir;

  model HalfCoilForcedConvection
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 0, G = 10, redeclare package Medium = Water1, P(displayUnit = "Pa") = 5e5, T(displayUnit = "K") = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-42, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Water1, G = 10) annotation(
      Placement(visible = true, transformation(origin = {42, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.HalfCoilForcedConvection HalfCoil(redeclare package Medium = Water1, PortB(H(start = 1e6)), basePipeDi = 0.05, num = 10, numActiveTubes = 4, numTubes = 4, path = 0.15, useThermalConnector = true) annotation(
      Placement(visible = true, transformation(origin = {-2, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 350, isTsource = true) annotation(
      Placement(visible = true, transformation(origin = {-2, -22}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  algorithm
    HalfCoil.HalfCoilDiam := 2;
  equation
    connect(ThSource.PortH, HalfCoil.PortH) annotation(
      Line(points = {{-2, -16}, {-2, -2}}, color = {255, 0, 0}));
    connect(HalfCoil.PortB, Sink.PortA) annotation(
      Line(points = {{8, 2}, {32, 2}}, color = {0, 127, 255}));
    connect(Source.PortB, HalfCoil.PortA) annotation(
      Line(points = {{-32, 2}, {-12, 2}}, color = {0, 127, 255}));
  end HalfCoilForcedConvection;

  model HalfCoilCondensing
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSource source(Elevation = 1, redeclare package Medium = WaterS, T = 373.15, isGsource = false, sourceOption = FreeFluids.Types.SourceOption.useSatGasT) annotation(
      Placement(visible = true, transformation(origin = {-42, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink sink(P = 99999.99999999999, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {40, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.HalfCoilCondensing HalfCoil(redeclare package Medium = WaterS, basePipeDi = 0.05, condensationOption = FreeFluids.Types.CondensationOption.totalCondensation, elevDifference = -1, num = 10, path = 0.15, thickness = 3e-3, useThermalConnector = false) annotation(
      Placement(visible = true, transformation(origin = {-2, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //basePipeDi = 0.05, num = 10, path = 0.15,
  algorithm
    HalfCoil.HalfCoilDiam := 2;
    HalfCoil.TwallC := 90 + 273.15;
  equation
    connect(HalfCoil.PortB, sink.PortA) annotation(
      Line(points = {{8, 4}, {30, 4}, {30, 4}, {30, 4}}, color = {0, 127, 255}));
    connect(source.PortB, HalfCoil.PortA) annotation(
      Line(points = {{-32, 4}, {-12, 4}, {-12, 4}, {-12, 4}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body>If used with FreeFluids.TMedia, run with the old frontend. this is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 1),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
  end HalfCoilCondensing;

  model CoilMarlothermSHThermal
    extends Modelica.Icons.Example;
    Interfaces.FlowSourceSP Source(Elevation = 0, G = 7.5, redeclare package Medium = MarlothermSH, P(displayUnit = "Pa") = 1e5, T(displayUnit = "K") = 443.15) annotation(
      Placement(visible = true, transformation(origin = {-78, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    CoilForcedConvection Coil(redeclare package Medium = MarlothermSH, coilDiam = 2.615, di = 0.064, num = 12, path = 0.14, rhoWall = 8000, thickness = 0.003, useThermalConnector = true) annotation(
      Placement(visible = true, transformation(origin = {-30, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Interfaces.FlowSink Sink(redeclare package Medium = MarlothermSH, G = 7.5) annotation(
      Placement(visible = true, transformation(origin = {38, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.ThermalSource ThSource(T = 483.15, isTsource = true) annotation(
      Placement(visible = true, transformation(origin = {-30, -18}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  equation
    connect(ThSource.PortH, Coil.PortH) annotation(
      Line(points = {{-30, -13}, {-30, 10}}, color = {255, 0, 0}));
    connect(Coil.PortB, Sink.PortA) annotation(
      Line(points = {{-20, 14}, {28, 14}}, color = {0, 127, 255}));
    connect(Source.PortB, Coil.PortA) annotation(
      Line(points = {{-68, 14}, {-40, 14}, {-40, 14}, {-40, 14}}, color = {0, 127, 255}));
  end CoilMarlothermSHThermal;

  model ModelicaPipeTestS
    extends Modelica.Fluid.Examples.IncompressibleFluidNetwork(valveOpening1.height = -0.99, each heat8.Q_flow = 0, redeclare package Medium = Modelica.Media.Water.StandardWater);
    annotation(
      experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-6, Interval = 0.4),
      Documentation(info = "<html><head></head><body>The Modelica.Fluid.Examples.IncompressibleFluidNetwork run with Modelica.Media.Water.StandardWater.</body></html>"));
  end ModelicaPipeTestS;

  model ModelicaPipeTestF
    extends Modelica.Fluid.Examples.IncompressibleFluidNetwork(valveOpening1.height = -0.99, each heat8.Q_flow = 0, redeclare package Medium = FreeFluids.LMedia.Fluids.Water(inputChoice = "ph"));
    annotation(
      experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-6, Interval = 0.4),
      Documentation(info = "<html><head></head><body>The Modelica IncompressibleFluidNetwor example ran with FreeFluids LMedia Water, in order to see the influence of the medium used.</body></html>"));
  end ModelicaPipeTestF;

  model ModelicaPipeTest
    extends Modelica.Fluid.Examples.IncompressibleFluidNetwork(valveOpening1.height = -0.99, each heat8.Q_flow = 0);
    annotation(
      experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-6, Interval = 0.4),
      Documentation(info = "<html><head></head><body>The Modelica.Fluid.Examples.IncompressibleFluidNetwork run with Modelica.Media.Water.StandardWater.</body></html>"));
  end ModelicaPipeTest;

  model ModelicaPipeTestFF
    extends Modelica.Icons.Example;
    FreeFluids.Interfaces.FlowSourceSP Source(Elevation = 0.0, redeclare package Medium = Glycol47, P = 5e+5, T(displayUnit = "K") = 300, isGsource = false) annotation(
      Placement(visible = true, transformation(origin = {-102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe1(redeclare package Medium = Glycol47, di = 0.025, lTube = 10, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe2(redeclare package Medium = Glycol47, di = 0.025, lTube = 0.5, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-48, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe3(redeclare package Medium = Glycol47, di = 0.025, lTube = 0.5, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-48, -22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    FreeFluids.Valves.ValveIncompressible Valve1(redeclare package Medium = Glycol47, aperture = 1, fixedKv = 17.53, useFixedAperture = false) annotation(
      Placement(visible = true, transformation(origin = {-40, 56}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Blocks.Sources.Step valveOpening1(height = -0.999, offset = 1, startTime = 50) annotation(
      Placement(visible = true, transformation(extent = {{-78, 90}, {-58, 70}}, rotation = 0)));
    FreeFluids.Valves.ValveIncompressible Valve2(redeclare package Medium = Glycol47, fixedKv = 17.53, useFixedAperture = false) annotation(
      Placement(visible = true, transformation(origin = {-40, -46}, extent = {{-7, 7}, {7, -7}}, rotation = 0)));
    Modelica.Blocks.Sources.Step valveOpening2(height = -0.5, offset = 1, startTime = 100) annotation(
      Placement(visible = true, transformation(extent = {{-80, -60}, {-60, -80}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Glycol47, fix = FreeFluids.Types.BoundaryOption.fixPressure) annotation(
      Placement(visible = true, transformation(origin = {98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe4(redeclare package Medium = Glycol47, calcEnthalpyDifference = false, di = 0.025, lTube = 2, passComposition = false, roughness = 2.5e-05, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-16, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe5(redeclare package Medium = Glycol47, calcEnthalpyDifference = false, di = 0.025, lTube = 20, passComposition = false, roughness = 2.5e-05, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {20, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe6(redeclare package Medium = Glycol47, calcEnthalpyDifference = false, di = 0.025, lTube = 20, passComposition = false, roughness = 2.5e-05, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {20, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe7(redeclare package Medium = Glycol47, di = 0.025, lTube = 10, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-16, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe8(redeclare package Medium = Glycol47, PLossFriction(displayUnit = "Pa"), calcEnthalpyDifference = true, di = 0.025, fixedW = 0, lTube = 10, passComposition = true, roughness = 2.5e-05, thermalType = FreeFluids.Types.ThermalType.fixedPower, useElevDifference = true, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {-6, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe9(redeclare package Medium = Glycol47, di = 0.025, lTube = 10, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {20, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe10(redeclare package Medium = Glycol47, calcEnthalpyDifference = false, di = 0.025, lTube = 10, passComposition = false, roughness = 2.5e-05, useElevDifference = false, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {20, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Pipes.PipeFlow1Ph Pipe11(redeclare package Medium = Glycol47, di = 0.025, lTube = 0.5, roughness = 2.5e-05, useTubeLength = true) annotation(
      Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Valves.ValveIncompressible Valve3(redeclare package Medium = Glycol47, Q(displayUnit = "m3/s"), fixedKv = 17.53, useFixedAperture = false) annotation(
      Placement(visible = true, transformation(origin = {74, 0}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Modelica.Blocks.Sources.Step valveOpening3(height = -0.5, offset = 1, startTime = 150) annotation(
      Placement(visible = true, transformation(extent = {{42, 88}, {62, 68}}, rotation = 0)));
  equation
//Pipe8.W = 16000;
    connect(valveOpening2.y, Valve2.Opening) annotation(
      Line(points = {{-58, -70}, {-40, -70}, {-40, -52}}, color = {0, 0, 127}));
    connect(Source.PortB, Pipe1.PortA) annotation(
      Line(points = {{-92, 0}, {-76, 0}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Pipe2.PortA) annotation(
      Line(points = {{-56, 0}, {-48, 0}, {-48, 20}, {-48, 20}}, color = {0, 127, 255}));
    connect(Pipe1.PortB, Pipe3.PortA) annotation(
      Line(points = {{-56, 0}, {-48, 0}, {-48, -12}, {-48, -12}}, color = {0, 127, 255}));
    connect(Pipe3.PortB, Valve2.PortA) annotation(
      Line(points = {{-48, -32}, {-47, -32}, {-47, -46}}, color = {0, 127, 255}));
    connect(Pipe2.PortB, Valve1.PortA) annotation(
      Line(points = {{-48, 40}, {-47.5, 40}, {-47.5, 56}, {-47, 56}}, color = {0, 127, 255}));
    connect(Valve1.PortB, Pipe7.PortA) annotation(
      Line(points = {{-32, 56}, {-26, 56}}, color = {0, 127, 255}));
    connect(Pipe7.PortB, Pipe9.PortA) annotation(
      Line(points = {{-6, 56}, {10, 56}}, color = {0, 127, 255}));
    connect(Valve2.PortB, Pipe4.PortA) annotation(
      Line(points = {{-32, -46}, {-27, -46}, {-27, -74}, {-26, -74}}, color = {0, 127, 255}));
    connect(Pipe7.PortB, Pipe8.PortA) annotation(
      Line(points = {{-6, 56}, {-6, 28}}, color = {0, 127, 255}));
    connect(Pipe9.PortB, Pipe11.PortA) annotation(
      Line(points = {{30, 56}, {40, 56}, {40, 0}}, color = {0, 127, 255}));
    connect(Pipe10.PortB, Pipe11.PortA) annotation(
      Line(points = {{30, -24}, {40, -24}, {40, 0}}, color = {0, 127, 255}));
    connect(Pipe6.PortB, Pipe11.PortA) annotation(
      Line(points = {{30, -46}, {40, -46}, {40, 0}}, color = {0, 127, 255}));
    connect(Pipe5.PortB, Pipe11.PortA) annotation(
      Line(points = {{30, -68}, {40, -68}, {40, 0}}, color = {0, 127, 255}));
    connect(Pipe11.PortB, Valve3.PortA) annotation(
      Line(points = {{60, 0}, {67, 0}}, color = {0, 127, 255}));
    connect(Valve3.PortB, Sink.PortA) annotation(
      Line(points = {{81, 0}, {88, 0}}, color = {0, 127, 255}));
    connect(valveOpening3.y, Valve3.Opening) annotation(
      Line(points = {{64, 78}, {74, 78}, {74, 6}, {74, 6}}, color = {0, 0, 127}));
    connect(valveOpening1.y, Valve1.Opening) annotation(
      Line(points = {{-56, 80}, {-40, 80}, {-40, 62}, {-40, 62}}, color = {0, 0, 127}));
    connect(Pipe8.PortB, Pipe10.PortA) annotation(
      Line(points = {{-6, 8}, {10, 8}, {10, -24}, {10, -24}}, color = {0, 127, 255}));
    connect(Pipe8.PortB, Pipe6.PortA) annotation(
      Line(points = {{-6, 8}, {10, 8}, {10, -46}, {10, -46}, {10, -46}}, color = {0, 127, 255}));
    connect(Pipe8.PortB, Pipe5.PortA) annotation(
      Line(points = {{-6, 8}, {10, 8}, {10, -68}, {10, -68}}, color = {0, 127, 255}));
    connect(Pipe4.PortB, Pipe5.PortA) annotation(
      Line(points = {{-6, -74}, {10, -74}, {10, -68}, {10, -68}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-6, Interval = 0.4),
      Documentation(info = "<html><head></head><body>An equivalent circuit, using FreeFluids models, of the Modelica.Fluid.Examples.IncompressibleFluidNetwork.</body></html>"));
  end ModelicaPipeTestFF;
end Examples;

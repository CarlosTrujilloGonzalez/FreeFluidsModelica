within FreeFluids.HeatExchangers;

package Examples
  extends Modelica.Icons.ExamplesPackage;

  model HEXsimpleWater
    extends Modelica.Icons.Example;
    FreeFluids.HeatExchangers.HEXsimple HEX(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, dP(displayUnit = "bar") = -12750, extCp = 1000, extTin = 383.15, extTout = 353.15, refG = 19.44444444444444, s = 31.1646, u = 400, useFixedDiffP = false) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, P(displayUnit = "bar") = 600000, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, G = 19.44444444444444) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, HEX.PortA) annotation(
      Line(points = {{-42, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(HEX.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {50, 0}, {50, 0}, {50, 0}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body><br></body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
  end HEXsimpleWater;

  model HEXgeneric1PhWater
    extends Modelica.Icons.Example;
    FreeFluids.HeatExchangers.HEXgeneric1Ph HEX(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, di = 0.02058, exchangerType = FreeFluids.Types.ExchangerType.shellAndTubes, extCp = 4194, extG(displayUnit = "kg/h"), extTin = 383.15, extTout = 353.15, fixExternalFlow = false, numPasses = 2, numTubesTotal = 80, thickness = 0.00211, tubeLength = 5, u = 400) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, P(displayUnit = "Pa") = 600000, T = 298.15) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.LMedia.Fluids.Water, G = 19.44444444444444) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, HEX.PortA) annotation(
      Line(points = {{-42, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(HEX.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {50, 0}, {50, 0}, {50, 0}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body><br></body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
  end HEXgeneric1PhWater;

  model HEXgeneric1PhThermalOil
    extends Modelica.Icons.Example;
    FreeFluids.HeatExchangers.HEXgeneric1Ph HEX(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, di = 0.02058, exchangerType = FreeFluids.Types.ExchangerType.crossflow, extCp = 1006.01, extG(displayUnit = "kg/h") = 22.3439, extSratio = 28.35, extTin = 305.15, extTout = 333.15, fixExternalFlow = true, numPasses = 3, numRows = 6, numTubesTotal = 132, numVelocityHeads = 4.5, thickness = 0.00241, tubeLength = 4, u = 16.155) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSourceSP Source(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, P(displayUnit = "Pa") = 499999.9999999999, T = 375.15) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = FreeFluids.LMedia.Fluids.MarlothermSH, G = 27.5278) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, HEX.PortA) annotation(
      Line(points = {{-42, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(HEX.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {50, 0}, {50, 0}, {50, 0}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body><br></body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
  end HEXgeneric1PhThermalOil;

  model HEXgeneric2PhSteamCondensing
    extends Modelica.Icons.Example;
    FreeFluids.HeatExchangers.HEXgeneric2Ph HEX(redeclare package Medium = Modelica.Media.Water.StandardWater, PLossFrictionG(displayUnit = ""), PLossFrictionL(displayUnit = ""), di = 0.02058, extCp = 1003, extG(displayUnit = "kg/s") = 22.3439, extSratio = 28, extTin = 305.15, extTout = 333.15, fixExternalFlow = true, isCircular = true, numPasses = 3, numRows = 6, numTubesTotal = 132, numVelocityHeads = 7, perimeter = 0.0735, section = 0.0002851, thickness = 0.00211, tubeLength = 4, u = 11.39, useDiameter = true, useSectionAndPerimeter = false) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSource Source(redeclare package Medium = Modelica.Media.Water.StandardWater, P = 200000, sourceOption = FreeFluids.Types.SourceOption.useSatGasP) annotation(
      Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FlowSink Sink(redeclare package Medium = Modelica.Media.Water.StandardWater, G = 0.5555555555555556) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Source.PortB, HEX.PortA) annotation(
      Line(points = {{-42, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(HEX.PortB, Sink.PortA) annotation(
      Line(points = {{10, 0}, {50, 0}, {50, 0}, {50, 0}}, color = {0, 127, 255}));
    annotation(
      Documentation(info = "<html><head></head><body><span style=\"font-family: 'Bitstream Vera Sans Mono'; font-size: 12px;\">If used with FreeFluids.TMedia, run with the old frontend. This is done by checking the box at \"Simulation Setup/Translation flags/Enable old frontend for code generation\"</span></body></html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 1),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
  end HEXgeneric2PhSteamCondensing;
end Examples;

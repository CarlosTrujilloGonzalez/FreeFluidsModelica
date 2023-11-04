within FreeFluids.ExternalPure;

package BuildingsTests
model Carnot_TCon
  "Test model for heat pump based on Carnot efficiency and condenser outlet temperature control signal"
  extends Modelica.Icons.Example;
  replaceable package Medium1 = Buildings.Media.Water "Medium model";
  replaceable package Medium2 = Buildings.Media.Water "Medium model";
  parameter Modelica.Units.SI.TemperatureDifference dTEva_nominal=-5
    "Temperature difference evaporator inlet-outlet";
  parameter Modelica.Units.SI.TemperatureDifference dTCon_nominal=10
    "Temperature difference condenser outlet-inlet";
  parameter Modelica.Units.SI.HeatFlowRate QCon_flow_nominal=100E3
    "Evaporator heat flow rate";
  parameter Modelica.Units.SI.MassFlowRate m1_flow_nominal=QCon_flow_nominal/
      dTCon_nominal/4200 "Nominal mass flow rate at condenser";

  Buildings.Fluid.HeatPumps.Carnot_TCon heaPum(
    redeclare package Medium1 = Medium1,
    redeclare package Medium2 = Medium2,
    dTEva_nominal=dTEva_nominal,
    dTCon_nominal=dTCon_nominal,
    m1_flow_nominal=m1_flow_nominal,
    show_T=true,
    allowFlowReversal1=false,
    allowFlowReversal2=false,
    use_eta_Carnot_nominal=true,
    QCon_flow_nominal=QCon_flow_nominal,
    dp1_nominal=6000,
    dp2_nominal=6000) "Heat pump"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Buildings.Fluid.Sources.MassFlowSource_T sou1(nPorts=1,
    redeclare package Medium = Medium1,
    m_flow=m1_flow_nominal,
    T=293.15)
    annotation (Placement(transformation(extent={{-60,-4},{-40,16}})));
  Buildings.Fluid.Sources.MassFlowSource_T sou2(nPorts=1,
    redeclare package Medium = Medium2,
    use_T_in=false,
    use_m_flow_in=true,
    T=288.15)
    annotation (Placement(transformation(extent={{60,-16},{40,4}})));
  Modelica.Blocks.Sources.Ramp TConLvg(
    duration=60,
    startTime=1800,
    height=15,
    offset=273.15 + 35) "Control signal for condenser leaving temperature"
    annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
  Modelica.Blocks.Math.Gain mEva_flow(k=-1/cp2_default/dTEva_nominal)
    "Evaporator mass flow rate"
    annotation (Placement(transformation(extent={{34,-88},{54,-68}})));
  Modelica.Blocks.Math.Add QEva_flow(k2=-1) "Evaporator heat flow rate"
    annotation (Placement(transformation(extent={{32,-48},{52,-28}})));
  Buildings.Fluid.Sources.Boundary_pT sin2(
    redeclare package Medium = Medium2,
    nPorts=1)
    annotation (Placement(transformation(extent={{-60,-40},{-40,-20}})));
  Buildings.Fluid.Sources.Boundary_pT sin1(
    redeclare package Medium = Medium1,
    nPorts=1)
    annotation (Placement(transformation(extent={{60,28},{40,48}})));

  final parameter Modelica.Units.SI.SpecificHeatCapacity cp2_default=
      Medium2.specificHeatCapacityCp(Medium2.setState_pTX(
      Medium2.p_default,
      Medium2.T_default,
      Medium2.X_default))
    "Specific heat capacity of medium 2 at default medium state";

equation
  connect(sou1.ports[1], heaPum.port_a1) annotation (Line(
      points={{-40,6},{-10,6}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sou2.ports[1], heaPum.port_a2) annotation (Line(
      points={{40,-6},{10,-6}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(QEva_flow.y,mEva_flow. u) annotation (Line(points={{53,-38},{64,-38},
          {64,-58},{24,-58},{24,-78},{32,-78}}, color={0,0,127}));
  connect(TConLvg.y, heaPum.TSet) annotation (Line(points={{-59,50},{-20,50},{
          -20,9},{-12,9}}, color={0,0,127}));
  connect(mEva_flow.y, sou2.m_flow_in) annotation (Line(points={{55,-78},{74,
          -78},{74,-10},{74,2},{62,2}}, color={0,0,127}));
  connect(QEva_flow.u1, heaPum.QCon_flow) annotation (Line(points={{30,-32},{20,
          -32},{20,9},{11,9}}, color={0,0,127}));
  connect(QEva_flow.u2, heaPum.P) annotation (Line(points={{30,-44},{16,-44},{16,
          0},{11,0}},    color={0,0,127}));
  connect(sin2.ports[1], heaPum.port_b2) annotation (Line(points={{-40,-30},{-20,
          -30},{-20,-6},{-10,-6}}, color={0,127,255}));
  connect(heaPum.port_b1, sin1.ports[1]) annotation (Line(points={{10,6},{30,6},
          {30,38},{40,38}}, color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatPumps/Examples/Carnot_TCon.mos"
        "Simulate and plot"),
    Documentation(
info= "<html><head></head><body><p>Example copied from the Buildings library 10.0.0</p><p>
Example that simulates a chiller whose efficiency is scaled based on the
Carnot cycle.
The chiller takes as an input the evaporator leaving water temperature.
The condenser mass flow rate is computed in such a way that it has
a temperature difference equal to <code>dTEva_nominal</code>.
</p>
</body></html>",
revisions="<html>
<ul>
<li>
February 10, 2023, by Michael Wetter:<br/>
Removed binding of parameter with same value as the default.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1692\">#1692</a>.
</li>
<li>
May 2, 2019, by Jianjun Hu:<br/>
Replaced fluid source. This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1072\"> #1072</a>.
</li>
<li>
November 25, 2015, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
  __OpenModelica_simulationFlags(jacobian = "numerical", lv = "LOG_STATS", s = "dassl", variableFilter = ".*"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=nonewInst");
end Carnot_TCon;

  model Carnot_TCon_TMedia
  extends Carnot_TCon(redeclare package Medium1 = FreeFluids.TMedia.Fluids.Water, redeclare package Medium2 = FreeFluids.TMedia.Fluids.Water);
    annotation(
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=nonewInst");
  end Carnot_TCon_TMedia;
  
  model ReciprocatingCompressor
  extends Modelica.Icons.Example;
  Buildings.Fluid.HeatPumps.Compressors.ReciprocatingCompressor com(
    redeclare replaceable package ref = Buildings.Media.Refrigerants.R410A,
    pisDis=0.00162,
    cleFac=0.069,
    etaEle=0.696,
    PLos=100,
    dTSup=9.82,
    pDro=99290) "Reciprocating compressor"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature eva(T=253.15)
    "Evaporating temprature"
    annotation (Placement(transformation(extent={{-82,-10},{-62,10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature con(T=323.15)
    "Condensing temperature"
    annotation (Placement(transformation(extent={{84,-10},{64,10}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heaFloEva
    "Evaporator heat flow rate sensor"
    annotation (Placement(transformation(extent={{-46,-10},{-26,10}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heaFloCon
    "Condenser heat flow rate sensor"
    annotation (Placement(transformation(extent={{28,-10},{48,10}})));
  Modelica.Blocks.Sources.Constant on(k=1) "Compressor control signal"
    annotation (Placement(transformation(extent={{-26,18},{-6,38}})));
equation
  connect(eva.port, heaFloEva.port_a)
    annotation (Line(points={{-62,0},{-46,0}}, color={191,0,0}));
  connect(heaFloEva.port_b, com.port_a)
    annotation (Line(points={{-26,0},{-10,0}}, color={191,0,0}));
  connect(com.port_b, heaFloCon.port_a)
    annotation (Line(points={{10,0},{19,0},{28,0}}, color={191,0,0}));
  connect(heaFloCon.port_b, con.port)
    annotation (Line(points={{48,0},{56,0},{64,0}}, color={191,0,0}));
  connect(on.y,com.y)
    annotation (Line(points={{-5,28},{6,28},{6,11}}, color={0,0,127}));
  annotation (__Dymola_Commands(file= "modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatPumps/Compressors/Validation/ReciprocatingCompressor.mos"
        "Simulate and plot"),
    experiment(
      Tolerance=1e-6, StopTime=100),
    Documentation(info= "<html><head></head><body><p>Copied from Buildings.Fluid.HeatPumps.Compressors.Validation.ReciprocatingCompressor</p><p>
Model that demonstrates the use of the ReciprocatingCompressor model.
</p>
<p>
The compressor power, condenser heat transfer rate and evaporator heat transfer rate are calculated for given refrigerant temperatures.
</p>
</body></html>", revisions="<html>
<ul>
<li>
October 17, 2016, by Massimo Cimmino:<br/>
First implementation.
</li>
</ul>
</html>"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end ReciprocatingCompressor;

  model ReciprocatingCompressor_ExternalPure
  extends ReciprocatingCompressor(com(redeclare package ref = FreeFluids.ExternalPure.Fluids.R410A(thermoModel=3, refState=2, ancillaries=3)));
    annotation(
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  Documentation);
  end ReciprocatingCompressor_ExternalPure;
  
  model ReciprocatingCompressor_Propane
  extends ReciprocatingCompressor(com(redeclare package ref = FreeFluids.ExternalPure.Fluids.Propane(thermoModel=3, refState=2, ancillaries=3)));
    annotation(
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  Documentation);
  end ReciprocatingCompressor_Propane;

  model RefrigerantsTest
    package Medium1 = Buildings.Media.Refrigerants.R410A;
    package Medium2 = FreeFluids.ExternalPure.Fluids.R410A;
    parameter Real v=1/10;
    parameter Real T=250;
    Real dP_dV1, dP_dV2;
    Real dP_dT1, dP_dT2;
  algorithm
    dP_dV1:=Medium1.dPressureVap_dSpecificVolume_Tv(T,v);
    dP_dV2:=Medium2.dPressureVap_dSpecificVolume_Tv(T,v);
    dP_dT1:=Medium1.dPressureVap_dTemperature_Tv(T,v);
    dP_dT2:=Medium2.dPressureVap_dTemperature_Tv(T,v);
  end RefrigerantsTest;
  annotation(
    Documentation(info = "<html><head></head><body>Needs Buldings library 10.0.0</body></html>"));
end BuildingsTests;

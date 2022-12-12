within FreeFluids.HeatExchangers;

package BaseClasses
  partial model HEXtubesDefinition
    replaceable package MediumI = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "inner pipe medium";
    replaceable package MediumO = FreeFluids.TMedia.Fluids.Water constrainedby Modelica.Media.Interfaces.PartialMedium "outer flow medium";
    FreeFluids.Interfaces.FluidPortA Iin annotation(
      Placement(visible = true, transformation(origin = {-92, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB Iout annotation(
      Placement(visible = true, transformation(origin = {92, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortA Oin annotation(
      Placement(visible = true, transformation(origin = {70, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    FreeFluids.Interfaces.FluidPortB Oout annotation(
      Placement(visible = true, transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter Modelica.Units.SI.Distance iElevDiff=0 "inner pipe differential elevation: out elev. - in elev." annotation(
      Dialog(tab = "Flow data"));
    parameter Modelica.Units.SI.Distance oElevDiff=0 "outer flow differential elevation: out elev. - in elev." annotation(
      Dialog(tab = "Flow data"));
    parameter Boolean iIsCircular=true "if false, the inner pipe is not circular" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iDi(displayUnit = "mm") = 0.0 "inner pipe internal diameter, if circular" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Area iSection = 0 "inner pipe flow section, if not circular" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iPerimeter = 0 "inner pipe perimeter, if not circular" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iThickness(displayUnit = "mm") = 1e-3 "inner pipe thickness" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iRoughness(displayUnit = "mm") = 1.5e-005 "inner pipe roughness. SS:1.5e-5, Steel new:4.6e-5, Steel old:2.0e-4, Concrete:1.5e-3" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Density iRhoWall(displayUnit = "kg/m3") = 8000 "inner pipe wall density" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.ThermalConductivity iKwall = 16 "inner pipe wall thermal conductivity. Al=210, Cu=390, Steel=50, SS=17" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.ThermalInsulance iFoulingF = 0.0002 "inner pipe side fouling factor.Typical: 0.00018 for thermal oil, or treated cooling water" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Distance iLTube = 0 "stright length of each individual inner pipe" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Integer iNumPipes = 1 "number of inner pipes that form a common pass" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Integer iNumParallel = 1 "number of parallel inner passes. Normally it is the number of exchangers" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Integer iNumSerial=2 "number of serial inner passes  (Each hairpin has two serial passes)" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Integer iNumUturns=1 "number of U turns in each exchanger" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Real iNumVelocityHeads=1.5 "inner pipe number of velocity heads to consider in pressure loss" annotation(
      Dialog(tab = "Flow data"));
    parameter Modelica.Units.SI.ThermalInsulance oFoulingF = 0.0002 "outer flow fouling factor.Typical: 0.00018 for thermal oil, or treated cooling water" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Boolean useHTcorr=true "if true, the heat transfer correction factor is used"  annotation(
      Dialog(tab = "Heat transfer"));
              
    parameter Boolean useFins=false "if true, fins are taken into account in calculations" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Density finRho = 0 "fin density. Al=2700, Cu=8930, Steel=7830, SS=8000" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.ThermalConductivity finK = 0 "fin thermal conductivity. Al=210, Cu=390, Steel=50, SS=17" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finThickness(displayUnit = "mm") = 0.0 "fin thickness" annotation(
      Dialog(tab = "Fins data"));
  
    parameter Boolean counterCurrent=true "true if flow is counter-current, otherwise it will be considered co-current" annotation(
      Dialog(tab = "Heat transfer"));
  
    Integer ItotalNumPipes "total number of individual inner pipes";
    Modelica.Units.SI.Distance IiPerimeter "inner pipe internal wetted perimeter"; 
    Modelica.Units.SI.Area IiSection "inner pipe internal flow section";
    Modelica.Units.SI.Distance IDh "internal pipe hydraulic diameter";
    Modelica.Units.SI.Area IiArea "total inner pipes internal heat transfer area";
    Modelica.Units.SI.Distance IoPerimeter "inner pipe external perimeter";  
    Modelica.Units.SI.Area IoSection "inner pipe external section";
    Modelica.Units.SI.Area IoArea "total inner pipes external heat transfer area";
    MediumI.ThermodynamicState IStateIn "inner flow inlet state";
    MediumI.ThermodynamicState IStateOut "inner flow outlet state";
    Modelica.Units.SI.Temperature ITin "inner flow inlet temperature";
    Modelica.Units.SI.Temperature ITout "iiner flow outlet temperature";
    Modelica.Units.SI.Temperature IT "inner pipe average temperature";  
    MediumI.ThermodynamicState IStateWall "inner pipe wall state for correction factors
     calculation";
    Modelica.Units.SI.Temperature ITwall "inner pipe wall temperature";
    Modelica.Units.SI.DynamicViscosity IMuWall(min = 1e-6, start = 1e-3, max = 1e6) "inner pipe wall dynamic viscosity";
    Real IPlossCorr "internal pipe pressure loss correction factor";
    Modelica.Units.SI.AbsolutePressure IPloss "internal pipe pressure loss by friction";
    Real IHTcorr "internal pipe heat transfer correction factor";
    Modelica.Units.SI.NusseltNumber INu "inner pipe Nusselt number";
    Modelica.Units.SI.CoefficientOfHeatTransfer IH(min = 1, start = 1000) "inner pipe average heat transfer coefficient"; 
      
  
    MediumO.ThermodynamicState OStateIn "outer flow inlet state";
    MediumO.ThermodynamicState OStateOut "outer flow outlet state";
    Modelica.Units.SI.Temperature OTin "outer flow inlet temperature";
    Modelica.Units.SI.Temperature OTout "outer flow outlet temperature";
    Modelica.Units.SI.Temperature OT "outer flow average temperature";  
    Real OPlossCorr "outer flow pressure loss correction factor";
    Modelica.Units.SI.AbsolutePressure OPloss "outer flow pressure loss by friction";
    Modelica.Units.SI.Area OiArea "total outer flow internal area for heat transfer";
    Fraction OiEfficiency "outer flow internal area efficiency";
    Modelica.Units.SI.NusseltNumber ONu "outer flow nusselt number for internal heat transfer";
    Modelica.Units.SI.CoefficientOfHeatTransfer OHi(min = 1, start = 1000) "outer flow average heat transfer coefficient to inner pipe";
  
    Modelica.Units.SI.ThermalResistance TRii "total resitance for the transfer between inner fluid and inner wall";
    Modelica.Units.SI.ThermalResistance TRio "total resitance for the transfer between outer fluid and inner wall";
    Modelica.Units.SI.ThermalResistance TRi "total resitance for the transfer between inner and outer fluids";
    Modelica.Units.SI.CoefficientOfHeatTransfer Ui "global heat transfer coefficient referenced to total inner flow area"; 
    Modelica.Units.SI.TemperatureDifference ILMTD "LMTD between inner and outer fluids";
  
    Modelica.Units.SI.Power IW "heat exchanged between inner and outer fluids. Positive if heat enters the inner fluid" ;
    
    Modelica.Units.SI.Length IEquivLength "internal pipe equivalent length for pressure drop";
    Fraction ILMTDcorr(min=0.5) "internal LMTD correction factor";
    Real Rntu(min=0.0);
    Real NTU(min=0.0);
  annotation(
      Documentation(info = "<html><head></head><body>Contains just the definitions of the inner pipe elements needed for the development of the outer flow models. With no equations, in order to allow multiple inheritance of the model without giving too many equations.</body></html>")); 
  end HEXtubesDefinition;

  partial model HEXtubes
    extends HEXtubesDefinition; 
  
  algorithm
    if iIsCircular==true then
      IiPerimeter:=pi*iDi;
      IiSection:=pi*iDi^2/4;
      IDh:=iDi;
      IoPerimeter:=pi*(iDi+2*iThickness);
      IoSection:=pi*(iDi+2*iThickness)^2/4;
    else
      IiPerimeter:=iPerimeter;
      IiSection:=iSection;
      IDh:=4*iSection/iPerimeter;
      IoPerimeter:=iPerimeter+8*iThickness;
      IoSection:=iSection+(iPerimeter+4*iThickness)*iThickness;
    end if;
    ItotalNumPipes:=iNumParallel*iNumSerial*iNumPipes;
    IEquivLength:=iNumSerial*iLTube;
    IiArea:=ItotalNumPipes*iLTube*IiPerimeter; 
  equation
    Iout.X=Iin.X;
    Iout.Elevation=Iin.Elevation+iElevDiff;
    Iout.G=-Iin.G;
    IStateIn = MediumI.setState_phX(Iin.P, Iin.H, Iin.X);
    IStateOut = MediumI.setState_phX(Iout.P, Iout.H, Iout.X); 
    ITin = MediumI.temperature(IStateIn);
    ITout = MediumI.temperature(IStateOut);
    ITwall=(TRio*IT+TRii*OT)/(TRio+TRii);
//ITwall=(IH*IiArea*IT+OHi*OiArea*OiEfficiency*OT)/(IH*IiArea+OHi*OiArea*OiEfficiency);
//ITwall =(ITin+ITout+OTin+OTout)/4;
    IStateWall = MediumI.setState_pTX((Iin.P + Iout.P) / 2, ITwall, Iin.X); 
    IMuWall = MediumI.dynamicViscosity(IStateWall);
  
    Oout.X=Oin.X;
    Oout.Elevation=Oin.Elevation+oElevDiff;
    Oout.G=-Oin.G;
    OStateIn = MediumO.setState_phX(Oin.P, Oin.H, Oin.X);
    OStateOut = MediumO.setState_phX(Oout.P, Oout.H, Oout.X);    
    OTin = MediumO.temperature(OStateIn);
    OTout = MediumO.temperature(OStateOut);
    
    if counterCurrent==true then  
      if noEvent((ITin - OTout) * (ITout - OTin) <= 0) then
        ILMTD = 0;
      elseif noEvent(abs(ITin - OTout - ITout + OTin)<0.001) then
        ILMTD = ITin - OTout;
      else
        ILMTD = (ITin - OTout - ITout + OTin) / log((ITin - OTout) / (ITout - OTin));
      end if;
    else
      if noEvent((ITin - OTin) * (ITout - OTout) <= 0) then
        ILMTD = 0;
      elseif noEvent(abs(ITin - OTin - ITout + OTout)<0.001) then
        ILMTD = ITin - OTout;
      else
        ILMTD = (ITin - OTin - ITout + OTout) / log((ITin - OTin) / (ITout - OTout));
      end if;  
    end if;
    TRii=1/(IH*IiArea)+iFoulingF/IiArea;
    TRio=oFoulingF/(OiArea*OiEfficiency)+1/(OHi*OiArea*OiEfficiency);
    TRi=1/(IH*IiArea)+iFoulingF/IiArea+log((IDh+2*iThickness)/IDh)/(iKwall*2*pi*ItotalNumPipes*iLTube)+oFoulingF/(OiArea*OiEfficiency)+1/(OHi*OiArea*OiEfficiency);
    IW=-ILMTD*ILMTDcorr/TRi;
    Iout.H=Iin.H+IW/Iin.G "kinetic and gravitational energy are not taken into account";
    Oout.H=Oin.H-IW/Oin.G "kinetic and gravitational energy are not taken into account";
    Ui=1/(TRi*IiArea);
  annotation(
      Documentation(info = "<html><head></head><body>Here come the equations missing in the HEXtubesDefinition model. The model has the elements and equations common to forced convection and gas condensation.<div>It is not clear how to calculate the average wall temperature used to introduce corrections in the calculated heat transfer coefficient. I see three possibilities:</div><div>· The average temperature of both streams at both ends.</div><div>· Split the temperature difference (between the two streams average temperatures) according to the respective heat trasfer resistance of each flow. It seems better, as it is clear that the wall temperature will be closer to the stream with less heat transfer resistance.</div><div>· Calculate a wall temperature that will produce the same heat transfer considering the LMTD between this temperature and one of the streams. This is a more complex calculation, and I have not checked that the obtained temperature is the same for both streams.</div><div>For now, split the temperature according to the resistance seems the best option.</div></body></html>"));end HEXtubes;

  partial model HEXtubesForcedConvection
    extends HEXtubes;
    parameter Boolean useTwistedTapeInserts = false "if true, twisted tape inserts will be used in inner pipes" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Length tapeWidth(displayUnit="mm") = iDi "only if twisted tape inserts are used" annotation(
      Dialog(tab = "Heat transfer"));
    parameter Modelica.Units.SI.Length tapeThickness(displayUnit="mm") = iThickness "only if twisted tape inserts are used" annotation(
      Dialog(tab = "Heat transfer"));    
    parameter Real twistRatio = 6 "Length of a 180 degrees twist divided by Di, if twisted tape inserts are used" annotation(
      Dialog(tab = "Heat transfer"));  
    MediumI.ThermodynamicState IStateAvg "inner pipe average state for physical properties calculation";
    Modelica.Units.SI.Density IRho(displayUnit = "kg/m3") "inner pipe average density";
    Modelica.Units.SI.DynamicViscosity IMu(min = 1e-6, start = 1e-3, max = 1e6) "inner pipe average dynamic viscosity";
    Modelica.Units.SI.SpecificHeatCapacity ICp(start = 2000.0) "inner pipe average heat capacity";
    Modelica.Units.SI.ThermalConductivity IK(start = 0.1) "inner pipe average thermal conductivity";
    Modelica.Units.SI.PrandtlNumber IPr "inner pipe Prandt number";
    Modelica.Units.SI.VolumeFlowRate IQ(displayUnit = "m3/h") "inner pipe total volume flow rate at average conditions";
    Modelica.Units.SI.Velocity IV(start = 1) "velocity at average conditions referenced to the empty pipe. Normally between 0.9 and 3.0 m/s for liquids";
    Modelica.Units.SI.ReynoldsNumber IRe(min = 0.01, start = 20000) "inner pipe Reynolds number at average conditions";
    Real IFl(start = 0.01) "internal flowlaminar Darcy's friction factor at average conditions";
    Real IFt(start = 0.01) "internal flow turbulent Darcy's friction factor at average conditions";
    Real IF(start = 0.01) "internal flow Darcy's friction factor at average conditions";
    Modelica.Units.SI.RayleighNumber IRa " inner flow Rayleigh number";
    Real ISw "inner flowSwirl number";
    Real IfSmooth;
  
  equation
    IStateAvg = MediumI.setState_phX((Iin.P + Iout.P) / 2, (Iin.H + Iout.H) / 2, Iin.X);
    IT = MediumI.temperature(IStateAvg);   
    IRho = abs(MediumI.density(IStateAvg));
    IMu = MediumI.dynamicViscosity(IStateAvg);
    ICp = MediumI.specificHeatCapacityCp(IStateAvg);
    IK = MediumI.thermalConductivity(IStateAvg);
    IPr = ICp * IMu / IK;
    IQ = Iin.G / IRho;
    IV = abs(IQ / (iNumParallel*iNumPipes*IiSection));
    if noEvent(Iin.G < (-1e-7) or Iin.G > 1e-7) then
      IRe = IDh * abs(Iin.G) / (iNumParallel*iNumPipes*IiSection * IMu) + 0.01 "in order to avoid division by 0 calculating F";
      
    else
      IRe = 0.01;
  
    end if;
  if useTwistedTapeInserts == true then
      ISw = IRe * pi / (twistRatio ^ 0.5 * (pi - 4 * tapeThickness / IDh)) * (1 + (pi / (2 * twistRatio)) ^ 2) ^ 0.5;
//IRa=g_n*IRho^2*IDh^3*MediumI.isobaricExpansionCoefficient(IStateAvg)*abs(Twall-Tsurf)*IPr/IMu^2;
      IRa = 0;
      IFl = 4 * 15.767 / IRe * (pi / (pi - 4 * tapeThickness / IDh)) * ((pi + 2 - 2 * tapeThickness / IDh) / (pi - 4 * tapeThickness / IDh)) ^ 2 * (1 + (pi / 2 / twistRatio) ^ 2) * (1 + 1e-6 * ISw ^ 2.55) ^ (1 / 6) "laminar flow friction factor";
      IFt = 4 * 0.0791 / IRe ^ 0.25 * (pi / (pi - 4 * tapeThickness / IDh)) * ((pi + 2 - 2 * tapeThickness / IDh) / (pi - 4 * tapeThickness / IDh)) ^ 1.25 * (1 + 2.752 / twistRatio ^ 1.29) "turbulent flow friction factor";
      IF = (IFl ^ 10 + IFt ^ 10) ^ 0.1;
    else
      ISw = 0;
      IRa = 0;
      IFl = 0;
      IFt = 0;
      IF = 8 * ((8 / IRe) ^ 12 + ((37530 / IRe) ^ 16 + (-2.457 * log((7 / IRe) ^ 0.9 + 0.27 * iRoughness / IDh)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
    end if;
    
    if IRho>= 300 then
      IPlossCorr=(IMuWall/IMu)^0.25 "friction correction factor for liquids. S.Kalkaç pag. 89 turbulent:0.25, laminar: 0.5. Serth: 0.14 for turbulent and 0.25 for laminar. D.G.Kroger turbulent:0.24 laminar: 0.56";
    else
      IPlossCorr=(IT/ITwall)^0.3 "friction correction factor for gases. Kalkaç: turbulent 0.3, laminar -0.9. D.G.Kroger turbulent:0.x laminar -0.9";
    end if;
    IPloss = homotopy(0.5 * IV ^ 2 * IRho * IPlossCorr * (IF * IEquivLength / IDh + iNumVelocityHeads), 0.5 * IV ^ 2 * IRho * (IF * IEquivLength / IDh + iNumVelocityHeads));
    Iout.P-Iin.P = (-sign(Iin.G) * IPloss) + (Iin.Elevation - Iout.Elevation + 1e-5) * g_n * IRho "momentum change is not taken into account.1 e-5 is to avoid division by 0";
    
     if useHTcorr == true then
      if useTwistedTapeInserts==true then
        if IRho>300 then
          if noEvent(ITwall>IT) then
            IHTcorr = homotopy((IMu / IMuWall) ^ 0.18,1) "liquid heating heat transfer correction factor with inserts";
          else
            IHTcorr = homotopy((IMu / IMuWall) ^ 0.3,1) "liquid cooling heat transfer correction factor with inserts";
          end if;
        else
          if noEvent(ITwall>IT) then
            IHTcorr = homotopy((IT / ITwall) ^ 0.45,1) "gas heating heat transfer correction factor with inserts";
          else
            IHTcorr = homotopy((IT / ITwall) ^ 0.15,1) "gas cooling heat transfer correction factor with inserts";
          end if;
        end if;
      else
        if IRho>300 then
          IHTcorr = homotopy((IMu / IMuWall) ^ 0.11,1) "liquids heat transfer correction factor without inserts";
        else
          IHTcorr = 1.0 "gases without inserts";
        end if;
      end if;
    else
      IHTcorr = 1;
    end if;
  if noEvent(IRe > 10000) then
      IfSmooth = (0.78173 * log(IRe) - 1.5) ^ (-2);
      INu = OpenModelica.Internal.realAbs(IfSmooth / 8 * IRe * IPr / (1 + 12.7 * (IfSmooth / 8) ^ 0.5 * (IPr ^ (2 / 3) - 1)) * (1 + (IDh / iLTube) ^ (2 / 3)) * IHTcorr) "Pethukov/Gnielinsky equation for smooth tubes, VDI mean";
      //INu=0.023*IRe^0.8*IPr^(1/3);
    elseif noEvent(IRe < 2100) then
      IfSmooth = 0;
      if noEvent(iLTube > 0.05 * IRe * IPr * IDh) then
        INu = abs((3.66 ^ 3 + 0.7 ^ 3 + (1.615 * (IRe * IPr * IDh / iLTube) ^ (1 / 3) - 0.7) ^ 3) ^ (1 / 3)) * IHTcorr "VDI Atlas";
//INu = abs((3.66 ^ 3 + 0.7 ^ 3 + (1.65 * (IRe * IPr * IDh / iLTube) ^ 0.333 - 0.7) ^ 3 + ((2 / (1 + 22 * IPr)) ^ (1 / 6) * (IRe * IPr * IDh / iLTube) ^ 0.5) ^ 3) ^ 0.333 * IHTcorr) "Gnielinsky-Martin correlation for constant wall temperature: VDI mean";
      else
        INu = abs(1.615 * (IRe * IPr * IDh / iLTube) ^ (1 / 3)) * IHTcorr;
//INu = abs(1.86 * (IRe * IPr * IDh / iLtube) ^ 0.333 * IHTcorr) "Sieder-Tate equation for laminar flow";
      end if;
    else
//interpolation between turbulent and laminar flow
      IfSmooth = (0.78173 * log(10000) - 1.5) ^ (-2);
      INu = (abs(IfSmooth / 8 * 10000 * IPr / (1 + 12.7 * (IfSmooth / 8) ^ 0.5 * (IPr ^ (2 / 3) - 1)) * (1 + (IDh / iLTube) ^ (2 / 3))) * (IRe - 2100) / 7900 + abs((3.66 ^ 3 + 0.7 ^ 3 + (1.65 * (2100 * IPr * IDh / iLTube) ^ (1 / 3) - 0.7) ^ 3) ^ (1 / 3)) * (10000 - IRe) / 7900) * IHTcorr "VDI G1 4.2";
    end if;
    IH=INu*IK/IDh;
    
    Rntu = (OTout - OTin) / (ITin - ITout);
    NTU=(ITout-ITin)/(TRi*IW);    
  annotation(
      Documentation(info = "<html><head></head><body>This partial model performs the calculation of the heat transfer coefficient, and the pressure loss, of tubes side of an exchanger.<div>Twistted tape inserts can be used optionally.</div></body></html>"));end HEXtubesForcedConvection;

  partial model DoublePipeHEX
  "double pipe heat exchanger with internal and external forced convection transfer"
    extends HEXtubesDefinition;
    parameter Boolean oIsCircular=true "if false, the outer pipe is not circular" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Distance oDi(displayUnit = "mm") = 0.0 "outer pipe internal diameter, if circular" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Area oSection = 0 "outer pipe raw section, if not circular. The section of the inner pipe will be substracted from it" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Distance oPerimeter = 0 "outer pipe internal perimeter, if not circular" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Distance oThickness(displayUnit = "mm") = 1e-3 "outer pipe thickness" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Distance oRoughness(displayUnit = "mm") = 1.5e-005 "outer pipe roughness. SS:1.5e-5, Steel new:4.6e-5, Steel old:2.0e-4, Concrete:1.5e-3" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Density oRhoWall(displayUnit = "kg/m3") = 8000 "outer pipe wall density" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.ThermalConductivity oKwall = 16 "outer pipe wall thermal conductivity. Al=210, Cu=390, Steel=50, SS=17" annotation(
      Dialog(tab = "Outer pipe data"));
    parameter Modelica.Units.SI.Distance thicknessInsul(displayUnit = "mm") = 0 "insulation thickness" annotation(
      Dialog(tab = "Outer pipe data"));    
    parameter Real oNumVelocityHeads=0 "outer flow number of velocity heads to consider in pressure loss" annotation(
      Dialog(tab = "Flow data"));    
    
    parameter Boolean oSerialFlow=false "if false, annulus flow is splitted in the same number of parallel flows than the inner flow" annotation(
      Dialog(tab = "Flow data"));  
  
    parameter Integer finNum = 1 "number of longitudinal fins on each individual inner pipe" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finLength = iLTube "length of each longitudinal fin" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finHeight(displayUnit = "mm") = 0.0 "fin height" annotation(
      Dialog(tab = "Fins data"));
  
    Modelica.Units.SI.Distance FinsPerimeter "wetted fins perimeter on each internal pipe";
    Modelica.Units.SI.Area FinsSection "fins section on each internal pipe";
    Modelica.Units.SI.Area FinsArea "total fins area for heat transfer";
  
    Modelica.Units.SI.Distance OoPerimeter "outer pipe external perimeter";  
    Modelica.Units.SI.Area OoSection "Outer pipe external section";
    Modelica.Units.SI.Distance OiPerimeter "Outer pipe internal wetted perimeter"; 
    Modelica.Units.SI.Area OiSection "Outer pipe total internal section";
    
    Modelica.Units.SI.Distance Operimeter "annulus wetted perimeter";
    Modelica.Units.SI.Area Osection "annulus nett flow section";  
    Modelica.Units.SI.Distance ODh "annulus hydraulic diameter";
    
    MediumO.ThermodynamicState OStateAvg "outer flow average state for physical properties calculation";
    MediumO.ThermodynamicState OiStateWall "annulus internal wall state for correction factors calculation";
    Modelica.Units.SI.Density ORho(displayUnit = "kg/m3") "annulus average density";
    Modelica.Units.SI.DynamicViscosity OMu(min = 1e-6, start = 1e-3, max = 1e6) "annulus average dynamic viscosity";
    Modelica.Units.SI.SpecificHeatCapacity OCp(start = 2000.0) "annulus average heat capacity";
    Modelica.Units.SI.ThermalConductivity OK(start = 0.1) "annulus average thermal conductivity";  
    Modelica.Units.SI.DynamicViscosity OiMuWall(min = 1e-6, start = 1e-3, max = 1e6) "annulus inner wall dynamic viscosity"; 
    Modelica.Units.SI.PrandtlNumber OPr "annulus Prandt number";
    Modelica.Units.SI.VolumeFlowRate OQ(displayUnit = "m3/h") "total annulus volume flow rate at average conditions";
    Modelica.Units.SI.ReynoldsNumber ORe(min = 0.01, start = 20000) "annulus average Reynolds number";
    Real OF(start = 0.01) "annulus Darcy's friction factor at average conditions";
  
    Real OfSmooth;
    Modelica.Units.SI.CoefficientOfHeatTransfer OHo(min = 1, start = 1000) "annulus average heat transfer coefficient to outer pipe";
    Modelica.Units.SI.Diameter ODei "annulus diameter for heat transfer to internal pipe";
    Modelica.Units.SI.Diameter ODeo "annulus diameter for heat transfer to external surface";
  
    Modelica.Units.SI.Velocity OV(start = 1) "annulus flow velocity at average conditions. Normally between 0.9 and 3.0 m/s for liquids";
    Real FinM;
  
    Modelica.Units.SI.Length OEquivLength "external pipe equivalent length for pressure drop";
  
    Real OiHTcorr "external flow internal heat trasfer correction factor. 0.873";
  
  algorithm
    if useFins==true then
      IoArea:=ItotalNumPipes*iLTube*(IoPerimeter-finNum*finThickness);
    else
      IoArea:=ItotalNumPipes*iLTube*IoPerimeter;      
    end if;
  
  
    FinsPerimeter:=finNum*(2*finHeight+finThickness);
    FinsSection:=finNum*finHeight*finThickness;
    FinsArea:=ItotalNumPipes*finLength*FinsPerimeter;
    
    if oIsCircular==true then
      OoPerimeter:=pi*(oDi+2*oThickness);
      OoSection:=pi*(oDi+2*oThickness)^2/4;  
      OiPerimeter:=pi*oDi;
      OiSection:=pi*oDi^2/4;
    else
      OoPerimeter:=oPerimeter+8*oThickness;
      OoSection:=oSection+(oPerimeter+4*oThickness)*oThickness;
      OiPerimeter:=oPerimeter;
      OiSection:=oSection;
    end if;
    if oSerialFlow==false then
      OEquivLength:=iNumSerial*iLTube;  
    else
      OEquivLength:=iNumParallel*iNumSerial*iLTube;
    end if;
    
    if useFins==true then
      Operimeter:=iNumPipes*(IoPerimeter-finNum*finThickness+FinsPerimeter)+OiPerimeter;
      Osection:=OiSection-iNumPipes*(IoSection+FinsSection);
      ODei:=4*Osection/(iNumPipes*(IoPerimeter-finNum*finThickness+FinsPerimeter));
      OiArea:=IoArea+FinsArea;
    else
      Operimeter:=iNumPipes*IoPerimeter+OiPerimeter;
      Osection:=OiSection-iNumPipes*IoSection;
      ODei:=4*Osection/(iNumPipes*IoPerimeter);
      OiArea:=IoArea;
    end if;
    ODh:=4*Osection/Operimeter;
    ODeo:=4*Osection/OiPerimeter; 
  
   
  equation
   
    OStateAvg = MediumO.setState_phX((Oin.P + Oout.P) / 2, (Oin.H + Oout.H) / 2, Oin.X);
    OiStateWall = MediumO.setState_pTX((Oin.P + Oout.P) / 2, ITwall, Oin.X);
    OT = MediumO.temperature(OStateAvg);   
    ORho = abs(MediumO.density(OStateAvg));
    OMu = MediumO.dynamicViscosity(OStateAvg);
    OCp = MediumO.specificHeatCapacityCp(OStateAvg);
    OK = MediumO.thermalConductivity(OStateAvg);
    OiMuWall = MediumO.dynamicViscosity(OiStateWall);
    OPr = OCp * OMu / OK;
    
    OQ = Oin.G / ORho;
    if oSerialFlow==false then
      OV = abs(OQ / (Osection*iNumParallel));  
    else
      OV = abs(OQ / (Osection));
    end if;
    if noEvent(Oin.G < (-1e-7) or Oin.G > 1e-7) then
      ORe = ODh * abs(Oin.G) / (Osection * OMu) + 0.01 "in order to avoid division by 0 calculating F";
      OF = 8 * ((8 / ORe) ^ 12 + ((37530 / ORe) ^ 16 + (-2.457 * log((7 / ORe) ^ 0.9 + 0.27 * oRoughness / ODh)) ^ 16) ^ (-1.5)) ^ (1 / 12) "Churchill equation for Darcy's friction factor";
    else
      ORe = 0.01;
      OF = 6400;
    end if;
    if ORho >= 300 then
      //OPlossCorr = 1.64;
      OPlossCorr = (OiMuWall/OMu)^0.25 "friction correction factor for liquids. In laminar flow is higher(0.5) according to S.Kalkaç pag. 89, 95, 105. Serth recommends 0.14 for turbulent and 0.25 for laminar";
      if useHTcorr then
        OiHTcorr=homotopy((OMu/OiMuWall)^0.11,1) "Kalkaç turbulent: 0.11 for heating, 0.25 for cooling, laminar: 0.14 for heating and cooling. VDI Atlas 0.11" ;
      else
        OiHTcorr=1;
      end if;
    else
      OPlossCorr = (OT/ITwall)^0.3; //(ITwall / OT) ^ m "laminar m=0.9, turbulent m=-0.3 friction correction factor for gases as per Kalkaç";
      OiHTcorr=1 "(OT/ITwall)^n Kalkaç turbulent: heating 0.47, cooling 0.36  ,laminar: gas 0. VDI Atlas very variable";
    end if;
    
    OPloss = homotopy(0.5 * OV ^ 2 * ORho * OPlossCorr * (OF * OEquivLength / ODh + oNumVelocityHeads), 0.5 * OV ^ 2 * ORho * (OF * OEquivLength / ODh + oNumVelocityHeads));
    Oout.P-Oin.P = (-sign(Oin.G) * OPloss) + (Oin.Elevation - Oout.Elevation + 1e-5) * g_n * ORho "momentum change is not taken into account.ForcedConvection1 e-5 is to avoid division by 0";
    
  
  if noEvent(ORe > 10000) then
      OfSmooth = (0.78173 * log(ORe) - 1.5) ^ (-2);
      ONu = OpenModelica.Internal.realAbs(OfSmooth / 8 * ORe * OPr / (1 + 12.7 * (OfSmooth / 8) ^ 0.5 * (OPr ^ (2 / 3) - 1)) * (1 + (ODh / iLTube) ^ (2 / 3)) * OiHTcorr) "Pethukov/Gnielinsky equation for smooth tubes, VDI Atlas mean";
    elseif noEvent(ORe < 2100) then
      OfSmooth = 0;
      if noEvent(iLTube > 0.05 * ORe * OPr * ODh) then
        ONu = abs((3.66 ^ 3 + 0.7 ^ 3 + (1.615 * (ORe * OPr * ODh / iLTube) ^ (1 / 3) - 0.7) ^ 3) ^ (1 / 3)) * OiHTcorr "VDI Atlas";
      else
        ONu = 1.615 * (ORe * OPr * ODh / iLTube) ^ (1 / 3) * OiHTcorr "VDI Atlas mean. Developing laminar flow, constant wall temperature";
  //ONu=1.86 * (ORe * OPr * ODh / iLTube) ^ (1/3)* OiHTcorr "Sieder-Tate. Recommended by Kalkaç. But this is for local Nusselt, not for average";
  //ONu=3.66+1.2*(4*OiSection/((IDh+2*iThickness)*OiPerimeter))^0.8+(0.19*(1+0.14*(4*OiSection/((IDh+2*iThickness)*OiPerimeter))^0.5)*(ORe*OPr*ODh/iLTube)^0.8)/(1+0.117*(ORe*OPr*ODh/iLTube)^0.467) "Gnielinski, Serth, pag.55";
      end if;
    else
  //interpolation between turbulent and laminar flow
      OfSmooth = (0.78173 * log(10000) - 1.5) ^ (-2);
      ONu = (abs(OfSmooth / 8 * 10000 * OPr / (1 + 12.7 * (OfSmooth / 8) ^ 0.5 * (OPr ^ (2 / 3) - 1)) * (1 + (ODh / iLTube) ^ 0.667)) * (ORe - 2100) / 7900 + abs((3.66 ^ 3 + 0.7 ^ 3 + (1.615 * (2100 * OPr * ODh / iLTube) ^ (1 / 3) - 0.7) ^ 3) ^ (1 / 3)) * (10000 - ORe) / 7900) * OiHTcorr "VDI G1 4.2";
    end if;
    OHi=ONu*OK/ODh;
    OHo=ONu*OK/ODeo;
    if useFins==true then
      FinM=(2*OHi/(finThickness*finK))^0.5;
      OiEfficiency=1-(1-tanh(FinM*finHeight)/(FinM*finHeight))*FinsArea/OiArea;
    else
      FinM=0;
      OiEfficiency=1;
    end if;
    
    ILMTDcorr=1;
  annotation(
      Icon(graphics = {Rectangle(origin = {0, 35}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 5}, {80, -15}}), Text(origin = {90, -8}, lineColor = {0, 0, 255}, extent = {{-150, 100}, {146, 42}}, textString = "%name"), Rectangle(origin = {0, -35}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 15}, {80, -5}}), Rectangle(origin = {-20, 14}, lineColor = {0, 85, 255}, extent = {{-60, 6}, {100, -34}}), Line(origin = {0, 60}, points = {{0, -20}, {0, 20}}, color = {32, 102, 241}, thickness = 0.5), Line(origin = {0, -60}, points = {{0, 20}, {0, -20}, {0, -20}}, color = {30, 81, 241}, thickness = 0.5)}),
      Documentation(info = "<html><head></head><body>This partial model performs the calculation of the heat transfer coefficient, and pressure drop, of the annulus of a doble pipe heat exchanger. Longitudinal fins can be optionally used.</body></html>"));       
  end DoublePipeHEX;
  
  partial model GasCooledHEX
    extends HEXtubesDefinition(redeclare replaceable package MediumO=Modelica.Media.Air.DryAirNasa, final counterCurrent=true);
    parameter Integer iNumRows(start = 1) "number of pipe rows passed by the gas flow" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Boolean iStaggered = true "if false, inline distribution" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iTubePitch(displayUnit = "mm")=2*iDi "distance between tubes" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance iDop(displayUnit = "mm")=iDi+2*iThickness "inner pipe outer diameter, or width, perpenticular to gas flow" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Modelica.Units.SI.Distance finDistance(displayUnit = "mm") "distance between fins" annotation(
      Dialog(tab = "Fins data"));
    parameter Boolean finIsCircular = true "if false rectangular fin is used" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finDiameter(displayUnit = "mm") "fin external diameter, if circular" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finWidth(displayUnit = "mm") "side length perpenticular to flow, if rectangular" annotation(
      Dialog(tab = "Fins data"));
    parameter Modelica.Units.SI.Distance finHeight(displayUnit = "mm") "side length along flow, if rectangular" annotation(
      Dialog(tab = "Fins data"));    
    parameter Boolean useBYnusselt = false "if true, the Briggs and Youg Nusselt calc. will be used, otherwise the Ganguli's one" annotation(
      Dialog(tab = "Heat transfer"));
    Modelica.Units.SI.Area IoAreaRaw(start = 1.0) "total external surface of tubes without the fins";
    Modelica.Units.SI.Area FinsArea(start = 1.0) "total external surface of fins";
    Modelica.Units.SI.Area OiArea(start = 1.0) "total external surface of tubes and fins";
    Modelica.Units.SI.Area Sface(start = 1.0) "total area perpenticular to gas flow";
    Modelica.Units.SI.Area Sflow(start = 1.0) "total free area available for gas flow";
    Real OFinRatio;
    Real OPhi;
    MediumO.ThermodynamicState OStateAvg "outer flow average state for physical properties calculation";
    Modelica.Units.SI.Density ORho(displayUnit = "kg/m3") "outer flow average density";
    Modelica.Units.SI.DynamicViscosity OMu(min = 1e-6, start = 1e-3, max = 1e6) "outer flow average dynamic viscosity";
    Modelica.Units.SI.SpecificHeatCapacity OCp(start = 2000.0) "outer flow average heat capacity";
    Modelica.Units.SI.ThermalConductivity OK(start = 0.1) "outer flow average thermal conductivity";
    Modelica.Units.SI.PrandtlNumber OPr "outer flow Prandt number";
    Modelica.Units.SI.VolumeFlowRate OQ(displayUnit = "m3/h") "total outer volume flow rate at average conditions";
    Modelica.Units.SI.Velocity Vmax "gas velocity at the free surface";
    Modelica.Units.SI.Velocity Vface "gas velocity as per Sface. Typical: 2.5-3.5 m/s";
    Modelica.Units.SI.ReynoldsNumber ORe(min = 0.01, start = 20000) "outer flow average Reynolds number";
    Real Oa;
    Modelica.Units.SI.ReynoldsNumber OReff "outer flow effective Reynolds number";
    Real OF(start = 0.01) "outer flow Fanning's friction factor at average conditions";
    Modelica.Units.SI.NusseltNumber ONuGan "gas side Nusselt number as per Ganguli";
    Modelica.Units.SI.NusseltNumber ONuBY "gas side Nusselt as for Briggs and Young"; 
    Real OX;
    
  algorithm
    IoAreaRaw:=ItotalNumPipes*iLTube*IoPerimeter;
    Sface := iLTube * iTubePitch * ItotalNumPipes / iNumRows;
  if finIsCircular == true then
      FinsArea := (0.5 * pi * finDiameter ^ 2 - 2 * IoSection) * iLTube / finDistance * ItotalNumPipes;
      Sflow := (iTubePitch - iDop - (finDiameter - iDop) * finThickness / finDistance) * iLTube * ItotalNumPipes / iNumRows;
      OFinRatio := finDiameter / iDop;
    else
      FinsArea := 2 * (finWidth * finHeight - IoSection) * iLTube / finDistance * ItotalNumPipes;
      Sflow := (iTubePitch - iDop - (finWidth - iDop) * finThickness / finDistance) * iLTube * ItotalNumPipes / iNumRows;
      if iStaggered == false then
        OFinRatio := 1.28 * min(finWidth, finHeight) / iDop * (max(finWidth, finHeight) / min(finWidth, finHeight) - 0.2) ^ 0.5 "as per Schmidt, Zabronsky, and
          Rich, Air Cooled Heat Exchangers, D.G.Kroger. But rectified according to VDI Atlas";
      else
        OFinRatio := 1.27 * min(finWidth, finHeight) / iDop * (max(finWidth, finHeight) / min(finWidth, finHeight) - 0.3) ^ 0.5 "as per Schmidt, Zabronsky, and
          Rich, Air Cooled Heat Exchangers, D.G.Kroger. According to VDI Atlas this is not totally correct";
      end if;
    end if;
//Sflow := ((iTubePitch - iDop) * (finDistance - finThickness) - (iTubePitch - finDiameter) * finThickness) * iLTube / finDistance * (ItotalNumPipes / iNumRows);
    OPhi := (OFinRatio - 1) * (1 + 0.35 * log(OFinRatio));
  
    if useFins==true then
      IoArea:=ItotalNumPipes*iLTube*IoPerimeter*(finDistance-finThickness)/finDistance;    
      OiArea := FinsArea + IoArea;
    else
      Sflow:=Sface;
      IoArea:=ItotalNumPipes*iLTube*IoPerimeter;
      OiArea := IoArea;         
    end if;
//missing quadrangular
  equation
    OStateAvg = MediumO.setState_phX((Oin.P + Oout.P) / 2, (Oin.H + Oout.H) / 2, Oin.X);
    OT = MediumO.temperature(OStateAvg);   
    ORho = abs(MediumO.density(OStateAvg));
    OMu = MediumO.dynamicViscosity(OStateAvg);
    OCp = MediumO.specificHeatCapacityCp(OStateAvg);
    OK = MediumO.thermalConductivity(OStateAvg);
    OPr = OCp * OMu / OK;
    OQ = Oin.G / ORho;
    Vmax = OQ / Sflow;
    Vface = OQ / Sface;
    ORe = iDop * Vmax * ORho / OMu;
    Oa = (iTubePitch - finDiameter) / iDop;
    OReff = 2 * ORe * (finDistance - finThickness) / (finDiameter - iDop);
    OF = (1 + 2 * exp(Oa / 4) / (1 + Oa)) * (0.021 + 27.2 / OReff + 0.29 / OReff ^ 0.2);
    OPlossCorr=1.0;
    OPloss = 2 * OF * iNumRows * Vmax ^ 2 * ORho;
    Oout.P-Oin.P = (-sign(Oin.G) * OPloss) + (Oin.Elevation - Oout.Elevation + 1e-5) * g_n * ORho "momentum change is not taken into account. 1 e-5 is to avoid division by 0";
    
    if iNumRows < 4 then
      if iStaggered == false then
        ONuGan = 0.2 * ORe ^ 0.6 * OPr ^ 0.33333 * (OiArea / IoAreaRaw) ^ (-0.15);
      else
        ONuGan = if iNumRows == 3 then 0.36 * ORe ^ 0.6 * OPr ^ 0.33333 * (OiArea / IoAreaRaw) ^ (-0.15) else 0.33 * ORe ^ 0.6 * OPr ^ 0.33333 * (OiArea / IoAreaRaw) ^ (-0.15);
      end if;
    else
      ONuGan = if iStaggered == true then 0.38 * ORe ^ 0.6 * OPr ^ 0.33333 * (OiArea / IoAreaRaw) ^ (-0.15) else 0.22 * ORe ^ 0.6 * OPr ^ 0.33333 * (OiArea / IoAreaRaw) ^ (-0.15);
    end if;
    ONuBY = 0.134 * ORe ^ 0.681 * OPr ^ 0.3333 * ((finDistance - finThickness) / (finDiameter - iDop)) ^ 0.2 * ((finDistance - finThickness) / finThickness) ^ 0.1134;
    if useBYnusselt == true then
      ONu = ONuBY;
    else
      ONu = ONuGan;
    end if;
      OHi = ONu * OK / iDop;
    if useFins==true then
      OX = 0.5 * OPhi * iDop * (2 * OHi / (finK * finThickness)) ^ 0.5;
      OiEfficiency=1 - (1 - tanh(OX) / OX) * FinsArea / OiArea;
    else
      OX=0;
      OiEfficiency=1;
    end if;
  
    ILMTDcorr=homotopy(FreeFluids.HeatExchangers.Functions.CrossLMTDfactor(iNumSerial, iNumRows, Rntu, NTU),1);
  annotation(
      Icon(graphics = {Line(origin = {-79.8217, -0.588363}, points = {{0, 20}, {0, -20}}), Line(origin = {79.0307, -0.479763}, points = {{0, 20}, {0, -20}, {0, -20}}), Line(origin = {-0.682365, -0.371163}, points = {{-80, 0}, {80, 0}, {80, 0}}), Line(origin = {-70.3095, -0.766643}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-60.705, -0.653653}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-40.818, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-51.2135, -0.653653}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-0.761468, -0.427663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-11.157, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-20.6485, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-30.253, -0.653653}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {68.4475, -0.427663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {58.956, -0.427663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {49.3515, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {38.7865, -0.427663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {28.391, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {18.8995, -0.540663}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {9.29505, -0.653653}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {68.5605, 19.0638}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-0.569375, 19.1203}, points = {{-80, 0}, {80, 0}, {80, 0}}), Line(origin = {59.069, 19.0638}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {49.4645, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {38.8995, 19.0638}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {28.504, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {19.0125, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {9.40803, 18.8378}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-0.648474, 19.0638}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-11.044, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-20.5355, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-30.14, 18.8378}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-40.705, 18.9508}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-51.1005, 18.8378}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-60.592, 18.8378}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-70.1965, 18.7248}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {68.5605, -20.8797}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-0.56937, -20.8232}, points = {{-80, 0}, {80, 0}, {80, 0}}), Line(origin = {59.069, -20.8797}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {49.4645, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {38.8995, -20.8797}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {28.504, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {19.0125, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {9.40803, -21.1057}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-0.648474, -20.8797}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-11.044, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-20.5355, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-30.14, -21.1057}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-40.705, -20.9927}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-51.1005, -21.1057}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-60.592, -21.1057}, points = {{0, 4}, {0, -4}, {0, -4}}), Line(origin = {-70.1965, -21.2187}, points = {{0, 4}, {0, -4}, {0, -4}}), Ellipse(origin = {-33, -74}, extent = {{-33, 6}, {33, -6}}), Ellipse(origin = {33, -74}, extent = {{-33, 6}, {33, -6}}), Line(origin = {-59.7131, -49.5658}, points = {{0, -14}, {0, 6}, {0, 6}}), Polygon(origin = {-61, -38}, points = {{1, 6}, {-5, -6}, {7, -6}, {1, 6}, {1, 6}}), Line(origin = {-29.5437, -49.4527}, points = {{0, -14}, {0, 6}, {0, 6}}), Polygon(origin = {-31, -38}, points = {{1, 6}, {-5, -6}, {7, -6}, {1, 6}, {1, 6}}), Line(origin = {30.9692, -49.5179}, points = {{0, -14}, {0, 6}, {0, 6}}), Polygon(origin = {29, -38}, points = {{1, 6}, {-5, -6}, {7, -6}, {1, 6}, {1, 6}}), Polygon(origin = {-1, -38}, points = {{1, 6}, {-5, -6}, {7, -6}, {1, 6}, {1, 6}}), Line(origin = {0.0476988, -49.9179}, points = {{0, -14}, {0, 6}, {0, 6}}), Polygon(origin = {59, -38}, points = {{1, 6}, {-5, -6}, {7, -6}, {1, 6}, {1, 6}}), Line(origin = {60.4952, -49.6221}, points = {{0, -14}, {0, 6}, {0, 6}}), Text(origin = {72, 53}, lineColor = {22, 41, 240}, extent = {{-58, -25}, {58, 25}}, textString = "%name")}),
      Documentation(info = "<html><head></head><body>This partial model performs the calculation of the heat exchange coefficient, and pressure drop, of the gas part (normally air) of an air cooled, finned tubes, exchanger.<div>The LMTD correction factor is also calculated.</div></body></html>"));
  end GasCooledHEX;
  
  model ShellAndTubesHEX
    extends HEXtubesDefinition(final iNumParallel=1, final iIsCircular=true, final iPerimeter=0, final iSection=0, final useFins=false, final finRho=0, final finK=0, final finThickness=0);
    parameter Real pitchRatio = 1.33 "ratio between tube pitch and tube diameter. Normally 1.25-1.5" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Boolean triangularPitch = true "if false, square layout" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Boolean rotatedPitch = false "if true, rotated 45º square layout" annotation(
      Dialog(tab = "Inner pipe data"));
    parameter Integer numShells = 1 "total number of serial shells" annotation(
      Dialog(tab = "Shell data"));
    parameter Types.TemaShell shellType=Types.TemaShell.E "shell type" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Diameter dShI(displayUnit = "mm") = 0 "internal shell diameter" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.ShearStress shellDTS = 60.0e6 "shell design tensile stress. AISI316L=60e6. AISI316=68.9e6" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Diameter shellNozzlesDi(displayUnit = "mm") "diameter of shell nozzles" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance cTuBa(displayUnit = "mm") = 0.8e-3 "tube-baffle hole diameter clearance. Typical: 0.8 mm" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance cBaSh(displayUnit = "mm") = 1.6e-3+0.004*dShI "baffle-shell diameter clearance. Typical 5 mm" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance cBuSh(displayUnit = "mm") = 35e-3 "tube bundle-shell diameter clearance. Typical: 35 mm" annotation(
      Dialog(tab = "Shell data"));
    parameter Integer numBa = 23 annotation(
      Dialog(tab = "Shell data"));
    parameter Boolean equalBaffleDistance = true "if true, all baffles are at equal distance. If false, inlet and outlet are different" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance inletBaffleDistance(displayUnit = "mm") = 0 "if equalBaffleDistance=false" annotation(
      Dialog(tab = "Shell data"));
    parameter SI.Distance outletBaffleDistance(displayUnit = "mm") = 0 "if equalBaffleDistance=false" annotation(
      Dialog(tab = "Shell data"));
    parameter Fraction BaCutRatio = 0.25 "between 0.15 and 0.45 of shell diameter. Normally 0.2-0.35" annotation(
      Dialog(tab = "Shell data"));
    parameter Integer numSS = 0 "number of sealing strips pairs" annotation(
      Dialog(tab = "Shell data"));
    Integer TubesInShell;
    Integer PassesInShell;
    Modelica.Units.SI.Diameter IDo "inner pipe external diameter";
    Real TubeCountFactor;
    Real TubeLayoutFactor;
    SI.Distance Pitch "should be at least 6 mm more than IDo";
    SI.Distance PitchP "pitch parallel to flow";
    SI.Distance PitchN "pitch normal to flow";
    SI.Diameter DshIProposal "internal shell diameter recommended from number of tubes";
    SI.Diameter DbuO "tube bundle external diameter";
    SI.Distance BaCut "max distance from baffle to shell";
    SI.Distance BaSep "inter baffles distance. Minimum 1/5 of shell diameter or 5 cms. max: 29.54*IDo^0.75. Optimum: 0.3-0.6 of dShI";
    SI.Angle ThetaBuC "angle of baffle cut intersection with outer tubes centers diameter";
    SI.Angle ThetaSh "angle of baffle cut intersection with shell";
    Fraction Fw "fraction of total tubes in each window";
    Fraction Fc "fraction of total tubes in full crossflow";
    SI.Area OSm "crossflow area at axis of exchanger";
    SI.Area OSBaSh "baffle to shell leakage area for each baffle";
    SI.Area OSTuBa "tube to baffle leakage area for each baffle";
    //For multiple pass tubes, OSb must be corrected if partition lane is parallel to flow
    SI.Area OSb "bundle bypass area";
    Real Nc "number of rows in full crossflow";
    Real Ncw "number of effective crossflow rows in each window";
    SI.Area OSw "free area for shell flow at windows";
  equation
    TubesInShell=ItotalNumPipes/numShells;
    PassesInShell=iNumSerial/numShells;
    IDo=iDi+2*iThickness;
    if iNumSerial == 1 then
      TubeCountFactor = 0.93;
    elseif iNumSerial == 2 then
      TubeCountFactor = 0.9;
    elseif iNumSerial == 3 then
      TubeCountFactor = 0.85;
    else
      TubeCountFactor = 0.825;
    end if;
    if triangularPitch == true then
      TubeLayoutFactor = 0.866;
    else
      TubeLayoutFactor = 1;
    end if;
    BaCut = dShI * BaCutRatio;
    if equalBaffleDistance == true then
      BaSep = iLTube / (numBa + 1);
    else
      BaSep = (iLTube - inletBaffleDistance - outletBaffleDistance) / (numBa - 1);
    end if;
    DbuO = dShI - cBuSh;
    ThetaBuC = pi - 2 * asin((dShI / 2 - BaCut) * 2 / (DbuO - IDo));
    ThetaSh = pi - 2 * asin((dShI / 2 - BaCut) * 2 / dShI);
    Fw = ThetaBuC / (2 * pi) - sin(ThetaBuC) / (2 * pi);
    Fc = 1 - 2 * Fw;
    DshIProposal = (ItotalNumPipes * Pitch ^ 2 * TubeLayoutFactor / TubeCountFactor * 4 / pi) ^ 0.5;
  //baffle leakage correction factor, streams between baffle and tubes, and between baffle and shell
    Pitch = pitchRatio * IDo;
    if triangularPitch == true then
      if rotatedPitch == false then
        PitchP = Pitch * 0.866;
        PitchN = Pitch * 0.5;
      else
        PitchP = Pitch * 0.5;
        PitchN = Pitch * 0.866;
      end if;
    else
      if rotatedPitch == false then
        PitchP = Pitch;
        PitchN = Pitch;
      else
        PitchP = Pitch * 0.7071;
        PitchN = Pitch * 0.7071;
      end if;
    end if;
    Nc = (dShI - 2 * BaCut) / PitchP "number of rows in full crossflow";
  //Fc=1.2-2.27*BaCut/dShI;
  //Fc=1/pi*(pi+2*(dShI-2*BaCut)/DbuO*sin(acos((dShI-2*BaCut)/DbuO))-2*acos((dShI-2*BaCut)/DbuO));
    Ncw = 0.8 * BaCut / PitchP "number of effective rows in window crossflow";
    if triangularPitch == false then
      OSm = BaSep * (cBuSh + (DbuO - IDo) / PitchN * (Pitch - IDo));
    else
      OSm = BaSep * (cBuSh + (DbuO - IDo) / Pitch * (Pitch - IDo));
    end if;
    OSBaSh = (dShI ^ 2 - (dShI - cBaSh) ^ 2) * (2 * pi - ThetaSh) / 8;
    OSTuBa = pi / 4 * ((IDo + cTuBa) ^ 2 - IDo ^ 2) * ItotalNumPipes * (1 - Fw);
    OSb = BaSep * cBuSh;
    OSw = dShI * dShI * (ThetaSh / 8 - abs(sin(ThetaSh / 2) * cos(ThetaSh / 2)) / 4) - ItotalNumPipes * Fw * pi * IDo ^ 2 / 4;
  
    IoArea=ItotalNumPipes*iLTube*IoPerimeter;  
    OiArea=IoArea;  
    OiEfficiency=1;
   
  annotation(
      Icon(graphics = {Rectangle(origin = {2, 0}, rotation = 90, fillColor = {85, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-30, 60}, {30, -60}}), Ellipse(origin = {60, 0}, rotation = 90, fillColor = {85, 170, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-30, 6}, {30, -10}}), Ellipse(origin = {-60, 0}, rotation = 90, fillColor = {85, 170, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-30, 6}, {30, -10}}), Line(origin = {0.02, 23.57}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {-58.5738, 0}, points = {{0, -30}, {0, 30}}), Line(origin = {60.7295, -0.0409836}, points = {{0, -30}, {0, 30}}), Line(origin = {0.27, 16.07}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {0.23, 0.53}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {-0.02, 8.32}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {0.56, -7.74}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {0.52, -15.82}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {0.52, -23.85}, rotation = -90, points = {{0, 60}, {0, -58}}, thickness = 1), Line(origin = {-72.5, 0}, points = {{6.5, 0}, {-7.5, 0}, {-5.5, 0}}), Line(origin = {75, 0}, points = {{-5, 0}, {5, 0}}), Line(origin = {0, -55}, points = {{0, -25}, {0, 25}}), Line(origin = {0, 55}, points = {{0, 25}, {0, -25}, {0, -25}}), Text(origin = {53, 49}, lineColor = {15, 10, 253}, extent = {{-43, 15}, {43, -15}}, textString = "%name")}));end ShellAndTubesHEX;

  partial model ShellAndTubesHEXfc
    extends FreeFluids.HeatExchangers.BaseClasses.ShellAndTubesHEX;
    Real Jc "bafflecut correction factor";
    Real Jl "baffle leakage correction factor";
    Real Jb "Bundle bypass correction factor";
    Real Cbh "empirical factor for bundle bypass correction, depending on Reynolds number";
    Real Js "unequal baffle spacing correction factor";
    Real Jm "wall viscosity correction factor";
    Real N "exponent for baffle spacing correction"; 
    
    Real A1, A2, A3, A4 "empirical constants for ideal heat transfer calculation";
    Real B1, B2, B3, B4 "empirical constants for pressure loss calculation"; 
  
    MediumO.ThermodynamicState OStateAvg "shell average state for physical properties calculation";
    MediumO.ThermodynamicState OiStateWall "shell internal wall state for correction factors calculation";
    Modelica.Units.SI.Density ORho(displayUnit = "kg/m3") "shell average density";
    Modelica.Units.SI.DynamicViscosity OMu(min = 1e-6, start = 1e-3, max = 1e6) "annulus average dynamic viscosity";
    Modelica.Units.SI.SpecificHeatCapacity OCp(start = 2000.0) "shell average heat capacity";
    Modelica.Units.SI.ThermalConductivity OK(start = 0.1) "shell average thermal conductivity";  
    Modelica.Units.SI.DynamicViscosity OiMuWall(min = 1e-6, start = 1e-3, max = 1e6) "shell inner wall dynamic viscosity"; 
    Modelica.Units.SI.PrandtlNumber OPr "shell Prandt number";
    Modelica.Units.SI.VolumeFlowRate OQ(displayUnit = "m3/h") "total shell volume flow rate at average conditions";
    Modelica.Units.SI.Velocity OV(start = 1) "shell flow velocity at average conditions. Normally between 0.9 and 3.0 m/s for liquids";
    Modelica.Units.SI.ReynoldsNumber ORe(min = 0.01, start = 20000) "shell average Reynolds number";
    Modelica.Units.SI.CoefficientOfHeatTransfer OHiIdeal(min = 1, start = 1000) "outer flow ideal heat transfer coefficient to inner pipe";
    Real OFcf "crossflow friction factor";    
    SI.Pressure OPlossCfId(start = 0) "ideal crossflow friction head loss in each baffle compartment";
    SI.Pressure OPlossCf(start = 0) "crossflow friction loss for all baffle compartments";
    Real Cbp "empirical factor for bundle bypass correction, depending on Reynolds number";
    Real Rb "bypass correction factor";
    Real Rl "leakage correction factor";
    SI.Pressure OPlossWin(start = 0) "windows friction loss for all baffle compartments";
    SI.Pressure OPlossEn(start = 0) "entrance friction loss";
    SI.Pressure OPlossNoz(start = 0) "nozzles friction loss";
    Modelica.Units.SI.Velocity OVnoz(start = 1) "nozzles velocity";
    Modelica.Units.SI.ReynoldsNumber OReNoz(min = 0.01, start = 20000) "nozzles average Reynolds number";
    
  algorithm
    if triangularPitch == true then
      if rotatedPitch == false then
        A3 := 1.45;
        A4 := 0.519;
        B3 := 7.0;
        B4 := 0.5;
        if ORe > 1e4 then
          A1 := 0.321;
          A2 := -0.388;
          B1 := 0.372;
          B2 := -0.123;
        elseif ORe > 1e3 then
          A1 := 0.321;
          A2 := -0.388;
          B1 := 0.486;
          B2 := -0.152;
        elseif ORe > 1e2 then
          A1 := 0.593;
          A2 := -0.477;
          B1 := 4.57;
          B2 := -0.476;
        elseif ORe > 1e1 then
          A1 := 1.36;
          A2 := -0.657;
          B1 := 45.1;
          B2 := -0.973;
        else
          A1 := 1.4;
          A2 := -0.667;
          B1 := 48.0;
          B2 := -1.0;
        end if;
      else
        A3 := 0;
        A4 := 0;
        B3 := 0;
        B4 := 0;
        A1 := 0;
        A2 := 0;
        B1 := 0;
        B2 := 0;
      end if;
    else
      if rotatedPitch == false then
        A3 := 1.187;
        A4 := 0.37;
        B3 := 6.3;
        B4 := 0.378;
        if ORe > 1e4 then
          A1 := 0.37;
          A2 := -0.395;
          B1 := 0.391;
          B2 := -0.148;
        elseif ORe > 1e3 then
          A1 := 0.107;
          A2 := -0.266;
          B1 := 0.0815;
          B2 := 0.022;
        elseif ORe > 1e2 then
          A1 := 0.408;
          A2 := -0.46;
          B1 := 6.09;
          B2 := -0.602;
        elseif ORe > 1e1 then
          A1 := 0.9;
          A2 := -0.631;
          B1 := 32.1;
          B2 := -0.963;
        else
          A1 := 0.97;
          A2 := -0.667;
          B1 := 35;
          B2 := -1;
        end if;
      else
        A3 := 1.93;
        A4 := 0.5;
        B3 := 6.59;
        B4 := 0.52;
        if ORe > 1e4 then
          A1 := 0.37;
          A2 := -0.396;
          B1 := 0.303;
          B2 := -0.126;
        elseif ORe > 1e3 then
          A1 := 0.37;
          A2 := -0.396;
          B1 := 0.333;
          B2 := -0.136;
        elseif ORe > 1e2 then
          A1 := 0.73;
          A2 := -0.5;
          B1 := 3.5;
          B2 := -0.476;
        elseif ORe > 1e1 then
          A1 := 0.498;
          A2 := -0.656;
          B1 := 26.2;
          B2 := -0.913;
        else
          A1 := 1.55;
          A2 := -0.667;
          B1 := 32;
          B2 := -1;
        end if;
      end if;
    end if;
  equation
    OStateAvg = MediumO.setState_phX((Oin.P + Oout.P) / 2, (Oin.H + Oout.H) / 2, Oin.X);
    OiStateWall = MediumO.setState_pTX((Oin.P + Oout.P) / 2, ITwall, Oin.X);
    OT = MediumO.temperature(OStateAvg);   
    ORho = abs(MediumO.density(OStateAvg));
    OMu = MediumO.dynamicViscosity(OStateAvg);
    OCp = MediumO.specificHeatCapacityCp(OStateAvg);
    OK = MediumO.thermalConductivity(OStateAvg);
    OiMuWall = MediumO.dynamicViscosity(OiStateWall);
    OPr = OCp * OMu / OK;
    OQ = Oin.G / ORho;
    OV = OQ/iNumParallel / OSm;
    ORe = IDo * abs(OV) * ORho / OMu;
    if ORe > 100 then
      Cbh = 1.25;
      N = 0.6;
    else
      Cbh = 1.35;
      N = 0.333;
    end if;
    Jb = exp(-Cbh * OSb / OSm * (1 - (2 * numSS / Nc) ^ 0.333));
    if equalBaffleDistance == true then
      Js = 1;
    else
      Js = (numBa - 1 + (inletBaffleDistance / BaSep) ^ (1 - N) + (outletBaffleDistance / BaSep) ^ (1 - N)) / (numBa - 1 + inletBaffleDistance / BaSep + outletBaffleDistance / BaSep);
    end if;
    if ORho > 300 then
      Jm = (OMu / OiMuWall) ^ 0.11;
    else
      if Oin.H > Oout.H then
        Jm = (OT / ITwall) ^ 0.25;
      else
        Jm = 1.0;
      end if;
    end if;
  //bafflecut correction factor
    Jc = 0.55 + 0.72 * Fc;
  //baffle leak correction factor
    Jl = 0.44 * (1 - OSBaSh / (OSBaSh + OSTuBa)) + (1 - 0.44 * (1 - OSBaSh / (OSBaSh + OSTuBa))) * exp(-2.2 * (OSBaSh + OSTuBa) / OSm);
    
    OHiIdeal = A1 * (1.33 * IDo / Pitch) ^ (A3 / (0.14 * ORe ^ A4)) * ORe ^ A2 * OCp * Oin.G / OSm * OPr ^ (-0.666667);
    
    OHi = OHiIdeal * Jc * Jl * Jb * Js * Jm;
  
  //Pressure loss calculation
    OFcf=B1 * (1.33 * IDo / Pitch) ^ (B3 / (1 + 0.14 * ORe ^ B4)) * ORe ^ B2;
    OPlossCfId = 2 * OFcf * Nc * (Oin.G / OSm) ^ 2/ ORho * OPlossCorr;
    if ORe > 100 then
      Cbp = 3.7;
    else
      Cbp = 4.5;
    end if;
    Rb = exp(-Cbp * OSb / OSm * (1 - (2 * numSS / Nc) ^ 0.333));
    Rl = exp(-1.33 * (1 + OSBaSh / (OSBaSh + OSTuBa)) * ((OSBaSh + OSTuBa) / OSm) ^ ((-0.15 * (1 + OSBaSh / (OSBaSh + OSTuBa))) + 0.8));
    OPlossCf = OPlossCfId * (numBa - 1) * Rb * Rl;
    OPlossWin = numBa * ((2 + 0.6 * Ncw) * Oin.G ^ 2 / (OSm * OSw * ORho * 2)) * Rl;
    OPlossEn = 2 * ORho * (OQ / (pi * shellNozzlesDi ^ 2 / 4)) ^ 2; 
    OVnoz=4*OQ/(numShells*pi*shellNozzlesDi^2);
    OReNoz=shellNozzlesDi*OVnoz*ORho/OMu;
    if OReNoz>1e4 then
      OPlossNoz=0.75*numShells*OVnoz^2*ORho;
    else
      OPlossNoz=1.5*numShells*OVnoz^2*ORho;  
    end if;        
    OPloss = OPlossCf + OPlossWin + OPlossEn+OPlossNoz;
    
  
    OPlossCorr=1;
    Oout.P-Oin.P = (-sign(Oin.G) * OPloss) + (Oin.Elevation - Oout.Elevation + 1e-5) * g_n * ORho "momentum change is not taken into account. 1 e-5 is to avoid division by 0";
    ONu=0;
    if (iNumSerial==1 and counterCurrent==true) then
      ILMTDcorr=1.0;
    else
      ILMTDcorr=homotopy(FreeFluids.HeatExchangers.Functions.ShellLMTDfactor(shellType, numShells, iNumSerial, Rntu, NTU),1);
    end if;
    //ILMTDcorr=1;
  annotation(
      Documentation(info = "<html><head></head><body>This partial model performs the calculation of the heat transfer coefficient and the pressure drop of a TEMA E shell, according to the Bell-Delaware methodology.<div>The LMTD correction factor is also calculated.</div></body></html>"));end ShellAndTubesHEXfc;
  annotation(
    Documentation(info = "<html><head></head><body>Partial classes for detailed heat exchangers.</body></html>"));
end BaseClasses;
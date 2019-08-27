within FreeFluids;

package MediaCommon
 "MediaCommon.mo by Carlos Trujillo
      This file is part of the Free Fluids application
      Copyright (C) 2008-2019  Carlos Trujillo Gonzalez
        
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

    record DataRecord
      String name="", description="", CAS="";
      Integer family=0;
      Real MW = 0.0, Tc = 0.0, Pc = 0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 6.667e-10, lnuA=0.0, lnuB=0.0;
      Integer Cp0Corr=0, VpCorr=0, BtCorr=0, HvCorr=0, lDensCorr=0, lCpCorr=0, lTfromHsatCorr=0, lViscCorr=0, lThCondCorr=0, lSurfTensCorr=0, lBulkModRCorr=0, gSatDensCorr=0, gViscCorr=0, gThCondCorr=0;
      Real Cp0Coef[13]={0,0,0,0,0,0,0,0,0,0,0,0,0}, VpCoef[6]={0,0,0,0,0,0}, BtCoef[6]={0,0,0,0,0,0}, HvCoef[6]={0,0,0,0,0,0}, lDensCoef[6]={0,0,0,0,0,0}, lCpCoef[6]={0,0,0,0,0,0}, lTfromHsatCoef[6]={0,0,0,0,0,0}, lViscCoef[6]={0,0,0,0,0,0}, lThCondCoef[6]={0,0,0,0,0,0}, lSurfTensCoef[6]={0,0,0,0,0,0}, lBulkModRCoef[6]={0,0,0,0,0,0}, gSatDensCoef[6]={0,0,0,0,0,0}, gViscCoef[6]={0,0,0,0,0,0}, gThCondCoef[6]={0,0,0,0,0,0};
      Real Cp0LimI=0, VpLimI=0, BtLimI=0, HvLimI=0, lDensLimI=0, lCpLimI=0, lTfromHsatLimI=0, lViscLimI=0, lThCondLimI=0, lSurfTensLimI=0, lBulkModRLimI=0, gSatDensLimI=0, gViscLimI=0, gThCondLimI=0;
      Real Cp0LimS=0, VpLimS=0, BtLimS=0, HvLimS=0, lDensLimS=0, lCpLimS=0, lTfromHsatLimS=0, lViscLimS=0, lThCondLimS=0, lSurfTensLimS=0, lBulkModRLimS=0, gSatDensLimS=0, gViscLimS=0, gThCondLimS=0;
    end DataRecord;
  
    constant DataRecord MediaDataTemplate( name="", description="", CAS="", family=0, MW = 0.0, Tc = 0.0, Pc = 0.0, Vc = 0.0, Zc = 0.0, w = 0.0, Tb = 0.0, mu = 0.0, IsothComp = 0.0, lnuA=0.0, lnuB=0.0,
	Cp0Corr = 0, Cp0Coef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, Cp0LimI = 0.0, Cp0LimS = 0.0,
	VpCorr = 0, VpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, VpLimI = 0.0, VpLimS = 0.0,
	BtCorr = 0, BtCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, BtLimI = 0.0, BtLimS = 0.0,
	HvCorr = 0, HvCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, HvLimI = 0.0, HvLimS = 0.0,
	lDensCorr = 0, lDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lDensLimI = 0.0, lDensLimS = 0.0,
	lCpCorr = 0, lCpCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lCpLimI = 0.0, lCpLimS = 0.0,
	lTfromHsatCorr = 0, lTfromHsatCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lTfromHsatLimI = 0.0, lTfromHsatLimS = 0.0,
	lViscCorr = 0, lViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lViscLimI = 0.0, lViscLimS = 0.0,
	lThCondCorr = 0, lThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lThCondLimI = 0.0, lThCondLimS = 0.0,
	lSurfTensCorr = 0, lSurfTensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lSurfTensLimI = 0.0, lSurfTensLimS = 0.0, lBulkModRCorr = 0, lBulkModRCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, lBulkModRLimI = 0.0, lBulkModRLimS = 0.0, gSatDensCorr = 0, gSatDensCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gSatDensLimI = 0.0, gSatDensLimS = 0.0, gViscCorr = 0, gViscCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gViscLimI = 0.0, gViscLimS = 0.0, gThCondCorr = 0, gThCondCoef = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, gThCondLimI = 0.0, gThCondLimS = 0.0);
   

  type InputChoice = enumeration(
    pT "(p,T) as inputs",
    ph "(p,h) as inputs",
    ps "(p,s) as inputs",    
    dT "(d,T) as inputs");
  
package Functions
  function PhysPropCorr "Calculates a physical property from a given correlation data and a given temperature, or pressure, using an external function. Should be the general form of calculation, but is not working in OpenModelica 1.14 when iterating in an array of DataRecord"
    input Integer corr;
    input Real coef[:];
    input Real MW;
    input Real x;
    output Real y;
  
    external "C" FF_PhysPropCorr(corr, coef, MW, x, y) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFphysprop.c\"");
  end PhysPropCorr;

  function Cp0Corr "Calculates ideal gas heat capacity from a given DataRecord and a given temperature, using an external function. The intermediate pass of the protected variable is necessary for OpenModelica 1.14"
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real y;
  protected
    Real coef[13]=data.Cp0Coef; 
    external "C" FF_PhysPropCorr(data.Cp0Corr, coef, data.MW, T, y) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFphysprop.c\"");
  end Cp0Corr;

  function FF_Cp0Corr "Calculates ideal gas heat capacity from a given DataRecord and a given temperature, using DIPPR107 correlation in J/(kg·K)"
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real y;
  algorithm
      y:=data.Cp0Coef[1]+data.Cp0Coef[2]*(data.Cp0Coef[3]/(T*sinh(data.Cp0Coef[3]/T)))^2+data.Cp0Coef[4]*(data.Cp0Coef[5]/(T*cosh(data.Cp0Coef[5]/T)))^2;
  end FF_Cp0Corr;

  function SpecificEnthalpyCorr2 "Calculates specific enthalpy from a given DataRecord at a given temperature, using an external function.  The intermediate pass of the protected variable is necessary for OpenModelica 1.14"
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real h;
  protected
    Real coef[13]=data.Cp0Coef;
      external "C" FF_SpecificEnthalpyCorr(data.Cp0Corr, coef, data.MW, T, h) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFphysprop.c\"");
  end SpecificEnthalpyCorr2;

  function SpecificEntropyCorr2 "Calculates specific enthalpy from a given DataRecord at a given temperature. As per specific enthalpy "
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real s;
  protected
    Real coef[13]=data.Cp0Coef;
    external "C" FF_SpecificEntropyCorr(data.Cp0Corr, coef, data.MW, T, s) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFphysprop.c\"");
  end SpecificEntropyCorr2;
  
    function gasViscCorr "Calculates low pressure gas viscosity from a given DataRecord and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gViscCoef directly, coming from an iteration."
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real eta;
    protected
      Real coef[6]= data.gViscCoef;
    external "C" FF_PhysPropCorr(data.gViscCorr, coef, data.MW, T, eta) annotation(
      IncludeDirectory = "modelica://FreeFluids/Resources",
      Include = "#include \"FFphysprop.c\"");
  end gasViscCorr;
  
  function FF_gasViscCorr "Calculates low pressure gas viscosity from a given DataRecord and a given temperature. Faster than gasViscCorr"
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real eta;
  algorithm
    eta:= if data.gViscCorr==111 then data.gViscCoef[1]+T*(data.gViscCoef[2]+T*(data.gViscCoef[3]+T*(data.gViscCoef[4]+T*data.gViscCoef[4]))) else data.gViscCoef[1]*T^data.gViscCoef[2]/(1+data.gViscCoef[3]/T+data.gViscCoef[4]/T^2);
  end FF_gasViscCorr;
    
  function gasViscLowPressureChung "Dynamic viscosity of a low pressure gas according to Chung"
    input FreeFluids.MediaCommon.DataRecord data;
    input SI.Temp_K T "Gas temperature";
    output SI.DynamicViscosity eta "Dynamic viscosity";
  protected
    Real Ta, sigma, omega,muR, muR4,Fc,k;
  algorithm
    if data.family==6 then
      k:=0.076 "water";
    elseif data.family==7 then
      k:=0.0682+4.74/data.MW "alcohol";
    elseif data.family==8 then
      k:=0.0682+2*4.74/data.MW "polyol";
    else
      k:=0.0;
    end if;
    Ta := 1.2593*T / data.Tc;
    sigma:= 0.809*data.Vc^0.33333;
    omega:= 1.16145*Ta^(- 0.14874)+0.52487*exp(- 0.7732*Ta)+2.16178*exp(- 2.43787*Ta);
    muR := if (data.mu < 999.0 and data.mu > 0.0) then 131.3*data.mu/(data.Vc*1e6*data.Tc)^0.5 else 0.0;
    muR4:=muR*muR*muR*muR;
    Fc:=1-0.2756*data.w+0.059035*muR4+k;
    eta := 26.692*Fc*(data.MW*T)^0.5*1e-11/(sigma*sigma*omega);
  end gasViscLowPressureChung;

  function gasThCondCorr "Calculates low pressure gas thermal conductivity from a given correlation data and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gThCondCoef directly, coming from an iteration."
  input FreeFluids.MediaCommon.DataRecord data;
  input Real T;
  output Real lambda;
  protected
    Real coef[6]= data.gThCondCoef;
  external "C" FF_PhysPropCorr(data.gThCondCorr, coef, data.MW, T, lambda) annotation(
    IncludeDirectory = "modelica://FreeFluids/Resources",
    Include = "#include \"FFphysprop.c\"");
  end gasThCondCorr;
  
  function FF_gasThCondCorr "Calculates low pressure gas thermal conductivity from a given correlation data and a given temperature. Needed due to bugs in OpenModelica 1.13, that doesn't work if you pass data.gThCondCoef directly, coming from an iteration."
    input FreeFluids.MediaCommon.DataRecord data;
    input Real T;
    output Real lambda;
  algorithm
    lambda:= if data.gThCondCorr==121 then data.gThCondCoef[1]+T*(data.gThCondCoef[2]+T*(data.gThCondCoef[3]+T*(data.gThCondCoef[4]+T*data.gThCondCoef[4]))) else data.gThCondCoef[1]*T^data.gThCondCoef[2]/(1+data.gThCondCoef[3]/T+data.gThCondCoef[4]/T^2);
  end FF_gasThCondCorr;

  function gasThCondLowPressureChung "Thermal conductivity of a low pressure gas according to Chung"
    input FreeFluids.MediaCommon.DataRecord data;
    input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp0;
    input SI.DynamicViscosity eta "low pressure dynamic viscosity";
    input SI.Temp_K T "Gas temperature";
    output Modelica.Media.Interfaces.Types.ThermalConductivity lambda;
  protected
    Real R,alpha,beta,Tr,Z,psi;
  algorithm
    R:=Modelica.Constants.R*1000/data.MW;
    alpha:=(Cp0-R)/R-1.5;
    beta:=0.7862 - 0.7109*data.w + 1.3168 *data.w*data.w;
    Tr:=T/data.Tc;
    Z:=2.0 + 10.5*Tr*Tr;
    psi:=1+alpha*((0.215 + 0.28288 *alpha - 1.061 *beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta ));
    lambda:=3.75*psi*eta*R;
  end gasThCondLowPressureChung;

  function gasMixViscosityWilke
    input SI.Temp_K T;
    input SI.MoleFraction X[:];
    input SI.DynamicViscosity etaX[:];
    input SI.MolarMass MW[:];
    output SI.DynamicViscosity eta;
  protected
    Real phi[size(MW,1),size(MW,1)];
    Real sigma;
  algorithm
    for i in 1:size(X,1) loop
      for j in 1:size(X,1) loop
        phi[i,j]:=(etaX[i]/etaX[j])^0.5*(MW[i]/MW[j])^0.25;
        phi[i,j]:=(1+phi[i,j])^2.0/ (8*(1+MW[i]/MW[j]))^0.5;
      end for;
    end for;
    eta:=0;
    for i in 1:size(X,1) loop
      sigma:=0;
      for j in 1:size(X,1) loop
        sigma:=sigma+X[j]*phi[i,j];
      end for;
      eta:=eta+X[i]*etaX[i]/sigma;
    end for;
  end gasMixViscosityWilke;
  
  function gasMixThCondMason
    input SI.Temp_K T;
    input SI.MoleFraction X[:];
    input SI.ThermalConductivity lambdaX[size(X,1)];
    input SI.MolarMass MW[size(X,1)];
    input SI.Temp_K Tc[size(X,1)];
    input SI.AbsolutePressure Pc[size(X,1)];
    output SI.ThermalConductivity lambda;
  protected
    Real Tr[size(X,1)];
    Real Gamma[size(X,1)];
    Real A[size(X,1),size(X,1)];
    Real ratio;
    Real sigma;
  algorithm
    for i in 1:size(X,1) loop
      Tr[i]:=T/Tc[i];
      Gamma[i]:=210*(Tc[i]*MW[i]*MW[i]*MW[i]/
                           (Pc[i]*Pc[i]*Pc[i]*Pc[i]*1e-20))^0.1666667;
    end for;
    lambda:=0;
    for i in 1:size(X,1) loop
      sigma:=0;
      for j in 1:size(X,1) loop
        ratio:=Gamma[j]*(exp(0.0464*Tr[i])-exp(-0.2412*Tr[i]))/(Gamma[i]*(exp(0.0464*Tr[j])-exp(-0.2412*Tr[j])));
        A[i,j]:=(1+ratio^0.5*(MW[i]/MW[j])^0.25)^2/(8*(1+MW[i]/MW[j]))^0.5;
        sigma:=sigma+X[j]*A[i,j];
      end for;
      lambda:=lambda+X[i]*lambdaX[i]/sigma;
    end for;
  end gasMixThCondMason;

  function liqThCondLatini
    input FreeFluids.MediaCommon.DataRecord data;
    input SI.Temp_K T;
    output SI.ThermalConductivity lambda;
  protected
      Real Tr=T/data.Tc,A,alpha,beta,gamma;
  algorithm
      if data.family == 1 then
          A :=0.00350 "Alkane";
          alpha :=1.2;
          beta :=0.5;
          gamma :=0.167;
      elseif data.family == 2 or data.family == 3 then
          A :=0.0361 "Alkene, Alkyne";
          alpha :=1.2;
          beta :=1.0;
          gamma :=0.167;
      elseif data.family == 4 then
          A :=0.031 "Cycloalkane";
          alpha :=1.2;
          beta :=1.0;
          gamma :=0.167;
      elseif data.family == 5 or data.family == 6 then
          A :=0.0346 "Aromatic, Water";
          alpha :=1.2;
          beta :=1.0;
          gamma :=0.167;
      elseif data.family == 7 or data.family == 8 then
          A :=0.00339 "Alcohol, Polyol";
          alpha :=1.2;
          beta :=0.5;
          gamma :=0.167;
      elseif data.family == 10 then
          A :=0.0385 "Ether";
          alpha :=1.2;
          beta :=1.0;
          gamma :=0.167;
      elseif data.family == 12 then
          A :=0.00383 "Ketone";
          alpha :=1.2;
          beta :=0.5;
          gamma :=0.167;
      elseif data.family == 13 then
          A :=0.00319 "Acid";
          alpha :=1.2;
          beta :=0.5;
          gamma :=0.167;
      elseif data.family == 14 then
          A :=0.0415 "Ester";
          alpha :=1.2;
          beta :=1.0;
          gamma :=0.167;
      else
          A :=0.494;
          alpha :=0.0;
          beta :=0.5;
          gamma :=-0.167;
      end if;
      lambda := A*data.Tb^alpha*(1-Tr)^0.38/(data.MW^beta*data.Tc^gamma*Tr^(1/6));
  end liqThCondLatini;
end Functions;

end MediaCommon;

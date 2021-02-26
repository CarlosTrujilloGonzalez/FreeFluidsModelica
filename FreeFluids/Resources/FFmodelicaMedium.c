/*
 * FFmodelicaMedium.c
 *
 *  Created on: 5/9/2016
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2020  Carlos Trujillo Gonzalez

 *This program is free software; you can redistribute it and/or
 *modify it under the terms of the GNU General Public License version 2
 *as published by the Free Software Foundation

 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.

 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

// contains mainly calculations for pure substances physical properties, by non-eos methods
//Single substance, physical properties correlations
//==================================================


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FFbasic.h"
#include "ModelicaUtilities.h"
#include "FFphysprop.c"

//#define WINDOWS
#if (defined(_WIN32) || defined(__WIN32__) || defined(WINDOWS))
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

int nPoints=1;


//HERE BEGINS THE MODELICA INTERFACE
//==================================

void *FF_createSubstanceData(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP) {
  double T,P,V,answerL[3],answerG[3],H,S;
  char option,state;
  FF_SubstanceData *data = (FF_SubstanceData*) calloc(1,sizeof(FF_SubstanceData));
  if (!data) {
    ModelicaError("Memory allocation error\n");
  }
  //char path[FILENAME_MAX]="C:/Users/Carlos/workspace/FreeFluids/Resources/";
  char path[FILENAME_MAX]="";
  //"C:/Users/CARLOS~1/workspaceNew/ExternalMedium/Resources/"
  //char path[FILENAME_MAX];
  //char *omlib=getenv("OPENMODELICALIBRARY");
  //ModelicaFormatError("%s\n",omlib);
  //strcpy(path,omlib);
  //strcat(path,"/FreeFluids/");
  strcat(path,resDir);
  strcat(path,"/");
  strcat(path,name);
  strcat(path,".sd");
  //printf("%s\n",path);
  //FILE * file= fopen("Z:/Datos/acetone.sd", "rb");
  //FILE * file= fopen("C:/Users/CARLOS~1/workspaceNew/ExternalMedium/acetone.sd", "rb");
  //FILE * file= fopen(name, "rb");
  FILE * file= fopen(path, "rb");
  if (file != NULL) {
    fread(data, sizeof(FF_SubstanceData), 1, file);
    fclose(file);
  }
  else ModelicaError("unable to charge the substance data\n");

  if (thermoModel==1) data->model=FF_CubicType;
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else ModelicaError("thermoModel out of range\n");

  if(!(((thermoModel==1)&&(data->cubicData.eos>0))||((thermoModel==2)&&(data->saftData.eos>0))||((thermoModel==3)&&(data->swData.eos>0)))) ModelicaError("Not valid EOS supplied\n");

  //reference enthalpy and entropy calculation
  if(refState==1){//ASHRAE
      T=233.15;
      FF_VpEosS(T,data,&P);
  }
  else if(refState==2){//IIR
      T=273.15;
      FF_VpEosS(T,data,&P);
  }
  else if(refState==3){//NBP
      P=101325;
      FF_TbEosS(P,data,&T);
  }
  else {//User defined as default
      T=refT;
      P=refP;
  }

  if((refState==1)||(refState==2)||(refState==3)){
      option='l';
      FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
      V=answerL[0];
  }
  else{
      option='s';
      FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
      if(state=='L') V=answerL[0];
      else V=answerG[0];
  }
  data->refT=0.0;//ideal gas calculation from 0K
  data->refP=101325;//ideal gas entropy referenced to 1 atm
  FF_HSfromTVPeosS(T,V,P,data,&data->refH,&data->refS);
  if(refState==2){//IIR
      data->refH=data->refH-2.0e2*data->baseProp.MW;
      data->refS=data->refS-1*data->baseProp.MW;
  }
  return data;
}

void FF_destroySubstanceData(void *object) {
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (data==NULL) return;
  free(data);
}

//Basic functions
//---------------

//Pressure from d,T by EOS
void FF_pressureEOS_dT(void *object, double d, double T, double *P){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double V=data->baseProp.MW/(d*1e3);
  FF_PfromTVeosS(T,V,data,P);
}

//Densities from T,P by EOS. Bellow the critical temperature
void FF_densitiesEOS_pT(FF_SubstanceData *data, double P, double T, int phase, double *ld, double *gd){
  //FF_SubstanceData *data = (FF_SubstanceData *)object;
  if((data->model==FF_SAFTtype)&&(T>(data->saftData.Tc-1.0))) *ld=*gd=data->saftData.Pc*data->saftData.MW/(1000*data->saftData.Zc*8.314472*data->saftData.Tc);
  else if((data->model==FF_SWtype)&&(T>(data->swData.Tc-1.0))) *ld=*gd=data->swData.Pc*data->swData.MW/(1000*data->swData.Zc*8.314472*data->swData.Tc);
  else{
      double answerL[3],answerG[3];//for cubic eos we can arrive to the critical point
      char option,state;
      if (phase==1) option='g';
      else if (phase==0) option='b';
      else if (phase==2) option='l';
      else ModelicaError("Phase requested for density calculation not adequate");
      FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
      if (!(phase==1)) *ld=data->baseProp.MW/answerL[0]/1e3;
      else *ld=0;
      if (!(phase==2)) *gd=data->baseProp.MW/answerG[0]/1e3;
      else *gd=0;
  }
}



//Dynamic viscosity
void FF_dynamicViscosity(FF_SubstanceData *data, double T, double P, double gf, double *eta){
  //FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (gf==0.0){//only liquid phase
    if (data->lViscCorr.id>0) FF_PhysPropCorrM(data->lViscCorr.form,data->lViscCorr.coef,data->baseProp.MW,T,eta);
    else ModelicaError("unable to compute liquid viscosity\n");
    if((P>20e5)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)){
        double Vp,satLiqVisc;
        satLiqVisc= * eta;
        FF_Vp(T,data,&Vp);
        if((P>Vp)&&(Vp>0)) FF_LiqViscPcorLucas(T,P,Vp,&data->baseProp,satLiqVisc,eta);
    }
  }
  else if (gf==1.0){//only gas phase
    double lpGasVisc=0;
    if ((data->gViscCorr.id>0)&&(T>=data->gViscCorr.limI)&&(T<=data->gViscCorr.limS)) FF_PhysPropCorrM(data->gViscCorr.form,data->gViscCorr.coef,data->baseProp.MW,T,&lpGasVisc);
    if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_GasViscTPcpLucas(T,P,&data->baseProp,&lpGasVisc,eta);
    else *eta=lpGasVisc;
    if(*eta==0) ModelicaError("unable to compute gas viscosity\n");
  }
  else ModelicaError("unable to compute two phases viscosity\n");
}



//Thermal conductivity
void FF_thermalConductivity(FF_SubstanceData *data, double T, double P, double gf, double *lambda){
  //FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (gf==0.0){
    if (data->lThCCorr.id>0) FF_PhysPropCorrM(data->lThCCorr.form,data->lThCCorr.coef,data->baseProp.MW,T,lambda);
    else if(data->baseProp.Tc>0)FF_LiquidThCondLatini(T,&data->baseProp,lambda);
    else ModelicaError("unable to compute liquid thermal conductivity\n");
  }
  else if (gf==1.0){
    if ((data->gThCCorr.id>0)&&(T>=data->gThCCorr.limI)&&(T<=data->gThCCorr.limS)&&!(data->gThCCorr.form==122)) FF_PhysPropCorrM(data->gThCCorr.form,data->gThCCorr.coef,data->baseProp.MW,T,lambda);
    else if(data->cp0Corr.form>0){
        double Cp0;
        FF_PhysPropCorrM(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,T,&Cp0);
        FF_GasLpThCondTCpChung(T,Cp0,&data->baseProp,lambda);
    }
    else ModelicaError("unable to compute gas thermal conductivity\n");
  }
  else ModelicaError("unable to compute two phases thermal conductivity\n");  
}

//Surface tension
void FF_surfaceTension(FF_SubstanceData *data, double T, double *sigma){
  //FF_SubstanceData *data = (FF_SubstanceData *)object;
  if ((data->lSurfTCorr.id>0)&&(T>=data->lSurfTCorr.limI)&&(T<=data->lSurfTCorr.limS)) FF_PhysPropCorrM(data->lSurfTCorr.form,data->lSurfTCorr.coef,data->baseProp.MW,T,sigma);
  else if ((data->baseProp.Pa>0)&&(data->lDensCorr.form>0)&&((data->vpCorr.form>0)||((data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)))){
  FF_SurfTensMcLeod(T,data,sigma);
   }
  else if((data->baseProp.Tb>0)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_SurfTensSastri(T,&data->baseProp,sigma);
  else ModelicaError("unable to compute liquid surface tension\n");
}






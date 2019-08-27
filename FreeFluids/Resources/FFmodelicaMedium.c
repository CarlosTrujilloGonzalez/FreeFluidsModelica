/*
 * FFphyspropModelica.c
 *
 *  Created on: 5/9/2016
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2019  Carlos Trujillo Gonzalez

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

#define WINDOWS  /* uncomment this line to use it for windows.*/
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



void *FF_createSubstanceData(const char *name, int thermoModel, double refT, double refP) {
  FF_SubstanceData *data = (FF_SubstanceData*) calloc(1,sizeof(FF_SubstanceData));
  if (!data) {
    ModelicaError("Memory allocation error\n");
  }
  char path[FILENAME_MAX]="C:/Users/Carlos/workspacePH/FreeFluids/Resources/";
  //"C:/Users/CARLOS~1/workspaceNew/ExternalMedium/Resources/"
  //char path[FILENAME_MAX];
  //char *omlib=getenv("OPENMODELICALIBRARY");
  //ModelicaFormatError("%s\n",omlib);
  //strcpy(path,omlib);
  //strcat(path,"/FreeFluids/");
  strcat(path,name);
  strcat(path,".sd");
  //FILE * file= fopen("Z:/Datos/acetone.sd", "rb");
  //FILE * file= fopen("C:/Users/CARLOS~1/workspaceNew/ExternalMedium/acetone.sd", "rb");
  //FILE * file= fopen(name, "rb");
  FILE * file= fopen(path, "rb");
  if (file != NULL) {
    fread(data, sizeof(FF_SubstanceData), 1, file);
    fclose(file);
  }
  else ModelicaError("unable to charge the substance data\n");
  data->model=thermoModel;
  data->refT=refT;
  data->refP=refP;
  return data;
}

void FF_destroySubstanceData(void *object) {
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (data==NULL) return;
  free(data);
}

//Basic functions
//---------------
//Saturation pressure
void FF_saturationPressure(void *object,double T, double *Vp){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  //ModelicaFormatError("saft eos:%i\n",data->saftData.eos);
  int eosType;
  if ((T>=data->baseProp.Tc)&&(data->baseProp.Tc>0)) *Vp=1e15;
  else if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
	FF_VpEOS(&eosType,&T,&data->cubicData,Vp);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
	eosType=FF_SAFTtype;
	FF_VpEOS(&eosType,&T,&data->saftData,Vp);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
	eosType=FF_SWtype;
	FF_VpEOS(&eosType,&T,&data->swData,Vp);
  }
  else ModelicaError("unable to compute vapor pressure with the selected EOS\n");
}

//Saturation temperature
void FF_saturationTemperature(void *object,double P, double *Tb){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  if (P>=data->baseProp.Pc) *Tb=data->baseProp.Tc;
  else if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
	FF_TbEOS(&eosType,&P,&data->cubicData,Tb);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
	eosType=FF_SAFTtype;
	FF_TbEOS(&eosType,&P,&data->saftData,Tb);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
        eosType=FF_SWtype;
	FF_TbEOS(&eosType,&P,&data->swData,Tb);
  }
  else ModelicaError("unable to compute boiling temperature with the selected EOS\n");
}

//Pressure from d,T by vp and lDens correlations
void FF_liquidPressureCorr_dT(void *object, double d, double T, double *P){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if ((T<data->baseProp.Tc)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)){//this grants the density correction by pressure
    double dSat,Vp,Tr,N,Zc;
    if ((data->lDensCorr.form>0)&&(T>=data->lDensCorr.limI)&&(T<=data->lDensCorr.limS)) FF_PhysPropCorr(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,&dSat);
    else if (data->lDens.y>0) FF_LiqDensSatRackett(&data->baseProp,&data->lDens.x,&data->lDens.y,&T,&dSat);
    else ModelicaError("unable to compute the saturated liquid density");
    FF_saturationPressure(data,T,&Vp);
    Tr=T/data->baseProp.Tc;
    N=(1-0.89*data->baseProp.w)*exp(6.9547-76.2853*Tr+191.306*Tr*Tr-203.5472*pow(Tr,3)+82.763*pow(Tr,4));
    if (data->baseProp.Zc==0) Zc=0.29056-0.08775*data->baseProp.w;
    else Zc=data->baseProp.Zc;
    *P=(pow(d/dSat,9)-1)* data->baseProp.Pc/(9*Zc*N)+Vp;
  }
  else ModelicaError("T over Tc, or not enough critical data to compute Chueh-Prausnitz correction");
}

//Pressure from d,T by EOS
void FF_pressureEOS_dT(void *object, double d, double T, double *P){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  double V=data->baseProp.MW/(d*1e3);
  if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
	FF_PfromTVeos(&eosType,&T,&V,&data->cubicData,P);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
	eosType=FF_SAFTtype;
	FF_PfromTVeos(&eosType,&T,&V,&data->saftData,P);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
	eosType=FF_SWtype;
	FF_PfromTVeos(&eosType,&T,&V,&data->swData,P);
  }
  else ModelicaError("unable to compute the pressure from d and T with the selected EOS\n");
}

//Densities from T,P by EOS
void FF_densitiesEOS_pT(void *object, double P, double T, int phase, double *ld, double *gd){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  double answerL[3],answerG[3];
  char option,state;
  if (phase==1) option='g';
  else if (phase==0) option='b';
  else if (phase==2) option='l';
  else ModelicaError("Phase requested for density calculation not adequate");
  if ((data->cubicData.eos>0)&&((data->model==1)||(data->model==11))){
    //añadir un else if para calcular PR in situ si no había datos
	//************************************************************
	eosType=FF_CubicType;
	FF_VfromTPeos(&eosType,&T,&P,&data->cubicData,&option,answerL,answerG,&state);
  }
  else if ((data->saftData.eos>0)&&((data->model==2)||(data->model==12))){
	eosType=FF_SAFTtype;
	FF_VfromTPeos(&eosType,&T,&P,&data->saftData,&option,answerL,answerG,&state);
  }
  else if ((data->swData.eos>0)&&((data->model==3)||(data->model==13))){
	eosType=FF_SWtype;
	FF_VfromTPeos(&eosType,&T,&P,&data->swData,&option,answerL,answerG,&state);
  }
  else ModelicaError("unable to compute the densities with the given EOS");
  if (answerL[0]>0) *ld=data->baseProp.MW/answerL[0]/1e3;
  else *ld=0;
  if (answerG[0]>0) *gd=data->baseProp.MW/answerG[0]/1e3;
  else *gd=0;
}


//Bubble enthalpy, from a SaturationProperties record (p and T)
void FF_bubbleEnthalpy(void *object, double P, double T,double *lh){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double ld,gd;
  FF_ThermoProperties th;
  int eosType;
  FF_densitiesEOS_pT(data,P,T,2,&ld,&gd);//we get the liquid density by EOS
  th.T=T;
  th.P=P;
  th.V=data->baseProp.MW/ld/1e3;
  if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
    FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
    eosType=FF_SAFTtype;
	FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
	eosType=FF_SWtype;
	FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  *lh=th.H/data->baseProp.MW*1e3;
}

//Dew enthalpy from a SaturationProperties record (p and T)
void FF_dewEnthalpy(void *object, double P, double T,double *gh){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double ld,gd,lh,hv;
  FF_ThermoProperties th;
  int eosType;
    FF_densitiesEOS_pT(data,P,T,1,&ld,&gd);//we get the gas density
    th.T=T;
    th.P=P;
    th.V=data->baseProp.MW/gd/1e3;
    if ((data->cubicData.eos>0)&&(data->model==1)){
          eosType=FF_CubicType;
          FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
    }
    else if ((data->saftData.eos>0)&&(data->model==2)){
          eosType=FF_SAFTtype;
          FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
    }
    else if ((data->swData.eos>0)&&(data->model==3)){
          eosType=FF_SWtype;
          FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
    }
    *gh=th.H/data->baseProp.MW*1e3;

}

//Thermodynamic properties from a thermodynamic record of T and densities
void FF_thermoPropertiesEOS_Td(void *object, double T, double ld, double gd, double gf, double *h, double *u, double *s, double *cp, double *cv, double *dP_dT, double *dP_dV, double *ss, double *jt, double *it){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  FF_ThermoProperties thl,thg;
  if (gf<1){
      thl.T=T;
      thl.V=data->baseProp.MW/ld/1e3;
      if ((data->cubicData.eos>0)&&(data->model==1)){
        eosType=FF_CubicType;
        FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
      else if ((data->saftData.eos>0)&&(data->model==2)){
        eosType=FF_SAFTtype;
        FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
      else if ((data->swData.eos>0)&&(data->model==3)){
        eosType=FF_SWtype;
        FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
  }
  if (gf>0){
      thg.T=T;
      thg.V=data->baseProp.MW/gd/1e3;
      if ((data->cubicData.eos>0)&&(data->model==1)){
            eosType=FF_CubicType;
            FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
      else if ((data->saftData.eos>0)&&(data->model==2)){
            eosType=FF_SAFTtype;
            FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
      else if ((data->swData.eos>0)&&(data->model==3)){
            eosType=FF_SWtype;
            FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
  }
  *h=(thl.H*(1-gf)+thg.H*gf)*1e3/data->baseProp.MW;
  *u=(thl.U*(1-gf)+thg.U*gf)*1e3/data->baseProp.MW;
  *s=(thl.S*(1-gf)+thg.S*gf)*1e3/data->baseProp.MW;
  if (gf==0.0){
	*cp=thl.Cp*1e3/data->baseProp.MW;
	*cv=thl.Cv*1e3/data->baseProp.MW;
    *ss=thl.SS;
    *dP_dV=thl.dP_dV*data->baseProp.MW/1e3;
	*dP_dT=thl.dP_dT;
  }
  else if (gf==1.0){
	*cp=thg.Cp*1e3/data->baseProp.MW;
	*cv=thg.Cv*1e3/data->baseProp.MW;
    *ss=thg.SS;
    *dP_dV=thg.dP_dV*data->baseProp.MW/1e3;
	*dP_dT=thg.dP_dT;
  }
  else{
	*cp=0.0;
	*cv=0.0;
	*ss=0.0;
	*dP_dV=0.0;
	*dP_dT=0.0;
  }
}

//Specific enthalpy from a thermodynamic record
void FF_specificEnthalpy(void *object, double T, double ld, double gd, double gf, double *h){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  double lh,gh;
  FF_ThermoProperties thl,thg;
  if (gf<1){
      thl.T=T;
      thl.V=data->baseProp.MW/ld/1e3;
      if ((data->cubicData.eos>0)&&(data->model==1)){
        eosType=FF_CubicType;
        FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
      else if ((data->saftData.eos>0)&&(data->model==2)){
        eosType=FF_SAFTtype;
        FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
      else if ((data->swData.eos>0)&&(data->model==3)){
        eosType=FF_SWtype;
        FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
      }
      lh=thl.H/data->baseProp.MW*1e3;
  }
  if (gf>0){
      thg.T=T;
      thg.V=data->baseProp.MW/gd/1e3;
      if ((data->cubicData.eos>0)&&(data->model==1)){
            eosType=FF_CubicType;
            FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
      else if ((data->saftData.eos>0)&&(data->model==2)){
            eosType=FF_SAFTtype;
            FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
      else if ((data->swData.eos>0)&&(data->model==3)){
            eosType=FF_SWtype;
            FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
      }
      gh=thg.H/data->baseProp.MW*1e3;
  }
  *h=lh*(1-gf)+gh*gf;
  }

//Specific heat capacity from a thermodynamic state record
void FF_specificHeatCapacityCp(void *object, double P, double T, double ld, double gd, double gf, double *cp){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  FF_ThermoProperties thl,thg;
  double cpl;
  if (gf<1){
    thl.T=T;
    thl.P=P;
    thl.V=data->baseProp.MW/ld/1e3;
    if ((data->cubicData.eos>0)&&(data->model==1)){
      eosType=FF_CubicType;
      FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
    else if ((data->saftData.eos>0)&&(data->model==2)){
      eosType=FF_SAFTtype;
      FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
    else if ((data->swData.eos>0)&&(data->model==3)){
      eosType=FF_SWtype;
      FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
    cpl=thl.Cp*1e3/data->baseProp.MW;
  }
  if (gf>0){
    thg.T=T;
    thg.P=P;
    thg.V=data->baseProp.MW/gd/1e3; 
    if ((data->cubicData.eos>0)&&(data->model==1)){
	  eosType=FF_CubicType;
	  FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
    else if ((data->saftData.eos>0)&&(data->model==2)){
	  eosType=FF_SAFTtype;
	  FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
    else if ((data->swData.eos>0)&&(data->model==3)){
	  eosType=FF_SWtype;
	  FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
  }
  *cp=cpl*(1-gf)+(thg.Cp*gf)*1e3/data->baseProp.MW;
}

//Bubble entropy, from a SaturationProperties record (p and T)
void FF_bubbleEntropy(void *object, double P, double T,double *ls){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double ld,gd;
  FF_ThermoProperties th;
  int eosType;
  FF_densitiesEOS_pT(data,P,T,2,&ld,&gd);
  th.T=T;
  th.P=P;
  th.V=data->baseProp.MW/ld/1e3;
  if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
	FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
	eosType=FF_SAFTtype;
	FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
	eosType=FF_SWtype;
	FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  *ls=th.S/data->baseProp.MW*1e3;
}

//Dew entropy from a SaturationProperties record (p and T)
void FF_dewEntropy(void *object, double P, double T,double *gs){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double ld,gd;
  FF_ThermoProperties th;
  int eosType;
  FF_densitiesEOS_pT(data,P,T,1,&ld,&gd);
  th.T=T;
  th.P=P;
  th.V=data->baseProp.MW/gd/1e3;
  if ((data->cubicData.eos>0)&&(data->model==1)){
	eosType=FF_CubicType;
	FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->saftData.eos>0)&&(data->model==2)){
	eosType=FF_SAFTtype;
	FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  else if ((data->swData.eos>0)&&(data->model==3)){
	eosType=FF_SWtype;
	FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &th);
  }
  *gs=th.S/data->baseProp.MW*1e3;
}

//Spefific entropy
void FF_specificEntropy(void *object, double P, double T, double ld, double gd, double gf, double *s){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  int eosType;
  FF_ThermoProperties thl,thg;
  if (gf<1){
    thl.T=T;
    thl.P=P;
    thl.V=data->baseProp.MW/ld/1e3; 
    if ((data->cubicData.eos>0)&&(data->model==1)){
	  eosType=FF_CubicType;
	  FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
    else if ((data->saftData.eos>0)&&(data->model==2)){
	  eosType=FF_SAFTtype;
	  FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
    else if ((data->swData.eos>0)&&(data->model==3)){
	  eosType=FF_SWtype;
	  FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thl);
    }
  }
  if (gf>0){
    thg.T=T;
    thg.P=P;
    thg.V=data->baseProp.MW/gd/1e3; 
    if ((data->cubicData.eos>0)&&((data->model==1)||(data->model==11)||(data->model==21))){
	  eosType=FF_CubicType;
	  FF_ThermoEOS(&eosType,&data->cubicData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
    else if ((data->saftData.eos>0)&&((data->model==2)||(data->model==12)||(data->model==22))){
	  eosType=FF_SAFTtype;
	  FF_ThermoEOS(&eosType,&data->saftData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
    else if ((data->swData.eos>0)&&((data->model==3)||(data->model==13)||(data->model==23))){
	  eosType=FF_SWtype;
          FF_ThermoEOS(&eosType,&data->swData,&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP, &thg);
    }
  }
  *s=(thl.S*(1-gf)+thg.S*gf)/data->baseProp.MW*1e3;  
}


//Dynamic viscosity
void FF_dynamicViscosity(void *object, double T, double P, double gf, double *eta){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (gf==0.0){//only liquid phase
    if (data->lViscCorr.id>0) FF_PhysPropCorr(data->lViscCorr.form,data->lViscCorr.coef,data->baseProp.MW,T,eta);
    else ModelicaError("unable to compute liquid viscosity\n");
    if((P>20e5)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)){
        double Vp,satLiqVisc;
        satLiqVisc= * eta;
        FF_Vp(T,data,&Vp);
        if((P>Vp)&&(Vp>0)) FF_LiqViscPcorLucas(&T,&P,&Vp,&data->baseProp,&satLiqVisc,eta);
    }
  }
  else if (gf==1.0){//only gas phase
    double lpGasVisc=0;
    if ((data->gViscCorr.id>0)&&(T>=data->gViscCorr.limI)&&(T<=data->gViscCorr.limS)) FF_PhysPropCorr(data->gViscCorr.form,data->gViscCorr.coef,data->baseProp.MW,T,&lpGasVisc);
    if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_GasViscTPcpLucas(&T,&P,&data->baseProp,&lpGasVisc,eta);
    else *eta=lpGasVisc;
    if(*eta==0) ModelicaError("unable to compute gas viscosity\n");
  }
  else ModelicaError("unable to compute two phases viscosity\n");
}

//Thermal conductivity
void FF_thermalConductivity(void *object, double T, double P, double gf, double *lambda){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (gf==0.0){
    if (data->lThCCorr.id>0) FF_PhysPropCorr(data->lThCCorr.form,data->lThCCorr.coef,data->baseProp.MW,T,lambda);
    else if(data->baseProp.Tc>0)FF_LiquidThCondLatini(&T,&data->baseProp,lambda);
    else ModelicaError("unable to compute liquid thermal conductivity\n");
  }
  else if (gf==1.0){
    if ((data->gThCCorr.id>0)&&(T>=data->gThCCorr.limI)&&(T<=data->gThCCorr.limS)) FF_PhysPropCorr(data->gThCCorr.form,data->gThCCorr.coef,data->baseProp.MW,T,lambda);
    else ModelicaError("unable to compute gas thermal conductivity\n");
  }
  else ModelicaError("unable to compute two phases thermal conductivity\n");  
}

//Surface tension
void FF_surfaceTension(void *object, double T, double P, double gf, double *sigma){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if (gf<1.0){
        if ((data->lSurfTCorr.id>0)&&(T>=data->lSurfTCorr.limI)&&(T<=data->lSurfTCorr.limS)) FF_PhysPropCorr(data->lSurfTCorr.form,data->lSurfTCorr.coef,data->baseProp.MW,T,sigma);
        else if ((data->baseProp.Pa>0)&&(data->lDensCorr.form>0)&&((data->vpCorr.form>0)||((data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)))){
            FF_SurfTensMcLeod(T,data,sigma);
        }
        else if((data->baseProp.Tb>0)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_SurfTensSastri(&T,&data->baseProp,sigma);
        else ModelicaError("unable to compute liquid surface tension\n");
  }
  else ModelicaError("unable to compute surface tension without liquid\n");  
}






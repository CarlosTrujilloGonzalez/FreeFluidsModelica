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



void *FF_createSubstanceData(const char *name, int thermoModel, int refState, double refT, double refP) {
  double T,P,V,answerL[3],answerG[3],H,S;
  char option,state;
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

  if (thermoModel==1) data->model=FF_CubicType;
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else ModelicaError("thermoModel out of range\n");

  if(!(((thermoModel==1)&&(data->cubicData.eos>0))||((thermoModel==2)&&(data->saftData.eos>0))||((thermoModel==3)&&(data->swData.eos>0)))) ModelicaError("Not valid EOS supplied\n");

  //reference enthalpy and entropy calculation
  if(refState==1){//ASHRAE
      T=233.15;
      FF_VpEOSs(&T,data,&P);
  }
  else if(refState==2){//IIR
      T=273.15;
      FF_VpEOSs(&T,data,&P);
  }
  else if(refState==3){//NBP
      P=101325;
      FF_TbEOSs(&P,data,&T);
  }
  else {//User defined as default
      T=refT;
      P=refP;
  }

  if((refState==1)||(refState==2)||(refState==3)){
      option='l';
      FF_VfromTPeosS(&T,&P,data,&option,answerL,answerG,&state);
      V=answerL[0];
  }
  else{
      option='s';
      FF_VfromTPeosS(&T,&P,data,&option,answerL,answerG,&state);
      if(state=='L') V=answerL[0];
      else V=answerG[0];
  }
  data->refT=0.0;//ideal gas calculation from 0K
  data->refP=101325;//ideal gas entropy referenced to 1 atm
  FF_HSfromTVPeosS(&T,&V,&P,data,&data->refH,&data->refS);
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

typedef struct{double a,av,avv,at,att,avt;} FF_HelmholtzDerivatives;
typedef struct{double cp,h,s;} FF_IdealThermoProp;


//Basic functions
//---------------
//Saturation pressure
void FF_saturationPressure(void *object,double T, double *Vp){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if((data->model==FF_CubicType)&&(data->cubicData.Tc>0)&&(T>=data->cubicData.Tc)) *Vp=1e15;
  else if((data->model==FF_SAFTtype)&&(data->saftData.Tc>0)&&(T>=data->saftData.Tc)) *Vp=1e15;
  else if((data->model==FF_SWtype)&&(data->swData.Tc>0)&&(T>=data->swData.Tc)) *Vp=1e15;
  else if((data->model==FF_SWtype)&&(data->vpCorr.form>0)) FF_PhysPropCorr(data->vpCorr.form,data->vpCorr.coef,data->swData.MW,T,Vp);
  else FF_VpEOSs(&T,data,Vp);
}


//Saturation temperature
void FF_saturationTemperature(void *object,double P, double *Tb){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  if((data->model==FF_CubicType)&&(data->cubicData.Pc>0)&&(P>=data->cubicData.Pc)) *Tb=data->cubicData.Tc;
  else if((data->model==FF_SAFTtype)&&(data->saftData.Pc>0)&&(P>=data->saftData.Pc)) *Tb=data->saftData.Tc;
  else if((data->model==FF_SWtype)&&(data->swData.Pc>0)&&(P>=data->swData.Pc)) *Tb=data->swData.Tc;
  else if((data->model==FF_SWtype)&&(data->btCorr.form>0)) FF_PhysPropCorr(data->btCorr.form,data->btCorr.coef,data->swData.MW,P,Tb);
  FF_TbEOSs(&P,data,Tb);
}


//Pressure from d,T by EOS
void FF_pressureEOS_dT(void *object, double d, double T, double *P){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double V=data->baseProp.MW/(d*1e3);
  FF_PfromTVeosS(&T,&V,data,P);
}

//Densities from T,P by EOS
void FF_densitiesEOS_pT(void *object, double P, double T, int phase, double *ld, double *gd){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double answerL[3],answerG[3];
  char option,state;
  if (phase==1) option='g';
  else if (phase==0) option='b';
  else if (phase==2) option='l';
  else ModelicaError("Phase requested for density calculation not adequate");
  FF_VfromTPeosS(&T,&P,data,&option,answerL,answerG,&state);
  if (!(phase==1)) *ld=data->baseProp.MW/answerL[0]/1e3;
  else *ld=0;
  if (!(phase==2)) *gd=data->baseProp.MW/answerG[0]/1e3;
  else *gd=0;
}

void FF_satProp_T(void *object, double T, double p, int phase, double *nMols, double *d, double deriv[6], double ideal[3]){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double answerL[3],answerG[3];
    char option,state;
    FF_CubicParam param;
    FF_ThermoProperties th0;
    *nMols=1/MW;
    if (phase==1){
        if((data->model==FF_SWtype)&&(data->gDensCorr.form>0)){
            FF_PhysPropCorr(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,d);
            th0.V=MW/ *d;
        }
        else{
            option='g';
            FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
            th0.V=answerG[0];
            *d=MW/answerG[0];
        }
    }
    else if(phase==2){
        if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)){
            FF_PhysPropCorr(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,d);
            th0.V=MW/ *d;
        }
        else{
            option='l';
            FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
            th0.V=answerL[0];
            *d=MW/answerL[0];
        }
    }
    th0.T=T;
    if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
    else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
    ideal[0]=th0.Cp;
    ideal[1]=th0.H-data->refH;
    ideal[2]=th0.S-data->refS;
    if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&th0.V,&data->swData,deriv);
    else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&th0.V,&data->saftData,deriv);
    else{
        FF_FixedParamCubic(&data->cubicData,&param);
        FF_ThetaDerivCubic(&T,&data->cubicData,&param);
        FF_ArrDerCubic(&T,&th0.V,&param,deriv);
    }
}

void FF_solveEOS_TpOrig(void *object, double T, double p, double *nMols, double *ld, double *gd, double lDer[6], double gDer[6], double lIdeal[3], double gIdeal[3]){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double answerL[3],answerG[3];
    char option,state;
    FF_CubicParam param;
    FF_ThermoProperties th0;
    *nMols=1/MW;
    //Volume calculation
    if (T>=data->baseProp.Tc) option='g';
    else if (p>data->baseProp.Pc) option='l';
    else option='s';
    FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
    th0.T=T;
    if(state=='L'){
        th0.V=answerL[0];
        *ld=MW/answerL[0];
        *gd=0;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        lIdeal[0]=th0.Cp;
        lIdeal[1]=th0.H-data->refH;
        lIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&answerL[0],&data->swData,lDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&answerL[0],&data->saftData,lDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&answerL[0],&param,lDer);
        }
    }
    else if(state=='G'){
        th0.V=answerG[0];
        *gd=MW/answerG[0];
        *ld=0;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        gIdeal[0]=th0.Cp;
        gIdeal[1]=th0.H-data->refH;
        gIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&answerG[0],&data->swData,gDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&answerG[0],&data->saftData,gDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&answerG[0],&param,gDer);
        }
    }
}

void FF_solveEOS_Tp(void *object, double T, double p, double *nMols, double *ld, double *gd, void *lDeriva, void *gDeriva, void *lId, void *gId){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    FF_HelmholtzDerivatives *lDer=(FF_HelmholtzDerivatives *)lDeriva;
    FF_HelmholtzDerivatives *gDer=(FF_HelmholtzDerivatives *)gDeriva;
    FF_IdealThermoProp *lIdeal=(FF_IdealThermoProp *)lId;
    FF_IdealThermoProp *gIdeal=(FF_IdealThermoProp *)gId;
    /*
    double lDeriv[6],gDeriv[6];
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double answerL[3],answerG[3];
    char option,state;
    FF_CubicParam param;
    FF_ThermoProperties th0;
    *nMols=1/MW;
    //Volume calculation
    if (T>=data->baseProp.Tc) option='g';
    else if (p>data->baseProp.Pc) option='l';
    else option='s';
    FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
    th0.T=T;
    if(state=='L'){
        th0.V=answerL[0];
        *ld=MW/answerL[0];
        *gd=0;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        lIdeal->cp=th0.Cp;
        lIdeal->h=th0.H-data->refH;
        lIdeal->s=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&answerL[0],&data->swData,lDeriv);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&answerL[0],&data->saftData,lDeriv);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&answerL[0],&param,lDeriv);
        }
        lDer->a=lDeriv[0];
        lDer->av=lDeriv[1];
        lDer->avv=lDeriv[2];
        lDer->at=lDeriv[3];
        lDer->att=lDeriv[4];
        lDer->avt=lDeriv[5];
    }
    else if(state=='G'){
        th0.V=answerG[0];
        *gd=MW/answerG[0];
        *ld=0;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        gIdeal->cp=th0.Cp;
        gIdeal->cp=th0.H-data->refH;
        gIdeal->cp=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&answerG[0],&data->swData,gDeriv);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&answerG[0],&data->saftData,gDeriv);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&answerG[0],&param,gDeriv);
        }
        gDer->a=gDeriv[0];
        gDer->av=gDeriv[1];
        gDer->avv=gDeriv[2];
        gDer->at=gDeriv[3];
        gDer->att=gDeriv[4];
        gDer->avt=gDeriv[5];
    }*/
}

void FF_solveEOS_Td(void *object, double T, double d, double *nMols, double *p, double *ld, double *gd, double lDer[6], double gDer[6], double lIdeal[3], double gIdeal[3]){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double V;
    double Vp;
    double answerL[3],answerG[3];
    char option,state;
    FF_CubicParam param;
    FF_ThermoProperties th0;
    *nMols=1/MW;
    V=MW/d;
    //Volume calculation
    if (T>data->baseProp.Tc) state='G';
    else{
        FF_VpEOSs(&T,data,&Vp);
        option='b';
        FF_VfromTPeosS(&T,&Vp,data,&option,answerL,answerG,&state);
        if(V<=answerL[0]) state='L';
        else if (V>=answerG[0]) state='G';
        else state='B';
    }
    if (state=='L'){
        *ld=d;
        *gd=0.0;
        FF_PfromTVeosS(&T,&V,data,p);
    }
    else if (state=='G'){
        *gd=d;
        *ld=0.0;
        FF_PfromTVeosS(&T,&V,data,p);
    }
    else{
        *ld=MW/answerL[0];
        *gd=MW/answerG[0];
        *p=Vp;
    }

    th0.T=T;
    if(!(state=='G')){
        th0.V=MW/ *ld;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        lIdeal[0]=th0.Cp;
        lIdeal[1]=th0.H-data->refH;
        lIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&th0.V,&data->swData,lDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&th0.V,&data->saftData,lDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&th0.V,&param,lDer);
            }
    }
    if(!(state=='L')){
        th0.V=MW/ *gd;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        gIdeal[0]=th0.Cp;
        gIdeal[1]=th0.H-data->refH;
        gIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(&T,&th0.V,&data->swData,gDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(&T,&th0.V,&data->saftData,gDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(&T,&data->cubicData,&param);
            FF_ArrDerCubic(&T,&th0.V,&param,gDer);
        }
    }
}

void FF_solveEOS_ph(void *object, double p, double h, double *nMols, double *T, double *ld, double *gd, double *gf, double lDer[6], double gDer[6], double lIdeal[3], double gIdeal[3]){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double Tb;
    double Hraw,Hl,Hh,Sl,Sh,Ht,Cp0;
    double answerL[3],answerG[3];
    char option,state;
    FF_CubicParam param;
    FF_ThermoProperties th0;
    *nMols=1/MW;
    Hraw=h*MW + data->refH;
    th0.P=p;
    if(p>data->baseProp.Pc) Tb=data->baseProp.Tc;
    else if((data->model==FF_SWtype)&&(data->btCorr.form>0)) FF_PhysPropCorr(data->btCorr.form,data->btCorr.coef,data->baseProp.MW,p,&Tb);
    else FF_TbEOSs(&p,data,&Tb);
    option='b';
    if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)&&(data->gDensCorr.form>0)){//if we can use the saturated densities correlations
        FF_PhysPropCorr(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,Tb,ld);
        FF_PhysPropCorr(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,Tb,gd);
        answerL[0]=MW/ *ld;
        answerG[0]=MW/ *gd;
    }
    else  FF_VfromTPeosS(&Tb,&p,data,&option,answerL,answerG,&state);
    FF_HSfromTVPeosS(&Tb,&answerL[0],&p,data,&Hl,&Sl);
    FF_HSfromTVPeosS(&Tb,&answerG[0],&p,data,&Hh,&Sh);
    //printf("Tb:%f Hraw:%f Hh:%f Hl:%f\n",Tb,Hraw,Hh,Hl);

    if((Hraw>=Hl)&&(Hraw<=Hh)){
        *T=Tb;
        *ld=MW/answerL[0];
        *gd=MW/answerG[0];
        *gf=(Hraw-Hl)/(Hh-Hl);
    }

    else if (Hraw>Hh){
        *gf=1.0;
        *ld=0.0;
        option='g';
        while (fabs((Hraw-Hh)/Hraw)>0.0001){
            FF_PhysPropCorr(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,Tb,&Cp0);
            Tb=Tb+(Hraw-Hh)/(1.5*MW*Cp0);
            FF_VfromTPeosS(&Tb,&p,data,&option,answerL,answerG,&state);
            FF_HSfromTVPeosS(&Tb,&answerG[0],&p,data,&Hh,&Sh);
            //printf("Tb:%f Hraw:%f, Hh:%f Cp0:%f\n",Tb,Hraw,Hh,Cp0);
        }
        *T=Tb;
        *gd=MW/answerG[0];
    }
    else{
        *gf=0.0;
        *gd=0.0;
        option='l';
        while (fabs((Hraw-Hl)/Hraw)>0.0001){
            FF_PhysPropCorr(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,Tb,&Cp0);
            Tb=Tb+(Hraw-Hl)/(1.5*MW*Cp0);
            FF_VfromTPeosS(&Tb,&p,data,&option,answerL,answerG,&state);
            FF_HSfromTVPeosS(&Tb,&answerL[0],&p,data,&Hl,&Sl);
            //printf("Tb:%f Hraw:%f, Hl:%f Cp0:%f\n",Tb,Hraw,Hl,Cp0);
        }
        *T=Tb;
        *ld=MW/answerL[0];
    }
    th0.T=*T;
    if(!(*gf==1.0)){
        th0.V=MW/ *ld;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        lIdeal[0]=th0.Cp;
        lIdeal[1]=th0.H-data->refH;
        lIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(T,&th0.V,&data->swData,lDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(T,&th0.V,&data->saftData,lDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(T,&data->cubicData,&param);
            FF_ArrDerCubic(T,&th0.V,&param,lDer);
            }
    }
    if(!(*gf==0.0)){
        th0.V=MW/ *gd;
        if((data->model==FF_SWtype)&&(data->swData.eos==FF_IAPWS95)) FF_IdealThermoWater(&th0);
        else FF_IdealThermoEOS(&data->cp0Corr.form,data->cp0Corr.coef,&data->refT,&data->refP,&th0);
        gIdeal[0]=th0.Cp;
        gIdeal[1]=th0.H-data->refH;
        gIdeal[2]=th0.S-data->refS;

        if(data->model==FF_SWtype) FF_ArrDerSWTV(T,&th0.V,&data->swData,gDer);
        else if(data->model==FF_SAFTtype) FF_ArrDerSAFT(T,&th0.V,&data->saftData,gDer);
        else{
            FF_FixedParamCubic(&data->cubicData,&param);
            FF_ThetaDerivCubic(T,&data->cubicData,&param);
            FF_ArrDerCubic(T,&th0.V,&param,gDer);
        }
    }


}

void FF_hsEOS_Tdp(void *object, double T, double ld, double gd, double gf, double p, double *h, double *s){
    FF_SubstanceData *data = (FF_SubstanceData *)object;
    double Vl,Vg,Hl,Hg,Sl,Sg;
    if (gf<1.0){
        Vl=data->baseProp.MW/ld/1e3;
        FF_HSfromTVPeosS(&T,&Vl,&p,data,&Hl,&Sl);
    }
    if (gf>0.0){
        Vg=data->baseProp.MW/gd/1e3;
        FF_HSfromTVPeosS(&T,&Vg,&p,data,&Hg,&Sg);
    }
    *h=((Hl*(1-gf)+Hg*gf))*1e3/data->baseProp.MW-data->refH;
    *s=((Sl*(1-gf)+Sg*gf))*1e3/data->baseProp.MW-data->refS;
}

//Thermodynamic properties from a thermodynamic record of T and densities
void FF_thermoPropertiesEOS_Td(void *object, double T, double ld, double gd, double gf, double *h, double *u, double *s, double *cp, double *cv, double *dP_dT, double *dP_dV, double *ss, double *jt, double *it){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  FF_ThermoProperties thl,thg;
  if (gf<1.0){
      thl.T=T;
      thl.V=data->baseProp.MW/ld/1e3;
      FF_ThermoEOSs(data,&thl);
  }
  if (gf>0.0){
      thg.T=T;
      thg.V=data->baseProp.MW/gd/1e3;
      FF_ThermoEOSs(data,&thg);
  }
  *h=((thl.H*(1-gf)+thg.H*gf))*1e3/data->baseProp.MW-data->refH;
  *u=*h-thl.P*((1-gf)*thl.V+gf*thg.V);
  //*u=(thl.U*(1-gf)+thg.U*gf)*1e3/data->baseProp.MW;
  *s=((thl.S*(1-gf)+thg.S*gf))*1e3/data->baseProp.MW-data->refS;
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


//Bubble enthalpy, from a SaturationProperties record (p and T)
void FF_bubbleEnthalpy(void *object, double p, double T,double *lh){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double MW;
  char option,state;
  double answerL[3],answerG[3];
  double ls;
  MW=data->baseProp.MW*0.001;
  option='l';
  FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
  FF_HSfromTVPeosS(&T,&answerL[0],&p,data,lh,&ls);
  *lh=(*lh-data->refH)/MW;
}

//Dew enthalpy from a SaturationProperties record (p and T)
void FF_dewEnthalpy(void *object, double p, double T,double *gh){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double MW;
  char option,state;
  double answerL[3],answerG[3];
  double gs;
  MW=data->baseProp.MW*0.001;
  option='g';
  FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
  FF_HSfromTVPeosS(&T,&answerG[0],&p,data,gh,&gs);
  *gh=(*gh-data->refH)/MW;
}


//Specific enthalpy from T and d
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
void FF_bubbleEntropy(void *object, double p, double T,double *ls){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double MW;
  char option,state;
  double answerL[3],answerG[3];
  double lh;
  MW=data->baseProp.MW*0.001;
  option='l';
  FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
  FF_HSfromTVPeosS(&T,&answerL[0],&p,data,&lh,ls);
  *ls=(*ls-data->refS)/MW;
}

//Dew entropy from a SaturationProperties record (p and T)
void FF_dewEntropy(void *object, double p, double T,double *gs){
  FF_SubstanceData *data = (FF_SubstanceData *)object;
  double MW;
  char option,state;
  double answerL[3],answerG[3];
  double gh;
  MW=data->baseProp.MW*0.001;
  option='g';
  FF_VfromTPeosS(&T,&p,data,&option,answerL,answerG,&state);
  FF_HSfromTVPeosS(&T,&answerG[0],&p,data,&gh,gs);
  *gs=(*gs-data->refS)/MW;
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

typedef struct {double a,b,c,d;} FF_RecordProvaExternal;

void FF_FuncioProvaExternal(void *st, void *st2){
   FF_RecordProvaExternal *rec=(FF_RecordProvaExternal*) st;
   FF_RecordProvaExternal *rec2=(FF_RecordProvaExternal*) st2;
   rec->a=1;
   rec->b=2;
   rec->c=3;
   rec->d=4;
   rec2->a=5;
   rec2->b=6;
   rec2->c=7;
   rec2->d=8;
}

void FF_FuncioProvaExternal2(void *lDeriva, void *gDeriva, void *lId, void *gId){
    FF_HelmholtzDerivatives *lDer=(FF_HelmholtzDerivatives *)lDeriva;
    FF_HelmholtzDerivatives *gDer=(FF_HelmholtzDerivatives *)gDeriva;
    FF_IdealThermoProp *lIdeal=(FF_IdealThermoProp *)lId;
    FF_IdealThermoProp *gIdeal=(FF_IdealThermoProp *)gId;
    lDer->a=1.0;
    gDer->a=2.0;
    lIdeal->h=3.0;
    gIdeal->h=4.0;
}







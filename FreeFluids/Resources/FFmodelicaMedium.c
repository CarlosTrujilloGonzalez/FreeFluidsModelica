/*
 * FFmodelicaMedium.c
 *
 *  Created on: 5/9/2016
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2023  Carlos Trujillo Gonzalez

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
//#include "ModelicaUtilities.h"
#include "FFeosPure.c"
#include "FFphysprop.c"
//#include "FFeosMix.h"

//#define WINDOWS
#if (defined(_WIN32) || defined(__WIN32__) || defined(WINDOWS))
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

int nPoints=1;


//HERE BEGINS THE MODELICA INTERFACE FOR PURE SUBSTANCES
//======================================================

//Creates a static substances array that is charged each time it is called
void *FF_createSubstanceData(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP) {
    double T,P,V,answerL[3],answerG[3],H,S;
    char option,state;
    int i;
    static int n;//defined substances counter
    static char substance[10][30];
    static FF_SubstanceData data[10];
    for (i=0;i<n;i++){
        if (strcmp(substance[i],name)==0){
            //printf("Found substance:%i \n",i);
            break;
        }
    }
    if (i==n){
        //i=n;
        n++;//we increment the number of defined substances
        strcpy(substance[i],name);
        //printf("New substance:%s %i\n",substance[i],i);
        char path[FILENAME_MAX]="";
        strcat(path,resDir);
        strcat(path,"/Fluids/");
        strcat(path,name);
        strcat(path,".sd");
        //printf("%s\n",path);
        FILE * file= fopen(path, "rb");
        if (file != NULL) {
          fread(&data[i], sizeof(FF_SubstanceData), 1, file);
          fclose(file);
        }
        else printf("unable to charge the substance data\n");

        if (thermoModel==1) data[i].model=FF_CubicType;
        else if (thermoModel==2) data[i].model=FF_SAFTtype;
        else if (thermoModel==3) data[i].model=FF_SWtype;
        else printf("thermoModel out of range\n");

        if(!(((thermoModel==1)&&(data[i].cubicData.eos>0))||((thermoModel==2)&&(data[i].saftData.eos>0))||((thermoModel==3)&&(data[i].swData.eos>0)))) printf("Not valid EOS supplied\n");

        //reference enthalpy and entropy calculation
        if(refState==1){//ASHRAE
            T=233.15;
            FF_VpEosS(T,&data[i],&P);
        }
        else if(refState==2){//IIR
            T=273.15;
            FF_VpEosS(T,&data[i],&P);
        }
        else if(refState==3){//NBP
            P=101325;
            FF_TbEosS(P,&data[i],&T);
        }
        else {//User defined as default
            T=refT;
            P=refP;
        }

        if((refState==1)||(refState==2)||(refState==3)){
            option='l';
            FF_VfromTPeosS(T,P,&data[i],option,answerL,answerG,&state);
            V=answerL[0];
        }
        else{
            option='s';
            FF_VfromTPeosS(T,P,&data[i],option,answerL,answerG,&state);
            if(state=='L') V=answerL[0];
            else V=answerG[0];
        }
        data[i].refT=0.0;//ideal gas calculation from 0K
        data[i].refP=101325;//ideal gas entropy referenced to 1 atm
        FF_HSfromTVPeosS(T,V,P,&data[i],&data[i].refH,&data[i].refS);
        if(refState==2){//IIR
            data[i].refH=data[i].refH-2.0e2*data[i].baseProp.MW;
            data[i].refS=data[i].refS-1*data[i].baseProp.MW;
        }
    }

    return &data[i];
}



//Basic functions
//---------------

void FF_eosCriticalConstantsM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, double *Tc, double *Pc){
    FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
    if (data->model==FF_CubicType){
        *Tc=data->cubicData.Tc;
        *Pc=data->cubicData.Pc;
    }
    else if (data->model==FF_SAFTtype){
        *Tc=data->saftData.Tc;
        *Pc=data->saftData.Pc;
    }
    else{
        *Tc=data->swData.Tc;
        *Pc=data->swData.Pc;
    }
}

//Pressure from d,T by EOS
void FF_pressureEOS_dT(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, double d, double T, double *P){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else printf("thermoModel out of range\n");
  double V=data->baseProp.MW/(d*1e3);
  FF_PfromTVeosS(T,V,data,P);
}

//Densities from T,P by EOS. Bellow the critical temperature
void FF_densities_pTM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, int aid, char *var, double P, double T, double *ld, double *gd){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else printf("thermoModel out of range\n");
  double MW=data->baseProp.MW*1e-3;//Mol weight in kg
  double Tc,Pc,Vc;
  //stablish the critical point
  if (data->model==FF_SWtype){
      Tc=data->swData.Tc;
      Pc=data->swData.Pc;
      Vc=1/data->swData.rhoRef;
  }
  else if(data->model==FF_SAFTtype){
      Tc=data->saftData.Tc;
      Pc=data->saftData.Pc;
      Vc=R*data->saftData.Zc*data->saftData.Tc/data->saftData.Pc;
  }
  else {
      Tc=data->cubicData.Tc;
      Pc=data->cubicData.Pc;
      Vc=R*data->cubicData.Zc*data->cubicData.Tc/data->cubicData.Pc;
  }
  if((var[0]=='e')||(var[0]=='E')){//Bubble state
      if (T>=Tc){
         *ld=MW/Vc;
      }
      else{
          if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)&&(aid==3)) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,ld);
          else if((data->model==FF_SAFTtype)&&(data->lDensCorr.form>0)&&(aid==2)) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,ld);
          else{
              double answerL[3],answerG[3];
              char option,state;
              option='l';
              FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
              *ld=MW/answerL[0];
              }
          *gd=0.0;
      }
      *gd=0.0;
  }
  else if((var[0]=='w')||(var[0]=='W')){//Dew state
      if (T>=Tc) *gd=MW/Vc;
      else{
          if ((data->model==FF_SWtype)&&(data->gDensCorr.form>0)&&(aid==3)) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,gd);
          else if ((data->model==FF_SAFTtype)&&(data->gDensCorr.form>0)&&(aid==2)) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,gd);
          else{
              double answerL[3],answerG[3];
              char option,state;
              option='g';
              FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
              *gd=MW/answerG[0];
          }
          *ld=0;
      }
      *ld=0;
  }
  else{//'l', 'g', or 'b' General case for nonsaturated conditions
      double answerL[3],answerG[3];//for cubic eos we can arrive to the critical point
      char state;
      FF_VfromTPeosS(T,P,data,var[0],answerL,answerG,&state);
      //printf("Modelica Vl:%f\n",answerL[0]);
      if (!(var[0]=='g')) *ld=data->baseProp.MW/answerL[0]/1e3;
      else *ld=0;
      if (!(var[0]=='l')) *gd=data->baseProp.MW/answerG[0]/1e3;
      else *gd=0;
  }

}

//Saturation pressure. For Modelica usage
//-----------------------------------------------------------------
void FF_saturationPressureM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, int aid, double T, double *Vp){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else printf("thermoModel out of range\n");
  if((data->model==FF_CubicType)&&(data->cubicData.Tc>0)&&(T>=data->cubicData.Tc)) *Vp=data->cubicData.Pc;
  else if((data->model==FF_SAFTtype)&&(data->saftData.Tc>0)&&(T>=data->saftData.Tc)) *Vp=data->saftData.Pc;
  else if((data->model==FF_SWtype)&&(data->swData.Tc>0)&&(T>=data->swData.Tc)) *Vp=data->swData.Pc;
  else if((data->model==FF_SWtype)&&(data->vpCorr.form>0) && !(data->swData.eos==FF_IAPWS95)&&(aid==3)) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->swData.MW,T,Vp);
  else if((data->model==FF_SAFTtype)&&(data->saftData.Tc>0)&&(T>=data->saftData.Tc)) *Vp=data->saftData.Pc;
  else if((data->model==FF_SAFTtype)&&(data->vpCorr.form>0) && (aid==2)) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,T,Vp);
  else FF_VpEosS(T,data,Vp);  //for water IAPWS uses a very accurate correlation
}


//Saturation temperature.  For Modelica usage
//--------------------------------------------------------------------
void FF_saturationTemperatureM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, int aid, double P, double *Tb){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else printf("thermoModel out of range\n");
  if((data->model==FF_CubicType)&&(data->cubicData.Pc>0)&&(P>=data->cubicData.Pc)) *Tb=data->cubicData.Tc;
  else if((data->model==FF_SAFTtype)&&(data->saftData.Pc>0)&&(P>=data->saftData.Pc)) *Tb=data->saftData.Tc;
  else if((data->model==FF_SWtype)&&(data->swData.Pc>0)&&(P>=data->swData.Pc)) *Tb=data->swData.Tc;
  else if((data->model==FF_SWtype)&&(data->btCorr.form>0) && !(data->swData.eos==FF_IAPWS95)&&(aid==3)) FF_PhysPropCorrM(data->btCorr.form,data->btCorr.coef,data->swData.MW,P,Tb);
  else if((data->model==FF_SAFTtype)&&(data->btCorr.form>0) && (aid==2)) FF_PhysPropCorrM(data->btCorr.form,data->btCorr.coef,data->swData.MW,P,Tb);
  else FF_TbEosS(P,data,Tb);
}

//all thermo properties from T/P or T/D or P/H or P/S, for Modelica usage
void CALLCONV FF_solveEosM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, int aid, char *variable, double x, double y, double *T, double *p, double *gd, double *gh, double *gs, double *gCv, double *gCp, double *gDvp,
                             double *gDvT, double *ld, double *lh, double *ls, double *lCv, double *lCp, double *lDvp, double *lDvT, double *gf){
    FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
    if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
    else if (thermoModel==2) data->model=FF_SAFTtype;
    else if (thermoModel==3) data->model=FF_SWtype;
    else printf("thermoModel out of range\n");
    FF_solveEos(variable,data,aid,x, y, T, p, gd, gh, gs, gCv, gCp, gDvp, gDvT, ld, lh, ls, lCv, lCp, lDvp, lDvT, gf);
}






//Densities, enthalpy and entropy from T,P by EOS.
void FF_dhsFrom_pTM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, int aid, char *var, double P, double T, double *ld, double *lh, double *ls, double *gd, double *gh, double *gs){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
  else if (thermoModel==2) data->model=FF_SAFTtype;
  else if (thermoModel==3) data->model=FF_SWtype;
  else printf("thermoModel out of range\n");
  double MW=data->baseProp.MW*1e-3;//Mol weight in kg
  double Tc,Pc,Vc;
  //stablish the critical point
  if (data->model==FF_SWtype){
      Tc=data->swData.Tc;
      Pc=data->swData.Pc;
      Vc=1/data->swData.rhoRef;
  }
  else if(data->model==FF_SAFTtype){
      Tc=data->saftData.Tc;
      Pc=data->saftData.Pc;
      Vc=R*data->saftData.Zc*data->saftData.Tc/data->saftData.Pc;
  }
  else {
      Tc=data->cubicData.Tc;
      Pc=data->cubicData.Pc;
      Vc=R*data->cubicData.Zc*data->cubicData.Tc/data->cubicData.Pc;
  }
  if((var[0]=='e')||(var[0]=='E')){//Bubble state
      if (T>=Tc){
          T=Tc;
         *ld=MW/Vc;
      }
      else{
          if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)&&(aid==3)) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,ld);
          else if((data->model==FF_SAFTtype)&&(data->lDensCorr.form>0)&&(aid==2)) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,ld);
          else{
              double answerL[3],answerG[3];
              char option,state;
              option='l';
              FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
              *ld=MW/answerL[0];
              }
      }
      *gd=0.0;
  }
  else if((var[0]=='w')||(var[0]=='W')){//Dew state
      if (T>=Tc){
          T=Tc;
          *gd=MW/Vc;
      }
      else{
          if ((data->model==FF_SWtype)&&(data->gDensCorr.form>0)&&(aid==3)) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,gd);
          else if ((data->model==FF_SAFTtype)&&(data->gDensCorr.form>0)&&(aid==2)) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,gd);
          else{
              double answerL[3],answerG[3];
              char option,state;
              option='g';
              FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
              *gd=MW/answerG[0];
          }
      }
      *ld=0;
  }
  else{//'l', 'g', or 'b' General case for nonsaturated conditions
      double answerL[3],answerG[3];//for cubic eos we can arrive to the critical point
      char state;
      FF_VfromTPeosS(T,P,data,var[0],answerL,answerG,&state);
      //printf("Modelica Vl:%f\n",answerL[0]);
      if (!(var[0]=='g')) *ld=data->baseProp.MW/answerL[0]/1e3;
      else *ld=0;
      if (!(var[0]=='l')) *gd=data->baseProp.MW/answerG[0]/1e3;
      else *gd=0;
  }
  if (*ld>0){
      FF_HSfromTVPeosS(T, MW / *ld, P, data, lh, ls);
      *lh=(*lh-data->refH)/MW;
      *ls=(*ls-data->refS)/MW;
  }
  else lh=ls=0;
  if (*gd>0){
      FF_HSfromTVPeosS(T, MW / *gd, P, data, gh, gs);
      *gh=(*gh-data->refH)/MW;
      *gs=(*gs-data->refS)/MW;
  }
  else gh=gs=0;
}

//Viscosity of a pure substance
void FF_ViscosityM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, const char *refName, double T,double dens, double P, double gf, double *visc){
    FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);//obtain the substance
    FF_SubstanceData *ref=FF_createSubstanceData(refName,resDir,thermoModel,refState,refT,refP);//obtain the reference substance
    if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
    else if (thermoModel==2) data->model=FF_SAFTtype;
    else if (thermoModel==3) data->model=FF_SWtype;
    else printf("thermoModel out of range\n");
    FF_Viscosity(data, ref, T, dens, P, gf, visc);
}

//Dynamic viscosity
void FF_dynamicViscosityM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, double T, double P, double gf, double *eta){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (gf==0.0){//only liquid phase
    if (data->lViscCorr.id>0) FF_PhysPropCorrM(data->lViscCorr.form,data->lViscCorr.coef,data->baseProp.MW,T,eta);
    else printf("unable to compute liquid viscosity\n");
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
    if(*eta==0) printf("unable to compute gas viscosity\n");
  }
  else printf("unable to compute two phases viscosity\n");
}



//Thermal conductivity
void FF_thermalConductivityM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, double T, double P, double gf, double *lambda){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if (gf==0.0){
    if (data->lThCCorr.id>0) FF_PhysPropCorrM(data->lThCCorr.form,data->lThCCorr.coef,data->baseProp.MW,T,lambda);
    else if(data->baseProp.Tc>0)FF_LiquidThCondLatini(T,&data->baseProp,lambda);
    else printf("unable to compute liquid thermal conductivity\n");
  }
  else if (gf==1.0){
    if ((data->gThCCorr.id>0)&&(T>=data->gThCCorr.limI)&&(T<=data->gThCCorr.limS)&&!(data->gThCCorr.form==122)) FF_PhysPropCorrM(data->gThCCorr.form,data->gThCCorr.coef,data->baseProp.MW,T,lambda);
    else if(data->cp0Corr.form>0){
        double Cp0;
        FF_PhysPropCorrM(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,T,&Cp0);
        FF_GasLpThCondTCpChung(T,Cp0,&data->baseProp,lambda);
    }
    else printf("unable to compute gas thermal conductivity\n");
  }
  else printf("unable to compute two phases thermal conductivity\n");
}

//Surface tension
void FF_surfaceTensionM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, double T, double *sigma){
  FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);
  if ((data->lSurfTCorr.id>0)&&(T>=data->lSurfTCorr.limI)&&(T<=data->lSurfTCorr.limS)) FF_PhysPropCorrM(data->lSurfTCorr.form,data->lSurfTCorr.coef,data->baseProp.MW,T,sigma);
  else if ((data->baseProp.Pa>0)&&(data->lDensCorr.form>0)&&((data->vpCorr.form>0)||((data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)))){
  FF_SurfTensMcLeod(T,data,sigma);
   }
  else if((data->baseProp.Tb>0)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_SurfTensSastri(T,&data->baseProp,sigma);
  else printf("unable to compute liquid surface tension\n");
}

//HERE BEGINS THE MODELICA INTERFACE FOR MIXTURES
//===============================================
/*
//Creates a static mixture array that is charged each time it is called
void *FF_createMixData(const char *name, int numSubs, const char subsNames[15][30], const char *resDir, int thermoModel, int eosType, int mixRule, int activityModel, int refCalc, double refT, double refP) {
    int i,j;
    static int n;//defined mixtures counter
    static char mixture[3][30];//stored mixture names
    static FF_MixData mixData[3];

    for (i=0;i<n;i++){
        if (strcmp(mixture[i],name)==0){
            printf("Found mixture:%i \n",i);
            break;
        }
    }
    if (i==n){
        n++;//we increment the number of defined mixtures
        strcpy(mixture[i],name);
        printf("New mixture:%s %i\n",mixture[i],i);
        FF_SubstanceData subsData[numSubs];
        for (j=0;j<numSubs;j++){
            char path[FILENAME_MAX]="";
            strcat(path,resDir);
            strcat(path,"/Fluids/");
            strcat(path,subsNames[j]);
            strcat(path,".sd");
            printf("%s\n",path);
            FILE * file= fopen(path, "rb");
            if (file != NULL) {
              fread(&subsData[j], sizeof(FF_SubstanceData), 1, file);
              fclose(file);
            }
            else printf("unable to charge the substance data\n");
            subsData[j].refT=0.0;//ideal gas calculation from 0K
            subsData[j].refP=101325;//ideal gas entropy referenced to 1 atm
        }
        FF_MixFillDataWithSubsData2(numSubs,subsData,&mixData[i]);

    }

    return &mixData[i];
}*/



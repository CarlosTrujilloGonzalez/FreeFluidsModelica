/*
 * FFmodelicaMedium.c
 *
 *  Created on: 5/9/2016
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2024  Carlos Trujillo Gonzalez

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
//#include "ModelicaUtilities.h"
#include "FFmodelicaMedium.h"
#include "FFbasic.h"


#include "FFeosPure.c"
#include "FFphysprop.c"
#include "FFeosMix.c"
#include "FFactivity.c"
#include "FFequilibrium.c"



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
      if (data->saftData.Vc>0) Vc=data->saftData.Vc;
      else Vc=R*data->saftData.Zc*data->saftData.Tc/data->saftData.Pc;
  }
  else {
      Tc=data->cubicData.Tc;
      Pc=data->cubicData.Pc;
      if (data->cubicData.Vc>0) Vc=data->cubicData.Vc;
      else Vc=R*data->cubicData.Zc*data->cubicData.Tc/data->cubicData.Pc;
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

//Dynamic viscosity by correlations or predictive, pressure correction for liquids
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
    if ((data->gViscCorr.id>0)&&(T>=data->gViscCorr.limI)&&(T<=data->gViscCorr.limS)&&!(data->gViscCorr.form==112)) FF_PhysPropCorrM(data->gViscCorr.form,data->gViscCorr.coef,data->baseProp.MW,T,&lpGasVisc);
    if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_GasViscTPcpLucas(T,P,&data->baseProp,&lpGasVisc,eta);
    else *eta=lpGasVisc;
    if(*eta==0) printf("unable to compute gas viscosity\n");
  }
  else printf("unable to compute two phases viscosity\n");
}

//thermal conductivity of a pure substance
void FF_ThCondM(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP, const char *refName, double T,double dens, double P, double gf, double *thCond){
    FF_SubstanceData *data=FF_createSubstanceData(name,resDir,thermoModel,refState,refT,refP);//obtain the substance
    FF_SubstanceData *ref=FF_createSubstanceData(refName,resDir,thermoModel,refState,refT,refP);//obtain the reference substance
    if (thermoModel==1) data->model=FF_CubicType;//because it could be that the substance has been previously created with another thermoModel
    else if (thermoModel==2) data->model=FF_SAFTtype;
    else if (thermoModel==3) data->model=FF_SWtype;
    else printf("thermoModel out of range\n");
    FF_ThCond(data, ref, T, dens, P, gf, thCond);
}



//Thermal conductivity by correlations or predictive, no density correction
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

//Creates a static mixture array that is charged each time it is called
//void *FF_createMixData(const char *name, int numSubs, const char subsNames[15][30], const char *resDir) {
void *FF_createMixData(const char *name, int numSubs, char *subsNamesOr, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel) {
    int i,j,k,count;
    static int n;//defined mixtures counter
    static char mixture[3][30];//stored mixture names
    static FF_MixData mixData[3];
    char subsNames[FILENAME_MAX]="";
    strcat(subsNames,subsNamesOr);
    //printf("numSubs:%i subsNames:%s resDir:%s\n",numSubs,subsNames,resDir);
    for (i=0;i<n;i++){
        if (strcmp(mixture[i],name)==0){
            //printf("Found mixture:%i \n",i);
            break;
        }
    }
    if (i==n){
        const char delimiters[] = ",";
        char *subsName;
        n++;//we increment the number of defined mixtures
        strcpy(mixture[i],name);
        //printf("New mixture:%s %i\n",mixture[i],i);
        FF_SubstanceData subsData[numSubs];
        subsName=strtok(subsNames,delimiters);

        for (j=0;j<numSubs;j++){
            char path[FILENAME_MAX]="";
            strcat(path,resDir);
            strcat(path,"/Fluids/");
            strcat(path,subsName);
            strcat(path,".sd");
            //printf("%s\n",path);
            FILE * file= fopen(path, "rb");
            if (file != NULL) {
              fread(&subsData[j], sizeof(FF_SubstanceData), 1, file);
              fclose(file);
            }
            else printf("unable to charge the substance data\n");
            subsData[j].refT=0.0;//ideal gas calculation from 0K
            subsData[j].refP=101325;//ideal gas entropy referenced to 1 atm
            subsName=strtok(NULL,delimiters);
        }
        char path[FILENAME_MAX]="";
        strcat(path,resDir);
        strcat(path,"/Extra/");
        FF_MixFillDataWithSubsData2(numSubs,subsData,path,&mixData[i]);
        mixData[i].numSubs=numSubs;
        mixData[i].thModelActEos=1;
        mixData[i].refVpEos=1;
        mixData[i].refT=0.0;
        mixData[i].refP=101325.0;
        mixData[i].eosType=FF_CubicPRtype;//EOS type by default
        if (strcmp(eosType,"SRK")==0) mixData[i].eosType=FF_CubicSRKtype;
        else if(strcmp(eosType,"PCSAFT")==0) mixData[i].eosType=FF_SAFTtype;

        if((mixData[i].eosType==FF_CubicPRtype)||(mixData[i].eosType==FF_CubicSRKtype)){
            mixData[i].actModel=FF_UNIFACDort;
            if(strcmp(activityModel,"None")==0) mixData[i].actModel=FF_NoModel;
            else if(strcmp(activityModel,"UNIFACstd")==0) mixData[i].actModel=FF_UNIFACStd;
            else if(strcmp(activityModel,"UNIFACpsrk")==0) mixData[i].actModel=FF_UNIFACPSRK;
            else if(strcmp(activityModel,"UNIQUAC")==0) mixData[i].actModel=FF_UNIQUAC;
            else if(strcmp(activityModel,"NRTL")==0) mixData[i].actModel=FF_NRTL;

            mixData[i].mixRule=FF_LCVM;//default mixing rule
            if(strcmp(cubicMixRule,"VdWnoInt")==0) mixData[i].mixRule=FF_VdWnoInt;
            else if(strcmp(cubicMixRule,"VdW")==0) mixData[i].mixRule=FF_VdW;
            else if (strcmp(cubicMixRule,"HV")==0) mixData[i].mixRule=FF_HV;
            else if (strcmp(cubicMixRule,"MHV1")==0) mixData[i].mixRule=FF_MHV1;
            else if (strcmp(cubicMixRule,"MHV2")==0) mixData[i].mixRule=FF_MHV2;
            else if (strcmp(cubicMixRule,"UMR")==0) mixData[i].mixRule=FF_UMR;
            else if (strcmp(cubicMixRule,"PSRK")==0) mixData[i].mixRule=FF_PSRK;
        }
        else {
            if(strcmp(cubicMixRule,"IndAssoc")==0) mixData[i].mixRule=FF_IndAssoc;
            //printf("mixrule:%i \n",mixData[i].mixRule);
        }

        //load of interaction parameters
        char path2[FILENAME_MAX]="";
        char cas1[22],cas2[22],form[20];
        float a,b,c,d,e,f,ai,bi,ci,di,ei,fi;
        strcat(path2,resDir);
        if (mixData[i].eosType==FF_SAFTtype){
            strcat(path2,"/Interactions/PCSAFT.txt");
            FILE * file= fopen(path2, "rb");
            if (file != NULL) {
              count=0;
              while (fscanf(file,"%s %s %f %f %f",cas1,cas2,&a,&b,&c)==5){
                  for (j=0;j<numSubs;j++){
                      for (k=0;k<numSubs;k++){
                          if((strcmp(cas1,mixData[i].CAS[j])==0)&&(strcmp(cas2,mixData[i].CAS[k])==0)){
                              count++;
                              mixData[i].intParam[j][k][0]=mixData[i].intParam[k][j][0]=a;
                              mixData[i].intParam[j][k][1]=mixData[i].intParam[k][j][1]=b;
                              mixData[i].intParam[j][k][2]=mixData[i].intParam[k][j][2]=c;
                          }
                      }

                  }
                  //if (count==(numSubs*(numSubs-1)/2)) break;//problematic if there are duplicate values
              }
              //printf("cas1:%s cas2:%s a:%f d:%f intParam[0][1][0]:%f intParam[0][1][3]:%f\n",cas1,cas2,a,d,mixData[i].intParam[0][1][0],mixData[i].intParam[0][1][3]);
              fclose(file);
            }
            else printf("unable to charge the interaction parameters\n");
        }
        else if (mixData[i].mixRule==FF_VdWnoInt);
        else if (mixData[i].mixRule==FF_VdW){
            if (mixData[i].eosType==FF_CubicPRtype) strcat(path2,"/Interactions/VdWPR.txt");
            if (mixData[i].eosType==FF_CubicSRKtype) strcat(path2,"/Interactions/VdWSRK.txt");
            FILE * file= fopen(path2, "rb");
            if (file != NULL) {
              count=0;
              while (fscanf(file,"%s %s %f %f %f %f",cas1,cas2,&a,&b,&c,&d)==6){
                  for (j=0;j<numSubs;j++){
                      for (k=0;k<numSubs;k++){
                          if((strcmp(cas1,mixData[i].CAS[j])==0)&&(strcmp(cas2,mixData[i].CAS[k])==0)){
                              count++;
                              mixData[i].intParam[j][k][0]=mixData[i].intParam[k][j][0]=a;
                              mixData[i].intParam[j][k][1]=mixData[i].intParam[k][j][1]=b;
                              mixData[i].intParam[j][k][2]=mixData[i].intParam[k][j][2]=c;
                              mixData[i].intParam[j][k][3]=mixData[i].intParam[k][j][3]=d;
                          }
                      }

                  }
              }
              //printf("cas1:%s cas2:%s a:%f d:%f intParam[0][1][0]:%f intParam[0][1][3]:%f\n",cas1,cas2,a,d,mixData[i].intParam[0][1][0],mixData[i].intParam[0][1][3]);
              fclose(file);
            }
            else printf("unable to charge the interaction parameters\n");
        }
        else if (mixData[i].actModel==FF_UNIQUAC){
            strcat(path2,"/Interactions/UNIQUAC.txt");
            FILE * file= fopen(path2, "rb");
            if (file != NULL) {
              count=0;
              while (fscanf(file,"%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f",cas1,cas2,form,&a,&b,&c,&d,&e,&f,&ai,&bi,&ci,&di,&ei,&fi)==15){
                  for (j=0;j<numSubs;j++){
                      for (k=0;k<numSubs;k++){
                          if((strcmp(cas1,mixData[i].CAS[j])==0)&&(strcmp(cas2,mixData[i].CAS[k])==0)){
                              count++;
                              if (mixData[i].intParam[j][k][0]==0){
                                  if (strcmp(form,"Pol1K")==0) mixData[i].intForm=FF_Pol1K;
                                  else if (strcmp(form,"Pol1J")==0) mixData[i].intForm=FF_Pol1J;
                                  else if (strcmp(form,"Pol1C")==0) mixData[i].intForm=FF_Pol1C;
                                  mixData[i].intParam[j][k][0]=a;
                                  mixData[i].intParam[j][k][1]=b;
                                  mixData[i].intParam[j][k][2]=c;
                                  mixData[i].intParam[j][k][3]=d;
                                  mixData[i].intParam[j][k][4]=e;
                                  mixData[i].intParam[j][k][5]=f;

                                  mixData[i].intParam[k][j][0]=ai;
                                  mixData[i].intParam[k][j][1]=bi;
                                  mixData[i].intParam[k][j][2]=ci;
                                  mixData[i].intParam[k][j][3]=di;
                                  mixData[i].intParam[k][j][4]=ei;
                                  mixData[i].intParam[k][j][5]=fi;
                              }
                          }
                      }

                  }
              }

              //printf("cas1:%s cas2:%s a:%f d:%f intParam[0][1][0]:%f intParam[0][1][3]:%f\n",cas1,cas2,a,d,mixData[i].intParam[0][1][0],mixData[i].intParam[0][1][3]);
              fclose(file);
            }
            else printf("unable to charge the interaction parameters\n");
        }
        else if (mixData[i].actModel==FF_NRTL){
            strcat(path2,"/Interactions/NRTL.txt");
            FILE * file= fopen(path2, "rb");
            if (file != NULL) {
              count=0;
              while (fscanf(file,"%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f",cas1,cas2,form,&a,&b,&c,&d,&e,&f,&ai,&bi,&ci,&di,&ei,&fi)==15){
                  for (j=0;j<numSubs;j++){
                      for (k=0;k<numSubs;k++){
                          if((strcmp(cas1,mixData[i].CAS[j])==0)&&(strcmp(cas2,mixData[i].CAS[k])==0)){
                              count++;
                              if (mixData[i].intParam[j][k][0]==0){
                                  if (strcmp(form,"Pol1K")==0) mixData[i].intForm=FF_Pol1K;
                                  else if (strcmp(form,"Pol1J")==0) mixData[i].intForm=FF_Pol1J;
                                  else if (strcmp(form,"Pol1C")==0) mixData[i].intForm=FF_Pol1C;
                                  mixData[i].intParam[j][k][0]=a;
                                  mixData[i].intParam[j][k][1]=b;
                                  mixData[i].intParam[j][k][2]=c;
                                  mixData[i].intParam[j][k][3]=d;
                                  mixData[i].intParam[j][k][4]=e;
                                  mixData[i].intParam[j][k][5]=f;

                                  mixData[i].intParam[k][j][0]=ai;
                                  mixData[i].intParam[k][j][1]=bi;
                                  mixData[i].intParam[k][j][2]=ci;
                                  mixData[i].intParam[k][j][3]=di;
                                  mixData[i].intParam[k][j][4]=ei;
                                  mixData[i].intParam[k][j][5]=fi;
                              }
                          }
                      }

                  }
              }

              //printf("cas1:%s cas2:%s a:%f d:%f intParam[0][1][0]:%f intParam[0][1][3]:%f\n",cas1,cas2,a,d,mixData[i].intParam[0][1][0],mixData[i].intParam[0][1][3]);
              fclose(file);
            }
            else printf("unable to charge the interaction parameters\n");
        }



    }
    return &mixData[i];
}

//molar masses of components
void FF_molarMassesM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double molarMass[15]){
  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  for(int i=0;i<numSubs;i++){
      molarMass[i]=0.001*mix->baseProp[i].MW;
      //printf("%f %f\n",mix->baseProp[i].MW),molarMass[i];
  }
}

//Pressure from T,D,X by EOS.
void FF_pressure_dTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double D, double T, double zMass[15], double *p, double *MW){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  *MW=0;
  double z[numSubs];
  double V;
  for(int i=0;i<numSubs;i++){
      nMols=nMols+zMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
      *MW=*MW+z[i]*mix->baseProp[i].MW;
      //printf("z[%i]:%f\n",i,z[i]);
  }
  *MW=*MW*0.001;
  //printf("MW in kgr/mol:%f\n",MW);
  V=*MW/D;

  FF_MixPfromTVeos(mix,&T,&V,z,p);
}

//Density from T,P,X by EOS.
void FF_density_pTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, char *var, double P, double T, double zMass[15], double *ld, double *gd, double *MW){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  double nMols=0;//number of moles in 1 gr
  *MW=0;
  double z[numSubs];
  for(int i=0;i<numSubs;i++){
      nMols=nMols+zMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
      //*MW=*MW+z[i]*mix->baseProp[i].MW;
      //printf("z[%i]:%f\n",i,z[i]);
  }
  //*MW=*MW*0.001;
  *MW=0.001/nMols;
  //printf("MW in kgr/mol:%f\n",MW);

  char state;
  double answerL[3],answerG[3];
  FF_MixVfromTPeos(mix,&T,&P,z,var,answerL,answerG,&state);

  //mix->eosType=FF_CubicPRtype;//type of EOS to use, at least for the gas phase. Also for the liquid if the EOS is cubic and thModelActEos==0
  //FF_MixVfromTPeos(mix,&T,&P,z,&option,answerL,answerG,&state);
  //*d=*MW/answerL[0];
  *ld=*MW/answerL[0];
  *gd=*MW/answerG[0];
}

//Thermodynamic properties from T,D,X by EOS.
void FF_thermo_dTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel,
                    double D, double T, double zMass[15], double *p, double *h, double *s, double *Cv, double *Cp, double *Dvp, double *DvT){
  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  FF_PhaseThermoProp th;
  double nMols=0;
  th.MW=0;
  for(int i=0;i<numSubs;i++){
      nMols=nMols+zMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      th.c[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
      th.MW=th.MW+th.c[i]*mix->baseProp[i].MW;
      //printf("z[%i]:%f\n",i,z[i]);
  }
  double MW=th.MW*0.001;
  th.T=T;
  th.V=MW/D;
  //printf("th.T:%f th.V:%f MW in gr/mol:%f\n",th.T,th.V,th.MW);

  FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
  *p=th.P;
  *h=th.H/MW;
  *s=th.S/MW;
  *Cv=th.Cv/MW;
  *Cp=th.Cp/MW;
  *Dvp=1/(th.dP_dV*MW);
  *DvT=-th.dP_dT/(th.dP_dV*MW);
  //printf("D:%f T:%f Cp:%f J/mol MW:%f kg/mol Cp:%f J/kg\n",D,T,th.Cp,MW,*Cp);
}

void FF_activity_TXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel,
                     double T, double xMass[15], double gamma[15], double*gE){
   FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
   double nMols=0;
   double x[numSubs];
   FF_SubsActivityData actData[numSubs];
   for(int i=0;i<numSubs;i++){
       nMols=nMols+xMass[i]/mix->baseProp[i].MW;
   }
   for(int i=0;i<numSubs;i++){
       x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
   }
   FF_Activity(mix,&T,x,actData);
   *gE=0;
   for(int i=0;i<numSubs;i++){
       gamma[i]=actData[i].gamma;
       *gE=*gE+x[i]*log(gamma[i]); //(J/mol)
   }
   *gE=*gE*1000*nMols;//(J/kg)

}

void FF_fugacity_pTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel,
                     double p, double T, double xMass[15], double phiL[15], double phiG[15], double kWilson[15]){
   FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
   double nMols=0;
   double x[numSubs];
   FF_SubsActivityData actData[numSubs];
   for(int i=0;i<numSubs;i++){
       nMols=nMols+xMass[i]/mix->baseProp[i].MW;
   }
   for(int i=0;i<numSubs;i++){
       x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
   }
   char option='l';
   FF_MixPhiEOS(mix,&T,&p,x,&option,phiL);
   option='g';
   FF_MixPhiEOS(mix,&T,&p,x,&option,phiG);
   for(int i=0;i<numSubs;i++){
       kWilson[i]=exp(log(mix->baseProp[i].Pc/ p)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
   }
}


//bubble pressure and composition from T,X by EOS.
void FF_bubble_TXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double T, double xMass[15], double *p, double yMass[15]){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  double nMols=0;
  double MW=0;
  double x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
  for(int i=0;i<numSubs;i++){
      nMols=nMols+xMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
      //printf("FF_bubble_TXM x[%i]:%f\n",i,x[i]);
  }

  double pGuess=0;
  FF_BubbleP(mix,&T,x,&pGuess,p,y,phiL,phiG);
  for(int i=0;i<numSubs;i++){
      MW=MW+y[i]*mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      yMass[i]=y[i]*mix->baseProp[i].MW/MW;
  }
}

//bubble temperature and composition from p,X by EOS.
void FF_bubble_pXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double xMass[15], double *T, double yMass[15]){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  double nMols=0;
  double MW=0;
  double x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];;
  for(int i=0;i<numSubs;i++){
      nMols=nMols+xMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
  }

  double TGuess=0;
  FF_BubbleT(mix,&p,x,&TGuess,T,y,phiL,phiG);
  for(int i=0;i<numSubs;i++){
      MW=MW+y[i]*mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      yMass[i]=y[i]*mix->baseProp[i].MW/MW;
  }
}

//dew pressure and composition from T,y by EOS.
void FF_dew_TXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double T, double yMass[15], double *p, double xMass[15]){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  double nMols=0;
  double MW=0;
  double x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];;
  for(int i=0;i<numSubs;i++){
      nMols=nMols+yMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      y[i]=yMass[i]/(mix->baseProp[i].MW*nMols);
  }

  double pGuess=0;
  FF_DewP(mix,&T,y,&pGuess,p,x,phiL,phiG);
  for(int i=0;i<numSubs;i++){
      MW=MW+x[i]*mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      xMass[i]=x[i]*mix->baseProp[i].MW/MW;
  }
}

//dew temperature and composition from p,y by EOS.
void FF_dew_pXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double yMass[15], double *T, double xMass[15]){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

  double nMols=0;
  double MW=0;
  double x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
  for(int i=0;i<numSubs;i++){
      nMols=nMols+yMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      y[i]=yMass[i]/(mix->baseProp[i].MW*nMols);
  }

  double TGuess=0;
  FF_DewT(mix,&p,y,&TGuess,T,x,phiL,phiG);
  for(int i=0;i<numSubs;i++){
      MW=MW+x[i]*mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      xMass[i]=x[i]*mix->baseProp[i].MW/MW;
  }
}

//VL flash calculation, given T, P, feed composition, eos and mixing rule
void FF_TwoPhasesFlashPTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double zMass[],
                           double xMass[],double yMass[],double *gfMass){
    FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

    double nMols=0;
    double lMW=0,gMW=0;
    double gf=0.33;//gas molar fraction
    double z[numSubs],x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
    int i;
    for(int i=0;i<numSubs;i++){
        nMols=nMols+zMass[i]/mix->baseProp[i].MW;
    }
    for(i=0;i<numSubs;i++){
        z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
    }

    FF_TwoPhasesPreFlashPT(mix,&T,&p,z,x,y,phiL,phiG,&gf);
    if (gf==0){
        *gfMass=0;
        for(i=0;i<numSubs;i++){
            xMass[i]=zMass[i];
            yMass[i]=0.0;
        }
    }
    else if (gf==1){
        *gfMass=1;
        for(i=0;i<numSubs;i++){
            yMass[i]=zMass[i];
            xMass[i]=0.0;
        }
    }
    else{
        for(int i=0;i<numSubs;i++){
            lMW=lMW+x[i]*mix->baseProp[i].MW;
            gMW=gMW+y[i]*mix->baseProp[i].MW;
        }
        for(int i=0;i<numSubs;i++){
            xMass[i]=x[i]*mix->baseProp[i].MW/lMW;
            yMass[i]=y[i]*mix->baseProp[i].MW/gMW;
            //xMass[i]=x[i];
            //yMass[i]=y[i];

        }
        *gfMass=gf*gMW/(gf*gMW+(1-gf)*lMW);
    }

}

//VL flash calculation, given T, P, feed composition, eos and mixing rule
void FF_TwoPhasesFlashPTM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double zMass[],
                           double xMass[],double yMass[],double *gfMass){
    FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);


    double nMols=0;
    double lMW=0,gMW=0;
    double gf=0;//gas molar fraction
    double z[numSubs],x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
    for(int i=0;i<numSubs;i++){
        nMols=nMols+zMass[i]/mix->baseProp[i].MW;
    }
    for(int i=0;i<numSubs;i++){
        z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
    }
    FF_TwoPhasesFlashPT(mix,&T,&p,z,x,y,phiL,phiG,&gf);

    for(int i=0;i<numSubs;i++){
        lMW=lMW+x[i]*mix->baseProp[i].MW;
        gMW=gMW+y[i]*mix->baseProp[i].MW;
    }
    for(int i=0;i<numSubs;i++){
        xMass[i]=x[i]*mix->baseProp[i].MW/lMW;
        yMass[i]=y[i]*mix->baseProp[i].MW/lMW;
    }
    *gfMass=gf*gMW/(gf*gMW+(1-gf)*lMW);
}

//VL flash calculation, given H or S, P, feed composition, eos and mixing rule
void FF_TwoPhasesFlashP_HS_XM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, char *election, double p, double e, double *T, double zMass[],
                           double xMass[],double yMass[],double *gfMass){
    FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

    double nMols=0;//number of moles in 1 gr
    double lMW=0,gMW=0;
    double gf=0.33;//gas molar fraction
    double z[numSubs],x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
    int i;
    double eMolar;
    for(int i=0;i<numSubs;i++){
        nMols=nMols+zMass[i]/mix->baseProp[i].MW;
    }
    for(i=0;i<numSubs;i++){
        z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
    }
    eMolar=e/(1000*nMols);
    //printf("z molar 0:%f eMolar:%f \n",z[0],eMolar);



    FF_TwoPhasesPreFlashP_HS(mix,*election,&eMolar,&p,z,T,x,y,phiL,phiG,&gf);
    if (gf==0){
        *gfMass=0;
        for(i=0;i<numSubs;i++){
            xMass[i]=zMass[i];
            yMass[i]=0.0;
        }
    }
    else if (gf==1){
        *gfMass=1;
        for(i=0;i<numSubs;i++){
            yMass[i]=zMass[i];
            xMass[i]=0.0;
        }
    }
    else{
        for(int i=0;i<numSubs;i++){
            lMW=lMW+x[i]*mix->baseProp[i].MW;
            gMW=gMW+y[i]*mix->baseProp[i].MW;
        }
        for(int i=0;i<numSubs;i++){
            xMass[i]=x[i]*mix->baseProp[i].MW/lMW;
            yMass[i]=y[i]*mix->baseProp[i].MW/gMW;
        }
        *gfMass=gf*gMW/(gf*gMW+(1-gf)*lMW);
    }

}

//VL flash calculation, given gas fraction, P, feed composition, eos and mixing rule
void FF_TwoPhasesFlashPThetaXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, char *opt, double p, double Theta, double *T, double zMass[],
                           double xMass[],double yMass[]){
    FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);

    double nMols=0;//number of moles in 1 gr
    double lMW=0,gMW=0;
    double gf=0.33;//gas molar fraction
    double z[numSubs],x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
    int i;
    for(int i=0;i<numSubs;i++){
        nMols=nMols+zMass[i]/mix->baseProp[i].MW;
    }
    for(i=0;i<numSubs;i++){
        z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
    }

    FF_TwoPhasesPreFlashPTheta(mix,opt, &Theta,&p,z,T,x,y,phiL,phiG);

    for(int i=0;i<numSubs;i++){
        lMW=lMW+x[i]*mix->baseProp[i].MW;
        gMW=gMW+y[i]*mix->baseProp[i].MW;
    }
    for(int i=0;i<numSubs;i++){
        xMass[i]=x[i]*mix->baseProp[i].MW/lMW;
        yMass[i]=y[i]*mix->baseProp[i].MW/gMW;
    }
    //*gfMass=gf*gMW/(gf*gMW+(1-gf)*lMW);

}


//VL flash calculation, given T, P, feed composition, eos and mixing rule
void FF_RachfordRiceSolverPTXM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double zMass[],
                           double kInit[], double xMass[],double yMass[],double *gfMass){
    FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);


    double nMols=0;
    double lMW=0,gMW=0;
    double gf=0;//gas molar fraction
    int numRep=16;
    double z[numSubs],x[numSubs],y[numSubs],phiL[numSubs],phiG[numSubs];
    for(int i=0;i<numSubs;i++){
        nMols=nMols+zMass[i]/mix->baseProp[i].MW;
    }
    for(int i=0;i<numSubs;i++){
        z[i]=zMass[i]/(mix->baseProp[i].MW*nMols);
    }
    FF_RachfordRiceSolver(mix,&T,&p,z,kInit,&numRep,x,y,phiL,phiG,&gf);

    for(int i=0;i<numSubs;i++){
        lMW=lMW+x[i]*mix->baseProp[i].MW;
        gMW=gMW+y[i]*mix->baseProp[i].MW;
    }
    for(int i=0;i<numSubs;i++){
        xMass[i]=x[i]*mix->baseProp[i].MW/lMW;
        yMass[i]=y[i]*mix->baseProp[i].MW/lMW;
    }
    *gfMass=gf*gMW/(gf*gMW+(1-gf)*lMW);
}

//Liquid mix viscosity from p,T,X by Teja-Rice method
void FF_mixLiquidViscosityM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double xMass[], double *eta){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  double x[numSubs]; //liquid molar fractions
  for(int i=0;i<numSubs;i++){
      nMols=nMols+xMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
  }

  FF_MixLiqViscTeja(mix,&T,&p,x,eta);
}

//Gas mix viscosity from p,T,X by Lucas method
void FF_mixGasViscosityM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double yMass[], double *eta){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  double y[numSubs]; //liquid molar fractions
  for(int i=0;i<numSubs;i++){
      nMols=nMols+yMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      y[i]=yMass[i]/(mix->baseProp[i].MW*nMols);
  }

  FF_MixGasViscTPcpLucas(mix,T,p,y,eta);
}

//Liquid mix thermal conductivity from p,T,X by Latini method
void FF_mixLiquidThCondM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double xMass[], double *lambda){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  double x[numSubs]; //liquid molar fractions
  for(int i=0;i<numSubs;i++){
      nMols=nMols+xMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
  }

  FF_MixLiqThCondLi(mix,&T,&p,x,lambda);
}

//Gas mix thermal conductivity from p,T,X by Mason method
void FF_mixGasThCondM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double yMass[], double *lambda){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  double y[numSubs]; //liquid molar fractions
  for(int i=0;i<numSubs;i++){
      nMols=nMols+yMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      y[i]=yMass[i]/(mix->baseProp[i].MW*nMols);
  }
  FF_MixLpGasThCondTpMason(mix,T,y,lambda);
}

//Liquid mix thermal conductivity from p,T,X by Latini method
void FF_mixSurfaceTensionM(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel, double p, double T, double xMass[], double *sigma){

  FF_MixData* mix = FF_createMixData(name,numSubs,subsNames,resDir,eosType,cubicMixRule,activityModel);
  double nMols=0;
  double x[numSubs]; //liquid molar fractions
  for(int i=0;i<numSubs;i++){
      nMols=nMols+xMass[i]/mix->baseProp[i].MW;
  }
  for(int i=0;i<numSubs;i++){
      x[i]=xMass[i]/(mix->baseProp[i].MW*nMols);
  }
 FF_MixLiqSurfTensWinterfeld(mix,T,p,x,sigma);
}


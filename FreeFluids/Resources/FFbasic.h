/*
 * FFbasic.h
 *
 *  Created on: 11/03/2018
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2019  Carlos Trujillo Gonzalez

 *This program is free software; you can redistribute it and/or
 *modify it under the terms of the GNU General Public License version 3
 *as published by the Free Software Foundation

 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.

 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

// contains basic definitions for FreeFluids

#ifndef FFBASIC_H
#define FFBASIC_H

#if (defined(_WIN32) || defined(__WIN32__))
  #define CALLCONV __cdecl//Definition of the calling convention, option is __stdcall
  /*** You should define EXPORTS at the gcc command line only when building the shared library ***/
  #ifdef EXPORTS
    #define EXP_IMP __declspec(dllexport)/*used to declare functions as export*/
  #else
    #ifdef IMPORTS
      #define EXP_IMP __declspec(dllimport)/*used to declare functions as import*/
    #else
      #define EXP_IMP
    #endif
  #endif
#else/*If we are not in Windows, we give EXP_IMP and CALLCONV null value*/
  #define EXP_IMP
  #define CALLCONV
#endif

#include <stdbool.h>

//Global variables
//----------------
extern const double Av; //molecules/mol
extern const double kb;//J/K
extern const double Pi;
extern const double R; //Pa·m3/(K·mol)
extern const double FF_PCSAFTap[7][3];
extern const double FF_PCSAFTbp[7][3];

//Enumerations
//------------
enum FF_EosType{FF_NoType,FF_IdealType,FF_CubicType,FF_CubicPRtype,FF_CubicSRKtype,FF_SAFTtype,FF_SWtype};
enum FF_EOS{FF_IdealGas,FF_PR76,FF_PR78,FF_PRSV1,FF_PRBM,FF_PRMELHEM,FF_PRSOF,FF_PRALMEIDA,FF_PRMC,FF_PRTWU91,FF_PRvTWU91,FF_PRTWU95,FF_PRPOL1,FF_PRFIT3,FF_PRFIT4,FF_PRFIT4B,
            FF_SRK,FF_SRKSOF,FF_SRKMC,FF_SRKTWU91,FF_SRKPOL2,
            FF_PCSAFT,FF_PCSAFT1A,FF_PCSAFT2A,FF_PCSAFT2B,FF_PCSAFT3B,FF_PCSAFT4C,FF_PPCSAFT_GV,FF_PPCSAFT_JC,
            FF_PPCSAFT1A_GV,FF_PPCSAFT2B_GV,FF_PPCSAFT2B_JC,FF_PPCSAFT3B_GV,FF_PPCSAFT4C_GV,FF_PCSAFTPOL1,
            FF_SAFTVRMie,FF_SAFTVRMie1A,FF_SAFTVRMie2B,FF_SAFTVRMie4C,FF_PSAFTVRMie_GV,FF_PSAFTVRMie_JC,FF_SAFTVRMie2B_GV,FF_SAFTVRMie4C_GV,
            FF_SW,FF_IAPWS95,IF97};
enum FF_MixingRule{FF_NoMixRul,FF_VdW,FF_PR,FF_MKP,FF_HV,FF_MHV1,FF_PSRK,FF_HVOS,FF_LCVM,FF_MHV2,FF_UMR,FF_OPTgE,FF_PSRKnew,FF_VTPR,FF_VdWnoInt,FF_BL};
enum FF_CorrEquation{FF_DIPPR100,FF_Polynomial,FF_Polynomial2,FF_DIPPR100Ld,FF_expDIPPR100,FF_DIPPR101,FF_DIPPR101Vp,FF_DIPPR101Lv,FF_logDIPPR101,
                  FF_DIPPR102,FF_DIPPR103,FF_DIPPR104,FF_DIPPR105,FF_DIPPR106,FF_DIPPR106Hv,FF_DIPPR106Ld,FF_DIPPR106SurfT,FF_DIPPR107,
                  FF_DIPPR107Cp,FF_DIPPR114,FF_DIPPR115,FF_DIPPR116,FF_DIPPR116Ld,FF_Wilhoit,FF_Cooper,FF_Jaechske,FF_ChemSep16,FF_Antoine1,
                  FF_Antoine2,FF_Wagner25,FF_Wagner36,FF_PPDS9,FF_PPDS10,FF_PCWIN,FF_Rackett,FF_ExtAndrade1,FF_ExtAndrade2,FF_ChericVisc,FF_WagnerGd,FF_Tait,FF_ExtWagner,FF_PPDS15};

enum FF_ActModel{FF_NoModel,FF_Wilson,FF_NRTL,FF_UNIQUAC,FF_UNIQUACFV,FF_UNIFACStd,FF_UNIFACPSRK,FF_UNIFACDort,FF_UNIFACNist,FF_EntropicFV,FF_UNIFACZM,
                 FF_Hildebrand,FF_Hansen,FF_Chi};
//Enumeration for binary interaction parameter calculation formulas
enum FF_IntParamForm{FF_NoForm,FF_Pol1,FF_Pol1K,FF_Pol1J,FF_Pol1C,FF_Pol2,FF_Pol2K,FF_Pol2J,FF_Pol2C,FF_Pol3,FF_Pol3K,FF_Pol3J,FF_Pol3C,FF_Ant1,FF_Ant2,FF_Ant3};
//FF_Pol1 corresponds to:a+b*T+c*T^2
//FF_Pol2: a+b*T+c/T^2
//FF_Pol3: a+b/T+c*T
//FF_Ant1: a+b/T+c*ln(T)+d*T
//FF_Ant2: a+b/T+c*ln(T)+d*T+e*T^2
//FF_Ant3: a+b/T+c*ln(T)+d*T+e/T^2.
//In case the result contains energy, K,J,C indicates if the energy is expresed in K, Joules o Calories
//Enumeration for units used in Flory-Huggins model
enum FF_Units{FF_cal_cm3_05,FF_MPa_05};
enum FF_SubstanceType{FF_NoFamily,FF_Alkane,FF_Alkene,FF_Alkyne,FF_Cycloalkane,FF_Aromatic,FF_Water,FF_Alcohol,FF_Polyol,FF_Phenol,FF_Ether,FF_Aldehyde,FF_Ketone,FF_Acid,FF_Ester,FF_Amine,FF_Polymer};
//Structures
//----------
//Basic properties for pure or pseudopure substances
typedef struct {int id,type,numMono,unifacPSRKSubg[10][2],unifacDortSubg[10][2];double MW,MWmono,Tc,Pc,Vc,Zc,w,Zra,r,q,qRes,VdWV,Hf0g,Gf0g,S0g,
                Pa,Vliq,FV,mu,Q,RadGyr,Tm,Hm,Tb,Hildebrand,HansenD,HansenP,HansenH,LnuA,LnuB;}FF_BaseProp;
//Data for cubic EOS. Critical properties used are reported in order to better fix the EOS
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,VdWV,c,k1,k2,k3,k4;}FF_CubicEOSdata;
//Cubic EOS parameters once given composition and T
typedef struct {double a,Theta,b,c,u,w,dTheta,d2Theta,Tc,Pc,Zc;}FF_CubicParam;
//Data for SAFT type EOS. Critical properties used are reported in order to better fix the EOS
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,sigma,m,epsilon,lr,la,chi,kAB,epsilonAB,mu,xp,Q;int nPos,nNeg,nAcid;}FF_SaftEOSdata;
//Data for Schmidt-Wagner type EOS
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,tRef,rhoRef,n[60],t[60],a[16],e[16],b[16],g[16],af[5],bf[5],
                Af[5],Bf[5],Cf[5],Df[5],betaf[5];int d[60],c[55],nPol,nExp,nSpec,nFinal;}FF_SWEOSdata;
//Data for physical properties correlation formulas
typedef struct {int id,form;double coef[14],limI,limS;}FF_Correlation;
//Data for single physical property
typedef struct {double x,y;}FF_SinglePointData;

//Data for a pure substance, including substructures for basic data, EOS and physical properties
typedef struct {char name[50],CAS[22],description[150];int id,model,UnifStdSubg[20][2],UnifPSRKSubg[20][2],UnifDortSubg[20][2],UnifNistSubg[20][2];
                double refT,refP,refH,refS;FF_BaseProp baseProp;FF_SinglePointData RI,cp0,vp,hVsat,lCp,lDens,lVisc,lThC,lSurfT,lIsothComp,gVisc,gThC,sDens,sCp;
                FF_CubicEOSdata cubicData;FF_SaftEOSdata saftData;FF_SWEOSdata swData;FF_Correlation cp0Corr,vpCorr,btCorr,hVsatCorr,lCpCorr,
                lTfromHCorr,lDensCorr,lViscCorr,lThCCorr,lSurfTCorr,lBulkModRCorr,gDensCorr,gTfromDcorr,gViscCorr,gThCCorr,sDensCorr,sCpCorr;}FF_SubstanceData;

//UNIFAC data for a mixture. Used to speed-up calculations
//Prepared for 20 substances and 30 subgroups. FV must be filled with the free volume of each substance if EntropicFV model is to be used
//subgroups:List of subgroups with its groups.
typedef struct {int model,numSubs,numSubg,subgroup[30][2],subsSubg[20][30];double subgData[30][2], subgInt[30][30][3],subsR[20],subsQ[20],
                FV[20];}FF_UnifacData;
//Flory-Huggins data for a binary mixture. alternatives for using Hildebrand,Hansen, or Chi
typedef struct {enum FF_ActModel model;enum FF_IntParamForm formula;enum FF_Units units;double deltaHil[2],deltaHan[2][3],chiData[3],V[2];}FF_FloryData;
//Activity data for a single substance
typedef struct{double lnGammaC,lnGammaSG,lnGammaR,gamma,dgEC,dgESG,dgER;}FF_SubsActivityData;
typedef struct{double gEC,gESG,gER,gE,dgE_dT;}FF_ExcessData;
//Thermodynamic model for a mixture. Not used till now
typedef struct{int numSubs,thModelActEos,actModel,refVpEos,eosType,mixRule;}FF_ThermoModel;
//Data for a mixture. Includes as arrays(fixed to 15 substances) the substances data, plus definition of the thermo model to use. BIP ar not included
//as they would change depending on the model choosen for calculations(activity/eos) and the EOS selected
//thModelActEos: 0(gamma-phi), 1(phi-phi), 2(gamma-gamma). If an  activity model is used for the liquid phase the BIP will be for this activity model.
typedef struct {char name[30],description[150],subsName[15][30],CAS[15][22];int model,numSubs,thModelActEos,actModel,refVpEos,eosType,mixRule,intForm,id[15];
                double refT,refP,intParam[15][15][6];FF_BaseProp baseProp[15];FF_UnifacData unifStdData,
                unifPSRKData,unifDortData,unifNistData;FF_SinglePointData RI[15],cp0[15],vp[15],hVsat[15],lCp[15],lDens[15],lVisc[15],lThC[15],lSurfT[15],
                gVisc[15],gThC[15],sDens[15],sCp[15];FF_CubicEOSdata cubicData[15];FF_SaftEOSdata saftData[15];FF_SWEOSdata swData[15];
                FF_Correlation cp0Corr[15],vpCorr[15],btCorr[15],hVsatCorr[15],lCpCorr[15],lDensCorr[15],lViscCorr[15],lThCCorr[15],lSurfTCorr[15],
                gDensCorr[15],gViscCorr[15],gThCCorr[15],sDensCorr[15],sCpCorr[15];}FF_MixData;
//Thermodynamic properties records
typedef struct {double MW,T,P,V,A,G,S,U,H,dP_dT,dP_dV,Cv,Cp,SS,JT,IT;}FF_ThermoProperties;//smaller record mainly for ideal gas calculations
typedef struct {double fraction,MW,T,P,V,A,G,S,U,H,dP_dT,dP_dV,Cv,Cp,SS,JT,IT,ArrDer[6],c[15],subsPhi[15];}FF_PhaseThermoProp;//fraction stands for (phase moles)/(total mix moles)
typedef struct {double MW,T,P;FF_PhaseThermoProp phase[4];}FF_MixThermoProp;
typedef struct {double uE,hE,gE;}ExcessProp;
//Data for a defined composition,T and P, with log of fugacity/pressure as an aid to speed up calculations , if needed
typedef struct {FF_MixData *mix;double T,P,H,S,z[15],logSubstFugacity[15];}FF_FeedData;

typedef struct {enum FF_CorrEquation eq ;double Tc,Pc,rhoC;unsigned nPoints;double x[40],y[40];}FF_CorrelationData;//The dimension of arrays must be declared in structures
typedef struct {enum FF_EOS eos ;double MW,Tc,Pc,Zc,w,VdWV,mu,xp,m,chi,ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints;double points[40][3];}FF_EOSPvRhoData;//The dimension of arrays must be declared in structures
typedef struct{int eosType;FF_SaftEOSdata *eos;double xp,ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints,nVpPoints,nLdPoints;double vpPoints[30][2],ldPoints[30][3];}FF_SAFTFitData;
typedef struct{int eosType;FF_CubicEOSdata *eos;double ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints;double points[40][3];}FF_CubicFitData;

#endif /* FFBASIC_H */

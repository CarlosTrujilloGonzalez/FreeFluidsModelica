/*
 * FFeosMix.h
 *
 *  Created on: 14/07/2013
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2018  Carlos Trujillo Gonzalez

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

// contains EOS definitions for pure substances
#ifndef FFEOSMIX_H
#define FFEOSMIX_H

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
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFactivity.h"


#ifdef __cplusplus
extern "C"
{
#endif

//Mixture calculations
//====================

//Get mix data from an exported file
EXP_IMP FF_MixData * CALLCONV FF_MixDataFromFile(const char *name);

//Create mixture data structure from an array of substance data structures
EXP_IMP FF_MixData * CALLCONV FF_MixDataFromSubsData(int numSubs,const FF_SubstanceData *subsData[]);


//Fill a Mixture data structure from an array of substance data structures
EXP_IMP void CALLCONV FF_MixFillDataWithSubsData2(int numSubs,FF_SubstanceData subsData[], const char *path, FF_MixData *mixData);


//Fill a Mixture data structure from an array of substance data structures
EXP_IMP void CALLCONV FF_MixFillDataWithSubsData(int *numSubs,FF_SubstanceData *subsData[],FF_MixData *mixData);

//Write a mixture data to a file. Adds ".md" extension
EXP_IMP void CALLCONV FF_MixDataToFile(const char *name,FF_MixData *mix);

//Mixture cubic EOS calculations
//-------------------------------

//Calculates Theta,b,c and their composition derivatives, given a cubic EOS,a mixing rule, composition, and pure substance parameters
void CALLCONV FF_MixParamXderCubicEOS(const int *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[15][15][6],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[]);

//Calculates Theta,b,delta and epsilon for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
EXP_IMP void CALLCONV FF_MixParamTderCubicEOS(const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
                                     const double pintParam[15][15][6],const double x[], FF_CubicParam *param);

//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i] and db/dX[i] for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
EXP_IMP void CALLCONV FF_MixParamCubicEOSOld(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[]);

//Calculates Theta,b and c, given a cubic EOS,excess G, composition, and pure substance parameters
//------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_MixParamCubicEOSgE(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param);

//Calculates Theta,b,c and their composition derivatives, given a cubic EOS,excess G, composition, and pure substance parameters
EXP_IMP void CALLCONV FF_MixParamNderCubicEOSgE(const FF_MixData *mix,const double *T,const double x[],FF_CubicParam *param,double dNalpha_dNi[],double dNb_dNi[]);

//Calculates Theta,b,c and their temperature derivatives, given a cubic EOS,excess G, composition, and pure substance parameters
EXP_IMP void CALLCONV FF_MixParamTderCubicEOSgE(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param);

//Calculates Theta,b,delta and epsilon for a mixture, given a cubic EOS,a gE mixing rule, composition,pure substance parameters, and gE
//EXP_IMP void CALLCONV calcMixParamCubicEOSgE(const enum FF_EOS *eos,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicParam sParam[],
//        const double *gE,const double x[], FF_CubicParam *param,double dNb_dNi[]);


//Mixture SAFT EOS calculation
//------------------------------
//Mixture Z and Arr calculation for a mixture, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixArrZfromTVSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                          FF_SaftEOSdata data[],const double pintParam[15][15][6],const double x[],double *Arr,double *Z);
//Mixture P calculation given T, V, and composition according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixPfromTVSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[15][15][6],const double x[],double *P);
//Mixture V,Arr and Z  calculation, given T and P and composition, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixVfromTPSAFTOld(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                                     const double pintParam[15][15][6],const double x[],const char *option,double resultL[3],double resultG[3],char *state);

//Mixture V,Arr and Z  calculation, given T and P and composition, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixVfromTPSAFT(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                                     const double pintParam[15][15][6],const double x[],const char *option,double resultL[3],double resultG[3],char *state);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a mixture, given T and V, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_MixArrDerSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                              const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double result[6]);

//Mixture common calculations
//---------------------------
//Mixture P calculation from T and V by eos
EXP_IMP void CALLCONV FF_MixPfromTVeos(const FF_MixData *mix,const double *T,const double *V,const double x[],double *P);
//Mixture V calculation from T and P by eos
EXP_IMP void CALLCONV FF_MixVfromTPeos(const FF_MixData *mix,const double *T,const double *P,const double x[],
                                       const char *option,double resultL[3],double resultG[3],char *state);
//Mixture Ideal gas thermodynamic properties calculation, from a reference state, specified by T and P, where H and S are 0
EXP_IMP void CALLCONV FF_MixIdealThermoEOS(const int *numSubs,const  FF_Correlation cp0[], const FF_BaseProp baseProp[],const double x[],double *refT,double *refP, FF_ThermoProperties *th0);

//Mixture Residual thermodynamic properties calculation from T and V, using EOS
EXP_IMP void CALLCONV FF_MixResidualThermoEOS(FF_MixData *mix,FF_PhaseThermoProp *thR);

//Mixture thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
EXP_IMP void CALLCONV FF_MixThermoEOS(FF_MixData *mix,double *refT,double *refP, FF_PhaseThermoProp *th);

//Mixture fugacity coeff. calculation from T and P by eos
EXP_IMP void CALLCONV FF_MixPhiEOS(const FF_MixData *mix,const double *T,const double *P,const double x[],const char *option,double phi[]);

//Fast Phi computation for a cubic EOS supplying parameters and derivatives
EXP_IMP void CALLCONV FF_MixPhiEOScubic(const FF_MixData *mix,FF_CubicParam *param,double da_di[],double db_di[],double dc_di[],
                                     const double *T,const double *P,const double x[],const char *option,double phi[]);


#ifdef __cplusplus
}
#endif

#endif // FFEOSMIX_H


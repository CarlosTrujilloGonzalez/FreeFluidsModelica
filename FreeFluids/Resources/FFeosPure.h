/*
 * FFeosPure.h
 *
 *  Created on: 14/07/2013
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2023  Carlos Trujillo Gonzalez

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

#ifndef FFEOSPURE_H
#define FFEOSPURE_H

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
#include "FFphysprop.h"

#ifdef __cplusplus
extern "C"
{
#endif


//Interface with permanent storage, plus reference set
//----------------------------------------------------
//Get substance data from an exported file
EXP_IMP FF_SubstanceData * CALLCONV FF_GetSubsDataFromFile(const char *name);
//Load substance data from an exported file
EXP_IMP int CALLCONV FF_LoadSubsFromFile(FF_SubstanceData *subsData, const char *name);
//Write a substance data to a file. Adds ".sd" extension
EXP_IMP void CALLCONV FF_SubsDataToFile(const char *name,FF_SubstanceData *subsData);
//Sets the reference enthalpy and pressure
EXP_IMP void CALLCONV FF_SetReference(int refState, double refT, double refP,FF_SubstanceData *subsData);

//Cubic EOS calculations
//----------------------
//Calculates u,w,b a cubic EOS. Using critical and others constants
EXP_IMP void CALLCONV FF_FixedParamCubic(const  FF_CubicEOSdata *data, FF_CubicParam *param);
//Calculates Theta and its derivatives, given a cubic EOS and T. Using critical and others constants
EXP_IMP void CALLCONV FF_ThetaDerivCubic(const double *T,const  FF_CubicEOSdata *data, FF_CubicParam *param);
//Arr and dArr/dV at constant T calculation for a pure substance, given T and V, according to cubic EOS
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to cubic EOS
EXP_IMP void CALLCONV FF_ArrDerCubic(double T,double V,const  FF_CubicParam *param,double result[6]);
//Arr and dArr/dT as for cubic EOS
EXP_IMP void CALLCONV FF_ArrDerCubic0T(double T,double V,const  FF_CubicParam *param,double result[2]);
EXP_IMP void CALLCONV FF_ArrZfromTVcubic(double T,double V,const  FF_CubicParam *param,double *Arr,double *Z);
//P calculation from T and V using cubic eos
EXP_IMP void CALLCONV FF_PfromTVcubic(double T,double V,const FF_CubicParam *param,double *P);
//V calculation for a pure substance, given T and P, according to cubic EOS
EXP_IMP void CALLCONV FF_VfromTPcubic(double T,double P,const  FF_CubicParam *param,char option,double resultL[3],double resultG[3],char *state);


//SAFT type EOS calculation
//-------------------------
//Auxiliary calculation for ArrZfromTVPCSAFT
EXP_IMP void CALLCONV FF_calcI1I2(double m,double eta,double I[4]);
//Arr and Z calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_ArrZfromTVSAFT(double T,double V,const  FF_SaftEOSdata *data,double *Arr,double *Z);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_ArrDerSAFT(double T,double V,const  FF_SaftEOSdata *data,double result[6]);
//Arr and its partial derivatives symbolic calculation for a pure substance, given T and V, according to PCSAFT EOS. To be implemented
EXP_IMP void CALLCONV FF_ArrDerPCSAFT(double T,double V,const  FF_SaftEOSdata *data,double result[6]);
//Arr and dArr/dT as per SAFT EOS
EXP_IMP void CALLCONV FF_ArrDerSAFT0T(double T,double V,const  FF_SaftEOSdata *data,double result[2]);
//P calculation from T and V using according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_PfromTVSAFT(double T,double V,const  FF_SaftEOSdata *data,double *P);
//V,Arr and Z calculation for a pure substance, given T and P, according to FF_PCSAFT EOS
//Volume solver from T and P using a SAFT or SW EOS. Regula Falsi solver, Anderson-Bjork modification.
EXP_IMP double CALLCONV FF_VsolverRegula(void *data, double T, double P, double a, double b);
EXP_IMP void CALLCONV FF_VfromTPsaft(double T,double P,const  FF_SaftEOSdata *data,char option,double resultL[3],double resultG[3],char *state);


//Span and Wagner EOS calculation
//-------------------------------
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given tau and delta, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSW(double tau,double delta,const  FF_SWEOSdata *data,double result[6]);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSWTV(double T,double V,const  FF_SWEOSdata *data,double result[6]);
//Arr and dArr/dT for a pure substance, given tau and delta, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSW0T(double tau,double delta,const  FF_SWEOSdata *data,double result[2]);
//Arr and Z calculation for a pure substance, given T and V, according to Span and Wagner EOS
EXP_IMP void CALLCONV FF_ArrZfromTVsw(double T,double V,const  FF_SWEOSdata *data,double *Arr,double *Z);
//P calculation from T and V using according to Span and Wagner EOS
EXP_IMP void CALLCONV FF_PfromTVsw(double T,double V,const  FF_SWEOSdata *data,double *P);
//P and dP_ddelta calculation for a pure substance, given reduced T and reduced rho, according to SW EOS
EXP_IMP void CALLCONV FF_PresDerSW(double tau,double delta,const  FF_SWEOSdata *data,double result[4]);
//V,Arr and Z calculation for a pure substance, given T and P, according to SW EOS
EXP_IMP void CALLCONV FF_VfromTPsw(double T,double P,const  FF_SWEOSdata *data,char option,double resultL[3],double resultG[3],char *state);
//V,Arr and Z calculation for a pure substance, given T and P, according to SW EOS
EXP_IMP void CALLCONV FF_VfromTPswS(double T,double P,const  FF_SubstanceData *data,char option,double resultL[3],double resultG[3],char *state);

//Common P,V,T calculations
//-------------------------

//P calculation from T and V by eos
EXP_IMP void CALLCONV FF_PfromTVeos(int eosType,double T,double V,void *data,double *P);
//P calculation from T and V by eos
EXP_IMP void CALLCONV FF_PfromTVeosS(double T,double V,const FF_SubstanceData *data,double *P);
//V,Arr and Z calculation for a pure substance, given T and P by eos
EXP_IMP void CALLCONV FF_VfromTPeos(int eosType,double T,double P,const void *data,char option,double resultL[3],double resultG[3],char *state);
//V,Arr and Z calculation for a pure substance, given T and P by eos.
EXP_IMP void CALLCONV FF_VfromTPeosS(double T,double P,const FF_SubstanceData *data,char option,double resultL[3],double resultG[3],char *state);

//Properties at saturation. We use p and T in order not to perform again the calculation of Vp
//------------------------
//Boiling point calculation
EXP_IMP void CALLCONV FF_TbEos(int eosType,double P,const void *data,double *Tb);
//Boiling point calculation
EXP_IMP void CALLCONV FF_TbEosS(double P,const FF_SubstanceData *data,double *Tb);
//Vapor pressure calculation
EXP_IMP void CALLCONV FF_VpEos(int eosType,double T,const void *data,double *Vp);
//Vapor pressure calculation
EXP_IMP void CALLCONV FF_VpEosS(double T,const FF_SubstanceData *data,double *Vp);

//Thermodynamic properties calculation from T and V
//-------------------------------------------------
//Thermodynamic properties calculation for a ideal gas at same T and V, from a reference state, specified by refT and refP, where H and S are 0
EXP_IMP void CALLCONV FF_IdealThermoEos(int equation,const double coef[],double refT,double refP, FF_ThermoProperties *th0);
//Ideal gas thermodynamic properties of water calculation , from a reference state specified by the triple point where H and S are 0
EXP_IMP void CALLCONV FF_IdealThermoWater( FF_ThermoProperties *th0);
//Enthalpy and entropy calculation from T,V and P using EOS
EXP_IMP void CALLCONV FF_HSfromTVPeosS(double T, double V, double P, const FF_SubstanceData *data, double *H, double *S);
//Residual extended thermodynamic properties calculation from T and V, using EOS
EXP_IMP void CALLCONV FF_ExtResidualThermoEosS(const FF_SubstanceData *data, FF_ThermoProperties *thR);
//Thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
EXP_IMP void CALLCONV FF_ThermoEosS(const FF_SubstanceData *data, FF_ThermoProperties *th);


//Phase solution and properties
//-----------------------------
//Phases and T,V solver from P and T or H or U or S
EXP_IMP void CALLCONV FF_TVsFromPXNewton(char var, FF_SubstanceData *data, int aid, double p, double x, double *T, double *Vg, double *Vl, double *gf);
//Phases and T,V solver from V and T
EXP_IMP void CALLCONV FF_TVsFromVX(char var, FF_SubstanceData *data, int aid, double v, double x, double *T, double *Vg, double *Vl, double *gf);
//Phases and T,D solver from P and H or S
EXP_IMP void CALLCONV FF_TVfromPX(char var, FF_SubstanceData *data, double p, double x, double *T, double *gd, double *ld);
//all thermo properties from T/P or T/D or P/H or P/S
EXP_IMP void CALLCONV FF_solveEos(char *variable, FF_SubstanceData *data,int aid, double x, double y, double *T, double *p, double *gd, double *gh, double *gs, double *gCv, double *gCp, double *gDvp,
                          double *gDvT, double *ld, double *lh, double *ls, double *lCv, double *lCp, double *lDvp, double *lDvT, double *gf);
//Phases thermo properties calculation
EXP_IMP void CALLCONV FF_solveEosOld(char *variable, FF_SubstanceData *data, double x, double y, double *T, double *p, double *gd, double *gh, double *gs, double *gCv, double *gCp, double *gDvp,
                          double *gDvT, double *ld, double *lh, double *ls, double *lCv, double *lCp, double *lDvp, double *lDvT);

#ifdef __cplusplus
}
#endif
#endif /* FFEOSPURE_H */

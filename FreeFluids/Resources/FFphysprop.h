/*
 * FFphysprop.h
 *
 *  Created on: 26/12/2015
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

// contains mainly definitions for pure substances physical properties calculations, by non-eos methods
#ifndef FFPHYSPROP
#define FFPHYSPROP
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

#include "FFbasic.h"

#ifdef __cplusplus
extern "C"
{
#endif
//Single substance, physical properties correlations
//--------------------------------------------------
//Calculates the result of the given equation
EXP_IMP void CALLCONV FF_CorrelationResult(int eq,const double coef[],int nPoints,double x[],double y[]);

//Calculates physical property, with input and output in SI units(kgr, not moles), using the given correlation, that may or not be in SI units
EXP_IMP void CALLCONV FF_PhysPropCorr(int cor,const double coef[],double MW,int nPoints,double x[],double y[]);

//Calculates the result of the given equation
//Same as previous, but now without an array of values
EXP_IMP void CALLCONV FF_CorrelationResultM(int eq,const double coef[],double x,double *y);

//Calculates physical property, with input and output in SI units(kgr, not moles), using the given correlation, that may or not be in SI units
//Same as previous, but now without an array of values
EXP_IMP void CALLCONV FF_PhysPropCorrM(int cor,const double coef[],double MW,double x,double *y);

//Finds independent variable of a physical property correlation
EXP_IMP int CALLCONV FF_CorrelationSolver(double y, FF_Correlation *corr, double MW, double *x, double *yFound);

//Calculates specific enthalpy from a Cp correlation, with reference T=0 K
EXP_IMP void CALLCONV FF_SpecificEnthalpyCorr(int cor,const double coef[],double MW, double x, double *H);

//Calculates specific entropy from a Cp correlation, with reference T=0 K
EXP_IMP void CALLCONV FF_SpecificEntropyCorr(int cor,const double coef[],double MW, double x, double *S);

//Calculates specific enthalpy and entropy from a Cp correlation, with reference T=0 K
EXP_IMP void CALLCONV FF_SpecificEnthalpyEntropyCorr(int cor,double coef[],double MW,int nPoints,double x[],double H[],double S[]);

//Vapor pressure related calculations
//-----------------------------------
//Acentric factor calculation by the definition equation
EXP_IMP void CALLCONV FF_WfromDefinition(double Pc,double Vp,double *w);

//Acentric factor calculation from one vapor pressure
EXP_IMP void CALLCONV FF_WfromOneVap(double Tc,double Pc,double T,double Vp,double *w);

//Vapor pressure using Ambrose-Walton equation. Needs only Tc,Pc and w
EXP_IMP void CALLCONV FF_VpAmbroseWalton(FF_BaseProp *baseProp,double T,double *Vp);

//Vapor pressure using Riedel-Vetere equation. Needs only Tc,Pc and one boiling point
EXP_IMP void CALLCONV FF_VpRiedelVetere(FF_BaseProp *baseProp,double Tref,double VpRef,double T,double *Vp);

//Vapor pressure calculation by correlations
EXP_IMP void CALLCONV FF_Vp(double T,FF_SubstanceData *data,double *Vp);

//Liquid density calculations
//---------------------------
//Rackett equation for satured liquid density. If you supply ref rho, Vc is not used. It is better to use w than Zra, and Zra than Zc
EXP_IMP void CALLCONV FF_LiqDensSatRackett(FF_BaseProp *baseProp,double Tref,double rhoRef,double T,double *rho);

//Chueh-Prausnitz pressure correction for liquid density
EXP_IMP void CALLCONV FF_LiqDensChuehPrausnitz(FF_BaseProp *baseProp, double T,double P,double Pin,double rhoIn, double *rho);

//Tait equation for polymer density, with T and P dependence
EXP_IMP void CALLCONV FF_LiqDensTait(int eq,double coef[],FF_BaseProp *baseProp,int nPoints,double T[],double P[], double rho[]);

//Liquid density from T,P
EXP_IMP void CALLCONV FF_LiqDensTP(double T,double P,FF_SubstanceData *data,double *liqDens);

//Others thermodynamic properties
//-------------------------------
//Liquid Cp. Bondi method
void CALLCONV FF_LiqCpBondi(FF_SubstanceData *data,double T,double *Cp);

//Transport properties of liquids
//-------------------------------
//Lucas liquid viscosity pressure correction.
EXP_IMP void CALLCONV FF_LiqViscPcorLucas(double T,double P,double Pref,FF_BaseProp *data,double rpVisc,double *visc);

//Liquid viscosity from T,P
EXP_IMP void CALLCONV FF_LiqViscTP(double T,double P,FF_SubstanceData *data,double *liqVisc);

//Grunberg-Nissan method for mixture liquid viscosity
//BIPs are missing
EXP_IMP void CALLCONV FF_MixLiqViscGrunberg(FF_MixData *mix,double T,double P,double x[],double *visc);

//Teja-Rice method for mixture liquid viscosity
//BIP is missing
EXP_IMP void CALLCONV FF_MixLiqViscTeja(FF_MixData *mix,double *T,double *P,double x[],double *visc);

//Andrade method for mixture liquid viscosity
EXP_IMP void CALLCONV FF_MixLiqViscAndrade(FF_MixData *mix,double T,double P,double x[],double *visc);

//Thermal conductivity of liquids. Latini method
EXP_IMP void CALLCONV FF_LiquidThCondLatini(double T,FF_BaseProp *data,double *thCond);

//Liquid thermal conductivity from T
EXP_IMP void CALLCONV FF_LiqThCondT(double T,FF_SubstanceData *data,double *liqThCond);

//Li method for mixture liquid thermal conductivity
EXP_IMP void CALLCONV FF_MixLiqThCondLi(FF_MixData *mix,double *T,double *P,double x[],double *thCond);

//SurfaceTension, MacLeod-Sugden method. Very sensible to Parachor value
EXP_IMP void CALLCONV FF_SurfTensMcLeod(double T,FF_SubstanceData *data,double *surfTens);

//Surface tension,Sastri-Rao method
EXP_IMP void CALLCONV FF_SurfTensSastri(double T,FF_BaseProp *data,double *surfTens);

//Surface tension from T
EXP_IMP void CALLCONV FF_SurfTensT(double T,FF_SubstanceData *data,double *surfTens);

//Linear method for mixture surface tension
EXP_IMP void CALLCONV FF_MixLiqSurfTensLinear(FF_MixData *mix,double T,double P,double x[],double *surfTens);

//Winterfeld method for mixture surface tension.
EXP_IMP void CALLCONV FF_MixLiqSurfTensWinterfeld(FF_MixData *mix,double T,double P,double x[],double *surfTens);

//McLeod-Sugden method for mixture surface tension.
EXP_IMP void CALLCONV FF_MixLiqSurfTensMcLeod(FF_MixData *mix,double rhoL,double rhoG,double x[],double y[],double *surfTens);

//Transport properties of gases
//-----------------------------
//Gas viscosity pressure prediction/correction. Chung method
EXP_IMP void CALLCONV FF_GasViscTVcpChung(double T,double V,FF_BaseProp *data,double *ldVisc,double *visc);

//Gas viscosity pressure prediction/correction. Lucas method
EXP_IMP void CALLCONV FF_GasViscTPcpLucas(const double T,const double P,const FF_BaseProp *data,double *lpVisc,double *visc);

//Gas viscosity from T,P
EXP_IMP void CALLCONV FF_GasViscTP(double T,double P,FF_SubstanceData *data,double *gasVisc);

//Viscosity of gas mixtures. Wilke method
EXP_IMP void CALLCONV FF_MixGasViscTPcpWilke(FF_MixData *mix,double T,double P,double y[],double *gVisc);

//Viscosity of gas mixtures. Lucas method
EXP_IMP void CALLCONV FF_MixGasViscTPcpLucas(FF_MixData *mix,double T,double P,double y[],double *gVisc);

//Gas low pressure thermal conductivity prediction
EXP_IMP void CALLCONV FF_GasLpThCondTCpChung(double T,double Cp0,FF_BaseProp *data,double *ldThCond);

//Gas thermal conductivity pressure correction
EXP_IMP void CALLCONV FF_GasThCondTVcorChung(double T,double V,FF_BaseProp *data,double *ldThCond,double *thCond);

//Gas thermal conductivity from T,V
EXP_IMP void CALLCONV FF_GasThCondTV(double T,double V,FF_SubstanceData *data,double *gasThCond);

//Thermal conductivity of low pressure gas mixtures. Mason and Saxena method
EXP_IMP void CALLCONV FF_MixLpGasThCondTpMason(FF_MixData *mix,double T,double y[],double *gThCond);

//Transport properties phase independent
//--------------------------------------
//Temperature and density dependent, phase independent, viscosity.
EXP_IMP void CALLCONV FF_ViscosityTDens(int subsRef, double T, double rhoMolar,double eta[3]);
//Viscosity of a pure substance

EXP_IMP void CALLCONV FF_Viscosity(FF_SubstanceData *subs,FF_SubstanceData *ref,double T,double dens, double P,double gf,double *visc);
//Temperature and density dependent, phase independent, thermal conductivity and viscosity.

EXP_IMP void CALLCONV FF_ThCondViscTDens(FF_SubstanceData *data, double T, double rhoMolar,double lambda[3],double eta[3]);

//Thermal conductivity of a pure substance
EXP_IMP void CALLCONV FF_ThCond(FF_SubstanceData *subs,FF_SubstanceData *ref,double T,double dens, double P,double gf,double *thCond);

//Conformal state calculation using SW EOS. Returns the reducing ratios f and h (h refered to molar units)
EXP_IMP void CALLCONV FF_ConformalStateSW(FF_SubstanceData *subs, FF_SubstanceData *ref, double T, double rhoMolar, double *f, double *h);

//Conformal state calculation of a pure subs. or a mixture, using SRK EOS. Returns the reducing ratios f and h (h refered to molar units)
EXP_IMP void CALLCONV  FF_ConformalStateSRK(FF_SubstanceData *subs,FF_SubstanceData *ref, double T,double *f, double *h);

//Conformal state calculation at saturation. Returns the reducing ratios f and h (h refered to molar units)
EXP_IMP void CALLCONV FF_ConformalStateSat(FF_SubstanceData *subs, FF_SubstanceData *ref, double T, char option, double *f, double *h);

//Will find the reducing ratios by the best method
EXP_IMP void CALLCONV FF_ConformalState(FF_SubstanceData *subs, FF_SubstanceData *ref, double T, double dens, double *f, double *h);

//Viscosity and thermal conductivity using extended corresponding states
EXP_IMP void CALLCONV FF_CorrespondingStatesTransport(FF_SubstanceData *subs, FF_SubstanceData *ref, double T, double dens, double f, double h,
                                                      double eta[3], double lambda[3]);

//Liquid viscosity and thermal conductivity at saturation using extended corresponding states
EXP_IMP void CALLCONV FF_CorrespondingStatesSat(FF_SubstanceData *subs, FF_SubstanceData *ref, double T, double f, double h,
                                                double *eta, double *lambda);

#ifdef __cplusplus
}
#endif
#endif // FFPHYSPROP


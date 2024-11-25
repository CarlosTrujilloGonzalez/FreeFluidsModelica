/*

 * FFequilibrium.h
 *
 *  Created on: 31/07/2018
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
#ifndef FFEQUILIBRIUM_H
#define FFEQUILIBRIUM_H

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
#include "FFeosMix.h"

#ifdef __cplusplus
extern "C"
{
#endif
//Calculation of the minimum tangent plane distance and the corresponding composition
EXP_IMP void CALLCONV FF_StabilityCheck(FF_FeedData *data,double *tpd,double tpdX[]);
//Calculation of the minimum tangent plane distance and the corresponding composition, by simulated annealing
void CALLCONV FF_StabilityCheckSA(FF_FeedData *data,double *tpd,double tpdX[]);

//Mixture bubble temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleT(FF_MixData *mix,const double *P,const double x[],const double *bTguess, double *bT,double y[],
                                 double substPhiL[],double substPhiG[]);
//Mixture dew temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewT(FF_MixData *mix,const double *P, const double y[],const double *dTguess, double *dT,double x[],
                              double substPhiL[],double substPhiG[]);
//Temperature envelope of a binary mixture
void CALLCONV FF_TemperatureEnvelope(FF_MixData *mix,const double *P, const int *nPoints, double c[],double bT[],double y[],double dT[],double x[]);
//Mixture bubble pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleP(FF_MixData *mix,const double *T, const double c[],const double *bPguess, double *bP,
                                 double y[],double substPhiL[],double substPhiG[]);
//Mixture dew pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewP(FF_MixData *mix,const double *T, const double y[],const double *dPguess, double *dP,double x[],double substPhiL[],double substPhiG[]);
//Pressure envelope of a binary mixture
EXP_IMP void CALLCONV FF_PressureEnvelope(FF_MixData *mix,const double *T, const int *nPoints, double c[], double bP[],double y[],double dP[],double x[]);

//VL flash calculation, given T, P, feed composition, eos and mixing rule, using bubble and dew pressure as help
EXP_IMP void CALLCONV FF_TwoPhasesPreFlashPT(FF_MixData *mix,const double *T,const double *P,const double f[],
                                     double x[],double y[],double substPhiL[],double substPhiG[],double *beta);

//VL flash calculation, given T, P, feed composition, eos and mixing rule
EXP_IMP void CALLCONV FF_TwoPhasesFlashPT(FF_MixData *mix,const double *T,const double *P,const double f[],
                                   double x[],double y[],double substPhiL[],double substPhiG[],double *beta);
//Mixture VL flash, given P,T, composition, and thermo model to use. By global optimization simulated annealing of residual Gibbs energy
EXP_IMP void CALLCONV FF_TwoPhasesFlashPTSA(FF_FeedData *data, double x[],double y[],double substPhiL[],double substPhiG[],double *beta,double *Gr);
//Mixture 2 phases flash, given P,T, composition, and thermo model to use. By differential evolution global minimization of the reduced Gibbs energy
void CALLCONV FF_TwoPhasesFlashPTDE(FF_FeedData *data, double x[],double y[],double substPhiB[],double substPhiA[],double *beta,double *Gr);

//VL flash calculation, given h, P, feed composition, eos and mixing rule, using bubble and dew pressure as help
EXP_IMP void CALLCONV FF_TwoPhasesPreFlashPH(FF_MixData *mix, const double *H,const double *P,const double f[], double *T,
                                     double x[],double y[],double substPhiL[],double substPhiG[],double *beta);

//Mixture VL flash, given P,T, composition, and thermo model to use. By simulated annealing global optimization of the reduced Gibbs energy
EXP_IMP void CALLCONV FF_ThreePhasesFlashPTSA(FF_FeedData *data, double x[],double y[],double z[],double substPhiA[],double substPhiB[],double substPhiC[],
                                      double *betaA,double *betaB,double *Gr);
#ifdef __cplusplus
}
#endif
#endif // FFEQUILIBRIUM_H

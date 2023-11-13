/*
 * FFactivity.h
 *
 *  Created on: 10/04/2017
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


#ifndef FFACTIVITY_H
#define FFACTIVITY_H

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFphysprop.h"


#ifdef __cplusplus
extern "C"
{
#endif
//Calculates activity coefficients according to Wilson (modified) equation, at given T and composition
EXP_IMP void CALLCONV FF_ActivityWilson(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                        const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]);

//Calculates activity coefficients according to NRTL equation, at given T and composition
EXP_IMP void CALLCONV FF_ActivityNRTL(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                      const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]);

//Calculates activity coefficients according to UNIQUAQ equation, at given T and composition
EXP_IMP void CALLCONV FF_ActivityUNIQUAC(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                 const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]);

//Calculates activity coefficients in polymer systems, using FV for combinatorial and segment UNIQUAQ for residual, at given T and composition
//For polymers x and FV are that of the monomer unit. In numMon we supply the number of monomeric units of the chain, must be <1 for solvents
EXP_IMP void CALLCONV FF_ActivityUNIQUACFV(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]);

//Flory-Huggins. Temperature dependent chi? Use FV in combinatorial part?
EXP_IMP void FF_ActivityFloryHuggins(const int *model,const  FF_BaseProp baseProp[2],const double chiData[15][15][6],const int *form,const double *T,const double x[2],FF_SubsActivityData actData[2]);

//Calculates the common data, independent of the molar fraction and temperature for the UNIFAC model.
void CALLCONV FF_UNIFACParams(int numData, const int data[][3], FF_UnifacData *uni);

//Calculates activity coefficients according to UNIFAC model, at given T and composition.
EXP_IMP void CALLCONV FF_ActivityUNIFAC(FF_UnifacData *data, const double *T, const double x[],FF_SubsActivityData actData[]);

//Calculates activity coefficents according any of the defined models
EXP_IMP void CALLCONV FF_Activity(FF_MixData *mix,double *T,double x[],FF_SubsActivityData actData[]);

//Calculates fugacity and activity coefficients, at given T and composition, from an activity model
EXP_IMP void CALLCONV FF_PhiAndActivity(FF_MixData *mix,const double *T,const double *P,const double x[],FF_SubsActivityData actData[],double phi[]);

//Calculates fugacity from given activity coefficients, at given T,P and composition, using an activity model
void CALLCONV FF_PhiFromActivity(FF_MixData *mix,const double *T,const double *P,const double x[],FF_SubsActivityData actData[],double phi[]);

//Calculates ln of activities and the derivatives of gE for UNIFAC
EXP_IMP void CALLCONV FF_UNIFACDerivatives(FF_UnifacData *data, const double *T, const double x[],FF_SubsActivityData actData[],FF_ExcessData *excData);

//Calculates ln of activities and the derivatives of gE for act.coef. models
EXP_IMP void CALLCONV FF_ActivityDerivatives(const int *actModel,const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],const int *form,
                                        const double *T,const double x[],FF_SubsActivityData actData[],FF_ExcessData *excData);

#ifdef __cplusplus
}
#endif

#endif // FFACTIVITY_H

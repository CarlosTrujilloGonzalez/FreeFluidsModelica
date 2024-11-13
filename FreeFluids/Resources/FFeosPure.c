/*

 * FFeosPure.c
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

// contains EOS calculations for pure substances
//==============================================

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include <string.h>
#include <unistd.h>

const double Av = 6.02214E+23; //molecules/mol
const double kb = 1.3806504E-023;//J/K
const double Pi = 3.141592;
const double R = 8.314472; //Pa·m3/(K·mol)
const double FF_PCSAFTap[7][3]={{0.910563145, -0.308401692, -0.090614835},{0.636128145, 0.186053116, 0.452784281},
        {2.686134789, -2.503004726, 0.596270073},{-26.54736249, 21.41979363, -1.724182913},
        {97.75920878, -65.25588533, -4.130211253},{-159.5915409, 83.31868048, 13.77663187},
        {91.29777408, -33.74692293, -8.672847037}};
const double FF_PCSAFTbp[7][3]={{0.724094694, -0.575549808, 0.097688312},{2.238279186, 0.699509552, -0.255757498},
        {-4.002584948, 3.892567339, -9.155856153},{-21.00357681, -17.21547165, 20.64207597},
        {26.85564136, 192.6722645, -38.80443005},{206.5513384, -161.8264616, 93.62677408},
        {-355.6023561, -165.2076935, -29.66690559}};

//Interface with permanent storage
//================================

//Get substance data from an exported file
EXP_IMP FF_SubstanceData * CALLCONV FF_GetSubsDataFromFile(const char *name){
    FF_SubstanceData *subsData = (FF_SubstanceData*) calloc(1,sizeof(FF_SubstanceData));
    char path[FILENAME_MAX]="";
    strcat(path,name);
    strcat(path,".sd");
    FILE * file= fopen(path, "rb");
    if (file != NULL) {
        fread(subsData, sizeof(FF_SubstanceData), 1, file);
        fclose(file);
    }
    else printf("Substance data file not found\n");
    if ((subsData->baseProp.FV==0)&&(subsData->baseProp.Vliq>0)&&(subsData->baseProp.VdWV>0)) subsData->baseProp.FV=subsData->baseProp.Vliq-
            1.2*subsData->baseProp.VdWV;
    //printf("Tc:%f\n",subsData->baseProp.Tc);
    return subsData;
}

//Load substance data from an exported file
EXP_IMP int CALLCONV FF_LoadSubsFromFile(FF_SubstanceData *subsData, const char *name){
    char path[FILENAME_MAX]="";
    strcat(path,name);
    strcat(path,".sd");
    FILE * file= fopen(path, "rb");
    if (file != NULL) {
        fread(subsData, sizeof(FF_SubstanceData), 1, file);
        fclose(file);
    }
    else{
        //printf("Substance data file not found\n");
        return 1;
    }
    if ((subsData->baseProp.FV==0)&&(subsData->baseProp.Vliq>0)&&(subsData->baseProp.VdWV>0)) subsData->baseProp.FV=subsData->baseProp.Vliq-
            1.2*subsData->baseProp.VdWV;
    //printf("Tc:%f\n",subsData->baseProp.Tc);
    return 0;
}

//Write a substance data to a file. Adds ".sd" extension
EXP_IMP void CALLCONV FF_SubsDataToFile(const char *name,FF_SubstanceData *subsData){
    char path[FILENAME_MAX]="Data/";
    strcat(path,name);
    strcat(path,".sd");
    FILE * file= fopen(path, "wb");
    if (file != NULL) {
        fwrite (subsData, sizeof(FF_SubstanceData), 1, file);
        fclose(file);
    }
    else printf("Error open file\n");
}

//Set of references
//=================

//Sets the reference enthalpy and pressure
void CALLCONV FF_SetReference(int refState, double refT, double refP,FF_SubstanceData *subsData){
    double T,P,V,answerL[3],answerG[3],H,S;
    char option,state;
    if(refState==1){//ASHRAE
        T=233.15;
        if(subsData->vpCorr.form>0) FF_PhysPropCorrM(subsData->vpCorr.form,subsData->vpCorr.coef,subsData->baseProp.MW,T,&P);
        else FF_VpEosS(T,subsData,&P);
    }
    else if(refState==2){//IIR
        T=273.15;
        if (subsData->id==399) P=611;
        else if(subsData->vpCorr.form>0) FF_PhysPropCorrM(subsData->vpCorr.form,subsData->vpCorr.coef,subsData->baseProp.MW,T,&P);
        else FF_VpEosS(T,subsData,&P);
    }
    else if(refState==3){//NBP
        P=101325;
        double pFound;
        if(subsData->vpCorr.form>0) FF_CorrelationSolver(P,&subsData->vpCorr,subsData->baseProp.MW,&T,&pFound);
        else FF_TbEosS(P,subsData,&T);
    }
    else {//User defined as default
        T=refT;
        P=refP;
    }
    subsData->refT=0.0;//ideal gas calculation will be always from 0K
    subsData->refP=101325;//ideal gas entropy referenced always to 1 atm
    if((refState==1)||(refState==2)||(refState==3)){
        option='l';
        FF_VfromTPeosS(T,P,subsData,option,answerL,answerG,&state);
        V=answerL[0];
        FF_HSfromTVPeosS(T,V,P,subsData,&subsData->refH,&subsData->refS);
        if(refState==2){//IIR
            subsData->refH=subsData->refH-2.0e2*subsData->baseProp.MW;
            subsData->refS=subsData->refS-1*subsData->baseProp.MW;
        }
    }
    else{
        if(refT==0){
            subsData->refH=0;
            subsData->refS=0;
        }
        else{
            option='s';
            FF_VfromTPeosS(T,P,subsData,option,answerL,answerG,&state);
            if(state=='L') V=answerL[0];
            else V=answerG[0];
            FF_HSfromTVPeosS(T,V,P,subsData,&subsData->refH,&subsData->refS);
        }
    }
}

//Cubic EOS calculations
//======================

//Calculates a,b,u,w for a cubic EOS. Using critical constants
//------------------------------------------------------------
void CALLCONV FF_FixedParamCubic(const  FF_CubicEOSdata *data, FF_CubicParam *param){
    switch (data->eos)
    {
    case FF_PR76:
    case FF_PR78://Peng-Robinson.
    case FF_PRSV1://Peng-Robinson Stryjek-Vera with 1 extra parameter.
    case FF_PRBM://Peng-Robinson Boston-Mathias with 1 extra parameter. Very similar to Stryjeck-Vera
    case FF_PRMELHEM://Peng-robinson Melhem with 2 extra parameters.
    case FF_PRSOF://Peng-Robinson Stamateria-Olivera-Fuentes with 2 extra parameters
    case FF_PRALMEIDA://Peng-Robinson Almeida with 3 extra parameters.
    case FF_PRMC://Peng-Robinson Mathias-Copeman with 3 extra parameters.
    case FF_PRTWU91://Peng-Robinson Twu with 3 extra parameters.
    case FF_PRTWU95://Peng-robinson Twu with parameters from Tc and w
    case FF_PRFIT4B://Peng Robinson Stryjeck- Vera with Tc,Pc,w,k1 fitted
        //printf("Tc:%f Pc:%f w:%f Zc%f:\n",data->Tc,data->Pc,data->w,data->Zc);
        param->a = 0.457235 * pow(R*data->Tc,2)/ data->Pc;
        param->b = 0.077796 * R * data->Tc / data->Pc;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;//if we supply a volume correction it is used
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.45)) param->c=R*data->Tc*(0.1014048-0.3892896*data->Zc)/data->Pc;//if volume correction is 0, we apply that of Peneloux
        else param->c=0.0;//a negative number will indicate not to use volume correction
        param->Zc=0.307;
        break;
    case FF_PRvTWU91:
        param->a = param->a = 0.42748 * pow(R*data->Tc,2)/ data->Pc*1.08;//0.599877;
        //param->b = data->r*15.17*1.6*1e-6;//would be an alternative using UNIQUAC r parameter
        param->b = data->VdWV*1.55;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        //if (data->c>0) param->c=data->c;
        //else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        //else param->c=0;
        param->c=data->c;
        param->Zc=0.307;
        break;
    case FF_PRFIT3://similar to FF_PRSV1, some parameters have no relationship with critical properties. There are 3 adjustable parameters:a, Tc, k1
        param->a = data->k1;
        param->b = 0.9*0.0778 * R * data->Tc / data->Pc;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        else param->c=0;
        param->Zc=0.307;
        break;
    case FF_PRFIT4://parameters have no relationship with critical properties. There are 4 adjustable parameters:b,a,w and Tc
        param->a = data->k2;
        param->b = data->k1;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        else param->c=0;
        param->Zc=0.307;
        break;
    case FF_SRK://Soaves-Redlich-Kwong.
    case FF_SRKSOF://Soaves-Redlich-Kwong Stamateria-Olivera-Fuentes with 2 extra parameters
    case FF_SRKMC://Soaves-Redlich-Kwong Mathias-Copeman, with 3 extra parameters.
    case FF_SRKTWU91:
        param->a = 0.42748 * pow(R*data->Tc,2)/ data->Pc;
        param->b = 0.08664 * R * data->Tc / data->Pc;
        param->u=1;
        param->w=0;
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->w>0)&&(data->w<1)) param->c=R*data->Tc/data->Pc*(0.001569568+0.03577392*data->w);
        else param->c=0;
        param->Zc=0.3333;
        break;
    case FF_PRPOL1://Peng-robinson for polymers b=k1*MW, a=1,Theta=k2*MW (a=1)
        param->a=1;
        param->b=data->MW*data->k1;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        break;
    case FF_SRKPOL2://Soave-Redlich-Kwong for polymers b=k1*MW, a=1
        param->a=1;
        param->b=data->MW*data->k1;
        param->u=1;
        param->w=0;
        break;
    default:
        break;
    }
    param->Tc=data->Tc;
    param->Pc=data->Pc;
}

//Calculates Theta and its derivatives, given a cubic EOS and T. Using critical and others constants
//--------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_ThetaDerivCubic(const double *T,const  FF_CubicEOSdata *data, FF_CubicParam *param)
{
    double Tr,Tx,dTx,d2Tx,fw,Alpha,dAlpha,d2Alpha;
    Tr = *T / data->Tc;
    Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
    dTx=-0.5*pow(Tr,-0.5)/data->Tc;//This is the 1sr derivative redarding T
    d2Tx=0.25*pow(Tr,-1.5)/pow(data->Tc,2);//and this the 2nd
    double d;//used in FF_PRBM
    switch (data->eos)
    {
    case FF_PR76:
        fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PR78://Peng-Robinson.
        if (data->w<=0.491) fw = 0.37464 + 1.54226 * data->w - 0.26992 * pow(data->w,2);//0.491
        else fw = 0.379642 + 1.487503 * data->w - 0.164423 * pow(data->w,2)+ 0.016666 * pow(data->w,3);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PRSV1://Peng-Robinson Stryjek-Vera with 1 extra parameter.
    case FF_PRFIT4B://Peng Robinson Stryjeck- Vera with Tc,Pc,w,k1 fitted
        fw = 0.378893 + 1.4897153 * data->w - 0.17131848 * pow(data->w,2) + 0.0196554 * pow(data->w,3);
        if (*T < data->Tc){
            Alpha = pow((1 + fw * Tx+ data->k1 * (1 - Tr) * (0.7 - Tr)),2);
            dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k1*(2*Tr-1.7)/data->Tc);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k1*2/pow(data->Tc,2));
        }
        else{
            Alpha = pow((1 + fw * Tx),2);
            dAlpha=2*(1 + fw * Tx)*fw*dTx;
            d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        }
        break;
    case FF_PRBM:
        fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
        if (*T < data->Tc){
            Alpha = pow((1 + fw * Tx- data->k1 * (1 - Tr) * (0.7 - Tr)),2);
            dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k1*(2*Tr-1.7)/data->Tc);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k1*2/pow(data->Tc,2));
        }
        else{
            d=1+fw/2+0.3*data->k1;
            Alpha = pow(exp((1-1/d)*(1-pow(Tr,d))),2);
            dAlpha=-2*Alpha*(d-1)*pow(Tr,d-1)/data->Tc;
            d2Alpha=-2*(d-1)/data->Tc*(dAlpha*pow(Tr,d-1)+Alpha*(d-1)*pow(Tr,d-2)/data->Tc);
        }

        break;
    case FF_PRMELHEM://Peng-robinson Melhem with 2 extra parameters.
        Alpha =exp(data->k1*(1-Tr)+data->k2*pow(Tx,2));
        dAlpha=Alpha*(-data->k1/data->Tc+2*data->k2*Tx*dTx);
        d2Alpha=dAlpha*(-data->k1/data->Tc+2*data->k2*Tx*dTx)+Alpha*2*data->k2*(dTx*dTx+Tx*d2Tx);
        break;
    case FF_PRSOF://Stamateria-Olivera-Fuentes
    case FF_SRKSOF:
        Alpha=Tr*(1+data->k1/(data->k2-1))-data->k1/(data->k2-1)*pow(Tr,2-data->k2);
        dAlpha=(1+data->k1/(data->k2-1))/data->Tc-(2-data->k2)*data->k1/(data->k2-1)*pow(Tr,1-data->k2)/data->Tc;
        d2Alpha=-(2-data->k2)*data->k1*(1-data->k2)/(data->k2-1)*pow(Tr,-data->k2)/data->Tc/data->Tc;
        break;
    case FF_PRMC://Peng-Robinson Mathias-Copeman with 3 extra parameters.
    case FF_SRKMC://Soaves-Redlich-Kwong Mathias-Copeman, with 3 extra parameters.
        if (*T < data->Tc)
        {
            Alpha = pow(1+data->k1*Tx+data->k2*pow(Tx,2)+data->k3*pow(Tx,3),2);
            dAlpha=2*pow(Alpha,0.5)*(data->k1*dTx+data->k2*2*Tx*dTx+data->k3*3*pow(Tx,2)*dTx);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(data->k1*d2Tx+data->k2*2*(pow(dTx,2)+Tx*d2Tx)+data->k3*3*(2*Tx*pow(dTx,2)+pow(Tx,2)*d2Tx));
        }
        else
        {
            Alpha =pow(1+data->k1*Tx,2);
            dAlpha=2*pow(Alpha,0.5)*data->k1*dTx;
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*data->k1*d2Tx;
        }
        //printf("T: %f Alpha: %f dAlpha: %f d2Alpha: %f\n",*T,Alpha,dAlpha,d2Alpha);
        break;
    case FF_PRALMEIDA:
        Alpha=exp(data->k1*(1-Tr)*pow(fabs(1-Tr),data->k3-1)+data->k2*(1/Tr-1));
        dAlpha=Alpha*((-data->k1*(pow(fabs(1-Tr),data->k3-1))-(data->k1*(1-Tr))*(data->k3-1)*pow(fabs(1-Tr),data->k3-2)-data->k2/Tr/Tr)/data->Tc);
        d2Alpha=dAlpha*dAlpha/Alpha+Alpha*(2*data->k1*(data->k3-1)*pow(fabs(1-Tr),data->k3-2)-data->k1*(1-Tr)*(data->k3-1)*(data->k3-2)*pow(fabs(1-Tr),data->k3-3)+
                2*data->k2/Tr/Tr/Tr)/data->Tc/data->Tc;
        break;
    case FF_PRTWU91://Peng-Robinson Twu with 3 extra parameters.
    case FF_PRvTWU91:
    case FF_SRKTWU91://Idem SRK
        Alpha =pow(Tr,data->k3*(data->k2-1))*exp(data->k1*(1-pow(Tr,(data->k2*data->k3))));
        dAlpha=Alpha*(data->k3*(data->k2-1)/Tr-data->k1*data->k2*data->k3*pow(Tr,data->k2*data->k3-1))/data->Tc;
        d2Alpha=dAlpha*dAlpha/Alpha+Alpha*(-data->k3*(data->k2-1)/Tr/Tr-data->k1*data->k2*data->k3*(data->k2*data->k3-1)*pow(Tr,data->k2*data->k3-2))/data->Tc/data->Tc;
        break;
    case FF_PRTWU95://Peng-robinson Twu with parameters from Tc and w
        if (*T <= data->Tc)
        {
            //double Alpha0=pow(Tr,1.19764*(0.924779-1))*exp(0.272838*(1-pow(Tr,1.19764*0.924779)));
            //double Alpha1=pow(Tr,2.46022*(0.792014-1))*exp(0.625701*(1-pow(Tr,2.46022*0.792014)));
            double Alpha0=pow(Tr,-0.0900876784)*exp(0.272838*(1-pow(Tr,1.1075523216)));
            double Alpha1=pow(Tr,-0.5116913169)*exp(0.625701*(1-pow(Tr,1.9485286831)));
            double dAlpha0=Alpha0*(-0.0900876784/Tr-0.3021823603*pow(Tr,0.1075523216))/data->Tc;
            double dAlpha1=Alpha1*(-0.5116913169/Tr-1.2191963455*pow(Tr,0.9485286831))/data->Tc;
            double d2Alpha0=dAlpha0*dAlpha0/Alpha0+Alpha0*(0.0900876784/Tr/Tr-0.0325004144*pow(Tr,-0.8924476784))/data->Tc/data->Tc;
            double d2Alpha1=dAlpha1*dAlpha1/Alpha1+Alpha1*(0.5116913169/Tr/Tr-1.1564427041*pow(Tr,-0.8924476784))/data->Tc/data->Tc;
            Alpha=Alpha0+data->w*(Alpha1-Alpha0);
            dAlpha=dAlpha0+data->w*(dAlpha1-dAlpha0);
            d2Alpha=d2Alpha0+data->w*(d2Alpha1-d2Alpha0);
        }
        else
        {
            double Alpha0=pow(Tr,-0.2*(4.73020-1))*exp(0.373949*(1-pow(Tr,-0.2*4.7302)));
            double Alpha1=pow(Tr,-8.0*(1.24615-1))*exp(0.0239035*(1-pow(Tr,-8.0*1.24615)));
            double dAlpha0=Alpha0*(-0.2*(4.73020-1)/Tr-0.373949*-0.2*4.7302*pow(Tr,-0.2*4.7302-1))/data->Tc;
            double dAlpha1=Alpha1*(-8.0*(1.24615-1)/Tr-0.0239035*-8.0*1.24615*pow(Tr,-8.0*1.24615-1))/data->Tc;
            double d2Alpha0=dAlpha0*dAlpha0/Alpha0+Alpha0*(0.2*(4.73020-1)/Tr/Tr-0.373949*-0.2*4.7302*(-0.2*4.7302-1)*pow(Tr,-0.2*4.7302-2))/data->Tc/data->Tc;
            double d2Alpha1=dAlpha1*dAlpha1/Alpha1+Alpha1*(8.0*(1.24615-1)/Tr/Tr-0.0239035*-8.0*1.24615*(-8.0*1.24615-1)*pow(Tr,-8.0*1.24615-2))/data->Tc/data->Tc;
            Alpha=Alpha0+data->w*(Alpha1-Alpha0);
            dAlpha=dAlpha0+data->w*(dAlpha1-dAlpha0);
            d2Alpha=d2Alpha0+data->w*(d2Alpha1-d2Alpha0);
        }
        break;
    case FF_PRFIT3://FF_PRSV1 with some parameters no related to critical properties, 3 adjustable parameters instead of a, Tc and k1 (now k3)
        fw = 0.378893 + 1.4897153 * data->w - 0.17131848 * pow(data->w,2) + 0.0196554 * pow(data->w,3);
        Tr = *T / data->k2;
        Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
        dTx=-0.5*pow(Tr,-0.5)/data->k2;//This is the 1sr derivative redarding T
        d2Tx=0.25*pow(Tr,-1.5)/data->k2/data->k2;//and this the 2ndAlpha = pow((1 + data->k2 * (1-pow(*T/data->k3,0.5)),2);
        Alpha = pow((1 + fw * Tx+ data->k3 * (1 - Tr) * (0.7 - Tr)),2);
        dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k3*(2*Tr-1.7)/data->k2);
        d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k3*2/pow(data->k2,2));
        break;
    case FF_PRFIT4://PR with parameters no related to critical properties. There are 4 adjustable parameters:b,a,w and Tc
        fw = 0.37464 + 1.54226 * data->k3 - 0.26992 * data->k3 * data->k3;
        Tr = *T / data->k4;
        Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
        dTx=-0.5*pow(Tr,-0.5)/data->k4;//This is the 1sr derivative redarding T
        d2Tx=0.25*pow(Tr,-1.5)/data->k4/data->k4;//and this the 2ndAlpha = pow((1 + data->k2 * (1-pow(*T/data->k3,0.5)),2);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_SRK://Soaves-Redlich-Kwong.
        fw = 0.48 + 1.574 * data->w - 0.176 * pow(data->w,2);
        //alfa = pow((1 + fw * (1 - pow(Tr,0.5))),2);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PRPOL1://Peng-robinson for polymers b=k1*MW, a=1,Theta=k2*MW (a=1)
        Alpha = data->MW*data->k2;
        dAlpha=0;
        d2Alpha=0;
        break;
    case FF_SRKPOL2://Soave-Redlich-Kwong for polymers b=k1*MW, a=1, Theta=k2*MW*exp(k3*T)
        Alpha = data->MW*data->k2*exp(data->k3* *T);
        dAlpha=Alpha*data->k3;
        d2Alpha=dAlpha*data->k3;
        break;
    default:
        Alpha = 0;
        dAlpha=0;
        d2Alpha=0;
        break;
    }
    param->Theta=param->a*Alpha;
    param->dTheta=param->a*dAlpha;
    param->d2Theta=param->a*d2Alpha;
}

//Calculation of Arr and its derivatives
//--------------------------------------
void CALLCONV FF_ArrDerCubic(double T,double V,const  FF_CubicParam *param,double result[6])
{
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double Veos=V+param->c;
    double T2=T * T;
    double logA=log((Veos+ub)/(Veos+wb));
    result[0]=param->Theta/(param->b*R* T*(param->w-param->u))*logA+log(V/(Veos-param->b));//This is Arr
    result[1]=1/ V- 1/(Veos-param->b)+param->Theta/(R * T * (Veos + ub)*(Veos + wb));//dArr/dV
    result[2]=-1/pow(V,2)+1/pow(Veos-param->b,2)-param->Theta * (Veos * 2+ub+wb)/(R * T *pow(Veos+ub,2)*pow(Veos+wb,2));//d2Arr/dV2
    result[3]=logA*(T * param->dTheta-param->Theta)/(R*param->b*(param->w-param->u)*T2);//dArr/dT
    result[4]=logA* (param->d2Theta*T2-2*(T * param->dTheta-param->Theta))/(R*param->b*(param->w-param->u)*T2* T);//d2Arr/dT2
    result[5]=(T * param->dTheta-param->Theta)/(R*T2*(Veos+ub)*(Veos+wb));//d2Arr/dVdT
}

//Arr and dArr/dT as for cubic EOS
//--------------------------------
void CALLCONV FF_ArrDerCubic0T(double T,double V,const  FF_CubicParam *param,double result[2])
{
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double Veos=V+param->c;
    double T2=T * T;
    double logA=log((Veos+ub)/(Veos+wb));
    result[0]=param->Theta/(param->b*R* T*(param->w-param->u))*logA+log(V/(Veos-param->b));//This is Arr
    result[1]=logA*(T * param->dTheta-param->Theta)/(R*param->b*(param->w-param->u)*T2);//dArr/dT
}

//Arr(reduced residual Helmholtz free energy)and Z calculation for a pure substance, given T and V, according to cubic EOS
//------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVcubic(double T,double V,const  FF_CubicParam *param,double *Arr,double *Z)
{
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double Veos=V+param->c;
    *Arr=param->Theta/(param->b*R* T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(V/(Veos-param->b));//This is Arr
    *Z=V/(Veos-param->b)-param->Theta* V/(R * T * (Veos + ub)*(Veos + wb));//Z
}

//P calculation from T and V using cubic eos
//-------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVcubic(double T,double V,const FF_CubicParam *param,double *P)
{
    double Veos=V+param->c;
    *P=R* T/(Veos-param->b)-param->Theta/((Veos+param->u*param->b)*(Veos+param->w*param->b));
}

//V calculation for a pure substance, given T and P, according to cubic EOS. Arr and Z are also given
//---------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPcubic(double T,double P,const  FF_CubicParam *param,char option,double resultL[3],double resultG[3],char *state)
{

    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double A,B,uw,a2,a1,a0,L,M,N,Nsqr,m,phi1,Root[4],Z[4],Veos;
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    A = param->Theta * P / (R * R * T * T);
    B = param->b * P / (R * T);
    uw=param->u*param->w;
    a2=(param->u+param->w-1)*B-1;
    a1=A+B*(B*(uw-param->u-param->w)-param->u-param->w);
    a0=-B*(A+B*uw*(1+B));
    //printf("a2:%f a1:%f a0:%f\n",a2,a1,a0);
    L = (3 * a1 - a2*a2) / 3;//P
    M = (2 * a2*a2*a2 - 9 * a2 * a1 + 27 * a0) / 27;//Q
    N = M*M / 4 + L*L*L / 27;//R
    //printf("P:%f Q:%f R:%f\n",L,M,N);
    if (N > 0){
        Nsqr=pow(N,0.5);
        Root[0] = (-M / 2 + Nsqr)/fabs((-M / 2 + Nsqr))*pow(fabs((-M / 2 + Nsqr)),0.33333333333) + (-M / 2 - Nsqr)/fabs(-M / 2 - Nsqr)*pow((fabs(-M / 2 - Nsqr)),0.33333333333);//complexity is due to grant a positive base
        Z[0] = Root[0] - a2 / 3;
        //printf("Z:%f\n",Z[0]);
        Veos=Z[0]* R * T / P;
        resultL[0]=resultG[0]=Veos-param->c;//this is V
        resultL[1]=resultG[1]=param->Theta/(param->b*R* T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultL[0]/(Veos-param->b));//This is Arr
        resultL[2]=resultG[2]=Z[0]*resultL[0]/Veos;//Z as per the translated volume
        if (T>=param->Tc) *state='G';
        else if (P>=param->Pc) *state='L';
        else if (Z[0]>=param->Zc) *state='G';
        else *state='L';
        //printf("One solution V:%f state:%c\n",resultL[0],*state);
    }
    else{
        m = 2 * pow(-L / 3,0.5);
        phi1 = acos(3 * M / L / m) / 3;
        Root[1] = m * cos(phi1);
        Root[2] = m * cos(phi1 + 4*M_PI/3);
        Root[3] = m * cos(phi1 + 2*M_PI/3);
        Z[1] = Root[1] - (a2 / 3);
        Z[2] = Root[2] - a2 / 3;
        Z[3] = Root[3] - a2 / 3;
        //printf("Z:%f %f %f\n",Z[1],Z[2],Z[3]);
        if ((Z[1] > Z[2])&& (Z[1] > Z[3])) resultG[2] = Z[1];
        else if (Z[2] > Z[3]) resultG[2] = Z[2];
        else resultG[2] = Z[3];//this is Z gas as P*Veos/(R*T)
        Veos=resultG[2]* R * T / P;
        resultG[0]=Veos-param->c;//this is V
        resultG[1]=param->Theta/(param->b*R* T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultG[0]/(Veos-param->b));//This is Arr
        resultG[2]=resultG[2]*resultG[0]/Veos;//This is Z gas as P*V/(R*T)
        if ((Z[1] < Z[2]) && (Z[1] < Z[3]) && (Z[1]>0)) resultL[2] = Z[1];
        else if ((Z[2] < Z[3]) && (Z[2]>0)) resultL[2] = Z[2];
        else if (Z[3]>0) resultL[2] = Z[3];
        else resultL[2]=Z[1];
        Veos=resultL[2]* R * T / P;
        //printf("Veos: %f c:%f\n",Veos,param->c);
        resultL[0]=Veos-param->c;//this is V
        resultL[1]=param->Theta/(param->b*R* T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultL[0]/(Veos-param->b));//This is Arr
        resultL[2]=resultL[2]*resultL[0]/Veos;
        //printf("V:%f Z:%f\n",resultL[0],resultL[2]);
        *state='b';
        //printf("Two solutions T:%f P:%f Vl:%f Vg:%f Zl:%f Zg:%f\n",T,P,resultL[0],resultG[0],resultL[2],resultG[2]);
    }

}


//SAFT type EOS calculation
//=========================

//Auxiliary calculation for FF_ArrZfromTVSAFT and calcMixPresFF_PCSAFT
//--------------------------------------------------------------------
void CALLCONV FF_calcI1I2(double m,double eta,double I[4])
{
    int i;
    double Aeta[7],Beta[7];
    for (i=0;i<7;i++)
    {
        Aeta[i]=(FF_PCSAFTap[i][0] + (m - 1) / m * FF_PCSAFTap[i][1] + (m - 1) * (m - 2) / pow(m,2) * FF_PCSAFTap[i][2])* pow(eta,i);
        I[0]=I[0]+Aeta[i];
        I[1]=I[1]+Aeta[i]*(i+1);
        Beta[i]=(FF_PCSAFTbp[i][0] + (m - 1) / m * FF_PCSAFTbp[i][1] + (m - 1) * (m - 2) / pow(m,2) * FF_PCSAFTbp[i][2])* pow(eta,i);
        I[2]=I[2]+Beta[i];
        I[3]=I[3]+Beta[i]*(i+1);
    }
}


//Z and Arr calculation for a pure substance, given T and V, according to PCSAFT EOS
//----------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVSAFT(double T,double V,const  FF_SaftEOSdata *data,double *Arr,double *Z)
{
    //static int counter=0;
    //counter++;
    int i;
    double sigma,epsilon,epsilon_kT;
    sigma=data->sigma*1e-10;//in SI units
    epsilon=data->epsilon*kb;//Energy depth in SI units
    epsilon_kT=data->epsilon/ T;
    double d;//Temperature corrected segment diameter
    double Vm,rhoM,rhoS,eta;//Molecular volumen and density, segment density and packing fraction. In SI units
    double Ahs,Zhs;//contribution by hard spheres
    double Amono,Zmono;//total monomer contribution
    double ghs,dLghs_dEta, dLghs_dRhoM,dLghs_dRhoS;//hard spheres radial distribution function, and its derivative regarding eta,rhoM and rhoS
    double Achain,Zchain;//Total chain contribution
    double Aassoc=0,Zassoc=0;//
    double Add=0,Zdd=0;

    double eta2;//Second power of eta
    Vm = V / Av;//molecular volume in m3
    rhoM = 1 / Vm;//number of molecules/m3
    rhoS=data->m*rhoM;//number of segments/m3

    //PCSAFT variations
    //-----------------
    double Ahchain,Zhchain;//contribution by hard chains formation
    double Adisp,Zdisp;//contribution by dispersion between chains
    d = sigma * (1 - 0.12 * exp(-3 * epsilon_kT)); //Hard sphere diameter in m, at given T
    eta = (Pi * pow(d,3) / 6) * rhoS; //Volume fraction filled with hard spheres.
    eta2=eta*eta;

    //Contribution by monomers
    Ahs=(4 * eta - 3 * eta2) / pow((1 - eta),2);
    Zhs=(4 * eta - 2 * eta2) / pow((1 - eta),3);
    Amono = data->m*Ahs;
    Zmono = data->m *Zhs;
    //printf("PCSAFT Amono:%f Zmono:%f\n",Amono,Zmono);
    //contribution by chain
    ghs = (1 - 0.5*eta) / pow((1 - eta),3);//radial distribution function for hard spheres
    dLghs_dEta =(2.5-eta)/((1-0.5*eta)*(1-eta));
    dLghs_dRhoM = dLghs_dEta*eta*Vm;
    Ahchain = -(data->m-1) * log(ghs);
    Zhchain = -(data->m-1) * eta*dLghs_dEta;
    //printf("ghs:%f dLghs_dEta:%f\n",ghs,dLghs_dEta);
    //printf("PCSAFT Ahchain:%f Zhchain:%f\n",Ahchain,Zhchain);
    //contribution by dispersion (attraction between chains)
    double Z1,C1,C2,Z2,I[4]={0.0,0.0,0.0,0.0};;
    FF_calcI1I2(data->m,eta,I);
    Z1 = -2 * Pi / Vm * I[1] * pow(data->m,2) * epsilon_kT * pow(sigma,3);
    C1 = 1/(1 + data->m * (8 * eta - 2 * eta2) / pow((1 - eta),4) + (1 - data->m) * (20 * eta - 27 * eta2
            + 12 * pow(eta,3) - 2 * pow(eta,4)) / pow(((1 - eta) * (2 - eta)),2));
    C2 = C1 * (data->m * (-4 * eta2 + 20 * eta + 8) / pow((1 - eta),5) + (1 - data->m) * (2 * pow(eta,3)
            + 12 * eta2 - 48 * eta + 40) / pow(((1 - eta) * (2 - eta)),3));
    Z2 = -Pi / Vm * data->m * C1 * (I[3] - C2 * eta * I[2])* pow(data->m,2) * pow(epsilon_kT,2) * pow(sigma,3);
    Adisp = -2 * Pi / Vm * I[0] * pow(data->m,2) * data->epsilon * pow(sigma,3) / T - Pi / Vm * data->m * C1
            * I[2] * pow(data->m,2) * pow(epsilon_kT,2) * pow(sigma,3);
    Zdisp = Z1 + Z2;
    //printf("PCSAFT Adisp:%f Zdisp:%f\n",Adisp,Zdisp);
    Achain=Ahchain+Adisp;
    Zchain=Zhchain+Zdisp;

    //Contribution by molecular association

    if ((data->kAB > 0) && (data->epsilonAB > 0)) //If the molecule has association parameters
    {
        double DeltaAB,X[data->nPos+data->nNeg+data->nAcid]; //X=[] is fraction of molecules not associated at site i
        double sum;
        //DeltaAB = pow(d,3) * ghs * data->kAB * (exp(data->epsilonAB / T) - 1);
        DeltaAB = pow(sigma,3) * ghs * data->kAB * (exp(data->epsilonAB / T) - 1);
        //Calculation taking account of number of association sites of the molecule(1=acids,2=alcohol,4=water or diols)
        if (data->nAcid==1){//1A
            X[0]=(-1 + pow((1 + 4 * rhoM * DeltaAB),0.5)) / (2 * rhoM * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if (data->nPos==1 && data->nNeg==1){//2B
            X[0]=X[1]=(-1 + pow((1 + 4 * rhoM * DeltaAB),0.5)) / (2 * rhoM * DeltaAB);
        }
        else if (data->nPos==2 && data->nNeg==2)//4C
        {
            X[0]=X[1]=X[2]=X[3]=(-1 + pow((1 + 8 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
        }
        else if ((data->nPos==2 && data->nNeg==1)||(data->nPos==1 && data->nNeg==2)){//3B
            X[0]=X[1]=(-(1 - rhoM * DeltaAB) + pow((pow((1 + rhoM * DeltaAB),2) + 4 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
            X[2]=(2 * X[0] - 1);
        }
        else if (data->nAcid==2){//2A
            X[0]=X[1]=(-1 + pow((1 + 8 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if ((data->nPos==3 && data->nNeg==1)||(data->nPos==1 && data->nNeg==3)){//4B
            X[0]=X[1]=X[2]=(-(1 - 2*rhoM * DeltaAB) + pow((pow((1 + 2*rhoM * DeltaAB),2) + 4 * rhoM * DeltaAB),0.5)) / (6 * rhoM * DeltaAB);
            X[3]=(3 * X[0] - 2);
        }
        else if (data->nAcid==3){//3A
            X[0]=X[1]=X[2]=(-1 + pow((1 + 12 * rhoM * DeltaAB),0.5)) / (6 * rhoM * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if (data->nAcid==4){//4A
            X[0]=X[1]=X[3]=X[4]=(-1 + pow((1 + 16 * rhoM * DeltaAB),0.5)) / (8 * rhoM * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else{
            for (i=0;i<data->nPos+data->nNeg+data->nAcid;i++) X[i]=1;
        }
        sum=0;
        for (i=0;i<(data->nPos+data->nNeg+data->nAcid);i++){
            //printf("i:%i Xi:%f\n",i,X[i]);
            sum = sum + 1 -X[i];
        }
        //printf("sum:%f dLghs_drhoS:%f\n",sum,dLghs_drhoS);
        Aassoc = (data->nPos+data->nNeg+data->nAcid)/ 2;
        for (i=0; i<(data->nPos+data->nNeg+data->nAcid);i++)
            Aassoc = Aassoc + (log(X[i]) - X[i] / 2);
        Zassoc=-0.5*(1+rhoM*dLghs_dRhoM)*sum;
        //printf("Aassoc:%f Zassoc:%f\n",Aassoc,Zassoc);
    }

    //contribution by polar forces
    //1Debbie=3.33564 e-30 C.m(SI units)
    //U=mu*mu/(4*pi*epsilon0*r^3)
    if(data->mu>0){
        if (data->xp>=1){  //Gross and Vrabeck model.
            double Add2,dAdd2_dRhoM,dAdd_dRhoM;
            double mEfec,muRed2;
            double ad[3][5]={{0.3043504,-0.1358588,1.4493329,0.3556977,-2.0653308},{0.9534641,-1.8396383,2.0131180,-7.3724958,8.2374135},{-1.1610080,4.5258607,0.9751222,-12.281038,5.9397575}};
            double bd[3][5]={{0.2187939,-1.1896431,1.1626889,0.0,0.0},{-0.5873164,1.2489132,-0.5085280,0.0,0.0},{3.4869576,-14.915974,15.372022,0.0,0.0}};
            double cd[3][5]={{-0.0646774,0.1975882,-0.8087562,0.6902849,0.0},{-0.9520876,2.9924258,-2.3802636,-0.2701261,0.0},{-0.6260979,1.2924686,1.6542783,-3.4396744,0.0}};

            double a[5],b[5],c[5];
            double Jdd2=0,Jdd3=0,dJdd2_dRhoM=0,dJdd3_dRhoM=0;
            double sigma3,d3,mEfecAux,aux1,aux2;
            sigma3=sigma*sigma*sigma;
            d3=d*d*d;
            //eta = (Pi * sigma3 / 6) * rhoS;
            if (data->m>2) mEfec=2;
            else mEfec=data->m;
            mEfecAux=(mEfec-1)*(mEfec-2)/(mEfec*mEfec);
            muRed2=pow(data->mu,2)/(data->m*sigma3*data->epsilon)*7.24311E-27;//muRed2=pow((data->mu*3.33564e-30),2)/(4*Pi*8.854e-12*m*d^3*epsilon);//epsilon=kb*data->epsilon
            //printf("muRed:%f\n",pow(muRed2,0.5));
            //printf("sigma3:%f epsilonT:%f muRed2:%f\n",sigma3*1e30,data->epsilon/ T,muRed2);
            //printf("eta:%f\n",eta);
            int n;
            for (n=0;n<5;n++)
            {
                a[n]=ad[0][n]+ad[1][n]*(mEfec-1)/mEfec+ad[2][n]*mEfecAux;
                b[n]=bd[0][n]+bd[1][n]*(mEfec-1)/mEfec+bd[2][n]*mEfecAux;
                c[n]=cd[0][n]+cd[1][n]*(mEfec-1)/mEfec+cd[2][n]*mEfecAux;
                Jdd2=Jdd2+(a[n]+b[n]*epsilon_kT)*pow(eta,n);
                dJdd2_dRhoM=dJdd2_dRhoM+(a[n]+b[n]*epsilon_kT)*n*pow(eta,n-1);//these are the derivatives regarding eta. Later we change to rhoM
                Jdd3=Jdd3+c[n]*pow(eta,n);
                dJdd3_dRhoM=dJdd3_dRhoM+c[n]*n*pow(eta,n-1);
            }
            //printf("Jdd2:%f dJdd2_Eta:%f\n",Jdd2,dJdd2_dRhoM);
            //printf("Jdd3:%f dJdd3_Eta:%f\n",Jdd3,dJdd3_dRhoM);
            dJdd2_dRhoM=dJdd2_dRhoM*Pi*data->m*d3/6;
            dJdd3_dRhoM=dJdd3_dRhoM*Pi*data->m*d3/6;
            aux1=-Pi*pow(epsilon_kT*data->xp*muRed2,2)*sigma3;
            Add2=aux1*rhoM*Jdd2;
            dAdd2_dRhoM=aux1*(Jdd2+rhoM*dJdd2_dRhoM);
            aux2=4*Pi*data->xp*muRed2*epsilon_kT*sigma3/3;
            Add=Add2/(1-aux2*rhoM*Jdd3/Jdd2);
            dAdd_dRhoM=(dAdd2_dRhoM*(1-aux2*rhoM*Jdd3/Jdd2)+Add2*aux2*(((Jdd3+rhoM*dJdd3_dRhoM)*Jdd2-rhoM*Jdd3*dJdd2_dRhoM)/(Jdd2*Jdd2)))/((1-aux2*rhoM*Jdd3/Jdd2)*(1-aux2*rhoM*Jdd3/Jdd2));
            //printf("Add2:%f Add3:%f\n",Add2,aux2*aux1*rhoM*rhoM*Jdd3);
            Zdd=rhoM*dAdd_dRhoM;
            //printf("AddGV:%f ZddGV:%f\n",Add,Zdd);
        }

        else if ((data->xp>0)&&(data->xp<1)){	//Jog and Chapman model
            //In the original work of Jog and Chapman in Add2 does not appear data->m in the multiplication,
            //But yes in literature. And results are better using Al Saifi coefficients, at least for methanol
            double Add2;
            double mu2,muRed2,rhoRed,I2,I3;
            double rhoRed2,rhoRed3,aux1,aux2;
            double dI2_dRhoRed,dI3_dRhoRed,dAdd2_dRhoRed,dAdd_dRhoRed;
            if (data->mu>0 && data->mu<10 && data->xp>0 && data->xp<1)
            {
                mu2=data->mu*data->mu*1.000021e-49;//mu2=pow((data->mu*3.33564e-30),2)/(4*Pi*8.854e-12);
                muRed2=mu2/(kb* T*d*d*d);
                //printf("MuRed:%f\n",pow(muRed2,2));
                rhoRed=rhoS*pow(d,3);
                rhoRed2=rhoRed*rhoRed;
                rhoRed3=rhoRed2*rhoRed;
                aux1=1-0.3618*rhoRed-0.3205*rhoRed2+0.1078*rhoRed3;
                aux2=1-0.5236*rhoRed;
                I2=aux1/(aux2*aux2);
                dI2_dRhoRed=((-0.3618-2*0.3205*rhoRed+3*0.1078*rhoRed2)*aux2+2*0.5236*aux1)/(aux2*aux2*aux2);
                aux1=(1+0.62378*rhoRed-0.11658*rhoRed2);
                aux2=(1-0.59056*rhoRed+0.20059*rhoRed2);
                I3=aux1/aux2;
                dI3_dRhoRed=((0.62378-2*0.11658*rhoRed)*aux2-aux1*(-0.59056+2*0.20059*rhoRed))/(aux2*aux2);
                //printf("I2:%f I3:%f\n",I2,I3);
                aux1=-2*Pi*data->m*pow(data->xp*muRed2,2)/9;
                Add2=aux1*rhoRed*I2;
                dAdd2_dRhoRed=aux1*(I2+rhoRed*dI2_dRhoRed);
                aux1=5*Pi*data->xp*muRed2/36;
                aux2=1+aux1*rhoRed*I3/I2;
                Add=Add2/aux2;
                dAdd_dRhoRed=(dAdd2_dRhoRed*aux2-Add2*(aux1*((I3+rhoRed*dI3_dRhoRed)*I2-rhoRed*I3*dI2_dRhoRed)/(I2*I2)))/(aux2*aux2);
                //printf("Add2:%f\n",Add2);
                Zdd=rhoRed*dAdd_dRhoRed;
                //printf("AddJC:%f ZddJC:%f\n",Add,Zdd);
            }
        }
    }
    if(data->Q>0){
        double Qred2;
        Qred2=7243*data->Q*data->Q/(data->epsilon*pow(data->sigma,5));
        //Aqq2=-14*Pi*Av*rhoM*Qred2*Qred2*J10/(5*kb* T*sigma3*sigma3*sigma);

    }


    *Arr= Amono + Achain + Aassoc+ Add;//Reduced residual Helmholz free energy
    *Z = 1 + Zmono + Zchain + Zassoc + Zdd;//Z
    //printf("Arr:%f Z:%f\n",*Arr,*Z);
}



//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to PCSAFT or SAFT-VRMie EOS
//----------------------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrDerSAFT(double T,double V,const  FF_SaftEOSdata *data,double result[6])
{
    double dV=V * 0.00001;//increments of V and T used to obtain dArr/dV,dArr/dT and d2Arr/dT2 in SAFT eos
    double Vplus=V + dV;
    double dT=0.01;
    double Tplus=T+dT;
    double Tminus=T-dT;
    double Arr,Z,ArrVplus,ZVplus,ArrTplus,ZTplus,ArrTminus,ZTminus;
    FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
    FF_ArrZfromTVSAFT(T,Vplus,data,&ArrVplus,&ZVplus);
    FF_ArrZfromTVSAFT(Tplus,V,data,&ArrTplus,&ZTplus);
    FF_ArrZfromTVSAFT(Tminus,V,data,&ArrTminus,&ZTminus);
    result[0]=Arr;//This is Arr
    result[1]=(1- Z)/ V;//dArr/dV at constant T
    result[2]=((1- ZVplus)/ Vplus-result[1])/dV;//d2Arr/dV2 at constant T
    result[3]=(ArrTplus-Arr)/dT;//dArr/T at constant V
    result[4]=(result[3]-(Arr-ArrTminus)/dT)/dT;//d2Arr/dT2 at constant V
    result[5]=((1- ZTplus)/ V-result[1])/dT;//d2Arr/dVdT
}

//Arr and its partial derivatives symbolic calculation for a pure substance, given T and V, according to PCSAFT EOS. To be implemented
//-----------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrDerPCSAFT(double T,double V,const  FF_SaftEOSdata *data,double result[6])
{
    int i;
    double sigma,epsilon,epsilon_kT;
    sigma=data->sigma*1e-10;//in SI units
    epsilon=data->epsilon*kb;//Energy depth in SI units
    epsilon_kT=data->epsilon/ T;
    double d;//Temperature corrected segment diameter
    double Vm,rhoM,rhoS,eta;//Molecular volumen and density, segment density and packing fraction. In SI units
    double Amono,Zmono;//total monomer contribution
    double ghs,dLghs_dEta, dLghs_dRhoM;//hard spheres radial distribution function, and its derivative regarding eta,rhoM and rhoS
    double Achain,Zchain;//Total chain contribution
    double Aassoc=0,Zassoc=0;//
    double Add=0,Zdd=0;

    double eta2;//Second power of eta
    Vm = V / Av;//molecular volume in m3
    rhoM = 1 / Vm;//number of molecules/m3
    rhoS=data->m*rhoM;//number of segments/m3
    double Ahchain,Zhchain;//contribution by hard chains formation
    double Adisp,Zdisp;//contribution by dispersion between chains
    d = sigma * (1 - 0.12 * exp(-3 * epsilon_kT)); //Hard sphere diameter in m, at given T
    eta = (Pi * pow(d,3) / 6) * rhoS; //Volume fraction filled with hard spheres.
    eta2=eta*eta;

    //Contribution by monomers
    Amono = data->m*(4 * eta - 3 * eta2) / pow((1 - eta),2);
    Zmono = data->m *(4 * eta - 2 * eta2) / pow((1 - eta),3);

    //contribution by chain
    ghs = (1 - 0.5*eta) / pow((1 - eta),3);//radial distribution function for hard spheres
    dLghs_dEta =(2.5-eta)/((1-0.5*eta)*(1-eta));
    dLghs_dRhoM = dLghs_dEta*eta*Vm;
    Ahchain = -(data->m-1) * log(ghs);
    Zhchain = -(data->m-1) * eta*dLghs_dEta;

    //contribution by dispersion (attraction between chains)
    double Z1,C1,C2,Z2,I[4]={0.0,0.0,0.0,0.0};;
    FF_calcI1I2(data->m,eta,I);
    Z1 = -2 * Pi / Vm * I[1] * pow(data->m,2) * epsilon_kT * pow(sigma,3);
    C1 = 1/(1 + data->m * (8 * eta - 2 * eta2) / pow((1 - eta),4) + (1 - data->m) * (20 * eta - 27 * eta2
            + 12 * pow(eta,3) - 2 * pow(eta,4)) / pow(((1 - eta) * (2 - eta)),2));
    C2 = C1 * (data->m * (-4 * eta2 + 20 * eta + 8) / pow((1 - eta),5) + (1 - data->m) * (2 * pow(eta,3)
            + 12 * eta2 - 48 * eta + 40) / pow(((1 - eta) * (2 - eta)),3));
    Z2 = -Pi / Vm * data->m * C1 * (I[3] - C2 * eta * I[2])* pow(data->m,2) * pow(epsilon_kT,2) * pow(sigma,3);
    Adisp = -2 * Pi / Vm * I[0] * pow(data->m,2) * data->epsilon * pow(sigma,3) / T - Pi / Vm * data->m * C1
            * I[2] * pow(data->m,2) * pow(epsilon_kT,2) * pow(sigma,3);
    Zdisp = Z1 + Z2;

    Achain=Ahchain+Adisp;
    Zchain=Zhchain+Zdisp;
    if ((data->kAB > 0) && (data->epsilonAB > 0)) //If the molecule has association parameters
    {
        double DeltaAB,X[data->nPos+data->nNeg+data->nAcid]; //X=[] is fraction of molecules not associated at site i
        double sum;
        //DeltaAB = pow(d,3) * ghs * data->kAB * (exp(data->epsilonAB / T) - 1);
        DeltaAB = pow(sigma,3) * ghs * data->kAB * (exp(data->epsilonAB / T) - 1);
        //Calculation taking account of number of association sites of the molecule(1=acids,2=alcohol,4=water or diols)
        if (data->nAcid==1) X[0]=(-1 + pow((1 + 4 * rhoM * DeltaAB),0.5)) / (2 * rhoM * DeltaAB);
        else if (data->nPos==1 && data->nNeg==1) X[0]=X[1]=(-1 + pow((1 + 4 * rhoM * DeltaAB),0.5)) / (2 * rhoM * DeltaAB);
        else if (data->nPos==2 && data->nNeg==2) X[0]=X[1]=X[2]=X[3]=(-1 + pow((1 + 8 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
        else if ((data->nPos==2 && data->nNeg==1)||(data->nPos==1 && data->nNeg==2)){//3B
            X[0]=X[1]=(-(1 - rhoM * DeltaAB) + pow((pow((1 + rhoM * DeltaAB),2) + 4 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
            X[2]=(2 * X[0] - 1);
        }
        else if (data->nAcid==2) X[0]=X[1]=(-1 + pow((1 + 8 * rhoM * DeltaAB),0.5)) / (4 * rhoM * DeltaAB);
        else if ((data->nPos==3 && data->nNeg==1)||(data->nPos==1 && data->nNeg==3)){//4B
            X[0]=X[1]=X[2]=(-(1 - 2*rhoM * DeltaAB) + pow((pow((1 + 2*rhoM * DeltaAB),2) + 4 * rhoM * DeltaAB),0.5)) / (6 * rhoM * DeltaAB);
            X[3]=(3 * X[0] - 2);
        }
        else if (data->nAcid==3) X[0]=X[1]=X[2]=(-1 + pow((1 + 12 * rhoM * DeltaAB),0.5)) / (6 * rhoM * DeltaAB);
        else if (data->nAcid==4) X[0]=X[1]=X[3]=X[4]=(-1 + pow((1 + 16 * rhoM * DeltaAB),0.5)) / (8 * rhoM * DeltaAB);
        else for (i=0;i<data->nPos+data->nNeg+data->nAcid;i++) X[i]=1;
        sum=0;
        for (i=0;i<(data->nPos+data->nNeg+data->nAcid);i++) sum = sum + 1 -X[i];
        Aassoc = (data->nPos+data->nNeg+data->nAcid)/ 2;
        for (i=0; i<(data->nPos+data->nNeg+data->nAcid);i++)  Aassoc = Aassoc + (log(X[i]) - X[i] / 2);
        Zassoc=-0.5*(1+rhoM*dLghs_dRhoM)*sum;
    }
/*
    result[0]=Arr;//This is Arr
    result[1]=(1- Z)/ V;//dArr/dV at constant T
    result[2]=((1- ZVplus)/ Vplus-result[1])/dV;//d2Arr/dV2 at constant T
    result[3]=(ArrTplus-Arr)/dT;//dArr/T at constant V
    result[4]=(result[3]-(Arr-ArrTminus)/dT)/dT;//d2Arr/dT2 at constant V
    result[5]=((1- ZTplus)/ V-result[1])/dT;//d2Arr/dVdT*/
}

//Arr and dArr/dT as per SAFT EOS
void CALLCONV FF_ArrDerSAFT0T(double T,double V,const  FF_SaftEOSdata *data,double result[2])
{
    double dT=0.01;
    double Tplus=T+dT;
    double Arr,Z,ArrTplus,ZTplus;
    FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
    FF_ArrZfromTVSAFT(Tplus,V,data,&ArrTplus,&ZTplus);
    result[0]=Arr;//This is Arr
    result[1]=(ArrTplus-Arr)/dT;//dArr/T at constant V
}


//P calculation from T and V using according to FF_PCSAFT EOS
//--------------------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVSAFT(double T,double V,const  FF_SaftEOSdata *data,double *P)
{
    double Arr,Z;
    FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
    *P=Z*R* T/ V;
}


//Volume solver from T and P using a SAFT or SW EOS. Regula Falsi solver, Anderson-Bjork modification.
double CALLCONV FF_VsolverRegula(void (*f)(double, double, void *, double *), void *data, double T, double P, double a, double b){//a and b define the interval
    double ftol=1.0e-5;
    int niter=15;
    double fx,fa,fb;
    double x=+HUGE_VALF;//The solution
    int n=0;
    int side=0;
    double ex,ea,eb;//errors at x, a and b points
    (*f)(T,a,data,&fa);
    ea=fa-P;
    (*f)(T,b,data,&fb);
    eb=fb-P;
    if (ea*eb>0) return x;//if both functions have the same sign, it is not sure that a root exists
    if ((ea/P<ftol)&&(ea/P>(-ftol))) x=a;//if one of the given value was already the solution
    else if ((eb/P<ftol)&&(eb/P>(-ftol))) x=b;
    else{
        while (n<niter){
            x=a+(b-a)*ea/(ea-eb);
            (*f)(T,x,data,&fx);
            n++;
            ex=fx-P;
            //printf("n:%i P:%f a:%f ea:%f b:%f eb:%f x:%e fx:%e ex:%f\n",n+1,P,a,ea,b,eb,x,fx,ex);
            if (((ex/ea)>0.8)||((ex/eb)>0.8)){
                x=(a+b)/2;
                (*f)(T,x,data,&fx);
                ex=fx-P;
                if ((eb*ex)>0){
                    b=x;
                    eb=ex;
                }
                else{
                    a=x;
                    ea=ex;
                }
                n++;
                //if (T>350) printf("New n:%i P:%f a:%f ea:%f b:%f eb:%f x:%e fx:%e ex:%f\n",n,P,a,ea,b,eb,x,fx,ex);
            }
            if ((ex/P<ftol)&&(ex/P>(-ftol))) break;
            if ((eb*ex)>0){//x and b are at the same side
                if (side==1){//if it happened also in the previous loop
                    if ((ex/eb)<1.0) ea=ea-ea*(ex/eb);//we decrease P-f(a) for the next loop calculation
                    else ea*=0.5;
                }
                side = 1;
                b=x;
                eb=ex;
            }
            else{
                if (side==-1){
                    if ((ex/ea)<1.0) eb=eb-eb*(ex/ea);
                    else eb *= 0.5;
                }
                side = -1;
                a=x;
                ea=ex;
            }
            //n=n+1;
            //if (T>350) printf("Iteracion n:%i \n",n);
        }
        if (n>=niter) x=HUGE_VALF;
    }
    return x;
}


//Z,Arr,V calculation for a pure substance, given T and P, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------
void CALLCONV FF_VfromTPsaft(double T,double P,const  FF_SaftEOSdata *data,char option,double resultL[3],double resultG[3],char *state)
{
    int ver=0;
    double V,Arr,Z,Pcalc,Vplus,ArrPlus,Zplus,error,dP_dV;
    int i;
    double maxError=0.0001;//maximum relative error accepted
    //double maxError=0.00005;//maximum relative error accepted
    double fw;
    //if (data->Vc==0) data->Vc=data->Zc*R*data->Tc/data->Pc;
    FF_CubicParam param;
    param.b = 0.077796 * R * data->Tc / data->Pc;
    param.a = 5.877359*param.b*R*data->Tc;//0.457235 * pow(R*data->Tc,2)/ data->Pc;
    param.u=2.414214;//1+2^0.5
    param.w=-0.414214;//1-2^0.5
    param.c=0.0;
    fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
    param.Theta=param.a*pow((1 + fw * (1-pow((T / data->Tc),0.5))),2);
    char stateCubic;
    //printf("T:%f P:%f eos:%i Vc:%f\n",T,P,data->eos,data->Vc);
    //printf("b:%f u:%f w:%f a:%f Theta:%f\n",param.b, param.u,param.w,param.a,param.Theta);
    FF_VfromTPcubic(T,P,&param,option,resultL,resultG,&stateCubic);
    if (ver==1) printf("V finding T:%f P to obtain:%f Cubic: Vl:%f Vg:%f\n",T,P,resultL[0],resultG[0]);
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    //if ((option!='g')||(!(stateCubic=='b'))){//if the cubic has only one solution it does not mean that SAFT also
    if (option!='g'){//we calculate the liquid phase if not only the gas phase has been asked, or there is only one phase
        if (data->Vc>0) V=fmin(resultL[0]*0.8,data->Vc*0.5);
        else V=resultL[0]*0.8;
        Vplus=V*1.00001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
        FF_ArrZfromTVSAFT(T,Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        Pcalc=Z*R* T/V;
        error=P-Pcalc;
        i=1;
        if (ver==1) printf("i= %d  Vl=%f Zl:%f error= %f  P= %f  dP_dV= %f\n",i,V,Z,error,Pcalc,dP_dV);
        while ((fabs(error/ P)>maxError) && (dP_dV < 0)&&(i<15))
        {
            if (T<0.97*data->Tc){
                if (error>0) V=V+0.5*error/dP_dV;//we need to decrease volume
                else V=V+error/dP_dV;//increase volume, should be the normal case
            }

            else if(T<1.1*data->Tc)
            {
                if (error>0) V=V+fmax(error/dP_dV,-0.1*V);//decrease volume
                else V=V+fmin(error/dP_dV,0.1*V);
            }
            else//over Tc
            {
                if (error>0) V=V+0.03*error/dP_dV;//decrease volume
                else V=V+0.7*error/dP_dV;
            }
            Vplus=V*1.00001;
            FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
            Pcalc=Z*R* T/V;
            error=P-Pcalc;
            FF_ArrZfromTVSAFT(T,Vplus,data,&ArrPlus,&Zplus);
            dP_dV=(R* T*Zplus/Vplus-Pcalc)/(Vplus-V);
            i=i+1;
            if (ver==1) printf("i= %d  Vl=%f Zl:%f error= %f  P= %f  dP_dV= %f\n",i,V,Z,error,Pcalc,dP_dV);
        }
        if ((dP_dV < 0)&&(fabs(error/ P)<=maxError)){
            resultL[0]=V;
            resultL[1]=Arr;
            resultL[2]=Z;
            *state='l';
           if (ver==1) printf("Final liquid Vl:%f Zl:%f Arrl:%f\n",resultL[0],resultL[2],resultL[1]);
        }
        else
        {
            if (ver==1) printf("fallo liquido\n");
            resultL[0]=resultL[1]=resultL[2]=0;
        }
    }

    if (option!='l')//and the gas phase if not only the liquid one has been asked, or single value and not found as liquid
    {   V=fmax(resultG[0]*1.4,data->Vc*3);
        //V=resultG[0];
        //Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        Vplus=V*0.999999; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
        FF_ArrZfromTVSAFT(T,Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        Pcalc=Z*R* T/V;
        error=P-Pcalc;
        i=1;
        if (ver==1) printf("Inicio gas state:%c Vl:%f i= %d  Vg:%f Pcalc:%f error= %f  dP_dV= %f\n",*state,resultL[0],i,V, Pcalc,error,dP_dV);
        while ((fabs(error/ P)>maxError) && (dP_dV < 0)&&(i<15)&&(!((*state=='l')&&(V<=resultL[0]))))
        {   if (T<0.8*data->Tc){
                if (error>0) V=V+0.7*error/dP_dV;//need to reduce volume
                else V=V+0.85*error/dP_dV;
            }
            else if(T<0.95*data->Tc)
            {
                if (error>0) V=V+fmax(error/dP_dV,-0.4*V);//reduce volume
                else V=V+fmin(error/dP_dV,0.4*V);
            }
            else if(T<=1.1*data->Tc)
            {
                if (error>0) V=V+fmax(error/dP_dV,-0.1*V);//reduce volume
                else V=V+fmin(error/dP_dV,0.1*V);
            }
            else V=V+error/dP_dV;
            //Vplus=V*1.000001;
            Vplus=V*0.999999;
            FF_ArrZfromTVSAFT(T,V,data,&Arr,&Z);
            Pcalc=Z*R* T/V;
            error=P-Pcalc;
            FF_ArrZfromTVSAFT(T,Vplus,data,&ArrPlus,&Zplus);
            dP_dV=(R* T*Zplus/Vplus-Pcalc)/(Vplus-V);
            i=i+1;
            if (ver==1) printf("i= %d  Vg= %f Pcalc=%f error= %f dP_dV= %f\n",i,V,Pcalc,error,dP_dV);
        }

        if ((dP_dV < 0)&&(fabs(error/ P)<=maxError))
        //if (fabs(error/ P)<=maxError)
        {
            resultG[0]=V; //This is V
            resultG[1]=Arr;
            resultG[2]=Z;
            if (*state=='l'){
                if (fabs((resultL[0]-resultG[0])/resultL[0])>0.01) *state='b';
                else *state='u';
            }
            else *state='g';
        }
        else
        {
            resultG[0]=resultG[1]=resultG[2]=0;
        }
    }
    if (((*state=='u')||(*state=='l'))&&((T>=data->Tc)||(resultL[0]>=R*data->Zc*data->Tc/data->Pc))){
            //printf("resultL[0]:%f\n",resultL[0]);
            resultG[0]=resultL[0];
            resultG[1]=resultL[1];
            resultG[2]=resultL[2];
            *state='G';
    }
    if (*state=='l') *state='L';
    else if  (*state=='g') *state='G';
    if (ver==1) printf("Final state:%c Vl:%f Vg:%f Zl:%f Zg:%f\n",*state,resultL[0],resultG[0],resultL[2],resultG[2]);
}



//SW EOS calculation
//==================

//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given tau and delta, according to SW EOS
//----------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrDerSW(double tau,double delta,const  FF_SWEOSdata *data,double result[6])
{   //It seems that there is a problem with the calculation of partial derivatives, at least regarding delta in the additional terms
    //The result of all derivatives are OK according to the IAPWS95 test, and are OK normally. But close to the critical point they differ from
    //numerical derivatives, and are really incorrect. The correct ones are the numerical. For IAPWS95 and CO2 perhaps is better to use numerical close to CP
    int i,j;
    double tau2=tau * tau;
    double delta2=delta* delta;
    double dt=delta * tau;
    double partial;
    double d12,D,Db,z,p;
    double dp_d, dz_d,dD_d,dDb_d;
    double d2z_d2,d2p_d2,d2D_d2,d2Db_d2;
    double dp_t,dz_t,dD_t,dDb_t,d2p_t2,d2Db_t2;
    double d2p_dt,d2D_dt,d2Db_dt;
    for (i=0;i<6;i++) result[i]=0;
    for (i=0;i<data->nPol;i++)//Polinomial terms
    {
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);//n*delta^d*tau^t
        result[0]=result[0]+partial;//This will be Arr
        result[1]=result[1]+data->d[i]*partial/ delta;//This will be dArr/ddelta at constant tau
        result[2]=result[2]+data->d[i]*(data->d[i]-1)*partial/ delta2;//This will be d2Arr/ddelta2 at constant tau
        result[3]=result[3]+data->t[i]*partial/ tau;//This will be dArr/dtau at constant delta
        result[4]=result[4]+data->t[i]*(data->t[i]-1)*partial/ tau2;//This will be d2Arr/dtau2 at constant delta
        result[5]=result[5]+data->t[i]*data->d[i]*partial/ dt;//This will be d2Arr/ddelta_dteta
        //result[6]=result[6]+data->d[i]*(data->d[i]-1)*(data->d[i]-2)*partial/ delta / delta/ delta;//This will be d3Arr/ddelta3 at constant tau
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)//exponential terms
    {
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]));//n*delta^d*tau^t*e^(-delta^c)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->d[i]-data->c[i]*pow(delta,data->c[i]))*partial/ delta;
        result[2]=result[2]+((data->d[i]-data->c[i]*pow(delta,data->c[i]))*(data->d[i]-1-data->c[i]*pow(delta,data->c[i]))-data->c[i]*data->c[i]*pow(delta,data->c[i]))
                *partial/ delta2;
        result[3]=result[3]+data->t[i]*partial/ tau;
        result[4]=result[4]+data->t[i]*(data->t[i]-1)*partial/ tau2;
        result[5]=result[5]+data->t[i]*(data->d[i]-data->c[i]*pow(delta,data->c[i]))*partial/ dt;
        //result[6]=result[6]+
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)//special exponential terms
    {
        j=i-data->nPol-data->nExp;
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        //n*delta^d*tau^t*e^(-a*(delta-e)^2-b*(tau-g)^2)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*partial/ delta;
        result[2]=result[2]+(pow(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]),2)-data->d[i]-2* delta* delta*data->a[j])*partial/ delta2;
        result[3]=result[3]+(data->t[i]-2* tau*data->b[j]*(tau-data->g[j]))*partial/ tau;
        result[4]=result[4]+(pow(data->t[i]-2* tau*data->b[j]*(tau-data->g[j]),2)-data->t[i]-2* tau* tau*data->b[j])*partial/ tau2;
        result[5]=result[5]+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*(data->t[i]-2* tau*data->b[j]*(tau-data->g[j]))*partial/ dt;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=(delta-1.0)*(delta-1.0);
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((tau-1),2));//equivalent to Psi
        z=(1-tau)+data->Af[j]*pow(d12,(0.5/data->betaf[j]));//equivalent to Theta
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);//equivalent to Dis
        Db=pow(D,data->bf[j]);//equivalent to Disb
        result[0]=result[0]+data->n[i]*Db* delta*p;//Arr

        dp_d=-2*data->Cf[j]*(delta-1)*p;
        if(!((delta-1)==0)){
            dz_d= data->Af[j]*pow(d12,0.5/data->betaf[j]-1)*(delta-1)/data->betaf[j];
            dD_d=2*z*dz_d+2*data->af[j]*data->Bf[j]*pow(d12,data->af[j]-1)*2*(delta-1);
            dDb_d=data->bf[j]*pow(D,(data->bf[j]-1))*dD_d;
        }
        else{
            dz_d=0.0;
            dD_d=0.0;
            dDb_d=0.0;
        }
        result[1]=result[1]+data->n[i]*(pow(D,data->bf[j])*(p+ delta*dp_d)+dDb_d* delta*p);//dArr/ddelta

        d2p_d2=(2*data->Cf[j]*d12-1)*2*data->Cf[j]*p;
        d2z_d2=data->Af[j]/data->betaf[j]*(1/data->betaf[j]-1)*pow(d12,0.5/data->betaf[j]-1);
        d2D_d2 = 2*(z*d2z_d2 + dz_d*dz_d + data->Bf[j]*(2*data->af[j]*data->af[j]-data->af[j])*pow(d12,data->af[j]-1));
        d2Db_d2 = data->bf[j]*(pow(D,data->bf[j]-1)*d2D_d2+(data->bf[j]-1)*pow(D,data->bf[j]-2.0)*pow(dD_d,2));
        result[2]=result[2]+data->n[i]*(Db*(2*dp_d+ delta*d2p_d2)+2*dDb_d*(p+ delta*dp_d)+d2Db_d2* delta*p);//d2Arr/ddelta2

        dp_t=-2*data->Df[j]*(tau-1)*p;
        dz_t=-1;
        dD_t=-2*z;
        dDb_t=-2*z*data->bf[j]*pow(D,(data->bf[j]-1));
        result[3]=result[3]+data->n[i]* delta*(dDb_t*p+pow(D,data->bf[j])*dp_t);//dArr/dtau

        d2p_t2=2*data->Df[j]*p*(2*data->Df[j]*pow((tau-1),2)-1);
        d2Db_t2=2*data->bf[j]*pow(D,(data->bf[j]-1))+4*pow(z,2)*data->bf[j]*(data->bf[j]-1)*pow(D,(data->bf[j]-2));
        result[4]=result[4]+data->n[i]* delta*(d2Db_t2*p+2*dDb_t*dp_t+pow(D,data->bf[j])*d2p_t2);//d2Arr/dtau2

        d2p_dt=data->Cf[j]*data->Df[j]*(2* delta-2)*(2* tau-2)*p;
        d2D_dt= -2* dz_d;
        if(!(D==0)) d2Db_dt=data->bf[j]* ((data->bf[j]-1)*pow(D,(data->bf[j]+1))*dD_d*dD_t + pow(D,(data->bf[j]+2))*d2D_dt)/pow(D,3);
        else d2Db_dt=0;
        result[5]=result[5]+ data->n[i]*(delta*Db*d2p_dt+ delta*p*d2Db_dt+ delta*dDb_d*dp_t+
                                         delta*dDb_t*dp_d + Db*dp_t + p*dDb_t);
    }
        //printf("%f  %f  %f %f %f %f\n",result[0],result[1],result[2],result[3],result[4],result[5]);
}

void CALLCONV FF_ArrDerSWTV(double T,double V,const  FF_SWEOSdata *data,double result[6])
{
    double delta,tau;
    double ArrDer[6];
    delta=1/(V * data->rhoRef);
    tau=data->tRef/ T;
    double V2=V * V;
    double T2=T * T;
    FF_ArrDerSW(tau,delta,data,ArrDer);
    result[0]=ArrDer[0];
    result[1]=-ArrDer[1]/(data->rhoRef* V2);
    result[2]=ArrDer[2]/(data->rhoRef* V2)/(data->rhoRef* V2)+ArrDer[1]*2/(data->rhoRef* V2 * V);
    result[3]=-ArrDer[3]*data->tRef/T2;
    result[4]=ArrDer[4]*data->tRef/T2*data->tRef/T2+ArrDer[3]*data->tRef*2/(T2 * T);
    result[5]=ArrDer[5]/(data->rhoRef* V2)*data->tRef/T2;
}

//Arr and dArr/dT for a pure substance, given tau and delta, according to SW EOS
void CALLCONV FF_ArrDerSW0T(double tau,double delta,const  FF_SWEOSdata *data,double result[2]){
    int i,j;
    double partial;
    double d12,D,Db,z,p;
    double dp_t,dDb_t;
    for (i=0;i<2;i++) result[i]=0;
    for (i=0;i<data->nPol;i++)//Polinomial terms
    {
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);//n*delta^d*tau^t
        result[0]=result[0]+partial;//This will be Arr
        result[1]=result[1]+data->t[i]*partial/ tau;//This will be dArr/dtau at constant delta
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)//exponential terms
    {
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]));//n*delta^d*tau^t*e^(-delta^c)
        result[0]=result[0]+partial;
        result[1]=result[1]+data->t[i]*partial/ tau;
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)//special exponential terms
    {
        j=i-data->nPol-data->nExp;
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        //n*delta^d*tau^t*e^(-a*(delta-e)^2-b*(tau-g)^2)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->t[i]-2* tau*data->b[j]*(tau-data->g[j]))*partial/ tau;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=pow((delta-1),2);
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((tau-1),2));//equivalent to Psi
        z=(1-tau)+data->Af[j]*pow(d12,(0.5/data->betaf[j]));//equivalent to Theta
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);//equivalent to Dis
        Db=pow(D,data->bf[j]);//equivalent to Disb
        result[0]=result[0]+data->n[i]*Db* delta*p;//Arr

        dp_t=-2*data->Df[j]*(tau-1)*p;
        dDb_t=-2*z*data->bf[j]*pow(D,(data->bf[j]-1));
        result[1]=result[1]+data->n[i]* delta*(dDb_t*p+pow(D,data->bf[j])*dp_t);//dArr/dtau
    }
        //printf("%f  %f  %f %f %f %f\n",result[0],result[1],result[2],result[3],result[4],result[5]);
}

//Arr and dArr/dV at constant T calculation for a pure substance, given T and V, according to Span and Wagner EOS
//---------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVsw(double T,double V,const  FF_SWEOSdata *data,double *Arr,double *Z)
{
    double dArr_d=0;
    *Arr=0;
    double tau=data->tRef/ T;
    double delta=1/(V*data->rhoRef);//rhoRef en mol/m3
    double powDeltaTau;
    double partial;
    int i,j;
    double d12,D,z,p;
    double dp_d,dz_d,dD_d,dDb_d;
    for (i=0;i<data->nPol;i++)
    {
        powDeltaTau=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);
        *Arr=*Arr+powDeltaTau;//This will be Arr
        dArr_d=dArr_d+data->d[i]*powDeltaTau/delta;//And this dArr/ddelta at constant T
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)
    {
        powDeltaTau=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]));

        *Arr=*Arr+powDeltaTau;
        dArr_d=dArr_d+powDeltaTau*(data->d[i]-data->c[i]*pow(delta,data->c[i]))/delta;
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)
    {
        j=i-data->nPol-data->nExp;
        //printf("j:%i\n",j);
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        *Arr=*Arr+partial;
        dArr_d=dArr_d+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*partial/delta;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=pow((delta-1),2);
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((tau-1),2));
        z=(1-tau)+data->Af[j]*pow(d12,(1/2/data->betaf[j]));
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);
        *Arr=*Arr+data->n[i]*pow(D,data->bf[j])* delta*p;//Arr
        dp_d=-2*data->Cf[j]*(delta-1)*p;
        dz_d= data->Af[j]*pow(d12,0.5/data->betaf[j]-1)*(delta-1)/data->betaf[j];
        dD_d=2*z*dz_d+2*data->af[j]*data->Bf[j]*pow(d12,data->af[j]-1)*2*(delta-1);
        dDb_d=data->bf[j]*pow(D,(data->bf[j]-1))*dD_d;
        dArr_d=dArr_d+data->n[i]*(pow(D,data->bf[j])*(p+ delta*dp_d)+dDb_d* delta*p);//dArr/ddelta
    }
    *Z=1+delta*dArr_d;//this is Z
}

//P calculation from T and V using according to Span and Wagner EOS
//-----------------------------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVsw(double T,double V,const  FF_SWEOSdata *data,double *P)
{
    double Z;
    double dArr_d=0;
    double tau=data->tRef/ T;
    double delta=1/(V*data->rhoRef);//rhoRef en mol/m3
    double powDeltaTau;
    double partial;
    unsigned int i,j;
    double d12,D,z,p;
    double dp_d,dz_d,dD_d,dDb_d;
    double nPol=data->nPol;
    double nExp=nPol+data->nExp;
    double nSpec=nExp+data->nSpec;
    double nFinal=nSpec+data->nFinal;
    for (i=0;i<nPol;i++)
    {
        powDeltaTau=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);
        dArr_d=dArr_d+data->d[i]*powDeltaTau/delta;//And this dArr/ddelta at constant T
    }
    for (i=nPol;i<nExp;i++)
    {
        powDeltaTau=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]));
        dArr_d=dArr_d+powDeltaTau*(data->d[i]-data->c[i]*pow(delta,data->c[i]))/delta;
    }
    for (i=nExp;i<nSpec;i++)
    {
        j=i-nExp;
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        dArr_d=dArr_d+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*partial/delta;
    }
    for (i=nSpec;i<nFinal;i++)//additional terms
    {
        j=i-nSpec;
        d12=(delta-1)*(delta-1);
        p=exp(-data->Cf[j]*d12-data->Df[j]*(tau-1)*(tau-1));
        z=(1-tau)+data->Af[j]*pow(d12,(0.5/data->betaf[j]));
        D=z*z+data->Bf[j]*pow(d12,data->af[j]);
        dp_d=-2*data->Cf[j]*(delta-1)*p;
        dz_d= data->Af[j]*pow(d12,0.5/data->betaf[j]-1)*(delta-1)/data->betaf[j];
        dD_d=2*z*dz_d+2*data->af[j]*data->Bf[j]*pow(d12,data->af[j]-1)*2*(delta-1);
        dDb_d=data->bf[j]*pow(D,(data->bf[j]-1))*dD_d;
        dArr_d=dArr_d+data->n[i]*(pow(D,data->bf[j])*(p+ delta*dp_d)+dDb_d* delta*p);//dArr/ddelta
        //printf("j:%i\n",j);
    }
    Z=1+delta*dArr_d;//this is Z
    *P=R* T*Z / V;//This is P in Pa
}

//P and dP/drho calculation for a pure substance, given reduced T and rho, according to SW EOS
//--------------------------------------------------------------------------------------------
void CALLCONV FF_PresDerSW(double tau,double delta,const  FF_SWEOSdata *data,double result[4])
{
    double delta2=delta*delta;
    double partial,aux,aux2;
    int i,j;
    double d12,D,Db,z,p;
    double dp_d,dz_d,dD_d,dDb_d;
    double d2z_d2,d2p_d2,d2D_d2,d2Db_d2;
    double nPol=data->nPol;
    double nExp=nPol+data->nExp;
    double nSpec=nExp+data->nSpec;
    double nFinal=nSpec+data->nFinal;
    for (i=0;i<4;i++) result[i]=0;
    for (i=0;i<nPol;i++)
    {
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);//n*delta^d*tau^t
        result[0]=result[0]+partial;//This will be Arr
        result[1]=result[1]+data->d[i]*partial/ delta;//This will be dArr/ddelta at constant tau
        result[2]=result[2]+data->d[i]*(data->d[i]-1)*partial/ delta2;//This will be d2Arr/ddelta2 at constant tau
    }
    for (i=nPol;i<nExp;i++)//ndelta^dtau^t*e^(-delta^c)
    {
        aux=pow(delta,data->c[i]);
        aux2=data->d[i]-data->c[i]*aux;
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-aux);//n*delta^d*tau^t*e^(-delta^c)
        result[0]=result[0]+partial;
        result[1]=result[1]+(aux2)*partial/ delta;
        result[2]=result[2]+((aux2)*(aux2-1)-data->c[i]*data->c[i]*aux)*partial/ delta2;
    }
    for (i=nExp;i<nSpec;i++)
    {
        j=i-nExp;
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        //n*delta^d*tau^t*e^(-a*(delta-e)^2-b*(tau-g)^2)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*partial/ delta;
        result[2]=result[2]+(pow(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]),2)-data->d[i]-2* delta* delta*data->a[j])*partial/ delta2;
    }
    for (i=nSpec;i<nFinal;i++)//additional terms
    {
        j=i-nSpec;
        d12=pow((delta-1),2);
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((tau-1),2));//equivalent to Psi
        z=(1-tau)+data->Af[j]*pow(d12,(0.5/data->betaf[j]));//equivalent to Theta
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);//equivalent to Dis
        Db=pow(D,data->bf[j]);//equivalent to Disb
        result[0]=result[0]+data->n[i]*Db* delta*p;//Arr

        dp_d=-2*data->Cf[j]*(delta-1)*p;
        if(!((delta-1)==0)){
            dz_d= data->Af[j]*pow(d12,0.5/data->betaf[j]-1)*(delta-1)/data->betaf[j];
            dD_d=2*z*dz_d+2*data->af[j]*data->Bf[j]*pow(d12,data->af[j]-1)*2*(delta-1);
            dDb_d=data->bf[j]*pow(D,(data->bf[j]-1))*dD_d;
        }
        else{
            dz_d=0.0;
            dD_d=0.0;
            dDb_d=0.0;
        }
        result[1]=result[1]+data->n[i]*(pow(D,data->bf[j])*(p+ delta*dp_d)+dDb_d* delta*p);//dArr/ddelta

        d2p_d2=(2*data->Cf[j]*d12-1)*2*data->Cf[j]*p;
        d2z_d2=data->Af[j]/data->betaf[j]*(1/data->betaf[j]-1)*pow(d12,0.5/data->betaf[j]-1);
        d2D_d2 = 2*(z*d2z_d2 + dz_d*dz_d + data->Bf[j]*(2*data->af[j]*data->af[j]-data->af[j])*pow(d12,data->af[j]-1));
        d2Db_d2 = data->bf[j]*(pow(D,data->bf[j]-1)*d2D_d2+(data->bf[j]-1)*pow(D,data->bf[j]-2.0)*pow(dD_d,2));
        result[2]=result[2]+data->n[i]*(Db*(2*dp_d+ delta*d2p_d2)+2*dDb_d*(p+ delta*dp_d)+d2Db_d2* delta*p);//d2Arr/ddelta2
    }
    //printf("Derivatives: %e %e %e\n",result[0],result[1],result[2]);
    result[3]=R*(data->tRef/ tau)*(1+2*delta*result[1]+delta2*result[2])*data->rhoRef;//this is dP/ddelta=R*T*(1+2*delta*dArr_d+delta*delta*d2Arr_d2
    result[2]=R*(data->tRef/ tau)*delta*data->rhoRef*(1+delta*result[1]);//This is P=R*T*rho*(1+delta*dArr_d)
    result[1]=1+delta*result[1]; //This is Z=1+delta*dArr_d
    //printf("Presion:%f derivada:%f\n",result[2],result[3]);
}


//Z,Arr,V and dP/dV at constant T calculation for a pure substance, given T and P, according to SW EOS
//----------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPsw(double T,double P,const  FF_SWEOSdata *data,char option,double resultL[3],double resultG[3],char *state)
{
    double Tc=data->Tc;
    double tau=data->tRef/ T;
    double delta,error;
    double deltaIni;
    double answer[4];
    int i;
    double maxError;
    if ((T>0.99*Tc)&&(T<1.01*Tc)) maxError=0.000001;
    else maxError=0.0001;

    double fw;
    FF_CubicParam param;
    param.b = 0.077796 * R * Tc / data->Pc;
    param.a = 5.877359*param.b*R*Tc;//0.457235 * pow(R*Tc,2)/ data->Pc;
    param.u=2.414214;//1+2^0.5
    param.w=-0.414214;//1-2^0.5
    param.c=0.0;
    fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
    param.Theta=param.a*pow((1 + fw * (1-pow((T / Tc),0.5))),2);
    char stateCubic;
    //printf("b:%f u:%f w:%f a:%f Theta:%f\n",param.b, param.u,param.w,param.a,param.Theta);
    FF_VfromTPcubic(T,P,&param,option,resultL,resultG,&stateCubic);

    *state='f';
    if ((option!='g')||(!(stateCubic=='b'))){
    //we calculate the liquid phase if not only the gas phase has been asked, or there is only one phase
        if ((T<=0.95*Tc)||(T>=1.05*Tc)) delta=1.35/(resultL[0]*data->rhoRef);//See how to use saturated density
        //else if ((T>=0.99*Tc)&&(T<1.01*Tc)&&(P>=0.9*data->Pc)&&(P<=1.1*data->Pc)) delta=2.0;
        //else if ((T>=0.97*Tc)&&(T<1.03*Tc)) delta=1.8/(resultL[0]*data->rhoRef);
        else if ((T>=0.95*Tc)&&(T<1.05*Tc)&&(P>=0.75*data->Pc)&&(P<=1.25*data->Pc)) delta=1.7;
        else delta=1.4/(resultL[0]*data->rhoRef);//we multiply by 1.4 the density obtained from PR EOS to begin the exploration
        FF_PresDerSW(tau,delta,data,answer);
        error =P-answer[2];
        i=1;
        //printf("i= %d  deltaL= %f rhoL:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*data->rhoRef*data->MW/1000,error,answer[0],answer[1]);
        while ((fabs(error/ P)>maxError) && (answer[3] > 0)&&(delta>1.0)&&(i<15))
        {   delta=delta+error/answer[3];
            FF_PresDerSW(tau,delta,data,answer);
            error =P-answer[2];
            i=i+1;
            //printf("i= %d  deltaL= %f rhoL:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*data->rhoRef*data->MW/1000,error,answer[0],answer[1]);
        }
        if ((answer[3]>=0)&&(fabs(error/ P)<=maxError))//if dP/ddelta >=0
        //if (((PtPlus-Pt)/(0.001*delta)>=0)&&(fabs(error/ P)<=maxError))//if dP/ddelta >=0
        {
            resultL[0]=1/(delta*data->rhoRef); //This is V
            resultL[1]=answer[0];
            resultL[2]=answer[1];
            //FF_ArrZfromTVsw(T,resultL[0],data,&resultL[1],&resultL[2]);
            *state='l';
        }
        else
        {
            //printf("fallo liquido\n");
            resultL[0]=resultL[1]=resultL[2]=0;
        }
        //printf("i= %d  deltaL= %f rhoL:%f error= %f  P= %f  dP_d= %f state:%c\n",i,delta,delta*data->rhoRef*data->MW/1000,error,answer[0],answer[1],*state);
    }

    if ((option!='l')||(!(stateCubic=='b')))//and the gas phase if not only the liquid one has been asked, or single value and not found as liquid
    {
        //if ((T>=0.99*Tc)&&(T<1.01*Tc)&&(P>=0.9*data->Pc)&&(P<=1.1*data->Pc)) deltaIni=0.4;
        if ((T>=0.99*Tc)&&(T<1.01*Tc)&&(P>=0.9*data->Pc)&&(P<=1.1*data->Pc)) deltaIni=0.5;
        //else deltaIni=0.8/(resultG[0]*data->rhoRef);//initial guess for gas volume
        else deltaIni=1.0/(resultG[0]*data->rhoRef);//initial guess for gas volume
        delta=deltaIni;
        //printf("Vg inicial:%f d:%f delta:%f\n",resultG[0],1/resultG[0],delta);
        //printf("tau= %f  delta= %f \n",tau,delta);
        FF_PresDerSW(tau,delta,data,answer);
        error =P-answer[2];
        i=1;
        //printf("i= %d  error= %f  dP_d= %f\n",i,error,answer[1]);
        while ((fabs(error/ P)>maxError) && (answer[3] >0)&&(i<55)&&(!((*state=='l')&&(1/(delta*data->rhoRef)<=resultL[0]))))
        {   delta=delta+error/answer[3];
            FF_PresDerSW(tau,delta,data,answer);
            error =P-answer[2];
            i=i+1;
            //printf("i= %d  deltaG= %f rhoG:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*data->rhoRef*data->MW/1000,error,answer[0],answer[1]);
        }
        if ((answer[3]>=0)&&(fabs(error/ P)<=maxError)&&(delta<5.0*deltaIni))
        {
            resultG[0]=1.0/(delta*data->rhoRef); //This is V
            resultG[1]=answer[0];
            resultG[2]=answer[1];
            //FF_ArrZfromTVsw(T,resultG[0],data,&resultG[1],&resultG[2]);
            if (*state=='l'){
                if (fabs((resultL[0]-resultG[0])/resultL[0])>0.01) *state='b';
                else *state='u';
            }
            else *state='g';
        }
        else
        {
            //printf("fallo gas\n");
            resultG[0]=resultG[1]=resultG[2]=0;
        }
        //printf("i= %d  deltaG= %f rhoG:%f error= %f  P= %f  dP_d= %f state:%c\n",i,delta,delta*data->rhoRef*data->MW/1000,error,answer[0],answer[1],*state);
    }
    if (((*state=='u')||(*state=='l'))&&((T>=Tc)||(1/resultG[0]<data->rhoRef))){
            resultG[0]=resultL[0];
            resultG[1]=resultL[1];
            resultG[2]=resultL[2];
            *state='G';
    }
    if (*state=='l') *state='L';
    else if  (*state=='g') *state='G';
    //printf("Final state:%c\n\n",*state);
}


//Z,Arr,V and dP/dV at constant T calculation for a pure substance, given T and P, according to SW EOS
//----------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPswS(double T,double P,const  FF_SubstanceData *data,char option,double resultL[3],double resultG[3],char *state)
{
    double Tc=data->swData.Tc-0.001;
    double Pc=data->swData.Pc-0.001;
    double rhoRef=data->swData.rhoRef;
    double tau=data->swData.tRef/ T;
    double delta,error;
    double gDeltaIni=0.0;
    double answer[4];
    int i;
    double maxError;//maximum relative error in pressure
    double Vp=0;//to detect its calculation
    double liqDensIni=0;
    double fw;
    char optionIn=option;
    FF_CubicParam param;
    char stateCubic;
    double rCubicL[3],rCubicG[3];

    //indicate situation of fail
    *state='f';
    resultL[0]=resultL[1]=resultL[2]=resultG[0]=resultG[1]=resultG[2]=0.0;

    //set maximum relative error for pressure
    if ((T>0.99*Tc)&&(T<1.01*Tc)) maxError=0.00005;
    else maxError=0.001;

    //set phase to search if possible
    if (T>=Tc) option='g';//over Tc only gas phase must be searched
    else if(P>=Pc) option='l';
    else if (option=='s'){//If the stable phase has been requested, this is the fastest method to stablish
        if (data->vpCorr.form>0){
            FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->swData.MW,T,&Vp);
                if (P<Vp*0.99) option='g';
                else if (P>Vp*1.01) option='l';
        }
    }
    //printf("vp:%i ld:%i lnuA:%f option:%c\n",data->vpCorr.form,data->lDensCorr.form,data->baseProp.LnuA,option);

    //set accurate initial liquid density if possible
    if(!(option=='g')&&(data->vpCorr.form>0)&&(data->lDensCorr.form>0)){
        double lDensSat;
        if (Vp==0) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->swData.MW,T,&Vp);
        FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->swData.MW,T,&lDensSat);
        if (data->baseProp.LnuA>0) liqDensIni=lDensSat+log(data->baseProp.LnuA*data->baseProp.MW*(P-Vp)/
                 (1000*R*T*exp(data->baseProp.LnuB+data->baseProp.LnuA*lDensSat))+1)/data->baseProp.LnuA;
        else FF_LiqDensChuehPrausnitz(&data->baseProp,T,P,Vp,lDensSat,&liqDensIni);
        //printf("lDensSat:%f liqDensIni:%f\n",lDensSat,liqDensIni);
    }
    //Cubic eos calculation to obtain start values, if not only liquid with better start value
    if (!((option=='l')&&(liqDensIni>0))){
        param.b = 0.077796 * R * Tc / data->swData.Pc;
        param.a = 5.877359*param.b*R*Tc;//0.457235 * pow(R*Tc,2)/ data->swData.Pc;
        param.u=2.414214;
        param.w=-0.414214;
        param.c=0.0;
        fw = 0.37464 + 1.54226 * data->swData.w - 0.26992 * data->swData.w * data->swData.w;
        param.Theta=param.a*pow((1 + fw * (1-pow((T / Tc),0.5))),2);
        //printf("b:%f u:%f w:%f a:%f Theta:%f\n",param.b, param.u,param.w,param.a,param.Theta);
        FF_VfromTPcubic(T,P,&param,option,rCubicL,rCubicG,&stateCubic);
        //param->c=R*data->Tc*(0.1014048-0.3892896*data->Zc)/data->Pc
    }

    //Zone where Newton method is slow
    if((T>=Tc)&&(T<Tc+2)&&(P>0.95*Pc)&&(P<1.05*Pc)){
        resultG[0]=FF_VsolverRegula(FF_PfromTVsw,&data->swData,T,P,0.6*rCubicG[0],1.2*rCubicG[0]);
        FF_ArrZfromTVsw(T,resultG[0],&data->swData,&resultG[1],&resultG[2]);//We need to fill Arr and Z
        *state='G';
        return;
    }

    //Newton method zone
    //we calculate the liquid phase if not only the gas phase has been asked, or there is only one phase
    if ((option!='g')){
        i=0;
        if (liqDensIni>0) delta=liqDensIni*1e3/(rhoRef*data->swData.MW);
        else if (T<=0.8*Tc) delta=1.2/(rCubicL[0]*rhoRef);//Bellow 0.8 Tc, we multiply by 1.25 the PR EOS density
        else if ((T>=0.95*Tc)) delta=2.0;
        else delta=1.3/(rCubicL[0]*rhoRef);//Between 0.8 and 0.95 Tc. We multiply by 1.4 the density obtained from PR EOS to begin the exploration
        FF_PresDerSW(tau,delta,&data->swData,answer);
        i++;
        if (answer[3]<0){//increase delta if it is too low
            delta*=1.1;
            FF_PresDerSW(tau,delta,&data->swData,answer);
            i++;
        }
        error =P-answer[2];
        //printf("i= %d  deltaL= %f rhoL:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*rhoRef*data->swData.MW/1000,error,answer[2],answer[3]);
        while ((fabs(error/ P)>maxError)&&(answer[3]>0)&&(delta>1.0)&&(i<10)){
            delta=delta+error/answer[3];
            FF_PresDerSW(tau,delta,&data->swData,answer);
            error =P-answer[2];
            i++;
            //printf("i= %d  deltaL= %f rhoL:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*rhoRef*data->swData.MW/1000,error,answer[2],answer[3]);
        }
        if ((fabs(error/ P)<=maxError)&&(answer[3]>0))//if dP/ddelta >=0
        {
            resultL[0]=1/(delta*rhoRef); //This is V
            resultL[1]=answer[0];//Arr
            resultL[2]=answer[1];//Z
            *state='l';
        }
        //if (*state=='f') printf("fallo liquido\n");
    }

    //Calculate gas phase if not only the liquid one has been asked
    if (option!='l'){
        i=0;
        /*if ((T>Tc)&&(P>Pc)) gDeltaIni=1.2/(rCubicG[0]*rhoRef);//initial guess for gas volume from cubic EOS, search from heavy side
        else gDeltaIni=0.9/(rCubicG[0]*rhoRef);//initial guess for gas volume from cubic EOS, search from light side*/
        if (T>Tc){//look for the best side to use in the one phase region
            if (P<=0.9*Pc) gDeltaIni=0.9/(rCubicG[0]*rhoRef);
            else if(P>1.2*Pc){
                if (data->id==399) gDeltaIni=1.4/(rCubicG[0]*rhoRef);
                else gDeltaIni=1.0/(rCubicG[0]*rhoRef);
            }
            else{
                double Paux;
                FF_PfromTVsw(T,1/data->swData.rhoRef,&data->swData,&Paux);//rhoRef to be sustitued by Vc
                if(P>Paux){
                    if (data->id==399) gDeltaIni=1.4/(rCubicG[0]*rhoRef);
                    else gDeltaIni=1.1/(rCubicG[0]*rhoRef);//initial guess for gas volume from cubic EOS, search from heavy side
                }
                else gDeltaIni=0.9/(rCubicG[0]*rhoRef);//initial guess for gas volume from cubic EOS, search from light side
            }
        }
        else gDeltaIni=0.9/(rCubicG[0]*rhoRef);
        delta=gDeltaIni;
        FF_PresDerSW(tau,delta,&data->swData,answer);
        i++;
        if (answer[3]<0){
            delta*=0.9;
            FF_PresDerSW(tau,delta,&data->swData,answer);
            i++;
        }
        error =P-answer[2];
        //printf("i= %d  deltaG= %f rhoG:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*rhoRef*data->swData.MW/1000,error,answer[2],answer[3]);
        while ((fabs(error/ P)>maxError) && (answer[3] >0)&&(i<10)&&(!((*state=='l')&&(1/(delta*rhoRef)<=resultL[0]))))
        {   delta=delta+error/answer[3];
            FF_PresDerSW(tau,delta,&data->swData,answer);
            error =P-answer[2];
            i++;
            //printf("i= %d  deltaG= %f rhoG:%f error= %f  P= %f  dP_d= %f\n",i,delta,delta*rhoRef*data->swData.MW/1000,error,answer[2],answer[3]);
        }
        if ((fabs(error/ P)<=maxError)&& (answer[3] >0)&&(delta<5.0*gDeltaIni)){
            resultG[0]=1.0/(delta*rhoRef); //This is V
            resultG[1]=answer[0];//Arr
            resultG[2]=answer[1];//Z
            if (*state=='l'){
                if (fabs((resultL[0]-resultG[0])/resultL[0])>0.01) *state='b';
                else *state='u';
            }
            else *state='g';
        }
        //if (!((*state=='g')||(*state=='u'))) (printf("fallo gas\n");
    }
    if(optionIn=='s'){
        if(*state=='l') *state='L';
        else if(*state=='g') *state='G';
    }

    //printf("Final state:%c\n\n",*state);
}


//Common P,V,T calculations
//=========================


//P calculation from T and V by eos
//---------------------------------
void CALLCONV FF_PfromTVeos(int eosType,double T,double V,void *data,double *P)
{
     FF_CubicParam param;
    switch (eosType)
    {
        case FF_IdealType:
            *P=R* T/ V;
            break;
        case FF_SAFTtype:
            //*( FF_SaftEOSdata*) data;
            FF_PfromTVSAFT(T,V,data,P);
            break;
        case FF_SWtype:
            FF_PfromTVsw(T,V,data,P);
            break;
        default://if we have a cubic eos the first step is to calculate its parameters
            //*( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(&T,data, &param);
            FF_PfromTVcubic(T,V,&param,P);
            break;
    }
}

//P calculation from T and V by eos, using a substance pointer
//------------------------------------------------------------
void CALLCONV FF_PfromTVeosS(double T,double V,const FF_SubstanceData *data,double *P)
{
     FF_CubicParam param;
    switch (data->model)
    {
        case FF_IdealType:
            *P=R* T/ V;
            break;
        case FF_SAFTtype:
            //*( FF_SaftEOSdata*) data;
            FF_PfromTVSAFT(T,V,&data->saftData,P);
            break;
        case FF_SWtype:
            FF_PfromTVsw(T,V,&data->swData,P);
            break;
        default://if we have a cubic eos the first step is to calculate its parameters
            //*( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(&data->cubicData, &param);
            FF_ThetaDerivCubic(&T,&data->cubicData, &param);
            FF_PfromTVcubic(T,V,&param,P);
            break;
    }
}


//V,Arr and Z calculation, given T and P by eos.
//----------------------------------------------
void CALLCONV FF_VfromTPeos(int eosType,double T,double P,const void *data,char option,double resultL[3],double resultG[3],char *state)
{
    //option= 'l':get liquid volume, 'g':gas, 'b':both, 's':determine also stable phase
    //state= 'u':just one solution found, 'l': solution found from liquid end, 'g':found from gas end, 'b':two solutions found, 'L':the stable phase is liquid, 'G': is gas
    //printf("P:%f",P);
     FF_CubicParam param;
    switch (eosType)//if we have a cubic eos the first step is to calculate its parameters
    {
        case FF_IdealType:
            resultL[0]=resultG[0]=R* T/ P;
            resultL[1]=resultG[1]=0;
            resultL[2]=resultG[2]=1;
            *state='u';
            break;
        case FF_SAFTtype:
            //( FF_SaftEOSdata*) data;
            FF_VfromTPsaft(T,P,data,option,resultL,resultG,state);
            break;
        case FF_SWtype:
            //*( FF_SWEOSdata*) data;
            FF_VfromTPsw(T,P,data,option,resultL,resultG,state);
            break;
        default:
            //( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(&T,data, &param);
            //printf("Hola soy V from TP eos: c:%f b:%f a:%f Theta:%f dTheta:%f d2Theta:%f\n",param.c,param.b,param.a,param.Theta,param.dTheta,param.d2Theta);
            FF_VfromTPcubic(T,P,&param,option,resultL,resultG,state);
            break;
    }
    //printf("T:%f  P:%f Vl:%f ArrL:%f Zl:%f\n",T,P,resultL[0],resultL[1],resultL[2]);
    //the answer from SAFT and SW is very elaborated to L or G, but from Cubic is u or b
    if (option=='s'){
        if (*state=='b'){
            if ((resultL[1]+resultL[2]-1-log(resultL[2]))<(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='L';//we compare Gdr
            else if ((resultL[1]+resultL[2]-1-log(resultL[2]))>(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='G';
            else *state='E';//if Gdr is the same we are in equilibrium
        }
    }
}

//V,Arr and Z calculation, given T and P by eos, using a substance pointer
//------------------------------------------------------------------------
void CALLCONV FF_VfromTPeosS(double T,double P,const FF_SubstanceData *data,char option,double resultL[3],double resultG[3],char *state)
{
    //option= 'l':get liquid volume, 'g':gas, 'b':both, 's':determine also stable phase
    //state= state='f': calculation fail, 'l': solution found from liquid end, 'g':found from gas end, 'b':two solutions found,
    //'u':just one solution found, 'L':the stable phase is liquid, 'G': is gas
    //printf("P:%f",P);
     FF_CubicParam param;
    switch (data->model)//if we have a cubic eos the first step is to calculate its parameters
    {
        case FF_IdealType:
            resultL[0]=resultG[0]=R* T/ P;
            resultL[1]=resultG[1]=0;
            resultL[2]=resultG[2]=1;
            *state='u';
            break;
        case FF_SAFTtype:
            //( FF_SaftEOSdata*) data;
            FF_VfromTPsaft(T,P, &data->saftData,option,resultL,resultG,state);
            break;
        case FF_SWtype:
            //*( FF_SWEOSdata*) data;
            FF_VfromTPswS(T,P,data,option,resultL,resultG,state);
            break;
        default:
            //( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(&data->cubicData, &param);
            FF_ThetaDerivCubic(&T,&data->cubicData, &param);
            //printf("Hola soy V from TP eos: c:%f b:%f a:%f Theta:%f dTheta:%f d2Theta:%f\n",param.c,param.b,param.a,param.Theta,param.dTheta,param.d2Theta);
            FF_VfromTPcubic(T,P,&param,option,resultL,resultG,state);
            break;
    }
    //printf("T:%f  P:%f Vl:%f ArrL:%f Zl:%f\n",T,P,resultL[0],resultL[1],resultL[2]);
    //the answer from SAFT and SW is very elaborated to L or G, but from Cubic is u or b
    if (option=='s'){
        if (*state=='b'){
            //printf("Calculating stable phase L:%f G:%f\n",(resultL[1]+resultL[2]-1-log(resultL[2])),(resultG[1]+resultG[2]-1-log(resultG[2])));
            if ((resultL[1]+resultL[2]-1-log(resultL[2]))<(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='L';//we compare Gdr
            else if ((resultL[1]+resultL[2]-1-log(resultL[2]))>(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='G';
            else *state='E';//if Gdr is the same we are in equilibrium
        }
    }
}

//Saturation point determination
//==============================


//Boiling point calculation
void CALLCONV FF_TbEos(int eosType,double P,const void *data,double *Tb)
{
    int n=0;//number of calculations done
    //printf("hola soy Tb\n");
    double phiMaxError;
    double Tc,Pc;
    int i;
    switch (eosType)//We read Tc,Pc and w. and in the case of IAPWS95 we calculate directly the Tb
    {
        case FF_SAFTtype:
            if ((( FF_SaftEOSdata*)data)->eos==FF_PCSAFTPOL1){
                *Tb=0;
                return;
            }
            Tc=(( FF_SaftEOSdata*)data)->Tc;
            Pc=(( FF_SaftEOSdata*)data)->Pc;
            break;
        case FF_SWtype:
            if (((( FF_SWEOSdata*)data)->eos==FF_IAPWS95)||((( FF_SWEOSdata*)data)->eos==IF97)){
                Pc=22.064e6;
                if (P>=Pc)
                {
                    *Tb=647.096;
                }
                else
                {
                    double beta=pow(P * 1e-6,0.25);
                    double B2=beta*beta;
                    double E=B2-0.17073846940092e2*beta+0.1491510861353e2;
                    double F=0.11670521452767e4*B2+0.1202082470247e5*beta-0.48232657361591e4;
                    double G=-0.72421316703206e6*B2-0.32325550322333e7*beta+0.40511340542057e6;
                    double D=2*G/(-F-pow(F*F-4*E*G,0.5));
                    *Tb=(0.65017534844798e3+D-pow((0.65017534844798e3+D)*(0.65017534844798e3+D)-4*(-0.23855557567849+0.65017534844798e3*D),0.5))/2;
                }
                return;
            }
            Tc=(( FF_SWEOSdata*)data)->Tc;
            Pc=(( FF_SWEOSdata*)data)->Pc;
            break;
        default://Cubic eos
            if (((( FF_CubicEOSdata*)data)->eos==FF_PRPOL1)||((( FF_CubicEOSdata*)data)->eos==FF_SRKPOL2)){
                *Tb=0;
                return;
            }
            Tc=(( FF_CubicEOSdata*)data)->Tc;
            Pc=(( FF_CubicEOSdata*)data)->Pc;
            break;
    }
    if ((P>=Pc)&&(!(Pc==0)))//If P>=Pc -> Tb=Tc
    {
        *Tb=Tc;
    }
    else//We need to calculate Tb
    {
        phiMaxError=0.0005;
        double T,phiL,phiG;//Temperature and fugacity coef.
        double answerL[3],answerG[3];
        char option='b',state;
        T=Tc/2;//Initial pressure guess
        i=4;//i will be used in the calculation of the new temperature
        //printf("Initial Tb guess:%f\n",T);
        FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
        //printf("T:%f P:%f  Vl:%f  Zl:%f  Vg:%f  Zg:%f state:%c\n",T,P,answerL[0],answerL[2],answerG[0],answerG[2],state);

        while (!(state=='b'))//Till we find a temperature with different possitive liquid and gas solutions
        {
            if ((state=='G')||(state=='g'))T=T-Tc/i;
            else if ((state=='L')||(state=='l'))T=T+Tc/i;
            FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
            i=i*2;
            n=n+1;
            //printf("Finding two phases n:%i i:%i T:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,T,answerL[0],answerL[2],answerG[0],answerG[2]);
            if ((n>15)||(state=='f')||(state=='u')){
                return;
            }
        }
        phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
        phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
        //printf("Initial phi calculation T:%f P:%f  phiL:%f  phiG:%f\n",*T,P,phiL,phiG);
        while (fabs(phiL-phiG)>phiMaxError)
        {
            if ((state=='G')||(state=='g'))T=T-Tc/i;
            else if ((state=='L')||(state=='l'))T=T+Tc/i;
            FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
            i=i*2;
            n=n+1;
            //printf("Finding two phases n:%i i:%i state:%c T:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,state,T,answerL[0],answerL[2],answerG[0],answerG[2]);
            if ((n>15)||(state=='f')||(state=='u')){
                return;
            }
        }
        *Tb=T;
    }
}

//Boiling point calculation
void CALLCONV FF_TbEosS(double P,const FF_SubstanceData *data,double *Tb)
{
    int ver=0;//controls the ejecution of the printf statements to provide information
    *Tb=0;
    int n=0;//number of calculations done
    //printf("hola soy Tb\n");
    double error, phiMaxError;
    double Tc,Pc;
    int i;
    switch (data->model)//We read Tc,Pc and w. and in the case of IAPWS95 we calculate directly the Tb
    {
        case FF_SAFTtype:
            if (data->saftData.eos==FF_PCSAFTPOL1){
                return;
            }
            Tc=data->saftData.Tc;//In order to allow for the EOS give a higher Tc than real
            Pc=data->saftData.Pc;
            break;
        case FF_SWtype:
            if ((data->swData.eos==FF_IAPWS95)||(data->swData.eos==IF97)){
                Pc=22.064e6;
                if (P>=Pc)
                {
                    *Tb=647.096;
                }
                else
                {
                    double beta=pow(P * 1e-6,0.25);
                    double B2=beta*beta;
                    double E=B2-0.17073846940092e2*beta+0.1491510861353e2;
                    double F=0.11670521452767e4*B2+0.1202082470247e5*beta-0.48232657361591e4;
                    double G=-0.72421316703206e6*B2-0.32325550322333e7*beta+0.40511340542057e6;
                    double D=2*G/(-F-pow(F*F-4*E*G,0.5));
                    *Tb=(0.65017534844798e3+D-pow((0.65017534844798e3+D)*(0.65017534844798e3+D)-4*(-0.23855557567849+0.65017534844798e3*D),0.5))/2;
                }
                return;
            }
            Tc=data->swData.Tc;
            Pc=data->swData.Pc;
            break;
        default://Cubic eos
            if ((data->cubicData.eos==FF_PRPOL1)||(data->cubicData.eos==FF_SRKPOL2)){
                return;
            }
            Tc=data->cubicData.Tc;
            Pc=data->cubicData.Pc;
            break;
    }
    if((P>=Pc)&&(!(Pc==0))){
        *Tb=Tc;
        return;
    }
    if (data->model==FF_SAFTtype) Tc=1.1*Tc;//In SAFT the Tc reported by the eos can be much higher than the real.

    //We continue with the calculation
    phiMaxError=0.0003;
    double T,phiL,phiG;//Temperature and fugacity coef.
    double answerL[3],answerG[3];
    char option='b',state;

    T=Tc/2;//Initial pressure guess
    i=2;//i will be used in the calculation of the new temperature
    if (ver==1) printf("P%f Initial Tb guess:%f \n",P,T);
    FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
    n++;
    if (ver==1) printf("T:%f P:%f  Vl:%f  Zl:%f  Vg:%f  Zg:%f state:%c \n",T,P,answerL[0],answerL[2],answerG[0],answerG[2],state);

    while (!(state=='b'))//Till we find a temperature with different possitive liquid and gas solutions
    {
        i=i*2;
        if ((state=='u')||(n>13)){
            *Tb=T;
            if (ver==1) printf("Tb calculation finished at n:%i state:%c T:%f \n",n,state,T);
            return;
        }
        else if ((state=='f')||(n>20)){
            if (ver==1) printf("Tb calculation failed n:%i state:%c \n",n,state);
            return;
        }
        else if ((state=='G')||(state=='g'))T=T-Tc/i;
        else if ((state=='L')||(state=='l'))T=T+Tc/i;
        FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
        n++;
        if (ver==1) printf("Tb Finding two phases n:%i i:%i state:%c T:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i,state,T,answerL[0],answerL[2],answerG[0],answerG[2]);

    }
    phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
    phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
    error=phiL-phiG;
    if (ver==1) printf("Initial phi calculation T:%f P:%f  phiL:%e  phiG:%e error:%f\n",T,P,phiL,phiG,error);
    if(fabs(error)>0.01){
        while (fabs(error)>phiMaxError)
        {
            T=T*pow((phiG / phiL),0.05);
            FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
            n++;
            phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
            phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
            error=phiL-phiG;
            if (ver==1) printf("n:%i Tb:%f state:%c phiL:%f  phiG:%f error:%f \n",n,T, state, phiL,phiG,error);
            if (n>50){
                return;
            }
        }
    }
    else{
        while (fabs(phiL-phiG)>phiMaxError)
        {
            i=i*2;
            if(phiL>phiG) T=T-Tc/i;
            else T=T+Tc/i;
            FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
            n++;
            while (!(state=='b'))//Till we find a temperature with different possitive liquid and gas solutions
            {
                i=i*2;
                if ((state=='G')||(state=='g'))T=T-Tc/i;
                else if ((state=='L')||(state=='l'))T=T+Tc/i;
                FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
                n++;
                if (ver==1) printf("Finding two phases again n:%i i:%i %c P:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,state,P,answerL[0],answerL[2],answerG[0],answerG[2]);
                if ((n>30)||(state=='f')||(state=='u')){
                    return;
                }
            }
            phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
            phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
            if (ver==1) printf("n:%i T:%f  phiL:%f  phiG:%f \n",n,T,phiL,phiG);
            if (n>20){
                return;
            }
        }
    }
    *Tb=T;
}


//Vapor pressure calculation
void CALLCONV FF_VpEos(int eosType,double T,const void *data,double *Vp)
{
    //printf("Hola, soy Vp, EOS type:%i\n",eosType);
    int n=0;//number of calculations done
    double error,phiMaxError;
     double Tc,Pc;
    int i;
    //printf("Eos type:%i\n",eosType);

    *Vp=0;
    switch (eosType)
    {
        case FF_SAFTtype:
            if ((( FF_SaftEOSdata*)data)->eos==FF_PCSAFTPOL1){
                return;
            }
            Tc=(( FF_SaftEOSdata*)data)->Tc;
            Pc=(( FF_SaftEOSdata*)data)->Pc;
            //printf("Tc:%f Pc:%f sigma:%f m:%f epsilon:%f \n",Tc,Pc,(( FF_SaftEOSdata*)data)->sigma,(( FF_SaftEOSdata*)data)->m,(( FF_SaftEOSdata*)data)->epsilon);
            break;
        case FF_SWtype:
            if (((( FF_SWEOSdata*)data)->eos==FF_IAPWS95)||((( FF_SWEOSdata*)data)->eos==IF97)){
                Tc=647.096;
                if (T>=Tc)
                {
                    *Vp=22.064e6;
                }
                else
                {
                    double tau=T-0.23855557567849/(T-0.65017534844798e3);
                    double tau2=tau*tau;
                    double A=tau2+0.11670521452767e4*tau-0.72421316703206e6;
                    double B=-0.17073846940092e2*tau2+0.1202082470247e5*tau-0.32325550322333e7;
                    double C=0.1491510861353e2*tau2-0.48232657361591e4*tau+0.40511340542057e6;
                    *Vp=pow(2*C/(-B+pow(B*B-4*A*C,0.5)),4)*1e6;
                }
                return;
            }
            Tc=(( FF_SWEOSdata*)data)->Tc;
            Pc=(( FF_SWEOSdata*)data)->Pc;
            //printf("EOS:%i\n",(( FF_SWEOSdata*)data)->eos);
            break;
        default://Cubic eos
        {   enum FF_EOS tipo=(( FF_CubicEOSdata*)data)->eos;
            //printf("estoy en Vp, a ver que esos es:%i\n",tipo);
            if (((( FF_CubicEOSdata*)data)->eos==FF_PRPOL1)||((( FF_CubicEOSdata*)data)->eos==FF_SRKPOL2)){
                return;
            }
            Tc=(( FF_CubicEOSdata*)data)->Tc;
            Pc=(( FF_CubicEOSdata*)data)->Pc;
            //printf("Hola, soy Vp cubic eos:%i T:%f Tc:%f Pc:%f w:%f k1:%f\n",eosType,T,Tc,Pc,(( FF_CubicEOSdata*)data)->w,(( FF_CubicEOSdata*)data)->k1);
            break;
        }
    }
    if ((T>=Tc)&&(!(Tc==0)))//If T> supplied Tc*1.01 no calculation is made
    {
        *Vp=Pc;//Temperature over Tc
        return;
    }
    else//We need to calculate Vp
    {
        if(T>Tc-0.5) phiMaxError=0.00003;
        //else if(T>0.99*Tc) phiMaxError=0.00001;
        else phiMaxError=0.00001;
        double P,phiL,phiG;//Pressure and fugacity coef.
        double answerL[3],answerG[3];
        char option='b',state;
        P=Pc*1.04;
        P=Pc/2;//Initial pressure guess
        i=4;//i will be used in the calculation of the new temperature
        //printf("Initial Vp guess:%f\n",P);
        FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
        //printf("T:%f P:%f  Vl:%f  Zl:%f  Vg:%f  Zg:%f state:%c\n",T,P,answerL[0],answerL[2],answerG[0],answerG[2],state);

        while (!(state=='b'))//Till we find a pressure with different possitive liquid and gas solutions
        {
            if ((state=='G')||(state=='g'))P=P+Pc/i;
            else if ((state=='L')||(state=='l'))P=P-Pc/i;
            FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
            i=i*2;
            n=n+1;
            //printf("Finding two phases n:%i i:%i P:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,P,answerL[0],answerL[2],answerG[0],answerG[2]);
            if ((n>20)||(state=='f')||(state=='u')){
                return;
            }
        }
        phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
        phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
        //printf("Initial phi calculation T:%f P:%f  phiL:%f  phiG:%f\n",T,P,phiL,phiG);
        error=phiL-phiG;
        if(T<Tc-10.0){
            while (fabs(error)>phiMaxError)
            {
                P=P * phiL / phiG;
                FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
                phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
                phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
                error=phiL-phiG;
                n++;
                if(n==15) phiMaxError*=2;
                if(n==19) phiMaxError*=2;
                if(n==23) phiMaxError*=2;
                if(n==27) phiMaxError*=2;
                if(n==31) phiMaxError*=2;
                //printf("%i Pb:%f  phiL:%f  phiG:%f error:%f\n",n-1,P,phiL,phiG,error);
                if (n>35) return;
            }
        }
        else{
            while (fabs(error)>phiMaxError)
            {
                if (phiL>phiG) P=P+Pc/i;
                else P=P-Pc/i;
                FF_VfromTPeos(eosType,T,P,data,option,answerL,answerG,&state);
                phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
                phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
                error=phiL-phiG;
                i=i*2;
                n++;
                if(n==15) phiMaxError*=2;
                if(n==19) phiMaxError*=2;
                if(n==23) phiMaxError*=2;
                if(n==27) phiMaxError*=2;
                if(n==31) phiMaxError*=2;
                //printf("n:%i Pb:%f  phiL:%f  phiG:%f error:%f\n",n,P,phiL,phiG,error);
                if(n>35) return;
            }
        }
        *Vp=P;
    }
}

//Vapor pressure calculation
void CALLCONV FF_VpEosS(double T,const FF_SubstanceData *data,double *Vp)
{
    int ver=0;//controls the ejecution of the printf statements to provide information
    *Vp=0;
    int n=0;//number of calculations done
    double error,phiMaxError;
    double Tc,Pc;
    int i;
    switch (data->model)
    {
        case FF_SAFTtype:
            if (data->saftData.eos==FF_PCSAFTPOL1){
                return;
            }
            Tc=data->saftData.Tc;
            Pc=data->saftData.Pc;
            break;
        case FF_SWtype:
            if ((data->swData.eos==FF_IAPWS95)||(data->swData.eos==IF97)){
                Tc=647.096;
                if (T>=Tc)
                {
                    *Vp=22.064e6;
                }
                else
                {
                    double tau=T-0.23855557567849/(T-0.65017534844798e3);
                    double tau2=tau*tau;
                    double A=tau2+0.11670521452767e4*tau-0.72421316703206e6;
                    double B=-0.17073846940092e2*tau2+0.1202082470247e5*tau-0.32325550322333e7;
                    double C=0.1491510861353e2*tau2-0.48232657361591e4*tau+0.40511340542057e6;
                    *Vp=pow(2*C/(-B+pow(B*B-4*A*C,0.5)),4)*1e6;
                }
                return;
            }
            Tc=data->swData.Tc;
            Pc=data->swData.Pc;
            break;
        default://Cubic eos
        {
            if ((data->cubicData.eos==FF_PRPOL1)||(data->cubicData.eos==FF_SRKPOL2)){
                return;
            }
            Tc=data->cubicData.Tc;
            Pc=data->cubicData.Pc;
            if (ver==1) printf("Tc cubic:%f \n",Tc);
            break;
        }
    }
    if ((T>=Tc)&&(!(Tc==0))){
        *Vp=Pc;
        return;
    }
    if (data->model==FF_SAFTtype) Pc=1.3*Pc;//In SAFT the Pc reported by the eos can be much higher than the real.
    //Vp calculation
    //if(T>Tc-0.5) phiMaxError=0.00005;
    //else phiMaxError=0.00001;
    if(T>0.85*Tc) phiMaxError=0.00006;
    else if(T>0.95*Tc) phiMaxError=0.0001;
    else phiMaxError=0.00001;
    double P,phiL,phiG;//Pressure and fugacity coef.
    double answerL[3],answerG[3];
    char option='b',state;
    P=Pc/2;//Initial pressure guess
    i=2;//i will be used in the calculation of the new temperature
    if (ver==1) printf("T:%f Initial Vp guess:%f\n",T,P);
    FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
    n++;
    if (ver==1) printf("T:%f P:%f  Vl:%f  Zl:%f  Vg:%f  Zg:%f state:%c\n",T,P,answerL[0],answerL[2],answerG[0],answerG[2],state);

    while (!(state=='b'))//Till we find a pressure with different possitive liquid and gas solutions
    {
        i=i*2;
        if (state=='u'){
            *Vp=P;
            return;
        }
        else if ((state=='G')||(state=='g'))P=P+Pc/i;
        else if ((state=='L')||(state=='l'))P=P-Pc/i;
        FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
        n=n+1;
        if (ver==1) printf("Vp Finding two phases n:%i i:%i %c T:%f P:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,state, T, P,answerL[0],answerL[2],answerG[0],answerG[2]);
        //if ((n>20)||(state=='f')||(state=='u')){
        if ((n>20)||(state=='f')){
            //printf("Vp calculation failed n:%i state:%c \n",n,state);
            return;
        }
    }
    phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
    phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
    error=phiL-phiG;
    if (ver==1) printf("Initial phi calculation T:%f P:%f  phiL:%f  phiG:%f error:%f\n",T,P,phiL,phiG,error);

    if(fabs(error)>0.01){
        while (fabs(error)>phiMaxError){
               P=P *fmax(fmin((phiL / phiG),5),0.2);
              FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
              n++;
              phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
              phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
              error=phiL-phiG;
              if(n==15) phiMaxError*=2;
              if(n==19) phiMaxError*=2;
              if(n==23) phiMaxError*=2;
              if(n==27) phiMaxError*=2;
              if(n==31) phiMaxError*=2;
              if (ver==1) printf("n:%i Pv:%f  state:%c Vl:%f phiL:%f  Vg:%f phiG:%f error:%f\n",n,P,state,answerL[0],phiL,answerG[0],phiG,error);
              if (n>35) return;
          }
    }
    else{
        while ((fabs(error)>phiMaxError)&&(n<15)){
            i=i*2;
            if (phiL>phiG) P=P+Pc/i;
            else P=P-Pc/i;
            FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
            n=n+1;
            while (!(state=='b'))//Till we find a pressure with different possitive liquid and gas solutions
            {
                i=i*2;
                if ((state=='G')||(state=='g'))P=P+Pc/i;
                else if ((state=='L')||(state=='l'))P=P-Pc/i;
                FF_VfromTPeosS(T,P,data,option,answerL,answerG,&state);
                n++;
                if (ver==1) printf("Finding two phases again n:%i i:%i %c P:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,state,P,answerL[0],answerL[2],answerG[0],answerG[2]);
                if (n>15) break;
                else if ((state=='f')||(state=='u')){
                    return;
                }
            }
        phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
        phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
        error=phiL-phiG;
        if (ver==1) printf("n:%i Pv:%f  state:%c phiL:%f  phiG:%f error:%f\n",n,P,state,phiL,phiG,error);
        }
    }
    if (ver==1) printf("Vp final:%f \n",P);
    *Vp=P;
}



//Thermodynamic properties calculation from T and V
//=================================================

//Thermodynamic properties calculation for a ideal gas at same T and V, from a reference state, specified by refT and refP, where H and S are 0
//---------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_IdealThermoEos(int equation,const double coef[],double refT,double refP, FF_ThermoProperties *th0)
{   //the result delivered is always in J/mol·K
    int i;
    th0->P=R*th0->T/th0->V;//this is the pressure of a ideal gas with the same volume
    //printf("T:%f P:%f refT:%f refP:%f\n",th0->T,th0->P,refT,refP);
    if(refT==0){//Direct integration without reference T
        switch (equation)
        {
            case 1://DIPPR 100 in KJ/kgr·K
            case 2://DIPPR 100 in J/mol·K
            case 6://DIPPR 100 in J/kgr·K
                th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4);
                th0->H=coef[0]*th0->T+coef[1]*(th0->T*th0->T)/2+coef[2]*(th0->T*th0->T*th0->T)/3+
                        coef[3]*(th0->T*th0->T*th0->T*th0->T)/4+
                        coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T)/5;//This is the integration at actual T at constant pressure
                        //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
                th0->S=coef[0]*log(th0->T)+coef[1]*(th0->T)+coef[2]*(th0->T*th0->T)/2+coef[3]*(th0->T*th0->T*th0->T)/3+coef[4]*(th0->T*th0->T*th0->T*th0->T)/4;
                break;
            case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 in cal/(mol·K)
                th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4)+coef[5]*pow(th0->T,5);
                th0->H=coef[0]*(th0->T)+coef[1]*(th0->T*th0->T)/2+coef[2]*(th0->T*th0->T*th0->T)/3+ coef[3]*(th0->T*th0->T*th0->T*th0->T)/4+
                    coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T)/5+
                    coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T*th0->T)/6;//This is the integration at actual T at constant pressure
                    //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
                th0->S=coef[0]*log(th0->T)+coef[1]*(th0->T)+coef[2]*(th0->T*th0->T)/2+coef[3]*(th0->T*th0->T*th0->T)/3+
                        coef[4]*(th0->T*th0->T*th0->T*th0->T)/4+
                    coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T)/5;
                break;
            case 3://DIPPR 107 correlation in calories/mol·K
            case 4://DIPPR 107 correlation in J/Kmol·K
            case 200://DIPPR 107 correlation in J/kg·K
                th0->Cp=(coef[0]+coef[1]*pow(coef[2]/th0->T/sinh(coef[2]/th0->T),2)+coef[3]*pow(coef[4]/th0->T/cosh(coef[4]/th0->T),2));
                th0->H=(coef[0]*(th0->T)+coef[1]*coef[2]*(1/tanh(coef[2]/th0->T)-1/tanh(coef[2]))+coef[3]*coef[4]*(tanh(coef[4])-tanh(coef[4]/th0->T)));

                th0->S=(coef[0]*log(th0->T)+coef[1]*(coef[2]/th0->T/tanh(coef[2]/th0->T)-log(sinh(coef[2]/th0->T)))-
                        coef[3]*(coef[4]/th0->T*tanh(coef[4]/th0->T)-log(cosh(coef[4]/th0->T))));
                break;
            case 9:{//ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
                th0->Cp=coef[0]+exp(coef[1]/th0->T+coef[2]+coef[3]*th0->T+coef[4]*pow(th0->T,2));
                th0->H=coef[0]*(th0->T- refT);
                th0->S=coef[0]*log(th0->T/ refT);
                int i=20;
                double interval=(th0->T- refT)/i;
                double T,Cp1,Cp2;
                Cp1=exp(coef[1]/ refT+coef[2]+coef[3]* refT+coef[4]*pow(refT,2));
                for (i=1;i<21;i++){
                    T=refT+i*interval;
                    Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                    th0->H=th0->H+(Cp1+Cp2)/2*interval;
                    th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                    Cp1=Cp2;
                    }
                }
                break;
            case 5:{//Wilhoit equation J/mol·K (8 coefficients)
                double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
                a7_a6=coef[7]/coef[6];
                a7_a6_2=a7_a6*a7_a6;
                a7_a6_4=a7_a6_2*a7_a6_2;
                x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);

                if (th0->T<=coef[7]) y=0;
                else y=(th0->T-coef[7])/(th0->T+coef[6]);
                y2=y*y;
                y4=y2*y2;
                th0->Cp=R*(coef[0] + coef[1]/pow(th0->T,2)*exp(-coef[2]/th0->T) + coef[3]*y2 + (coef[4] - coef[5] /pow((th0->T -coef[7]),2))*y4*y4);

                if (th0->T<=coef[7]) h=0;
                else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                        coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
                th0->H= R*th0->T*(coef[0]+coef[1]*exp(-coef[2]/th0->T)/(coef[2]*th0->T))+R*h;

                if (th0->T<=coef[7]) s=0;
                else{
                    z = th0->T*(coef[7] + coef[6])/((th0->T + coef[6])*coef[7]);
                    w=0;
                    for (i=1;i<8;i++) w=w+(x1*pow(-a7_a6,6-i) - coef[4])*pow(y,i)/i;
                    s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                        (coef[3] + coef[4])*log((th0->T + coef[6])/(coef[6] + coef[7]))-
                        (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
                }
                th0->S=R*(coef[0]*log(th0->T)+coef[1]*(1+coef[2]/th0->T)*exp(-coef[2]/th0->T)/(coef[2]*coef[2])+s);

                }
                break;
            case 7://Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
                th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
                for (i=3;i<13;i=i+2){
                    if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
                }
                th0->Cp=th0->Cp*R;
                th0->H=coef[0]*(th0->T)+coef[1]/(coef[2]+1)*(pow(th0->T,(coef[2]+1)));
                for (i=3;i<13;i=i+2) if (coef[i]>0) th0->H=th0->H+coef[i]*coef[i+1]*(1/(exp(coef[i+1]/th0->T)-1));
                th0->H=th0->H*R;
                if (coef[1]>0) th0->S=coef[0]*log(th0->T)+coef[1]/coef[2]*(pow(th0->T,coef[2]));
                else th0->S=coef[0]*log(th0->T);
                //printf("T:%f S0(0):%f\n",th0->T,th0->S*R/th0->MW);
                for (i=3;i<13;i=i+2){
                    if (coef[i]>0) th0->S=th0->S+coef[i]*coef[i+1]*(exp(coef[i+1]/th0->T)/(th0->T*(exp(coef[i+1]/th0->T)-1))-
                        log(exp(coef[i+1]/th0->T)-1)/coef[i+1]);
                    //printf("S0(%i):%f\n",i,th0->S*R/th0->MW);
                }
                th0->S=th0->S*R;
                break;
            case 8://Jaeschke and Schley equation (9 coefficients). Used by GERG2004
                th0->Cp=1+coef[0];
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow(coef[i+1]/th0->T/sinh(coef[i+1]/th0->T),2)+coef[i+2]*pow(coef[i+3]/th0->T/cosh(coef[i+3]/th0->T),2);
                th0->Cp=th0->Cp*R;
                th0->H=(1+coef[0])*(th0->T);
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->H=th0->H+2*coef[i]*coef[i+1]*(1/(exp(2*coef[i+1]/th0->T)-1))+
                        +2*coef[i+2]*coef[i+3]*(1/(exp(2*coef[i+3]/th0->T)+1));
                th0->H=th0->H*R;
                th0->S=(1+coef[0])*log(th0->T);
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->S=th0->S+coef[i]*(coef[i+1]/th0->T/tanh(coef[i+1]/th0->T)-log(sinh(coef[i+1]/th0->T)))-
                        coef[i+2]*(coef[i+3]/th0->T*tanh(coef[i+3]/th0->T)-log(cosh(coef[i+3]/th0->T)));
                th0->S=th0->S*R;
                break;
            default:
            th0->Cp=th0->H=th0->S=0;
        }
    }
    else{//Integration from a reference T
        switch (equation)
        {
            case 1://DIPPR 100 in KJ/kgr·K
            case 2://DIPPR 100 in J/mol·K
            case 6://DIPPR 100 in J/kgr·K
                th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4);
                th0->H=coef[0]*(th0->T-refT)+coef[1]*(th0->T*th0->T-refT* refT)/2+coef[2]*(th0->T*th0->T*th0->T-refT* refT* refT)/3+
                        coef[3]*(th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT)/4+
                        coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT* refT)/5;//This is the integration from reference T to actual T at constant pressure
                        //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
                th0->S=coef[0]*log(th0->T/ refT)+coef[1]*(th0->T- refT)+coef[2]*(th0->T*th0->T- refT* refT)/2+coef[3]*(th0->T*th0->T*th0->T-
                        refT* refT* refT)/3+coef[4]*(th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT)/4;
                break;
            case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 in cal/(mol·K)
                th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4)+coef[5]*pow(th0->T,5);
                th0->H=coef[0]*(th0->T-refT)+coef[1]*(th0->T*th0->T-refT* refT)/2+coef[2]*(th0->T*th0->T*th0->T-refT* refT* refT)/3+
                    coef[3]*(th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT)/4+
                    coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT* refT)/5+
                    coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT* refT* refT)/6;//This is the integration from reference T to actual T at constant pressure
                    //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
                th0->S=coef[0]*log(th0->T/ refT)+coef[1]*(th0->T- refT)+coef[2]*(th0->T*th0->T- refT* refT)/2+coef[3]*(th0->T*th0->T*th0->T-
                    refT* refT* refT)/3+coef[4]*(th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT)/4+
                    coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T-refT* refT* refT* refT* refT)/5;
                break;
            case 3://DIPPR 107 correlation in calories/mol·K
            case 4://DIPPR 107 correlation in J/Kmol·K
            case 200://DIPPR 107 correlation in J/kg·K
                th0->Cp=(coef[0]+coef[1]*pow(coef[2]/th0->T/sinh(coef[2]/th0->T),2)+coef[3]*pow(coef[4]/th0->T/cosh(coef[4]/th0->T),2));
                th0->H=(coef[0]*(th0->T-refT)+coef[1]*coef[2]*(1/tanh(coef[2]/th0->T)-1/tanh(coef[2]/ refT))+coef[3]*coef[4]*(tanh(coef[4]/ refT)-tanh(coef[4]/th0->T)));
                th0->S=(coef[0]*log(th0->T/ refT)+coef[1]*(coef[2]/th0->T/tanh(coef[2]/th0->T)-log(sinh(coef[2]/th0->T)))-coef[1]*(coef[2]/
                        refT/tanh(coef[2]/ refT)-log(sinh(coef[2]/ refT)))-coef[3]*(coef[4]/th0->T*tanh(coef[4]/th0->T)-log(cosh(coef[4]/th0->T)))+
                        coef[3]*(coef[4]/ refT*tanh(coef[4]/ refT)-log(cosh(coef[4]/ refT))));
                break;
            case 9:{//ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
                th0->Cp=coef[0]+exp(coef[1]/th0->T+coef[2]+coef[3]*th0->T+coef[4]*pow(th0->T,2));
                th0->H=coef[0]*(th0->T- refT);
                th0->S=coef[0]*log(th0->T/ refT);
                int i=20;
                double interval=(th0->T- refT)/i;
                double T,Cp1,Cp2;
                Cp1=exp(coef[1]/ refT+coef[2]+coef[3]* refT+coef[4]*pow(refT,2));
                for (i=1;i<21;i++){
                    T=refT+i*interval;
                    Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                    th0->H=th0->H+(Cp1+Cp2)/2*interval;
                    th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                    Cp1=Cp2;
                    }
                }
                break;
            case 5:{//Wilhoit equation J/mol·K (8 coefficients)
                double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
                a7_a6=coef[7]/coef[6];
                a7_a6_2=a7_a6*a7_a6;
                a7_a6_4=a7_a6_2*a7_a6_2;
                x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);

                if (th0->T<=coef[7]) y=0;
                else y=(th0->T-coef[7])/(th0->T+coef[6]);
                y2=y*y;
                y4=y2*y2;
                th0->Cp=R*(coef[0] + coef[1]/pow(th0->T,2)*exp(-coef[2]/th0->T) + coef[3]*y2 + (coef[4] - coef[5] /pow((th0->T -coef[7]),2))*y4*y4);

                if (th0->T<=coef[7]) h=0;
                else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                        coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
                th0->H= R*th0->T*(coef[0]+coef[1]*exp(-coef[2]/th0->T)/(coef[2]*th0->T))+R*h;

                if (th0->T<=coef[7]) s=0;
                else{
                    z = th0->T*(coef[7] + coef[6])/((th0->T + coef[6])*coef[7]);
                    w=0;
                    for (i=1;i<8;i++) w=w+(x1*pow(-a7_a6,6-i) - coef[4])*pow(y,i)/i;
                    s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                        (coef[3] + coef[4])*log((th0->T + coef[6])/(coef[6] + coef[7]))-
                        (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
                }
                th0->S=R*(coef[0]*log(th0->T)+coef[1]*(1+coef[2]/th0->T)*exp(-coef[2]/th0->T)/(coef[2]*coef[2])+s);

                if (refT<coef[7]) y=0;
                else y=(refT-coef[7])/(refT+coef[6]);
                y2=y*y;
                y4=y2*y2;
                if (refT<=coef[7]) h=0;
                else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                        coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
                th0->H=th0->H-R* refT*(coef[0]+coef[1]*exp(-coef[2]/ refT)/(coef[2]* refT))-R*h;

                if (refT<=coef[7]) s=0;
                else{
                    z = refT*(coef[7] + coef[6])/((refT + coef[6])*coef[7]);
                    w=0;
                    for (i=1;i<8;i++) w=w+(x1*pow(-a7_a6,6-i) - coef[4])*pow(y,i)/i;
                    s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                            (coef[3] + coef[4])*log((refT + coef[6])/(coef[6] + coef[7]))-
                            (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
                }
                th0->S=th0->S-R*(coef[0]*log(refT)+coef[1]*(1+coef[2]/ refT)*exp(-coef[2]/ refT)/(coef[2]*coef[2])+s);
                }
                break;
            case 7://Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
                th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
                for (i=3;i<13;i=i+2){
                    if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
                }
                th0->Cp=th0->Cp*R;
                th0->H=coef[0]*(th0->T- refT)+coef[1]/(coef[2]+1)*(pow(th0->T,(coef[2]+1))-pow(refT,(coef[2]+1)));
                for (i=3;i<13;i=i+2) if (coef[i]>0) th0->H=th0->H+coef[i]*coef[i+1]*(1/(exp(coef[i+1]/th0->T)-1)-1/(exp(coef[i+1]/ refT)-1));
                th0->H=th0->H*R;
                if (coef[1]>0) th0->S=coef[0]*log(th0->T/ refT)+coef[1]/coef[2]*(pow(th0->T,coef[2])-pow(refT,coef[2]));
                else th0->S=coef[0]*log(th0->T/ refT);
                //printf("T:%f S0(0):%f\n",th0->T,th0->S*R/th0->MW);
                for (i=3;i<13;i=i+2){
                    if (coef[i]>0) th0->S=th0->S+coef[i]*coef[i+1]*(exp(coef[i+1]/th0->T)/(th0->T*(exp(coef[i+1]/th0->T)-1))-
                        log(exp(coef[i+1]/th0->T)-1)/coef[i+1]-exp(coef[i+1]/ refT)/(refT*(exp(coef[i+1]/ refT)-1))+log(exp(coef[i+1]/ refT)-1)/coef[i+1]);
                    //printf("S0(%i):%f\n",i,th0->S*R/th0->MW);
                }
                th0->S=th0->S*R;
                break;
            case 8://Jaeschke and Schley equation (9 coefficients). Used by GERG2004
                th0->Cp=1+coef[0];
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow(coef[i+1]/th0->T/sinh(coef[i+1]/th0->T),2)+coef[i+2]*pow(coef[i+3]/th0->T/cosh(coef[i+3]/th0->T),2);
                th0->Cp=th0->Cp*R;
                th0->H=(1+coef[0])*(th0->T- refT);
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->H=th0->H+2*coef[i]*coef[i+1]*(1/(exp(2*coef[i+1]/th0->T)-1)-1/(exp(2*coef[i+1]/ refT)-1))+
                        +2*coef[i+2]*coef[i+3]*(1/(exp(2*coef[i+3]/th0->T)+1)-1/(exp(2*coef[i+3]/ refT)+1));
                th0->H=th0->H*R;
                th0->S=(1+coef[0])*log(th0->T/ refT);
                for (i=1;i<9;i=i+4) if (coef[i]>0) th0->S=th0->S+coef[i]*(coef[i+1]/th0->T/tanh(coef[i+1]/th0->T)-log(sinh(coef[i+1]/th0->T)))-coef[i]*(coef[i+1]/
                        refT/tanh(coef[i+1]/ refT)-log(sinh(coef[i+1]/ refT)))-coef[i+2]*(coef[i+3]/th0->T*tanh(coef[i+3]/th0->T)-log(cosh(coef[i+3]/th0->T)))+
                        coef[i+2]*(coef[i+3]/ refT*tanh(coef[i+3]/ refT)-log(cosh(coef[i+3]/ refT)));
                th0->S=th0->S*R;
                break;
            default:
            th0->Cp=th0->H=th0->S=0;
        }
    }

    switch (equation)//units conversion if necessary
    {
    case 1://DIPPR 100 in KJ/kgr·K
        th0->Cp=th0->Cp*th0->MW;
        th0->H=th0->H*th0->MW;
        th0->S=th0->S*th0->MW;
        break;
    case 6://DIPPR 100 in J/kgr·K
    case 200://DIPPR 107 correlation in J/kg·K
        th0->Cp=th0->Cp*th0->MW*1e-3;
        th0->H=th0->H*th0->MW*1e-3;
        th0->S=th0->S*th0->MW*1e-3;
        break;
    case 3://DIPPR 107 in calories/mol·K
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 in cal/(mol·K)
        th0->Cp=th0->Cp*4.1868;
        th0->H=th0->H*4.1868;
        th0->S=th0->S*4.1868;
        break;
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K
        th0->Cp=th0->Cp*1e-3;
        th0->H=th0->H*1e-3;
        th0->S=th0->S*1e-3;
        break;
    }

    if (th0->Cp>0){
        th0->S=th0->S-R*log(th0->P/ refP);//We change the reference of entropy from the same P to the reference P
        th0->Cv=th0->Cp-R;
        th0->U=th0->H-R*th0->T;
        th0->A=th0->U-th0->T*th0->S;
        th0->G=th0->H-th0->T*th0->S;
    }
    else{
        th0->Cv=th0->U=th0->A=th0->G=0;
    }


    //printf("T:%f MW:%f Cp0:%f  G0:%f  H0:%f  S0:%f\n",th0->T,th0->MW,th0->Cp,th0->G,th0->H,th0->S);
}


//Ideal gas thermodynamic properties of water calculation , from a reference state specified by the triple point where H and S are 0
//----------------------------------------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_IdealThermoWater( FF_ThermoProperties *th0)
{
    double delta=1/(th0->V*1.787371609e4);
    double tau=647.096/th0->T;
    double A0r,dt,dt2;//A0r and derivatives
    //double dd,dd2;
    //double n[8];//={–8.3204464837497,6.6832105275932,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873};

    A0r=(log(delta)-8.3204464837497+6.6832105275932*tau+3.00632*log(tau)+0.012436*log(1-exp(-1.28728967*tau))+0.97315*log(1-exp(-3.53734222*tau))+
            1.27950*log(1-exp(-7.74073708*tau))+0.96956*log(1-exp(-9.24437796*tau))+0.24873*log(1-exp(-27.5075105*tau)));
    //dd=1/delta;
    //dd2=-1/pow(delta,2);
    dt=6.6832105275932+3.00632/tau+0.012436*1.28728967*(pow((1-exp(-1.28728967*tau)),-1)-1)+0.97315*3.53734222*(pow((1-exp(-3.53734222*tau)),-1)-1)+
            1.27950*7.74073708*(pow((1-exp(-7.74073708*tau)),-1)-1)+0.96956*9.24437796*(pow((1-exp(-9.24437796*tau)),-1)-1)+
            0.24873*27.5075105*(pow((1-exp(-27.5075105*tau)),-1)-1);
    dt2=-3.00632/pow(tau,2)-0.012436*pow(1.28728967,2)*exp(-1.28728967*tau)*pow((1-exp(-1.28728967*tau)),-2)-
            0.97315*pow(3.53734222,2)*exp(-3.53734222*tau)*pow((1-exp(-3.53734222*tau)),-2)-
            1.27950*pow(7.74073708,2)*exp(-7.74073708*tau)*pow((1-exp(-7.74073708*tau)),-2)-
            0.96956*pow(9.24437796,2)*exp(-9.24437796*tau)*pow((1-exp(-9.24437796*tau)),-2)-
            0.24873*pow(27.5075105,2)*exp(-27.5075105*tau)*pow((1-exp(-27.5075105*tau)),-2);
    th0->A=R*th0->T*A0r;
    th0->S=R*(tau*dt-A0r);
    th0->U=th0->A+th0->S*th0->T;
    th0->G=th0->A+R*th0->T;
    th0->H=th0->U+R*th0->T;
    th0->Cv=-R*pow(tau,2)*dt2;
    th0->Cp=th0->Cv+R;
}

//Enthalpy and entropy calculation from T,V and P using EOS
//---------------------------------------------------------
void CALLCONV FF_HSfromTVPeosS(double T, double V, double P, const FF_SubstanceData *data, double *H, double *S)
{
    FF_ThermoProperties th0;
    double ArrDer[2];
     FF_CubicParam param;
    double delta,tau;
    bool water=false;
    th0.T= T;
    th0.V= V;

    switch (data->model)
    {
        case FF_SAFTtype:
            FF_ArrDerSAFT0T(T,V,&data->saftData,ArrDer);
            th0.MW=data->saftData.MW;
            break;
        case FF_SWtype:
            delta=1/(V*data->swData.rhoRef);
            tau=data->swData.tRef/ T;
            FF_ArrDerSW0T(tau,delta,&data->swData,ArrDer);
            //printf("HS calc Arr[0]:%f Arr[1]:%f\n",ArrDer[0],ArrDer[1]);
            th0.MW=data->swData.MW;
            if (data->swData.eos==FF_IAPWS95)water=true;
            break;
        default:
            FF_FixedParamCubic(&data->cubicData, &param);
            FF_ThetaDerivCubic(&T,&data->cubicData, &param);
            FF_ArrDerCubic0T(T,V,&param,ArrDer);
            th0.MW=data->cubicData.MW;
            break;
    }
    if (water==true){
        FF_IdealThermoWater(&th0);
    }
    else{
        FF_IdealThermoEos(data->cp0Corr.form,data->cp0Corr.coef,data->refT,data->refP,&th0);
        //printf("refT:%f refP:%f th0.H:%f th0.S:%f\n",data->refT,data->refP,th0.H,th0.S);
    }
    if (data->model==FF_SWtype) *S=R*(tau*ArrDer[1] - ArrDer[0]) + th0.S;
    else *S=-R*(T*ArrDer[1] + ArrDer[0]) + th0.S;
    *H=R* T*ArrDer[0] + P * V - R* T+ T * (*S-th0.S) + th0.H;
    //printf("H:%f\n",*H);
}


//Residual extended thermodynamic properties calculation from T and V, using EOS
//------------------------------------------------------------------------------
void CALLCONV FF_ExtResidualThermoEosS(const FF_SubstanceData *data, FF_ThermoProperties *thR)
{
    double ArrDer[6];
     FF_CubicParam param;
    double delta,tau;
    switch (data->model)
    {
        case FF_SAFTtype:
            FF_ArrDerSAFT(thR->T,thR->V,&data->saftData,ArrDer);
            break;
        case FF_SWtype:
            delta=1/(thR->V*data->swData.rhoRef);
            tau=data->swData.tRef/thR->T;
            //printf("delta:%f  tau:%f\n",delta,tau);
            FF_ArrDerSW(tau,delta,&data->swData,ArrDer);
            break;
        default:
            FF_FixedParamCubic(&data->cubicData, &param);
            FF_ThetaDerivCubic(&thR->T,&data->cubicData, &param);
            FF_ArrDerCubic(thR->T,thR->V,&param,ArrDer);
            break;
    }
    //printf("Arr:%f dArr/dV:%f d2Arr/dV2:%f dArr/dT:%f d2Arr/dT2:%f d2Arr/dVdT:%f\n",ArrDer[0],ArrDer[1],ArrDer[2],ArrDer[3],ArrDer[4],ArrDer[5]);
    if (data->model==FF_SWtype)
    {
        thR->P=R*thR->T*(1+delta*ArrDer[1])/thR->V;
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=R*(tau*ArrDer[3]-ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        thR->H=thR->G+thR->T*thR->S;
        thR->dP_dV=-R*thR->T*(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])/thR->V/thR->V;
        thR->dP_dT=R*(1+delta*ArrDer[1]-delta*tau*ArrDer[5])/thR->V;
        thR->Cv=-R*tau*tau*ArrDer[4];
        thR->Cp=thR->Cv+R*pow((1+delta*ArrDer[1]-delta*tau*ArrDer[5]),2)/(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])-R;
        //thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;
    }
    else
    {
        thR->P=R*thR->T*(1/thR->V-ArrDer[1]);
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=-R*(thR->T*ArrDer[3]+ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        //thR->U=-R*thR->T*thR->T*ArrDer[3];
        //printf("%f %f\n",thR->U,-R*thR->T*thR->T*ArrDer[3]);
        thR->H=thR->G+thR->T*thR->S;
        //thR->H=thR->U-R*thR->T*thR->V*ArrDer[1];
        //thR->H=-R*thR->T*(thR->V*ArrDer[1]+thR->T*ArrDer[3]);
        thR->Cv=-2*R*thR->T*ArrDer[3]-R*thR->T*thR->T*ArrDer[4];
        thR->dP_dV=R*thR->T*(-1/(thR->V*thR->V)-ArrDer[2]);
        thR->dP_dT=R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]);
        //thR->Cp=thR->Cv-R*pow((1/thR->V-ArrDer[1]-thR->T*ArrDer[5]),2)/(-1/thR->V/thR->V-ArrDer[2])-R;
        thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;

    }
}


//Thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0,using FF_SubstanceData
//-------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ThermoEosS(const FF_SubstanceData *data, FF_ThermoProperties *th)
{
     FF_ThermoProperties th0,thR;
    bool water=false;
    switch (data->model)
    {
        case FF_SAFTtype:
            th->MW=data->saftData.MW;
            break;
        case FF_SWtype:
            th->MW=data->swData.MW;
            if (data->swData.eos==FF_IAPWS95)water=true;
            break;
        default://Cubic EOS
            th->MW=data->cubicData.MW;
            break;
    }
    th0.MW=thR.MW=th->MW;//Perhaps not necessary
    th0.T=thR.T=th->T;
    th0.V=thR.V=th->V;
    if (water==true){
        FF_IdealThermoWater(&th0);
    }
    else{
        FF_IdealThermoEos(data->cp0Corr.form,data->cp0Corr.coef,data->refT,data->refP,&th0);
    }
    if (data->model==FF_IdealType)
    {
        th->P=R*th->T/th->V;
        th->A=th0.A;
        th->G=th0.G;
        th->S=th0.S;
        th->U=th0.U;
        th->H=th0.H;
        th->dP_dV=-R*th->T/(th->V*th->V);
        th->dP_dT=R/th->V;
        th->Cv=th0.Cv;
        th->Cp=th0.Cp;
    }
    else
    {
        FF_ExtResidualThermoEosS(data,&thR);
        th->P=thR.P;
        th->A=th0.A+thR.A;
        th->G=th0.G+thR.G;
        th->S=th0.S+thR.S;
        th->U=th0.U+thR.U;
        th->H=th0.H+thR.H;
        th->dP_dV=thR.dP_dV;
        th->dP_dT=thR.dP_dT;
        th->Cv=th0.Cv+thR.Cv;
        th->Cp=th0.Cp+thR.Cp;
    }
    th->SS=th->V*pow(-th->Cp*th->dP_dV/th->Cv/th->MW*1e3,0.5);
    th->JT=-(th->T*th->dP_dT/th->dP_dV+th->V)/th->Cp;//Joule Thomson coefficient =(dT/dP) at constant H
    th->IT=-th->JT*th->Cp;//Isothermal throttling coefficient = (dH/dP) at constant T
    //printf("T:%f V:%f P:%f Cv:%f Cp:%f H:%f U:%f S:%f G:%f A:%f dP_dV:%f dP_dT:%f SS:%f JT:%f IT:%f\n",th->T,th->V,th->P,th->Cv,th->Cp,th->H,th->U,th->S,th->G,th->A,th->dP_dV,th->dP_dT,th->SS,th->JT,th->IT);
    //printf("Ideal T:%f V:%f P:%f Cv0:%f Cp0:%f H0:%f U0:%f S0:%f G:%f A:%f\n",th0.T,th0.V,th0.P,th0.Cv,th0.Cp,th0.H,th0.U,th0.S,th0.G,th0.A);
    //printf("Residual T:%f V:%f P:%f Cvr:%f Cpr:%f Hr:%f Ur:%f Sr:%f Gr:%f Ar:%f\n",thR.T,thR.V,thR.P,thR.Cv,thR.Cp,thR.H,thR.U,thR.S,thR.G,thR.A);
}





//Multiphase state resolution and thermodynamic properties calculation
//====================================================================



//Phases and T,V solver from P and T or H or U or S
void CALLCONV FF_TVsFromPXNewton(char var, FF_SubstanceData *data, int aid, double p, double x, double *T, double *Vg, double *Vl, double *gf){
    int ver=0;//controls the ejecution of the printf statements to provide information
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double Tc,Pc,Vc,Tb=0;
    double Xraw = 0.0,Xh,Xl;
    double answerL[3],answerG[3];
    char option,state;
    double ld,gd;
    int i;
    double error;
    double grad;//T derivative of the state variable at constant pressure
    FF_ThermoProperties th;//where to store a phase properties
    if ((var=='h')||(var=='H')) Xraw=x + data->refH;
    else if((var=='s')||(var=='S')) Xraw=x + data->refS;
    else if((var=='u')||(var=='U')) Xraw=x;
    th.P=p;
    if (ver==1) printf("Pressure:%f, looking for %c, Objective value, original:%f raw:%f\n",p,var,x,Xraw);

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
        Vc=R*data->cubicData.Zc*data->cubicData.Tc/data->cubicData.Pc;
    }

    if ((var=='t')||(var=='T')){
        if((p>0.9999*Pc)&&(p<1.0001*Pc)&&((data->model==FF_SWtype)||(data->model==FF_CubicType))){//SAFT data may contain, or not, the true critical constants
            *Vg=Vc;
            *Vl=0;
            *gf=1;
        }
        else{
            option='s';
            FF_VfromTPeosS(x,p,data,option,answerL,answerG,&state);
            if((state=='l')||(state=='L')){
                *Vl=answerL[0];
                *Vg=0;
                *gf=0;
            }
            else if((state=='g')||(state=='G')){
                *Vl=0;
                *Vg=answerG[0];
                *gf=1;
            }
        }
        error=0;
    }
    else{//we need to find T
        //Stablish boiling point
        if(p>0.9999*Pc) Tb=Tc;
        else if ((((data->model==FF_SWtype)&&(aid==3))||((data->model==FF_SAFTtype)&&(aid==2)))&&data->vpCorr.form>0){//Bajo Pc hay que buscar primero el punto de ebullición
            double pResult;//The returned pressure, to check
            FF_CorrelationSolver(p,&data->vpCorr,data->swData.MW,&Tb,&pResult);
            if (ver==1) printf("Tb found for SW by correlation:%f \n",Tb);
        }
        else FF_TbEosS(p,data,&Tb);//If pressure bellow Pc we find the Tb
        if (ver==1) printf("Tb asigned:%f\n",Tb);

        //Calculate saturated volumes/densities
        if((Tb==Tc)&&(p>0.9999*Pc)&&(p<1.0001*Pc)&&((data->model==FF_SWtype)||(data->model==FF_CubicType))) *Vl=*Vg=Vc;
        else if ((data->model==FF_SWtype)&&(p<0.9999*Pc)&&(data->lDensCorr.form>0)&&(data->gDensCorr.form>0)&&(aid==3)){
            FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,Tb,&ld);
            FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,Tb,&gd);
            *Vl=MW/ ld;
            *Vg=MW/ gd;
            if (ver==1) printf("Saturated densities found for SW by correlations ld:%f gd:%f \n",ld,gd);
        }
        else{
            option='b';
            FF_VfromTPeosS(Tb,p,data,option,answerL,answerG,&state);
            if(answerL[0]>0) *Vl=*Vg=answerL[0];
            if(answerG[0]>0){
                if(*Vl==0) *Vl=*Vg=answerG[0];
                else if(!(answerG[0]==*Vl)) *Vg=answerG[0];
             }
        }
        if (ver==1) printf("Initial P:%f Tb:%f volumes found Vl:%f Vg:%f state:%c\n",p,Tb,*Vl,*Vg,state);
        //determination of liquid and gas values at boiling point
        th.T=Tb;
        if (*Vl>0){
            th.V=*Vl;
            FF_ThermoEosS(data,&th);
            //FF_HSfromTVPeosS(Tb,Vl,p,data,&Hl,&Sl);
            if((var=='h')||(var=='H')) Xh=Xl=th.H;
            else if((var=='s')||(var=='S')) Xh=Xl=th.S;
            else if((var=='u')||(var=='U')) Xh=Xl=th.U;

        }
        if ((*Vg>0)&&!(*Vg==*Vl)){
            th.V=*Vg;
            FF_ThermoEosS(data,&th);
            //FF_HSfromTVPeosS(Tb,Vg,p,data,&Hh,&Sh);
            if(*Vl==0){
                if((var=='h')||(var=='H')) Xh=Xl=th.H;
                else if((var=='s')||(var=='S')) Xh=Xl=th.S;
                else if((var=='u')||(var=='U')) Xh=Xl=th.U;
            }
            else{
                if((var=='h')||(var=='H')) Xh=th.H;
                else if((var=='s')||(var=='S')) Xh=th.S;
                else if((var=='u')||(var=='U')) Xh=th.U;
            }
        }
        if (ver==1) printf("Initial T:%f P:%f Vg:%f Vl:%f Xraw:%f Xh:%f Xl:%f\n",th.T,th.P,*Vg,*Vl,Xraw,Xh,Xl);


       if((Xraw>Xl)&&(Xraw<Xh)){
           //two phases situation. We already have the solution for T and D
           *T=Tb;
           //*Vl=answerL[0];
           //*Vg=answerG[0];
           *gf=(Xraw-Xl)/(Xh-Xl);
           error=0;
           if (ver==1)printf("Two phases\n");
       }
       else if (Xraw==Xl){
           *T=Tb;
           //*Vl=answerL[0];
           *Vg=0;
           *gf=0;
           error=0;
       }
       else if (Xraw==Xh){
           *T=Tb;
           *Vl=0;
           //*Vg=answerG[0];
           *gf=1;
           error=0;
       }
       else{//Begin search of T and d
           if (Xraw>Xh){
               //option='g';
               error=Xraw-Xh;
           }
           else{
               //option='l';
               error=Xraw-Xl;
           }
           th.T=Tb;
           if((var=='h')||(var=='H')) grad=th.Cp;
           else if((var=='s')||(var=='S')) grad=th.Cp/th.T;
           else if((var=='u')||(var=='U')) grad=th.Cp+th.P*th.dP_dT/th.dP_dV;

           i=1;
           while (((fabs(error/Xraw)>0.000001))&&(i<10)){
               if (ver==1) printf("i:%i T:%f P:%f V%f Xraw:%f error:%f grad:%f \n",i,th.T,th.P, th.V, Xraw,error,grad);
               *T=th.T+error/grad;
               FF_VfromTPeosS(*T,p,data,option,answerL,answerG,&state);
               th.T=*T;
               if (Xraw<Xl) th.V=answerL[0];
               else th.V=answerG[0];
               FF_ThermoEosS(data,&th);
               i++;
               if((var=='h')||(var=='H')){
                   error=Xraw-th.H;
                   grad=th.Cp;
               }
               else if((var=='s')||(var=='S')){
                   error=Xraw-th.S;
                   grad=th.Cp/th.T;
               }
               else if((var=='u')||(var=='U')){
                   error=Xraw-th.U;
                   grad=th.Cp+th.P*th.dP_dT/th.dP_dV;
               }
           }
           if (Xraw<Xl){
               *Vl=answerL[0];
               *Vg=0;
               *gf=0;
           }
           else{
               *Vl=0;
               *Vg=answerG[0];
               *gf=1;
           }

       }
    }
   if (ver==1) printf("T:%f error:%f Vl:%f Vg:%f gf:%f\n",*T, error,*Vl,*Vg,*gf);
}

//Phases and T,D solver from P and H or S. Old function
void CALLCONV FF_TVfromPX(char var, FF_SubstanceData *data, double p, double x, double *T, double *gd, double *ld){
    int searchDir=0;//1= searh up, -1=search interval
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double Tc,Pc,Vc,Tb=0,Th,Tl,Tn;
    double Xraw = 0.0,Hl,Hh,Sl,Sh,Xh,Xl,Cp0,Hn,Sn,Xn;
    double answerL[3],answerG[3];
    char option,state;
    double V=0,Vl=0,Vg=0;
    int i;
    FF_ThermoProperties th;//where to store a phase properties
    if ((var=='h')||(var=='H')) Xraw=x*MW + data->refH;//enthalpy in J/mol
    else if((var=='s')||(var=='S')) Xraw=x*MW + data->refS;//Entropy in In J/(mol·K)
    th.P=p;

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
        Vc=R*data->cubicData.Zc*data->cubicData.Tc/data->cubicData.Pc;
    }

    //Stablish boiling point
    if(p>0.9999*Pc) Tb=Tc;
    else if ((data->model==FF_SWtype)&&(data->vpCorr.form>0)){//Bajo Pc hay que buscar primero el punto de ebullición
        double pResult;//The returned pressure, to check
        FF_CorrelationSolver(p,&data->vpCorr,data->swData.MW,&Tb,&pResult);
    }
    else FF_TbEosS(p,data,&Tb);//If pressure bellow Pc we find the Tb
    //printf("Tb asigned:%f\n",Tb);

    //Calculate saturated volumes/densities
    if((Tb==Tc)&&(p>0.9999*Pc)&&(p<1.0001*Pc)&&((data->model==FF_SWtype)||(data->model==FF_CubicType))) Vl=Vg=Vc;
    else if((data->model==FF_SWtype)&&(p<0.9999*Pc)&&(data->lDensCorr.form>0)&&(data->gDensCorr.form>0)){
        FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,Tb,ld);
        FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,Tb,gd);
        Vl=MW/ *ld;
        Vg=MW/ *gd;
    }
    else{
        option='b';
        FF_VfromTPeosS(Tb,p,data,option,answerL,answerG,&state);
        if(answerL[0]>0) Vl=Vg=answerL[0];
        if(answerG[0]>0){
            if(Vl==0) Vl=Vg=answerG[0];
            else if(!(answerG[0]==Vl)) Vg=answerG[0];
         }
    }
    //printf("Initial P:%f Tb:%f volumes found Vl:%f Vg:%f state:%c\n",p,Tb,Vl,Vg,state);
    //determination of liquid and gas enthalpies at boiling point
    if (Vl>0){
        FF_HSfromTVPeosS(Tb,Vl,p,data,&Hl,&Sl);
        if((var=='h')||(var=='H')) Xh=Xl=Hl;
        else if((var=='s')||(var=='S')) Xh=Xl=Sl;

    }
    if ((Vg>0)&&!(Vg==Vl)){
        FF_HSfromTVPeosS(Tb,Vg,p,data,&Hh,&Sh);
        if(Vl==0){
            if((var=='h')||(var=='H')) Xh=Xl=Hh;
            else if((var=='s')||(var=='S')) Xh=Xl=Sh;
        }
        else{
            if((var=='h')||(var=='H')) Xh=Hh;
            else if((var=='s')||(var=='S')) Xh=Sh;
        }
    }
    //printf("Initial T:%f Xraw:%f Xh:%f Xl:%f\n",Tb,Xraw,Xh,Xl);

    //Begin search of T and d
    if (Xraw>=Xh){
        Tl=Th=Tb;
        Xl=Xh;
        searchDir=1;//search upwards, no limit. Gas
        //printf("Gas, search:%i\n",searchDir);
    }
    else if (Xraw<=Xl){
        Th=Tb;
        Xh=Xl;
        if (data->lCpCorr.limI>0) Tl=data->lCpCorr.limI;
        else if (data->baseProp.Tm>0) Tl=data->baseProp.Tm;
        else if (data->cp0Corr.limI>0) Tl=data->cp0Corr.limI;
        else Tl=273.15;
        option='l';
        FF_VfromTPeosS(Tl,p,data,option,answerL,answerG,&state);
        FF_HSfromTVPeosS(Tl,answerL[0],p,data,&Hl,&Sl);
        if((var=='h')||(var=='H')) Xl=Hl;
        else if((var=='s')||(var=='S')) Xl=Sl;
        searchDir=-1;//search interval. Liquid
        //printf("Liquid, search:%i\n",searchDir);
    }
    else if((Xraw>Xl)&&(Xraw<Xh)){//two phases situation. We already have the solution for T and D
        *T=Tb;
        *ld=MW/answerL[0];
        *gd=MW/answerG[0];
        //printf("Two phases\n");
    }

    if (searchDir==1){//search up. Only SAFT could be liquid for some degrees, and vaporize later. The others are always gas
        *ld=0.0;
        option='b';
        FF_PhysPropCorrM(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,Th,&Cp0);
        if((var=='h')||(var=='H')) Tn=Th+(Xraw-Xh)/Cp0;
        else if((var=='s')||(var=='S')) Tn=Th+(Xraw-Xh)*Th/Cp0;
        FF_VfromTPeosS(Tn,p,data,option,answerL,answerG,&state);
        if (state=='L') V=answerL[0];
        else V=answerG[0];
        FF_HSfromTVPeosS(Tn,V,p,data,&Hn,&Sn);
        Th=Tn;
        if((var=='h')||(var=='H')) Xh=Xn=Hn;
        else if((var=='s')||(var=='S')) Xh=Xn=Sn;
        i=1;
        //printf("Initial search Cp0:%f Xraw:%f Tl:%f Hl:%f, Th:%f Hh:%f\n",Cp0,Xraw,Tl,Hl,Th,Hh);
        while (fabs((Xraw-Xn)/Xraw)>0.00001){//Regula Falsi search, try to improve to Illinois or Anderson
            Tn=Tn+(Xraw-Xn)/(Xh-Xl)*(Th-Tl);
            FF_VfromTPeosS(Tn,p,data,option,answerL,answerG,&state);
            if (state=='L') V=answerL[0];
            else V=answerG[0];
            FF_HSfromTVPeosS(Tn,V,p,data,&Hn,&Sn);
            if((var=='h')||(var=='H')) Xn=Hn;
            else if((var=='s')||(var=='S')) Xn=Sn;
            if(Xn>Xraw){
                Th=Tn;
                Xh=Xn;
            }
            else{
                Tl=Tn;
                Xl=Xn;
            }
            //printf("Seeking gas temperature i:%i Tn:%f Hn:%f Xraw:%f, Xl:%f Xh:%f\n",i,Tn,Hn,Xraw,Xl,Xh);
            i++;
            if (i>20) break;
        }
        *T=Tn;
        *gd=MW/V;
    }

    else if (searchDir==-1){//search interval. Liquid phase
        option='l';
        i=1;
        Tn=Tl+(Xraw-Xl)/(Xh-Xl)*(Th-Tl);
        FF_VfromTPeosS(Tn,p,data,option,answerL,answerG,&state);
        if (answerL[0]>0) V=answerL[0];
        else V=answerG[0];
        FF_HSfromTVPeosS(Tn,V,p,data,&Hn,&Sn);
        if((var=='h')||(var=='H')) Xn=Hn;
        else if((var=='s')||(var=='S')) Xn=Sn;
        //printf("Initial interval search Tl:%f Xl:%f, Th:%f Xh:%f Tn:%f Xn:%f\n",Tl,Xl,Th,Xh,Tn,Xn);
        while (fabs((Xraw-Xn)/Xraw)>0.00001){
            if (Xraw>Xn){
                Xl=Xn;
                Tl=Tn;
            }
            else{
                Xh=Xn;
                Th=Tn;
            }
            Tn=Tl+(Xraw-Xl)/(Xh-Xl)*(Th-Tl);
            FF_VfromTPeosS(Tn,p,data,option,answerL,answerG,&state);
            if (answerL[0]>0) V=answerL[0];
            else V=answerG[0];
            FF_HSfromTVPeosS(Tn,V,p,data,&Hn,&Sn);
            if((var=='h')||(var=='H')) Xn=Hn;
            else if((var=='s')||(var=='S')) Xn=Sn;
            //printf("Seeking interval i:%i Tn:%f Xn:%f Xraw:%f, Xl:%f\n",i,Tn,Xn,Xraw,Xl);
            i++;
            if (i>20) break;
        }
        *T=Tn;
        if((state=='l')||(state=='L')){
            *ld=MW/V;
            *gd=0;
        }
        else{
            *gd=MW/V;
            *ld=0;
        }
    }
    //printf("T:%f ld:%f gd:%f\n",*T,*ld,*gd);
}

//Phases and T,V solver from V and T
void CALLCONV FF_TVsFromVX(char var, FF_SubstanceData *data, int aid, double v, double x, double *T, double *Vg, double *Vl, double *gf){
    int ver=0;//controls the ejecution of the printf statements to provide information
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    if((var=='t')||(var=='T')){//T and overall D are fixed
        double Tc,V,dh,dl;
        double Vp;
        double answerL[3],answerG[3];
        char option,state;

        *T=x;
        if (data->model==FF_SWtype) Tc=data->swData.Tc;
        else if(data->model==FF_SAFTtype) Tc=data->saftData.Tc;
        else Tc=data->cubicData.Tc;

        if (x>=Tc){//over Tc only gas phase exists
            *Vg=v;
            *Vl=0.0;
            *gf=1;
        }
        else{
            int end=0;
            if((((data->model==FF_SWtype)&&(aid==3))||((data->model==FF_SAFTtype)&&(aid==2)))&&data->vpCorr.form>0) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,x,&Vp);
            else FF_VpEosS(x,data,&Vp);
            FF_VpEosS(x,data,&Vp);

            if((((data->model==FF_SWtype)&&(aid==3))||((data->model==FF_SAFTtype)&&(aid==2)))&&(data->lDensCorr.form>0)&&(data->gDensCorr.form>0)){
                FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,x,&dh);
                FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,x,&dl);
                if(v<=0.95*MW/dh){
                    *Vl=v;
                    *Vg=0;
                    *gf=0;
                    end=1;
                }
                else if(v>=1.05*MW/dl){
                    *Vg=v;
                    *Vl=0.0;
                    *gf=1;
                    end=1;
                }
            }
            if (end==0){//for cubic and SAFT type eos, and not clearly monophasic SW eos
                //if((data->model==FF_SWtype)&&(data->vpCorr.form>0)) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,x,&Vp);
                //else FF_VpEosS(x,data,&Vp);
                option='b';
                FF_VfromTPeosS(x,Vp,data,option,answerL,answerG,&state);//calculation of liquid and gas volumes at saturation
                if(v<=answerL[0]){
                    *Vl=v;
                    *Vg=0;
                    *gf=0;
                }
                else if (v>=answerG[0]){
                    *Vg=v;
                    *Vl=0.0;
                    *gf=1;
                }
                else{
                     *Vl=answerL[0];
                     *Vg=answerG[0];
                     *gf=(v-*Vl)/(*Vg-*Vl);
                }
            }
        }
        if (ver==1) printf("Two phases T:%f V:%f d:%f Vp;%f Vl:%f Vg:%f Dl:%f Dg:%f gf:%f\n",x,v,MW/v,Vp,*Vl, *Vg, MW/ *Vl,MW/ *Vg,*gf);
    }
}




//all thermo properties from T/P or T/D or P/H or P/S
void CALLCONV FF_solveEos(char *variable, FF_SubstanceData *data, int aid, double x, double y, double *T, double *p, double *gd, double *gh, double *gs, double *gCv, double *gCp, double *gDvp,
                          double *gDvT, double *ld, double *lh, double *ls, double *lCv, double *lCp, double *lDvp, double *lDvT, double *gf){
    int ver=0;
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double nMols=1/MW;//number of moles in a kg
    double Tc,Pc,Vc;
    double Vl,Vg;
    char var;
    var=variable[0];
    FF_ThermoProperties th;

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
    if((var=='e')||(var=='E')){//Bubble state
        if (x>=Tc){
           *T=Tc;//T is fixed
           *p=Pc;//P is fixed
           *ld=MW/Vc;
        }
        else{
            *T=x;//T is fixed
            *p=y;//P is fixed
            if((((data->model==FF_SWtype)&&(aid==3))||((data->model==FF_SAFTtype)&&(aid==2)))&&data->lDensCorr.form>0) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,x,ld);
            else{
                double answerL[3],answerG[3];
                char option,state;
                option='l';
                FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
                *ld=MW/answerL[0];
                }
            *gd=0.0;
        }
        *gd=0.0;
        *gf=0;
    }
    if((var=='w')||(var=='W')){//Dew state
        if (x>=Tc){
           *T=Tc;//T is fixed
           *p=Pc;//P is fixed
           *gd=MW/Vc;
        }
        else{
            *T=x;//T is fixed
            *p=y;//P is fixed
            if ((((data->model==FF_SWtype)&&(aid==3))||((data->model==FF_SAFTtype)&&(aid==2)))&&data->gDensCorr.form>0) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,x,gd);
            else{
                double answerL[3],answerG[3];
                char option,state;
                option='g';
                FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
                *gd=MW/answerG[0];
            }
            *ld=0;
        }
        *ld=0;
        *gf=1;
    }
    else if((var=='p')||(var=='P')){//T and P fixed
        double answerL[3],answerG[3];
        char option,state;
        *T=x;//T is fixed
        *p=y;//P is fixed
        option='s';
        FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
        if((state=='l')||(state=='L')){
            *ld=MW/answerL[0];
            *gd=0.0;
            *gf=0;
        }
        else if((state=='g')||(state=='G')){
            *gd=MW/answerG[0];
            *ld=0.0;
            *gf=1;
        }
    }
    else if((var=='d')||(var=='D')){
        *T=x;//T is fixed
        if (ver==1) printf("D T:%f V:%f D:%f \n",x,MW/y,y);
        FF_TVsFromVX('T',data,aid,MW/y,x,T,&Vg,&Vl,gf);

        if (*gf>0.0) *gd=MW/Vg;
        else *gd=0;
        if (*gf<1.0) *ld=MW/Vl;
        else *ld=0;

    }
    else if((var=='h')||(var=='H')||(var=='s')||(var=='S')||((var=='u')||(var=='U'))){
        *p=x;
        y=y*MW;//pas to molar quantities
        FF_TVsFromPXNewton(var,data,aid,x,y,T,&Vg,&Vl,gf);
        if (*gf>0.0) *gd=MW/Vg;
        else *gd=0;
        if (*gf<1.0) *ld=MW/Vl;
        else *ld=0;

    }
    //thermodynamic calculation once solved Vg,Vl and gf
    th.T=*T;
    if (*gf==1.0){
        *lh= *ls= *lCv= *lCp= *lDvp= *lDvT= 0;
    }
    else{
        th.V=MW/ *ld;
        FF_ThermoEosS(data,&th);
        *lh=(th.H-data->refH)*nMols;
        *ls=(th.S-data->refS)*nMols;
        *lCv=th.Cv*nMols;
        *lCp=th.Cp*nMols;
        *lDvp=nMols/th.dP_dV;
        *lDvT=-nMols*th.dP_dT/th.dP_dV;
        if((var=='d')||(var=='D')) *p=th.P;
    }
    if (*gf==0.0){
        *gh= *gs= *gCv= *gCp= *gDvp= *gDvT= 0;
    }
    else{
        th.V=MW/ *gd;
        FF_ThermoEosS(data,&th);
        *gh=(th.H-data->refH)*nMols;
        *gs=(th.S-data->refS)*nMols;
        *gCv=th.Cv*nMols;
        *gCp=th.Cp*nMols;
        *gDvp=nMols/th.dP_dV;
        *gDvT=-nMols*th.dP_dT/th.dP_dV;
        if((var=='d')||(var=='D')) *p=th.P;
    }
}


//all thermo properties from T/P or T/D or P/H or P/S
void CALLCONV FF_solveEosOld(char *variable, FF_SubstanceData *data, double x, double y, double *T, double *p, double *gd, double *gh, double *gs, double *gCv, double *gCp, double *gDvp,
                             double *gDvT, double *ld, double *lh, double *ls, double *lCv, double *lCp, double *lDvp, double *lDvT){
    double MW=data->baseProp.MW*1e-3;//Mol weight in kg
    double nMols=1/MW;//number of moles in a kg
    char var;
    var=variable[0];
    FF_ThermoProperties th;
    if((var=='e')||(var=='E')){//Bubble state
        *T=x;//T is fixed
        *p=y;//P is fixed
        if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)) FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,x,ld);
        else{
            double answerL[3],answerG[3];
            char option,state;
            option='l';
            FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
            *ld=MW/answerL[0];
            }
        *gd=0.0;
    }
    if((var=='w')||(var=='W')){//Dew state  
        *T=x;//T is fixed
        *p=y;//P is fixed
        if ((data->model==FF_SWtype)&&(data->gDensCorr.form>0)) FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,x,gd);
        else{
            double answerL[3],answerG[3];
            char option,state;
            option='g';
            FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
            *gd=MW/answerG[0];
        }
        *ld=0;

    }
    else if((var=='p')||(var=='P')){//T and P fixed
        double answerL[3],answerG[3];
        char option,state;
        *T=x;//T is fixed
        *p=y;//P is fixed
        option='s';
        FF_VfromTPeosS(x,y,data,option,answerL,answerG,&state);
        if((state=='l')||(state=='L')){
            *ld=MW/answerL[0];
            *gd=0.0;
        }
        else if((state=='g')||(state=='G')){
            *gd=MW/answerG[0];
            *ld=0.0;
        }
    }
    else if((var=='d')||(var=='D')){//T and overall D fixed
        double Tc,V,dh,dl;
        double Vp;
        double answerL[3],answerG[3];
        char option,state;
        if (data->model==FF_SWtype) Tc=data->swData.Tc;
        else if(data->model==FF_SAFTtype) Tc=data->saftData.Tc;
        else Tc=data->cubicData.Tc;
        *T=x;//T is fixed
        V=MW/y;
        if (x>=Tc){//over Tc only gas phase exists
            *gd=y;
            *ld=0.0;
        }
        else{
            if((data->model==FF_SWtype)&&(data->vpCorr.form>0)) FF_PhysPropCorrM(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,x,&Vp);
            else FF_VpEosS(x,data,&Vp);
            FF_VpEosS(x,data,&Vp);
            int end=0;
            if((data->model==FF_SWtype)&&(data->lDensCorr.form>0)&&(data->gDensCorr.form>0)){
                FF_PhysPropCorrM(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,x,&dh);
                FF_PhysPropCorrM(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,x,&dl);
                if(y>=1.05*dh){
                    *ld=y;
                    *gd=0.0;
                    end=1;
                }
                else if(y<=0.95*dl){
                    *gd=y;
                    *ld=0.0;
                    end=1;
                }
            }
            if (end==0){//for cubic and SAFT type eos, and not clearly monophasic SW
                option='b';
                FF_VfromTPeosS(x,Vp,data,option,answerL,answerG,&state);//calculation of liquid and gas volumes at saturation
                if(V<=answerL[0]){
                    *ld=y;
                    *gd=0.0;
                }
                else if (V>=answerG[0]){
                    *gd=y;
                    *ld=0.0;
                }
                else{
                     *ld=MW/answerL[0];
                     *gd=MW/answerG[0];

                    //printf("Vp:%f \n",*p);
                }
            }
        }
        /*if ((*ld==0)||(*gd==0)) FF_PfromTVeosS(x,V,data,p);
        else *p=Vp;*/
        //printf("P obtenida:%f\n",*p);
    }
    else if((var=='h')||(var=='H')||(var=='s')||(var=='S')){
        *p=x;
        FF_TVfromPX(var,data,x,y,T,gd,ld);
    }
    th.T=*T;
    if (*ld==0.0){
        *lh= *ls= *lCv= *lCp= *lDvp= *lDvT= 0;
    }
    else{
        th.V=MW/ *ld;
        FF_ThermoEosS(data,&th);
        *lh=(th.H-data->refH)*nMols;
        *ls=(th.S-data->refS)*nMols;
        *lCv=th.Cv*nMols;
        *lCp=th.Cp*nMols;
        *lDvp=nMols/th.dP_dV;
        *lDvT=-nMols*th.dP_dT/th.dP_dV;
        if((var=='d')||(var=='D')) *p=th.P;
    }
    if (*gd==0.0){
        *gh= *gs= *gCv= *gCp= *gDvp= *gDvT= 0;
    }
    else{
        th.V=MW/ *gd;
        FF_ThermoEosS(data,&th);
        *gh=(th.H-data->refH)*nMols;
        *gs=(th.S-data->refS)*nMols;
        *gCv=th.Cv*nMols;
        *gCp=th.Cp*nMols;
        *gDvp=nMols/th.dP_dV;
        *gDvT=-nMols*th.dP_dT/th.dP_dV;
        if((var=='d')||(var=='D')) *p=th.P;
    }
}



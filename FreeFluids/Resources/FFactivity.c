/*
 * FFactivity.c
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

// contains mainly calculations for activity coefficients
//=======================================================

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFactivity.h"

//Calculates activity coefficients according to Wilson (modified) equation, at given T and composition
void CALLCONV FF_ActivityWilson(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]){
    int i,j;//the loop variables
    double V[*numSubs],Zc;//will contain the saturated liquid volume for each substance
    double Ru=0;//The gas constant in the necessary units
    double RuTinv;//1/(Ru*T) to speed calculation
    double Lambda[*numSubs][*numSubs];
    double LambdaX[*numSubs];
    if (*form!=FF_Ant3) for(i=0;i<*numSubs;i++){
        if(baseProp[i].Zra>0) V[i]=R*baseProp[i].Tc/baseProp[i].Pc*pow(baseProp[i].Zra,(1+pow(1-*T/baseProp[i].Tc,0.2857)));
        else if(baseProp[i].Zc>0) V[i]=R*baseProp[i].Tc/baseProp[i].Pc*pow(baseProp[i].Zc,(1+pow(1-*T/baseProp[i].Tc,0.2857)));
        else if(baseProp[i].Vc>0){
            Zc=baseProp[i].Pc*baseProp[i].Vc/(R*baseProp[i].Tc);
            V[i]=R*baseProp[i].Tc/baseProp[i].Pc*pow(baseProp[i].Zc,(1+pow(1-*T/baseProp[i].Tc,0.2857)));
        }
        else if(baseProp[i].Vliq>0) V[i]=baseProp[i].Vliq;
        else return;
    }
    //printf("V:%f %f\n",V[0],V[1]);
    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    for (i=0;i<*numSubs;i++){
        for(j=0;j<*numSubs;j++){
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C:
                Lambda[i][j]=V[j]/V[i]*exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+
                            pintParam[i][j][2]* *T * *T)*RuTinv);
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C:
                Lambda[i][j]=V[j]/V[i]*exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+
                            pintParam[i][j][2]/ (*T * *T))*RuTinv);
                break;
            case FF_Ant3:
                Lambda[i][j]=exp(pintParam[i][j][0]+pintParam[i][j][1]/ *T+pintParam[i][j][2]*log(*T)+
                             pintParam[i][j][3]* *T+ pintParam[i][j][4]/ (*T * *T));//extended Antoine form
                break;
            }
            //printf("Lambda:%f\n",Lambda[i][j]);
        }
    }
    for (i=0;i<*numSubs;i++){
        LambdaX[i]=0;
        for(j=0;j<*numSubs;j++){
            LambdaX[i]=LambdaX[i]+Lambda[i][j]*x[j];
        }
    }
    for (i=0;i<*numSubs;i++){
        actData[i].lnGammaC=0;
        actData[i].lnGammaSG=0;
        actData[i].lnGammaR=1-log(LambdaX[i]);
        for(j=0;j<*numSubs;j++){
            actData[i].lnGammaR=actData[i].lnGammaR-x[j]*Lambda[j][i]/LambdaX[j];
        }
    }
}

//Calculates activity coefficients according to NRTL equation, at given T and composition
//the first four parameters are for the calculation of tau, the 2 last ones for the calculation of alpha
void CALLCONV FF_ActivityNRTL(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],const enum FF_IntParamForm *form,const double *T,const double x[],
                              FF_SubsActivityData actData[]){
    double Ru=0;//The gas constant in different units
    double RuTinv;//1/(Ru*T) to speed calculation
    double tau[*numSubs][*numSubs];
    double alpha[*numSubs][*numSubs];
    double G[*numSubs][*numSubs];
    double xG[*numSubs],xGtau[*numSubs];
    int i,j;//the loop variables
    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    for (i=0;i<*numSubs;i++){
        for(j=0;j<*numSubs;j++){//we switch between the different calculations of tau, and alpha
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C://polynomial 1 form: a+b*T+c*T^2 for the interaction energy increment tau
                tau[i][j]=(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]* *T * *T)*RuTinv;//Interaction energy differential
                alpha[i][j]=pintParam[i][j][4];
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C://polynomial 2 form: a+b*T+c/T^2
                tau[i][j]=(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]/ (*T* *T))*RuTinv;//polynomial 2 form
                alpha[i][j]=pintParam[i][j][4];
                break;
            case FF_Ant1://Antoine form: a +b/T+c*ln(T)+d*T
                tau[i][j]=pintParam[i][j][0]+pintParam[i][j][1]/ *T+pintParam[i][j][2]*log(*T)+
                                                        pintParam[i][j][3]* *T;//extended Antoine form
                alpha[i][j]=pintParam[i][j][4]+pintParam[i][j][5]* (*T-273.15);
                break;
            }
            G[i][j]=exp(-alpha[i][j]*tau[i][j]);
            //printf("alpha:%f tau:%f\n",alpha[i][j],tau[i][j]);
        }
    }
    for (i=0;i<*numSubs;i++){
        xG[i]=0;
        xGtau[i]=0;
        for(j=0;j<*numSubs;j++){
            xG[i]=xG[i]+x[j]*G[j][i];
            xGtau[i]=xGtau[i]+x[j]*G[j][i]*tau[j][i];
        }
    }
    for (i=0;i<*numSubs;i++){
        actData[i].lnGammaC=0;
        actData[i].lnGammaSG=0;
        actData[i].lnGammaR=xGtau[i]/xG[i];
        for(j=0;j<*numSubs;j++){
            actData[i].lnGammaR=actData[i].lnGammaR+x[j]*G[i][j]/xG[j]*(tau[i][j]-xGtau[j]/xG[j]);
        }
    }


}

//Calculates activity coefficients according to UNIQUAQ equation, at given T and composition
void CALLCONV FF_ActivityUNIQUAC(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],
                                 const enum FF_IntParamForm *form,const double *T,const double x[],FF_SubsActivityData actData[]){
    double Ru=0;//The gas constant in different units
    double RuTinv;//1/(Ru*T) to speed calculation
    double R=0,Q=0,Qres=0;//mean molecular volume and surface
    double Phi[*numSubs],Theta[*numSubs],ThetaRes[*numSubs];//the volume and surface fractions of all individual substances
    int i,j;//the loop variables

    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
    case FF_Pol3K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
    case FF_Pol3J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
    case FF_Pol3C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    //combinatorial part calculation
    for (i=0;i<*numSubs;i++){
        R=R+x[i]*baseProp[i].r;
        Q=Q+x[i]*baseProp[i].q;
        Qres=Qres+x[i]*baseProp[i].qRes;
    }
    for (i=0;i<*numSubs;i++){
        Phi[i]=baseProp[i].r/R;
        Theta[i]=baseProp[i].q/Q;
        actData[i].lnGammaC=1-Phi[i]+log(Phi[i]);
        actData[i].lnGammaSG=-5*baseProp[i].q*(1-Phi[i]/Theta[i]+log(Phi[i]/Theta[i]));
        //printf("GammaC: %f\n",exp(actData[i].lnGammaC));
    }

    //residual part
    double tau[*numSubs][*numSubs], S[*numSubs];

    for (i=0;i<*numSubs;i++){
        ThetaRes[i]=x[i]*baseProp[i].qRes/Qres;
        for(j=0;j<*numSubs;j++){
            switch(*form){
            case FF_Pol1K://polynomial 1 form: a+b*T+c*T^2
            case FF_Pol1J:
            case FF_Pol1C:
                tau[i][j]=exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]* *T * *T)*RuTinv);
                break;
            case FF_Pol2K://polynomial 2 form: a+b*T+c/T^2
            case FF_Pol2J:
            case FF_Pol2C:
                tau[i][j]=exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]/ (*T* *T))*RuTinv);
                break;
            case FF_Pol3K://polynomial 3 form: a+b/T+c*T
            case FF_Pol3J:
            case FF_Pol3C:
                tau[i][j]=exp(-(pintParam[i][j][0]+pintParam[i][j][1]/ *T+pintParam[i][j][2]* *T)*RuTinv);
                break;
            case FF_Ant3:
                tau[i][j]=exp(pintParam[i][j][0]+pintParam[i][j][1]/ *T+pintParam[i][j][2]*log(*T)+
                                    pintParam[i][j][3]* *T+pintParam[i][j][4]/ (*T * *T));//extended Antoine form
                break;
            }
        }
    }
    for (i=0;i<*numSubs;i++){
        S[i]=0;
        for(j=0;j<*numSubs;j++){
            S[i]=S[i]+ThetaRes[j]*tau[j][i];
            //S[i]=S[i]+Theta[j]*tau[j][i];
        }
    }
    for (i=0;i<*numSubs;i++){
        actData[i].lnGammaR=1-log(S[i]);
        for(j=0;j<*numSubs;j++){
            actData[i].lnGammaR=actData[i].lnGammaR-tau[i][j]*ThetaRes[j]/S[j];
            //lnGammaR[i]=lnGammaR[i]-tau[i][j]*Theta[j]/S[j];
        }
        actData[i].lnGammaR=baseProp[i].q*actData[i].lnGammaR;
        //printf("GammaR: %f\n",exp(actData[i].lnGammaR));
    }
}

//Calculates activity coefficients in polymer systems, using FV for combinatorial and segment UNIQUAQ for residual, at given T and composition
//For polymers x is that of the polymer. FV is for the monomeric unit. In baseProp.numMono we supply the number of monomeric units of the chain, must be =1 for solvents
void CALLCONV FF_ActivityUNIQUACFV(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],const enum FF_IntParamForm *form,
                                   const double *T,const double x[],FF_SubsActivityData actData[]){
    double Ru=0;//The gas constant in different units
    double RuTinv;//1/(Ru*T) to speed calculation
    double meanR=0,meanQres=0;//mean molecular surface for residual part calculation
    double Phi[*numSubs],ThetaRes[*numSubs];//the mol fraction considering the polymer, and volume and surface fractions of all individual substances
    double tau[*numSubs][*numSubs];
    double S[*numSubs];
    int i,j;//the loop variables

    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);

    //common and combinatorial part calculation according to FV
    /*
    for (i=0;i<*numSubs;i++) xp[i]=x[i];//if there is no polymer the molar concentration is equal to the solvent or monomer concentrations
    for (i=0;i<*numSubs;i++){//Calculation of molar fractions passing the monomers to polymer
        if (baseProp[i].numMono>1){//If we find a polymer we recalculate the molar fractions
            for (j=0;j<*numSubs;j++){
                if (i!=j) xp[j]=xp[j]/((1-xp[i])+xp[i]/baseProp[i].numMono);//recalculation of the others molar fractions
            }
            xp[i]=xp[i]/((1-xp[i])*baseProp[i].numMono+xp[i]);//recalculation of the molar fraction for the polymer
        }
    }
    //printf("xp[0]:%f xp[1]:%f\n",xp[0],xp[1]);
    */
    for (i=0;i<*numSubs;i++){
        meanR=meanR+x[i]*baseProp[i].numMono*baseProp[i].FV;
        meanQres=meanQres+x[i]*baseProp[i].numMono*baseProp[i].qRes;
    }

    for (i=0;i<*numSubs;i++){
        Phi[i]=baseProp[i].numMono*baseProp[i].FV/meanR;//Free Volume ratio
        actData[i].lnGammaC=1-Phi[i]+log(Phi[i]);
        actData[i].lnGammaSG=0;
        //printf("GammaC: %f\n",exp(actData[i].lnGammaC));
    }

    //residual part
    for (i=0;i<*numSubs;i++){
        ThetaRes[i]=x[i]*baseProp[i].numMono*baseProp[i].qRes/meanQres;
        for(j=0;j<*numSubs;j++){
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C:
                tau[i][j]=exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]* *T * *T)*RuTinv);//polynomial form
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C:
                tau[i][j]=exp(-(pintParam[i][j][0]+pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T)*RuTinv);//polynomial form 2
                break;
            case FF_Ant3:
                tau[i][j]=exp(pintParam[i][j][0]+pintParam[i][j][1]/ *T+pintParam[i][j][2]*log(*T)+
                                    pintParam[i][j][3]* *T+pintParam[i][j][4]/ (*T * *T));//extended Antoine form
                break;
            }
        }
    }
    for (i=0;i<*numSubs;i++){
        S[i]=0;
        for(j=0;j<*numSubs;j++){
            S[i]=S[i]+ThetaRes[j]*tau[j][i];
            //S[i]=S[i]+Theta[j]*tau[j][i];
        }
    }
    for (i=0;i<*numSubs;i++){
        actData[i].lnGammaR=1-log(S[i]);
        for(j=0;j<*numSubs;j++){
            actData[i].lnGammaR=actData[i].lnGammaR-tau[i][j]*ThetaRes[j]/S[j];
            //lnGammaR[i]=lnGammaR[i]-tau[i][j]*Theta[j]/S[j];
        }
        actData[i].lnGammaR=baseProp[i].numMono*baseProp[i].q*actData[i].lnGammaR;
        //printf("GammaR: %f\n",exp(actData[i].lnGammaR));
    }
}

//Flory-Huggins. Temperature dependent chi? Use FV in combinatorial part? component 0 must be the solvent
void FF_ActivityFloryHuggins(const int *model,const  FF_BaseProp baseProp[],const double chiData[15][15][6],const int *form,const double *T,const double x[],FF_SubsActivityData actData[2]){
    double meanR;//mean molecular volume
    meanR=x[0]* baseProp[0].Vliq* baseProp[0].numMono+x[1]* baseProp[1].numMono*baseProp[1].Vliq;
    //printf("%f %f %f\n",baseProp[0].Vliq* baseProp[0].numMono,baseProp[1].numMono*baseProp[1].Vliq,meanR);
    double phi[2],chi;//the volume ratio of the solvent and the interaction parameter
    phi[0]=baseProp[0].Vliq* baseProp[0].numMono/meanR;
    phi[1]=baseProp[1].Vliq* baseProp[1].numMono/meanR;
    actData[0].lnGammaC=1-phi[0]+log(phi[0]);
    actData[1].lnGammaC=1-phi[1]+log(phi[1]);
    if(*model==FF_Chi) chi=chiData[0][1][0]+chiData[0][1][1]/ *T;
    //if sol.param. is in (cal/cm3)^0.5 multiply by 2045.5 to pass to Pa^0.5
    else if (*model==FF_Hildebrand) chi=0.35+baseProp[0].Vliq*pow(baseProp[0].Hildebrand-baseProp[1].Hildebrand,2)/(R* *T);
    else if (*model==FF_Hansen) chi=0.6*baseProp[0].Vliq*(pow(baseProp[0].HansenD-baseProp[1].HansenD,2)+0.25*pow(baseProp[0].HansenP-baseProp[1].HansenP,2)+
            0.25*pow(baseProp[0].HansenH-baseProp[1].HansenH,2))/(R* *T);
    //printf("chi:%f\n",chi);
    actData[0].lnGammaSG=0;
    actData[0].lnGammaR=chi*pow((x[1]*baseProp[1].numMono*baseProp[1].Vliq/meanR),2);
    actData[1].lnGammaSG=0;
    actData[1].lnGammaR=chi*pow((x[0]*baseProp[0].numMono*baseProp[0].Vliq/meanR),2);
}

/*
//Flory-Huggins. Use FV in combinatorial part?. Component 0 must be the solvent
void FF_ActivityFloryHugginsM(const int *numSubs,const  FF_BaseProp baseProp[],const double chiData[],const double *T,const double x[],FF_SubsActivityData actData[]){
    double meanR;//mean molecular volume

    meanR=x[0]* data->V[0]+x[1]* data->V[1];
    double phi,chi;//the volume ratio of the solvent
    phi=data->V[0]/meanR;
    *lnGammaC=1-phi+log(phi);
    if (data->model==FF_Hildebrand) *lnGammaR=(0.35+data->V[0]*pow(data->deltaHil[0]-data->deltaHil[1],2))*pow(data->V[1]/meanR,2);
    else if(data->model==FF_Hansen) *lnGammaR=data->V[0]*(pow(data->deltaHan[0][0]-data->deltaHan[1][0],2)+
            0.25*pow(data->deltaHan[0][1]-data->deltaHan[1][1],2)+0.25*pow(data->deltaHan[0][2]-data->deltaHan[1][2],2))*pow(data->V[1]/meanR,2)/1000;
    else if (data->model==FF_Chi){
        chi=data->chiData[0]+chiData[1]/ *T;
        *lnGammaR=chi*pow(x[1]*data->V[1]/meanR,2);
    }
}
*/


//Calculates the common data, independent of the molar fraction and temperature for the UNIFAC model.
//In data, we receive for each substance, its composition in number of each subgroup. So each register is: id substance, id subgroup, num of ocurrences
//for polymers the num of subgroups in data must be the multiplication of the monomer composition by the number of monomer units
//In uni we answer all the calculations
void CALLCONV FF_UNIFACParams(int numData, const int data[][3], const char *resDir, FF_UnifacData *uni){
    int ver=0;
    char path1[FILENAME_MAX]="";
    char path2[FILENAME_MAX]="";
    strcat(path1,resDir);
    strcat(path2,resDir);
    int i,j,k,newSubs,newSubg;
    //We need to count the different subgroups used, and to store their id
    int numSubgroups=0;
    //We need also to count the different substances, storing it original id
    int numSubs=0;
    int substance[20];//Will contain the original id(number of order) of the substances
    int numLines;//number of lines of interaction parameters file
    for (i=0;i<20;i++) for(j=0;j<30;j++) uni->subsSubg[i][j]=0;//We fill with 0 the conmposition of the substances in subgroups
    //We fill the first register with the first data
    uni->subgroup[0][0]=data[0][1];//this will contain the list of subgroups and their corresponding group
    numSubgroups=1;
    substance[0]=data[0][0];
    numSubs=1;
    uni->subsSubg[0][0]=data[0][2];//this will contain the ocurrence of subgroups in the substances

    for (i=1; i<numData;i++){//We follow all the descriptors
        newSubs=1;//We say that the substance is new
        for (j=0;j<numSubs;j++){//We follow all the already defined substances
            if (data[i][0]==substance[j]){
                newSubs=0;//We say that the substance is not new
                break;
            }
        }
        if (newSubs){//If not found we add the substance
            j=numSubs;//j must contain the index of the substance
            substance[j]=data[i][0];
            numSubs++;
        }
        newSubg=1;//we say the subgroup is new
        for (k=0; k<numSubgroups;k++){//We follow the already defined subgroups
            if (data[i][1]==uni->subgroup[k][0]){
                newSubg=0;//we say that the subgroup is not new
                break;
            }
        }
        if (newSubg){//if not found we add the subgroups
            k=numSubgroups;
            uni->subgroup[k][0]=data[i][1];
            numSubgroups++;
        }
        uni->subsSubg[j][k]=data[i][2];
        if (ver==1) printf("Substance:%i Subgroup:%i Number:%i\n",j,k,uni->subsSubg[j][k]);
    }
    //we recover from the corresponding file the subgroups information, and fill the subgroup and subgData table
    FILE *f;
    unsigned sg,g,g1,g2;
    double A12,B12,C12,A21,B21,C21;
    double r,q;

    /*if (uni->model==FF_UNIFACStd) f=fopen("Data/UnifacSubgStd.txt","r");
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_EntropicFV)||(uni->model==FF_UNIFACZM)) f=fopen("Data/UnifacSubgPSRK.txt","r");
    else if (uni->model==FF_UNIFACDort) f=fopen("Data/UnifacSubgDort.txt","r");
    else if (uni->model==FF_UNIFACNist) f=fopen("Data/UnifacSubgNist.txt","r");*/
    if (uni->model==FF_UNIFACStd) strcat(path1,"Data/UnifacSubgStd.txt");
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_EntropicFV)||(uni->model==FF_UNIFACZM)) strcat(path1,"Data/UnifacSubgPSRK.txt");
    else if (uni->model==FF_UNIFACDort) strcat(path1,"Data/UnifacSubgDort.txt");
    else if (uni->model==FF_UNIFACNist) strcat(path1,"Data/UnifacSubgNist.txt");
    f=fopen(path1,"r");
    for (i=0;i<numSubgroups;i++){
        do{
            fscanf(f,"%3lu%3lu%6lf%6lf\n",&sg,&g,&r,&q);
            if(uni->subgroup[i][0]==sg){
                uni->subgroup[i][1]=g;
                uni->subgData[i][0]=r;
                uni->subgData[i][1]=q;
                break;
            }
        } while (sg<300);
        if (ver==1) printf("Index:%i Subgroup:%i Group:%i r:%f q:%f\n",i,uni->subgroup[i][0],uni->subgroup[i][1],uni->subgData[i][0],uni->subgData[i][1]);
        rewind(f);
    }
    fclose(f);

    //we recover now the subgroup interaction and fill the subgroupInter table
    for (i=0;i<numSubgroups;i++){
        for (j=0;j<numSubgroups;j++){
            uni->subgInt[i][j][0]=0;
            uni->subgInt[i][j][1]=0;
            uni->subgInt[i][j][2]=0;
        }
    }
    if (uni->model==FF_UNIFACStd){
        numLines=635;
        strcat(path2,"Data/UnifacInterStd.txt");
        f=fopen(path2,"r");
        if (f==NULL) printf("Error opening Data/UnifacInterStd.txt\n");
    }
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_EntropicFV)||(uni->model==FF_UNIFACZM)){
        numLines=956;
        strcat(path2,"Data/UnifacInterPSRK.txt");
        f=fopen(path2,"r");
        if (f==NULL) printf("Error opening Data/UnifacInterPSRK.txt\n");
    }
    else if ((uni->model==FF_UNIFACDort)){
        numLines=756;
        strcat(path2,"Data/UnifacInterDort.txt");
        f=fopen(path2,"r");
        if (f==NULL) printf("Error opening data/UnifacInterDort.txt\n");
    }
    else if (uni->model==FF_UNIFACNist){
        numLines=1969;
        strcat(path2,"Data/UnifacInterNist.txt");
        f=fopen(path2,"r");
        if (f==NULL) printf("Error opening data/UnifacInterNist.txt\n");
    }
    if (f==NULL) printf("Error\n");
    if (uni->model==FF_UNIFACStd){
        for (i=0;i<numSubgroups;i++){
            for (j=0;j<numSubgroups;j++){
                if (uni->subgroup[i][1]<uni->subgroup[j][1]){
                    for(k=1;k<numLines;k++){
                        fscanf(f,"%03lu%03lu%10lf%10lf\n",&g1,&g2,&A12,&A21);
                        //if (ver==1) printf("%i %i %f %f\n",g1,g2,A12,A21);
                        if ((g1==uni->subgroup[i][1])&& (g2==uni->subgroup[j][1])){
                            uni->subgInt[i][j][0]=A12;
                            uni->subgInt[j][i][0]=A21;
                            if (ver==1) printf("Unifac Std GroupA:%i GroupB:%i %f %f\n",g1,g2,uni->subgInt[i][j][0],uni->subgInt[j][i][0]);
                            break;
                        }
                    }
                    rewind(f);
                }
            }
        }
    }
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_UNIFACDort)||(uni->model==FF_UNIFACNist)||(uni->model==FF_UNIFACZM)||(uni->model==FF_EntropicFV)){
        for (i=0;i<numSubgroups;i++){
            for (j=0;j<numSubgroups;j++){
                if (uni->subgroup[i][1]<uni->subgroup[j][1]){
                    for(k=1;k<numLines;k++){
                        fscanf(f,"%03lu%03lu%10lf%10lf%10lf%10lf%10lf%10lf\n",&g1,&g2,&A12,&B12,&C12,&A21,&B21,&C21);
                        //if (ver==1) printf("%i %i %f %f\n",g1,g2,A12,A21);
                        if ((g1==uni->subgroup[i][1])&& (g2==uni->subgroup[j][1])){
                            uni->subgInt[i][j][0]=A12;
                            uni->subgInt[i][j][1]=B12;
                            uni->subgInt[i][j][2]=C12;
                            uni->subgInt[j][i][0]=A21;
                            uni->subgInt[j][i][1]=B21;
                            uni->subgInt[j][i][2]=C21;
                            if (ver==1) printf("Unifac thermic GroupA:%i GroupB:%i %f %f %f %f %f %f\n",g1,g2,uni->subgInt[i][j][0],uni->subgInt[i][j][1],uni->subgInt[i][j][2],
                                    uni->subgInt[j][i][0],uni->subgInt[j][i][1],uni->subgInt[j][i][2]);
                            break;
                        }
                    }
                    rewind(f);
                }
            }
        }
    }
    fclose(f);

    for (i=0; i<numSubs;i++){
        uni->subsR[i]=0;
        uni->subsQ[i]=0;
        for (j=0;j<numSubgroups;j++){
            uni->subsR[i]=uni->subsR[i]+uni->subsSubg[i][j]*uni->subgData[j][0];
            uni->subsQ[i]=uni->subsQ[i]+uni->subsSubg[i][j]*uni->subgData[j][1];
        }
    }
    uni->numSubs=numSubs;
    uni->numSubg=numSubgroups;
}

//Calculates activity coefficients according to UNIFAC models, at given T and composition.
//For polymers x q,r,FV are those of the polymer.
void CALLCONV FF_ActivityUNIFAC(FF_UnifacData *data, const double *T, const double x[],FF_SubsActivityData actData[]){
    int ver=0;
    int i,j,k;
    double aux;
    //Combinatorial part calculation
    double meanR=0, meanQ=0,theta[data->numSubs],phi[data->numSubs];
    if (data->model==FF_EntropicFV){
        for (i=0;i<data->numSubs;i++){
            meanR=meanR+data->FV[i]*x[i];
            meanQ=meanQ+data->subsQ[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            phi[i]=data->FV[i]/meanR;//Free Volume ratio
            theta[i]=data->subsQ[i]/meanQ;//Area fraction
        }
    }
    else{
        for (i=0;i<data->numSubs;i++){
            meanR=meanR+data->subsR[i]*x[i];
            meanQ=meanQ+data->subsQ[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            phi[i]=data->subsR[i]/meanR;//Volume fraction/x[i]
            theta[i]=data->subsQ[i]/meanQ;//Area fraction
        }
    }
    if ((data->model==FF_UNIFACStd)||(data->model==FF_UNIFACPSRK)){
            for (i=0;i<data->numSubs;i++){
                actData[i].lnGammaC=1-phi[i]+log(phi[i]);
                //lnGammaC[i]=log(phi[i])+5*data->subsQ[i]*log(theta[i]/phi[i])+data->subsL[i]-phi[i]*meanL;//alternative from Fredenslund
                actData[i].lnGammaSG=-5*data->subsQ[i]*(1-phi[i]/theta[i]+log(phi[i]/theta[i]));
                if (ver==1) printf("Unifac Std/PSRK lnGammaC[%i]:%f d\n",i,actData[i].lnGammaC);
                if (ver==1) printf("Unifac Std/PSRK lnGammaSG[%i]:%f dlnGammaSG:%f\n",i,actData[i].lnGammaSG);
            }
    }
    else if((data->model==FF_EntropicFV)){
        for (i=0;i<data->numSubs;i++){
            actData[i].lnGammaC=1-phi[i]+log(phi[i]);
            actData[i].lnGammaSG=0;
            if (ver==1) printf("Entropic FV gammaC[%i]:%f\n",i,exp(actData[i].lnGammaC));
        }
    }
    else if ((data->model==FF_UNIFACDort)||(data->model==FF_UNIFACNist)){
        double meanR2=0, phi2[data->numSubs];
        for (i=0;i<data->numSubs;i++) meanR2=meanR2+pow(data->subsR[i],0.75)*x[i];
        for (i=0;i<data->numSubs;i++){
            phi2[i]=pow(data->subsR[i],0.75)/meanR2;//Volume fraction /x[i]
            actData[i].lnGammaC=1-phi2[i]+log(phi2[i]);
            actData[i].lnGammaSG=-5*data->subsQ[i]*(1-phi[i]/theta[i]+log(phi[i]/theta[i]));
            if (ver==1) printf("Unifac Dortmund/Nist lnGammaC[%i]:%f \n",i,actData[i].lnGammaC);
            if (ver==1) printf("Unifac Dortmund/Nist llnGammaSG[%i]:%f \n",i,actData[i].lnGammaSG);
        }
    }
    else if (data->model==FF_UNIFACZM){
        double meanR2=0, phi2[data->numSubs];
        for (i=0;i<data->numSubs;i++){
            if (data->subsR[i]>65) meanR2=meanR2+data->subsR[i]*0.6583*x[i];//we consider a polymer if R is over 65
            else meanR2=meanR2+data->subsR[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            if (data->subsR[i]>65) phi2[i]=data->subsR[i]*0.6583/meanR2;//Volume fraction /x[i]
            else phi2[i]=data->subsR[i]/meanR2;
            actData[i].lnGammaC=1-phi2[i]+log(phi2[i]);
            actData[i].lnGammaSG=-5*data->subsQ[i]*(1-phi[i]/theta[i]+log(phi[i]/theta[i]));
            if (ver==1) printf("Unifac ZM lnGammaC[%i]:%f\n",i,actData[i].lnGammaC);
        }
    }
    //Residual part calculation
    //first we calculate the interaction between subgroups, as a temperature function
    double psi[data->numSubg][data->numSubg], substSubgrTheta[data->numSubs][data->numSubg],substSubgrLambda[data->numSubs][data->numSubg],
            substSubgrSum[data->numSubs][data->numSubg];
    if (data->model==FF_UNIFACStd){
        for (i=0;i<data->numSubg;i++) for (j=0;j<data->numSubg;j++) psi[i][j]=exp(-data->subgInt[i][j][0]/ *T);
        if (ver==1) printf("Unifac standard subgroupA index:%i subgroupB index:%i psi:%f\n",i,j,psi[i][j]);
    }
    else if ((data->model==FF_UNIFACPSRK)||(data->model==FF_UNIFACDort)||(data->model==FF_UNIFACNist)||(data->model==FF_EntropicFV)||(data->model==FF_UNIFACZM)){
        for (i=0;i<data->numSubg;i++) for (j=0;j<data->numSubg;j++){
            psi[i][j]=exp(-data->subgInt[i][j][0]/ *T - data->subgInt[i][j][1] - data->subgInt[i][j][2]* *T);
            if (ver==1) printf("Unifac complex subgroupA index:%i subgroupB index:%i a:%f b:%f c:%f\n",i,j,data->subgInt[i][j][0],data->subgInt[i][j][1],
                    data->subgInt[i][j][2]);
            if (ver==1) printf("Unifac complex subgroupA index:%i subgroupB index:%i psi:%f\n",i,j,psi[i][j]);
        }
    }
    //now the calculation of theta and lambda por each subgroup in each substance
    for (i=0;i<data->numSubs;i++){
        for (j=0;j<data->numSubg;j++) substSubgrTheta[i][j]=data->subsSubg[i][j]*data->subgData[j][1]/data->subsQ[i];
    }
    for (i=0;i<data->numSubs;i++) for (j=0;j<data->numSubg;j++){
        substSubgrSum[i][j]=0;
        for (k=0;k<data->numSubg;k++) substSubgrSum[i][j]=substSubgrSum[i][j]+substSubgrTheta[i][k]*psi[k][j];
    }
    for (i=0;i<data->numSubs;i++) for (j=0;j<data->numSubg;j++){
        substSubgrLambda[i][j]=0;
        if (data->subsSubg[i][j]>0){
            for (k=0;k<data->numSubg;k++) substSubgrLambda[i][j]=substSubgrLambda[i][j]+substSubgrTheta[i][k]*psi[k][j];
            substSubgrLambda[i][j]=1-log(substSubgrLambda[i][j]);
            for (k=0;k<data->numSubg;k++) substSubgrLambda[i][j]=substSubgrLambda[i][j]-substSubgrTheta[i][k]*psi[j][k]/substSubgrSum[i][k];
             substSubgrLambda[i][j]=data->subgData[j][1]* substSubgrLambda[i][j];
        }
    }
    //Now we calculate theta and lambda of the subgroups in the whole mixture
    double subgrTheta[data->numSubg],subgrSum[data->numSubg],subgrLambda[data->numSubg];
    for (j=0;j<data->numSubg;j++){
        subgrTheta[j]=0;
        for (i=0;i<data->numSubs;i++){
            aux=data->subgData[j][1]*x[i]*data->subsSubg[i][j]/meanQ;
            subgrTheta[j]=subgrTheta[j]+aux;
        }
    }
    for (j=0;j<data->numSubg;j++){
        subgrSum[j]=0;
        for (k=0;k<data->numSubg;k++) subgrSum[j]=subgrSum[j]+subgrTheta[k]*psi[k][j];
    }
    for (j=0;j<data->numSubg;j++){
        //subgrLambda[j]=0;
        //for (k=0;k<data->numSubg;k++) subgrLambda[j]=subgrLambda[j]+ subgrTheta[k]*psi[k][j];
        //subgrLambda[j]=1-log(subgrLambda[j]);
        subgrLambda[j]=1-log(subgrSum[j]);
        for (k=0;k<data->numSubg;k++) subgrLambda[j]=subgrLambda[j]-subgrTheta[k]*psi[j][k]/subgrSum[k];
        subgrLambda[j]=data->subgData[j][1]*subgrLambda[j];
    }
    for (i=0;i<data->numSubs;i++){
        actData[i].lnGammaR=0;
        for (j=0;j<data->numSubg;j++) actData[i].lnGammaR=actData[i].lnGammaR+data->subsSubg[i][j]*(subgrLambda[j]-substSubgrLambda[i][j]);
        actData[i].gamma=exp(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR);
        if (ver==1) printf("lnGammaR[%i]:%f\n",i,actData[i].lnGammaR);
    }
}

//Calculates activity coefficents according any of the defined models from a FF_MixData structure
EXP_IMP void CALLCONV FF_Activity(FF_MixData *mix,const double *T,double x[],FF_SubsActivityData actData[]){
    int i;
    switch(mix->actModel){
    case FF_Wilson:
        FF_ActivityWilson(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_NRTL:
        FF_ActivityNRTL(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIQUAC:
        FF_ActivityUNIQUAC(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIQUACFV:
        FF_ActivityUNIQUACFV(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_Hildebrand:
    case FF_Hansen:
    case FF_Chi:
        FF_ActivityFloryHuggins(&mix->actModel,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIFACStd:
        FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
        break;
    case FF_UNIFACPSRK://As the data is shared it is necessary to specify the calculation model
        mix->unifPSRKData.model=FF_UNIFACPSRK;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_UNIFACZM:
        mix->unifPSRKData.model=FF_UNIFACZM;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_EntropicFV:
        mix->unifPSRKData.model=FF_EntropicFV;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_UNIFACDort:
        FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
        break;
    case FF_UNIFACNist:
        FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
        break;
    }
    for(i=0;i<mix->numSubs;i++)actData[i].gamma=exp(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR);
}


//Calculates fugacity and activity coefficients, at given T,P and composition, from an activity model
void CALLCONV FF_PhiAndActivity(FF_MixData *mix,const double *T,const double *P,const double x[],FF_SubsActivityData actData[],double phi[]){
    int i;
    for(i=0;i<mix->numSubs;i++){
        if(*T>=mix->baseProp[i].Tc){
            phi[i]=0;
            return;
        }
    }
    //First the calculation of the activity coefficient
    switch(mix->actModel){
    case FF_Wilson:
        FF_ActivityWilson(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_NRTL:
        FF_ActivityNRTL(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIQUAC:
        FF_ActivityUNIQUAC(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIQUACFV:
        FF_ActivityUNIQUACFV(&mix->numSubs,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_Hildebrand:
    case FF_Hansen:
    case FF_Chi:
        FF_ActivityFloryHuggins(&mix->actModel,mix->baseProp,mix->intParam,&mix->intForm,T,x,actData);
        break;
    case FF_UNIFACStd:
        FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
        break;
    case FF_UNIFACPSRK://As the data is shared it is necessary to specify the calculation model
        mix->unifPSRKData.model=FF_UNIFACPSRK;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_UNIFACZM:
        mix->unifPSRKData.model=FF_UNIFACZM;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_EntropicFV:
        mix->unifPSRKData.model=FF_EntropicFV;
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        break;
    case FF_UNIFACDort:
        FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
        break;
    case FF_UNIFACNist:
        FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
        break;
    }
    //Calculation of gamma
    for(i=0;i<mix->numSubs;i++) actData[i].gamma=exp(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR);

    //Calculation of the reference fugacity at same T and P for each component
    if (mix->refVpEos==1){//Reference fugacity by Eos
        char option='l',state;
        double answerL[3],answerG[3];
        double phiL0[mix->numSubs];
        switch (mix->eosType)
        {
        case FF_SAFTtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->saftData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                }
                else phi[i]=0;
            }
            break;
        /*case FF_SWtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->swData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                }
                else phi[i]=0;
            }
            break;*/
        case FF_CubicType:
        case FF_CubicPRtype:
        case FF_CubicSRKtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->cubicData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                    //printf("phiRef[%i]: %f\n",i,phiL0[i]);
                }
                else phi[i]=0;
            }
            break;
        }
    }
    //In case of using Vp, it makes no sense to calculate the reference fugacity at Vp, and later the Poynting factor, to change to actual pressure
    //If we are able to calculate the reference fugacity at Vp we can do directly its calculation at actual pressure.Although this is perhaps not always true.
    //It makes sense only if assuming that reference fugacity is Vp and phi at saturation =1
    else{
        double Vp,lDens,Tref,rhoRef;
        int nPoints=1;
        for(i=0;i<mix->numSubs;i++){
            if(mix->baseProp[i].numMono<20){
                FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,*T,&Vp);
                if(mix->lDensCorr[i].form>0) FF_PhysPropCorrM(mix->lDensCorr[i].form,mix->lDensCorr[i].coef,mix->baseProp[i].MW,*T,&lDens);
                else FF_LiqDensSatRackett(&mix->baseProp[i],Tref,rhoRef,*T,&lDens);
              //printf("Vp:%f\n",Vp);
                phi[i]=actData[i].gamma*Vp/ *P*exp(mix->baseProp[i].MW*0.001*(*P-Vp)/(lDens*R* *T));
            }
            else phi[i]=0;
        }
    }
}


//Calculates fugacity from given activity coefficients, at given T,P and composition, using an activity model
void CALLCONV FF_PhiFromActivity(FF_MixData *mix,const double *T,const double *P,const double x[],FF_SubsActivityData actData[],double phi[]){
    int i;
    for(i=0;i<mix->numSubs;i++){
        if(*T>=mix->baseProp[i].Tc){
            phi[i]=0;
            return;
        }
    }
    //Calculation of the reference fugacity at same T and P for each component
    if (mix->refVpEos==1){//Reference fugacity by Eos
        char option='l',state;
        double answerL[3],answerG[3];
        double phiL0[mix->numSubs];
        switch (mix->eosType)
        {
        case FF_SAFTtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->saftData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                }
                else phi[i]=0;
            }
            break;
        /*case FF_SWtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->swData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                }
                else phi[i]=0;
            }
            break;*/
        case FF_CubicType:
        case FF_CubicPRtype:
        case FF_CubicSRKtype:
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].numMono<20){
                    FF_VfromTPeos(mix->eosType,*T,*P,&mix->cubicData[i],option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*actData[i].gamma;
                    //printf("phiRef[%i]: %f\n",i,phiL0[i]);
                }
                else phi[i]=0;
            }
            break;
        }
    }
    //In case of using Vp, it makes no sense to calculate the reference fugacity at Vp, and later the Poynting factor, to change to actual pressure
    //If we are able to calculate the reference fugacity at vp we can do directly its calculation at actual pressure.Although this is perhaps not always true.
    //It makes sense only if assuming that reference fugacity is Vp
    else{
        double Vp,lDens,Tref,rhoRef;
        for(i=0;i<mix->numSubs;i++){
            if(mix->baseProp[i].numMono<20){
                FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,*T,&Vp);
                if(mix->lDensCorr[i].form>0){
                    FF_PhysPropCorrM(mix->lDensCorr[i].form,mix->lDensCorr[i].coef,mix->baseProp[i].MW,*T,&lDens);
                }
                else FF_LiqDensSatRackett(&mix->baseProp[i],Tref,rhoRef,*T,&lDens);
              //printf("Vp:%f\n",Vp);
                phi[i]=actData[i].gamma*Vp/ *P*exp(mix->baseProp[i].MW*0.001*(*P-Vp)/(lDens*R* *T));
            }
            else phi[i]=0;
        }
    }
}

//Calculates ln of activities and the derivatives of gE
void CALLCONV FF_UNIFACDerivatives(FF_UnifacData *data, const double *T, const double x[],FF_SubsActivityData actData[],FF_ExcessData *excData){
    int i,j;
    excData->gEC=0;excData->gESG=0;excData->gER=0;
    double gECplus=0,gESGplus=0,gERplus=0,gEplus=0;
    double xPlus[data->numSubs];
    double Tplus=*T+0.01;
    double hE;
    FF_SubsActivityData actDataPlus[data->numSubs];
    //Calculation of the ln of activity coefficients combinatorial, SG and residual
    FF_ActivityUNIFAC(data,T,x,actData);
    //Calculation of excess g for combinatorial, Sg and residual parts
    for(i=0;i<data->numSubs;i++){
        excData->gEC=excData->gEC+x[i]*actData[i].lnGammaC;
        excData->gESG=excData->gESG+x[i]*actData[i].lnGammaSG;
        excData->gER=excData->gER+x[i]*actData[i].lnGammaR;
        xPlus[i]=x[i];
        actData[i].dgEC=actData[i].dgESG=actData[i].dgER=0;
    }
    excData->gE=excData->gEC+excData->gESG+excData->gER;// Calculate the global gE
    //Numeric partial derivatives of gE combinatorial, SG and residual regarding substances
    for (i=0;i<data->numSubs;i++){
        xPlus[i]=x[i]*1.001;
        FF_ActivityUNIFAC(data,T,xPlus,actDataPlus);
        gECplus=gESGplus=gERplus=0;
        for (j=0;j<data->numSubs;j++){
            gECplus=gECplus+xPlus[j]*actDataPlus[j].lnGammaC;
            gESGplus=gESGplus+xPlus[j]*actDataPlus[j].lnGammaSG;
            gERplus=gERplus+xPlus[j]*actDataPlus[j].lnGammaR;
        }
        actData[i].dgEC=(gECplus-excData->gEC)/(xPlus[i]-x[i]);
        actData[i].dgESG=(gESGplus-excData->gESG)/(xPlus[i]-x[i]);
        actData[i].dgER=(gERplus-excData->gER)/(xPlus[i]-x[i]);
        xPlus[i]=x[i];
        //printf("gamma[%i]: %f\n",i,exp(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR));
    }
    //Calculation of ln of act. coef. a somewhat higher temp.
    FF_ActivityUNIFAC(data,&Tplus,x,actDataPlus);
    gEplus=0;
    for (i=0;i<data->numSubs;i++) gEplus=gEplus+x[i]*(actDataPlus[i].lnGammaC+actDataPlus[i].lnGammaSG+actDataPlus[i].lnGammaR);
    excData->dgE_dT=((gEplus)-(excData->gE))/(Tplus- *T);
    //printf("hE:%f\n",-R* *T* *T* excData->dgE_dT);
}



//Calculates ln of activities and the derivatives of gE
void CALLCONV FF_ActivityDerivatives(const int *actModel,const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[15][15][6],const int *form,
                                        const double *T,const double x[],FF_SubsActivityData actData[],FF_ExcessData *excData){
    int i,j;
    excData->gEC=0;excData->gESG=0;excData->gER=0;//Excess data for reduced G
    double gECplus=0,gESGplus=0,gERplus,gEplus=0;//Values after increment of x[i] or T
    double xPlus[*numSubs];//Values after increment of x[i]
    double Tplus=*T+0.01;//Values after increment of T
    double hE;
    FF_SubsActivityData actDataPlus[*numSubs];//Where to store the activity coeff. after incrementing x[i] or T
    /*void (*ActCalc)(int *,FF_BaseProp *,double *,int *,double *,double *,FF_SubsActivityData *);
    if (*actModel==FF_UNIQUAC) ActCalc=FF_ActivityUNIQUAC;
    else if (*actModel==FF_NRTL) ActCalc=&FF_ActivityNRTL;
    else if (*actModel==FF_Wilson) ActCalc=&FF_ActivityWilson;*/

    //Calculation of the ln of activity coefficients combinatorial, SG and residual
    //ActCalc(numSubs,baseProp,pintParam,form,T,x,actData);
    if (*actModel==FF_UNIQUAC) FF_ActivityUNIQUAC(numSubs,baseProp,pintParam,form,T,x,actData);
    else if (*actModel==FF_NRTL) FF_ActivityNRTL(numSubs,baseProp,pintParam,form,T,x,actData);
    else if (*actModel==FF_Wilson)FF_ActivityWilson(numSubs,baseProp,pintParam,form,T,x,actData);
    //Calculation of excess g for combinatorial, Sg and residual parts
    for(i=0;i<*numSubs;i++){
        excData->gEC=excData->gEC+x[i]*actData[i].lnGammaC;
        excData->gESG=excData->gESG+x[i]*actData[i].lnGammaSG;
        excData->gER=excData->gER+x[i]*actData[i].lnGammaR;
        xPlus[i]=x[i];//Copy composition
        actData[i].dgEC=actData[i].dgESG=actData[i].dgER=0;//Put also all partial derivatives at 0
    }
    excData->gE=excData->gEC+excData->gESG+excData->gER;// Calculate the global gE
    //Numeric calculation of partial derivatives of gE combinatorial, SG and residual
    for (i=0;i<*numSubs;i++){
        gECplus=gESGplus=gERplus=0;
        xPlus[i]=x[i]*1.001;
        //ActCalc(numSubs,baseProp,pintParam,form,T,xPlus,actDataPlus);
        if (*actModel==FF_UNIQUAC) FF_ActivityUNIQUAC(numSubs,baseProp,pintParam,form,T,xPlus,actDataPlus);
        else if (*actModel==FF_NRTL) FF_ActivityNRTL(numSubs,baseProp,pintParam,form,T,xPlus,actDataPlus);
        else if (*actModel==FF_Wilson)FF_ActivityWilson(numSubs,baseProp,pintParam,form,T,xPlus,actDataPlus);
        for (j=0;j<*numSubs;j++){
            gECplus=gECplus+xPlus[j]*actDataPlus[j].lnGammaC;
            gESGplus=gESGplus+xPlus[j]*actDataPlus[j].lnGammaSG;
            gERplus=gERplus+xPlus[j]*actDataPlus[j].lnGammaR;
        }
        actData[i].dgEC=(gECplus-excData->gEC)/(xPlus[i]-x[i]);
        actData[i].dgESG=(gESGplus-excData->gESG)/(xPlus[i]-x[i]);
        actData[i].dgER=(gERplus-excData->gER)/(xPlus[i]-x[i]);
        xPlus[i]=x[i];
    }
    //Calculation of act. coef. a somewhat higher temp.
    //ActCalc(numSubs,baseProp,pintParam,form,&Tplus,x,actDataPlus);
    if (*actModel==FF_UNIQUAC) FF_ActivityUNIQUAC(numSubs,baseProp,pintParam,form,&Tplus,x,actDataPlus);
    else if (*actModel==FF_NRTL) FF_ActivityNRTL(numSubs,baseProp,pintParam,form,&Tplus,x,actDataPlus);
    else if (*actModel==FF_Wilson)FF_ActivityWilson(numSubs,baseProp,pintParam,form,&Tplus,x,actDataPlus);
    gEplus=0;
    for (i=0;i<*numSubs;i++) gEplus=gEplus+x[i]*(actDataPlus[i].lnGammaC+actDataPlus[i].lnGammaSG+actDataPlus[i].lnGammaR);
    excData->dgE_dT=((gEplus)-(excData->gE))/(Tplus- *T);
    //printf("hE:%f\n",-R* *T* excData->dgE_dT);
}




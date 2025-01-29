/*

 * FFeosMix.c
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

// contains EOS calculations for pure substances
//==============================================

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFeosMix.h"
#include "FFactivity.h"

//Get mix data from an exported file
EXP_IMP FF_MixData * CALLCONV FF_MixDataFromFile(const char *name){
    FF_MixData *mixData = (FF_MixData*) calloc(1,sizeof(FF_MixData));
    char path[FILENAME_MAX]="Data/";
    strcat(path,name);
    strcat(path,".md");
    FILE * file= fopen(path, "rb");
    if (file != NULL) {
        fread(mixData, sizeof(FF_MixData), 1, file);
        fclose(file);
    }
    else printf("Mixture data file not found\n");
    //if ((subsData->baseProp.FV==0)&&(subsData->baseProp.Vliq>0)&&(subsData->baseProp.VdWV>0)) subsData->baseProp.FV=subsData->baseProp.Vliq-
    //        1.2*subsData->baseProp.VdWV;
    return mixData;
}

//Create a mixture data structure from an array of substance data structures
EXP_IMP FF_MixData * CALLCONV FF_MixDataFromSubsData(int numSubs,const FF_SubstanceData *subsData[]){
    if (numSubs>15) numSubs=15;
    //FF_MixData *mixData =(FF_MixData*) calloc(1,sizeof(FF_MixData));
    static FF_MixData mixData;
    int h,i,j=0,k=0,l=0,m=0,n=0;
    int unifacCompStd[numSubs*10][3];//To contain the Unifac description for all substances
    int unifacCompPSRK[numSubs*10][3];
    int unifacCompDort[numSubs*10][3];
    int unifacCompNist[numSubs*10][3];
    int useUnifacStd=1,useUnifacPSRK=1,useUnifacDort=1,useUnifacNist=1;//Indicates that can be used in the mix. If not, we will put to 0

    /*for(i=0;i<numSubs*10;i++){
        unifacCompStd[i][0]=unifacCompStd[i][1]=unifacCompStd[i][2]=0;
    }*/

    for(h=0;h<numSubs*10;h++) for(i=0;i<3;i++) {
        unifacCompStd[h][i]=0;
        unifacCompPSRK[h][i]=0;
        unifacCompDort[h][i]=0;
        unifacCompNist[h][i]=0;
    }


    for(h=0;h<15;h++)for(i=0;i<15;i++)for(j=0;j<6;j++) mixData.intParam[h][i][j]=0;
    mixData.numSubs=numSubs;
    for (i=0;i<numSubs;i++){
        mixData.id[i]=subsData[i]->id;
        strcpy(mixData.subsName[i],subsData[i]->name);
        strncpy(mixData.CAS[i],subsData[i]->CAS,22);
        mixData.baseProp[i]=subsData[i]->baseProp;
        if((mixData.baseProp[i].MWmono>0)&&(mixData.baseProp[i].MW>0)) mixData.baseProp[i].numMono=557;//mixData.baseProp[i].numMono=mixData.baseProp[i].MW/mixData.baseProp[i].MWmono;
        else mixData.baseProp[i].numMono=1;
        //In FF_UnifacData we will have the FV of the polymer and in baseProp that of the monomer
        mixData.unifStdData.FV[i]=mixData.baseProp[i].numMono*subsData[i]->baseProp.FV;
        mixData.unifPSRKData.FV[i]=mixData.baseProp[i].numMono*subsData[i]->baseProp.FV;
        mixData.unifDortData.FV[i]=mixData.baseProp[i].numMono*subsData[i]->baseProp.FV;
        mixData.unifNistData.FV[i]=mixData.baseProp[i].numMono*subsData[i]->baseProp.FV;
        mixData.cubicData[i]=subsData[i]->cubicData;
        mixData.saftData[i]=subsData[i]->saftData;
        mixData.RI[i]=subsData[i]->RI;
        mixData.cp0[i]=subsData[i]->cp0;
        mixData.vp[i]=subsData[i]->vp;
        mixData.hVsat[i]=subsData[i]->hVsat;
        mixData.lCp[i]=subsData[i]->lCp;
        mixData.lDens[i]=subsData[i]->lDens;
        mixData.lVisc[i]=subsData[i]->lVisc;
        mixData.lThC[i]=subsData[i]->lThC;
        mixData.lSurfT[i]=subsData[i]->lSurfT;
        mixData.cp0Corr[i]=subsData[i]->cp0Corr;
        mixData.vpCorr[i]=subsData[i]->vpCorr;
        mixData.btCorr[i]=subsData[i]->btCorr;
        mixData.hVsatCorr[i]=subsData[i]->hVsatCorr;
        mixData.lCpCorr[i]=subsData[i]->lCpCorr;
        mixData.lDensCorr[i]=subsData[i]->lDensCorr;
        mixData.lViscCorr[i]=subsData[i]->lViscCorr;
        mixData.lThCCorr[i]=subsData[i]->lThCCorr;
        mixData.lSurfTCorr[i]=subsData[i]->lSurfTCorr;
        mixData.gDensCorr[i]=subsData[i]->gDensCorr;
        mixData.gViscCorr[i]=subsData[i]->gViscCorr;
        mixData.gThCCorr[i]=subsData[i]->gThCCorr;
        mixData.unifStdData.model=FF_UNIFACStd;
        mixData.unifPSRKData.model=FF_UNIFACPSRK;
        mixData.unifDortData.model=FF_UNIFACDort;
        mixData.unifNistData.model=FF_UNIFACNist;

        j=0;
        if(subsData[i]->UnifStdSubg[j][1]<=0) useUnifacStd=0;
        while(subsData[i]->UnifStdSubg[j][1]>0){//Put the subgroups description of all substances in a single array
            unifacCompStd[k][0]=i;
            unifacCompStd[k][1]=subsData[i]->UnifStdSubg[j][0];
            if(mixData.baseProp[i].MWmono>1) unifacCompStd[k][2]=mixData.baseProp[i].numMono*subsData[i]->UnifStdSubg[j][1];
            else unifacCompStd[k][2]=subsData[i]->UnifStdSubg[j][1];
            //printf("Std %i %i %i\n",unifacCompStd[k][0],unifacCompStd[k][1],unifacCompStd[k][2]);
            j++;
            k++;
        }
        j=0;
        if(subsData[i]->UnifPSRKSubg[j][1]<=0) useUnifacPSRK=0;
        while(subsData[i]->UnifPSRKSubg[j][1]>0){//Put the subgroups description of all substances in a single array
            unifacCompPSRK[l][0]=i;
            unifacCompPSRK[l][1]=subsData[i]->UnifPSRKSubg[j][0];
            if(mixData.baseProp[i].MWmono>1) unifacCompPSRK[l][2]=mixData.baseProp[i].numMono*subsData[i]->UnifPSRKSubg[j][1];
            else unifacCompPSRK[l][2]=subsData[i]->UnifPSRKSubg[j][1];
            //printf("PSRK %i %i %i\n",unifacCompPSRK[l][0],unifacCompPSRK[l][1],unifacCompPSRK[l][2]);
            j++;
            l++;
        }
        j=0;
        if(subsData[i]->UnifDortSubg[j][0]<=0) useUnifacDort=0;
        while(subsData[i]->UnifDortSubg[j][0]>0){
            unifacCompDort[m][0]=i;
            unifacCompDort[m][1]=subsData[i]->UnifDortSubg[j][0];
            if(mixData.baseProp[i].MWmono>1) unifacCompDort[m][2]=mixData.baseProp[i].numMono*subsData[i]->UnifDortSubg[j][1];
            else unifacCompDort[m][2]=subsData[i]->UnifDortSubg[j][1];
            //printf("Dortmund %i %i %i\n",unifacCompDort[m][0],unifacCompDort[m][1],unifacCompDort[m][2]);
            j++;
            m++;
        }
        j=0;
        if(subsData[i]->UnifNistSubg[j][0]<=0)useUnifacNist=0;
        while(subsData[i]->UnifNistSubg[j][0]>0){
            unifacCompNist[n][0]=i;
            unifacCompNist[n][1]=subsData[i]->UnifNistSubg[j][0];
            if(mixData.baseProp[i].MWmono>1) unifacCompNist[n][2]=mixData.baseProp[i].numMono*subsData[i]->UnifNistSubg[j][1];
            else unifacCompNist[n][2]=subsData[i]->UnifNistSubg[j][1];
            //printf("Nist %i %i %i\n",unifacCompNist[n][0],unifacCompNist[n][1],unifacCompNist[n][2]);
            j++;
            n++;
        }
    }
    //printf(" %i %i %i %i\n",k,unifacCompStd[k][0],unifacCompStd[k][1],unifacCompStd[k][2]);
    printf("Std:%i PSRK:%i Dort:%i Nist:%i\n",useUnifacStd,useUnifacPSRK,useUnifacDort,useUnifacNist==1);
    char path[FILENAME_MAX]="";
    if(useUnifacStd==1) FF_UNIFACParams(k,unifacCompStd,path,&mixData.unifStdData);
    if(useUnifacPSRK==1) FF_UNIFACParams(l,unifacCompPSRK,path,&mixData.unifPSRKData);
    if(useUnifacDort==1) FF_UNIFACParams(m,unifacCompDort,path,&mixData.unifDortData);
    if(useUnifacNist==1) FF_UNIFACParams(n,unifacCompNist,path,&mixData.unifNistData);
    return &mixData;

}

//Fill a Mixture data structure from an array of substance data structures
EXP_IMP void CALLCONV FF_MixFillDataWithSubsData(int numSubs,FF_SubstanceData subsData[], const char *path, FF_MixData *mixData){
    int ver=0;
    if (numSubs>15) numSubs=15;
    int h,i,j=0,k=0,l=0,m=0,n=0;
    int unifacCompStd[numSubs*10][3];//To contain the Unifac description for all substances, one after another:substance,subgroup,number
    int unifacCompPSRK[numSubs*10][3];
    int unifacCompDort[numSubs*10][3];
    int unifacCompNist[numSubs*10][3];
    int useUnifacStd=1,useUnifacPSRK=1,useUnifacDort=1,useUnifacNist=1;//Indicates that can be used in the mix. If not, we will put to 0

    for(i=0;i<numSubs*10;i++){
        unifacCompStd[i][0]=unifacCompStd[i][1]=unifacCompStd[i][2]=0;
    }

    for(h=0;h<numSubs*10;h++) for(i=0;i<3;i++) {
        unifacCompStd[h][i]=0;
        unifacCompPSRK[h][i]=0;
        unifacCompDort[h][i]=0;
        unifacCompNist[h][i]=0;
    }

    mixData->numSubs=numSubs;

    for (i=0;i<numSubs;i++){
        mixData->id[i]=subsData[i].id;
        strcpy(mixData->subsName[i],subsData[i].name);
        strcpy(mixData->CAS[i],subsData[i].CAS);
        mixData->baseProp[i]=subsData[i].baseProp;
        if((mixData->baseProp[i].MWmono>0)&&(mixData->baseProp[i].MW>0)) mixData->baseProp[i].numMono=557;//mixData->baseProp[i].numMono=mixData->baseProp[i].MW/mixData->baseProp[i].MWmono;
        else mixData->baseProp[i].numMono=1;
        //In FF_UnifacData we will have the FV of the polymer and in baseProp that of the monomer

        mixData->unifStdData.FV[i]=mixData->baseProp[i].numMono*subsData[i].baseProp.FV;
        mixData->unifPSRKData.FV[i]=mixData->baseProp[i].numMono*subsData[i].baseProp.FV;
        mixData->unifDortData.FV[i]=mixData->baseProp[i].numMono*subsData[i].baseProp.FV;
        mixData->unifNistData.FV[i]=mixData->baseProp[i].numMono*subsData[i].baseProp.FV;
        mixData->cubicData[i]=subsData[i].cubicData;
        mixData->saftData[i]=subsData[i].saftData;

        //Now preparation for using induced association
        if ((subsData[i].saftData.epsilonAB==0)&&(subsData[i].baseProp.mu>0.8)){
            mixData->saftData[i].nPos=1;
            mixData->saftData[i].nNeg=2;
            //printf("nPos:%i \n",mixData->saftData[i].nPos);
        }

        mixData->cp0[i]=subsData[i].cp0;
        mixData->vp[i]=subsData[i].vp;
        mixData->hVsat[i]=subsData[i].hVsat;
        mixData->lCp[i]=subsData[i].lCp;
        mixData->lDens[i]=subsData[i].lDens;
        mixData->lVisc[i]=subsData[i].lVisc;
        mixData->lThC[i]=subsData[i].lThC;
        mixData->lSurfT[i]=subsData[i].lSurfT;
        mixData->cp0Corr[i]=subsData[i].cp0Corr;
        mixData->vpCorr[i]=subsData[i].vpCorr;
        mixData->btCorr[i]=subsData[i].btCorr;
        mixData->hVsatCorr[i]=subsData[i].hVsatCorr;
        mixData->lCpCorr[i]=subsData[i].lCpCorr;
        mixData->lDensCorr[i]=subsData[i].lDensCorr;
        mixData->lViscCorr[i]=subsData[i].lViscCorr;
        mixData->lThCCorr[i]=subsData[i].lThCCorr;
        mixData->lSurfTCorr[i]=subsData[i].lSurfTCorr;
        mixData->gDensCorr[i]=subsData[i].gDensCorr;
        mixData->gViscCorr[i]=subsData[i].gViscCorr;
        mixData->gThCCorr[i]=subsData[i].gThCCorr;

        mixData->unifStdData.model=FF_UNIFACStd;
        mixData->unifPSRKData.model=FF_UNIFACPSRK;
        mixData->unifDortData.model=FF_UNIFACDort;
        mixData->unifNistData.model=FF_UNIFACNist;

        j=0;
        if(subsData[i].UnifStdSubg[j][1]<=0) useUnifacStd=0;//if one substance doesn't have this Unifac description no calculation can be made
        while(subsData[i].UnifStdSubg[j][1]>0){//Put the subgroups description of all substances in a single array of substance/subgroup/quantity
            unifacCompStd[k][0]=i;
            unifacCompStd[k][1]=subsData[i].UnifStdSubg[j][0];
            if(mixData->baseProp[i].MWmono>1) unifacCompStd[k][2]=mixData->baseProp[i].numMono*subsData[i].UnifStdSubg[j][1];
            else unifacCompStd[k][2]=subsData[i].UnifStdSubg[j][1];
            if (ver==1) printf("Temp Unifac Std: Substance:%i Subgroup:%i Quantity:%i\n",unifacCompStd[k][0],unifacCompStd[k][1],unifacCompStd[k][2]);
            j++;
            k++;
        }

        j=0;
        if(subsData[i].UnifPSRKSubg[j][1]<=0) useUnifacPSRK=0;
        while(subsData[i].UnifPSRKSubg[j][1]>0){
            unifacCompPSRK[l][0]=i;
            unifacCompPSRK[l][1]=subsData[i].UnifPSRKSubg[j][0];
            if(mixData->baseProp[i].MWmono>1) unifacCompPSRK[l][2]=mixData->baseProp[i].numMono*subsData[i].UnifPSRKSubg[j][1];
            else unifacCompPSRK[l][2]=subsData[i].UnifPSRKSubg[j][1];
            if (ver==1) printf("Temp Unifac PSRK: %i %i %i\n",unifacCompPSRK[l][0],unifacCompPSRK[l][1],unifacCompPSRK[l][2]);
            j++;
            l++;
        }

        j=0;
        if(subsData[i].UnifDortSubg[j][0]<=0) useUnifacDort=0;
        while(subsData[i].UnifDortSubg[j][0]>0){
            unifacCompDort[m][0]=i;
            unifacCompDort[m][1]=subsData[i].UnifDortSubg[j][0];
            if(mixData->baseProp[i].MWmono>1) unifacCompDort[m][2]=mixData->baseProp[i].numMono*subsData[i].UnifDortSubg[j][1];
            else unifacCompDort[m][2]=subsData[i].UnifDortSubg[j][1];
            if (ver==1) printf("Temp Unifac Dortmund: %i %i %i\n",unifacCompDort[m][0],unifacCompDort[m][1],unifacCompDort[m][2]);
            j++;
            m++;
        }

        j=0;
        if(subsData[i].UnifNistSubg[j][0]<=0)useUnifacNist=0;
        while(subsData[i].UnifNistSubg[j][0]>0){
            unifacCompNist[n][0]=i;
            unifacCompNist[n][1]=subsData[i].UnifNistSubg[j][0];
            if(mixData->baseProp[i].MWmono>1) unifacCompNist[n][2]=mixData->baseProp[i].numMono*subsData[i].UnifNistSubg[j][1];
            else unifacCompNist[n][2]=subsData[i].UnifNistSubg[j][1];
            if (ver==1) printf("Temp Unifac Nist: %i %i %i\n",unifacCompNist[n][0],unifacCompNist[n][1],unifacCompNist[n][2]);
            j++;
            n++;
        }

    }

    if(useUnifacStd==1)FF_UNIFACParams(k,unifacCompStd,path,&mixData->unifStdData);
    if(useUnifacPSRK==1)FF_UNIFACParams(l,unifacCompPSRK,path,&mixData->unifPSRKData);
    if(useUnifacDort==1)FF_UNIFACParams(m,unifacCompDort,path,&mixData->unifDortData);
    if(useUnifacNist==1)FF_UNIFACParams(n,unifacCompNist,path,&mixData->unifNistData);
    //return;
}

//Write a mixture data to a file. Adds ".md" extension
EXP_IMP void CALLCONV FF_MixDataToFile(const char *name,FF_MixData *mix){
    char path[FILENAME_MAX]="Data/";
    strcat(path,name);
    strcat(path,".md");
    FILE * file= fopen(path, "wb");
    if (file != NULL) {
        fwrite (mix, sizeof(FF_MixData), 1, file);
        fclose(file);
    }
    else printf("Error open file\n");
}
//Mixture Cubic EOS calculations
//==============================

//Calculates Theta,b,c and their composition derivatives, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//-----------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamXderCubicEOS(const int *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[15][15][6],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[]){
    int ver=0;
    int i,j;
    double k,kInv,l,a,aInv,b,Co0,Co;
    FF_CubicParam sParam[*numSubs];
    for (i=0;i< *numSubs;i++)//First we get the parameters of the individual substances
    {
        FF_FixedParamCubic(&data[i],&sParam[i]);
        FF_ThetaDerivCubic(T,&data[i],&sParam[i]);
        if (ver==1) printf(" i a,b,Theta,dTheta,d2Theta,x: %i %f %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta,sParam[i].dTheta,sParam[i].d2Theta,x[i]);
    }
    param->Theta=0;
    param->b=0;
    param->c=0;
    param->u=sParam[0].u;
    param->w=sParam[0].w;
    for (i=0;i< *numSubs;i++){
        if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))) return;//We check that all cubic eos are of the same type
        param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
        dc_dXi[i]=sParam[i].c;
    }
    switch (*rule)//six values are passed by every substances pair. Not all of them are used in every mixing rule.
    {
    case FF_VdW://van der Waals quadratic mixing rule
    case FF_VdWnoInt://If the BIP are for other model not use them
        for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
        {
            dTheta_dXi[i]=0;
            db_dXi[i]=0;
            for (j=0;j<*numSubs;j++)
            {
                //Three values are used for the calculation of the mixture parameter k[i,j]=k[j,i], the forth for l[i,j]=l[j,i], used in the calculation of b[i,j]
                if(*rule==FF_VdW){
                    k=pintParam[i][j][0] + pintParam[i][j][1]* *T + pintParam[i][j][2]/ *T ;
                    l=pintParam[i][j][3];
                }
                else{//If the BIP are for other model
                    k=0;
                    l=0;
                }

                a=fabs(pow(sParam[i].Theta*sParam[j].Theta,0.5))*(1-k);
                param->Theta=param->Theta+x[i]*x[j]*a;
                b=(sParam[i].b+sParam[j].b)/2*(1-l);
                param->b=param->b+x[i]*x[j]*b;
                dTheta_dXi[i]=dTheta_dXi[i]+x[j]*2*a;//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                db_dXi[i]=db_dXi[i]+x[j]*2*b;
            }
            //printf("dTheta_dXi:%f db_dXi:%f dc_dXi:%f\n",dTheta_dXi[i],db_dXi[i],dc_dXi[i]);
        }
        //printf("Theta:%f b:%f c:%f\n",param->Theta,param->b,param->c);
        break;
    case FF_PR://Panagiotopoulos and Reid composition dependent mixing rule
        for (i=0;i<*numSubs;i++) dTheta_dXi[i]=0;
        for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
        {
            //dTheta_dXi[i]=0;
            db_dXi[i]=0;
            for (j=0;j<*numSubs;j++)
            {
                Co=fabs(pow(sParam[i].Theta*sParam[j].Theta,0.5));
                k=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;
                kInv=pintParam[j][i][0] + pintParam[j][i][1]* *T+pintParam[j][i][2]/ *T ;
                a=Co*(1-k+(k-kInv)*x[i]);//this is theta[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5
                param->Theta=param->Theta+x[i]*x[j]*a;
                l=pintParam[i][j][3];
                b=(1-l)*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                param->b=param->b+x[i]*x[j]*b;
                aInv=Co*(1-kInv+(kInv-k)*x[j]);
                //dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a+x[i]*Co*(k-kInv)+aInv);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a+x[i]*Co*(k-kInv));
                dTheta_dXi[j]=dTheta_dXi[j]+x[i]*a;
                db_dXi[i]=db_dXi[i]+x[j]*2*b;
            }

        }
        break;
    case FF_MKP://In MKP th 3 first parameters are for k, forth and fith for lambda and the sixth for l (covolumen correction)
        {
        double kij,kji,lambda;
        double a2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
        double a3[*numSubs];//Is the summatory of second part of Theta calculation for each substance, and it derivative regarding T
        /**/
        for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
        {
            a3[i]=0;
            dTheta_dXi[i]=0;
            db_dXi[i]=0;
             for (j=0;j<*numSubs;j++)
            {
                //Three values are used for the calculation of the interaction parameters, one is k[i,j]=k[j,i], the second lambda[i,j]=-lambda[j,i], and l[i,j]=l[j,i]
                //First the part that is equal to VdW
                Co=fabs(pow(sParam[i].Theta*sParam[j].Theta,0.5));
                //printf("i,j k: %i %i %f\n",i,j,pintParam[i][j][0]);
                //Three values are used for the calculation of the mixture parameter k[i,j]=k[j,i], the forth for l[i,j]=l[j,i], used in the calculation of b[i,j]
                if(i!=j){
                    kij=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;//using same coef than Panagiotopoulos-Reid
                    kji=pintParam[j][i][0] + pintParam[j][i][1]* *T+pintParam[j][i][2]/ *T ;
                }
                else kij=kji=0;
                k=0.5*(kij+kji);
                //k=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;
                //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                a=Co*(1-k);
                param->Theta=param->Theta+x[i]*x[j]*a;
                dTheta_dXi[i]=dTheta_dXi[i]+x[j]*2*a;//This is the 1st part composition derivative (VdW) of Theta regarding x[i].
                //and taking into account that when calculate the term [j][i] we duplicate Theta and its derivatives
                //printf("a,da,d2a[i][j]: %i %i %f %f %f\n",i,j,a*1e3,da*1e3,d2a*1e3);
                l=pintParam[i][j][3];
                b=(1-l)*(sParam[i].b+sParam[j].b)/2;
                param->b=param->b+x[i]*x[j]*b;
                db_dXi[i]=db_dXi[i]+x[j]*2*b;

                //Now the additional part of Theta
                if(i!=j){
                    lambda=kji-kij;
                    //lambda=pintParam[j][i][3] + pintParam[j][i][4]* *T + pintParam[j][i][5]/ *T ;
                    a2[i][j]= x[j]*pow(fabs(Co*lambda),0.3333333);
                    if (lambda<0) a2[i][j]=-a2[i][j];
                }
                else{
                    a2[i][j]=0;
                }

                //printf("Co:%f lambda:%f a2:%f da2:%f d2a2:%f\n",Co,lambda,a2[i][j],da2[i][j],d2a2[i][j]);
                a3[i]=a3[i]+a2[i][j];
            }
            param->Theta=param->Theta+x[i]*pow(a3[i],3);//We increase Theta and its derivatives in the part special for MKP
        }
        for (i=0;i<*numSubs;i++)//And now is necessary to finish dTheta/dX[i]
        {
            for (j=0;j<*numSubs;j++)
            {
                dTheta_dXi[i]=dTheta_dXi[i]+3*x[j]*a3[j]*a3[j]*a2[j][i];//This is the derivative of the second part for the second component of the pair
            }
            dTheta_dXi[i]=dTheta_dXi[i]+pow(a3[i],3);//this is the derivative of the second part for the first component of the pair
        }

        break;
        }
    }
}


//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i], db/dX[i], dTheta_dXi for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamTderCubicEOS(const int *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[15][15][6],const double x[], FF_CubicParam *param)
{
    int i,j;
    FF_CubicParam sParam[*numSubs];
    double Pr,dPr,d2Pr;//Will be theta[i]*theta[j], and their derivatives
    double Co,dCo,d2Co;//This will be (theta[i]*theta[j])^0.5 and their derivatives regarding T
    double k,dk,d2k;//Will be the interaction coefficient for a[i][j], and their derivatives regarding T
    double kInv,dkInv,d2kInv,aInv,bInv;//Will be the corresponding [j][i] values
    double kCo,dkCo,d2kCo;//will be the combined interaction in PR and MKP rules
    double l;//the interaction coefficient for b[i][j]in VdW and PR rules, and for a2[i][j] for MKP rule
    double a,da,d2a;//Is the combined Theta[i,j] used in VdW and PR for the calculation of mixture Theta. It is also the first part of MKP Theta calculation
    double b;


    //double sumdTheta_dx=0,sumdb_dx=0;
    //First we get the parameters of the individual substances
    //later, calculation of the interaction coefficient and of its derivatives
    for (i=0;i< *numSubs;i++)//First we get the parameters of the individual substances
    {
        //printf("T: %f\n",*T);
        FF_FixedParamCubic(&data[i],&sParam[i]);
        FF_ThetaDerivCubic(T,&data[i],&sParam[i]);
        //printf("i:%i eos:%i Tc:%f Pc:%f w:%f k1:%f k2:%f k3:%f\n",i,data[i].eos,data[i].Tc,data[i].Pc,data[i].w,data[i].k1,data[i].k2,data[i].k3);
        //printf("i:%i a:%f b:%f Theta:%f dTheta*1e3:%f d2Theta*1e3:%f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
    }
    param->Theta=0;
    param->dTheta=0;
    param->d2Theta=0;
    param->b=0;
    param->c=0;
    param->u=sParam[0].u;
    param->w=sParam[0].w;
    for (i=0;i< *numSubs;i++){
        if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))) return;//We check that all cubic eos are of the same type
        param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
    }
    switch (*rule)//six values are passed by every substances pair. Not all of them are used in every mixing rule.
    {
        case FF_VdW://van der Waals quadratic mixing rule
        case FF_VdWnoInt:
            //printf("Using VdW mix rule\n");
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    Pr=sParam[i].Theta*sParam[j].Theta;
                    dPr=sParam[i].dTheta*sParam[j].Theta+sParam[i].Theta*sParam[j].dTheta;
                    d2Pr=sParam[i].d2Theta*sParam[j].Theta+2*sParam[i].dTheta*sParam[j].dTheta+sParam[i].Theta*sParam[j].d2Theta;
                    Co=fabs(pow(Pr,0.5));
                    dCo=0.5*pow(Pr,-0.5)*dPr;
                    d2Co=-0.25*pow(Pr,-1.5)*dPr*dPr+0.5*pow(Pr,-0.5)*d2Pr;
                    //printf("i,j k: %i %i %f\n",i,j,pintParam[i][j][0]);
                    //Three values are used for the calculation of the mixture parameter k[i,j]=k[j,i], the forth for l[i,j]=l[j,i], used in the calculation of b[i,j]
                    if(*rule==FF_VdW){
                        k=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;
                        dk=pintParam[i][j][1]-pintParam[i][j][2]/(*T * *T);
                        d2k=2* pintParam[i][j][2]/(*T * *T * *T);
                    }
                    else{
                        k=0;
                        dk=0;
                        d2k=0;
                    }
                    //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                    a=Co*(1-k);
                    da=dCo*(1-k)-Co*dk;
                    d2a=d2Co*(1-k)-2*dCo*dk-Co*d2k;
                    param->Theta=param->Theta+x[i]*x[j]*a;
                    param->dTheta=param->dTheta+x[i]*x[j]*da;
                    param->d2Theta=param->d2Theta+x[i]*x[j]*d2a;
                    //printf("a,da,d2a[i][j]: %i %i %f %f %f\n",i,j,a*1e3,da*1e3,d2a*1e3);
                    l=pintParam[i][j][3];
                    b=(sParam[i].b+sParam[j].b)/2*(1-l);//this is b[i,j]=(b[i]+b[j])/2*(1-l)
                    param->b=param->b+x[i]*x[j]*b;
                    //kInv=pintParam[(j* *numSubs+i)*6] + pintParam[(j* *numSubs+i)*6+1]* *T+pintParam[(j* *numSubs+i)*6+2]/ *T ;
                    //aInv=Co*(1-kInv);
                }
            }
            break;


        case FF_PR://Panagiotopoulos and Reid composition dependent mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    Pr=sParam[i].Theta*sParam[j].Theta;
                    dPr=sParam[i].dTheta*sParam[j].Theta+sParam[i].Theta*sParam[j].dTheta;
                    d2Pr=sParam[i].d2Theta*sParam[j].Theta+2*sParam[i].dTheta*sParam[j].dTheta+sParam[i].Theta*sParam[j].d2Theta;
                    Co=fabs(pow(Pr,0.5));
                    dCo=0.5*pow(Pr,-0.5)*dPr;
                    d2Co=-0.25*pow(Pr,-1.5)*dPr*dPr+0.5*pow(Pr,-0.5)*d2Pr;
                    k=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;
                    dk=pintParam[i][j][1]-pintParam[i][j][2]/(*T * *T);
                    d2k=2* pintParam[i][j][2]/(*T * *T * *T);
                    kInv=pintParam[j][i][0] + pintParam[j][i][1]* *T+pintParam[j][i][2]/ *T ;
                    dkInv=pintParam[j][i][1]-pintParam[j][i][2]/(*T * *T);
                    d2kInv=2* pintParam[j][i][2]/(*T * *T * *T);
                    //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                    kCo=1-k+(k-kInv)*x[i];
                    dkCo=-dk+(dk-dkInv)*x[i];
                    d2kCo=-d2k+(d2k-d2kInv)*x[i];
                    a=Co*kCo;//this is theta[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5
                    da=dCo*kCo+Co*dkCo;//this is the derivarive of a[i][j] regarding T
                    d2a=d2Co*kCo+2*dCo*dkCo+Co*d2kCo;//this is the second derivarive of a[i][j] regarding T
                    param->Theta=param->Theta+x[i]*x[j]*a;
                    param->dTheta=param->dTheta+x[i]*x[j]*da;
                    param->d2Theta=param->d2Theta+x[i]*x[j]*d2a;
                    l=pintParam[i][j][3];
                    b=(1-l)*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b;
                    //aInv=Co*(1-kInv+(kInv-k)*x[j]);
                }

            }
            break;

        case FF_MKP://Mathias, Klotz and Prausnitz composition dependent mixing rule
        {
            double kij,kji,lambda,dLambda,d2Lambda;
            double a2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double da2[*numSubs][*numSubs];//Its derivative
            double d2a2[*numSubs][*numSubs];//Its 2nd derivative
            double a3[*numSubs];//Is the summatory of second part of Theta calculation for each substance, and it derivative regarding T
            double da3[*numSubs];
            double d2a3[*numSubs];

            /**/
            param->dTheta=0;
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                a3[i]=0;
                da3[i]=0;
                d2a3[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j]=k[j,i], the second lambda[i,j]=-lambda[j,i], and l[i,j]=l[j,i]
                    //First the part that is equal to VdW
                    Pr=sParam[i].Theta*sParam[j].Theta;
                    dPr=sParam[i].dTheta*sParam[j].Theta+sParam[i].Theta*sParam[j].dTheta;
                    d2Pr=sParam[i].d2Theta*sParam[j].Theta+2*sParam[i].dTheta*sParam[j].dTheta+sParam[i].Theta*sParam[j].d2Theta;
                    Co=fabs(pow(Pr,0.5));
                    dCo=0.5*pow(Pr,-0.5)*dPr;
                    d2Co=-0.25*pow(Pr,-1.5)*dPr*dPr+0.5*pow(Pr,-0.5)*d2Pr;
                    //printf("i,j k: %i %i %f\n",i,j,pintParam[i][j][0]);
                    //Three values are used for the calculation of the mixture parameter k[i,j]=k[j,i], the forth for l[i,j]=l[j,i], used in the calculation of b[i,j]
                    kij=pintParam[i][j][0] + pintParam[i][j][1]* *T+pintParam[i][j][2]/ *T ;//using same coef than Panagiotopoulos-Reid
                    kji=pintParam[j][i][0] + pintParam[j][i][1]* *T+pintParam[j][i][2]/ *T ;
                    k=0.5*(kij+kji);
                    //k=pintParam[i][j][0] + pintParam[i][j][1]* *T + pintParam[i][j][2]/ *T ;
                    dk=0.5*(pintParam[i][j][1]+pintParam[j][i][1]-(pintParam[i][j][2]+pintParam[j][i][2])/(*T * *T));
                    d2k=(pintParam[i][j][2]+pintParam[i][j][2])/(*T * *T * *T);
                    //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                    a=Co*(1-k);
                    da=dCo*(1-k)-Co*dk;
                    d2a=d2Co*(1-k)-2*dCo*dk-Co*d2k;
                    param->Theta=param->Theta+x[i]*x[j]*a;
                    param->dTheta=param->dTheta+x[i]*x[j]*da;
                    param->d2Theta=param->d2Theta+x[i]*x[j]*d2a;
                    //printf("a,da,d2a[i][j]: %i %i %f %f %f\n",i,j,a*1e3,da*1e3,d2a*1e3);
                    l=pintParam[i][j][3];
                    b=(1-l)*(sParam[i].b+sParam[j].b)/2;
                    param->b=param->b+x[i]*x[j]*b;

                    //Now the additional part
                    if(i!=j){
                        lambda=kji-kij;
                        dLambda=pintParam[j][i][1]-pintParam[i][j][1]+(pintParam[i][j][2]-pintParam[j][i][2])/(*T * *T);
                        d2Lambda=2*(pintParam[j][i][2]-pintParam[i][j][2])/(*T * *T * *T);
                        //lambda=pintParam[j][i][3] + pintParam[j][i][4]* *T + pintParam[j][i][5]/ *T ;
                        //dLambda=pintParam[j][i][4] -pintParam[j][i][5]/(*T * *T);
                        //d2Lambda=2* pintParam[j][i][5]/(*T * *T * *T);
                        a2[i][j]= x[j]*pow(fabs(Co*lambda),(1.0/3.0));
                        //da2[i][j]=x[j]*(dCo*lambda+Co*dLambda)/(pow(fabs(Co*lambda),0.6666667)*3);
                        da2[i][j]=x[j]*pow(fabs(Co*lambda),(-2.0/3.0))*(dCo*lambda+Co*dLambda)/3;
                        d2a2[i][j]=-2*da2[i][j]*da2[i][j]/a2[i][j]+a2[i][j]*(d2Co*lambda+2*dCo*dLambda+Co*d2Lambda)/(3*Co*lambda);
                        if (lambda<0){
                            a2[i][j]=-a2[i][j];
                            //da2[i][j]=-da2[i][j];
                            //d2a2[i][j]=-d2a2[i][j];
                        }
                    }
                    else{
                        a2[i][j]=0;
                        da2[i][j]=0;
                        d2a2[i][j]=0;
                    }



                    //printf("Co:%f lambda:%f a2:%f da2:%f d2a2:%f\n",Co,lambda,a2[i][j],da2[i][j],d2a2[i][j]);
                    a3[i]=a3[i]+a2[i][j];
                    da3[i]=da3[i]+da2[i][j];
                    d2a3[i]=d2a3[i]+d2a2[i][j];

                }

                param->Theta=param->Theta+x[i]*pow(a3[i],3);//We increase Theta and its derivatives in the part special for MKP
                param->dTheta=param->dTheta+x[i]*3*pow(a3[i],2)*da3[i];
                param->d2Theta=param->d2Theta+x[i]*3*(2*a3[i]*da3[i]*da3[i]+a3[i]*a3[i]*d2a3[i]);
                //printf("a3:%f da3:%f d2a3:%f Theta:%f dTheta_dT:%f d2Theta_dT2:%f\n",a3[i],da3[i],d2a3[i],param->Theta,param->dTheta,param->d2Theta);
            }
            break;
        }
        default:
            param->Theta=0;
            param->b=0;
            param->c=0;
            break;
    }
    //printf("Theta:%f dTheta:%f d2Theta:%f b:%f\n",param->Theta,param->dTheta,param->d2Theta,param->b);
}

//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i] and db/dX[i] for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamCubicEOSOld(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[])
{
    int i,j;
     FF_CubicParam sParam[*numSubs];
    double rootTij;//This will be (theta[i]*theta[j])^0.5
    double dRootTij;//Will be the derivative of (theta[i]*theta[j])^0.5 regarding T
    double k[*numSubs][*numSubs];//the interaction coefficient for a[i][j]
    double dk[*numSubs][*numSubs];//Derivative of the interaction coefficient for a[i][j] regarding T
    double d2k[*numSubs][*numSubs];//Second derivative of the interaction coefficient for a[i][j] regarding T
    double l[*numSubs][*numSubs];//the interaction coefficient for b[i][j]in VdW and PR rules, and for b[i][j]
    double dl[*numSubs][*numSubs];//Derivative of the interaction coefficient for b[i][j]in VdW and PR rules, and of a2[i][j] for MKP rule
    double a[*numSubs][*numSubs];//Is the combined Theta[i,j] used in VdW and PR for the calculation of mixture Theta. It is also the first part of MKP Theta calculation
    double da[*numSubs][*numSubs];//Is the derivative of Theta[i,j] regarding T
    double b[*numSubs][*numSubs];

    //double sumdTheta_dx=0,sumdb_dx=0;
    //First we get the parameters of the individual substances
    //later, calculation of the interaction coefficient and of its derivatives
    for (i=0;i< *numSubs;i++)//First we get the parameters of the individual substances
    {
        FF_FixedParamCubic(&data[i],&sParam[i]);
        FF_ThetaDerivCubic(T,&data[i],&sParam[i]);
        for (j=0;j<*numSubs;j++){
            k[i][j]=pintParam[(i* *numSubs+j)*6] + pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T ;
            dk[i][j]=pintParam[(i* *numSubs+j)*6+1]-pintParam[(i* *numSubs+j)*6+2]/(*T * *T);
            d2k[i][j]=2* pintParam[(i* *numSubs+j)*6+2]/(*T * *T * *T);
            l[i][j]=pintParam[(i* *numSubs+j)*6+3] + pintParam[(i* *numSubs+j)*6+4]* *T+pintParam[(i* *numSubs+j)*6+5]/ *T;
            dl[i][j]=pintParam[(i* *numSubs+j)*6+4]-pintParam[(i* *numSubs+j)*6+5]/(*T * *T);
        }
    }
    param->Theta=0;
    param->dTheta=0;
    param->d2Theta=0;
    param->b=0;
    param->c=0;
    param->u=sParam[0].u;
    param->w=sParam[0].w;
    for (i=0;i< *numSubs;i++){
        param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
        dc_dXi[i]=sParam[i].c;
    }
    switch (*rule)//six values are passed by every substances pair. Not all of them are used in every mixing rule.
    {
        case FF_VdW://van der Waals quadratic mixing rule
        case FF_VdWnoInt:
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    //Two values are used for the calculation of the mixture parameters, one for k[i,j]=k[j,i], the other for l[i,j]=l[j,i], used in the calculation of b[i,j]
                    rootTij=fabs(pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    dRootTij=0.5*rootTij*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    a[i][j]=(1-k[i][j])*rootTij;//this is theta[i,j]=(1-k[i,j])*(theta[i]*theta[j])^0.5
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];
                    //da[i][j]=dRootTij*(1-k[i][j])-rootTij*dk[i][j];//this is the derivarive of a[i][j] regarding T
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta)-dk[i][j]*rootTij;
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(0.5*da[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta)+0.5*a[i][j]*
                                   ((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*
                                   sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta)-d2k[i][j]*rootTij-dk[i][j]*dRootTij);
                    //param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                                //sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                }
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+a[j][i]);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j] and , al
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
            }
            break;
        case FF_PR://Panagiotopoulos and Reid composition dependent mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j] not equal to k[j,i], the other is l[i,j]=l[j,i]
                    a[i][j]=fabs((1-*(pintParam+(i* *numSubs+j)*6)+x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6)))*pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    //this is theta[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);//this is the derivarive of a[i][j] regarding T
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                                sParam[i].Theta/sParam[i].Theta)+((sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(1-*(pintParam+(i* *numSubs+j)*6+3))*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                }
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+pow(sParam[i].Theta*sParam[j].Theta,0.5)*x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6))+
                                                      a[j][i]);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
            }
            break;
        case FF_MKP://Mathias, Klotz and Prausnitz composition dependent mixing rule
        {
            double a2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double da2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double a3[*numSubs];//Is the summatory of second part of Theta calculation for each substance, and it derivative regarding T
            double da3[*numSubs];//The same for second part in MKP
            double d2a3[*numSubs];//The same for second part in MKP
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                a3[i]=0;
                da3[i]=0;
                d2a3[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j]=k[j,i], the second lambda[i,j]=-lambda[j,i], and l[i,j]=l[j,i]
                    a[i][j]=fabs((1-*(pintParam+(i* *numSubs+j)*6)+x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6)))*pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    //this is Theta1[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5. It is the first part of Theta calculation
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];//We increase theta in the part that is as in VdW mixing rule
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                            sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(1-*(pintParam+(i* *numSubs+j)*6+2))*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                    a2[i][j]= x[j]*pow(sParam[i].Theta*sParam[j].Theta,0.1666667)* pow(fabs(*(pintParam+(j* *numSubs+i)*6+1)),0.333333);
                    if (*(pintParam+(j* *numSubs+i)*6+1)<0) a2[i][j]=-a2[i][j];
                    a3[i]=a3[i]+a2[i][j];
                    //printf("%f ",*(pintParam+(j* *numSubs+i)*3+1));
                    //printf("a2[i,j]:%f\n",a2[i][j]);
                    da2[i][j]=a2[i][j]/6*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    da3[i]=da3[i]+da2[i][j];
                    if (i!=i)
                    d2a3[i]=d2a3[i]+da2[i][j]*da2[i][j]/a2[i][j]+a2[i][j]*(((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                        sParam[i].Theta/sParam[i].Theta)+((sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta))/6;
                    //printf("d2a3[i]%f\n",d2a3[i]);
                }
                param->Theta=param->Theta+x[i]*pow(a3[i],3);//We increase Theta and its derivatives in the part special for MKP
                param->dTheta=param->dTheta+x[i]*3*pow(a3[i],2)*da3[i];
                param->d2Theta=param->d2Theta+x[i]*3*(2*a3[i]*da3[i]*da3[i]+a3[i]*a3[i]*d2a3[i]);
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+a[j][i]);//This is the first part (VdW) of the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j]
                    dTheta_dXi[i]=dTheta_dXi[i]+3*x[j]*a3[j]*a3[j]*a2[j][i]/x[i];//This is the derivative for the second part when [i] is not the first component of the pair
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
                dTheta_dXi[i]=dTheta_dXi[i]+pow(a3[i],3);//this is the derivative for the second part when [i] is the first component of the pair
            }
            break;
        }
        default:
            param->Theta=0;
            param->b=0;
            param->c=0;
            break;
    }
    //printf("Theta:%f dTheta:%f d2Theta:%f b:%f\n",param->Theta,param->dTheta,param->d2Theta,param->b);
}


//Calculates Theta,b and c, given a cubic EOS,excess G, composition, and pure substance parameters
//------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamCubicEOSgE(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param){
        int ver=0;
        int i,j;
        FF_SubsActivityData actData[mix->numSubs];
        FF_ExcessData excData={0,0,0,0};
        FF_CubicParam sParam[mix->numSubs];
        double RT=R* *T;//to speed up
        double bij;//used to calculate the combined b for substances i,j

        //First we need to get gE. The data needed is inside the mix structure
        switch(mix->actModel){
        case FF_UNIFACStd:
            FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
            break;
        case FF_UNIFACPSRK:
            FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
            break;
        case FF_UNIFACDort:
            FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
            break;
        case FF_UNIFACNist:
            FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
            break;
        default:
            FF_Activity(mix,T,x,actData);
            break;
        }
        //Now we obtain excess g
        for(i=0;i<mix->numSubs;i++){
            excData.gEC=excData.gEC+x[i]*actData[i].lnGammaC;
            excData.gESG=excData.gESG+x[i]*actData[i].lnGammaSG;
            excData.gER=excData.gER+x[i]*actData[i].lnGammaR;
        }
        excData.gE=excData.gEC+excData.gESG+excData.gER;
        if (ver==1) printf("gE: %f %f %f %f\n",excData.gEC,excData.gESG,excData.gER,excData.gE);

        //Next we get the parameters of the individual substances
        for (i=0;i< mix->numSubs;i++)//First we get the parameters of the individual substances
        {
            FF_FixedParamCubic(&mix->cubicData[i],&sParam[i]);
            FF_ThetaDerivCubic(T,&mix->cubicData[i],&sParam[i]);
            if (ver==1) printf("i a,b,Theta,dTheta,d2Theta: %i %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta*1e3,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
        }
        param->Theta=0;
        param->b=0;
        param->c=0;
        param->u=sParam[0].u;
        param->w=sParam[0].w;
        for (i=0;i< mix->numSubs;i++){
            if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))) {
                printf("Different type of cubic EOS\n");
                return;//We check that all cubic eos are of the same type
            }
            if (sParam[i].c>0) param->c=param->c+x[i]*sParam[i].c;//Calculation of mix volume translation
        }

        double sPartA=0,sPartB=0;//sum of individual substances contribution
        double q1=0,q2=0,lambda=0;//Parameters for  alpha (and Theta) calculation, lambda is for LCVM mixing rule

        switch (mix->mixRule)
        {
        case FF_HV:
            if (mix->eosType==FF_CubicPRtype) q1=-0.623;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.693;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gESG+excData.gER)/q1+sPartB);
            }
            break;
        case FF_MHV1:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
            }
            break;
        case FF_PSRK:
            if (mix->eosType==FF_CubicSRKtype){
                q1=-0.64663;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
            }
            break;
        case FF_LCVM:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.623;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.693;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if ((mix->actModel==FF_UNIFACStd)||(mix->actModel==FF_UNIFACPSRK)) lambda=0.36;
            else if ((mix->actModel==FF_UNIFACDort)||(mix->actModel==FF_UNIFACNist)) lambda=0.65;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(lambda*(RT*excData.gER/q2+sPartB)+(1-lambda)*(RT*(excData.gE+sPartA)/q1+sPartB));
            }
            break;
        case FF_MHV2:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.4347;
                q2=-0.003654;
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.4783;
                q2=-0.0047;
            }
            if (!(q2==0)){
                double alphai;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;;
                }
                for (i=0;i<mix->numSubs;i++){
                    alphai=sParam[i].Theta/(sParam[i].b*RT);
                    sPartB=sPartB+x[i]*(log(param->b/sParam[i].b)+q1*alphai+q2*pow(alphai,2));
                }
                param->Theta=param->b*RT*(-q1-pow(q1*q1+4*q2*(sPartB+excData.gE),0.5))/(2*q2);
            }
            break;
        case FF_UMR:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for(i=0;i<mix->numSubs;i++){
                    for(j=0;j<mix->numSubs;j++){
                        bij=pow((pow(sParam[i].b,0.5)+pow(sParam[j].b,0.5))/2,2);
                        param->b=param->b+x[i]*x[j]*bij;
                    }
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gER)/q1+sPartB);
            }
            break;
        default:
            param->c=0;
            param->b=0;
            param->Theta=0;
            break;
        }

    if (ver==1) printf("Theta:%f b:%f c:%f\n",param->Theta,param->b,param->c);
}

//Calculates Theta,b,c and their composition derivatives, given a cubic EOS,excess G, composition, and pure substance parameters
//------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamNderCubicEOSgE(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param,double dNalpha_dNi[],double dNb_dNi[]){
        int i,j;
        FF_SubsActivityData actData[mix->numSubs];
        FF_ExcessData excData={0,0,0,0};
        FF_CubicParam sParam[mix->numSubs];
        double sumdb_dXi;
        double RT=R* *T;//to speed up
        double bij;//used to calculate the combined b for substances i,j
        double alpha,alpha2,alphai;

        //First we need to get gE. The data needed is inside the mix structure
        switch(mix->actModel){
        case FF_UNIFACStd:
            FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
            break;
        case FF_UNIFACPSRK:
            FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
            break;
        case FF_UNIFACDort:
            FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
            break;
        case FF_UNIFACNist:
            FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
            break;
        default:
            FF_Activity(mix,T,x,actData);
            break;
        }
        //Now we obtain excess g
        for(i=0;i<mix->numSubs;i++){
            excData.gEC=excData.gEC+x[i]*actData[i].lnGammaC;
            excData.gESG=excData.gESG+x[i]*actData[i].lnGammaSG;
            excData.gER=excData.gER+x[i]*actData[i].lnGammaR;
        }
        excData.gE=excData.gEC+excData.gESG+excData.gER;
        //printf("gE: %f %f %f %f\n",excData.gEC,excData.gESG,excData.gER,excData.gE);

        //Next we get the parameters of the individual substances
        for (i=0;i< mix->numSubs;i++)//First we get the parameters of the individual substances
        {
            FF_FixedParamCubic(&mix->cubicData[i],&sParam[i]);
            FF_ThetaDerivCubic(T,&mix->cubicData[i],&sParam[i]);
            //printf("i a,b,Theta,dTheta,d2Theta: %i %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta*1e3,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
        }
        param->Theta=0;
        param->b=0;
        param->c=0;
        param->u=sParam[0].u;
        param->w=sParam[0].w;
        for (i=0;i< mix->numSubs;i++){
            if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))){
                printf("Cubic EOS are of different types\n");
                return;//We check that all cubic eos are of the same type
            }
        }

        double sPartA=0,sPartB=0,dPartB_dXi[mix->numSubs];//sum of individual substances contribution
        double q1=0,q2=0,lambda=0;//Parameters for  alpha (and Theta) calculation, lambda is for LCVM mixing rule
        double gE=0;

        switch (mix->mixRule)
        {
        case FF_HV:
            if (mix->eosType==FF_CubicPRtype) q1=-0.623;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.693;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    dNb_dNi[i]=sParam[i].b;
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gESG+excData.gER)/q1+sPartB);
                for (i=0;i<mix->numSubs;i++){
                    dNalpha_dNi[i]=sParam[i].Theta/(sParam[i].b*RT)+(actData[i].lnGammaSG+actData[i].lnGammaR)/q1;
                }
            }
            break;
        case FF_MHV1:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    dNb_dNi[i]=sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                for (i=0;i<mix->numSubs;i++){
                    dNalpha_dNi[i]=sParam[i].Theta/(sParam[i].b*RT)+(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR+log(param->b/sParam[i].b)+
                                      sParam[i].b/param->b-1)/q1;
                }
            }
            break;
        case FF_PSRK:
            if (mix->eosType==FF_CubicSRKtype){
                q1=-0.64663;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    dNb_dNi[i]=sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                for (i=0;i<mix->numSubs;i++){
                    dNalpha_dNi[i]=sParam[i].Theta/(sParam[i].b*RT)+(actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR+log(param->b/sParam[i].b)+
                                      sParam[i].b/param->b-1)/q1;
                }
            }
            break;
        case FF_LCVM:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.623;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.693;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if ((mix->actModel==FF_UNIFACStd)||(mix->actModel==FF_UNIFACPSRK)) lambda=0.36;
            else if ((mix->actModel==FF_UNIFACDort)||(mix->actModel==FF_UNIFACNist)) lambda=0.65;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    dNb_dNi[i]=sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                }
                param->Theta=param->b*(lambda*(RT*excData.gER/q2+sPartB)+(1-lambda)*(RT*(excData.gE+sPartA)/q1+sPartB));
                for (i=0;i<mix->numSubs;i++){
                    dNalpha_dNi[i]=sParam[i].Theta/(sParam[i].b*RT)+lambda*(actData[i].lnGammaSG+actData[i].lnGammaR)/q2+(1-lambda)*
                            (actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR+log(param->b/sParam[i].b)+sParam[i].b/param->b-1)/q1;
                }
            }
            break;
        case FF_MHV2:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.4347;
                q2=-0.003654;
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.4783;
                q2=-0.0047;
            }
            if (!(q2==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    dNb_dNi[i]=sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    alphai=sParam[i].Theta/(sParam[i].b*RT);
                    sPartB=sPartB+x[i]*(log(param->b/sParam[i].b)+q1*alphai+q2*pow(alphai,2));
                }
                param->Theta=param->b*RT*(-q1-pow(q1*q1+4*q2*(sPartB+excData.gE),0.5))/(2*q2);
                alpha=param->Theta/(RT*param->b);
                alpha2=alpha*alpha;
                for (i=0;i<mix->numSubs;i++){
                    alphai=sParam[i].Theta/(sParam[i].b*RT);
                    dNalpha_dNi[i]=(q1*alphai+actData[i].lnGammaC+actData[i].lnGammaSG+actData[i].lnGammaR+q2*(alphai*alphai+alpha2)+log(param->b/sParam[i].b)+
                                    sParam[i].b/param->b-1)/(q1+2*alpha*q2);
                }
            }
            break;
        case FF_UMR:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                sumdb_dXi=0;
                for(i=0;i<mix->numSubs;i++){
                    dNb_dNi[i]=0;
                    for(j=0;j<mix->numSubs;j++){
                        bij=pow((pow(sParam[i].b,0.5)+pow(sParam[j].b,0.5))/2,2);
                        param->b=param->b+x[i]*x[j]*bij;
                        dNb_dNi[i]=dNb_dNi[i]+2*x[j]*bij;//this is db_dXi[i]
                    }
                    sumdb_dXi=sumdb_dXi+x[i]*dNb_dNi[i];
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dPartB_dXi[i]=sParam[i].Theta/sParam[i].b;
                }
                param->Theta=param->b*(RT*(excData.gER)/q1+sPartB);
                for (i=0;i<mix->numSubs;i++){
                    dNb_dNi[i]=dNb_dNi[i]+param->b-sumdb_dXi;
                    dNalpha_dNi[i]=sParam[i].Theta/(sParam[i].b*RT)+(actData[i].lnGammaR+log(param->b/sParam[i].b)+
                                      sParam[i].b/param->b-1)/q1;
                }
            }
            break;
        default:
            param->c=0;
            param->b=0;
            param->Theta=0;
            break;
        }

    //printf("Theta:%f b:%f c:%f\n",param->Theta,param->b,param->c);
}


//Calculates Theta,b,c and their temperature derivatives, given a cubic EOS,excess G, composition, and pure substance parameters
//------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamTderCubicEOSgE(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param){
    int ver=0;
    int i,j;
    double dT=0.01;
    double Tplus,Tminus;//Temperature variation to obtain numeric temperature derivative of theta
    Tplus=*T + dT;
    Tminus=*T-dT;
    FF_SubsActivityData actData[mix->numSubs],actDataPlus[mix->numSubs],actDataMinus[mix->numSubs];
    //FF_SubsActivityData actDataPlus[mix->numSubs];
    FF_ExcessData excData={0,0,0,0},excDataPlus={0,0,0,0},excDataMinus={0,0,0,0};
    //FF_ExcessData excDataPlus={0,0,0,0};
    FF_CubicParam sParam[mix->numSubs];
    double RT=R* *T;//to speed up
    double bij;//used to calculate the combined b for substances i,j
    double dgEC,dgESG,dgER,dgE;//derivatives og gE regarding T
    double d2gEC,d2gESG,d2gER,d2gE;//second derivatives og gE regarding T
    //First we need to get activity. The data needed is inside the mix structure
    //printf("Act model: %i\n",mix->actModel);
    switch(mix->actModel){
    case FF_UNIFACStd:
        FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifStdData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifStdData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACPSRK:
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifPSRKData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifPSRKData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACDort:
        FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifDortData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifDortData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACNist:
        FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifNistData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifNistData,&Tminus,x,actDataMinus);
        break;
    default:
        FF_Activity(mix,T,x,actData);
        FF_Activity(mix,&Tplus,x,actDataPlus);
        FF_Activity(mix,&Tminus,x,actDataMinus);
        break;
    }
        //Now we obtain excess gE
        for(i=0;i<mix->numSubs;i++){
            excData.gEC=excData.gEC+x[i]*actData[i].lnGammaC;
            excData.gESG=excData.gESG+x[i]*actData[i].lnGammaSG;
            excData.gER=excData.gER+x[i]*actData[i].lnGammaR;
            excDataPlus.gEC=excDataPlus.gEC+x[i]*actDataPlus[i].lnGammaC;
            excDataPlus.gESG=excDataPlus.gESG+x[i]*actDataPlus[i].lnGammaSG;
            excDataPlus.gER=excDataPlus.gER+x[i]*actDataPlus[i].lnGammaR;
            excDataMinus.gEC=excDataMinus.gEC+x[i]*actDataMinus[i].lnGammaC;
            excDataMinus.gESG=excDataMinus.gESG+x[i]*actDataMinus[i].lnGammaSG;
            excDataMinus.gER=excDataMinus.gER+x[i]*actDataMinus[i].lnGammaR;
        }
        excData.gE=excData.gEC+excData.gESG+excData.gER;
        excDataPlus.gE=excDataPlus.gEC+excDataPlus.gESG+excDataPlus.gER;
        excDataMinus.gE=excDataMinus.gEC+excDataMinus.gESG+excDataMinus.gER;
        dgEC=(excDataPlus.gEC-excData.gEC)/dT;
        d2gEC=(dgEC-((excData.gEC-excDataMinus.gEC)/dT))/dT;
        dgESG=(excDataPlus.gESG-excData.gESG)/dT;
        d2gESG=(dgESG-((excData.gESG-excDataMinus.gESG)/dT))/dT;
        dgER=(excDataPlus.gER-excData.gER)/dT;
        d2gER=(dgER-((excData.gER-excDataMinus.gER)/dT))/dT;
        dgE=dgEC+dgESG+dgER;
        d2gE=d2gEC+d2gESG+d2gER;
        if (ver==1) printf("gE: %f %f %f %f\n",excData.gEC,excData.gESG,excData.gER,excData.gE);

        //Next we get the parameters of the individual substances
        for (i=0;i< mix->numSubs;i++)//First we get the parameters of the individual substances
        {
            FF_FixedParamCubic(&mix->cubicData[i],&sParam[i]);
            FF_ThetaDerivCubic(T,&mix->cubicData[i],&sParam[i]);
            if (ver==1) printf("i a,b,Theta,dTheta,d2Theta: %i %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta*1e3,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
        }
        param->Theta=0;
        param->b=0;
        param->c=0;
        param->u=sParam[0].u;
        param->w=sParam[0].w;
        for (i=0;i< mix->numSubs;i++){
            if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))){
                printf("Cubic EOS are of different types\n");
                return;//We check that all cubic eos are of the same type
            }
            if (sParam[i].c>0) param->c=param->c+x[i]*sParam[i].c;//Calculation of mix volume translation
        }
        double sPartA=0,sPartB=0;//sum of individual substances contribution
        double dsPartB=0,d2sPartB=0;//Its derivatives regarding T
        double q1=0,q2=0,lambda=0;//Parameters for  alpha (and Theta) calculation, lambda is for LCVM mixing rule

        switch (mix->mixRule)
        {
        case FF_HV:
            if (mix->eosType==FF_CubicPRtype) q1=-0.623;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.693;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gESG+excData.gER)/q1+sPartB);
                param->dTheta=param->b*(R*(excData.gESG+excData.gER)/q1+RT*(dgESG+dgER)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*(dgESG+dgER)/q1+RT*(d2gESG+d2gER)/q1+d2sPartB);
            }
            break;
        case FF_MHV1:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                //printf("gE:%f Theta:%f\n",excData.gE,param->Theta);
                param->dTheta=param->b*(R*(excData.gE+sPartA)/q1+RT*(dgE)/q1+dsPartB);
                //printf("dTheta:%f\n",param->dTheta);
                param->d2Theta=param->b*(2*R*dgE/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        case FF_PSRK:
            if (mix->eosType==FF_CubicSRKtype){
                q1=-0.64663;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                param->dTheta=param->b*(R*(excData.gE+sPartA)/q1+RT*(dgE)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*dgE/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        case FF_LCVM:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.623;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.693;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if ((mix->actModel==FF_UNIFACStd)||(mix->actModel==FF_UNIFACPSRK)) lambda=0.36;
            else if ((mix->actModel==FF_UNIFACDort)||(mix->actModel==FF_UNIFACNist)) lambda=0.65;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                //param->Theta=param->b*(lambda*(RT*excData.gER/q2+sPartB)+(1-lambda)*(RT*(excData.gE+sPartA)/q1+sPartB));
                param->Theta=param->b*(lambda*(RT*excData.gER/q2)+(1-lambda)*(RT*(excData.gE+sPartA)/q1)+sPartB);
                param->dTheta=param->b*(lambda*(R*excData.gER+RT*dgER)/q2+(1-lambda)*(R*(excData.gE+sPartA)+RT*dgE)/q1+dsPartB);
                param->d2Theta=param->b*(lambda*(2*R*dgER+RT*d2gER)/q2+(1-lambda)*(2*R*dgE+RT*d2gE)/q1+d2sPartB);
            }
            break;
        case FF_MHV2:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.4347;
                q2=-0.003654;
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.4783;
                q2=-0.0047;
            }
            if (!(q2==0)){
                double alphai,dalphai,d2alphai;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    alphai=sParam[i].Theta/(sParam[i].b*RT);
                    dalphai=sParam[i].dTheta/(sParam[i].b*RT)-alphai/ *T;
                    //d2alphai=sParam[i].d2Theta/(sParam[i].b*RT)-dalphai/ *T-(dalphai* *T-alphai)/(*T * *T);
                    d2alphai=sParam[i].d2Theta/(sParam[i].b*RT)-2*sParam[i].dTheta/(sParam[i].b*RT* *T)+2*sParam[i].Theta/(sParam[i].b*RT* *T * *T);
                    sPartB=sPartB+x[i]*(log(param->b/sParam[i].b)+q1*alphai+q2*pow(alphai,2));
                    dsPartB=dsPartB+x[i]*(q1*dalphai+2*q2*alphai*dalphai);
                    d2sPartB=d2sPartB+x[i]*(q1*d2alphai+2*q2*(dalphai*dalphai+alphai*d2alphai));

                }
                double aux=q1*q1+4*q2*(sPartB+excData.gE);
                double daux=4*q2*(dsPartB+dgE);
                param->Theta=param->b*RT*(-q1-pow(aux,0.5))/(2*q2);
                //printf("gE:%f Theta:%f\n",excData.gE,param->Theta);
                param->dTheta=param->Theta/ *T-param->b*RT*pow(aux,-0.5)*(dsPartB+dgE);
                //printf("dTheta:%f\n",param->dTheta);
                param->d2Theta=(param->dTheta* *T-param->Theta)/(*T * *T)-param->b*R*(pow(aux,-0.5)*(dsPartB+dgE)+
                               *T*(-0.5*pow(aux,-1.5)*daux*(dsPartB+dgE)+pow(aux,-0.5)*(d2sPartB+d2gE)));
            }
            break;
        case FF_UMR:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for(i=0;i<mix->numSubs;i++){
                    for(j=0;j<mix->numSubs;j++){
                        bij=pow((pow(sParam[i].b,0.5)+pow(sParam[j].b,0.5))/2,2);
                        param->b=param->b+x[i]*x[j]*bij;
                    }
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gER)/q1+sPartB);
                param->dTheta=param->b*(R*(excData.gER)/q1+RT*(dgER)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*dgER/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        default:
            param->c=0;
            param->b=0;
            param->Theta=0;
            param->dTheta=0;
            param->d2Theta=0;
            break;
        }
    if (ver==1) printf("Theta:%f dTheta:%f d2Theta:%f b:%f c:%f\n",param->Theta,param->dTheta, param->d2Theta, param->b,param->c);
}


//Calculates Theta,b,c and their temperature derivatives, given a cubic EOS,excess G, composition, and pure substance parameters
//------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamTderCubicEOSgE2(const FF_MixData *mix,const double *T,const double x[], FF_CubicParam *param){
    int ver=0;
    int i,j;
    double dT=0.01;
    double Tplus,Tminus,ThetaPlus,ThetaMinus;//Temperature variation to obtain numeric temperature derivative of theta
    Tplus=*T + dT;
    Tminus=*T-dT;
    FF_SubsActivityData actData[mix->numSubs],actDataPlus[mix->numSubs],actDataMinus[mix->numSubs];
    //FF_SubsActivityData actDataPlus[mix->numSubs];
    FF_ExcessData excData={0,0,0,0},excDataPlus={0,0,0,0},excDataMinus={0,0,0,0};
    //FF_ExcessData excDataPlus={0,0,0,0};
    FF_CubicParam sParam[mix->numSubs],sParamPlus[mix->numSubs],sParamMinus[mix->numSubs];
    double RT=R* *T;//to speed up
    double bij;//used to calculate the combined b for substances i,j
    double dgEC,dgESG,dgER,dgE;//derivatives og gE regarding T
    double d2gEC,d2gESG,d2gER,d2gE;//second derivatives og gE regarding T
    //First we need to get activity. The data needed is inside the mix structure
    //printf("Act model: %i\n",mix->actModel);
    switch(mix->actModel){
    case FF_UNIFACStd:
        FF_ActivityUNIFAC(&mix->unifStdData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifStdData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifStdData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACPSRK:
        FF_ActivityUNIFAC(&mix->unifPSRKData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifPSRKData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifPSRKData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACDort:
        FF_ActivityUNIFAC(&mix->unifDortData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifDortData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifDortData,&Tminus,x,actDataMinus);
        break;
    case FF_UNIFACNist:
        FF_ActivityUNIFAC(&mix->unifNistData,T,x,actData);
        FF_ActivityUNIFAC(&mix->unifNistData,&Tplus,x,actDataPlus);
        FF_ActivityUNIFAC(&mix->unifNistData,&Tminus,x,actDataMinus);
        break;
    default:
        FF_Activity(mix,T,x,actData);
        FF_Activity(mix,&Tplus,x,actDataPlus);
        FF_Activity(mix,&Tminus,x,actDataMinus);
        break;
    }
        //Now we obtain excess gE
        for(i=0;i<mix->numSubs;i++){
            excData.gEC=excData.gEC+x[i]*actData[i].lnGammaC;
            excData.gESG=excData.gESG+x[i]*actData[i].lnGammaSG;
            excData.gER=excData.gER+x[i]*actData[i].lnGammaR;
            excDataPlus.gEC=excDataPlus.gEC+x[i]*actDataPlus[i].lnGammaC;
            excDataPlus.gESG=excDataPlus.gESG+x[i]*actDataPlus[i].lnGammaSG;
            excDataPlus.gER=excDataPlus.gER+x[i]*actDataPlus[i].lnGammaR;
            excDataMinus.gEC=excDataMinus.gEC+x[i]*actDataMinus[i].lnGammaC;
            excDataMinus.gESG=excDataMinus.gESG+x[i]*actDataMinus[i].lnGammaSG;
            excDataMinus.gER=excDataMinus.gER+x[i]*actDataMinus[i].lnGammaR;
        }
        excData.gE=excData.gEC+excData.gESG+excData.gER;
        excDataPlus.gE=excDataPlus.gEC+excDataPlus.gESG+excDataPlus.gER;
        excDataMinus.gE=excDataMinus.gEC+excDataMinus.gESG+excDataMinus.gER;
        /*dgEC=(excDataPlus.gEC-excData.gEC)/dT;
        d2gEC=(dgEC-((excData.gEC-excDataMinus.gEC)/dT))/dT;
        dgESG=(excDataPlus.gESG-excData.gESG)/dT;
        d2gESG=(dgESG-((excData.gESG-excDataMinus.gESG)/dT))/dT;
        dgER=(excDataPlus.gER-excData.gER)/dT;
        d2gER=(dgER-((excData.gER-excDataMinus.gER)/dT))/dT;
        dgE=dgEC+dgESG+dgER;
        d2gE=d2gEC+d2gESG+d2gER;
        if (ver==1) printf("gE: %f %f %f %f\n",excData.gEC,excData.gESG,excData.gER,excData.gE);*/

        //Next we get the parameters of the individual substances
        for (i=0;i< mix->numSubs;i++)//First we get the parameters of the individual substances
        {
            FF_FixedParamCubic(&mix->cubicData[i],&sParam[i]);
            FF_ThetaDerivCubic(T,&mix->cubicData[i],&sParam[i]);
            FF_ThetaDerivCubic(&Tplus,&mix->cubicData[i],&sParamPlus[i]);
            FF_ThetaDerivCubic(&Tminus,&mix->cubicData[i],&sParamMinus[i]);
            if (ver==1) printf("i a,b,Theta,dTheta,d2Theta: %i %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta*1e3,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
        }
        param->Theta=0;
        param->b=0;
        param->c=0;
        param->u=sParam[0].u;
        param->w=sParam[0].w;
        for (i=0;i< mix->numSubs;i++){
            if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))){
                printf("Cubic EOS are of different types\n");
                return;//We check that all cubic eos are of the same type
            }
            if (sParam[i].c>0) param->c=param->c+x[i]*sParam[i].c;//Calculation of mix volume translation
        }
        double sPartA=0,sPartB=0,sPartBplus=0,sPartBminus=0;//sum of individual substances contribution
        double dsPartB=0,d2sPartB=0;//Its derivatives regarding T
        double q1=0,q2=0,lambda=0;//Parameters for  alpha (and Theta) calculation, lambda is for LCVM mixing rule

        switch (mix->mixRule)
        {
        case FF_HV:
            if (mix->eosType==FF_CubicPRtype) q1=-0.623;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.693;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gESG+excData.gER)/q1+sPartB);
                param->dTheta=param->b*(R*(excData.gESG+excData.gER)/q1+RT*(dgESG+dgER)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*(dgESG+dgER)/q1+RT*(d2gESG+d2gER)/q1+d2sPartB);
            }
            break;
        case FF_MHV1:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    sPartBplus=sPartBplus+x[i]*(sParamPlus[i].Theta/sParam[i].b);
                    sPartBminus=sPartBminus+x[i]*(sParamMinus[i].Theta/sParam[i].b);
                    //dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    //d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                ThetaPlus=param->b*(R*Tplus*(excDataPlus.gE+sPartA)/q1+sPartBplus);
                ThetaMinus=param->b*(R*Tminus*(excDataMinus.gE+sPartA)/q1+sPartBminus);
                //printf("gE:%f Theta:%f\n",excData.gE,param->Theta);
                param->dTheta=(ThetaPlus-param->Theta)/dT;
                param->d2Theta=(ThetaPlus-2*param->Theta+ThetaMinus)/(dT*dT);
                //param->dTheta=param->b*(R*(excData.gE+sPartA)/q1+RT*(dgE)/q1+dsPartB);
                //printf("dTheta:%f\n",param->dTheta);
                //param->d2Theta=param->b*(2*R*dgE/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        case FF_PSRK:
            if (mix->eosType==FF_CubicSRKtype){
                q1=-0.64663;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    sPartBplus=sPartBplus+x[i]*(sParamPlus[i].Theta/sParam[i].b);
                    sPartBminus=sPartBminus+x[i]*(sParamMinus[i].Theta/sParam[i].b);
                    //dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    //d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gE+sPartA)/q1+sPartB);
                ThetaPlus=param->b*(R*Tplus*(excDataPlus.gE+sPartA)/q1+sPartBplus);
                ThetaMinus=param->b*(R*Tminus*(excDataMinus.gE+sPartA)/q1+sPartBminus);
                param->dTheta=param->b*(R*(excData.gE+sPartA)/q1+RT*(dgE)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*dgE/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        case FF_LCVM:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.623;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.693;//Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if ((mix->actModel==FF_UNIFACStd)||(mix->actModel==FF_UNIFACPSRK)) lambda=0.36;
            else if ((mix->actModel==FF_UNIFACDort)||(mix->actModel==FF_UNIFACNist)) lambda=0.65;
            if (!(q1==0)){
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartA=sPartA+x[i]*log(param->b/sParam[i].b);
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                //param->Theta=param->b*(lambda*(RT*excData.gER/q2+sPartB)+(1-lambda)*(RT*(excData.gE+sPartA)/q1+sPartB));
                param->Theta=param->b*(lambda*(RT*excData.gER/q2)+(1-lambda)*(RT*(excData.gE+sPartA)/q1)+sPartB);
                param->dTheta=param->b*(lambda*(R*excData.gER+RT*dgER)/q2+(1-lambda)*(R*(excData.gE+sPartA)+RT*dgE)/q1+dsPartB);
                param->d2Theta=param->b*(lambda*(2*R*dgER+RT*d2gER)/q2+(1-lambda)*(2*R*dgE+RT*d2gE)/q1+d2sPartB);
            }
            break;
        case FF_MHV2:
            if (mix->eosType==FF_CubicPRtype)
            {
                q1=-0.4347;
                q2=-0.003654;
            }
            if (mix->eosType==FF_CubicSRKtype)
            {
                q1=-0.4783;
                q2=-0.0047;
            }
            if (!(q2==0)){
                double alphai,dalphai,d2alphai;
                for (i=0;i<mix->numSubs;i++){
                    param->b=param->b+x[i]*sParam[i].b;
                }
                for (i=0;i<mix->numSubs;i++){
                    alphai=sParam[i].Theta/sParam[i].b/RT;
                    dalphai=sParam[i].dTheta/sParam[i].b/RT-alphai/ *T;
                    d2alphai=sParam[i].d2Theta/sParam[i].b/RT-dalphai/ *T-(dalphai-alphai)/(*T * *T);
                    sPartB=sPartB+x[i]*(log(param->b/sParam[i].b)+q1*alphai+q2*pow(alphai,2));
                    dsPartB=dsPartB+x[i]*(q1*dalphai+2*q2*alphai*dalphai);
                    d2sPartB=d2sPartB+x[i]*(q1*d2alphai+2*q2*(dalphai*dalphai+alphai*d2alphai));

                }
                double aux=q1*q1+4*q2*(sPartB+excData.gE);
                double daux=4*q2*(dsPartB+d2gE);
                param->Theta=param->b*RT*(-q1-pow(aux,0.5))/(2*q2);
                //printf("gE:%f Theta:%f\n",excData.gE,param->Theta);
                param->dTheta=param->Theta/ *T-param->b*RT*pow(aux,-0.5)*(dsPartB+dgE);
                //printf("dTheta:%f\n",param->dTheta);
                param->d2Theta=(param->dTheta* *T-param->Theta)/(*T * *T)-param->b*R*(pow(aux,-0.5)*(dsPartB+dgE)+
                               *T*(-0.5*pow(aux,-1.5)*daux*(dsPartB+dgE)+pow(aux,-0.5)*(d2sPartB+d2gE)));
            }
            break;
        case FF_UMR:
            if (mix->eosType==FF_CubicPRtype) q1=-0.53;
            else if (mix->eosType==FF_CubicSRKtype) q1=-0.593;
            if (!(q1==0)){
                for(i=0;i<mix->numSubs;i++){
                    for(j=0;j<mix->numSubs;j++){
                        bij=pow((pow(sParam[i].b,0.5)+pow(sParam[j].b,0.5))/2,2);
                        param->b=param->b+x[i]*x[j]*bij;
                    }
                }
                for (i=0;i<mix->numSubs;i++){
                    sPartB=sPartB+x[i]*(sParam[i].Theta/sParam[i].b);
                    dsPartB=dsPartB+x[i]*(sParam[i].dTheta/sParam[i].b);
                    d2sPartB=d2sPartB+x[i]*(sParam[i].d2Theta/sParam[i].b);
                }
                param->Theta=param->b*(RT*(excData.gER)/q1+sPartB);
                param->dTheta=param->b*(R*(excData.gER)/q1+RT*(dgER)/q1+dsPartB);
                param->d2Theta=param->b*(2*R*dgER/q1+RT*d2gE/q1+d2sPartB);
            }
            break;
        default:
            param->c=0;
            param->b=0;
            param->Theta=0;
            param->dTheta=0;
            param->d2Theta=0;
            break;
        }
    if (ver==1) printf("Theta:%f dTheta:%f d2Theta:%f b:%f c:%f\n",param->Theta,param->dTheta, param->d2Theta, param->b,param->c);
}

//Mixture SAFT EOS calculations
//==============================

//Z and Arr calculation for a mixture, given T and V, according to FF_PCSAFT EOS. combRul=2
//------------------------------------------------------------------------------
void CALLCONV FF_MixArrZfromTVSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                          FF_SaftEOSdata data[],const double pintParam[15][15][6],const double x[],double *Arr,double *Z)
{
    //Initial calculations: molecular volume, molecular density, and hard molecule volume fraction, etc.
    int i,j,k;
    double Vmolecular,rhoM,mM=0.0,substRho[*numSubs],d[*numSubs],pairSigma[*numSubs][*numSubs],pairEpsilon[*numSubs][*numSubs],
            pairD[*numSubs][*numSubs],pairKAB[*numSubs][*numSubs],pairEpsilonAB[*numSubs][*numSubs];
    double Add=0,Zdd=0;
    double acidAux;
    Vmolecular = *V / Av;
    rhoM = 1 / Vmolecular;
    for (i=0;i<*numSubs;i++)
    {
        mM=mM+x[i]*data[i].m;
        substRho[i]=x[i]*rhoM;
        d[i] = data[i].sigma * (1 - 0.12 * exp(-3 * data[i].epsilon / *T))*1e-10;//Diameter of each segment of the chain in m
    }
    //************
    //int combRul=2;
    //************

    for (i=0;i<*numSubs;i++)
        for (j=0;j<*numSubs;j++)
        {
            pairSigma[i][j]=(data[i].sigma+data[j].sigma)*1e-10/2;
            //pairEpsilon[i][j]=pow((data[i].epsilon*data[j].epsilon),0.5);
            pairEpsilon[i][j]=(1-pintParam[i][j][0]-pintParam[i][j][1]* *T -pintParam[i][j][2]/ *T)*pow((data[i].epsilon*data[j].epsilon),0.5);
            pairD[i][j]=(d[i]+d[j])/2;

            if ((*rule==FF_IndAssoc)&&(data[i].nPos==1)&&(data[i].kAB==0)) data[i].kAB=data[j].kAB;
            else if ((*rule==FF_IndAssoc)&&(data[j].nPos==1)&&(data[j].kAB==0)) data[j].kAB=data[i].kAB;
            //printf("kAB:%f \n",data[i].kAB);

            //if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
            //else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
            pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5)*pow(2*pow(data[i].sigma*data[j].sigma,0.5)/(data[i].sigma+data[j].sigma),3);
            pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;


            /*switch (combRul)
            {
            case 1://CR-1 combining rule: (epsilonAB+epsilonAB)/2, kAB=(kAB*kAB)^0.5
                if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
                else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
                else pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
                break;
            case 2://CR-1 modified by Wolbach and Sandler
                //if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
                //else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
                pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5)*pow(2*pow(data[i].sigma*data[j].sigma,0.5)/(data[i].sigma+data[j].sigma),3);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
                break;

            case 6://CR-1 modified by Wolbach and Sandler, to use with cross-association.
                if ((data[i].kAB==0)&&(1==1)){
                    data[i].kAB=data[j].kAB;
                    data[i].epsilonAB=data[j].epsilonAB;
                }
                pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5)*pow(2*pow(data[i].sigma*data[j].sigma,0.5)/(data[i].sigma+data[j].sigma),3);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
                break;

            case 3://Gross and Sadowski
                if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
                else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
                else pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5)*2*pow(data[i].sigma*data[j].sigma,0.5)/(data[i].sigma+data[j].sigma);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
            case 4://combRul (epsilon*epsilon)**0.5, kAB=(kAB+kAB)/2. Not used
                pairKAB[i][j]=(data[i].kAB+data[j].kAB)/2;
                pairEpsilonAB[i][j]=pow((data[i].epsilonAB*data[j].epsilonAB),0.5);
                break;
            case 5://combRul (epsilon*epsilon)**0.5 , kAB=(...)**3. Not used
                pairKAB[i][j]=pow((pow(data[i].kAB,1/3)+pow(data[j].kAB,1/3))/2,3);
                pairEpsilonAB[i][j]=pow((data[i].epsilonAB*data[j].epsilonAB),0.5);
                break;

            }*/
        }
    //Contribution by hard sphere
    double dseta[4],auxDseta,Zhs,Ahs;
        for (i=0;i<4;i++)
        {
            auxDseta=0.0;
            for (j=0;j<*numSubs;j++)
                auxDseta=auxDseta+x[j]*data[j].m*pow(d[j],i)/6;
            dseta[i]=Pi*rhoM*auxDseta;
        }
        Zhs=mM*(1/(1-dseta[3])+3*dseta[1]*dseta[2]/dseta[0]/pow((1-dseta[3]),2)+(3*pow(dseta[2],3)-dseta[3]*pow(dseta[2],3))/dseta[0]/pow((1-dseta[3]),3)-1);
        Ahs=mM*1/dseta[0]*(3*dseta[1]*dseta[2]/(1-dseta[3])+pow(dseta[2],3)/dseta[3]/pow((1-dseta[3]),2)+(pow(dseta[2],3)/pow(dseta[3],2)-dseta[0])*log(1-dseta[3]));
    //Contribution by chain
    double pairGhs[*numSubs][*numSubs],pairDelta[*numSubs][*numSubs],dLghs_drho[*numSubs],Achain=0.0,Zchain=0.0,etaM,dLghsM_drho;
    for (i=0;i<*numSubs;i++)
    {
        for (j=0;j<*numSubs;j++)
        {
            pairGhs[i][j]=1/(1-dseta[3])+(d[i]*d[j]/(d[i]+d[j]))*3*dseta[2]/pow((1-dseta[3]),2)+
            pow((d[i]*d[j]/(d[i]+d[j])),2)*2*pow(dseta[2],2)/pow((1-dseta[3]),3);
            pairDelta[i][j]=pow(pairSigma[i][j],3)*pairGhs[i][j]*pairKAB[i][j]*(exp(pairEpsilonAB[i][j]/ *T)-1);//This will be used in contribution by molecular association
        }
        //Now we calculate d(ln(ghs(di))/d(rhoM) for each substance
        dLghs_drho[i] =(dseta[3]/pow((1-dseta[3]),2)+3*d[i]*dseta[2]/2/pow((1-dseta[3]),2)+3*d[i]*dseta[2]*dseta[3]/pow((1-dseta[3]),3)+
            pow(d[i],2)*pow(dseta[2],2)/pow((1-dseta[3]),3)+3*pow(d[i],2)*pow(dseta[2],2)*dseta[3]/2/pow((1-dseta[3]),4))*Vmolecular/pairGhs[i][i];
        Achain = Achain + x[i]*(1-data[i].m)* log(pairGhs[i][i]);
        Zchain= Zchain + x[i]*(1-data[i].m)* dLghs_drho[i];
    }
    Zchain= rhoM*Zchain;
    etaM=dseta[3];
    dLghsM_drho=(5*etaM/2-(pow(etaM,2)))/(1-etaM)/(1-0.5*etaM)*Vmolecular; //Derivative of ghs(mean diameter). Used later in association

    //contribution by dispersion
    double I[4]={0.0,0.0,0.0,0.0};
    FF_calcI1I2(mM,etaM,I);
    double sum1=0.0,sum2=0.0,Z1,C1,C2,Z2,Zdisp,Adisp;
    for (i=0;i<*numSubs;i++)
        for (j=0;j<*numSubs;j++)
        {	sum1=sum1+x[i]*x[j]*data[i].m*data[j].m*pairEpsilon[i][j]*pow(pairSigma[i][j],3);
            sum2=sum2+x[i]*x[j]*data[i].m*data[j].m*pow(pairEpsilon[i][j],2)*pow(pairSigma[i][j],3);
        }
    sum1=sum1/ *T;
    sum2=sum2/pow(*T,2);
    Z1 = -2 * Pi / Vmolecular * I[1] * sum1;
    C1 = 1/(1 + mM * (8 * etaM - 2 * pow(etaM,2)) / pow((1 - etaM),4) + (1 - mM) *
            (20 * etaM - 27 * pow(etaM,2)+ 12 * pow(etaM,3) - 2 * pow(etaM,4)) / pow(((1 - etaM) * (2 - etaM)),2));
    C2 = C1 * (mM * (-4 * pow(etaM,2) + 20 * etaM + 8) / pow((1 - etaM),5) + (1 - mM) * (2 * pow(etaM,3)+
            12 * pow(etaM,2) - 48 * etaM + 40) / pow(((1 - etaM) * (2 - etaM)),3));
    Z2 = -Pi / Vmolecular * mM * C1 * (I[3] - C2 * etaM * I[2])* sum2;
    Zdisp = Z1 + Z2;
    Adisp = - Pi * rhoM *(2* I[0] * sum1 +  mM * C1 * I[2] * sum2);

    //contribution by molecular association
    double xPos[*numSubs],xNeg[*numSubs],xAcid[*numSubs],Aassoc=0.0,Zassoc=0.0;//Fraction of non associated sites in each molecule
    for (i=0;i<*numSubs;i++) //initial values for the non-associated fraction
    {
        if (data[i].nPos>0) xPos[i]=xNeg[i]=0.5; else xPos[i]=xNeg[i]=1.0;
        if (data[i].nAcid==0) xAcid[i]=1.0;
        else if (data[i].nAcid==1) xAcid[i]=(-1 + pow((1 + 4 * substRho[i] * pairDelta[i][i]),0.5)) / (2 * substRho[i] * pairDelta[i][i]);
        else if (data[i].nAcid==2) xAcid[i]=(-1 + pow((1 + 8 * substRho[i] * pairDelta[i][i]),0.5)) / (4 * substRho[i] * pairDelta[i][i]);
        else if (data[i].nAcid==3) xAcid[i]=(-1 + pow((1 + 12 * substRho[i] * pairDelta[i][i]),0.5)) / (6 * substRho[i] * pairDelta[i][i]);
        else if (data[i].nAcid==4) xAcid[i]=(-1 + pow((1 + 16 * substRho[i] * pairDelta[i][i]),0.5)) / (8 * substRho[i] * pairDelta[i][i]);
        //printf("xAcid[%i]:%f \n",i,xAcid[i]);
    }
    for (k=0;k<13;k++){ //This is a iteration to approximate non associated fraction
        for (i=0;i<*numSubs;i++){
            if (data[i].nPos>0)
            {
                xPos[i]=0.0;
                for (j=0;j<*numSubs;j++)
                    xPos[i]=xPos[i]+substRho[j]*(data[j].nNeg*xNeg[j])*pairDelta[i][j];
                xPos[i]=1/(1+xPos[i]);
            }
            if (data[i].nNeg>0)
            {
                xNeg[i]=0.0;
                for (j=0;j<*numSubs;j++)
                    xNeg[i]=xNeg[i]+substRho[j]*(data[j].nPos*xPos[j])*pairDelta[i][j];
                xNeg[i]=1/(1+xNeg[i]);
            }

            if (data[i].nAcid>0)
            {
                acidAux=0;
                for (j=0;j<*numSubs;j++)
                    acidAux=acidAux+substRho[j]*data[j].nAcid*xAcid[j]*pairDelta[i][j];
                    xAcid[i]=1/(1+acidAux);
                                }
            //printf("K:%i xPos[%i]:%f xNeg[%i]:%f xAcid:%f \n",k,i,xPos[i],i,xNeg[i],xAcid[i]);
        }
    }
    //for (i=0;i<*numSubs;i++) printf("xPos[%i]:%f xNeg[%i]:%f xAcid[%i]:%f substRho[%i]:%e pairDelta[%i][%i]:%e \n",i,xPos[i],i,xNeg[i],i,xAcid[i],i,substRho[i],i,i,pairDelta[1][1]);
    for (i=0;i<*numSubs;i++){

        Aassoc=Aassoc+x[i]*(((data[i].nPos+data[i].nNeg+data[i].nAcid)/2)+data[i].nPos*(log(xPos[i])-xPos[i]/2)+data[i].nNeg*(log(xNeg[i])-xNeg[i]/2)+
                            data[i].nAcid*(log(xAcid[i])-xAcid[i]/2));
        Zassoc=Zassoc+x[i]*(data[i].nPos*(1-xPos[i])+data[i].nNeg*(1-xNeg[i])+data[i].nAcid*(1-xAcid[i]));
    }
        Zassoc=-0.5*(1+rhoM*dLghsM_drho)*Zassoc;


        //contribution by polar forces
        //1Debbie=3.33564 e-30 C.m(SI units)
        //U=mu*mu/(4*pi*epsilon0*r^3)
        int polar=0;
        for(i=0;i<*numSubs;i++){
            if(data[i].xp>0){
                polar=1;
                break;
            }
        }
        if(polar==1){
            int n;
            double AddSum=0,dAddSum_dEta=0;
            double Add2,dAdd2_dRhoM,Add3,dAdd3_dRhoM;
            double mEfec;
            double ad[3][5]={{0.3043504,-0.1358588,1.4493329,0.3556977,-2.0653308},{0.9534641,-1.8396383,2.0131180,-7.3724958,8.2374135},{-1.1610080,4.5258607,0.9751222,-12.281038,5.9397575}};
            double bd[3][5]={{0.2187939,-1.1896431,1.1626889,0.0,0.0},{-0.5873164,1.2489132,-0.5085280,0.0,0.0},{3.4869576,-14.915974,15.372022,0.0,0.0}};
            double cd[3][5]={{-0.0646774,0.1975882,-0.8087562,0.6902849,0.0},{-0.9520876,2.9924258,-2.3802636,-0.2701261,0.0},{-0.6260979,1.2924686,1.6542783,-3.4396744,0.0}};
            double sigmaS3[*numSubs],epsilonST[*numSubs],muRedS2[*numSubs];
            double a,b,c;
            double Jdd,dJdd_eta;
            double d3,mEfecAux,aux1;

            //Previous preparation
            d3=6*etaM/(Pi*rhoM*mM);
            for(i=0;i<*numSubs;i++){//New
                sigmaS3[i]=data[i].sigma*data[i].sigma*data[i].sigma*1e-30;
                epsilonST[i]=data[i].epsilon/ *T;
                muRedS2[i]=pow(data[i].mu,2)/(data[i].m*sigmaS3[i]*data[i].epsilon)*7.24311E-27;
                //printf("i:%i sigmaS3:%f epsilonST:%f muRedS2:%f\n",i,sigmaS3[i]*1e30,epsilonST[i],muRedS2[i]);
            }
            //printf("etaM:%f\n",etaM);
            //if (data->m>2) mEfec=2;
            //else mEfec=data->m;
            mEfec=2;
            mEfecAux=(mEfec-1)*(mEfec-2)/(mEfec*mEfec);

            //printf("muRed:%f\n",pow(muRed2,0.5));
            AddSum=0;
            dAddSum_dEta=0;
            for(i=0;i<*numSubs;i++){//Simplifications: taken m=2
                for(j=0;j<*numSubs;j++){
                    Jdd=0;
                    dJdd_eta=0;
                    for (n=0;n<5;n++){
                        a=ad[0][n]+ad[1][n]*(mEfec-1)/mEfec+ad[2][n]*mEfecAux;
                        b=bd[0][n]+bd[1][n]*(mEfec-1)/mEfec+bd[2][n]*mEfecAux;
                        aux1=(a+b*pairEpsilon[i][j]/ *T)*pow(etaM,n);
                        Jdd=Jdd+aux1;
                        dJdd_eta=dJdd_eta+n*aux1/etaM;
                    }
                    //printf("Jdd2:%f dJdd2_Eta:%f\n",Jdd,dJdd_eta);
                    aux1=x[i]*x[j]*epsilonST[i]*epsilonST[j]*sigmaS3[i]*sigmaS3[j]*data[i].xp*data[j].xp*muRedS2[i]*muRedS2[j]/
                            (pairSigma[i][j]*pairSigma[i][j]*pairSigma[i][j]);
                    AddSum=AddSum+aux1*Jdd;
                    dAddSum_dEta=dAddSum_dEta+aux1*dJdd_eta;

                }
            }

            Add2=-Pi*rhoM*AddSum;
            dAdd2_dRhoM=-Pi*(AddSum+rhoM*dAddSum_dEta*Pi*mM*d3/6);//Here we change the derivative from eta to rhoM

            AddSum=0;
            dAddSum_dEta=0;
            for(i=0;i<*numSubs;i++){//Simplifications: taken m=2
                for(j=0;j<*numSubs;j++){
                    for(k=0;k<*numSubs;k++){
                        Jdd=0;
                        dJdd_eta=0;
                        for (n=0;n<5;n++){
                            c=cd[0][n]+cd[1][n]*(mEfec-1)/mEfec+cd[2][n]*mEfecAux;
                            aux1=c*pow(etaM,n);
                            Jdd=Jdd+aux1;
                            dJdd_eta=dJdd_eta+n*aux1/etaM;
                        }
                        //printf("Jdd3:%f dJdd3_Eta:%f\n",Jdd,dJdd_eta);
                        aux1=x[i]*x[j]*x[k]*epsilonST[i]*epsilonST[j]*epsilonST[k]*sigmaS3[i]*sigmaS3[j]*sigmaS3[k]*data[i].xp*data[j].xp*data[k].xp*
                                muRedS2[i]*muRedS2[j]*muRedS2[k]/(pairSigma[i][j]*pairSigma[i][k]*pairSigma[j][k]);
                        AddSum=AddSum+aux1*Jdd;
                        dAddSum_dEta=dAddSum_dEta+aux1*dJdd_eta;
                    }
                }
            }
            Add3=-4*pow(Pi*rhoM,2)*AddSum/3;
            dAdd3_dRhoM=-4*(2*Pi*Pi*rhoM*AddSum+pow(Pi*rhoM,2)*dAddSum_dEta*Pi*mM*d3/6)/3;//Here we change the derivative from eta to rhoM
            aux1=(1-Add3/Add2);
            Add=Add2/aux1;
            Zdd=rhoM*(dAdd2_dRhoM*aux1+dAdd3_dRhoM-Add3*dAdd2_dRhoM/Add2)/(aux1*aux1);
            //printf("Add2:%f Add3:%f Zdd2:%f Add:%f Zdd:%f\n",Add2,Add3,rhoM*dAdd2_dRhoM,Add,Zdd);
            //printf("a:%f b*1e3:%f c*1e3:%f\n",a,b*1e3,c*1e3);
        }

    *Arr=Ahs+Achain+Adisp+Aassoc+Add;//Arr
    *Z=1+Zhs+Zchain+Zdisp+Zassoc+Zdd;//Z
    //printf("Arr:%f Z:%f\n",*Arr,*Z);

}

//Mixture P calculation given T and V, according to FF_PCSAFT EOS
//---------------------------------------------------------------
void CALLCONV FF_MixPfromTVSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[15][15][6],const double x[],double *P)
{
    double Arr,Z;
    FF_MixArrZfromTVSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
    //printf("T:%f V:%f Z:%f\n",*T,*V,Z);
    *P=Z*R* *T/ *V;
}


//Mixture V,Arr and Z  calculation, given T and P and composition, according to FF_PCSAFT EOS
//-------------------------------------------------------------------------------------------
void CALLCONV FF_MixVfromTPSAFTOld(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                              const double pintParam[15][15][6],const double x[],const char *option,double resultL[3],double resultG[3],char *state)
{
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double V,Arr,Z,Vplus,ArrPlus,Zplus,error,errorPrev,errorNew,errorRel,dP_dV,dP_dVnew;
    double maxError=0.000001,Vprev,Vnew,dP_dVprev;
    int i;
    double dM=0,mM=0;
    double eta,Vmolecular;
    for (i=0;i<*numSubs;i++)
    {
        dM = dM + x[i]*data[i].sigma * (1 - 0.12 * exp(-3 * data[i].epsilon / *T));
        mM = mM + x[i]*data[i].m;
    }
    eta = 0.6;  //As initial guess, we suppose the volume fraction occupied by the hard molecule is 0.4
    Vmolecular = mM * (Pi * pow(dM,3) / 6) / eta * Av / 1E+30; //This will be the initial guess for the molecular volume in m3/mol
    if (*option!='g'){//we calculate the liquid phase if not only the gas phase has been asked for
        V = Vmolecular; //initial guess for mole volume in m3
            //Vplus=V; //Vplus will be a volume that produces a P lower than that wanted
        Vplus=V*1.0000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        while (dP_dV>=0){//Initial loop to find a negative derivative of P regarding V
            //printf("finding a liquid negative derivative. V:%f\n",V);
            V = V*1.5;
            Vplus=V*1.0000001;
            FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        }
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        Vprev=V;
        errorPrev=error;
        dP_dVprev=dP_dV;
        i=1;
        //printf("Initial liquid:T:%f V:%f dP/dV:%f err:%f\n",*T,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(i<41))
        {
            if(error*errorPrev<0){
                Vnew=(V+Vprev)/2;
                FF_MixArrZfromTVSAFT(rule,T,&Vnew,numSubs,data,pintParam,x,&Arr,&Z);
                errorNew =*P-R* *T*Z/Vnew;
                if(errorNew*errorPrev>0){
                    Vprev=Vnew;
                    errorPrev=errorNew;
                }
                else{
                    V=Vnew;
                    error=errorNew;
                }
                errorRel=error/ *P;
                i++;

            }
            else if(dP_dV>0){//It is necessary to determinate the minimum
                Vnew=(V+Vprev)/2;
                Vplus=Vnew*1.0000001;//we calculate dP/dV numerically
                FF_MixArrZfromTVSAFT(rule,T,&Vnew,numSubs,data,pintParam,x,&Arr,&Z);
                FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
                dP_dVnew=R* *T*(Zplus/Vplus-Z/ Vnew)/(Vplus-Vnew);
                errorNew =*P-R* *T*Z/V;
                if(dP_dVnew<0){
                    Vprev=Vnew;
                    errorPrev=errorNew;
                    dP_dVprev=dP_dVnew;
                }
                else{
                    V=Vnew;
                    error=errorNew;
                    dP_dV=dP_dVnew;
                }
                errorRel=error/ *P;
                i++;
            }
            else{
                V=V+error/dP_dV;//Newton method for root finding
                Vplus=V*1.0000001;//we calculate dP/dV numerically
                FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
                FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
                dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
                error =*P-R* *T*Z/V;
                errorRel=error/ *P;
                i++;
            }
            //printf("i:%d Vl:%f Vplus:%f Z:%f Zplus:%fdP_dV:%f error:%f\n",i,V,Vplus,Z,Zplus,dP_dV,error);
        }
        if ((fabs(errorRel)<maxError)){
            resultL[0]=V;
            resultL[1]=Arr;
            resultL[2]=Z;
            *state='l';//We inform that we have done the calculation from the liquid end
        }
        else resultL[0]=resultL[1]=resultL[2]=0;
    }
    if (*option!='l')//and the gas phase if not only the liquid one has been asked for
    {
        V = R * *T / *P + Vmolecular;//initial guess for gas volume
        Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        i=1;
        //printf("%d V:%f dP/dV:%f error:%f\n",i,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(dP_dV <0)&&(V>0)&&(i<51))
        {
            V=V+error/dP_dV;//Newton method for root finding
            if ((V<=0)||(dP_dV>0))//We slow the Newton method if V is negative, or dP/dV is positive
            {
               V=V-0.9*error/dP_dV;//Newton method for root finding, slowed
            }

            Vplus=V*1.000001;
            FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            error =*P-R* *T*Z/V;
            errorRel=error/ *P;
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus- V);
            i=i+1;
            //printf("i:%d V:%f dP_dV:%f error:%f\n",i,V,dP_dV,error);
        }
        if (fabs(errorRel)<maxError){
            resultG[0]=V;
            resultG[1]=Arr;
            resultG[2]=Z;
            if (*state=='l') *state='b';
            else *state='g';
            //printf("hola\n");
        }
        else resultG[0]=resultG[1]=resultG[2]=0;
    }
}


//Z,Arr,V calculation for a pure substance, given T and P, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------
void CALLCONV FF_MixVfromTPSAFT(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                                   const double pintParam[15][15][6],const double x[],const char *option,double resultL[3],double resultG[3],char *state)
{
    //Interesting to determine max after finding a positive dP/dV from gas. Probably not so interesting minimum from liquid. All by regula falsi
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double V,Arr,Z,Pcalc,Vplus,ArrPlus,Zplus,error,errorRel,dP_dV,Zprev,Vmax,Vmin,Emax,Emin,dP_dVprev;
    int i,ok=0,side=0;
    double maxError=0.0001;//maximum relative error accepted
    double eta,Vinit;

    double sigmaM=0,mM=0;
    for (i=0;i<*numSubs;i++)
    {
        sigmaM = sigmaM + x[i]*data[i].sigma;
        mM = mM + x[i]*data[i].m;
    }

    eta = 0.55;  //As initial guess, we suppose the volume fraction occupied by the hard molecule is 0.55
    Vinit = mM * (Pi * pow(sigmaM,3) / 6) / eta * Av / 1E+30; //This will be the initial guess for the molecular volume in m3/mol
    //printf("Initial Vl guess:%f\n", Vinit);
    if (*option!='g')//we calculate the liquid phase if not only the gas phase has been asked
    {
        V = Vinit; //initial guess for liquid mole volume in m3
        Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        //FF_ArrZfromTVSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        while (dP_dV>=0){//Initial loop to find a negative derivative of P regarding V. Normally not necessary
            //printf("finding a liquid negative derivative V:%f\n",V);
            V = V*1.2;
            Vplus=V*1.000001;
            FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
            //FF_ArrZfromTVSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        }
        Pcalc=Z*R* *T/V;
        error=*P-Pcalc;
        Vmin=V;//We store a volume that gives a higher pressure
        Emin=error;
        Zprev=Z;
        for(i=0;i<20;i++){//Newton method to find a volume with pressure lower than objective,or a possitive dP/dV
            V=V+error/dP_dV;
            Zprev=Z;
            FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
            Pcalc=Z*R* *T/V;
            error=*P-Pcalc;
            //Calculation of the new dP/dV derivative
            if(i>6) Vplus=V*1.00001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
            else Vplus=V*1.000001;//Close to the critical point the increment must be larger
            FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            //FF_ArrZfromTVSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
            //printf("Looking for liq.root V:%f P:%f dP_dV:%f Z:%f\n",V,Pcalc,dP_dV,Z);
            if((dP_dV>0)||((Z>Zprev)&&(Z>0.4))) break;//there is no liquid solution
            if ((fabs(error/ *P)<maxError)||(error>0)){//We stop when we get a lower pressure
                Vmax=V;//we store a volume that gives a lower pressure
                Emax=error;
                ok=1;
                break;
            }
            else{

                Vmin=V;//We store a volume that gives a higher pressure
                Emin=error;
            }
        }
        if(ok==1){//We refine only if a lower pressure has been found
            i=0;
            while ((i<25)&&(fabs(error/ *P)>maxError)){//we begin to refine the root
                i++;
                V=(Vmin*Emax-Vmax*Emin)/(Emax-Emin);//This is the regula falsi method, Anderson-Bjork modified
                FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
                //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
                Pcalc=Z*R* *T/V;
                error=*P-Pcalc;
                //printf("Finding liquid root V:%f P:%f\n",V,Pcalc);
                if (fabs(error/ *P)<0.0001) break;//if we have arrived to the solution exit
                if ((Emax * error)>0){//if the proposed solution is of the same sign than error(Vmax)
                    if (side==-1){//if it happened also in the previous loop
                        if ((1-error/Emax)>0) Emin *= (1-error/Emax);//we decrease f(xmin) for the next loop calculation
                        else Emin *= 0.5;
                    }
                    Vmax=V;
                    Emax=error;
                    side = -1;//we register than the solution was of the same sign than previous xmax
                }
                else{//If the prosed solution is of the same sign than f(xmin) we apply the same technic
                    if (side==1){
                        if ((1-error/Emin)>0) Emax *= (1-error/Emin);
                        else Emax *= 0.5;
                    }
                    Vmin=V;
                    Emin=error;
                    side = 1;
                }
            }
        }
        if ((ok==1)&&(fabs(errorRel)<maxError)){
            resultL[0]=V;
            resultL[1]=Arr;
            resultL[2]=Z;
            *state='l';//We inform that we have done the calculation from the liquid end
        }
        else resultL[0]=resultL[1]=resultL[2]=0;
    }

    if (*option!='l'){//and the gas phase if not only the liquid one has been asked for
        ok=0;
        V = R * *T / *P + Vinit;//initial guess for gas volume
        Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
        //FF_ArrZfromTVSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        Pcalc=Z*R* *T/V;
        error=*P-Pcalc;
        //printf("Initial gas root V:%f P:%f dP_dV:%f Z:%f\n",V,Pcalc,dP_dV,Z);
        Vmax=V;//We store a volume that gives a lower pressure
        Emax=error;
        dP_dVprev=dP_dV;
        double Vprev=V;
        int finish=0;
        for(i=0;i<20;i++){//Newton method to find a volume with pressure higher than objective,or a possitive dP/dV
            if(dP_dV>dP_dVprev){//We notice the pass by the transition zone
                V=0.9*Vprev;
                Vprev=V;
                finish=1;
            }
            else{
                Vprev=V;
                dP_dVprev=dP_dV;
                V=V+0.5*error/dP_dV;//As pressure increase much faster than the prediction of the derivative it is necessary to go slowly
                if(V<0.7*Vprev) V=0.7*Vprev;
            }
            FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
            Pcalc=Z*R* *T/V;
            error=*P-Pcalc;
            //Calculation of the new dP/dV derivative
            if(i>6) Vplus=V*1.0001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
            else Vplus=V*1.00001;//Close to the critical point the increment must be larger
            FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            //FF_ArrZfromTVSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
            //printf("Looking for gas root V:%f P:%f dP_dV:%f Z:%f\n",V,Pcalc,dP_dV,Z);
            if((dP_dV>0)&&(finish==1)){
               break;//there is no gas solution
            }
            if ((fabs(error/ *P)<maxError)||(error<0)){//We stop when we get a higher pressure
                Vmin=V;//we store a volume that gives a lower pressure
                Emin=error;
                ok=1;
                break;
            }
            else{
                Vmax=V;//We store a volume that gives a lower pressure
                Emax=error;
            }
        }
        if(ok==1){//We refine only if a higher pressure has been found
            i=0;
            while ((i<25)&&(fabs(error/ *P)>maxError)){//we begin to refine the root
                i++;
                V=(Vmin*Emax-Vmax*Emin)/(Emax-Emin);//This is the regula falsi method, Anderson-Bjork modified
                FF_MixArrZfromTVSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
                //FF_ArrZfromTVSAFT(T,&V,data,&Arr,&Z);
                Pcalc=Z*R* *T/V;
                error=*P-Pcalc;
                //printf("Finding gas root V:%f P:%f\n",V,Pcalc);
                if (fabs(error/ *P)<0.0001) break;//if we have arrived to the solution exit
                if ((Emax * error)>0){//if the proposed solution is of the same sign than error(Vmax)
                    if (side==-1){//if it happened also in the previous loop
                        if ((1-error/Emax)>0) Emin *= (1-error/Emax);//we decrease f(xmin) for the next loop calculation
                        else Emin *= 0.5;
                    }
                    Vmax=V;
                    Emax=error;
                    side = -1;//we register than the solution was of the same sign than previous xmax
                }
                else{//If the prosed solution is of the same sign than f(xmin) we apply the same technic
                    if (side==1){
                        if ((1-error/Emin)>0) Emax *= (1-error/Emin);
                        else Emax *= 0.5;
                    }
                    Vmin=V;
                    Emin=error;
                    side = 1;
                }
            }
        }
        if ((ok==1)&&(fabs(errorRel)<maxError)){
            resultG[0]=V;
            resultG[1]=Arr;
            resultG[2]=Z;
            if (*state=='l'){
                if (fabs((resultL[0]-resultG[0])/resultL[0])>0.001) *state='b';
                else *state='u';
            }
            else *state='g';
            //printf("hola\n");
        }
        else resultG[0]=resultG[1]=resultG[2]=0;
    }
    //printf("state:%c\n",*state);
}




//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a mixture, given T and V, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixArrDerSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                              const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double result[6])
{
    double dV=*V * 0.00001;//increments of V and T used to obtain dArr/dV,dArr/dT and d2Arr/dT2 in SAFT eos
    double Vplus=*V + dV;
    double dT=0.01;
    double Tplus=*T+dT;
    double Tminus=*T-dT;
    double Arr,Z,ArrVplus,ZVplus,ArrTplus,ZTplus,ArrTminus,ZTminus;
    FF_MixArrZfromTVSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
    FF_MixArrZfromTVSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrVplus,&ZVplus);
    FF_MixArrZfromTVSAFT(rule,&Tplus,V,numSubs,data,pintParam,x,&ArrTplus,&ZTplus);
    FF_MixArrZfromTVSAFT(rule,&Tminus,V,numSubs,data,pintParam,x,&ArrTminus,&ZTminus);
    result[0]=Arr;//This is Arr
    result[1]=(1- Z)/ *V;//dArr/dV at constant T
    result[2]=((1- ZVplus)/ Vplus-result[1])/dV;//d2Arr/dV2 at constant T
    result[3]=(ArrTplus-Arr)/dT;//dArr/T at constant V
    result[4]=(result[3]-(Arr-ArrTminus)/dT)/dT;//d2Arr/dT2 at constant V
    result[5]=((1- ZTplus)/ *V-result[1])/dT;//d2Arr/dVdT
}


//Mixtures common calculations
//============================

//Mixture P calculation from T and V by eos
//-----------------------------------------
void CALLCONV FF_MixPfromTVeos(const FF_MixData *mix,const double *T,const double *V,const double x[],double *P)
{
     FF_CubicParam param;
    double Arr,Z;
    double da_di[mix->numSubs],db_di[mix->numSubs],dc_di[mix->numSubs];//just for compatibility of function call
    switch (mix->eosType)
    {
        //case FF_IdealGas:
        //    *P=R* *T/ *V;
        //    break;
        case FF_SAFTtype:
            //*( FF_SaftEOSdata*) data;
            FF_MixArrZfromTVSAFT(&mix->mixRule,T,V,&mix->numSubs,mix->saftData,mix->intParam,x,&Arr,&Z);
            //printf("T:%f V:%f Z:%f\n",*T,*V,Z);
            *P=Z*R* *T/ *V;
            break;
        case FF_SWtype:
            //FF_PfromTVsw(*T,*V,data,P);
            break;
        default://if we have a cubic eos the first step is to calculate its parameters
            if((mix->mixRule==FF_VdW)||(mix->mixRule==FF_VdWnoInt)||(mix->mixRule==FF_PR)||(mix->mixRule==FF_MKP)){
                FF_MixParamXderCubicEOS(&mix->mixRule,T,&mix->numSubs,mix->cubicData,mix->intParam,x,&param,da_di,db_di,dc_di);
            }
            else{
                FF_MixParamCubicEOSgE(mix,T,x,&param);
                //FF_MixParamTderCubicEOSgE(mix,T,x,&param);
            }
            FF_PfromTVcubic(*T,*V,&param,P);
            break;
    }
}

//Mixture V calculation from T and P by eos
//-----------------------------------------
void CALLCONV FF_MixVfromTPeos(const FF_MixData *mix,const double *T,const double *P,const double x[],
                               const char *option,double resultL[3],double resultG[3],char *state)
{
     FF_CubicParam param;//if we have a cubic eos the first step is to calculate its parameters
    double da_di[mix->numSubs],db_di[mix->numSubs],dc_di[mix->numSubs];//just for compatibility with the called function
    switch (mix->eosType)
    {
        case FF_IdealType:
            resultL[0]=resultG[0]=R* *T/ *P;
            resultL[1]=resultG[1]=0;
            resultL[2]=resultG[2]=1;
            *state='u';
            break;
        case FF_SAFTtype:
            //( FF_SaftEOSdata*) data;
            FF_MixVfromTPSAFT(&mix->mixRule,T,P,&mix->numSubs,mix->saftData,mix->intParam,x,option,resultL,resultG,state);
            break;
        case FF_SWtype:
            //*( FF_SWEOSdata*) data;
            //FF_VfromTPsw(eos,*T,*P,data,*option,resultL,resultG,state);
            break;
        default://Cubic eos
            //For a cubic EOs is always necessary to obtain first its parameters
            if((mix->mixRule==FF_VdW)||(mix->mixRule==FF_VdWnoInt)||(mix->mixRule==FF_PR)||(mix->mixRule==FF_MKP)){
                FF_MixParamXderCubicEOS(&mix->mixRule,T,&mix->numSubs,mix->cubicData,mix->intParam,x,&param,da_di,db_di,dc_di);
            }
            else{//gE EOS
                FF_MixParamCubicEOSgE(mix,T,x,&param);
            }
            //printf("Theta:%f, b:%f\n",param.Theta,param.b);
            FF_VfromTPcubic(*T,*P,&param,*option,resultL,resultG,state);
            break;
    }

    //printf("T:%f  P:%f Vl:%f ArrL:%f Zl:%f\n",*T,*P,resultL[0],resultL[1],resultL[2]);
    if (*option=='s'){
        if (*state=='b'){
            if (fabs((resultL[0]-resultG[0])/resultL[0])>0.001){
                if ((resultL[1]+resultL[2]-1-log(resultL[2]))<(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='L';//we compare Gdr
                else if ((resultL[1]+resultL[2]-1-log(resultL[2]))>(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='G';
                else *state='E';//if Gdr is the same we are in equilibrium
            }
            else *state='U';//We got the same result from both extremes. 'l' and 'g' indicates just one result found from just one extreme
        }
    }
}


EXP_IMP void CALLCONV FF_MixPhiEOS(const FF_MixData *mix,const double *T,const double *P,const double x[],const char *option,double phi[])
{//const FF_MixData *mix,const double *T,const double *P,const double x[],const char *option,double phi[]
    FF_CubicParam param;
    char state;
    double resultL[3],resultG[3];
    double V,Arr,Z;
    double xPlus[mix->numSubs],ArrPlus,Zplus;
    double da_di[mix->numSubs],db_di[mix->numSubs],dc_di[mix->numSubs],dArr_dXi[mix->numSubs],sumdArr_dXi=0;
    int i;

    switch (mix->eosType){//Mixture parameters calculation, and initial solution of the EOS for the given T and P
    case FF_SAFTtype:
        FF_MixVfromTPSAFT(&mix->mixRule,T,P,&mix->numSubs,mix->saftData,mix->intParam,x,option,resultL,resultG,&state);
        if(*option=='l'){
            V=resultL[0];
            Arr=resultL[1];
            Z=resultL[2];
        }
        else{
            V=resultG[0];
            Arr=resultG[1];
            Z=resultG[2];
        }
        //printf("option:%c state:%c V:%f Arr:%f Z:%f\n",*option,state,V,Arr,Z);
        for (i=0;i<mix->numSubs;i++) xPlus[i]=x[i];
        for (i=0;i<mix->numSubs;i++)
        {
            xPlus[i]=x[i]*1.00001; //This is the composition increment for the calculation of the variation of red.res.Helmholtz energy
            //printf("i:%i xPlus:%f\n",i,x[i]);
            FF_MixArrZfromTVSAFT(&mix->mixRule,T,&V,&mix->numSubs,mix->saftData,mix->intParam,xPlus,&ArrPlus,&Zplus);
            dArr_dXi[i]=(ArrPlus-Arr)/(xPlus[i]-x[i]);
            sumdArr_dXi=sumdArr_dXi+x[i]*dArr_dXi[i];
            xPlus[i]=x[i];
        }
        for (i=0;i<mix->numSubs;i++)
        {
            phi[i]=exp(dArr_dXi[i]+Arr-sumdArr_dXi+Z-1-log(Z));
            //printf("dArr_dXi:%f\n",dArr_dXi[i]);
            //printf("dArr_dNi:%f\n",dArr_dXi[i]-sumdArr_dXi);
            //printf("phi[%i]: %f\n",i,phi[i]);
        }
        break;
    default:
        if ((mix->mixRule==FF_VdW)||(mix->mixRule==FF_VdWnoInt)||(mix->mixRule==FF_PR)||(mix->mixRule==FF_MKP)){
            FF_MixParamXderCubicEOS(&mix->mixRule,T,&mix->numSubs,mix->cubicData,mix->intParam,x,&param,da_di,db_di,dc_di);
            FF_VfromTPcubic(*T,*P,&param,*option,resultL,resultG,&state);
            if(*option=='l'){
                V=resultL[0];
                Arr=resultL[1];
                Z=resultL[2];
            }
            else{
                V=resultG[0];
                Arr=resultG[1];
                Z=resultG[2];
            }
            for (i=0;i<mix->numSubs;i++)
            {
                //dArr_dXi[i]=param.Theta/(param.b*R* *T*(param.w-param.u))*((da_di[i]/param.Theta-db_di[i]/param.b)*log((V+param.u*param.b)/(V+param.w*param.b))+
                //            V*db_di[i]*(param.u-param.w)/(V+param.u*param.b)/(V+param.w*param.b))+db_di[i]/(V-param.b);//This is without volume traslation
                dArr_dXi[i]=param.Theta/(param.b*R* *T*(param.w-param.u))*((da_di[i]/param.Theta-db_di[i]/param.b)*log((V+param.c+param.u*param.b)/
                           (V+param.c+param.w*param.b))+((V+param.c)*db_di[i]-param.b*dc_di[i])*(param.u-param.w)/
                           ((V+param.c+param.u*param.b)*(V+param.c+param.w*param.b)))+(db_di[i]-dc_di[i])/(V+param.c-param.b);
                sumdArr_dXi=sumdArr_dXi+x[i]*dArr_dXi[i];
            }

            for (i=0;i<mix->numSubs;i++)
            {
                phi[i]=exp(dArr_dXi[i]+Arr-sumdArr_dXi+Z-1-log(Z));
                //printf("dArr_dXi:%f\n",dArr_dXi[i]);
                //printf("dArr_dNi:%f\n",dArr_dXi[i]-sumdArr_dXi);
                //printf("phi[%i]: %f\n",i,phi[i]);
            }
        }
        else{//gE mixing rule
            FF_MixParamNderCubicEOSgE(mix,T,x,&param,da_di,db_di);
            double B=*P*param.b/(R* *T);
            FF_VfromTPcubic(*T,*P,&param,*option,resultL,resultG,&state);
            if(*option=='l'){
                V=resultL[0];
                Arr=resultL[1];
                Z=resultL[2];
            }
            else{
                V=resultG[0];
                Arr=resultG[1];
                Z=resultG[2];
            }

            //printf("V:%f Arr;%f Z;%f\n",V,Arr,Z);
            for (i=0;i<mix->numSubs;i++){
                phi[i]=exp(db_di[i]/param.b*(Z-1)-log(Z-B)- da_di[i] * log((Z+param.u*B)/(Z+param.w*B))/(param.u-param.w));
            }
            //printf("Theta:%f b:%f c:%f\n",param.Theta,param.b,param.c);
        }


        break;
    }

}

//Fast Phi computation for a cubic EOS supplying parameters and derivatives
EXP_IMP void CALLCONV FF_MixPhiEOScubic(const FF_MixData *mix,FF_CubicParam *param,double da_di[],double db_di[],double dc_di[],
                                     const double *T,const double *P,const double x[],const char *option,double phi[]){
    int i;
    char state;
    double V,Arr,Z,resultL[3],resultG[3],dArr_dXi[mix->numSubs],sumdArr_dXi=0;
    double B=*P*param->b/(R* *T);
    if ((mix->mixRule==FF_VdW)||(mix->mixRule==FF_VdWnoInt)||(mix->mixRule==FF_PR)||(mix->mixRule==FF_MKP)){
        FF_VfromTPcubic(*T,*P,param,*option,resultL,resultG,&state);
        if(*option=='l'){
            V=resultL[0];
            Arr=resultL[1];
            Z=resultL[2];
        }
        else{
            V=resultG[0];
            Arr=resultG[1];
            Z=resultG[2];
        }
        for (i=0;i<mix->numSubs;i++)
        {
            //dArr_dXi[i]=param.Theta/(param.b*R* *T*(param.w-param.u))*((da_di[i]/param.Theta-db_di[i]/param.b)*log((V+param.u*param.b)/(V+param.w*param.b))+
            //            V*db_di[i]*(param.u-param.w)/(V+param.u*param.b)/(V+param.w*param.b))+db_di[i]/(V-param.b);//This is without volume traslation
            dArr_dXi[i]=param->Theta/(param->b*R* *T*(param->w-param->u))*((da_di[i]/param->Theta-db_di[i]/param->b)*log((V+param->c+param->u*param->b)/
                       (V+param->c+param->w*param->b))+((V+param->c)*db_di[i]-param->b*dc_di[i])*(param->u-param->w)/
                       ((V+param->c+param->u*param->b)*(V+param->c+param->w*param->b)))+(db_di[i]-dc_di[i])/(V+param->c-param->b);
            sumdArr_dXi=sumdArr_dXi+x[i]*dArr_dXi[i];
        }

        for (i=0;i<mix->numSubs;i++)
        {
            phi[i]=exp(dArr_dXi[i]+Arr-sumdArr_dXi+Z-1-log(Z));
            //printf("dArr_dXi:%f\n",dArr_dXi[i]);
            //printf("dArr_dNi:%f\n",dArr_dXi[i]-sumdArr_dXi);
            //printf("phi[%i]: %f\n",i,phi[i]);
        }
    }
    else{
        FF_VfromTPcubic(*T,*P,param,*option,resultL,resultG,&state);
        if(*option=='l'){
            V=resultL[0];
            Arr=resultL[1];
            Z=resultL[2];
        }
        else{
            V=resultG[0];
            Arr=resultG[1];
            Z=resultG[2];
        }

        //printf("V:%f Arr;%f Z;%f\n",V,Arr,Z);
        for (i=0;i<mix->numSubs;i++){
            phi[i]=exp(db_di[i]/param->b*(Z-1)-log(Z-B)- da_di[i] * log((Z+param->u*B)/(Z+param->w*B))/(param->u-param->w));
        }
        //printf("Theta:%f b:%f c:%f\n",param.Theta,param.b,param.c);
    }

}

//Mixture Ideal gas thermodynamic properties calculation, from a reference state, specified by T and P, where H and S are 0
void CALLCONV FF_MixIdealThermoEOS(const int *numSubs,const  FF_Correlation cp0[], const FF_BaseProp baseProp[], const double x[],double *refT,double *refP, FF_ThermoProperties *th0)
{
     FF_ThermoProperties th0S;//Ideal thermodynamic properties of the selected substance
    int i;
    th0S.T=th0->T;
    th0S.V=th0->V;
    th0->Cp=th0->H=th0->S=0;
    //printf("I am going for individual calculation. Num subs: %i\n",*numSubs);
    //printf("%i %f %f %f %f %f\n",equation[0],coef[0][0],coef[0][1],coef[0][2],coef[0][3],coef[0][4]);
    for (i=0;i<*numSubs;i++){
        th0S.MW=baseProp[i].MW;
        FF_IdealThermoEos(cp0[i].form,cp0[i].coef,*refT,*refP,&th0S);
        th0->Cp=th0->Cp+x[i]*th0S.Cp;
        th0->H=th0->H+x[i]*th0S.H;
        th0->S=th0->S+x[i]*(th0S.S-log(x[i]));
    }
    //printf("Ideal cp mix:%f\n",th0->Cp);
    th0->Cv=th0->Cp-R;
    th0->U=th0->H-R*(th0->T- *refT);//We need to substrat the integration of d(P*V)=d(R*T)=R*dT from reference T to actual T
    th0->A=th0->U-th0->T*th0->S;
    th0->G=th0->H-th0->T*th0->S;
}


//Mixture Residual thermodynamic properties calculation from T and V, using EOS
void CALLCONV FF_MixResidualThermoEOS(FF_MixData *mix,FF_PhaseThermoProp *thR)
{
    FF_CubicParam param;
    double delta,tau;
    switch (mix->eosType)
    {
        case FF_SAFTtype:
            FF_MixArrDerSAFT(&mix->mixRule,&thR->T,&thR->V,&mix->numSubs,mix->saftData,mix->intParam,thR->c,thR->ArrDer);
            break;
        /*case FF_SWtype:
            delta=1/(thR->V*((( FF_SWEOSdata*)data)->rhoRef));
            tau=((( FF_SWEOSdata*)data)->tRef)/thR->T;
            //printf("delta:%f  tau:%f\n",delta,tau);
            //FF_ArrDerSW(tau,delta,data,ArrDer);
            break;*/
        default://Cubic EOS
            switch(mix->mixRule){
            case FF_VdW:
            case FF_VdWnoInt:
            case FF_PR:
            case FF_MKP:
                //printf("Param calculation\n");
                FF_MixParamTderCubicEOS(&mix->mixRule,&thR->T,&mix->numSubs,mix->cubicData,mix->intParam,thR->c,&param);
                //printf("T: %f Theta*1e3: %f dTheta*1e3:%f d2theta*1e3: %f b:%f c:%f\n",thR->T,param.Theta,param.dTheta*1e3,param.d2Theta*1e3,param.b,param.c);
                break;
            default://gE mixing rules
                FF_MixParamTderCubicEOSgE(mix,&thR->T,thR->c,&param);
                break;
            }
        //else//gE mixing rules
            //printf("thR.T: %f thR.V: %f\n",thR->T,thR->V);
            FF_ArrDerCubic(thR->T,thR->V,&param,thR->ArrDer);
            break;
    }
    //printf("Arr:%f dArr/dV:%f d2Arr/dV2:%f dArr/dT:%f d2Arr/dT2:%f d2Arr/dVdT:%f\n",thR->ArrDer[0],thR->ArrDer[1],thR->ArrDer[2],
    //        thR->ArrDer[3],thR->ArrDer[4],thR->ArrDer[5]);
    if (mix->eosType==FF_SWtype)
    {
        /*thR->P=R*thR->T*(1+delta*ArrDer[1])/thR->V;
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=R*(tau*ArrDer[3]-ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        thR->H=thR->G+thR->T*thR->S;
        thR->dP_dV=-R*thR->T*(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])/thR->V/thR->V;
        thR->dP_dT=R*(1+delta*ArrDer[1]-delta*tau*ArrDer[5])/thR->V;
        thR->Cv=-R*tau*tau*ArrDer[4];
        thR->Cp=thR->Cv+R*pow((1+delta*ArrDer[1]-delta*tau*ArrDer[5]),2)/(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2]);*/

    }
    else
    {
        thR->P=R*thR->T*(1/thR->V-thR->ArrDer[1]);
        thR->A=R*thR->T*thR->ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=-R*(thR->T*thR->ArrDer[3]+thR->ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        //thR->U=-R*thR->T*thR->T*ArrDer[3];
        //printf("%f %f\n",thR->U,-R*thR->T*thR->T*ArrDer[3]);
        thR->H=thR->G+thR->T*thR->S;
        //thR->H=thR->U-R*thR->T*thR->V*ArrDer[1];
        //thR->H=-R*thR->T*(thR->V*ArrDer[1]+thR->T*ArrDer[3]);
        thR->Cv=-2*R*thR->T*thR->ArrDer[3]-R*thR->T*thR->T*thR->ArrDer[4];
        thR->dP_dV=R*thR->T*(-1/thR->V/thR->V-thR->ArrDer[2]);
        thR->dP_dT=R*(1/thR->V-thR->ArrDer[1]-thR->T*thR->ArrDer[5]);
        //thR->Cp=thR->Cv-thR->T*pow(R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]),2)/(R*thR->T*(-1/thR->V/thR->V-ArrDer[2]))-R;
        thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;

    }
}

//Mixture thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
void CALLCONV FF_MixThermoEOS(FF_MixData *mix,double *refT,double *refP, FF_PhaseThermoProp *th)
{
    //printf("Hi, I am FF_MixThermoEOS\n");
    int i;
    FF_ThermoProperties th0;
    FF_PhaseThermoProp thR;
    th->MW=0;
    for (i=0;i<mix->numSubs;i++) th->MW=th->MW+th->c[i]*mix->baseProp[i].MW;//MW should arrive already calculated
    for (i=0;i<mix->numSubs;i++)thR.c[i]=th->c[i];
    th0.MW=thR.MW=th->MW;//Perhaps not necessary
    th0.T=thR.T=th->T;
    th0.V=thR.V=th->V;
    //printf("Now going to calculate Ideal part\n");
    FF_MixIdealThermoEOS(&mix->numSubs,mix->cp0Corr,mix->baseProp,th->c,refT,refP,&th0);
    if (mix->eosType==FF_IdealType)
    {
        th->P=R*th->T/th->V;
        th->A=th0.A;
        th->G=th0.G;
        th->S=th0.S;
        th->U=th0.U;
        th->H=th0.H;
        th->dP_dV=-R*th->T/th->V/th->V;
        th->dP_dT=R/th->V;
        th->Cv=th0.Cv;
        th->Cp=th0.Cp;
    }
    else
    {
        //printf("Now going for residual part T:%f\n",thR.T);
        FF_MixResidualThermoEOS(mix,&thR);
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
    for(i=0;i<6;i++)th->ArrDer[i]=thR.ArrDer[i];
    //printf("MW:%f T:%f V:%f P:%f Cv:%f Cp:%f H:%f S:%f dP_dV:%f dP_dT:%f SS:%f JT:%f IT:%f\n",th->MW,th->T,th->V,th->P,th->Cv,th->Cp,th->H,th->S,th->dP_dV,th->dP_dT,th->SS,th->JT,th->IT);
    //printf("Ideal T:%f V:%f P:%f Cv0:%f Cp0:%f H0:%f S0:%f\n",th0.T,th0.V,th0.P,th0.Cv,th0.Cp,th0.H,th0.S);
    //printf("Residual T:%f V:%f P:%f Cvr:%f Cpr:%f Ur:%f Hr:%f Sr:%f\n",thR.T,thR.V,thR.P,thR.Cv,thR.Cp,thR.U,thR.H,thR.S);
}

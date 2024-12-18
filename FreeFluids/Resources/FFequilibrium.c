/*

 * FFequilibrium.c
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

// contains equilibrium calculations for mixtures
//===============================================
#include <math.h>
#include <stdio.h>
//#include <nlopt.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFeosMix.h"
#include "FFequilibrium.h"

//#include "FFtools.h"

//Calculation of the tangent plane distance for a test composition,regarding a base one. data.z in mol fraction and coef in mole number
double CALLCONV FF_TPD(unsigned nVar, const double coef[], double grad[], FF_FeedData *data){
    unsigned int i,equal;
    double xTotal,x[nVar],y[nVar],answerL[3],answerG[3],gA,gE;
    double substPhiA[nVar],substPhiE[nVar],tpd;
    char option,state;
    FF_SubsActivityData actData[nVar];
    xTotal=0;
    for(i=0;i<nVar;i++) xTotal=xTotal+coef[i];//count of the number of moles
    equal=1;
    for(i=0;i<nVar;i++){//We convert the mole number to fraction
        x[i]=coef[i]/xTotal;//Normalization of the test phase
        if(fabs(x[i]-data->z[i])>0.001) equal=0;//Detect that the two phases are of different composition
    }
    if(equal==1) return 1e5;
    else{
        gA=+HUGE_VALF;
        if(data->mix->thModelActEos!=1){//Not phi-phi
            FF_PhiAndActivity(data->mix,&data->T,&data->P,x,actData,substPhiA);
            gA=0;
            for(i=0;i<nVar;i++) gA=gA+x[i]*(log(x[i]*substPhiA[i]));
        }
        gE=+HUGE_VALF;
        if(data->mix->thModelActEos!=2){
            if(data->mix->thModelActEos==1){
                option='s';
                FF_MixVfromTPeos(data->mix,&data->T,&data->P,x,&option,answerL,answerG,&state);
                if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
                else if((state=='G')||(state=='g')) option='g';
                else return 0.0;
            }
            else option='g';
            FF_MixPhiEOS(data->mix,&data->T,&data->P,x,&option,substPhiE);
            gE=0;
            for(i=0;i<nVar;i++) gE=gE+log(data->z[i]*substPhiE[i]);
        }
        tpd=0;
        if(gA<gE){
            for(i=0;i<nVar;i++){
                tpd=tpd+x[i]*(log(x[i]*substPhiA[i])-data->logSubstFugacity[i]);
            }
        }
        else{
            for(i=0;i<nVar;i++){
                tpd=tpd+x[i]*(log(x[i]*substPhiE[i])-data->logSubstFugacity[i]);
            }
        }
        //printf("tpd found:%f\n",tpd);
        return tpd;
    }
}

//Calculation of the minimum tangent plane distance and the corresponding composition
void CALLCONV FF_StabilityCheck(FF_FeedData *data,double *tpd,double tpdX[]){
    unsigned nVar;//number of components
    char option,state;
    nVar=data->mix->numSubs;
    double answerL[3],answerG[3],substPhi[nVar],substPhiA[nVar],substPhiE[nVar],xTotal,gA,gE;
    //printf("nVar:%i\n",nVar);
    unsigned int i, code;
    FF_SubsActivityData actData[nVar];

    //Calculation of the reference Gibbs energy
    //printf("state:%c\n",state);
    //printf("Vl:%f Vg:%f\n",answerL[0],answerG[0]);
    gA=+HUGE_VALF;
    if(data->mix->thModelActEos!=1){//Not phi-phi. It is necessary to calculate activity
        FF_PhiAndActivity(data->mix,&data->T,&data->P,data->z,actData,substPhiA);
        gA=0;
        for(i=0;i<nVar;i++) gA=gA+data->z[i]*(log(data->z[i]*substPhiA[i]));
    }
    gE=+HUGE_VALF;
    if(data->mix->thModelActEos!=2){//Not gamma-gamma. It is necessary to calculate EOS fugacity
        if(data->mix->thModelActEos==1){//If phi-phi EOS is for liquid or gas phase
            option='s';
            FF_MixVfromTPeos(data->mix,&data->T,&data->P,data->z,&option,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
            else if((state=='G')||(state=='g')) option='g';
            else return;
        }
        else option='g';//EOS is only for gas phase
        FF_MixPhiEOS(data->mix,&data->T,&data->P,data->z,&option,substPhiE);
        gE=0;
        for(i=0;i<nVar;i++) gE=gE+log(data->z[i]*substPhiE[i]);
    }
    if(gA<gE){
        for(i=0;i<nVar;i++){
            substPhi[i]=substPhiA[i];
            data->logSubstFugacity[i]=log(data->z[i]*substPhiA[i]);
        }
    }
    else{
        for(i=0;i<nVar;i++){
            substPhi[i]=substPhiA[i];
            data->logSubstFugacity[i]=log(data->z[i]*substPhiE[i]);
        }
    }
    //printf("g feed, activity:%f EOS:%f\n",gA,gE);

    //Calculation of the tangent plane distance
    double k[nVar],xTest[nVar],substPhiTest[nVar],testTPD,testTPDA,testTPDE;
    int j;
    *tpd=+HUGE_VALF;
    //Loops of k calculation starting with original composition as gas phase
    for(j=0;j<8;j++){
        xTotal=0;
        for (i=0;i<nVar;i++){
            if(j==0){
                //Probably it is necessary a better initialization for LLE
                k[i]=exp(log(data->mix->baseProp[i].Pc/ data->P)+5.373*(1-data->mix->baseProp[i].w)*(1-data->mix->baseProp[i].Tc/ data->T));
                xTest[i]=data->z[i]/k[i];
            }
            else{
                xTest[i]=data->z[i]*substPhi[i]/substPhiTest[i];
            }
            xTotal=xTotal+xTest[i];
        }
        for(i=0;i<nVar;i++)xTest[i]=xTest[i]/xTotal;//normalization of xTest[i]

        testTPDA=+HUGE_VALF;
        if(data->mix->thModelActEos!=1){//Not phi-phi. It is necessary to calculate activity
            FF_PhiAndActivity(data->mix,&data->T,&data->P,xTest,actData,substPhiA);
            testTPDA=0;
            for(i=0;i<nVar;i++) testTPDA=testTPDA+xTest[i]*(log(xTest[i]*substPhiA[i])-data->logSubstFugacity[i]);
        }
        testTPDE=+HUGE_VALF;
        if(data->mix->thModelActEos!=2){//Not gamma-gamma. It is necessary to calculate EOS fugacity
            if(data->mix->thModelActEos==1){//If phi-phi EOS is for liquid or gas phase
                option='s';
                FF_MixVfromTPeos(data->mix,&data->T,&data->P,xTest,&option,answerL,answerG,&state);
                if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
                else if((state=='G')||(state=='g')) option='g';
                else return;
            }
            else option='g';//EOS is only for gas phase
            FF_MixPhiEOS(data->mix,&data->T,&data->P,xTest,&option,substPhiE);
            testTPDE=0;
            for(i=0;i<nVar;i++) testTPDE=testTPDE+xTest[i]*(log(xTest[i]*substPhiE[i])-data->logSubstFugacity[i]);
            }
        if(testTPDA<testTPDE){
            testTPD=testTPDA;
            for(i=0;i<nVar;i++) substPhiTest[i]=substPhiA[i];
        }
        else{
            testTPD=testTPDE;
            for(i=0;i<nVar;i++) substPhiTest[i]=substPhiE[i];
        }
        if(testTPD<*tpd){
            *tpd=testTPD;
            for(i=0;i<nVar;i++)tpdX[i]=xTest[i];
        }
    }
    //Now another loops as liquid phase   
    for(j=0;j<8;j++){
        xTotal=0;
        for (i=0;i<nVar;i++){
            if(j==0){
                //Probably it is necessary a better initialization for LLE
                k[i]=exp(log(data->mix->baseProp[i].Pc/ data->P)+5.373*(1-data->mix->baseProp[i].w)*(1-data->mix->baseProp[i].Tc/ data->T));
                xTest[i]=data->z[i]*k[i];
            }
            else{
                xTest[i]=data->z[i]*substPhi[i]/substPhiTest[i];
            }
            xTotal=xTotal+xTest[i];
        }
        for(i=0;i<nVar;i++)xTest[i]=xTest[i]/xTotal;//normalization of xTest[i]

        testTPDA=+HUGE_VALF;
        if(data->mix->thModelActEos!=1){//Not phi-phi. It is necessary to calculate activity
            FF_PhiAndActivity(data->mix,&data->T,&data->P,xTest,actData,substPhiA);
            testTPDA=0;
            for(i=0;i<nVar;i++) testTPDA=testTPDA+xTest[i]*(log(xTest[i]*substPhiA[i])-data->logSubstFugacity[i]);
        }
        testTPDE=+HUGE_VALF;
        if(data->mix->thModelActEos!=2){//Not gamma-gamma. It is necessary to calculate EOS fugacity
            if(data->mix->thModelActEos==1){//If phi-phi EOS is for liquid or gas phase
                option='s';
                FF_MixVfromTPeos(data->mix,&data->T,&data->P,xTest,&option,answerL,answerG,&state);
                if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
                else if((state=='G')||(state=='g')) option='g';
                else return;
            }
            else option='g';//EOS is only for gas phase
            FF_MixPhiEOS(data->mix,&data->T,&data->P,xTest,&option,substPhiE);
            testTPDE=0;
            for(i=0;i<nVar;i++) testTPDE=testTPDE+xTest[i]*(log(xTest[i]*substPhiE[i])-data->logSubstFugacity[i]);
            }
        if(testTPDA<testTPDE){
            testTPD=testTPDA;
            for(i=0;i<nVar;i++) substPhiTest[i]=substPhiA[i];
        }
        else{
            testTPD=testTPDE;
            for(i=0;i<nVar;i++) substPhiTest[i]=substPhiE[i];
        }
        if(testTPD<*tpd){
            *tpd=testTPD;
            for(i=0;i<nVar;i++)tpdX[i]=xTest[i];
        }
    }
}

//Calculation of the minimum tangent plane distance and the corresponding composition, by simulated annealing
void CALLCONV FF_StabilityCheckSA(FF_FeedData *data,double *tpd,double tpdX[]){
    int nVar,i,j,k,l,m;
    int na;//number of accepted points
    int ns=30;//number of cycles before step adjustement//20
    int nt=5;//number of iterations before temperature is reduced//4
    int ne=6;//number of temperature reductions//14
    float r;//randon value between -1 and 1
    float rm;//randon value between 0 and 1
    nVar=data->mix->numSubs;
    FF_SubsActivityData actData[nVar];

    double xTotal,answerL[3],answerG[3],dif,Vf,V,tpdLast,xLast[nVar],substPhiA[nVar],substPhiE[nVar],gA,gE,tpdA,tpdE;
    double tpdBest,xBest[nVar];
    double step;//step vector for next composition calculation
    double T0=5/(R*data->T);//Initial temperature//100
    double T;//working temperature of the optimization
    double accept,stepCorr,p1,p2;//p1 and p2 are the maximum and minimum accept rate desired
    char option,state;

    p1=0.3;
    p2=0.2;
    tpdLast=+HUGE_VALF;
    srand(time(NULL));
    T=T0;

    //Calculation of feed base Gibbs
    gA=+HUGE_VALF;
    if(data->mix->thModelActEos!=1){//Not phi-phi
        FF_PhiAndActivity(data->mix,&data->T,&data->P,data->z,actData,substPhiA);
        gA=0;
        for(i=0;i<nVar;i++){
            gA=gA+data->z[i]*(log(data->z[i]*substPhiA[i]));
        }
    }
    gE=+HUGE_VALF;
    if(data->mix->thModelActEos!=2) {//Not gamma-gamma
        option='s';
        FF_MixVfromTPeos(data->mix,&data->T,&data->P,data->z,&option,answerL,answerG,&state);
        if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
            option='l';
            Vf=answerL[0];
        }
        else if((state=='G')||(state=='g')){
            option='g';
            Vf=answerG[0];
        }
        else{
            *tpd=1e6;
            return;
        }
        FF_MixPhiEOS(data->mix,&data->T,&data->P,data->z,&option,substPhiE);
        gE=0;
        for(i=0;i<nVar;i++){
            gE=gE+data->z[i]*(log(data->z[i]*substPhiE[i]));
        }
    }
    if(gA<gE){
        for(i=0;i<nVar;i++){
            data->logSubstFugacity[i]=log(data->z[i]*substPhiA[i]);
        }
    }
    else{
        for(i=0;i<nVar;i++){
            data->logSubstFugacity[i]=log(data->z[i]*substPhiE[i]);
        }
    }
    //printf("g feed, activity:%f EOS:%f\n",gA,gE);

    //Initialization at feed composition, and the step vector
    for(i=0;i<nVar;i++){
        tpdX[i]=1.0/nVar;
    }
    step=1.0;
    for(m=0;m<ne;m++){
        for(l=0;l<nt;l++){//loops at fixed T
            na=0;
            for(k=0;k<ns;k++){//loops at fixed step vector
                r=(float) 2*rand()/RAND_MAX-1;//randon number between -1 and 1
                //printf("r:%f\n",r);
                for(j=0;j<nVar;j++){//loop along composition
                    for(i=0;i<nVar;i++) xLast[i]=tpdX[i];
                    if(r>0) tpdX[j]=tpdX[j]+r*step*(1-tpdX[j]);//New point to test
                    else tpdX[j]=tpdX[j]+r*step*tpdX[j];
                    xTotal=0;
                    for(i=0;i<nVar;i++) xTotal=xTotal+tpdX[i];//count of the number of moles
                    dif=0;
                    for(i=0;i<nVar;i++){
                        tpdX[i]=tpdX[i]/xTotal;//Normalization of concentrations
                        dif=dif+fabs(tpdX[i]-data->z[i]);
                    }
                    if(dif>0.01*nVar){
                        //Calculation of test tpd
                        tpdA=+HUGE_VALF;
                        tpdE=+HUGE_VALF;
                        if(data->mix->thModelActEos!=1){//Not phi-phi
                            FF_PhiAndActivity(data->mix,&data->T,&data->P,tpdX,actData,substPhiA);
                            tpdA=0;
                            for(i=0;i<nVar;i++){
                                tpdA=tpdA+tpdX[i]*(log(tpdX[i]*substPhiA[i])-data->logSubstFugacity[i]);
                            }
                        }
                        if(data->mix->thModelActEos!=2) {//Not gamma-gamma
                            option='s';
                            FF_MixVfromTPeos(data->mix,&data->T,&data->P,tpdX,&option,answerL,answerG,&state);
                            if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                                option='l';
                                V=answerL[0];
                            }
                            else if((state=='G')||(state=='g')){
                                option='g';
                                V=answerG[0];
                            }
                            //printf("state:%c V:%f\n",state,V);
                            FF_MixPhiEOS(data->mix,&data->T,&data->P,tpdX,&option,substPhiE);
                            tpdE=0;
                            for(i=0;i<nVar;i++){
                                tpdE=tpdE+tpdX[i]*(log(tpdX[i]*substPhiE[i])-data->logSubstFugacity[i]);
                            }
                        }
                        if(tpdA<tpdE) *tpd=tpdA;
                        else *tpd=tpdE;
                    }
                    else *tpd=+HUGE_VALF;
                    //printf("tpdX[0]:%f tpdX[1]:%f tpdX[2]:%f tpdA:%f tpdE:%f tpd:%f\n",tpdX[0],tpdX[1],tpdX[2],tpdA,tpdE,*tpd);
                    if(*tpd>=tpdLast){
                        rm=(float) rand()/RAND_MAX;//randon number between 0 and 1 for Metropolis criteria
                        if(exp((tpdLast-*tpd)/T)>rm){
                            //printf("R\n");
                            tpdLast=*tpd;
                            na=na+1;
                        }
                        else{
                            for(i=0;i<nVar;i++) tpdX[i]=xLast[i];
                            //printf("M rm:%f exp:%f tpd:%f\n",rm,exp((tpdLast-tpd)/T),*tpd);
                        }
                    }
                    else{
                        //printf("state:%c Va:%f Vb:%f\n",state,Va,Vb);
                        //printf("B\n");
                        tpdLast=*tpd;
                        if(*tpd<tpdBest){
                            tpdBest=*tpd;
                            for(i=0;i<nVar;i++)xBest[i]=tpdX[i];
                        }
                        na=na+1;
                    }
                }
            }
            accept=(double) na/(ns*nVar);
            //printf("accepted:%i tests:%i acceptance:%f\n",na,ns*nVar,accept);
            //printf("----------------------\n");
            if(accept>p1){
                stepCorr=1+2*(accept-p1)/p2;
                step=step*stepCorr;
                if (step>1) step=1;
                //for(i=0;i<nVar;i++) step[i]=step[i]*stepCorr;
                //printf("New step:%f stepCorr mult:%f\n",step,stepCorr);
            }
            else if(accept<p2){
                stepCorr=1+2*(p2-accept)/p2;
                step=step/stepCorr;
                if (step>1) step=1;
                //for(i=0;i<nVar;i++) step[i]=step[i]/stepCorr;
                //printf("New step:%f stepCorr divis:%f\n",step,stepCorr);
            }
        }

        if(m==ne-2){
            for(i=0;i<nVar;i++)tpdX[i]=xBest[i];
            T=1e-10;
            step=0.1*step;
        }
        else T=T*0.85;//0.85
        //printf("New T:%f\n",T);
        //printf("-------------------\n");
        //T=T*0.8;//0.8
        //printf("New T:%f\n",T);
        //printf("-------------------\n");
    }
    *tpd=tpdLast;
}

//Solves the Rachford-Rice equation from an initial k value. Returns phases concentrations, fugacity coefficients and gas fraction
void CALLCONV FF_RachfordRiceSolver(FF_MixData *mix,const double *T,const double *P,const double f[],const double kInit[],int *numRep,
                           double x[],double y[],double substPhiB[],double substPhiA[],double *a){
    FF_SubsActivityData actData[mix->numSubs];
    char option,state;
    int i,j,n=0,side=0;
    int nFails=0;//number of not successful iterations
    double answerL[3],answerG[3];
    double k[mix->numSubs];
    double xTotal,yTotal,aPrev;
    double xmin,xmax,fxmin,fxmax,fa;//minimum, maximum and test gas fractions, and its Rachford-Rice function results

    for(i=0;i<mix->numSubs;i++) k[i]=kInit[i];
    aPrev=-1;
    for(j=0;j<*numRep;j++){//We perform n times the solution of the Rachford-Rice equation
        //printf("loop: %i\n",j);
        xTotal=0;
        yTotal=0;
        for (i=0;i<mix->numSubs;i++){
            //printf("k[%i]: %f\n",i,k[i]);
            y[i]=f[i]*k[i];
            yTotal=yTotal+y[i];
        }
        //printf("yTotal: %f\n",yTotal);
        if(yTotal<=1){//Probably the mixture is at or below the bubble temperature
            for(i=0;i<mix->numSubs;i++){
                x[i]=f[i];
                y[i]=y[i]/yTotal;//normalization of y[i]
            }
            *a=0;
            nFails++;
        }
        else{
            for(i=0;i<mix->numSubs;i++){
                x[i]=f[i]/k[i];
                xTotal=xTotal+x[i];
            }
            //printf("xTotal: %f\n",xTotal);
            if(xTotal<=1){//probably the mixture is over the dew temperature
                for(i=0;i<mix->numSubs;i++){
                    y[i]=f[i];
                    x[i]=x[i]/xTotal;//normalization of x[i]
                }
                *a=1;
                nFails++;
            }
            else{//we need to solve the Rachford-Rice equation using k[i]
                xmin=0;
                xmax=1;
                fxmin=yTotal-1;//this must be positive
                fxmax=1-xTotal;//this must be negative
                fa=1;
                n=0;
                side=0;
                while ((n<20)&&(fabs(fa)>0.001)){//we begin to seek for the root
                    n=n+1;
                    *a=(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin);//This is the regula falsi method, Anderson-Bjork modified
                    fa=0;
                    for(i=0;i<mix->numSubs;i++) fa=fa+f[i]*((k[i]-1)/(1-*a*(1-k[i])));//Rachford-Rice equation
                    //printf("n:%i a:%f fa;%f\n",n,*a,fa);
                    //y=(*f)(*x);
                    if ((fa<0.001)&&(fa>-0.001)) break;//if we have arrived to the solution exit
                    if ((fxmax * fa)>0){//if the proposed solution is of the same sign than f(xmax)
                        if (side==-1){//if it happened also in the previous loop
                            if ((1-fa/fxmax)>0) fxmin *= (1-fa/fxmax);//we decrease f(xmin) for the next loop calculation
                            else fxmin *= 0.5;
                        }
                        xmax=*a;
                        fxmax=fa;
                        side = -1;//we register than the solution was of the same sign than previous xmax
                    }
                    else{//If the prosed solution is of the same sign than f(xmin) we apply the same technic
                        if (side==1){
                            if ((1-fa/fxmin)>0) fxmax *= (1-fa/fxmin);
                            else fxmax *= 0.5;
                        }
                        xmin=*a;
                        fxmin=fa;
                        side = 1;
                    }
                }
                //printf("j:%i Gas fraction:%f Sumatory:%f side:%i\n",j,*a,fa,side);
                //Now that we have the gas fraction we calculate the composition of the phases
                for(i=0;i<mix->numSubs;i++){
                    x[i]=f[i]/(1+*a*(k[i]-1));
                    y[i]=x[i]*k[i];
                    //printf("a:%f i:%i x:%f y:%f\n",*a,i,x[i],y[i]);
                }
            }
        }
        if(nFails>=3) break;//If we are not inside a plausible value of gas fraction after 3 tries
        //It follows the calculation of the new k[i]
        if(mix->thModelActEos!=1) FF_PhiAndActivity(mix,T,P,x,actData,substPhiB);
        else{
            option='s';
            FF_MixVfromTPeos(mix,T,P,x,&option,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                option='l';
            }
            else if((state=='G')||(state=='g')){
                option='g';
            }
            FF_MixPhiEOS(mix,T,P,x,&option,substPhiB);
        }
        if(mix->thModelActEos==2) FF_PhiAndActivity(mix,T,P,y,actData,substPhiA);
        else{
            option='s';
            FF_MixVfromTPeos(mix,T,P,y,&option,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                option='l';
            }
            else if((state=='G')||(state=='g')){
                option='g';
            }
            FF_MixPhiEOS(mix,T,P,y,&option,substPhiA);
        }

        for(i=0;i<mix->numSubs;i++) k[i]=substPhiB[i]/substPhiA[i];
        if((*a>0)&&(*a<1)&&(fabs(aPrev- *a)<0.0001)) break;
        else aPrev= *a;
    }
}


//Mixture bubble pressure calculation, given T, composition, and thermo model to use
void CALLCONV FF_BubblePtpdg(FF_MixData *mix,const double *T, const double x[],const double *bPguess, double *bP,double y[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *bP=0;
    double Pmax=400e5;//This is the maximum computable pressure
    double P;
    double k[mix->numSubs],newY[mix->numSubs],yTotal,dyTotal;
    int i=0,j,n=4,counter=0;
    char option,state;
    //If the liquid model is based only on activity all components must be below the critical temperature, and we need activity data
    if(mix->thModelActEos==0){
        for(i=0;i<mix->numSubs;i++){
            if(*T>mix->baseProp[i].Tc){
                *bP=0;
                return;
            }
        }
        //We get activity data
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
        //As BIP are for activity not use them in the eos
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;

        //now the approximation of molar fractions and P
        P=0;
        double Vp;
        for(i=0;i<mix->numSubs;i++){
            FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,*T,&Vp);
            //printf("Vp[%i]: form:%i %f\n",i,mix->vpCorr[i].form,Vp);
            P=P+x[i]*actData[i].gamma*Vp;//This is the initial P.
        }
        //printf("Approximate P from activity: %f\n",P);
        yTotal=0;
        for(i=0;i<mix->numSubs;i++){
            y[i]=x[i]*actData[i].gamma*Vp/P;//considering phi of gas=1 and not using Poynting factor
            yTotal=yTotal+y[i];
            //printf("y[%i] :%f\n",i,y[i]);
        }
        //printf("yTotal: %f\n",yTotal);
        //y[i] is already normalized
        if (*bPguess!=0) P=*bPguess;//if a P guess is given we take it

    }
    //As a comment: if cubic EOS, we can get the parameters for the liquid phase and use them in all the calculation. We can gain some speed here
    else if(mix->thModelActEos==1){  //If phi-phi
        if (*bPguess==0){
            P=Pmax*0.5;  //initialize P
            double yLow,yHigh,tpdg,tpdgMin,tpdgMax;
            yLow=0.64;//0.75;//0.65;
            yHigh=0.9999;
            tpdgMin=0.0005;//0.0005;
            tpdgMax=0.002;//0.002;
            do{  //bisection approximation using k values from liquid fugacity calculation
                //printf("P:%f\n",P);
                option='l';
                FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);  //phi of liquid phase
                for (i=0;i<mix->numSubs;i++){ //Initialization of phi of gas phase
                    substPhiG[i]=1.0;
                }
                for(j=0;j<5;j++){
                    yTotal=0;
                    for (i=0;i<mix->numSubs;i++){  //ks and concentrations calculation
                        //k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                        k[i]=substPhiL[i]/substPhiG[i];
                        y[i]=x[i]*k[i];
                        yTotal=yTotal+y[i];
                    }
                    for (i=0;i<mix->numSubs;i++){  //Normalization
                        y[i]=y[i]/yTotal;
                    }
                    option='g';  //Phi gas and tpd calculation
                    FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);
                    tpdg=0;
                    for(i=0;i<mix->numSubs;i++){
                        tpdg=tpdg+y[i]*(log(y[i])+log(substPhiG[i])-log(x[i])-log(substPhiL[i]));
                    }
                    //printf("x[0]:%f P:%f yTotal:%f y[0]:%f tpdg:%f\n",x[0],P,yTotal,y[0],tpdg);
                    if(yTotal<yLow) break; //Trying to avoid to find another liquid root that would bring yTotal=1
                    if (tpdg<=0) break;  //can save some calculations
                 }
                if((tpdg>tpdgMax)||(yTotal<yLow)){//probably we are too much over the bubble pressure
                    P=P-Pmax/n;
                    n=n*2;
                }
                else if((tpdg<tpdgMin)||(yTotal>yHigh)){
                    P=P+Pmax/n;
                    n=n*2;
                }
                counter++;
                //printf("x[0]:%f Counter:%i n:%i yTotal:%f New P:%f tpdg:%f\n",x[0],counter,n,yTotal,P,tpdg);
                if(counter>15) break;
            } while((tpdg<tpdgMin)||(tpdg>tpdgMax)||(yTotal<yLow)||(yTotal>yHigh));
            //y[i] could be not normalized
            //printf("x[0]:%f Counter:%i n:%i yTotal:%f P:%f tpdg:%f\n",x[0],counter,n,yTotal,P,tpdg);

        }
        else{ //P guessed and phi-phi
            P=*bPguess;
            option='l';
            FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);
            yTotal=0;
            for (i=0;i<mix->numSubs;i++){
                //y[i]=x[i]*exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                y[i]=x[i]*substPhiL[i];  //Alternative for phi. alternative for gamma needed
                yTotal=yTotal+y[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            //y[i] is not normalized
            //printf("P guess:%f yTotal:%f\n",P,yTotal);
            //If we got a guess from the user and we are using phi-phi
        }

        for (i=0;i<mix->numSubs;i++){
            y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
            //printf("initial composition: i[%i] y:%f yTotal: %f\n",i,y[i],yTotal);
        }
    }


    do{
        P = P *(1+yTotal)*0.5;
        if (mix->thModelActEos==0)FF_PhiFromActivity(mix,T,&P,x,actData,substPhiL);
        else{
            option='l';
            FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
        yTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<mix->numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        //printf("Counter:%i P:%f yTotal:%f\n",counter,P,yTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dyTotal > 0.001)&&(fabs(yTotal - 1) > 0.0015)){
            FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
            yTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            //printf("Adj.gas: P:%f substPhiG:%f,%f yTotal:%f\n",P,substPhiG[0],substPhiG[1],yTotal);
            for (i=0;i<mix->numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>100) return;
        }
        //printf("counter:%i P:%f yTotal:%f dyTotal:%f\n",counter,P,yTotal,dyTotal);
    } while (fabs(yTotal - 1) > 0.0002);
    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,T,&P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,T,&P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    //printf("P:%f Vl:%f Vdiff:%f Ydiff:%f\n",P,Vl,Vdiff,Ydiff);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        //printf("BubbleP:%f\n",P);
        *bP=P;
    }
}


//Mixture bubble pressure calculation, given T, composition, and thermo model to use. Fast but not as good as FF_BubbleP
void CALLCONV FF_BubbleP(FF_MixData *mix,const double *T, const double x[],const double *bPguess, double *bP,double y[],double substPhiL[],double substPhiG[]){
    int ver=0;
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *bP=0;
    double P,xTotal;
    double k[mix->numSubs],newY[mix->numSubs],yTotal,dyTotal;
    int i=0,counter=0;
    char option,state;

    if(mix->thModelActEos==0){//If the liquid model is based only on activity all components must be below the critical temperature, and we need activity data
        for(i=0;i<mix->numSubs;i++){
            if(*T>mix->baseProp[i].Tc){
                *bP=0;
                return;
            }
        }
        //We get activity data
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
        P=0;
        double Vp;
        for(i=0;i<mix->numSubs;i++){
            FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,*T,&Vp);
            //printf("Vp[%i]: form:%i %f\n",i,mix->vpCorr[i].form,Vp);
            P=P+x[i]*actData[i].gamma*Vp;//This is the initial P.
        }
        //printf("Approximate P from activity: %f\n",P);
        yTotal=0;
        for(i=0;i<mix->numSubs;i++){
            y[i]=x[i]*actData[i].gamma*Vp/P;//considering phi of gas=1 and not using Poynting factor
            yTotal=yTotal+y[i];
            //printf("y[%i] :%f\n",i,y[i]);
        }
        //printf("yTotal: %f\n",yTotal);
        //y[i] is already normalized
        if (*bPguess!=0) P=*bPguess;//if a P guess is given we take it
        //As BIP are for activity not use them in the eos
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;
    }
    //As a comment: if cubic EOS, we can get the parameters for the liquid phase and use them in all the calculation. We can gain some speed here

    else {//using and EOS for the liquid phase
        if(*bPguess==0){//We need to approximate the bubble pressure and gas composition
            //initialize P
            P=0.0;
            for (i=0;i<mix->numSubs;i++) P=P+x[i]*mix->baseProp[i].Pc;
            P=P*1.3;
            double Plow,Phigh,sumLow,sumHigh,Pnew;
            double yTarget;
            Plow=0;
            Phigh=0;
            xTotal=1.0;
            if(ver==11) printf("Initial P:%f \n",P);
            option='l';
            do{//secant approximation using k values from Wilson equation and/or liquid fugacity calculation
                if (P>60e5) yTarget=0.65;
                else if (P>40e5) yTarget=0.7;
                else yTarget=0.85;//0.775
                //yTarget=1.0;
                yTotal=0;
                FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//phi as liquid phase, but it could be gas at this pressure
                for (i=0;i<mix->numSubs;i++){
                    k[i]=substPhiL[i];
                    //k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                    y[i]=x[i]*k[i];
                    yTotal=yTotal+y[i];
                }
                if(yTotal>yTarget){
                    Plow=P;
                    sumHigh=yTotal-yTarget;
                    Pnew=P*1.2;
                }
                else if(yTotal<yTarget){
                    Phigh=P;
                    sumLow=yTotal-yTarget;
                    Pnew=P*0.8;
                }
                if (Phigh*Plow>0.0){
                    Pnew=(Phigh*sumHigh-Plow*sumLow)/(sumHigh-sumLow);
                }
                if (fabs(Pnew-P)<1.0) break;
                P=Pnew;
                counter++;
                if(ver==1) printf("Approach counter:%i xTotal:%f yTotal:%f New P:%f k[0]:%f k[1]:%f)\n",counter,xTotal,yTotal,P,k[0],k[1]);
                if(counter>25) break;
            } while(fabs(yTotal-yTarget)>0.00001);
            //y[i] is not necessarily normalized
        }

        else{   //If we got a bubble P guess from the user,
            option='l';
            FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//Speed up for cubics
            yTotal=0;
            for (i=0;i<mix->numSubs;i++){
                y[i]=x[i]*substPhiL[i];
                yTotal=yTotal+y[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            //y[i] is not normalized
        }

        for (i=0;i<mix->numSubs;i++){
            y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
            //printf("initial composition: i[%i] y:%f yTotal: %f\n",i,y[i],yTotal);
        }
    }

    if(ver==2) printf("Approximate bubble:%f xTotal:%f yTotal:%f approach counter:%i\n",P,xTotal,yTotal,counter);

    //now the main loop
    do{
        //P = P *(1+yTotal)*0.5;
        P=P*(1+(yTotal-1)*0.5);
        if (mix->thModelActEos==0)FF_PhiFromActivity(mix,T,&P,x,actData,substPhiL);
        else{
            option='l';
            FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
        yTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<mix->numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        if(ver==3) printf("Counter:%i P:%f yTotal:%f\n",counter,P,yTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dyTotal > 0.001)&&(fabs(yTotal - 1) > 0.00001)){
            FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
            yTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            if(ver==4) printf("Adj.gas: P:%f substPhiG:%f,%f yTotal:%f\n",P,substPhiG[0],substPhiG[1],yTotal);
            for (i=0;i<mix->numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>200) return;
        }
        //printf("counter:%i dyTotal:%f\n",counter,dyTotal);
    } while (fabs(yTotal - 1) > 0.0015);
    if (ver==5) printf("BubbleP:%f counter:%i x[0]:%f y[0]:%f\n",P,counter,x[0],y[0]);
    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,T,&P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,T,&P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl>0.02;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))){
        if (ver==6) printf("BubbleP:%f n:%i\n",P,counter);
        *bP=P;
    }
    //*bP=P;
}


//Mixture dew pressure calculation, given T, composition, and thermo model to use
void CALLCONV FF_DewPoriginal(FF_MixData *mix,const double *T, const double y[],const double *dPguess, double *dP,double x[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *dP=0;
    double Pmax=400e5;//This is the maximum computable pressure
    double P;
    double k[mix->numSubs],newX[mix->numSubs],xTotal,dxTotal;
    int i=0,j,n=4,counter=0;
    char option,state;
    //If the liquid model is based only on activity all components must be below the critical temperature, and we need activity data
    if(mix->thModelActEos==0){
        for(i=0;i<mix->numSubs;i++){
            if(*T>mix->baseProp[i].Tc) return;
        }
        //As BIP are for activity not use them in the eos
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;

    }

    if((mix->thModelActEos==1)&&(*dPguess==0)){  //If phi-phi without guess
        P=Pmax*0.5;  //initialize P
        double xLow,xHigh,tpdl;
        double answerL[3],answerG[3],VyL,VyG,VxL,VxG;
        double VyPrev,VxPrev,xTotalPrev,Pprev;
        char yPhase,xPhase;
        xLow=0.999;
        xHigh=1.0;
        option='l';
        FF_MixVfromTPeos(mix,T,&Pmax,y,&option,answerL,answerG,&state);
        VyPrev=answerL[0];//register the initial liquid volume
        VxPrev=answerL[0];
        Pprev=Pmax;
        yPhase='l';  //We begin saying that we are in liquid phase. We need to know when we change to have a gas root
        xPhase='l';
        xTotal=0;
        //printf("Initial V:%f\n",VyPrev);
        do{  //bisection approximation using k values from liquid fugacity calculation
            option='b';  //We ask for the calculations being made from liquid and gas sides (of course not necessary for cubics)
            FF_MixVfromTPeos(mix,T,&P,y,&option,answerL,answerG,&state);
            if (answerL[0]==0){  //If only one value is returned we need to work with it
                if(answerG[0]==0) return;
                else answerL[0]=answerG[0];
            }
            else if(answerG[0]==0) answerG[0]=answerL[0];
            VyL=answerL[0];
            VyG=answerG[0];
            if(VyL>(VyPrev*1.7)) yPhase='g';
            else if(VyL<(VyPrev/1.7)) yPhase='l';//Very important. Determination of the point where all roots are liquid
            if(fabs((VyG-VyL)/VyL)>0.05) yPhase='b';  //If we have two different roots we overwrite the dituation
            VyPrev=VyL;
            //printf("y[0]:%f P:%f VyL:%f VyG:%f yPhase:%c\n",y[0],P,VyL,VyG,yPhase);
            if(yPhase=='l'){  //still liquid
                P=P-Pmax/n;
                n=n*2;
                VxPrev=VyL;
            }
            else{  //there is a gas root
                option='l';
                FF_MixPhiEOS(mix,T,&P,y,&option,substPhiL);//phi as liquid just for initialization
                option='g';
                FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas for calculations
                xTotalPrev=0;
                for(j=0;j<10;j++){  //loop for improve the k values
                    xTotal=0;
                    for (i=0;i<mix->numSubs;i++){  //ks and concentrations calculation
                        if(j==0){  //You can't use the fugacity coef. from gas phase. In this case use Wilson. Perhaps is better to begin always with Wilson
                            if (yPhase=='g') k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                            else k[i]=substPhiL[i];
                            k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                        }
                        else k[i]=substPhiL[i]/substPhiG[i];
                        x[i]=y[i]/k[i];
                        xTotal=xTotal+x[i];
                    }
                    for (i=0;i<mix->numSubs;i++){  //Normalization
                        x[i]=x[i]/xTotal;
                    }
                    option='l';  //Phi liquid and tpd calculation
                    FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);
                    tpdl=0;
                    for(i=0;i<mix->numSubs;i++){
                        tpdl=tpdl+x[i]*(log(x[i])+log(substPhiL[i])-log(y[i])-log(substPhiG[i]));
                    }
                    if(fabs(xTotal-xTotalPrev)<0.0001) break;
                    xTotalPrev=xTotal;
                    //printf("xTotal:%f x[0]:%f tpdl:%f\n",xTotal,x[0],tpdl);
                }
                option='b';
                FF_MixVfromTPeos(mix,T,&P,x,&option,answerL,answerG,&state);
                if (answerL[0]==0){
                    if(answerG[0]==0) return;
                    else answerL[0]=answerG[0];
                }
                else if(answerG[0]==0) answerG[0]=answerL[0];
                VxL=answerL[0];
                VxG=answerG[0];
                if(VxL>VxPrev*1.7) xPhase='g';
                else if(VxL<VxPrev/1.7) xPhase='l';
                if(fabs((VxG-VxL)/VxL)>0.05) xPhase='b';
                VxPrev=VxL;
                if(xPhase=='g'){  //just a gas root for x. The condensate must be liquid!!
                    xTotal=0;
                    P=P+Pmax/n;
                    n=n*2;
                }
                else{  //everything is OK but it is necessary to adjust xTotal to 1
                    //if(((xPhase=='l')||(xPhase=='b'))&&(counter>14)) break;
                    if(xTotal>xHigh){
                        P=P-Pmax/n;
                        n=n*2;
                    }
                    else if(xTotal<xLow){
                        P=P+Pmax/n;
                        n=n*2;
                    }
                    else break;
                }
            }
            counter++;
            //printf("Counter:%i xTotal:%f x[0]:%f VyL:%f VyG:%f VxL:%f VxG:%f New P:%f\n",counter,xTotal,x[0],VyL,VyG,VxL,VxG,P);
            if(counter>20) break;
        } while((xTotal<xLow)||(xTotal>xHigh));//((tpdl<tpdlMin)||(tpdl>tpdlMax)||(xTotal<xLow)||(xTotal>xHigh))
    //printf("y[0]:%f Counter:%i n:%i xTotal:%f P:%f tpdl:%f\n",y[0],counter,n,xTotal,P,tpdl);
    }

    else{   //If we got a guess from the user or we are using gamma-phi
        if(*dPguess!=0){  //we got a guess
            P=*dPguess;
            xTotal=0;
            for (i=0;i<mix->numSubs;i++){
                x[i]=y[i]/exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                xTotal=xTotal+x[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            for (i=0;i<mix->numSubs;i++){
                x[i]=x[i]/xTotal; //  normalization
                //printf("Pguess:%f x[%i];%f\n",P,i,x[i]);
            }
            //printf("P guess:%f xTotal: %f\n",P,xTotal);
        }
        else{  //we are using gamma-phi
            double Vp[mix->numSubs];
            P=0;
            for(i=0;i<mix->numSubs;i++){
                if(mix->refVpEos==0) FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,*T,&Vp[i]);
                else{
                        if(mix->eosType==FF_SAFTtype) FF_VpEos(mix->eosType,*T,&mix->saftData[i],&Vp[i]);
                        else FF_VpEos(mix->eosType,*T,&mix->cubicData[i],&Vp[i]);
                    }
                P=P+y[i]/Vp[i];//Inverse of P assuming Raoult law and gas fugacity coef=1
                //printf("Activity Vp[%i]:%f\n",i,Vp[i]);
            }
            P=1/P;  //This will be the initial P
            for(i=0;i<mix->numSubs;i++){
                x[i]=y[i]*P/Vp[i];
                //printf("Activity x[%i]:%f\n",i,x[i]);
            }
            xTotal=1;
            //printf("Activity P:%f xTotal:%f\n",P,xTotal);
        }


        do{
            P =P /xTotal;
            if (mix->thModelActEos==0)FF_PhiAndActivity(mix,T,&P,x,actData,substPhiL);
            else{
                option='l';
                FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//phi as liquid phase
            }
            option='g';
            FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
            xTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                x[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+x[i];
            }
            for (i=0;i<mix->numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=1.0;
            //printf("Counter:%i new P:%f new xTotal:%f\n",counter,P,xTotal);
            //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiG:%f x:%f phiL:%f\n",i,substPhiG[i],x[i],substPhiL[i]);
            while ((dxTotal > 0.001)&&(fabs(xTotal - 1) > 0.001)){
                if(mix->thModelActEos==0) FF_PhiAndActivity(mix,T,&P,x,actData,substPhiL);
                else{
                    option='l';
                    FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);
                }

                xTotal=0;
                for (i=0;i<mix->numSubs;i++)
                {
                    newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                    xTotal=xTotal+newX[i];
                }
                for (i=0;i<mix->numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
                dxTotal=0;
                for (i=0;i<mix->numSubs;i++)
                {
                    dxTotal=dxTotal+fabs(newX[i]-x[i]);
                    x[i]=newX[i];
                }
                //printf("Adj.liquid: P:%f substPhiL:%f,%f xTotal:%f dxTotal:%f\n",P,substPhiL[0],substPhiL[1],xTotal,dxTotal);
                counter=counter+1;
                if (counter>100) return;
            }
            //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
        } while (fabs(xTotal - 1) > 0.0002);
    }


    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,T,&P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,T,&P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    //printf("P:%f Vl:%f Vdiff:%f Ydiff:%f\n",P,Vl,Vdiff,Ydiff);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        //printf("DewP:%f\n",P);
        *dP=P;
    }
}


//Mixture dew pressure calculation, given T, composition, and thermo model to use
void CALLCONV FF_DewP(FF_MixData *mix,const double *T, const double y[],const double *dPguess, double *dP,double x[],double substPhiL[],double substPhiG[]){
    int ver=0;
    FF_SubsActivityData actData[mix->numSubs];
    *dP=0;
    double Pmax=400e5;
    double P,yTotal=1.0;
    double k[mix->numSubs],newX[mix->numSubs],xTotal,dxTotal;
    int i=0,j,n=4,counter=0;
    char option,state;
    //If the liquid model is based only on activity all components must be below the critical temperature
    if(mix->thModelActEos==0){
        for(i=0;i<mix->numSubs;i++){
            if(*T>mix->baseProp[i].Tc){
                *dP=0;
                return;
            }
        }
        //As BIP are for activity not use them in the gas eos
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;
    }

    if (*dPguess==0) P=Pmax*0.5;
    else P=*dPguess;
    //now approximation by secant method
    //initialize P
    P=0.0;
    for (i=0;i<mix->numSubs;i++) P=P+y[i]*mix->baseProp[i].Pc;
    P=P*0.5;
    double Plow,Phigh,sumLow,sumHigh,Pnew;
    double xTarget=1.0;
    Plow=0;
    Phigh=0;

    do{//secant approximation using k values from Wilson equation and/or liquid fugacity calculation
        xTotal=0;
        for (i=0;i<mix->numSubs;i++){
            k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
            x[i]=y[i]/k[i];
            xTotal=xTotal+x[i];
        }
        if(xTotal>xTarget){
            Phigh=P;
            sumHigh=xTotal-xTarget;
            Pnew=P*0.8;
        }
        else if(xTotal<xTarget){
            Plow=P;
            sumLow=xTotal-xTarget;
            Pnew=P*1.2;
        }
        if (Phigh*Plow>0.0){
            Pnew=(Plow*sumHigh-Phigh*sumLow)/(sumHigh-sumLow);
        }
        if (fabs(Pnew-P)<1.0) break;
        P=Pnew;
        counter++;
        if(ver==1) printf("Approach counter:%i xTotal:%f yTotal:%f New P:%f k[0]:%f k[1]:%f)\n",counter,xTotal,yTotal,P,k[0],k[1]);
        if(counter>15) break;
    } while(fabs(xTotal-xTarget)>0.00001);

    for (i=0;i<mix->numSubs;i++){
        x[i]=x[i]/xTotal; //we normalize x in order to obtain molar fractions.
        //printf("initial composition: i[%i] x:%f xTotal: %f\n",i,x[i],xTotal);
    }

    //Now we have an initial P and can begin the calculation of the liquid and gas fugacity coefficients
    double answerL[3],answerG[3];
    if (mix->thModelActEos==0){
        FF_PhiAndActivity(mix,T,&P,x,actData,substPhiL);
        option='g';
        FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);
        //for(i=0;i<mix->numSubs;i++)printf("Initial. Subst[%i] x and y:%f phiL:%f phiG:%f\n",i,y[i],substPhiL[i],substPhiG[i]);
    }
    else {
        option='l';
        //FF_MixPhiEOS(mix,T,&P,y,&option,substPhiL);
        FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);
        option='g';
        FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
        //for(i=0;i<mix->numSubs;i++)printf("Initial. Subst[%i] x and y:%f phiL:%f phiG:%f\n",i,y[i],substPhiL[i],substPhiG[i]);
    }

    xTotal=0;
    for (i=0;i<mix->numSubs;i++){
        x[i]=y[i]*substPhiG[i]/substPhiL[i];//Initial composition of liquid phase
        xTotal=xTotal+x[i];
    }
    for (i=0;i<mix->numSubs;i++){
        x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        if(ver==1) printf("Enter in the loop with P:%f x[%i]: %f xTotal: %f\n",P,i,x[i],xTotal);
    }
    do{
        P =P /xTotal;
        if (mix->thModelActEos==0)FF_PhiFromActivity(mix,T,&P,x,actData,substPhiL);
        else{
            option='l';
            FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,T,&P,y,&option,substPhiG);//phi as gas phase
        xTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<mix->numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("Counter:%i new P:%f new xTotal:%f\n",counter,P,xTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiG:%f x:%f phiL:%f\n",i,substPhiG[i],x[i],substPhiL[i]);
        while ((dxTotal > 0.001)&&(fabs(xTotal - 1) > 0.001)){
            if(mix->thModelActEos==0) FF_PhiAndActivity(mix,T,&P,x,actData,substPhiL);
            else{
                option='l';
                FF_MixPhiEOS(mix,T,&P,x,&option,substPhiL);
            }

            xTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            for (i=0;i<mix->numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            //printf("Adj.liquid: P:%f substPhiL:%f,%f xTotal:%f dxTotal:%f\n",P,substPhiL[0],substPhiL[1],xTotal,dxTotal);
            counter=counter+1;
            if (counter>200) return;
        }
        //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
    } while (fabs(xTotal - 1) > 0.001);
    //check that the phases are different
    double Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,T,&P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,T,&P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    //printf("P:%f Vl:%f Vdiff:%f Ydiff:%f\n",P,Vl,Vdiff,Ydiff);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        //printf("DewP:%f\n",P);
        *dP=P;
    }
    //printf("BubbleP:%f n:%i\n",P,counter);
}



//Pressure envelope of a binary mixture
void CALLCONV FF_PressureEnvelope(FF_MixData *mix,const double *T, const int *nPoints, double c[],double bP[],double y[],double dP[],double x[]){
    if(mix->numSubs!=2) return;
    int i;
    double interval;
    double Vp;
    double in[2],out[2];
    double substPhiL[2];
    double substPhiG[2];
    double guess=0;
    interval=1.0 /(*nPoints - 1);
    for(i=0;i<*nPoints;i++){
        c[i]=i*interval;
    }
    if(mix->eosType==FF_SAFTtype)FF_VpEos(mix->eosType,*T,&mix->saftData[0],&Vp);
    else FF_VpEos(mix->eosType,*T,&mix->cubicData[0],&Vp);
    bP[*nPoints - 1]=dP[*nPoints - 1]=Vp;
    x[*nPoints - 1]=y[*nPoints - 1]=1;
    //printf("bP:%f\n",bP[*numPoints - 1]);
    for(i=*nPoints-2;i>0;i--){
        in[0]=c[i];
        in[1]=1-c[i];
        FF_BubbleP(mix,T,in,&guess,&bP[i],out,substPhiL,substPhiG);
        y[i]=out[0];
        FF_DewP(mix,T,in,&guess,&dP[i],out,substPhiL,substPhiG);
        x[i]=out[0];
        //printf("c:%f bP:%f y:%f dP:%f x:%f\n",in[0],bP[i],y[i],dP[i],x[i]);
    }
    if(mix->eosType==FF_SAFTtype)FF_VpEos(mix->eosType,*T,&mix->saftData[1],&Vp);
    else FF_VpEos(mix->eosType,*T,&mix->cubicData[1],&Vp);
    bP[0]=dP[0]=Vp;
    y[0]=x[0]=0;
}


//Mixture bubble temperature calculation, given P, composition, and thermo model to use
void CALLCONV FF_BubbleTtpdg(FF_MixData *mix,const double *P, const double x[],const double *bTguess, double *bT,double y[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *bT=0;
    double Tmax=800;//This is the maximum computable absolute temperature
    double T;//to use in the computations
    double k[mix->numSubs],newY[mix->numSubs],yTotal=0,dyTotal;
    int i=0,j,n=4,counter=0;
    char option,state;
    //If the liquid model is based on activity and Vp all components must be below the critical pressure
    if((mix->thModelActEos==0)&&(mix->refVpEos==0)){
        for(i=0;i<mix->numSubs;i++){
            if(*P>mix->baseProp[i].Pc) return;
        }
    }
    if(mix->thModelActEos==0){
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;//As BIP are for activity not use them in the gas pahse eos
    }

    if((mix->thModelActEos==1)&&(*bTguess==0)){  //If phi-phi without guess
        double yLow,yHigh,tpdg,tpdMin,tpdMax;
        T=Tmax*0.5;
        yLow=0.65;//0.65;
        yHigh=0.999;
        tpdMin=0.0005;
        tpdMax=0.002;
        do{//bisection approximation using k values from liquid fugacity calculation. To study to change to regula falsi
            //printf("P:%f\n",P);

            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
            for (i=0;i<mix->numSubs;i++){ //Initialization of phi of gas phase
                substPhiG[i]=1;
            }
            for(j=0;j<3;j++){
                yTotal=0;
                for (i=0;i<mix->numSubs;i++){  //Initial k and concentrations calculation
                    k[i]=substPhiL[i]/substPhiG[i];
                    //k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                    y[i]=x[i]*k[i];
                    yTotal=yTotal+y[i];
                }

                for (i=0;i<mix->numSubs;i++){  //Normalization
                    y[i]=y[i]/yTotal;
                }
                option='g';  //fugacity and tpd calculation
                FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);
                tpdg=0;
                for(i=0;i<mix->numSubs;i++){
                    tpdg=tpdg+y[i]*(log(y[i])+log(substPhiG[i])-log(x[i])-log(substPhiL[i]));
                }
                //printf("x[0]:%f T:%f yTotal:%f y[0]:%f tpdg:%f\n",x[0],T,yTotal,y[0],tpdg);
                if(yTotal<yLow) break; //Trying to avoid to find another liquid root that would bring yTotal=1
                //if(tpdg<0) break;  //doesn't work fine
            }
            if((tpdg>tpdMax)||(yTotal<yLow)){//probably we are too much below the bubble temperature
                T=T+Tmax/n;
                n=n*2;
            }
            else if((tpdg<tpdMin)||(yTotal>yHigh)){
                T=T-Tmax/n;
                n=n*2;
            }


            counter++;
            //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
            if(counter>15) break;
        } while((tpdg<tpdMin)||(tpdg>tpdMax)||(yTotal<yLow)||(yTotal>yHigh));
    }
    else{   //If we got a guess from the user or we are using gamma-phi
        if(*bTguess!=0){  //we got a guess
            T=*bTguess;
            yTotal=0;
            for (i=0;i<mix->numSubs;i++){
                y[i]=x[i]*exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                //y[i]=x[i]*substPhiL[i];  //Alternative for phi. alternative for gamma needed
                yTotal=yTotal+y[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            for (i=0;i<mix->numSubs;i++){
                y[i]=y[i]/yTotal; //  normalization
                //printf("Pguess:%f y[%i];%f\n",P,i,y[i]);
            }
        }
        else{  //gamma-phi model
            double Tb[mix->numSubs];
            T=0;
            for(i=0;i<mix->numSubs;i++){
                if(mix->eosType==FF_SAFTtype) FF_TbEos(mix->eosType,*P,&mix->saftData[i],&Tb[i]);
                else FF_TbEos(mix->eosType,*P,&mix->cubicData[i],&Tb[i]);
                T=T+x[i]*Tb[i];
            }
            //printf("Gamma-Phi T initial:%f\n",T);
            yTotal=0;
            for (i=0;i<mix->numSubs;i++){
                y[i]=x[i]*exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                //y[i]=x[i]*substPhiL[i];  //Alternative for phi. alternative for gamma needed
                yTotal=yTotal+y[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            for (i=0;i<mix->numSubs;i++){
                y[i]=y[i]/yTotal; //  normalization
                //printf("T init:%f y[%i];%f\n",T,i,y[i]);
            }

        }
    }

    //final loop to obtain bubble T
    do{
        T=T/pow(yTotal,0.08);//0.06
        if (mix->thModelActEos==0){
            FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
        }
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        yTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<mix->numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        //printf("Counter:%i T:%f yTotal:%f\n",counter,T,yTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dyTotal > 0.001)&&(fabs(yTotal - 1) > 0.0015)){
            FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
            yTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            //printf("Adj.gas: T:%f substPhiG:%f,%f yTotal:%f\n",T,substPhiG[0],substPhiG[1],yTotal);
            for (i=0;i<mix->numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>100) return;
        }
        //printf("counter:%i T:%f yTotal:%f\n",counter,T,yTotal);
        if(counter>100) break;
    } while (fabs(yTotal - 1) > 0.001);

    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl>0.02;
    //printf("T:%f Vl:%f Vdiff:%f Ydiff:%f\n",T,Vl,Vdiff,Ydiff);
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        *bT=T;
    }
    *bT=T;
}


//Mixture bubble temperature calculation, given P, composition, and thermo model to use
void CALLCONV FF_BubbleT1bona(FF_MixData *mix,const double *P, const double x[],const double *bTguess, double *bT,double y[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *bT=0;
    double Tmax=800;//This is the maximum computable absolute temperature
    double T,xTotal;//to use in the computations
    double k[mix->numSubs],newY[mix->numSubs],yTotal,dyTotal;
    int i=0,j,n=4,counter=0;
    char option,state;

    if(mix->thModelActEos==0){
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;//As BIP are for activity not use them in the gas pahse eos
    }
    if(*bTguess==0){  //We need to approximate the bubble temperature
        if((mix->thModelActEos==0)&&(mix->refVpEos==0)){//If activity and Vp are used for the liquid it is easier
            double Thigh=0,Tlow=0;
            double Vp,Pcalc=0,Phigh=0,Plow=0;;
            T=1000;
            for(i=0;i<mix->numSubs;i++){
                if(mix->baseProp[i].Tc<T) T=mix->baseProp[i].Tc;//We get the minimum Tc. Over this the method is not applicable
            }
            //We need to get the temperature using bisection
            //We get activity data
            switch(mix->actModel){
            case FF_UNIFACStd:
                FF_ActivityUNIFAC(&mix->unifStdData,&T,x,actData);
                break;
            case FF_UNIFACPSRK:
                FF_ActivityUNIFAC(&mix->unifPSRKData,&T,x,actData);
                break;
            case FF_UNIFACDort:
                FF_ActivityUNIFAC(&mix->unifDortData,&T,x,actData);
                break;
            case FF_UNIFACNist:
                FF_ActivityUNIFAC(&mix->unifNistData,&T,x,actData);
                break;
            default:
                FF_Activity(mix,&T,x,actData);
                break;
            }
            //and now vapor pressure and pressure
             for(i=0;i<mix->numSubs;i++){
                FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,T,&Vp);
                //printf("Vp[%i]: form:%i %f\n",i,mix->vpCorr[i].form,Vp);
                Pcalc=Pcalc+x[i]*actData[i].gamma*Vp;//This is the pressure at the maximum achivable temperature
            }
            if(*P>Pcalc) return;//there is no solution

            while(fabs(Pcalc-*P)/ *P>0.0001){
                if(Pcalc>*P){
                    Thigh=T;
                    Phigh=Pcalc;
                }
                else{
                    Tlow=T;
                    Plow=Pcalc;
                }
                T=(Thigh+Tlow)/2;
                //We get activity data
                switch(mix->actModel){
                case FF_UNIFACStd:
                    FF_ActivityUNIFAC(&mix->unifStdData,&T,x,actData);
                    break;
                case FF_UNIFACPSRK:
                    FF_ActivityUNIFAC(&mix->unifPSRKData,&T,x,actData);
                    break;
                case FF_UNIFACDort:
                    FF_ActivityUNIFAC(&mix->unifDortData,&T,x,actData);
                    break;
                case FF_UNIFACNist:
                    FF_ActivityUNIFAC(&mix->unifNistData,&T,x,actData);
                    break;
                default:
                    FF_Activity(mix,&T,x,actData);
                    break;
                }
                //and now vapor pressure and pressure
                 for(i=0;i<mix->numSubs;i++){
                    FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,T,&Vp);
                    //printf("Vp[%i]: form:%i %f\n",i,mix->vpCorr[i].form,Vp);
                    Pcalc=Pcalc+x[i]*actData[i].gamma*Vp;//This is the pressure at the maximum achivable temperature
                }
            }
            //printf("Approximate P from activity: %f\n",P);
            yTotal=0;
            for(i=0;i<mix->numSubs;i++){
                y[i]=x[i]*actData[i].gamma*Vp/ *P;//considering phi of gas=1 and not using Poynting factor
                yTotal=yTotal+y[i];
                //printf("y[%i] :%f\n",i,y[i]);
            }
            //printf("yTotal: %f\n",yTotal);
        }
        else{//If an EOS is used
            double yLow,yHigh,tpdg;
            T=Tmax*0.5;
            do{//bisection approximation using k values from liquid fugacity calculation
                //printf("P:%f\n",P);
                yLow=0.5;
                yHigh=0.85;
                xTotal=0;
                yTotal=0;
                option='l';
                FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
                for (i=0;i<mix->numSubs;i++){
                    k[i]=substPhiL[i];
                    //k[i]=exp(log(mix->baseProp[i].Pc/ P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
                    y[i]=x[i]/k[i];//Notation seems not OK, but we are obtaining liquid phase concentration
                    xTotal=xTotal+y[i];
                    y[i]=x[i]*k[i];
                    yTotal=yTotal+y[i];
                }
//xTotal<3.5*pow(yTotal,1.5)
                if((xTotal>90)||(yTotal<yLow)){//probably we are too much over the bubble pressure
                    T=T+Tmax/n;
                    n=n*2;
                }
                else if((xTotal<1.8)||(yTotal>yHigh)){
                    T=T-Tmax/n;
                    n=n*2;
                }

                counter++;
                //printf("Approach counter:%i n:%i xTotal:%f yTotal:%f New T:%f tpdg:%f\n",counter,n,xTotal,yTotal,T,tpdg);
                if(counter>15) break;
            } while((xTotal<1.8)||(xTotal>90)||(yTotal<yLow)||(yTotal>yHigh));
        }
    }
    else{   //If we got a guess from the user
        T=*bTguess;
        option='l';
        FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//Speed up for cubics
        yTotal=0;
        for (i=0;i<mix->numSubs;i++){
            y[i]=x[i]*substPhiL[i];
            yTotal=yTotal+y[i];
            //printf("k[%i]:%f\n",i,k[i]);
        }
    }


    for (i=0;i<mix->numSubs;i++){
        y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        //printf("initial composition: i[%i] y:%f yTotal: %f\n",i,y[i],yTotal);
    }

    //final loop to obtain bubble T
    do{
        T=T/pow(yTotal,0.06);//0.06 original
        if (mix->thModelActEos==0){
            switch(mix->actModel){
            case FF_UNIFACStd:
                FF_ActivityUNIFAC(&mix->unifStdData,&T,x,actData);
                break;
            case FF_UNIFACPSRK:
                FF_ActivityUNIFAC(&mix->unifPSRKData,&T,x,actData);
                break;
            case FF_UNIFACDort:
                FF_ActivityUNIFAC(&mix->unifDortData,&T,x,actData);
                break;
            case FF_UNIFACNist:
                FF_ActivityUNIFAC(&mix->unifNistData,&T,x,actData);
                break;
            default:
                FF_Activity(mix,&T,x,actData);
                break;
            }
            FF_PhiFromActivity(mix,&T,P,x,actData,substPhiL);
        }
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        yTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<mix->numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        //printf("Counter:%i T:%f yTotal:%f\n",counter,T,yTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dyTotal > 0.001)&&(fabs(yTotal - 1) > 0.0015)){
            FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
            yTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            //printf("Adj.gas: T:%f substPhiG:%f,%f yTotal:%f\n",T,substPhiG[0],substPhiG[1],yTotal);
            for (i=0;i<mix->numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>200) return;
        }
        //printf("counter:%i dyTotal:%f\n",counter,dyTotal);
    } while (fabs(yTotal - 1) > 0.0015);

    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl>0.02;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))){
        *bT=T;
    }
}


//Mixture bubble temperature calculation, given P, composition, and thermo model to use
void CALLCONV FF_BubbleT(FF_MixData *mix,const double *P, const double x[],const double *bTguess, double *bT,double y[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *bT=0;
    double Tl=0,Th=0,Tnew,yLow,yHigh;
    double T;//to use in the computations
    double k[mix->numSubs],newY[mix->numSubs],yTotal,dyTotal;
    int i=0,j,n=4,counter=0;
    char option,state;

    if(mix->thModelActEos==0){//gamma-phi model previous considerations
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;//As BIP are for activity not use them in the gas pahse eos
        //If the liquid model is based on activity and Vp all components must be below the critical pressure
        if(mix->refVpEos==0){
            for(i=0;i<mix->numSubs;i++){
                if(*P>mix->baseProp[i].Pc) return;
            }
        }
    }

    if(*bTguess==0){  // without guess
        double yTarget=1.0;  //0.68;
        T=0;
        for (i=0;i<mix->numSubs;i++) T=T+x[i]*mix->baseProp[i].Tc;
        T=0.7*T;
        option='l';
        do{
            //printf("P:%f\n",P);
            if (mix->thModelActEos==1){//phi-phi
                FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase, but it could be gas at this pressure
            }

            yTotal=0;
            for (i=0;i<mix->numSubs;i++){  //Initial k and concentrations calculation
                if (mix->thModelActEos==1) k[i]=1.47*substPhiL[i];//phi-phi
                else k[i]=(pow(mix->baseProp[i].Pc,((1/T-1/mix->baseProp[i].Tb)/(1/mix->baseProp[i].Tc-1/mix->baseProp[i].Tb)))/ *P);//gamma-phi. An alternative is to use gamma
                y[i]=x[i]*k[i];
                yTotal=yTotal+y[i];
            }
            if (yTotal<yTarget){
                Tl=T;
                yLow=yTotal-yTarget;
                Tnew=T*1.1;
            }
            else if(yTotal>yTarget){
                Th=T;
                yHigh=yTotal-yTarget;
                Tnew=T/1.1;
            }
            else break;
            if (Tl*Th>0) Tnew=(yHigh*Tl-yLow*Th)/(yHigh-yLow);//when we have already calculated Tl and Th
            if ((fabs(T-Tnew)<0.001)||(fabs(yTotal-yTarget)<0.00001)) break;
            else T=Tnew;
            counter++;
            //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
        } while(counter<10);

    }
    else{   //If we got a guess from the user
        T=*bTguess;
        if (mix->thModelActEos==1){//phi-phi
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase, but it could be gas at this pressure
        }
        yTotal=0;
        for (i=0;i<mix->numSubs;i++){
            if (mix->thModelActEos==1) k[i]=substPhiL[i];//phi-phi
            else k[i]=pow(mix->baseProp[i].Pc,((1/T-1/mix->baseProp[i].Tb)/(1/mix->baseProp[i].Tc-1/mix->baseProp[i].Tb)))/ *P;//gamma-phi.
            //An alternative is to use a calculated gamma
            y[i]=x[i]*k[i];
            yTotal=yTotal+y[i];
            //printf("k[%i]:%f\n",i,k[i]);
        }
    }

    //Now the final search
    for (i=0;i<mix->numSubs;i++){
        y[i]=y[i]/yTotal; //  normalization
        //printf("T init:%f y[%i];%f\n",T,i,y[i]);
    }
    Tl=0;
    Th=0;
    do{
        //printf("P:%f\n",P);
        //FF_MixPhiEOS(mix,&TbEst,P,x,&option,substPhiL);//phi as liquid phase, but it could be gas at this pressure
        if (mix->thModelActEos==0){
            FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
        }
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        yTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<mix->numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        //printf("Counter:%i T:%f yTotal:%f\n",counter,T,yTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dyTotal > 0.001)&&(fabs(yTotal - 1) > 0.00001)){
            FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
            yTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            //printf("Adj.gas: T:%f substPhiG:%f,%f yTotal:%f\n",T,substPhiG[0],substPhiG[1],yTotal);
            for (i=0;i<mix->numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>100) return;
        }

        if (yTotal<1){
            Tl=T;
            yLow=yTotal-1;
            Tnew=T*1.01;
        }
        else if(yTotal>1){
            Th=T;
            yHigh=yTotal-1;
            Tnew=T/1.01;
        }
        else break;
        if (Tl*Th>0) Tnew=(yHigh*Tl-yLow*Th)/(yHigh-yLow);//when we have already calculated Tl and Th
        if ((fabs(T-Tnew)<0.001)||(fabs(yTotal-1)<0.00001)) break;
        else T=Tnew;
        counter++;
        //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
    } while(counter<100);

    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl>0.02;
    //printf("T:%f Vl:%f Vdiff:%f Ydiff:%f\n",T,Vl,Vdiff,Ydiff);
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        *bT=T;
    }
}


//Mixture dew temperature calculation, given P, composition, and thermo model to use
void CALLCONV FF_DewTOriginal(FF_MixData *mix,const double *P, const double y[],const double *dTguess, double *dT,double x[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *dT=0;
    double Tmax=800;//This is the maximum computable temperature
    double T;
    double k[mix->numSubs],newX[mix->numSubs],xTotal,dxTotal;
    int i=0,j,n=4,counter=0;
    char option,state;
    //If the liquid model is based only on activity all components must be below the critical temperature, and we need activity data
    if((mix->thModelActEos==0)&&(mix->refVpEos==0)){
        for(i=0;i<mix->numSubs;i++){
            if(*P>mix->baseProp[i].Pc) return;
        }
        //As BIP are for activity not use them in the eos
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;

    }

    if((mix->thModelActEos==1)&&(*dTguess==0)){  //If phi-phi without guess
        T=Tmax*0.5;  //initialize P
        double xLow,xHigh,tpdl;
        double answerL[3],answerG[3],VyL,VyG,VxL,VxG;
        double VyPrev,VxPrev,xTotalPrev,Tprev;
        char yPhase,xPhase;
        xLow=0.999;
        xHigh=1.0;
        option='g';
        FF_MixVfromTPeos(mix,&Tmax,P,y,&option,answerL,answerG,&state);
        VyPrev=answerG[0];//register the initial liquid volume
        VxPrev=answerG[0];
        Tprev=Tmax;
        yPhase='g';  //We begin saying that we are in gas phase. We need to know when we change to have a liquid root
        xPhase='g';
        xTotal=0;
        //printf("Initial V:%f\n",VyPrev);
        do{  //bisection approximation using k values from liquid fugacity calculation
            option='b';  //We ask for the calculations being made from liquid and gas sides (of course not necessary for cubics)
            FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
            if (answerL[0]==0){  //If only one value is returned we need to work with it
                if(answerG[0]==0) return;
                else answerL[0]=answerG[0];
            }
            else if(answerG[0]==0) answerG[0]=answerL[0];
            VyL=answerL[0];
            VyG=answerG[0];
            if(VyL>(1.4*VyPrev*T/Tprev)) yPhase='g';
            else if(VyL<(0.57*VyPrev*T/Tprev)) yPhase='l';//Very important. Determination of the point where all roots are liquid
            if(fabs((VyG-VyL)/VyL)>0.05) yPhase='b';  //If we have two different roots we overwrite the dituation
            VyPrev=VyL;
            //printf("y[0]:%f T:%f VyL:%f VyG:%f yPhase:%c\n",y[0],T,VyL,VyG,yPhase);
            if(yPhase=='l'){  //necesary gas
                Tprev=T;
                T=T+Tmax/n;
                n=n*2;
                VxPrev=VyL;
            }
            else{  //there is a gas root
                option='l';
                FF_MixPhiEOS(mix,&T,P,y,&option,substPhiL);//phi as liquid just for initialization
                option='g';
                FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas for calculations
                xTotalPrev=0;
                for(j=0;j<10;j++){  //loop for improve the k values
                    xTotal=0;
                    for (i=0;i<mix->numSubs;i++){  //ks and concentrations calculation
                        if(j==0){  //You can't use the fugacity coef. from gas phase. In this case use Wilson. Perhaps is better to begin always with Wilson
                            k[i]=exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                            /*
                            if (yPhase=='g') k[i]=exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                            else k[i]=substPhiL[i];*/
                        }
                        else k[i]=substPhiL[i]/substPhiG[i];
                        x[i]=y[i]/k[i];
                        xTotal=xTotal+x[i];
                    }
                    for (i=0;i<mix->numSubs;i++){  //Normalization
                        x[i]=x[i]/xTotal;
                    }
                    option='l';  //Phi liquid and tpd calculation
                    FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);
                    tpdl=0;
                    for(i=0;i<mix->numSubs;i++){
                        tpdl=tpdl+x[i]*(log(x[i])+log(substPhiL[i])-log(y[i])-log(substPhiG[i]));
                    }
                    if(fabs(xTotal-xTotalPrev)<0.0001) break;
                    xTotalPrev=xTotal;
                    //printf("xTotal:%f x[0]:%f tpdl:%f\n",xTotal,x[0],tpdl);
                }
                option='b';
                FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
                if (answerL[0]==0){
                    if(answerG[0]==0) return;
                    else answerL[0]=answerG[0];
                }
                else if(answerG[0]==0) answerG[0]=answerL[0];
                VxL=answerL[0];
                VxG=answerG[0];
                if(VxL>1.4*VxPrev*T/Tprev) xPhase='g';
                else if(VxL<0.57*VxPrev*T/Tprev) xPhase='l';
                if(fabs((VxG-VxL)/VxL)>0.05) xPhase='b';
                VxPrev=VxL;
                if(xPhase=='g'){  //just a gas root for x. The condensate must be liquid!!
                    Tprev=T;
                    xTotal=0;
                    T=T-Tmax/n;
                    n=n*2;
                }
                else{  //everything is OK but it is necessary to adjust xTotal to 1
                    //if(((xPhase=='l')||(xPhase=='b'))&&(counter>14)) break;
                    if(xTotal>xHigh){
                        Tprev=T;
                        T=T+Tmax/n;
                        n=n*2;
                    }
                    else if(xTotal<xLow){
                        Tprev=T;
                        T=T-Tmax/n;
                        n=n*2;
                    }
                    else break;
                }
            }
            counter++;
            //printf("Counter:%i xTotal:%f x[0]:%f VyL:%f VyG:%f VxL:%f VxG:%f New T:%f\n",counter,xTotal,x[0],VyL,VyG,VxL,VxG,T);
            if(counter>20) break;
        } while((xTotal<xLow)||(xTotal>xHigh));//((tpdl<tpdlMin)||(tpdl>tpdlMax)||(xTotal<xLow)||(xTotal>xHigh))
    //printf("y[0]:%f Counter:%i n:%i xTotal:%f P:%f tpdl:%f\n",y[0],counter,n,xTotal,P,tpdl);
    }

    else{   //If we got a guess from the user or we are using gamma-phi
        if(*dTguess!=0){  //we got a guess
            T=*dTguess;
            xTotal=0;
            for (i=0;i<mix->numSubs;i++){
                x[i]=y[i]/exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                xTotal=xTotal+x[i];
                //printf("k[%i]:%f\n",i,k[i]);
            }
            for (i=0;i<mix->numSubs;i++){
                x[i]=x[i]/xTotal; //  normalization
                //printf("Pguess:%f x[%i];%f\n",P,i,x[i]);
            }
            //printf("P guess:%f xTotal: %f\n",P,xTotal);
        }
        else{  //we are using gamma-phi
            double Vp[mix->numSubs];
            P=0;
            for(i=0;i<mix->numSubs;i++){
                if(mix->refVpEos==0) FF_PhysPropCorrM(mix->vpCorr[i].form,mix->vpCorr[i].coef,mix->baseProp[i].MW,T,&Vp[i]);
                else{
                        if(mix->eosType==FF_SAFTtype) FF_VpEos(mix->eosType,T,&mix->saftData[i],&Vp[i]);
                        else FF_VpEos(mix->eosType,T,&mix->cubicData[i],&Vp[i]);
                    }
                T=T+y[i]/Vp[i];//Inverse of P assuming Raoult law and gas fugacity coef=1
                //printf("Activity Vp[%i]:%f\n",i,Vp[i]);
            }
            T=1/T;  //This will be the initial P
            for(i=0;i<mix->numSubs;i++){
                x[i]=y[i]*T/Vp[i];
                //printf("Activity x[%i]:%f\n",i,x[i]);
            }
            xTotal=1;
            //printf("Activity P:%f xTotal:%f\n",P,xTotal);
        }
    }
    do{
        T=T*pow(xTotal,0.06);
        if (mix->thModelActEos==0)FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        xTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<mix->numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("Counter:%i new P:%f new xTotal:%f\n",counter,P,xTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiG:%f x:%f phiL:%f\n",i,substPhiG[i],x[i],substPhiL[i]);
        while ((dxTotal > 0.001)&&(fabs(xTotal - 1) > 0.001)){
            if(mix->thModelActEos==0) FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
            else{
                option='l';
                FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);
            }

            xTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            for (i=0;i<mix->numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            //printf("Adj.liquid: P:%f substPhiL:%f,%f xTotal:%f dxTotal:%f\n",P,substPhiL[0],substPhiL[1],xTotal,dxTotal);
            counter=counter+1;
            if (counter>100) return;
        }
        //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
    } while (fabs(xTotal - 1) > 0.001);

    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    //printf("P:%f Vl:%f Vdiff:%f Ydiff:%f\n",P,Vl,Vdiff,Ydiff);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        //printf("DewP:%f\n",P);
        *dT=T;
    }
    //*dT=T;
}



//Mixture dew temperature calculation, given P, composition, and thermo model to use
void CALLCONV FF_DewT(FF_MixData *mix,const double *P, const double y[],const double *dTguess, double *dT,double x[],double substPhiL[],double substPhiG[]){
    FF_CubicParam param;
    FF_SubsActivityData actData[mix->numSubs];
    *dT=0;
    double Tl=0,Th=0,Tnew,xLow,xHigh;
    double T=0;//to use in the computations
    double k[mix->numSubs],newX[mix->numSubs],xTotal=0,dxTotal;
    int i=0,j,n=4,counter=0;
    char option,state;

    if(mix->thModelActEos==0){//gamma-phi model previous considerations
        if(mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;//As BIP are for activity not use them in the gas pahse eos
        //If the liquid model is based on activity and Vp all components must be below the critical pressure
        if(mix->refVpEos==0){
            for(i=0;i<mix->numSubs;i++){
                if(*P>mix->baseProp[i].Pc) return;
            }
        }
    }

    if(*dTguess==0){  // without guess, equal for phi-phi and gamma-phi
        double xTarget;  //1.9;//Entre 1.5 y 2.0. a mas presión mas alto
        if (*P>20e5) xTarget=2.0;
        else if (*P>10e5) xTarget=1.5;
        else xTarget=1.0;
        T=0;
        for (i=0;i<mix->numSubs;i++) T=T+y[i]*mix->baseProp[i].Tc;
        T=0.77*T;
        option='g';

        do{
            //printf("P:%f\n",P);
            //FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
            xTotal=0;
            for (i=0;i<mix->numSubs;i++){  //Initial k and concentrations calculation
                //x[i]=y[i]/(pow(0.000145*mix->baseProp[i].Pc,((1/(1.8*T)-1/(mix->baseProp[i].Tb*1.8))/(1/(mix->baseProp[i].Tc*1.8)-1/(mix->baseProp[i].Tb*1.8))))/ (0.000145 * *P));
                x[i]=y[i]/exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1+mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ T));
                if (T>mix->baseProp[i].Tc) x[i]=x[i];
                else x[i]=x[i];
                xTotal=xTotal+x[i];
            }
            if (xTotal<xTarget){
                Tl=T;
                xLow=xTotal-xTarget;
                Tnew=T/1.1;
            }
            else if(xTotal>xTarget){
                Th=T;
                xHigh=xTotal-xTarget;
                Tnew=T*1.1;
            }
            //else break;//xTotal==1
            if (Tl*Th>0) Tnew=Th-xLow*(Tl-Th)/(xHigh-xLow);//Tnew=(xHigh*Tl-xLow*Th)/(xHigh-xLow);
            if ((fabs(T-Tnew)<0.001)||(fabs(xTotal-xTarget)<0.00010)) break;
            else T=Tnew;
            //printf("Counter:%i xTotal:%f New T:%f \n",counter,xTotal,T);
            counter++;

        } while(counter<20);
    }
    else{   //If we got a guess from the user
        T=*dTguess;
        option='g';
        FF_MixPhiEOS(mix,&T,P,x,&option,substPhiG);//phi as gas phase
        xTotal=0;
        for (i=0;i<mix->numSubs;i++){
            x[i]=y[i]*substPhiG[i];
            xTotal=xTotal+x[i];
            //printf("k[%i]:%f\n",i,k[i]);
        }
    }

    //Now the final search
    for (i=0;i<mix->numSubs;i++){
        x[i]=x[i]/xTotal; //  normalization
        //printf("T init:%f y[%i];%f\n",T,i,y[i]);
    }
    Tl=0;
    Th=0;

    /*
    //alternative that works worse
    do{
        //printf("Initial T:%f\n",T);
        //FF_MixPhiEOS(mix,&TbEst,P,x,&option,substPhiL);//phi as liquid phase, but it could be gas at this pressure
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        if (mix->thModelActEos==0){
            FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
        }
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }

        xTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<mix->numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("Counter:%i T:%f xTotal:%f\n",counter,T,xTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiL:%f y:%f phiG:%f\n",i,substPhiL[i],y[i],substPhiG[i]);
        while ((dxTotal > 0.0001)&&(fabs(xTotal - 1) > 0.0001)){
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as gas phase
            xTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            //printf("Adj.gas: T:%f substPhiG:%f,%f xTotal:%f\n",T,substPhiG[0],substPhiG[1],xTotal);
            for (i=0;i<mix->numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            counter=counter+1;
            if (counter>500) return;
        }

        if (xTotal<1){
            Tl=T;
            xLow=xTotal-1;
            Tnew=T/1.01;
        }
        else if(xTotal>1){
            Th=T;
            xHigh=xTotal-1;
            Tnew=T*1.01;
        }
        else break;
        if (Tl*Th>0) Tnew=Th-xLow*(Tl-Th)/(xHigh-xLow);//Tnew=(xHigh*Tl-xLow*Th)/(xHigh-xLow);//when we have already calculated Tl and Th
        if ((fabs(T-Tnew)<0.001)||(fabs(xTotal-1)<0.0001)) break;
        else T=Tnew;
        counter++;
        //printf("Counter:%i n:%i xTotal:%f New T:%f \n",counter,n,xTotal,T);
    } while(counter<100);*/

    do{
        T=T*pow(xTotal,0.06);
        if (mix->thModelActEos==0)FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
        else{
            option='l';
            FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);//phi as liquid phase
        }
        option='g';
        FF_MixPhiEOS(mix,&T,P,y,&option,substPhiG);//phi as gas phase
        xTotal=0;
        for (i=0;i<mix->numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<mix->numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("Counter:%i new P:%f new xTotal:%f\n",counter,P,xTotal);
        //for(i=0;i<mix->numSubs;i++)printf("Subst[%i] phiG:%f x:%f phiL:%f\n",i,substPhiG[i],x[i],substPhiL[i]);
        while ((dxTotal > 0.001)&&(fabs(xTotal - 1) > 0.001)){
            if(mix->thModelActEos==0) FF_PhiAndActivity(mix,&T,P,x,actData,substPhiL);
            else{
                option='l';
                FF_MixPhiEOS(mix,&T,P,x,&option,substPhiL);
            }

            xTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            for (i=0;i<mix->numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<mix->numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            //printf("Adj.liquid: P:%f substPhiL:%f,%f xTotal:%f dxTotal:%f\n",P,substPhiL[0],substPhiL[1],xTotal,dxTotal);
            counter=counter+1;
            if (counter>100) return;
        }
        //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
    } while (fabs(xTotal - 1) > 0.001);

    //check that the phases are different
    double answerL[3],answerG[3],Vl,Vdiff,Ydiff=0;
    option='l';
    FF_MixVfromTPeos(mix,&T,P,x,&option,answerL,answerG,&state);
    Vl=answerL[0];
    option='g';
    FF_MixVfromTPeos(mix,&T,P,y,&option,answerL,answerG,&state);
    Vdiff=fabs(Vl-answerG[0])/Vl;
    for(i=0;i<mix->numSubs;i++) Ydiff=Ydiff+fabs(x[i]-y[i]);
    //printf("P:%f Vl:%f Vdiff:%f Ydiff:%f\n",P,Vl,Vdiff,Ydiff);
    if((Vdiff>0.1)||((Vdiff>0.03)&&(Ydiff>0.02))||((Vdiff>0.015)&&(Ydiff>0.08))||((Vdiff>0.005)&&(Ydiff>0.16))){
        //printf("DewP:%f\n",P);
        *dT=T;
    }
    //*dT=T;
}

//Temperature envelope of a binary mixture
void CALLCONV FF_TemperatureEnvelope(FF_MixData *mix,const double *P, const int *nPoints, double c[],double bT[],double y[],double dT[],double x[]){
    mix->numSubs=2;
    int i;
    double interval;
    double Tb;
    double in[2],out[2];
    double substPhiL[2];
    double substPhiG[2];
    double guess=0;
    interval=1.0 /(*nPoints - 1);
    for(i=0;i<*nPoints;i++){
        c[i]=i*interval;
    }
    if(mix->eosType==FF_SAFTtype)FF_TbEos(mix->eosType,*P,&mix->saftData[0],&Tb);
    else FF_TbEos(mix->eosType,*P,&mix->cubicData[0],&Tb);
    bT[*nPoints - 1]=dT[*nPoints - 1]=Tb;
    x[*nPoints - 1]=y[*nPoints - 1]=1;
    //printf("bP:%f\n",bP[*numPoints - 1]);
    for(i=*nPoints-2;i>0;i--){
        in[0]=c[i];
        in[1]=1-c[i];
        FF_BubbleT(mix,P,in,&guess,&bT[i],out,substPhiL,substPhiG);
        y[i]=out[0];
        FF_DewT(mix,P,in,&guess,&dT[i],out,substPhiL,substPhiG);
        x[i]=out[0];
        //printf("c:%f bP:%f y:%f dP:%f x:%f\n",in[0],bP[i],y[i],dP[i],x[i]);
    }
    if(mix->eosType==FF_SAFTtype)FF_TbEos(mix->eosType,*P,&mix->saftData[1],&Tb);
    else FF_TbEos(mix->eosType,*P,&mix->cubicData[1],&Tb);
    bT[0]=dT[0]=Tb;
    y[0]=x[0]=0;
}

//VL flash calculation, given T, P, feed composition, eos and mixing rule, using bubble and dew pressure as help
void CALLCONV FF_TwoPhasesPreFlashPT(FF_MixData *mix, const double *T,const double *P,const double f[],
                                     double x[],double y[],double substPhiL[],double substPhiG[],double *beta){
    int i,j,n;
    double guess=0,bP=0,dP=0, pAux;
    FF_BubbleP(mix,T,f,&guess,&bP,y,substPhiL,substPhiG);
    FF_DewP(mix,T,f,&guess,&dP,x,substPhiL,substPhiG);
    /*
    if ((mix->eosType!=FF_SAFTtype)&&((bP==0)||(dP==0))){
        FF_FeedData data;
        double Gr;
        data.mix=mix;
        data.P=*P;
        data.T=*T;
        for(i=0;i<mix->numSubs;i++) data.z[i]=f[i];
        FF_TwoPhasesFlashPTSA(&data,x,y,substPhiL,substPhiG,beta,&Gr);

        char option,state;
        double answerL[3],answerG[3],temp;
        option='s';
        FF_MixVfromTPeos(mix,T,P,x,&option,answerL,answerG,&state);
        if((state=='G')||(state=='g')){
            *beta=1-*beta;
            for (i=0;i<mix->numSubs;i++){
                temp=x[i];
                x[i]=y[i];
                y[i]=temp;
            }
        }

        return;
    }*/

    if (dP>bP){
        pAux=bP;
        bP=dP;
        dP=pAux;
    }

    if ((*P>=bP)&&(bP>0)){
        *beta=0;
        for(i=0;i<mix->numSubs;i++){
            x[i]=f[i];
            y[i]=0.0;
        }
    }
    else{
        if ((*P<=dP)||(bP==0)||(dP==0)){
            *beta=1;
            for(i=0;i<mix->numSubs;i++){
                y[i]=f[i];
                x[i]=0.0;
            }
        }
        else{
            double k[mix->numSubs];
            int numRep=10;
            double a, fa, aNew, aOld=1e6, aMin,aMax,faMin,faMax;
            double xTotal,yTotal;
            char option;
            for(i=0;i<mix->numSubs;i++){
                k[i]=(*P-dP)*(y[i]/f[i]-f[i]/x[i])/(bP-dP)+f[i]/x[i];//we initialize k
            }
            a=0.5;//and initialize a
            for(j=0;j<numRep;j++){//We will perform n times the solution of the Rachford-Rice equation
                aMin=0;
                aMax=0;
                faMin=0;
                faMax=0;
                n=0;
                do{ //find a for the given k
                  fa=0;
                  for(i=0;i<mix->numSubs;i++) fa=fa+(k[i]-1)*f[i]/(1+a*(k[i]-1));

                  if (fa>0){
                      aMax=a;
                      faMax=fa;
                      aNew=a*1.2;
                  }
                  else if (fa<0){
                      aMin=a;
                      faMin=fa;
                      aNew=a*0.8;
                  }
                  if (!((aMin*aMax)==0)) aNew=aMin-faMin*(aMax-aMin)/(faMax-faMin);
                  n++;
                  //printf("j:%i n:%i a:%f fa:%f \n",j,n,a,fa);//printf("a:%f fa:%f faMin:%f faMax:%f \n",a,fa,faMin,faMax);
                  if ((fabs(a-aNew)<0.0001)||(fabs(fa)<0.0001)) break;
                  a=aNew;

                }while (n<15);

                xTotal=0;
                yTotal=0;
                for(i=0;i<mix->numSubs;i++){
                    x[i]=f[i]/(1+a*(k[i]-1));
                    xTotal=xTotal+x[i];
                    y[i]=f[i]*k[i]/(1+a*(k[i]-1));
                    yTotal=yTotal+y[i];
                }
                for(i=0;i<mix->numSubs;i++){
                    x[i]=x[i]/xTotal;
                    y[i]=y[i]/xTotal;
                }
                option='l';
                FF_MixPhiEOS(mix,T,P,x,&option,substPhiL);
                option='g';
                FF_MixPhiEOS(mix,T,P,y,&option,substPhiG);

                if (fabs(aOld-a)<0.0001) break;
                aOld=a;
                for(i=0;i<mix->numSubs;i++) k[i]=substPhiL[i]/substPhiG[i];

            }
            *beta=a;
            //printf("x[0]:%f \n",x[0]);
        }
    }
}


//VL flash calculation, given T, P, feed composition, eos and mixing rule
void CALLCONV FF_TwoPhasesFlashPT(FF_MixData *mix,const double *T,const double *P,const double f[],
                           double x[],double y[],double substPhiL[],double substPhiG[],double *beta){//f is feed composition, beta is the gas fraction
    FF_SubsActivityData actData[mix->numSubs],actData2[mix->numSubs];
    double k[mix->numSubs],xTotal,yTotal,substPhiF[mix->numSubs],logFeedFugacity[mix->numSubs];
    char option,state;
    int i,j,numRep;
    double answerL[3],answerG[3];
    double kInit[mix->numSubs];//initial k
    double vL,vG;//volume of heavy and light trial phases
    double a;// gas fraction
    double tpdl,tpdg,tpd,tpdlmin,tpdgmin;//tangent plane distances of liquid and gas phases, total, and minimum of trials
    double substPhiTest[mix->numSubs];//to hold calculation values
    //printf("Arrived to flash\n");
    *beta=HUGE_VALF;

    if(mix->thModelActEos==2){//gamma-gamma. LLE. A modified Trebble(1989) method is used for k initialization.
        FF_PhiAndActivity(mix,T,P,f,actData,substPhiF);//one phase is assumed equal to feed
        yTotal=0;
        for(i=0;i<mix->numSubs;i++){
            y[i]=f[i]*actData[i].gamma;
            yTotal=yTotal+y[i];
        }
        for(i=0;i<mix->numSubs;i++)y[i]=y[i]/yTotal;//normalization
        FF_Activity(mix,T,y,actData2);
        for(i=0;i<mix->numSubs;i++){
            kInit[i]=actData[i].gamma/actData2[i].gamma;//initialization of k
            //printf("kInit[%i]:%f\n",i,kInit[i]);
            //if((kInit[i]>0.5)&&(kInit[i]<=1.0))kInit[i]=0.5;
            if((kInit[i]>1.0)&&(kInit[i]<3.0))kInit[i]=3.0;
            //printf("kInit[%i]:%f\n",i,kInit[i]);
        }
    }
    else if(mix->thModelActEos==0){//gamma-phi. Almost sure VLE
        //printf("Activity model\n");
        //Perhaps this alternative can be avoided and unify with the phi-phi approach for kInit, using the gas EOS
        if (mix->mixRule==FF_VdW) mix->mixRule=FF_VdWnoInt;//The BIP are for the activity model, not for the EOS
        FF_PhiAndActivity(mix,T,P,f,actData,substPhiL);//fugacity coef. of the feed as liquid
        option='g';
        FF_MixVfromTPeos(mix,T,P,f,&option,answerL,answerG,&state);
        tpdl=0;
        if(answerG[2]>0.3){//If we can grant a gas phase for the eos calculation
            FF_MixPhiEOS(mix,T,P,f,&option,substPhiG);
            for(i=0;i<mix->numSubs;i++){
                tpdl=tpdl+f[i]*(log(substPhiL[i])-log(substPhiG[i]));//checks if gas phase is inestable
                kInit[i]=substPhiL[i]/substPhiG[i];
            }
        }
        else for(i=0;i<mix->numSubs;i++){//If we do not find a gas phase with the eos we asume that phis are 1 in gas phase
            tpdl=tpdl+f[i]*log(substPhiL[i]);
            kInit[i]=substPhiL[i];
        }
        if(tpdl<0)for(i=0;i<mix->numSubs;i++) substPhiF[i]=substPhiL[i];//determine the more stable phase for the feed
        else for(i=0;i<mix->numSubs;i++) substPhiF[i]=substPhiG[i];
    }
    else{//phi-phi approach. Could be VLE or LLE
        //Calculate stable phase for the feed and its fugacity coefficients
        //printf("fugacity model\n");
        option='s';
        FF_MixVfromTPeos(mix,T,P,f,&option,answerL,answerG,&state);
        //printf("state:%c\n",state);
        if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
        else if((state=='G')||(state=='g')) option='g';
        else return;
        FF_MixPhiEOS(mix,T,P,f,&option,substPhiF);//Calculation of feed fugacity coef
        //printf("Feed state: %c \n",state);

        //initialize kInit[i] according to Wilson equation and check that the generated phases are not changed
        xTotal=0;
        yTotal=0;
        for (i=0;i<mix->numSubs;i++){
            kInit[i]=exp(log(mix->baseProp[i].Pc/ *P)+5.373*(1-mix->baseProp[i].w)*(1-mix->baseProp[i].Tc/ *T));
            x[i]=f[i]/kInit[i];
            xTotal=xTotal+x[i];
            y[i]=f[i]*kInit[i];
            yTotal=yTotal+y[i];
        }
        for(i=0;i<mix->numSubs;i++){
            x[i]=x[i]/xTotal;//normalization of x[i]
            y[i]=y[i]/yTotal;//normalization of y[i]
        }
        option='s';
        FF_MixVfromTPeos(mix,T,P,x,&option,answerL,answerG,&state);
        if((state=='L')||(state=='l')||(state=='U')||(state=='E')) vL=answerL[0];
        else vL=answerG[0];
        FF_MixVfromTPeos(mix,T,P,y,&option,answerL,answerG,&state);
        if((state=='G')||(state=='g')||(state=='U')||(state=='E')) vG=answerG[0];
        else vG=answerL[0];
        if(vL>vG) for(i=0;i<mix->numSubs;i++) kInit[i]=1/kInit[i];//If phases are changed we invert them. compare densities instead of volumes?
    }

    //Here begins calculation of the possible phase split
    //we calculate once the log of the fugacity/pressure for the feed
    for(i=0;i<mix->numSubs;i++)logFeedFugacity[i]=log(f[i]*substPhiF[i]);
    numRep=16;
    FF_RachfordRiceSolver(mix,T,P,f,kInit,&numRep,x,y,substPhiL,substPhiG,&a);

    //In case a good solution for beta has been found we check the unstability of the feed
    if((a>0)&&(a<1)){
        tpdl=0;
        tpdg=0;
        for(i=0;i<mix->numSubs;i++){
            tpdl=tpdl+x[i]*(log(x[i])+log(substPhiL[i])-logFeedFugacity[i]);
            tpdg=tpdg+y[i]*(log(y[i])+log(substPhiG[i])-logFeedFugacity[i]);
        }
        tpd=(1-a)*tpdl+a*tpdg;
        //printf("a:%f tpdl:%f tpdg:%f tpd:%f\n",a,tpdl,tpdg,tpd);
        if(tpd<0){
            *beta=a;
            return;
        }
    }

    //If feed seems stable we do the full stability analysis and, if necessary, return to the Rachford-Rice solver with better start k.
    else{
        tpdlmin=+HUGE_VALF;
        for(j=0;j<5;j++){
            xTotal=0;
            for (i=0;i<mix->numSubs;i++){
                if(j==0){
                    k[i]=kInit[i];
                    x[i]=f[i]/k[i];
                }
                else{
                    x[i]=f[i]*substPhiF[i]/substPhiTest[i];
                }
                xTotal=xTotal+x[i];
            }
            for(i=0;i<mix->numSubs;i++)x[i]=x[i]/xTotal;//normalization of x[i]
            option='s';
            FF_MixVfromTPeos(mix,T,P,x,&option,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
            else if((state=='G')||(state=='g')) option='g';
            FF_MixPhiEOS(mix,T,P,x,&option,substPhiTest);
            tpdl=0;
            for(i=0;i<mix->numSubs;i++) tpdl=tpdl+x[i]*(log(x[i]*substPhiTest[i])-logFeedFugacity[i]);
            //printf("Second round tpdl:%f\n",tpdl);
            if(tpdl<tpdlmin){
                tpdlmin=tpdl;
                for(i=0;i<mix->numSubs;i++)substPhiL[i]=substPhiTest[i];
            }
        }
        tpdgmin=+HUGE_VALF;
        for(j=0;j<5;j++){
            yTotal=0;
            for (i=0;i<mix->numSubs;i++){
                if(j==0){
                    k[i]=kInit[i];
                    y[i]=f[i]*k[i];
                }
                else{
                    y[i]=f[i]*substPhiF[i]/substPhiG[i];
                }
                yTotal=yTotal+y[i];
            }
            for(i=0;i<mix->numSubs;i++)y[i]=y[i]/yTotal;//normalization of y[i]
            option='s';
            FF_MixVfromTPeos(mix,T,P,y,&option,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')) option='l';
            else if((state=='G')||(state=='g')) option='g';
            FF_MixPhiEOS(mix,T,P,y,&option,substPhiTest);
            tpdg=0;
            for(i=0;i<mix->numSubs;i++) tpdg=tpdg+y[i]*(log(y[i]*substPhiTest[i])-logFeedFugacity[i]);
            //printf("Second round tpdg:%f\n",tpdg);
            if(tpdg<tpdgmin){
                tpdgmin=tpdg;
                for(i=0;i<mix->numSubs;i++)substPhiG[i]=substPhiTest[i];
            }
        }
        if((tpdlmin<0)||(tpdgmin<0)){
            for(i=0;i<mix->numSubs;i++) kInit[i]=substPhiL[i]/substPhiG[i];
            FF_RachfordRiceSolver(mix,T,P,f,kInit,&numRep,x,y,substPhiL,substPhiG,&a);
            *beta=a;
            return;
        }
        else{
            for(i=0;i<mix->numSubs;i++){
                x[i]=y[i]=f[i];
            }
        }
    }
}

//Mixture 2 phases flash, given P,T, composition, and thermo model to use. By simulated annealing global minimization of the reduced Gibbs energy
void CALLCONV FF_TwoPhasesFlashPTSA(FF_FeedData *data, double x[],double y[],double substPhiB[],double substPhiA[],double *beta,double *Gr){
    //printf("2 phases Gibbs min flash\n\n");
    int nVar,i,j,k,l,m;
    int na;//number of accepted points
    int ns=30;//number of cycles before step adjustement//20
    int nt=5;//number of iterations before temperature is reduced//5
    int ne=6;//number of temperature reductions//4
    float r;//randon value between -1 and 1
    float rm;//randon value between 0 and 1
    nVar=data->mix->numSubs;
    FF_SubsActivityData actDataB[nVar];
    FF_SubsActivityData actDataA[nVar];
    double a[nVar],aTotal,bTotal,answerL[3],answerG[3],dif,Va,Vb,GrLast,aLast,GrBest,aBest[nVar];
    double step;//step vector for next composition calculation
    double T0=5/(R*data->T);//Initial temperature
    double T;//working temperature of the optimization
    double accept,stepCorr,p1,p2;//p1 and p2 are the maximum and minimum accept rate desired
    char optionA,optionB,state;
    p1=0.6;//Maximum probability of acceptance wanted//0.6
    p2=0.4;//Minimum probability of acceptance wanted//0.4
    GrLast=1e6;
    GrBest=1e6;
    srand(time(NULL));
    T=T0;
    for(i=0;i<nVar;i++){//Initialization of both phases at feed composition, and the step vector
        a[i]=0.3*data->z[i];
        x[i]=data->z[i];
        y[i]=data->z[i];
    }
    step=1.0;
    for(m=0;m<ne;m++){
        for(l=0;l<nt;l++){//loops at fixed T
            na=0;
            for(k=0;k<ns;k++){//loops at fixed step vector
                r=(float) 2*rand()/RAND_MAX-1;//randon number between -1 and 1
                //printf("r:%f\n",r);
                for(j=0;j<nVar;j++){//loop along composition
                    aLast=a[j];
                    //a[j]=r*step[j]+a[j]*(1-step[j]/data->z[j]);//alternative with randon [0,1]
                    if(r>0) a[j]=a[j]+r*step*(data->z[j]-a[j]);//New point to test
                    else a[j]=a[j]+r*step*a[j];
                    aTotal=0;
                    for(i=0;i<nVar;i++) aTotal=aTotal+a[i];//count of the number of moles in the fraction 1
                    bTotal=1-aTotal;//number of moles in fraction 2
                    dif=0;
                    for(i=0;i<nVar;i++){//We convert the mole number to fraction
                        x[i]=(data->z[i]-a[i])/bTotal;//Normalization of the liquid phase to test
                        y[i]=a[i]/aTotal;
                        dif=dif+fabs(x[i]-y[i]);
                    }
                    //Calculation of V for the phases computed by EOS
                    if(data->mix->thModelActEos!=2){//Si no se usa gamma-gamma (LLE), la fase A se calcula mediante eos
                        optionA='s';
                        FF_MixVfromTPeos(data->mix,&data->T,&data->P,y,&optionA,answerL,answerG,&state);
                        if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                            optionA='l';
                            Va=answerL[0];
                        }
                        else if((state=='G')||(state=='g')){
                            optionA='g';
                            Va=answerG[0];
                        }
                        else *Gr=1e6;

                    }
                    if(data->mix->thModelActEos==1){//phi-phi. La fase B se calcula mediante EOS
                        optionB='s';
                        FF_MixVfromTPeos(data->mix,&data->T,&data->P,x,&optionB,answerL,answerG,&state);
                        if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                            optionB='l';
                            Vb=answerL[0];
                        }
                        else if((state=='G')||(state=='g')){
                            optionB='g';
                            Vb=answerG[0];
                        }
                        else *Gr=1e6;
                    }

                    if((dif>0.02*nVar)||((data->mix->thModelActEos==1)&&(fabs(Va-Vb)/Va)>0.05)){
                        if(data->mix->thModelActEos==2){//LLE mediante gamma-gamma. No hay fase gas
                            //printf("1 actibity\n");
                            FF_PhiAndActivity(data->mix,&data->T,&data->P,y,actDataA,substPhiA);
                        }
                        else FF_MixPhiEOS(data->mix,&data->T,&data->P,y,&optionA,substPhiA);

                        if(data->mix->thModelActEos==1) FF_MixPhiEOS(data->mix,&data->T,&data->P,x,&optionB,substPhiB);
                        else{//gamma-phi o gamma-gamma. Calculo mediante actividad
                            //printf("2 actibity\n");
                            FF_PhiAndActivity(data->mix,&data->T,&data->P,x,actDataB,substPhiB);
                        }
                        *Gr=0;//Residual Gibbs energy
                        for(i=0;i<nVar;i++){
                            *Gr=*Gr+bTotal*x[i]*(log(x[i]*substPhiB[i]))+aTotal*y[i]*(log(y[i]*substPhiA[i]));
                        }
                       //printf("Gas fract.:%f x[0]:%f x[1]:%f y[0]:%f y[1]:%f Gr:%f\n",bTotal,x[0],x[1],y[0],y[1],*Gr);
                    }
                    else *Gr=1e6;
                    if(*Gr>=GrLast){
                        rm=(float) rand()/RAND_MAX;//randon number between 0 and 1 for Metropolis criteria
                        if(exp((GrLast-*Gr)/T)>rm){
                            //printf("R rm:%f exp:%f Gas fract.:%f a[0]:%f a[1]:%f GrLast:%f Gr:%f\n",rm,exp((GrLast-*Gr)/T),bTotal,a[0],a[1],GrLast,*Gr);
                            //printf("R\n");
                            GrLast=*Gr;
                            na=na+1;
                        }
                        else{
                            a[j]=aLast;
                            //printf("M rm:%f exp:%f Gr:%f\n",rm,exp((GrLast-Gr)/T),Gr);
                        }
                    }
                    else{
                        //printf("state:%c Va:%f Vb:%f\n",state,Va,Vb);
                        //printf("B Afract:%f Bfract:%f TOTAL:%f x[0]:%f y[0]:%f Gr:%f\n",aTotal,bTotal,aTotal+bTotal,x[0],y[0],*Gr);
                        //printf("B\n");
                        GrLast=*Gr;
                        if(*Gr<GrBest){
                            GrBest=*Gr;
                            for(i=0;i<nVar;i++)aBest[i]=a[i];
                        }
                        na=na+1;
                    }
                }
            }
            accept=(double) na/(ns*nVar);
            //printf("\n accepted:%i tests:%i acceptance:%f\n",na,ns*nVar,accept);
            //printf("----------------------\n");
            if(accept>p1){
                stepCorr=1+2*(accept-p1)/p2;
                step=step*stepCorr;
                if (step>1) step=1;
                //printf("New step:%f stepCorr mult:%f\n",step,stepCorr);
            }
            else if(accept<p2){
                stepCorr=1+2*(p2-accept)/p2;
                step=step/stepCorr;
                if (step>1) step=1;
                //printf("New step:%f stepCorr div:%f\n",step,stepCorr);
            }
        }
        if(m==ne-2){
            for(i=0;i<nVar;i++)a[i]=aBest[i];
            T=1e-10;
            step=0.1*step;
        }
        else T=T*0.85;//0.85
        //printf("New T:%f\n",T);
        //printf("-------------------\n");
    }
    *beta=aTotal;

}

//Mixture 2 phases flash, given P,T, composition, and thermo model to use. By differential evolution global minimization of the reduced Gibbs energy
void CALLCONV FF_TwoPhasesFlashPTDE(FF_FeedData *data, double x[],double y[],double substPhiB[],double substPhiA[],double *beta,double *Gr){
    //printf("2 phases Gibbs min flash\n\n");
    int nVar,i,j,k,l,m;
    int nPop;//number of population used
    int best;//index of the best solution
    float rm;//randon value between 0 and 1
    int rm1,rm2,rm3;//integer random values
    nVar=data->mix->numSubs;
    FF_SubsActivityData actDataB[nVar];
    FF_SubsActivityData actDataA[nVar];
    double xTest[nVar],yTest[nVar],aTotal,bTotal,answerL[3],answerG[3],dif,Va,Vb,gBest;
    double F,Cr;
    double auxG,auxComp[nVar];
    char optionA,optionB,state;
    if(nVar<4) nPop=pow(5,nVar);
    else nPop=20*nVar;
    double actPop[nPop][nVar],test[nVar],g[nPop],gTest;//actual points and new points, with its number of moles composition
    srand(time(NULL));
    F=0.5;
    Cr=0.5;
    gBest=+HUGE_VALF;
    if(nVar==2){//Creation of an uniform grid of points
        for(i=0;i<5;i++){
            for(j=0;j<5;j++){
                    actPop[i*5+j][0]=i*data->z[0]/4;
                    actPop[i*5+j][1]=j*data->z[1]/4;
                }
            }
        }
    else if(nVar==3){
        for(i=0;i<5;i++){
            for(j=0;j<5;j++){
                for(k=0;k<5;k++){
                    actPop[i*25+j*5+k][0]=i*data->z[0]/4;
                    actPop[i*25+j*5+k][1]=j*data->z[1]/4;
                    actPop[i*25+j*5+k][2]=k*data->z[2]/4;
                }
            }
        }
    }
    else{
        for(i=0;i<nPop;i++){//Initial random population generation, with its mole content
            for(j=0;j<nVar;j++){
                rm=(float) rand()/RAND_MAX;//randon number between 0 and 1
                actPop[i][j]=rm*data->z[j];
            }
        }
    }
    //for(j=0;j<nPop;j++) printf("pop[%i] n[0]:%f n[1]:%f\n",j,actPop[j][0],actPop[j][1]);

    //Calculation of Gibbs energy for the initial population
    for(j=0;j<nPop;j++){
        aTotal=0;
        for(i=0;i<nVar;i++) aTotal=aTotal+actPop[j][i];//count of the number of moles in the fraction A of element j
        bTotal=1-aTotal;//number of moles in fraction 2
        dif=0;
        for(i=0;i<nVar;i++){//We convert the mole number to fraction
            yTest[i]=actPop[j][i]/aTotal;//Normalization of the liquid phase to test
            xTest[i]=(data->z[i]-actPop[j][i])/bTotal;
            dif=dif+fabs(x[i]-y[i]);
        }
        if(data->mix->thModelActEos!=2){//Si no se usa gamma-gamma (LLE), la fase A se calcula mediante eos
            optionA='s';
            FF_MixVfromTPeos(data->mix,&data->T,&data->P,yTest,&optionA,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                optionA='l';
                Va=answerL[0];
            }
            else if((state=='G')||(state=='g')){
                optionA='g';
                Va=answerG[0];
            }
            else g[j]=1e6;

        }
        if(data->mix->thModelActEos==1){//phi-phi. La fase B se calcula mediante EOS
            optionB='s';
            FF_MixVfromTPeos(data->mix,&data->T,&data->P,xTest,&optionB,answerL,answerG,&state);
            if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                optionB='l';
                Vb=answerL[0];
            }
            else if((state=='G')||(state=='g')){
                optionB='g';
                Vb=answerG[0];
            }
            else g[j]=1e6;
        }

        if((dif>0.02*nVar)||((data->mix->thModelActEos==1)&&(fabs(Va-Vb)/Va)>0.05)){
            if(data->mix->thModelActEos==2){//LLE mediante gamma-gamma. No hay fase gas
                //printf("2 gamma-gamma \n");
                FF_PhiAndActivity(data->mix,&data->T,&data->P,yTest,actDataA,substPhiA);
            }
            else FF_MixPhiEOS(data->mix,&data->T,&data->P,yTest,&optionA,substPhiA);

            if(data->mix->thModelActEos==1) FF_MixPhiEOS(data->mix,&data->T,&data->P,xTest,&optionB,substPhiB);
            else{//gamma-phi o gamma-gamma. Calculo mediante actividad
                //printf("2 actibity\n");
                FF_PhiAndActivity(data->mix,&data->T,&data->P,xTest,actDataB,substPhiB);
            }
            g[j]=0;//Residual Gibbs energy
            for(i=0;i<nVar;i++){
                g[j]=g[j]+bTotal*xTest[i]*(log(xTest[i]*substPhiB[i]))+aTotal*yTest[i]*(log(yTest[i]*substPhiA[i]));
            }
        }
        //printf("actPop:%i Gas fract.:%f y[0]:%f y[1]:%f x[0]:%f x[1]:%f Gr:%f\n",j,aTotal,yTest[0],yTest[1],xTest[0],xTest[1],g[j]);
        if(g[j]<gBest){
            best=j;
            gBest=g[j];
        }
    }
    //printf("best element:%i g:%f\n",best,gBest);

    //Initialization of the evolution
    for(k=0;k<400;k++){
        for(j=0;j<nPop;j++){
            rm1=nPop*rand()/RAND_MAX;//randon number
            rm2=nPop*rand()/RAND_MAX;
            rm3=nPop*rand()/RAND_MAX;
            for(i=0;i<nVar;i++){//Crossover plus mutation
                rm=(float) rand()/RAND_MAX;//randon number between 0 and 1
                if(rm>Cr) test[i]=actPop[j][i];
                else{
                    test[i]=actPop[rm1][i]+F*(actPop[rm2][i]-actPop[rm3][i]);
                    if(test[i]<0) test[i]=0;
                    else if(test[i]>data->z[i]) test[i]=data->z[i];
                }
            }
            //Calculation og g of the test
            aTotal=0;
            for(i=0;i<nVar;i++) aTotal=aTotal+test[i];//count of the number of moles in the fraction A of the test element
            bTotal=1-aTotal;//number of moles in fraction 2
            dif=0;
            for(i=0;i<nVar;i++){//We convert the mole number to fraction
                yTest[i]=test[i]/aTotal;//Normalization of the liquid phase to test
                xTest[i]=(data->z[i]-test[i])/bTotal;
                dif=dif+fabs(xTest[i]-yTest[i]);
            }
            if(data->mix->thModelActEos!=2){//Si no se usa gamma-gamma (LLE), la fase A se calcula mediante eos
                optionA='s';
                FF_MixVfromTPeos(data->mix,&data->T,&data->P,yTest,&optionA,answerL,answerG,&state);
                if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                    optionA='l';
                    Va=answerL[0];
                }
                else if((state=='G')||(state=='g')){
                    optionA='g';
                    Va=answerG[0];
                }
                else gTest=1e6;
            }
            if(data->mix->thModelActEos==1){//phi-phi. Ambas fases se calculan mediante EOS
                optionB='s';
                FF_MixVfromTPeos(data->mix,&data->T,&data->P,xTest,&optionB,answerL,answerG,&state);
                if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                    optionB='l';
                    Vb=answerL[0];
                }
                else if((state=='G')||(state=='g')){
                    optionB='g';
                    Vb=answerG[0];
                }
                else gTest=1e6;
            }

            if((dif>0.02*nVar)||((data->mix->thModelActEos==1)&&(fabs(Va-Vb)/Va)>0.05)){
                if(data->mix->thModelActEos==2){//LLE mediante gamma-gamma. No hay fase gas
                    //printf("1 actibity\n");
                    FF_PhiAndActivity(data->mix,&data->T,&data->P,yTest,actDataA,substPhiA);
                }
                else FF_MixPhiEOS(data->mix,&data->T,&data->P,yTest,&optionA,substPhiA);

                if(data->mix->thModelActEos==1) FF_MixPhiEOS(data->mix,&data->T,&data->P,xTest,&optionB,substPhiB);
                else{//gamma-phi o gamma-gamma. Calculo mediante actividad
                    //printf("2 actibity\n");
                    FF_PhiAndActivity(data->mix,&data->T,&data->P,xTest,actDataB,substPhiB);
                }
                gTest=0;//Residual Gibbs energy
                for(i=0;i<nVar;i++){
                    gTest=gTest+bTotal*xTest[i]*(log(xTest[i]*substPhiB[i]))+aTotal*yTest[i]*(log(yTest[i]*substPhiA[i]));
                }
            }
            if(gTest<g[j]){
                for(i=0;i<nVar;i++)actPop[j][i]=test[i];
                g[j]=gTest;
                if(g[j]<gBest){
                    best=j;
                    gBest=g[j];
                    *beta=aTotal;
                    for(i=0;i<nVar;i++){
                        x[i]=xTest[i];
                        y[i]=yTest[i];
                    }
                    //printf("generation:%i best element:%i g:%f\n",k,best,gBest);
                }
            }
        }
        if(nPop>40){
            if((k==4)||(k==16)||(k==32)||(k==64)){//Population reduction after n evolutions
                for(i=0;i<nPop;i++){
                    for(l=i+1;l<nPop;l++){
                        if(g[i]>g[l]){
                            auxG=g[i];
                            g[i]=g[l];
                            g[l]=auxG;
                            for(m=0;m<nVar;m++){
                                auxComp[m]=actPop[i][m];
                                actPop[i][m]=actPop[l][m];
                                actPop[l][m]=auxComp[m];
                            }
                        }
                    }
                }
                nPop=nPop/2;
            }
        }

    }
    *Gr=gBest;

}

//VL flash calculation, given h, P, feed composition, eos and mixing rule, using bubble and dew temperature as help
void CALLCONV FF_TwoPhasesPreFlashP_HS(FF_MixData *mix, char election, const double *E,const double *P,const double f[], double *T,
                                     double x[],double y[],double substPhiL[],double substPhiG[],double *beta){
    int i;
    double guess=0,bT=0,dT=0, tAux;
    double answerL[3],answerG[3];
    char option,state;
    FF_PhaseThermoProp th;
    th.P=*P;
    double Tl=0,Th=0,Tnew,Ed,Eb,Et,Elow,Ehigh;
    int counter=0;

    FF_DewT(mix,P,f,&guess,&dT,x,substPhiL,substPhiG);
    FF_BubbleT(mix,P,f,&guess,&bT,y,substPhiL,substPhiG);

    if (dT<bT){
        tAux=bT;
        bT=dT;
        dT=tAux;
    }

    option='g';
    FF_MixVfromTPeos(mix,&dT,P,f,&option,answerL,answerG,&state);
    th.T=dT;
    th.V=answerG[0];
    for (i=0;i<mix->numSubs;i++)th.c[i]=f[i];
    FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
    if (election=='h')  Ed=th.H;
    else if (election=='s') Ed=th.S;
    //printf("H:%f dew T:%f dew V:%f dew H:%f  \n",*E,th.T,th.V,th.H);
    if (*E>=Ed){
        option='g';
        *beta=1.0;
        for(i=0;i<mix->numSubs;i++){
            y[i]=f[i];
            x[i]=0.0;
        }
        *T=dT;
        Et=Ed;
        do{
            if (*E>Et){
                Tl=*T;
                Elow=Et;
                Tnew=*T*1.08;
            }
            else if(*E<Et){
                Th=*T;
                Ehigh=Et;
                Tnew=*T/1.08;
            }
            else break;
            if (Tl*Th>0) Tnew=Tl+(*E-Elow)*(Th-Tl)/(Ehigh-Elow);//Tnew=(Ehigh*Tl-Elow*Th)/(Ehigh-Elow);
            if ((fabs(*T-Tnew)<0.001)||(fabs(Et-*E)<0.01)) break;
            else *T=Tnew;
            FF_MixVfromTPeos(mix,T,P,f,&option,answerL,answerG,&state);
            th.T=*T;
            th.V=answerG[0];
            FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
            if (election=='h')  Et=th.H;
            else if (election=='s')Et=th.S;
            counter++;
            //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
        } while(counter<15);
    }
    else {
        option='l';
        FF_MixVfromTPeos(mix,&bT,P,f,&option,answerL,answerG,&state);
        th.T=bT;
        th.V=answerL[0];
        for (i=0;i<mix->numSubs;i++)th.c[i]=f[i];
        FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
        //printf("H:%f bub T:%f bub V:%f bub H:%f  \n",*E,th.T,th.V,th.H);
        if (election=='h') Eb=th.H;
        else if (election=='s') Eb=th.S;
        if (*E<=Eb){
            option='l';
            *beta=0.0;
            for(i=0;i<mix->numSubs;i++){
                x[i]=f[i];
                y[i]=0.0;
            }
            *T=bT;
            Et=Eb;
            do{
                if (*E>Et){
                    Tl=*T;
                    Elow=Et;
                    Tnew=*T*1.08;
                }
                else if(*E<Et){
                    Th=*T;
                    Ehigh=Et;
                    Tnew=*T/1.08;
                }
                else break;
                if (Tl*Th>0) Tnew= Tl+(*E-Elow)*(Th-Tl)/(Ehigh-Elow);//(Ehigh*Tl-Elow*Th)/(Ehigh-Elow);
                if ((fabs(*T-Tnew)<0.001)||(fabs(Et-*E)<1)) break;
                else *T=Tnew;
                FF_MixVfromTPeos(mix,T,P,f,&option,answerL,answerG,&state);
                th.T=*T;
                th.V=answerL[0];
                FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
                if (election=='h')  Et=th.H;
                else if (election=='s') Et=th.S;
                counter++;
                //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
            } while(counter<15);
        }
        else{
            double Eliq,Egas;
            *T=bT;//just in order to give some value
            Et=Eb;
            Th=dT;
            Ehigh=Ed;
            Tl=bT;
            Elow=Eb;
            do{ if (*E>Et){
                    Tl=*T;
                    Elow=Et;
                    //Tnew=*T*1.1;
                }
                else if(*E<Et){
                    Th=*T;
                    Ehigh=Et;
                    //Tnew=*T/1.1;
                }
                if (Tl*Th>0) Tnew=Tl+(Th-Tl)*(*E-Elow)/(Ehigh-Elow);//(Ehigh*Tl-Elow*Th)/(Ehigh-Elow);//(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin)
                if ((fabs(*T-Tnew)<0.001)||(fabs(Et-*E)<1.0)) break;
                else *T=Tnew;
                FF_TwoPhasesPreFlashPT(mix,T,P,f,x,y,substPhiL,substPhiG,beta);
                option='l';
                FF_MixVfromTPeos(mix,T,P,x,&option,answerL,answerG,&state);
                th.T=*T;
                th.V=answerL[0];
                for (i=0;i<mix->numSubs;i++)th.c[i]=x[i];
                FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
                if (election=='h')  Eliq=th.H;
                else if (election=='s') Eliq=th.S;
                option='g';
                FF_MixVfromTPeos(mix,T,P,y,&option,answerL,answerG,&state);
                th.T=*T;
                th.V=answerG[0];
                for (i=0;i<mix->numSubs;i++)th.c[i]=y[i];
                FF_MixThermoEOS(mix,&mix->refT,&mix->refP,&th);
                if (election=='h')  Egas=th.H;
                else if (election=='s') Egas=th.S;
                Et = *beta * Egas +(1.0- *beta)*Eliq;
                counter++;
                //printf("H:%f T:%f Et:%f  \n",*E,*T,Et);
                //printf("Counter:%i n:%i yTotal:%f New T:%f tpdg:%f\n",counter,n,yTotal,T,tpdg);
            } while(counter<15);

        }
    }
}

//VL flash calculation, given gas fraction, P, feed composition, eos and mixing rule, using bubble and dew temperature as help
void CALLCONV FF_TwoPhasesPreFlashPTheta(FF_MixData *mix, const char *opt, const double *Theta,const double *P,const double f[], double *T,
                                     double x[],double y[],double substPhiL[],double substPhiG[]){
    double lMW, gMW;
    double guess=0,bT=0,dT=0, tAux;
    double Tl=0,Th=0,Tnew,ThetaT,ThetaLow,ThetaHigh;
    int counter=0;

    FF_DewT(mix,P,f,&guess,&dT,x,substPhiL,substPhiG);
    FF_BubbleT(mix,P,f,&guess,&bT,y,substPhiL,substPhiG);

    if (dT<bT){
        tAux=bT;
        bT=dT;
        dT=tAux;
    }

    *T=bT;//just in order to give some value
    ThetaT=0;
    Th=dT;
    ThetaHigh=1;
    Tl=bT;
    ThetaLow=0;
    do{ if (*Theta>ThetaT){
            Tl=*T;
            ThetaLow=ThetaT;
        }
        else if(*Theta<ThetaT){
            Th=*T;
            ThetaHigh=ThetaT;
        }
        if (Tl*Th>0) Tnew=Tl+(Th-Tl)*(*Theta-ThetaLow)/(ThetaHigh-ThetaLow);//(Ehigh*Tl-Elow*Th)/(Ehigh-Elow);//(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin)
        if ((fabs(*T-Tnew)<0.001)||(fabs(ThetaT-*Theta)<0.0001)) break;
        else *T=Tnew;
        FF_TwoPhasesPreFlashPT(mix,T,P,f,x,y,substPhiL,substPhiG,&ThetaT);
        if (strcmp(opt,"mass")==0){//converts ThetaT to mass fraction if needed
            lMW=0;
            gMW=0;
            for(int i=0;i<mix->numSubs;i++){
                lMW=lMW+x[i]*mix->baseProp[i].MW;
                gMW=gMW+y[i]*mix->baseProp[i].MW;
            }
            ThetaT=ThetaT*gMW/(ThetaT*gMW+(1- ThetaT)*lMW);
        }

        counter++;
    } while(counter<30);

}

//Mixture 2 phases flash, given P,H, composition, and thermo model to use.
void CALLCONV FF_TwoPhasesFlashPH(FF_FeedData *data, FF_PhaseThermoProp phase[2],double *minT, double *maxT,double *refT,double *refP){
    double Tlow,Thigh,Hlow,Hhigh,answerL[3],answerG[3],V,Tplus,Hplus;
    char option,state;
    int i,j;
    if(*minT>0) Tlow=*minT;
    else Tlow=10;
    if(*maxT>0) Thigh=*maxT;
    else Thigh=800;
    option='s';//Ask for the stable phase
    phase[0].T=Tlow;
    FF_MixVfromTPeos(data->mix,&phase[0].T,&data->P,data->z,&option,answerL,answerG,&state);
    if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
        phase[0].V=answerL[0];
    }
    else if((state=='G')||(state=='g')){
        phase[0].V=answerG[0];
    }
    FF_MixThermoEOS(data->mix,refT,refP,&phase[0]);
    Hlow=phase[0].H;
    phase[0].T=Thigh;
    FF_MixVfromTPeos(data->mix,&phase[0].T,&data->P,data->z,&option,answerL,answerG,&state);
    if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
        phase[0].V=answerL[0];
    }
    else if((state=='G')||(state=='g')){
        phase[0].V=answerG[0];
    }
    FF_MixThermoEOS(data->mix,refT,refP,&phase[0]);
    Hhigh=phase[0].H;

    for(j=0;j<10;j++){//Regula falsi loops to bracket the T solution
        phase[0].T=Tlow+(data->H-Hlow)*(Thigh-Tlow)/(Hhigh-Hlow);
        FF_MixVfromTPeos(data->mix,&phase[0].T,&data->P,data->z,&option,answerL,answerG,state);
        if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
            phase[0].V=answerL[0];
        }
        else if((state=='G')||(state=='g')){
            phase[0].V=answerG[0];
        }
        FF_MixThermoEOS(data->mix,refT,refP,&phase[0]);
        if(fabs(phase[0].H-data->H)<0.001){//We have the solution in just one phase. Finish
            return;
        }
        else{
            if(phase[0].H>data->H) Thigh=phase[0].T;
            else Tlow=phase[0].T;
        }
        if((Thigh-Tlow)<0.01) break;//We have arrived to the discontinuity due to change of phase, and the solution will be biphasic
    }

    for(j=0;j<10;j++){
        FF_TwoPhasesFlashPT(data->mix,&phase[0].T,&data->P,data->z,phase[1].c,phase[0].c,phase[1].subsPhi,phase[0].subsPhi,&phase[0].fraction);
        FF_MixVfromTPeos(data->mix,&phase[0].T,&data->P,phase[0].c,&option,answerL,answerG,state);
        if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
            phase[0].V=answerL[0];
        }
        else if((state=='G')||(state=='g')){
            phase[0].V=answerG[0];
        }
        FF_MixThermoEOS(data->mix,refT,refP,&phase[0]);
        Tplus=phase[0].T+0.01;

    }



}


//Mixture 3 phases flash, given P,T, composition, and thermo model to use. By simulated annealing global minimization of the reduced Gibbs energy
void CALLCONV FF_ThreePhasesFlashPTSA(FF_FeedData *data, double x[],double y[],double z[],double substPhiA[],double substPhiB[],double substPhiC[],
                                      double *betaA,double *betaB,double *Gr){
    printf("3 phases flash\n");
    int nSubs,nVar,i,j,k,l,m;
    int na;//number of accepted points
    int ns=30;//number of cycles before step adjustement
    int nt=10;//number of iterations before temperature is reduced
    int ne=5;//number of temperature reductions
    float r;//randon value between 0 and 1
    float rm;//randon value between 0 and 1
    nSubs=data->mix->numSubs;
    nVar=nSubs*2;
    double a[nSubs],b[nSubs],c[nSubs],aTotal,bTotal,cTotal,GrBest,aBest[nVar],fractBest[nVar];//Number of moles of each substance(and its total) in each phase
    double fract[nSubs],answerL[3],answerG[3],difA,difB,difC,Va,Vb,Vc,GrLast,vLast;
    double step;//step vector for next composition calculation
    double T0=5/(R*data->T);//Initial temperature
    double T;//working temperature of the optimization
    double accept,stepCorr,p1,p2;//p1 and p2 are the maximum and minimum accept rate desired
    char optionA,optionB,optionC,state;
    GrLast=1e6;
    GrBest=1e6;
    p1=0.3;
    p2=0.2;

    srand(time(NULL));
    T=T0;
    aTotal=bTotal=cTotal=0;
    for(i=0;i<nSubs;i++){//Initialization of all phases, and the step vector
        a[i]=0.5*data->z[i];
        fract[i]=(double)1/(2+i);
        b[i]=(data->z[i]-a[i])*fract[i];
        c[i]=data->z[i]-a[i]-b[i];
        aTotal=aTotal+a[i];
        bTotal=bTotal+b[i];
        cTotal=cTotal+c[i];
        printf("Initial: a[%i]:%f fract[%i]:%f b[%i]:%f c[%i]:%f\n",i,a[i],i,fract[i],i,b[i],i,c[i]);
    }
    step=1.0;
    printf("Initial: aTotal:%f bTotal:%f cTotal:%f\n",aTotal,bTotal,cTotal);
    for(i=0;i<nSubs;i++){
        x[i]=a[i]/aTotal;
        y[i]=b[i]/bTotal;
        z[i]=c[i]/cTotal;
        printf("Initial: x[%i]:%f y[%i]:%f z[%i]:%f\n",i,x[i],i,y[i],i,z[i]);
    }
    for(m=0;m<ne;m++){//Loops decreasing T
        for(l=0;l<nt;l++){//loops at fixed T
            na=0;//number of accepted points
            printf("Step:%f\n",step);
            for(k=0;k<ns;k++){//loops at fixed step vector
                r=(float) 2*rand()/RAND_MAX-1;//random number between -1 and 1
                printf("Random:%f\n",r);
                for(j=0;j<nVar;j++){//loop along composition
                    if(j<nSubs){//If we change the composition of A phase
                        vLast=a[j];//save the actual value
                        //a[j]=r*step[j]+a[j]*(1-step[j]/data->z[j]);
                        if(r>0) a[j]=a[j]+r*step*(data->z[j]-a[j]);//New point to test
                        else a[j]=a[j]+r*step*a[j];
                        b[j]=(data->z[j]-a[j])*fract[j];
                        c[j]=data->z[j]-a[j]-b[j];
                        printf("change A phase a[%i]:%f b[%i]:%f c[%i]:%f total:%f\n",j,a[j],j,b[j],j,c[j],a[j]+b[j]+c[j]);
                    }
                    else{//If we change the distribution between b and c phases
                        vLast=fract[j-nSubs];//save the actual value
                        if(r>0) fract[j-nSubs]=fract[j-nSubs]+r*step*(1-fract[j-nSubs]);
                        else fract[j-nSubs]=fract[j-nSubs]+r*step*fract[j-nSubs];
                        b[j-nSubs]=(data->z[j-nSubs]-a[j-nSubs])*fract[j-nSubs];
                        c[j-nSubs]=data->z[j-nSubs]-a[j-nSubs]-b[j-nSubs];
                        printf("change ratio a[%i]:%f b[%i]:%f c[%i]:%f total:%f\n",j-nSubs,a[j-nSubs],j-nSubs,b[j-nSubs],j-nSubs,c[j-nSubs],a[j-nSubs]+b[j-nSubs]+c[j-nSubs]);
                    }
                    aTotal=bTotal=cTotal=0;
                    for(i=0;i<nSubs;i++){
                        aTotal=aTotal+a[i];//count of the number of moles in each fraction
                        bTotal=bTotal+b[i];
                        cTotal=cTotal+c[i];
                    }
                    //cTotal=1-aTotal-bTotal;//number of moles in fraction 2
                    difA=0;
                    difB=0;
                    difC=0;
                    for(i=0;i<nSubs;i++){//We convert the mole number to mole fraction
                        x[i]=a[i]/aTotal;//Normalization of the liquid phase to test
                        y[i]=b[i]/bTotal;
                        z[i]=c[i]/cTotal;
                        difA=difA+fabs(x[i]-y[i]);
                        difB=difB+fabs(x[i]-z[i]);
                        difC=difC+fabs(y[i]-z[i]);
                    }
                    //printf("aTotal:%f bTotal:%f cTotal:%f TOTAL:%f x[0]:%f x[1]:%f y[0]:%f y[1]:%f z[0]:%f z[1]:%f\n",aTotal,bTotal,cTotal,aTotal+bTotal+cTotal,x[0],x[1],y[0],y[1],z[0],z[1]);
                    optionA='s';
                    FF_MixVfromTPeos(data->mix,&data->T,&data->P,x,&optionA,answerL,answerG,&state);
                    if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                        optionA='l';
                        Va=answerL[0];
                    }
                    else if((state=='G')||(state=='g')){
                        optionA='g';
                        Va=answerG[0];
                    }
                    else *Gr=1e6;
                    //printf("state:%c Va:%f\n",state,Va);
                    optionB='s';
                    FF_MixVfromTPeos(data->mix,&data->T,&data->P,y,&optionB,answerL,answerG,&state);
                    if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                        optionB='l';
                        Vb=answerL[0];
                    }
                    else if((state=='G')||(state=='g')){
                        optionB='g';
                        Vb=answerG[0];
                    }
                    else *Gr=1e6;
                    optionC='s';
                    FF_MixVfromTPeos(data->mix,&data->T,&data->P,z,&optionC,answerL,answerG,&state);
                    if((state=='L')||(state=='l')||(state=='U')||(state=='u')){
                        optionB='l';
                        Vc=answerL[0];
                    }
                    else if((state=='G')||(state=='g')){
                        optionB='g';
                        Vc=answerG[0];
                    }
                    else *Gr=1e6;
                    //printf("state:%c Va:%f Vb:%f Vc:%f\n",state,Va,Vb,Vc);
                    if(((difA>0.01*nVar)&&(difB>0.01*nVar)&&(difC>0.01*nVar))||(((fabs(Va-Vb)/Va)>0.05)&&((fabs(Vb-Vc)/Vb)>0.05))){
                        //printf("Entrando en calculo de phi\n");
                        FF_MixPhiEOS(data->mix,&data->T,&data->P,x,&optionA,substPhiA);
                        FF_MixPhiEOS(data->mix,&data->T,&data->P,y,&optionB,substPhiB);
                        FF_MixPhiEOS(data->mix,&data->T,&data->P,z,&optionC,substPhiC);
                        *Gr=0;//Residual Gibbs energy
                        for(i=0;i<nSubs;i++){
                            *Gr=*Gr+aTotal*x[i]*log(x[i]*substPhiA[i])+bTotal*y[i]*log(y[i]*substPhiB[i])+cTotal*z[i]*log(z[i]*substPhiC[i]);
                        }
                        //printf("Calculo de Gr Afract:%f Bfract:%f x[0]:%f y[0]:%f z[0]:%f Gr:%f\n",aTotal,bTotal,x[0],y[0],z[0],*Gr);
                    }
                    else *Gr=1e6;

                    if(*Gr>=GrLast){
                        rm=(float) rand()/RAND_MAX;//randon number between 0 and 1 for Metropolis criteria
                        if(exp((GrLast-*Gr)/T)>rm){//
                            printf("R Afract:%f Bfract:%f Cfract:%f TOTAL:%f x[0]:%f y[0]:%f z[0]:%f Gr:%f\n",aTotal,bTotal,cTotal,aTotal+bTotal+cTotal,x[0],y[0],z[0],*Gr);
                            GrLast=*Gr;
                            na=na+1;
                        }
                        else{
                            if(j<nSubs){
                                a[j]=vLast;
                                b[j]=(data->z[j]-a[j])*fract[j];
                                c[j]=data->z[j]-a[j]-b[j];
                            }
                            else{
                                fract[j-nSubs]=vLast;
                                b[j-nSubs]=(data->z[j-nSubs]-a[j-nSubs])*fract[j-nSubs];
                                c[j-nSubs]=data->z[j-nSubs]-a[j-nSubs]-b[j-nSubs];
                            }
                            //printf("M Afract:%f Bfract:%f a[0]:%f b[0]:%f Gr:%f\n",aTotal,bTotal,a[0],b[0],*Gr);
                        }
                    }
                    else{
                        printf("state:%c Va:%f Vb:%f Vc:%f\n",state,Va,Vb,Vc);
                        printf("B Afract:%f Bfract:%f Cfract:%f TOTAL:%f x[0]:%f y[0]:%f z[0]:%f Gr:%f\n",aTotal,bTotal,cTotal,aTotal+bTotal+cTotal,x[0],y[0],z[0],*Gr);
                        GrLast=*Gr;
                        if(*Gr<GrBest){
                            GrBest=*Gr;
                            for(i=0;i<nVar;i++){
                                aBest[i]=a[i];
                                fractBest[i]=fract[i];
                            }
                        }
                        na=na+1;
                    }
                }
            }
            accept=(double) na/(ns*nVar);
            printf("accepted:%i tests:%i acceptance:%f\n",na,ns*nVar,accept);
            //printf("----------------------\n");
            if(accept>p1){
                stepCorr=1+2*(accept-p1)/p2;
                //stepCorr=1+(accept-0.6)/0.2;
                step=step*stepCorr;
                //for(i=0;i<nVar;i++) step[i]=step[i]*stepCorr;
                //printf("accept:%f stepCorr mult:%f\n",accept,stepCorr);
            }
            else if(accept<p2){
                stepCorr=1+2*(p2-accept)/p2;
                //stepCorr=1+(0.4-accept)/0.2;
                step=step/stepCorr;
                //for(i=0;i<nVar;i++) step[i]=step[i]/stepCorr;
                //printf("accept:%f stepCorr divis:%f\n",accept,stepCorr);
            }
            if (step>1) step=1;
        }
        if(m==ne-2){
            for(i=0;i<nVar;i++){
                a[i]=aBest[i];
                fract[i]=fractBest[i];
            }
            T=1e-10;
            step=0.1*step;
        }
        else T=T*0.85;//0.85
        //printf("New T:%f\n",T);
        //printf("-------------------\n");
        /*T=T*0.85;*/
        //printf("New T:%f\n",T);
        //printf("-------------------\n");
    }
    *betaA=aTotal;
    *betaB=bTotal;

}


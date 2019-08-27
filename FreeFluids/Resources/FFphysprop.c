/*
 * FFphysprop.c
 *
 *  Created on: 26/12/2015
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

// contains mainly calculations for pure substances physical properties, by non-eos methods
//Single substance, physical properties correlations
//==================================================

#ifndef FFPHYSPROP_C
#define FFPHYSPROP_C

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
//#include "FFphysprop.h"
#include "FFeosPure.c"
//Calculates the result of the given equation
void CALLCONV FF_CorrelationResult(const int *eq,const double coef[],double x,double *y){//x contains the in variable, and y the out variable
    int i,j;
    double Tm;
    switch (*eq)
    {
    case FF_DIPPR100://DIPPR-100. Polynomial a+b*T+c*T^2+d*T^3+e*T^4. Liquid density,heat capacity,thermal conductivity. Vapor viscosity. Ideal gas heat capacity.Dielectric constant
    case FF_DIPPR100Ld:
        //*y=coef[0]+coef[1]*x+coef[2]*pow(x,2)+coef[3]*pow(x,3)+coef[4]*pow(x,4);
		*y=coef[0]+x*(coef[1]+x*(coef[2]+x*(coef[3]+x*coef[4])));
        break;
    case FF_Polynomial://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5.Ideal gas heat capacity
        *y=coef[0]+coef[1]*x+coef[2]*pow(x,2)+coef[3]*pow(x,3)+coef[4]*pow(x,4)+coef[5]*pow(x,2);
        //printf("%f\n",y[0]);
        break;
    case FF_Polynomial2://Polynomial a+b*T^0.25+c*T^0.5+d*T+e*T^2+f*T^3.
        *y=coef[0]+coef[1]*pow(x,0.125)+coef[2]*pow(x,0.25)+coef[3]*pow(x,0.5)+coef[4]*x;
        //printf("%f\n",y[0]);
        break;
	case FF_expDIPPR100://exp(a + b*T + c*T^2 + d*T^3 + e*T^4). Used for liquid viscosity
                *y=exp(coef[0]+coef[1]*x+coef[2]*pow(x,2)+coef[3]*pow(x,3)+coef[4]*pow(x,4));
        break;
    case FF_DIPPR101://DIPPR-101. exp{a+b/T+c*ln(T)+d*T^e} Liquid viscosity. Vapor pressure.
    case FF_DIPPR101Lv:
    case FF_DIPPR101Vp:
        *y=exp(coef[0]+coef[1]/x+coef[2]*log(x)+coef[3]*pow(x,coef[4]));
          //printf("*y:%f\n",*y);
        break;
    case FF_logDIPPR101:
        *y=coef[0]+coef[1]/x+coef[2]*log(x)+coef[3]*pow(x,coef[4]);
        break;
    case FF_DIPPR102://DIPPR-102. a*T^b/(1+c/T+d/T^2). Vapor viscosity. solid heat capacity
        *y=coef[0]*pow(x,coef[1])/(1+coef[2]/x+coef[3]/pow(x,2));
        break;
    case FF_DIPPR103://DIPPR-103. a + b*exp(-c/T^d)
        *y=coef[0]+coef[1]*exp(-coef[2]/pow(x,coef[3]));
        break;
    case FF_DIPPR104://DIPPR-104 Polynomial. a+b/T+c/T^3+d/T^8+e/T^9. Second virial coefficient
        *y=coef[0]+coef[1]/x+coef[2]/pow(x,3)+coef[3]/pow(x,8)+coef[4]/pow(x,9);
        break;
    case FF_DIPPR105://DIPPR-105. a/b^{1+(1-T/c)^d}. Liquid density
        *y=coef[0]/pow(coef[1],1+pow(1-x/coef[2],coef[3]));
        break;
    case FF_DIPPR106://DIPPR-106. a*(1-Tr)^(b+c*Tr+d*Tr^2+e*Tr^3). Surface tension. Heat of vaporization. Liquid density. Independent variable is always Tr
    case FF_DIPPR106Hv:
    case FF_DIPPR106Ld:
    case FF_DIPPR106SurfT:
        Tm=x/coef[5];
        *y=coef[0]*pow((1-Tm),(coef[1]+coef[2]*Tm+coef[3]*pow(Tm,2)+coef[4]*pow(Tm,3)));
        break;
    case FF_DIPPR107://DIPPR-107. Aly & Lee. a+b*{(c/T)/sinh(c/T)}^2+d*{(e/T)/cosh(e/T)}^2. Ideal gas heat capacity
    case FF_DIPPR107Cp:
        *y=coef[0]+coef[1]*pow(coef[2]/x/sinh(coef[2]/x),2)+coef[3]*pow(coef[4]/x/cosh(coef[4]/x),2);
        break;
    case FF_DIPPR115://DIPPR115. exp{a+b/T+c*ln(T)+d*T^2+e/T^2}. Vapor pressure
        *y=exp(coef[0]+coef[1]/x+coef[2]*log(x)+coef[3]*pow(x,2)+coef[4]*pow(x,-2));
        break;
    case FF_DIPPR116://DIPPR116. a+b*Tr^0.35+c*Tr^(2/3)+ d*Tr + e*Tr^(4/3) with a=Rho crit and Tr=T/f with f=Tc. Liquid density
    case FF_DIPPR116Ld:
    case FF_PPDS10:
        Tm=1-x/coef[5];
        *y=coef[0]+coef[1]*pow(Tm,0.35)+coef[2]*pow(Tm,0.666667)+coef[3]*Tm+coef[4]*pow(Tm,1.333333);
        break;
    case FF_PPDS9://e*exp[a*((c-T)/(T-d))^1/3+b*((c-T)/(T-d))^4/3]
        *y=coef[4]*exp(coef[0]*pow((coef[2]-x)/(x-coef[3]),0.33333)+coef[1]*pow((coef[2]-x)/(x-coef[3]),1.33333));
        break;
    case FF_PPDS15://a/Tm+b+c*Tm+d*Tm^2+e*Tm^3
            Tm=1-x/coef[5];
            *y=coef[0]/Tm+coef[1]+coef[2]*Tm+coef[3]*Tm*Tm+coef[4]*Tm*Tm*Tm;
        break;
    case FF_Wilhoit://Wilhoit equation Cp0 J/mol·K (8 coefficients)
        //if (x>coef[7]) Tm=(x-coef[7])/(x+coef[6]);
        //else Tm=0;
        Tm=(x-coef[7])/(x+coef[6]);
        *y=R*(coef[0] + coef[1]/pow(x,2)*exp(-coef[2]/x) + coef[3]*pow(Tm,2) + (coef[4] - coef[5] /pow((x -coef[7]),2))*pow(Tm,8));
        break;
    case FF_Cooper://Cooper (11 coefficients used in IAPWS95 and CO2),with potential term is used in short fundamental equations with 11 coefficients also(lacks last exp terms)
        *y=coef[0]+coef[1]*pow(x,coef[2]);
        for (j=3;j<13;j=j+2){
            if (coef[j]>0) *y=*y+coef[j]*pow((coef[j+1]/x),2)*exp(coef[j+1]/x)/pow((exp(coef[j+1]/x)-1),2);
        }
        *y=*y*R;
        break;
    case FF_Jaechske://Jaeschke and Schley equation (9 coefficients). Used by GERG2004
        *y=1+coef[0];
        for (j=1;j<9;j=j+4) if (coef[j]>0) *y=*y+coef[j]*pow(coef[j+1]/x/sinh(coef[j+1]/x),2)+coef[j+2]*pow(coef[j+3]/x/cosh(coef[j+3]/x),2);
        *y=*y*R;
        break;
    case FF_ChemSep16://ChemSep nº16. a+exp(b/T+c+d*T+e*T^2). Ideal gas heat capacity. Liquid heat capacity,thermal conductivity,surface tension
        *y=coef[0]+exp(coef[1]/x+coef[2]+coef[3]*x+coef[4]*pow(x,2));
        break;
    case FF_Antoine1://Antoine. 10^{a-b/(T+c)}. Vapor pressure
        *y=pow(10,(coef[0]-coef[1]/(x+coef[2])));
        break;
    case FF_Antoine2://Antoine. exp{a-b/(T+c)}. Vapor pressure
        *y=exp(coef[0]-coef[1]/(x+coef[2]));
        break;
    case FF_Wagner25://Modified Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^2.5+e*Tm^5)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        Tm=1-x/coef[5];
        *y=coef[0]*exp((coef[1]*Tm+coef[2]*pow(Tm,1.5)+coef[3]*pow(Tm,2.5)+coef[4]*pow(Tm,5))/(1-Tm));
    case FF_Wagner36://Original Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^3+e*Tm^6)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        Tm=1-x/coef[5];
        *y=coef[0]*exp((coef[1]*Tm+coef[2]*pow(Tm,1.5)+coef[3]*pow(Tm,3)+coef[4]*pow(Tm,6))/(1-Tm));
        break;
    case FF_PCWIN://PCWIN. a+b*(1-T/e)+c*ln(1-T/e)+d*(1-T/e)^3 with e=Tc, liquid density
        *y=coef[0]+coef[1]*(1-x/coef[4])+coef[2]*log(1-x/coef[4])+coef[3]*pow((1-x/coef[4]),3);
        break;
    case FF_Rackett://Rackett. a/b^{(1-T/c)^d}. Liquid density
        *y=coef[0]/pow(coef[1],pow(1-x/coef[2],coef[3]));
        break;
    case FF_ExtAndrade1://Extended Andrade 1. exp(a+b/T+c*T+d*T^2). Surface tension. Liquid viscosity
        *y=exp(coef[0]+coef[1]/x+coef[2]*x+coef[3]*pow(x,2)+coef[4]*pow(x,3));
        break;
    case FF_ExtAndrade2://Extended Andrade 2. 10^(a+b/T+c*T+d*T^2). Surface tension
        *y=pow(10,(coef[0]+coef[1]/x+coef[2]*x+coef[3]*pow(x,2)+coef[4]*pow(x,3)));
        break;
    case FF_WagnerGd://Original Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^3+e*Tm^6)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        Tm=1-x/coef[5];
        *y=coef[0]*exp(coef[1]*pow(Tm,0.4)+coef[2]*Tm+coef[3]*pow(Tm,2.1)+coef[4]*pow(Tm,5.6));
        break;
    case FF_Tait://Tait equation with P=1e5, V=(a+b*T+c*T^2)*(1-0.0894*ln(1+P/(d*exp(-e*T))). Liquid volume
        *y=(coef[0]+coef[1]*x+coef[2]*pow(x,2))*(1-0.0894*log(1+1e5/(coef[3]*exp(-coef[4]*x))));
        break;
    case FF_ExtWagner://Wagner equation with 6 terms
        Tm=1-x/coef[1];
        *y=coef[0]*exp(coef[1]*(coef[2]*pow(Tm,coef[3])+coef[4]*pow(Tm,coef[5])+coef[6]*pow(Tm,coef[7])+coef[8]*pow(Tm,coef[9])+
                coef[10]*pow(Tm,coef[11])+coef[12]*pow(Tm,coef[13]))/x);
        break;
    }
}

//Calculates physical property, with input and output in SI units(kgr, not moles), using the given correlation, that may or not be in SI units
EXP_IMP void CALLCONV FF_PhysPropCorr(const int cor,const double coef[],const double MW,double x,double *y){
    //printf("%i %f %f %f %f %f %f %f\n",*cor,MW,coef[0],coef[1],coef[2],coef[3],coef[4],coef[5]);
    int i;
    int eq;
    //First we convert the input variable if necessary
    if ((cor==22)||(cor==63)||(cor==240)) x=x-273.15;
    //Second we chose the calculation equation to use
    switch (cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0 in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
    case 40://DIPPR 100 Liquid density in kgr/m3
    case 49://DIPPR 100 Liquid density in cm3/mol
    case 50://DIPPR 100 Liquid thermal conductivity W/(m·K)
    case 60://DIPPR 100 Liquid surface tension N/m
    case 63://DIPPR 100 Liquid surface tension dyna/cm
    case 70://DIPPR 100 Solid density in Kmol/m3
    case 71://DIPPR 100 Solid density in kgr/m3
    case 80://DIPPR 100 Solid Cp in J/Kmol·K
    case 82://DIPPR 100 Solid Cp in J/(kgr·K)
    case 100://DIPPR 100 Vapor density in kgr/m3
    case 111://DIPPR100 Gas viscosity in Pa·s
    case 121://DIPPR100 Gas thermal conductivity in W/(m·K)
    case 140://DIPPR100 T in K from liquid enthalpy in J/kgr
	case 150://DIPPR100 Isothermal compressibility adimensional
        eq=FF_DIPPR100;
        break;
    case 10://Polynomial Cp0 en cal/(mol·K)
        eq=FF_Polynomial;
        //printf("Equation:%i\n",eq);
        break;
    case 130://Polynomial2 boil temp in K
        eq=FF_Polynomial2;
        break;
    case 20://DIPPR101 Vp Pa
    case 21://DIPPR101 Vp KPa
    case 30://DIPPR101 Liquid viscosity Pa·s
        eq=FF_DIPPR101;
        break;
    case 81://DIPPR 102 Solid Cp in J/Kmol·K
    case 110://DIPPR 102 Gas viscosity in Pa·s
    case 120://DIPPR 102 Gas thermal conductivity in W/(m·K)
        eq=FF_DIPPR102;
        break;
    case 41://DIPPR105 Liquid density in mol/dm3
    case 42://DIPPR105 Liquid density in kgr/m3
        eq=FF_DIPPR105;
        break;
    case 45://DIPPR106 Liquid density in mol/dm3
    case 47://DIPPR106 Liquid density in kg/m3
    case 61://DIPPR 106 Liquid surface tension N/m
    case 90://DIPPR 106 HvSat J/Kmol
    case 91://DIPPR 106 HvSat J/kgr
        eq=FF_DIPPR106;
        break;
    case 3://DIPPR107 Cp in cal/mol·K
    case 4://DIPPR107 Cp in J/Kmol·K
	case 200://DIPPR107 Cp in J/kg·K
        eq=FF_DIPPR107;
        break;
    case 5://Wilhoit Cp0 J/mol·K
        eq=FF_Wilhoit;
        break;
    case 7://Cooper Cp0 J/mol·K
        eq=FF_Cooper;
        break;
    case 8://Jaechske Cp0 J/mol·K
        eq=FF_Jaechske;
        break;
    case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
    case 51://ChemSep nº16 Liquid thermal conductivity W/(m·K)
    case 62://chemSep nº16 Liquid surface tension N/m
        eq=FF_ChemSep16;
        break;
    case 22: //Antoine base 10 in C and mmHg
    case 36: //Antoine base 10 in K and cP
        eq=FF_Antoine1;
        break;
    case 23://Antoine base e in K and Pa
        eq=FF_Antoine2;
        break;
    case 24://Wagner 25 in K and Pa
        eq=FF_Wagner25;
        break;
    case 25://Wagner 36 in K and Pa
        eq=FF_Wagner36;
        break;
    case 26://Wagner 36 in K and Pa
        eq=FF_ExtWagner;
        break;
    case 31://Extended Andrade 1 in cP
    case 33://Extended Andrade 1 in cP/MW
        eq=FF_ExtAndrade1;
        break;
    case 32://Extended Andrade 2 in cP
        eq=FF_ExtAndrade2;
        break;
    case 34://Cheric viscosity in cP
        eq=FF_ChericVisc;
        break;
    case 37://PPDS9 in Pa·s
        eq=FF_PPDS9;
        break;
    case 19://PPDS15 Cp liquid in J/kg·K
    case 151://PPDS15 liquid isothermal compressibility adimensional
        eq=FF_PPDS15;
        break;
    case 43://PCWIN liquid density in kgr/m3
        eq=FF_PCWIN;
        break;
    case 48://Rackett liquid density in kgr/dm3
        eq=FF_Rackett;
        break;
    case 44://DIPPR116 liquid density in mol/dm3
    case 46://DIPPR116 liquid density in kg/m3
        eq=FF_DIPPR116;
        break;
    case 101:
        eq=FF_WagnerGd;
        break;
    case 240://Tait equation for polymer density in m3/kg at 1 bar
        eq=FF_Tait;
        break;
    }

    FF_CorrelationResult(&eq,coef,x,y);//We make the calculations

    switch (cor){//Now we convert to SI units
    case 1://DIPPR 100 Cp in KJ/kgr·K
    case 21://DIPPR101 Vp KPa
    case 48://Rackett liquid density in kgr/dm3
        *y=*y*1e3;
        break;
    case 2://DIPPR 100 Cp in J/mol·K
    case 5://Wilhoit Cp0 J/mol·K
    case 7://Cooper Cp0 J/mol·K
    case 8://Jaechske Cp0 J/mol·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
        *y=*y/ MW*1e3;
        break;
    case 3://DIPPR107 Cp in cal/mol·K
    case 10://Polynomial Cp0 en cal/(mol·K)
        *y=*y/ MW*1e3*4.1868;
        break;
    case 4://DIPPR107 Cp in J/Kmol·K
    case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
    case 80://DIPPR 100 Solid Cp in J/Kmol·K
    case 81://DIPPR 102 Solid Cp in J/Kmol·K
    case 90://DIPPR 106 HvSat J/Kmol
        *y=*y/ MW;
        break;
    case 41://DIPPR105 Liquid density in mol/dm3
    case 44://DIPPR116 liquid density in mol/dm3
    case 45://DIPPR106 liquid density in mol/dm3
    case 70://DIPPR 100 Solid density in Kmol/m3
        *y=*y* MW;
        break;
    case 22://Antoine base 10 in C and mmHg
        *y=*y*133.32239;
        break;
    case 31://Extended Andrade 1 in cP
    case 32://Extended Andrade 2 in cP
    case 34://Cheric viscosity in cP
    case 36://Antoine1 viscosity in cP
    case 63://DIPPR 100 Liquid surface tension dyna/cm
        *y=*y*1e-3;
        break;
    case 33://Extended Andrade 1 in cP/MW
        *y=*y* MW*1e-3;
        break;
    case 49://DIPPR 100 Liquid density in cm3/mol
        *y=MW*1e3/ *y;
        break;
    case 240://Tait equation for polymer density in m3/kg at 1 bar
        *y=1/ *y;
    }
    //Last we reconvert the input variable if necessary
    if ((cor==22)||(cor==63)||(cor==240)) x=x+273.15;
}

//Calculates specific enthalpy from a Cp correlation, with reference T=0 K
void FF_SpecificEnthalpyCorr(int cor,const double coef[],double MW, double x, double *H){
    switch (cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
        *H=x*(coef[0]+x*(coef[1]/2+x*(coef[2]/3+x*(coef[3]/4+x*coef[4]/5))));//This is the integration from Cp0
        break;
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        *H=x*(coef[0]+x*(coef[1]/2+x*(coef[2]/3+x*(coef[3]/4+x*(coef[4]/5+x*coef[5]/6)))));//This is the integration from Cp0
        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 4://DIPPR 107 correlation in J/Kmol·K
	case 200://DIPPR107 Cp in J/kg·K
        *H=coef[0]*x+coef[1]*coef[2]*(1/tanh(coef[2]/x))-coef[3]*coef[4]*tanh(coef[4]/x);
        break;
    /*
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        int j=20;
        double interval;
        double T,Cp1,Cp2;
        for (i=0;i<*nPoints;i++){
            interval=x/i;
            Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
            for (j=1;j<21;j++){
                T=*refT+j*interval;
                Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                th0->H=th0->H+(Cp1+Cp2)/2*interval;
                th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                Cp1=Cp2;
                }
        }
        break;
    */
    case 5:{//Wilhoit equation J/mol·K (8 coefficients)
        int j;
        double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
        a7_a6=coef[7]/coef[6];
        a7_a6_2=a7_a6*a7_a6;
        a7_a6_4=a7_a6_2*a7_a6_2;
        x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);
        if (x<=coef[7]) y=0;
        else y=(x-coef[7])/(x+coef[6]);
        y2=y*y;
        y4=y2*y2;
        if (x<=coef[7]) h=0;
        else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
        *H= R*x*(coef[0]+coef[1]*exp(-coef[2]/x)/(coef[2]*x))+R*h;
    }
        break;
    case 7:{//Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
        int j;
        /*
        th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
        for (i=3;i<13;i=i+2){
            if (coef>0) th0->Cp=th0->Cp+coef*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
        }
        th0->Cp=th0->Cp*R;*/
        *H=coef[0]*x+coef[1]/(coef[2]+1)*pow(x,(coef[2]+1));
        for (j=3;j<13;j=j+2) if (coef[j]>0) *H =*H+coef[j]*coef[j+1] / (exp(coef[j+1]/x)-1);
        *H=*H * R;
    }
        break;
    case 8:{//Jaeschke and Schley equation (9 coefficients). Used by GERG2004
        int j;
        *H=(1+coef[0])*x;
        for (j=1;j<9;j=j+4) if (coef[j]>0) *H = *H+2*coef[j]*coef[j+1] / (exp(2*coef[j+1]/x)-1)+
                +2*coef[j+2]*coef[j+3] / (exp(2*coef[j+3]/x)+1);
        *H = *H*R;
    }
        break;
    case 19:{//a/Tm+b+c*Tm+d*Tm^2+e*Tm^3 PPDS15 Cpl in J/(kgr·K)
        double Tm,Tm2;
        Tm=1-x/coef[5];
        Tm2=Tm*Tm;
        *H=-coef[5]*(coef[0]*log(Tm)+coef[1]*Tm+0.5*coef[2]*Tm2+coef[3]*Tm2*Tm/3+0.25*coef[4]*Tm2*Tm2);
      }
        break;
    }
    //now is necessary to pass to J/kg·K
    switch (cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
        *H=*H*1e3;
        break;
    case 2://DIPPR 100 Cp0in J/mol·K
    case 5://Wilhoit equation J/mol·K (8 coefficients)
    case 7://Cooper J/mol·K (11 coefficients)
    case 8://Jaeschke and Schley equation J/mol·K (9 coefficients)
    case 15://DIPPR 100 Liquid Cp in J/mol·K
        *H=*H*1e3/ MW;
        break;
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        *H=*H/ MW;
        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        *H=*H*4.1868*1e3/ MW;
        break;
    }
}

//Calculates specific entropy from a Cp correlation, with reference T=0 K
void FF_SpecificEntropyCorr(int cor,const double coef[],double MW, double x, double *S){
    switch (cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
        *S=coef[0]*log(x)+x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*coef[4]/4)));//This is the integration from Cp0
        break;
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        *S=coef[0]*log(x)+x*(coef[1]+x*(coef[2]/2+x*(coef[3]/3+x*(coef[4]/4+x*coef[5]/5))));
        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 4://DIPPR 107 correlation in J/Kmol·K
	case 200://DIPPR107 Cp in J/kg·K
        *S=coef[0]*log(x)+coef[1]*(coef[2]/x/tanh(coef[2]/x)-log(sinh(coef[2]/x)))-coef[3]*(coef[4]/x*tanh(coef[4]/x)-log(cosh(coef[4]/x)));
        break;
    /*
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        int j=20;
        double interval;
        double T,Cp1,Cp2;
        for (i=0;i<*nPoints;i++){
            interval=x/i;
            Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
            for (j=1;j<21;j++){
                T=*refT+j*interval;
                Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                th0->H=th0->H+(Cp1+Cp2)/2*interval;
                th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                Cp1=Cp2;
                }
        }
        break;
    */
    case 5:{//Wilhoit equation J/mol·K (8 coefficients)
        int j;
        double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
        a7_a6=coef[7]/coef[6];
        a7_a6_2=a7_a6*a7_a6;
        a7_a6_4=a7_a6_2*a7_a6_2;
        x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);
        if (x<=coef[7]) y=0;
        else y=(x-coef[7])/(x+coef[6]);
        y2=y*y;
        y4=y2*y2;
        if (x<=coef[7]) s=0;
        else{
            z = x*(coef[7] + coef[6])/((x + coef[6])*coef[7]);
            w=0;
            for (j=1;j<8;j++) w=w+(x1*pow(-a7_a6,6-j) - coef[4])*pow(y,j)/j;
            s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                (coef[3] + coef[4])*log((x + coef[6])/(coef[6] + coef[7]))-
                (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
            }
        *S=R*(coef[0]*log(x)+coef[1]*(1+coef[2]/x)*exp(-coef[2]/x)/(coef[2]*coef[2])+s);
    }
        break;
    case 7:{//Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
        int j;
        /*
        th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
        for (i=3;i<13;i=i+2){
            if (coef>0) th0->Cp=th0->Cp+coef*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
        }
        th0->Cp=th0->Cp*R;*/
        if (coef[1]>0) *S=coef[0]*log(x)+coef[1]/coef[2]*pow(x,coef[2]);
        else *S=coef[0]*log(x);
        //printf("T:%f S0(0):%f\n",x,S*R/th0->MW);
        for (j=3;j<13;j=j+2){
            if (coef[j]>0) *S = *S+coef[j]*coef[j+1]*(exp(coef[j+1]/x)/(x*(exp(coef[j+1]/x)-1))-
                log(exp(coef[j+1]/x)-1)/coef[j+1]);
            //printf("S0(%i):%f\n",i,S*R/th0->MW);
        }
        *S=*S * R;
    }
        break;
    case 8:{//Jaeschke and Schley equation (9 coefficients). Used by GERG2004
        int j;
        *S=(1+coef[0])*log(x);
        for (j=1;j<9;j=j+4) if (coef[j]>0) *S = *S+coef[j]*(coef[j+1]/x/tanh(coef[j+1]/x)-log(sinh(coef[j+1]/x)))-coef[j+2]*(coef[j+3]/x*tanh(coef[j+3]/x)-log(cosh(coef[j+3]/x)));
        *S = *S * R;
    }
        break;
    case 19:{//PPDS15 equation a/Tm+b+c*Tm+d*Tm^2+e*Tm^3 Cpl in J/(kgr·K)
        double Tc2;
        Tc2=coef[5]*coef[5];
        *S=-(6*Tc2*x*(coef[2]+2*coef[3]+3*coef[4])-3*coef[5]*x*x*(3*coef[4]+coef[3])+2*x*x*x*coef[4])/(6*Tc2*coef[5])+
                (coef[0]+coef[1]+coef[2]+coef[3]+coef[4])*log(x)-coef[0]*log(abs(x-coef[5]));
      }
        break;
		
    }
    //now is necessary to pass to J/kg·K
    switch (cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
        *S=*S*1e3;
        break;
    case 2://DIPPR 100 Cp0in J/mol·K
    case 5://Wilhoit equation J/mol·K (8 coefficients)
    case 7://Cooper J/mol·K (11 coefficients)
    case 8://Jaeschke and Schley equation J/mol·K (9 coefficients)
    case 15://DIPPR 100 Liquid Cp in J/mol·K
        *S=*S*1e3/ MW;
        break;
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        *S=*S/ MW;
        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        *S=*S*4.1868*1e3/ MW;
        break;
    }
}

//Vapor pressure related calculations
//-----------------------------------
//Acentric factor calculation by the definition equation
void CALLCONV FF_WfromDefinition(const double *Pc,const double *Vp,double *w){
    *w=-log(*Vp/ *Pc)/log(10)-1;
}

//Acentric factor calculation from one vapor pressure
void CALLCONV FF_WfromOneVap(const double *Tc,const double *Pc,const double *T,const double *Vp,double *w){
    double Tr,tau,f0,f1,f2;
    Tr=*T/ *Tc;
    tau=1-Tr;
    f0=(-5.97616*tau+1.29874*pow(tau,1.5)-0.60394*pow(tau,2.5)-1.06841*pow(tau,5))/Tr;
    f1=(-5.03365*tau+1.11505*pow(tau,1.5)-5.41217*pow(tau,2.5)-7.46628*pow(tau,5))/Tr;
    *w=(f0-log(*Pc/ *Vp))/f1;
}

//Vapor pressure using Ambrose-Walton equation. Needs only Tc,Pc and w
void CALLCONV FF_VpAmbroseWalton(const FF_BaseProp *baseProp,const double *T,double *Vp){
    int i;
    double Tr,tau,f0,f1,f2;
    *Vp=0;
    if((baseProp->Tc>0)&&(baseProp->Pc>0)&&(baseProp->w>0)){
        if((*T>0)&&(*T<baseProp->Tc)){
            Tr=*T/ baseProp->Tc;
            tau=1-Tr;
            f0=(-5.97616*tau+1.29874*pow(tau,1.5)-0.60394*pow(tau,2.5)-1.06841*pow(tau,5))/Tr;
            f1=(-5.03365*tau+1.11505*pow(tau,1.5)-5.41217*pow(tau,2.5)-7.46628*pow(tau,5))/Tr;
            f2=(-0.64771*tau+2.41539*pow(tau,1.5)-4.26979*pow(tau,2.5)+3.25259*pow(tau,5))/Tr;
            *Vp=baseProp->Pc*exp(f0+f1* baseProp->w+f2*pow(baseProp->w,2));
        }
    }
}

//Vapor pressure using Riedel-Vetere equation. Needs only Tc,Pc and one boiling point
void CALLCONV FF_VpRiedelVetere(const  FF_BaseProp *baseProp,const double *Tref,const double *VpRef,const double *T,double *Vp){
    double Tr,TrefR,PrefR,h,K,psi,alpha,Q;
    int i;
    *Vp=0;
    if((baseProp->Tc>0)&&(baseProp->Pc>0)&&(*Tref>0)&&(*VpRef>0)){
        TrefR=*Tref/ baseProp->Tc;
        PrefR=*VpRef/ baseProp->Pc;
        //first we proceed to calculate h and K
        //first we proceed to calculate h and K
        h=-TrefR*log(PrefR)/(1-TrefR);
        if(baseProp->type==FF_Alcohol)K=0.373-0.03*h;
        else if(baseProp->type==FF_Acid)K=-0.12+0.025*h;
        else K=0.0838;
        //Now psi,alpha
        psi=-35+36/TrefR+42*log(TrefR)-pow(TrefR,6);
        alpha=(3.758*K*psi-log(PrefR)) /(K*psi-log(TrefR));
        Q=K*(3.758-alpha);
        //printf("K,psi,alpha,Q: %f %f %f %f\n",K,psi,alpha,Q);
        //And now Vp
        if((*T>0)&&(*T<baseProp->Tc)){
            Tr=*T/ baseProp->Tc;
            *Vp=baseProp->Pc*exp(-35*Q+36*Q/Tr+(42*Q+alpha)*log(Tr)-Q*pow(Tr,6));
        }
    }
}

//Vapor pressure calculation by correlations
void CALLCONV FF_Vp(double T,FF_SubstanceData *data,double *Vp){
    *Vp=0;
    if(data->vpCorr.form>0) FF_PhysPropCorr(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,T,Vp);
    else if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0)){
        if((data->vp.x>0)&&(data->vp.y>0)) FF_VpRiedelVetere(&data->baseProp,&data->vp.x,&data->vp.y,&T,Vp);
        else if(data->baseProp.w>0) FF_VpAmbroseWalton(&data->baseProp,&T,Vp);
    }
}

//Density equations
//-----------------
//Rackett equation for saturated liquid density. If you supply ref rho, Vc is not used. It is better to use w than Zra, and Zra than Zc
void CALLCONV FF_LiqDensSatRackett(const  FF_BaseProp *baseProp,const double *Tref,const double *rhoRef,const double *T,double *rho){
    double Vc,Zc;//We define local variables in order to use them in the calculation
    int i;
    if (baseProp->Vc==0) Vc=baseProp->Zc*R*baseProp->Tc/baseProp->Pc;
    else Vc=baseProp->Vc;
    if (baseProp->w>0) Zc=0.29056-0.08775*baseProp->w;
    else if (baseProp->Zra>0) Zc=baseProp->Zra;
    else Zc=baseProp->Zc;
    //printf("%f %f %f %f %f %f\n",MW,*Tc,*Pc,*Zc,*Vc,*rhoRef);
        if ((*Tref>0)&&(*rhoRef>0)) *rho= *rhoRef * pow(Zc,(pow((1- *Tref / baseProp->Tc),0.28571)-pow((1- *T/ baseProp->Tc),0.28571)));//Rackett equation with reference data
        else *rho= baseProp->MW/(Vc* pow(Zc,pow((1- *T/ baseProp->Tc),0.28571)))/1000;//Rackett equation with no reference data in kgr/m3
}

//Chueh-Prausnitz pressure correction for liquid density
void CALLCONV FF_LiqDensChuehPrausnitz(const  FF_BaseProp *baseProp, const double *T,const double *P,const double *Pin,double *rhoIn, double *rho){
    int i;
    double Tr,N;
    Tr=*T/baseProp->Tc;
    N=(1-0.89*baseProp->w)*exp(6.9547-76.2853*Tr+191.306*Tr*Tr-203.5472*pow(Tr,3)+82.763*pow(Tr,4));
    *rho=pow((1+9*baseProp->Zc*N*(*P- *Pin)/baseProp->Pc),0.11111)* *rhoIn;
}

//Tait equation for polymer density, with T and P dependence
void CALLCONV FF_LiqDensTait(const int *eq,const double coef[],const  FF_BaseProp *baseProp,const double *T,const double *P, double *rho){
    int i;
    double Tu;
    if (*eq==240) Tu=*T-273.15;
    else Tu=*T;
    rho[i]=(coef[0]+coef[1]*Tu+coef[2]*pow(Tu,2))*(1-0.0894*log(1+P[i]/(coef[3]*exp(-coef[4]*Tu))));
    if (*eq==240) rho[i]=1/rho[i];//As Tait equation gives volume
}

//Liquid density calculation
void CALLCONV FF_LiqDensTP(double T,double *P,FF_SubstanceData *data,double *liqDens){
    double satLiqDens,Vp;
    *liqDens=0;
    if(data->lDensCorr.form>0){
        if(data->baseProp.type==FF_Polymer)FF_LiqDensTait(&data->lDensCorr.form,data->lDensCorr.coef,&data->baseProp,&T,P,&satLiqDens);
        else FF_PhysPropCorr(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,&satLiqDens);
    }
    else if((data->lDens.y>0)&&(data->baseProp.Tc>0)&&((data->baseProp.Vc>0)||((data->baseProp.Pc>0)&&(data->baseProp.Zc>0)))&&((data->baseProp.Zc>0)||(data->baseProp.w>0)))
         FF_LiqDensSatRackett(&data->baseProp,&data->lDens.x,&data->lDens.y,&T,&satLiqDens);
    if((*P>20e5)&&(satLiqDens>0)&&(data->baseProp.type!=FF_Polymer)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)){
        FF_Vp(T,data,&Vp);
        if((*P>Vp)&&(Vp>0)) FF_LiqDensChuehPrausnitz(&data->baseProp,&T,P,&Vp,&satLiqDens,liqDens);
        else *liqDens=satLiqDens;
    }
    else *liqDens=satLiqDens;
}

//Liquid Cp. Bondi method
void CALLCONV FF_LiqCpBondi(const  FF_SubstanceData *data,const double T,double *Cp){
    double Cp0,Tr;
    if(data->cp0Corr.form>0){
        Tr=T/data->baseProp.Tc;
        FF_PhysPropCorr(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,T,&Cp0);
        *Cp=R*1000*(1.586+0.49/(1-Tr)+data->baseProp.w*(4.2775+6.3*pow((1-Tr),0.33333)/Tr+0.4355/(1-Tr)))/data->baseProp.MW+Cp0;
    }
}


//Transport properties of liquids
//-------------------------------

//Lucas liquid viscosity pressure correction.
void CALLCONV FF_LiqViscPcorLucas(double *T,double *P,double *Pref,FF_BaseProp *data,double *rpVisc,double *visc){
    double Tr,Tr2,Tr3,Tr4,A,C,D,dPr;
    Tr=*T/data->Tc;
    Tr2=Tr*Tr;
    Tr3=Tr2*Tr;
    Tr4=Tr3*Tr;
    dPr=(*P- *Pref)/data->Pc;
    A=0.9991-(4.674e-4/(1.0523*pow(Tr,- 0.03877)-1.0513));
    C=-0.07921 + 2.1616*Tr - 13.404*Tr2 + 44.1706*Tr3 - 84.8291*Tr4 + 96.1209*Tr3*Tr2 - 59.8127*Tr3*Tr3 + 15.6719*Tr4*Tr3;
    D=(0.3257/pow(1.0039 - pow(Tr,2.573),0.2906)) - 0.2086;
    *visc=*rpVisc*(1+D*pow((dPr/2.118),A)/(1+C*data->w*dPr));
    //printf("A:%f C:%f D:%f ratio:%f\n",A,C,D,*visc/ *rpVisc);
}

//Liquid viscosity from T,P
void CALLCONV FF_LiqViscTP(double T,double *P,FF_SubstanceData *data,double *liqVisc){
    double satLiqVisc,Vp;
    *liqVisc=0;
    if(data->lViscCorr.form>0){
        FF_PhysPropCorr(data->lViscCorr.form,data->lViscCorr.coef,data->baseProp.MW,T,&satLiqVisc);
        if((*P>10e5)&&(data->baseProp.Tc>0)&&(data->baseProp.Pc>0)&&(data->baseProp.w>0)){
            FF_Vp(T,data,&Vp);
            if((*P>Vp)&&(Vp>0)) FF_LiqViscPcorLucas(&T,P,&Vp,&data->baseProp,&satLiqVisc,liqVisc);
            else *liqVisc=satLiqVisc;
        }
        else *liqVisc=satLiqVisc;
    }
}


//Thermal conductivity of liquids. Latini method
void CALLCONV FF_LiquidThCondLatini(double *T,FF_BaseProp *data,double *thCond){
    double Tr,A,alpha,beta,gamma;
    Tr=*T/data->Tc;
    switch(data->type){
    case FF_Alkane:
        A=0.00350;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Alkene:
    case FF_Alkyne:
        A=0.0361;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Cycloalkane:
        A=0.031;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Aromatic:
        A=0.0346;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Alcohol:
    case FF_Polyol:
        A=0.00339;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Ether:
        A=0.0385;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Ketone:
        A=0.00383;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Acid:
        A=0.00319;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Ester:
        A=0.0415;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    default:
        A=0.00383;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    }
    *thCond=A*pow(data->Tb,alpha)*pow((1-Tr),0.38)/(pow(data->MW,beta)*pow(data->Tc,gamma)*pow(Tr,(1/6)));
}

//Liquid thermal conductivity from T
void CALLCONV FF_LiqThCondT(double T,FF_SubstanceData *data,double *liqThCond){
    *liqThCond=0;
    if(data->lThCCorr.form>0){
        FF_PhysPropCorr(data->lThCCorr.form,data->lThCCorr.coef,data->baseProp.MW,T,liqThCond);
    }
    else FF_LiquidThCondLatini(&T,&data->baseProp,liqThCond);
}


//SurfaceTension, MacLeod-Sugden method. Very sensible to Parachor value
void CALLCONV FF_SurfTensMcLeod(double T,FF_SubstanceData *data,double *surfTens){
    double lDens,Vp,gDens;
    *surfTens=0;
    if(data->baseProp.Pa>0){
        if(data->lDensCorr.form>0){
            FF_PhysPropCorr(data->lDensCorr.form,data->lDensCorr.coef,data->baseProp.MW,T,&lDens);
        }
        else if((data->lDens.x>0)&&(data->lDens.y>0)) FF_LiqDensSatRackett(&data->baseProp,&data->lDens.x,&data->lDens.y,&T,&lDens);
        else return;

        if(data->gDensCorr.form>0){
            FF_PhysPropCorr(data->gDensCorr.form,data->gDensCorr.coef,data->baseProp.MW,T,&gDens);
        }
        else if(data->vpCorr.form>0){
            FF_PhysPropCorr(data->vpCorr.form,data->vpCorr.coef,data->baseProp.MW,T,&Vp);
        }
        else FF_VpAmbroseWalton(&data->baseProp,&T,&Vp);
        gDens=0.3*data->baseProp.MW*Vp/(R* T)*1e-3;
        *surfTens=pow((data->baseProp.Pa*(lDens-gDens)*1000/data->baseProp.MW),4);
        //printf("lDens:%f gDens:%f\n",lDens,gDens);
    }
}

//SurfaceTension, Sastri-Rao method. Very good approximation is obtained
void CALLCONV FF_SurfTensSastri(double *T,FF_BaseProp *data,double *surfTens){
    //double Tm=1- *T/data->Tc;
    //*surfTens=kb*data->Tc*pow((Av/data->Vc),0.666667)*(4.35+4.14*data->w)*pow(Tm,1.26)*(1+0.19*pow(Tm,0.5)-0.25*Tm);//Not good for alcohols
    //Use Sastri-Rao:
    double k,x,y,z,m;
    switch(data->type){
    case FF_Alcohol:
    case FF_Polyol:
        k=2.28;
        x=0.25;
        y=0.175;
        z=0;
        m=0.8;
        break;
    case FF_Acid:
        k=0.125;
        x=0.5;
        y=-1.5;
        z=1.85;
        m=1.222;
        break;
    default:
        k=0.158;
        x=0.5;
        y=-1.5;
        z=1.85;
        m=1.222;
        break;
    }
    *surfTens=k*pow(data->Pc*1e-5,x)*pow(data->Tb,y)*pow(data->Tc,z)*pow(((1- *T/data->Tc)/(1-data->Tb/data->Tc)),m)*1e-3;
}

//Surface tension from T
void CALLCONV FF_SurfTensT(double T,FF_SubstanceData *data,double *surfTens){
    *surfTens=0;
    if (data->lSurfTCorr.form>0){
        FF_PhysPropCorr(data->lSurfTCorr.form,data->lSurfTCorr.coef,data->baseProp.MW,T,surfTens);
    }
    else if (data->lDensCorr.form>0){
        FF_SurfTensMcLeod(T,data,surfTens);
    }
    else FF_SurfTensSastri(&T,&data->baseProp,surfTens);
}

//Transport properties of gases
//-----------------------------

//Gas viscosity pressure prediction/correction. Chung method
void CALLCONV FF_GasViscTVcpChung(double *T,double *V,FF_BaseProp *data,double *ldVisc,double *visc){
    double omega,Ta,sigma,Fc,muR,muR4,k;
    double E[10],y,G1,G2,eta2,eta1;
    int i;
    //Chung method
    Ta=1.2593* *T/data->Tc;
    sigma= 0.809*pow(data->Vc,0.33333);
    omega= 1.16145*pow(Ta,- 0.14874)+0.52487*exp(- 0.7732*Ta)+2.16178*exp(- 2.43787*Ta);
    if ((data->mu<99)&&(data->mu>0)) muR=131.3*data->mu/pow(data->Vc*1e6*data->Tc,0.5);
    else muR=0;
    muR4=muR*muR*muR*muR;
    if(data->type==FF_Water) k=0.076;//Water
    else if(data->type==FF_Alcohol) k=0.0682+4.74/data->MW;//alcohol
    else if(data->type==FF_Polyol) k=0.0682+4.74*2/data->MW;//Polyol
    else k=0;
    Fc=1-0.2756*data->w+0.059035*muR4+k;
    if(*ldVisc<=0){
        *ldVisc=26.692*Fc*pow((data->MW* *T),0.5)*1e-11/(sigma*sigma*omega);
    }

    double coef[10][4]={{6.324,50.412,-51.680,1189.0},{1.210e-3,-1.154e-3,-6.257e-3,0.03728},{5.283,254.209,-168.48,3898},{6.623,38.096,-8.464,31.42},
                        {19.745,7.63,-14.354,31.53},{-1.9,-12.537,4.985,-18.15},{24.275,3.45,-11.291,69.35},{0.7972,1.117,0.01235,-4.117},
                        {-0.2382,0.0677,-0.8163,4.025},{0.06863,0.3479,0.5926,-0.727}};
    for(i=0;i<10;i++) E[i]=coef[i][0]+coef[i][1]*data->w+coef[i][2]*muR4+coef[i][3]*k;
    y=data->Vc/(6* *V);
    G1=(1-0.5*y)/(pow((1-y),3));
    G2=(E[0]*((1-exp(-E[3]*y))/y)+E[1]*G1*exp(E[4]*y)+E[2]*G1)/(E[0]*E[3]+E[1]+E[2]);
    eta2=E[6]*y*y*G2*exp(E[7]+E[8]/Ta+E[9]/(Ta*Ta));
    eta1=pow(Ta,0.5)*(Fc*(1/G2+E[5]*y))/omega+eta2;
    if(eta1>1) *visc=*ldVisc*eta1;
    else *visc=*ldVisc;
    //printf("Ta:%f omega:%f muR:%f Fc:%f ldVisc:%f y:%f G1:%f G2:%f eta2:%f eta1:%f\n",Ta,omega,muR,Fc,*ldVisc,y,G1,G2,eta2,eta1);
    //for(i=0;i<10;i++) printf("E[%i]:%f\n",i,E[i]);
}

//Gas viscosity pressure prediction/correction. Lucas method
void CALLCONV FF_GasViscTPcpLucas(const double *T,const double *P,const FF_BaseProp *data,double *lpVisc,double *visc){
    double Z1,mur,Q,xi,Tr,Fp0,Fq0,Z2,Pr,a,b,c,d,e,f,Y,Fp,Fq;
    Tr=*T/data->Tc;
    Pr=*P/data->Pc;
    if((data->mu<999)&&(data->mu>0)) mur=52.46*data->mu*data->mu*data->Pc*1e-5/(data->Tc*data->Tc);
    else mur=0;
    if(mur<0.022) Fp0=1;
    else if(mur<0.075) Fp0=1+30.55*pow((0.292-data->Zc),1.72);
    else Fp0=1+30.55*pow((0.292-data->Zc),1.72)*fabs(0.96+0.1*(Tr-0.7));
    if(data->MW<4.04){
        if(data->MW<3) Q=0.76;//Hydrogen
        else if(data->MW<4.02) Q=1.38;//He
        Fq0=1.22*pow(Q,0.15);
    }
    else Fq0=1;
    xi=0.176*pow(data->Tc/(data->MW*data->MW*data->MW*data->Pc*data->Pc*data->Pc*data->Pc*1e-20),0.1666667)*1e-7;
    if(*lpVisc<=0){//If low pressure density is not supplied it is calculated
        *lpVisc=Fp0*Fq0*(0.807*pow(Tr,0.618)-0.357*exp(-0.449*Tr)+0.34*exp(-4.058*Tr)+0.018)*1e-14/xi;
    }
    //printf("lpVisc:%f\n",*lpVisc);
    Z1=*lpVisc*xi;
    if(Tr<=1){//The original method has been modified acording to Chemsep book
        Tr=1;
        //alpha=3.262+14.98*pow(Pr,5.508);
        //beta=1.39+5.746*Pr;
        //Z2=0.6+0.76*pow(Pr,alpha)+(6.99*pow(Pr,beta)-0.6)*(1-Tr);
    }
    a=1.245e-3*exp(5.1726*pow(Tr,-0.3286))/Tr;
    b=a*(1.6553*Tr-1.2723);
    c=0.4489*exp(3.0578*pow(Tr,-37.7332))/Tr;
    d=1.7368*exp(2.231*pow(Tr,-7.6351))/Tr;
    e=1.3088;
    f=0.9425*exp(-0.1853*pow(Tr,0.4489));
    Z2=Z1*(1+(a*pow(Pr,e))/(b*pow(Pr,f)+1/(1+c*pow(Pr,d))));

    Y=Z2/Z1;
    Fp=(1+(Fp0-1)*pow(Y,-3))/Fp0;
    Fq=(1+(Fq0-1)*(1/Y-0.007*pow(log(Y),4)))/Fq0;
    *visc=Z2*Fp*Fq/xi;
    //printf("lpVisc:%f visc:%f\n",*lpVisc,*visc);
}

//Gas viscosity from T,P
void CALLCONV FF_GasViscTP(double T,double *P,FF_SubstanceData *data,double *gasVisc){
    double lpGasVisc=0;
    *gasVisc=0;
    if(data->gViscCorr.form>0){
        FF_PhysPropCorr(data->gViscCorr.form,data->gViscCorr.coef,data->baseProp.MW,T,&lpGasVisc);
    }
    if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0))FF_GasViscTPcpLucas(&T,P,&data->baseProp,&lpGasVisc,gasVisc);
}

//Gas low pressure thermal conductivity prediction
void CALLCONV FF_GasLpThCondTCpChung(double *T,double *Cp0,FF_BaseProp *data,double *ldThCond){
    double ldVisc,visc,alpha,Tr,beta,Z,psi,P;
    //Chung method
    ldVisc=0;
    P=1e5;
    FF_GasViscTPcpLucas(T,&P,data,&ldVisc,&visc);
    alpha=(*Cp0-R)/R-1.5;
    beta=0.7862 - 0.7109*data->w + 1.3168 *data->w*data->w;
    Tr=*T/data->Tc;
    Z=2.0 + 10.5*Tr*Tr;
    psi=1+alpha*((0.215 + 0.28288 *alpha - 1.061 *beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta ));
    *ldThCond=3.75*psi*ldVisc*R*1000/data->MW;
}

//Gas thermal conductivity pressure correction
void CALLCONV FF_GasThCondTVcorChung(double *T,double *V,FF_BaseProp *data,double *ldThCond,double *thCond){
    int i;
    double Tr,y,B[7],G1,muR,muR4,k,G2;
    double coef[7][4]={{2.4166,7.4824e-1,-9.1858e-1,1.2172e2},{-5.0924e-1,-1.5094,-4.9991e1,6.9983e1},{6.6107,5.6207,6.4760e1,2.7039e1},{1.4543e1,-8.9139,-5.6379,7.4344e1},
                        {7.9274e-1,8.2019e-1,-6.9369e-1,6.3173},{-5.8634,1.2801e1,9.5893,6.5529e1},{9.1089e1,1.2811e2,-5.4217e1,5.2381e2}};
    Tr=*T/data->Tc;
    y=data->Vc/(6* *V);
    G1=(1-0.5*y)/(pow((1-y),3));
    if ((data->mu<99)&&(data->mu>0)) muR=131.3*data->mu/pow(data->Vc*1e6*data->Tc,0.5);
    else muR=0;
    muR4=muR*muR*muR*muR;
    if(data->type==FF_Water) k=0.076;//Water
    else if(data->type==FF_Alcohol) k=0.0682+4.74/data->MW;//alcohol
    else if(data->type==FF_Polyol) k=0.0682+4.74*2/data->MW;//Polyol
    else k=0;
    for(i=0;i<7;i++) B[i]=coef[i][0]+coef[i][1]*data->w+coef[i][2]*muR4+coef[i][3]*k;
    G2=(B[0]*(1-exp(-B[3]*y))/y+B[1]*G1*exp(B[4]*y)+B[2]*G1)/(B[0]*B[3]+B[1]+B[2]);
    *thCond=*ldThCond*(1/G2+B[5]*y)+3.586e-3*pow(data->Tc*1000/data->MW,0.5)*B[6]*y*y*pow(Tr,0.5)*G2/pow(data->Vc*1e6,0.6667);
    if(*thCond<*ldThCond) *thCond=*ldThCond;

    //printf("Ta:%f omega:%f muR:%f Fc:%f ldVisc:%f y:%f G1:%f G2:%f eta2:%f eta1:%f\n",Ta,omega,muR,Fc,*ldVisc,y,G1,G2,eta2,eta1);
    //for(i=0;i<10;i++) printf("E[%i]:%f\n",i,E[i]);
}

//Gas thermal conductivity from T,V
void CALLCONV FF_GasThCondTV(double T,double *V,FF_SubstanceData *data,double *gasThCond){
    double Cp0,ldGasThCond=0;
    *gasThCond=0;
    if(data->gThCCorr.form>0){
        FF_PhysPropCorr(data->gThCCorr.form,data->gThCCorr.coef,data->baseProp.MW,T,&ldGasThCond);
    }
    else{
        if(data->cp0Corr.form>0){
            FF_PhysPropCorr(data->cp0Corr.form,data->cp0Corr.coef,data->baseProp.MW,T,&Cp0);
            FF_GasLpThCondTCpChung(&T,&Cp0,&data->baseProp,&ldGasThCond);
        }
        else return;
    }
    FF_GasThCondTVcorChung(&T,V,&data->baseProp,&ldGasThCond,gasThCond);
}
#endif /* FFPHYSPROP_C */

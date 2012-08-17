//
//  derivatives.cpp
//  finance_analyze
//
//  Created by Thorsten Kurth on 17.08.12.
//
//

#include "derivatives.h"

BlackScholes::BlackScholes(double rr, double ssigma, double TT, double EE, double DD) : r(rr), sigma(ssigma), T(TT), E(EE), D(DD) {
    optiontype="none";
    set_modelid();
}

BlackScholes::BlackScholes(double rr, double ssigma, double TT, double EE, std::string ooptiontype, double DD) : r(rr), sigma(ssigma), T(TT), E(EE), D(DD), optiontype(ooptiontype) {
    set_modelid();
}

void BlackScholes::set_modelid(){
    if(optiontype.compare("eurocall")){
        modelid=1;
    }
    else if(optiontype.compare("europut")){
        modelid=2;
    }
    else if(optiontype.compare("bincall")){
        modelid=3;
    }
    else if(optiontype.compare("binput")){
        modelid=4;
    }
    else{
        modelid=-1;
    }
}

//Values of options:
double BlackScholes::Eurocall(double S, int flag){
    double result=0.;
    Normaldist N;

    switch(flag){
        //value:
        case 0:
            if(expired){
                result=max(S-E,0.);
            }
            else{
                result=S*exp(-D*tau)*N.cdf(d1)-E*exp(-r*tau)*N.cdf(d2);
            }
            break;
        //delta:
        case 1:
            if(!expired){
                result=exp(-D*tau)*N.cdf(d1);
            }
            break;
        //gamma:
        case 2:
            if(!expired){
                result=exp(-D*tau)*N.pdf(d1)/(sigma*S*sqrt(tau));
            }
            break;
        //Theta:
        case 3:
            if(!expired){
                result=-sigma*S*exp(-D*tau)*N.pdf(d1)/(2.*sqrt(tau))+D*S*N.cdf(d1)*exp(-D*tau)-r*E*exp(-r*tau)*N.cdf(d2);
            }
            break;
        //Speed:
        case 4:
            if(!expired){
                result=-exp(-D*tau)*N.pdf(d1)/(sigma*sigma*S*S*tau)*(d1+sigma*sqrt(tau));
            }
            break;
        //Vega:
        case 5:
            if(!expired){
                result=S*sqrt(tau)*exp(-D*tau)*N.pdf(d1);
            }
            break;
        //Rho:
        case 6:
            if(!expired){
                result=E*tau*exp(-r*tau)*N.cdf(d2);
            }
            break;
    }
    return result;
}

double BlackScholes::Europut(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(flag){
        //value:
        case 0:
            if(expired){
                result=max(E-S,0.);
            }
            else{
                result=-S*exp(-D*tau)*N.cdf(-d1)+E*exp(-r*tau)*N.cdf(-d2);
            }
            break;
        //delta:
        case 1:
            if(!expired){
                result=exp(-D*tau)*(N.cdf(d1)-1.);
            }
            break;
        //gamma:
        case 2:
            if(!expired){
                result=exp(-D*tau)*N.pdf(d1)/(sigma*S*sqrt(tau));
            }
            break;
        //Theta:
        case 3:
            if(!expired){
                result=sigma*S*exp(-D*tau)*N.pdf(-d1)/(2.*sqrt(tau))-D*S*N.cdf(-d1)*exp(-D*tau)+r*E*exp(-r*tau)*N.cdf(-d2);
            }
            break;
       //speed:
        case 4:
            if(!expired){
                result=-exp(-D*tau)*N.pdf(d1)/(sigma*sigma*S*S*tau)*(d1+sigma*sqrt(tau));
            }
            break;
        //Vega:
        case 5:
            if(!expired){
               result=S*sqrt(tau)*exp(-D*tau)*N.pdf(d1); 
            }
            break;
        //Rho:
        case 6:
            if(!expired){
                result=-E*tau*exp(-r*tau)*N.cdf(-d2);
            }
            break;
    }
    return result;
}

double BlackScholes::Bincall(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(flag){
        //value:
        case 0:
            if(expired){
                result=( S>E ? 1. : 0.);
            }
            else{
                result=exp(-r*tau)*N.cdf(d2);
            }
            break;
        //delta:
        case 1:
            if(!expired){
                result=exp(-r*tau*N.pdf(d2))/(sigma*S*sqrt(tau));
            }
            break;
        //gamma
        case 2:
            if(!expired){
                result=-exp(-r*tau)*d1*N.pdf(d2)/(sigma*sigma*S*S*tau);
            }
            break;
        //Theta:
        case 3:
            if(!expired){
                result=r*exp(-r*tau)*N.cdf(d2)+exp(-r*tau)*N.pdf(d2)*(d1/(2.*tau)-(r-D)/(sigma*sqrt(tau)));
            }
            break;
        //Speed:
        case 4:
            if(!expired){
                result=-exp(-r*tau)*N.pdf(d2)/(sigma*sigma*S*S*S*tau)*(-2.*d1+(1.-d1*d2)/(sigma*sqrt(tau)));
            }
            break;
        //Vega:
        case 5:
            if(!expired){
                result=-exp(-r*tau)*N.pdf(d2)*(sqrt(tau)+d2/sigma);
            }
            break;
        //Rho:
        case 6:
            if(!expired){
                result=-sqrt(tau)/sigma*exp(-r*tau)*N.pdf(d2);
            }
            break;
    }
    return result;
}

double BlackScholes::Binput(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(flag){
        //value:
        case 0:
            if(expired){
                result=( S<E ? 1. : 0.);
            }
            else{
                result=exp(-r*tau)*N.cdf(d2);
            }
            break;
        //delta:
        case 1:
            if(!expired){
                result=-exp(-r*tau)*N.pdf(d2)/(sigma*S*sqrt(tau));
            }
            break;
        //gamma:
        case 2:
            if(!expired){
                result=exp(-r*tau)*d1*N.pdf(d2)/(sigma*sigma*S*S*(tau));
            }
            break;
        //Theta:
        case 3:
            if(!expired){
                result=r*exp(-r*tau)*(1.-N.cdf(d2))-exp(-r*tau)*N.pdf(d2)*(d1/(2.*tau)-(r-D)/(sigma*sqrt(tau)));
            }
            break;
        //Speed:
        case 4:
            if(!expired){
                result=exp(-r*tau)*N.pdf(d2)/(sigma*sigma*S*S*S*tau)*(-2.*d1+(1.-d1*d2)/(sigma*sqrt(tau)));
            }
            break;
        //Vega:
        case 5:
            if(!expired){
                result=exp(-r*tau)*N.pdf(d2)*(sqrt(tau)+d2/sigma);
            }
            break;
        //Rho:
        case 6:
            if(!expired){
                result=sqrt(tau)/sigma*exp(-r*tau)*N.pdf(d2);
            }
            break;

    }
    return result;
}

double BlackScholes::operator()(double S, double t){
    return Value(S,t);
}

//generalized function of options:
double BlackScholes::Calc(double S, double t, int flag){
    double result;
    if( (fabs(t-T)<1.e-4) || t>T){
        tau=0.;
        expired=true;
    }
    else{
        tau=T-t;
        d1=(log(S/E)+(r-D+0.5*sigma*sigma)*tau)/sqrt(sigma*sigma*tau);
        d2=d1-sigma*sqrt(tau);
        expired=false;
    }
    
    switch(modelid){
        case 1:
            result=Eurocall(S,flag);
            break;
        case 2:
            result=Europut(S,flag);
            break;
        case 3:
            result=Bincall(S,flag);
            break;
        case 4:
            result=Binput(S,flag);
            break;
        default:
            result=0.;
            break;
    }
    return result;
}

//value of options:
double BlackScholes::Value(double S, double t){
    return Calc(S,t,0);
}

//Deltas of options:
double BlackScholes::Delta(double S, double t){
    return Calc(S,t,1);
}

//Gammas of options:
double BlackScholes::Gamma(double S, double t){
    return Calc(S,t,2);
}

//Thetas of options:
double BlackScholes::Theta(double S, double t){
    return Calc(S,t,3);
}

//Speed of options:
double BlackScholes::Speed(double S, double t){
    return Calc(S,t,4);
}

//Vega of options:
double BlackScholes::Vega(double S, double t){
    return Calc(S,t,5);
}

//Rho of options:
double BlackScholes::Rho(double S, double t){
    return Calc(S,t,6);
}


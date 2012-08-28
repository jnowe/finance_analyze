//
//  derivatives.cpp
//  finance_analyze
//
//  Created by Thorsten Kurth on 17.08.12.
//
//

#include "derivatives.h"

BlackScholes::BlackScholes(double rr, double ssigma, double TT, double EE, double DD) : r(rr), sigma(ssigma), T(TT), E(EE), D(DD) {
    kerneltype="vanilla";
    payofftype="none";
    set_kernelid();
    set_payoffid();
}

BlackScholes::BlackScholes(double rr, double ssigma, double TT, double EE, std::string ppayofftype, double DD) : r(rr), sigma(ssigma), T(TT), E(EE), D(DD), payofftype(ppayofftype) {
    kerneltype="vanilla";
    set_kernelid();
    set_payoffid();
}

BlackScholes::BlackScholes(double rr, double ssigma, double TT, double EEhigh, double EElow, double DD) : r(rr), sigma(ssigma), T(TT), E(EEhigh), Ehigh(EEhigh), Elow(EElow), D(DD) {
    kerneltype="inline";
    payofftype="none";
    set_kernelid();
    set_payoffid();
}

void BlackScholes::set_kernelid(){
    if(kerneltype.compare("vanilla")){
        kernelid=1;
    }
    else if(payofftype.compare("inline")){
        kernelid=2;
    }
    else{
        kernelid=-1;
    }
}

void BlackScholes::set_payoffid(){
    if(payofftype.compare("eurocall")){
        payoffid=1;
    }
    else if(payofftype.compare("europut")){
        payoffid=2;
    }
    else if(payofftype.compare("bincall")){
        payoffid=3;
    }
    else if(payofftype.compare("binput")){
        payoffid=4;
    }
    else{
        payoffid=-1;
    }
}

//Values of options:
double BlackScholes::Eurocall(double S, int flag){
    double result=0.;
    Normaldist N;

    switch(kernelid){
        //vanilla:
        case 1:
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
                    //sensitivity to dividend yield:
                case 7:
                    if(!expired){
                        result=-tau*S*exp(-D*tau)*N.cdf(d1);
                    }
                    break;
            }
            break;
        //inline:
        case 2:
            break;
    }
    return result;
}

double BlackScholes::Europut(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(kernelid){
        //vanilla:
        case 1:
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
                    //sensitivity to dividend yield:
                case 7:
                    if(!expired){
                        result=tau*S*exp(-D*tau)*N.cdf(-d1);
                    }
                    break;
            }
            break;
        //inline
        case 2:
            break;
    }
    return result;
}

double BlackScholes::Bincall(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(kernelid){
        //vanilla
        case 1:
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
                    //sensitivity to dividend yield:
                case 7:
                    if(!expired){
                        result=-tau/sigma*exp(-r*tau)*N.pdf(d2);
                    }
                    break;
            }
            break;
        //inline:
        case 2:
            break;
    }
    return result;
}

double BlackScholes::Binput(double S, int flag){
    double result=0.;
    Normaldist N;
    
    switch(kernelid){
        //vanilla
        case 1:
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
                    //sensitivity to dividend yield:
                case 7:
                    if(!expired){
                        result=tau/sigma*exp(-r*tau)*N.pdf(d2);
                    }
                    break;
            }
            break;
        //inline
        case 2:
            break;
    }
    return result;
}

//generalized function of options:
double BlackScholes::Calc(double S, double t, int flag){
    double result;
    tau=T-t;
    if( (fabs(tau)<1.e-4) || tau<0.){
        expired=true;
    }
    else{
        d1=(log(S/E)+(r-D+0.5*sigma*sigma)*tau)/sqrt(sigma*sigma*tau);
        d2=d1-sigma*sqrt(tau);
        expired=false;
    }
    
    switch(payoffid){
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

//kernels:
double BlackScholes::Kernel(double S, double Sprime, double t){
    double result=0.;
    double xdiff=log(S/Sprime);
    double xlow, xhigh, x, xprime;
    tau=T-t;
    
    if(fabs(tau)<1.e-4 || tau<0.){
        expired=true;
        result=0.;
    }
    else{
        switch(kernelid){
            //vanilla kernel:
            case 1:
                result=exp(-r*tau)/(sigma*sqrt(2.*pimath*tau))*exp(-sqr(xdiff+tau*(r-D-0.5*sigma*sigma))/(2.*tau*sigma*sigma));
                break;
            //inline kernel:
            case 2:
                alpha=-0.5*(2.*(r-D)/(sigma*sigma)-1.);
                beta=-0.5*(sqr(2.*r/(sigma*sigma)+1.)+sqr(2.*D/(sigma*sigma)+1.)-8.*r*D/(sqr(sigma*sigma))-1.)/sqr(sigma);
                xlow=log(Elow);
                xhigh=log(Ehigh);
                x=log(S);
                xprime=log(Sprime);
                
                if(x<xlow || xprime<xlow || x>xhigh || xprime>xhigh) result=0.;
                else{
                    //first partial sums:
                    double sn=exp(-sqr(x-xprime)/(2.*tau*sigma*sigma))-exp(-sqr(x+xprime-2.*xlow)/(2.*tau*sigma*sigma));
                    int nup=3;
                    for(int n=1; n<=nup; n++){
                        sn+=exp(-sqr(x-xprime+2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))+exp(-sqr(x-xprime-2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))
                            -exp(-sqr(x+xprime-2.*xlow-2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))-exp(-sqr(x+xprime-2.*xlow+2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma));
                    }
                    
                    //Levin u transformation:
                    Levin acc(250,1.e-10);
                    double omegan, an;
                    int nrun=nup+1;
                    while( !acc.cnvgd ){
                        an=exp(-sqr(x-xprime+2.*nrun*(xhigh-xlow))/(2.*tau*sigma*sigma))+exp(-sqr(x-xprime-2.*nrun*(xhigh-xlow))/(2.*tau*sigma*sigma))
                            -exp(-sqr(x+xprime-2.*xlow-2.*nrun*(xhigh-xlow))/(2.*tau*sigma*sigma))-exp(-sqr(x+xprime-2.*xlow+2.*nrun*(xhigh-xlow))/(2.*tau*sigma*sigma));
                        sn+=an;
                        omegan=(1.+(double)nrun)*an;
                        result=acc.next(sn,omegan);
                        nrun++;
                    }
                    //switch back to direct resummation if Levin fails:
                    if(isnan(result)){
                        result=exp(-sqr(x-xprime)/(2.*tau*sigma*sigma))-exp(-sqr(x+xprime-2.*xlow)/(2.*tau*sigma*sigma));
                        for(int n=1; n<=20; n++){
                            result+=exp(-sqr(x-xprime+2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))+exp(-sqr(x-xprime-2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))
                                    -exp(-sqr(x+xprime-2.*xlow-2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma))-exp(-sqr(x+xprime-2.*xlow+2.*n*(xhigh-xlow))/(2.*tau*sigma*sigma));
                        }
                    }
                    result*=exp(beta*tau+alpha*xdiff)/sqrt(2.*pimath*tau*sigma*sigma);
                }
                break;
        }
    }
    return result;
}

//value of options:
double BlackScholes::operator()(double S, double t){
    return Value(S,t);
}

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


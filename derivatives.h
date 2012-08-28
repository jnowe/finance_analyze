#include "mathutils.hpp"

#ifndef DERIVATIVES
#define DERIVATIVES

class BlackScholes{
public:
    BlackScholes(double rr, double ssigma, double TT, double EE, double DD=0.);
    BlackScholes(double rr, double ssigma, double TT, double EE, std::string ppayofftype, double DD=0.);
    BlackScholes(double rr, double ssigma, double TT, double EEhigh, double EElow, double DD=0.);
    double Kernel(double S, double Sprime, double t);
    double operator()(double S, double t);
    double Calc(double S, double t, int flag=0);
    double Value(double S, double t);
    double Delta(double S, double t);
    double Gamma(double S, double t);
    double Theta(double S, double t);
    double Speed(double S, double t);
    double Vega(double S, double t);
    double Rho(double S, double t);
    double ImpliedVol(double V, double S, double t);
    
private:
    bool expired;
    int payoffid, kernelid;
    double r; //risk-free rate of return
    double sigma; //volatility
    double T; //maturity-time
    double D; //continuous dividend yield or foreign currency interest rate
    double E; //strike price
    double Ehigh, Elow; //boundaries for double barrier options
    double d1,d2,tau; //helper values, combination of parameters
    double alpha,beta; //parameters used if BS equation is substituted
    std::string payofftype, kerneltype; //type of option
    void set_payoffid();
    void set_kernelid();
    
    //routines for specific options:
    double Eurocall(double S, int flag=0);
    double Europut(double S, int flag=0);
    double Bincall(double S, int flag=0);
    double Binput(double S, int flag=0);
};

#endif

#include "mathutils.hpp"

#ifndef DERIVATIVES
#define DERIVATIVES

class BlackScholes{
public:
    BlackScholes(double rr, double ssigma, double TT, double EE, double DD=0.);
    BlackScholes(double rr, double ssigma, double TT, double EE, std::string oopiontype, double DD=0.);
    double operator()(double S, double t);
    double BlackScholes::Calc(double S, double t, int flag=0);
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
    int modelid;
    double r; //risk-free rate of return
    double sigma; //volatility
    double T; //maturity-time
    double D; //continuous dividend yield or foreign currency interest rate
    double E; //strike price
    double d1,d2,tau; //helper values, combination of parameters
    std::string optiontype; //type of option
    void set_modelid();
    
    //routines for specific options:
    double Eurocall(double S, int flag=0);
    double Europut(double S, int flag=0);
    double Bincall(double S, int flag=0);
    double Binput(double S, int flag=0);
};

#endif

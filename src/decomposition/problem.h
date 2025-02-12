#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "cec09.h"
#include "RealLife-MOPs.h"

//// Toolkit includes. //////////////////////////////////////////////////////

#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;




// *********************** CEC 2009 ************************************


void CEC09_F1(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF1(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F2(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF2(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F3(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	CEC09::UF3(&(*(X.begin())), &(*(F.begin())), X.size());
}
void CEC09_F4(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF4(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF5(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F6(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF6(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F7(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF7(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F8(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF8(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F9(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF9(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F10(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF10(&(*(XX.begin())), &(*(F.begin())), X.size());
}

// ---------------- 5 objective test instances ------------------

void CEC09_R2_DTLZ2_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ2_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_R2_DTLZ3_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ3_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_WFG1_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	CEC09::WFG1_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

///////////////////////////////////////wfg problems...
void wfg1(std::vector<double> &F, std::vector<double> &X)
{
	 //---- Get the function name.

	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG1( X, k, M );
}
void wfg2(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG2( X, k, M );
}
void wfg3(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG3( X, k, M );
}
void wfg4(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG4( X, k, M );
}
void wfg5(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG5( X, k, M );
}
void wfg6(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG6( X, k, M );
}
void wfg7(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG7( X, k, M );
}
void wfg8(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG8( X, k, M );
}
void wfg9(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG9( X, k, M );
}
void dtlz1(std::vector<double> &F, std::vector<double> &X)
{

  int k = X.size()- F.size() + 1;

  double g = 0.0 ;

  for (int i = X.size() - k; i < X.size(); i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100 * (k + g);
  for (int i = 0; i < F.size(); i++)
    F[i] = (1.0 + g) * 0.5;

  for (int i = 0; i < F.size(); i++){
    for (int j = 0; j < F.size() - (i + 1); j++)
      F[i] *= X[j];
      if (i != 0){
        int aux = F.size() - (i + 1);
        F[i] *= 1 - X[aux];
      } //if
  }//for

}
void dtlz2(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
   double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } //if
  } // for
}
void dtlz3(std::vector<double> &F, std::vector<double> &X)
{
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  int k = numberOfVariables_ - numberOfObjectives_ + 1;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100.0 * (k + g);
  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } // if
  } //for

}
void dtlz4(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++) {
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(pow(X[j],alpha)*(M_PI/2.0));
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(pow(X[aux],alpha)*(M_PI/2.0));
      } //if
  } // for

}
void dtlz5(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz6(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += pow(X[i],0.1);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz7(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += X[i] ;

  g = 1 + (9.0 * g)/k ;


  for (int i = 0; i < numberOfObjectives_ - 1; i++)
    F[i] = X[i] ;

  double h = 0.0 ;
  for (int i = 0; i < numberOfObjectives_ - 1; i++){
    h+=(F[i]/(1.0+g))*(1 + sin(3.0*M_PI*F[i])) ;
  } //for

  h = numberOfObjectives_ - h ;

  F[numberOfObjectives_ - 1] = (1+g)*h ;	
}
void RWP1(std::vector<double> &F, std::vector<double> &X)
{
        RealLife_MOPs::INJ_4OBJ( &(*(X.begin())),  &(*(F.begin())));
}
void RWP2(std::vector<double> &F, std::vector<double> &X)
{
        RealLife_MOPs::CWD( &(*(X.begin())),  &(*(F.begin())));
}
/*
  IMB1
  two-objective
  recommended: n=10;
  PF
   0 <= f1 <= 0.2
   f2= 1- sqrt(f1)
    0 <= xi <= 1
  
*/
void imb1(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0] > 0.2)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-sqrt(x[0]));
}
/*
  IMB2
  two-objective
  recommended: n=10;
  PF
   0.4 <= f1 <= 0.6
   f2= 1- sqrt(f1)
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x_1)
  
*/
void imb2(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if( x[0]< 0.4 || x[0]>0.6)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB3
  two-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 = 1
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x_1)
  
*/
void imb3(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]< 0.8)
  {
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - sin(0.5*M_PI*x[0]);
	gx += 0.5*(-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5);
  f[1] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB4
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb4(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]<(2.0/3.0))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB5
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb5(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]>(0.5))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5)*cos(M_PI*x[1]*0.5);
  f[1] = (1.0+gx)*cos(M_PI*x[0]*0.5)*sin(M_PI*x[1]*0.5);
  f[2] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB6
  three-objective
  recommended: n=10;
  PF
   f1^2 + f2^2 + f3^2 = 1
    0 <= xi <= 1
  PS: x_j = (x1+x2)/2
  
*/
void imb6(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  if(x[0]>(0.75))
  {
     for(int i = 2; i < n; i++)
     {
	double ti = x[i] - ((x[0]+x[1])*0.5);
	gx += (-0.9*ti*ti + pow(fabs(ti), 0.6)); 
     }
     gx *= 2.0*cos(M_PI*x[0]*0.5);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0 - x[0]);
}
/*
  IMB7
  two-objective
  recommended: n=10;
  PF:  f2=1-sqrt(f1)
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb7(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 1; i < n; i++)
  {
     double ti = x[i] - 0.5, si=x[i]-sin(0.5*M_PI*x[0]);
     if(x[0]>=0.5 && x[0]<=0.8)
       gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-sqrt(x[0]));
}
/*
  IMB8
  two-objective
  recommended: n=10;
  PF:  f2=1-f1
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb8(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
     for(int i = 1; i < n; i++)
     {
	double ti = x[i] - 0.5, si=x[i]-sin(0.5*M_PI*x[0]);
	if(x[0]>=0.5 && x[0]<=0.8)
	  gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
	else 
	  gx += pow(fabs(ti), 0.6);
     }
  f[0] = (1.0+gx)*x[0];
  f[1] = (1.0+gx)*(1.0-x[0]);
}
/*
  IMB9
  two-objective
  recommended: n=10;
  PF:  f2=sqrt(1-f1)
    0 <= xi <= 1
  PS: x_j = 0.5
  
*/
void imb9(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 1; i < n; i++)
  {
     double ti = x[i] - 0.5, si = x[i]-sin(0.5*M_PI*x[0]);
     if(x[0]>=0.5 && x[0]<=0.8)
       gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*cos(M_PI*x[0]*0.5);
  f[1] = (1.0+gx)*sin(M_PI*x[0]*0.5);
}
/*
  IMB10
  three-objective
  recommended: n=10;
  PF:  f1+f2+f3=1
    0 <= xi <= 1
  PS: x_j = sin(0.5*M_PI*x1)
  
*/
void imb10(std::vector<double> &f, std::vector<double> &x)
{
  int n = x.size();
  double gx = 0.0;
  for(int i = 2; i < n; i++)
  {
     double ti = x[i] - x[0]*x[1], si=x[i]-((x[0]+x[1])/2.0);
     if((x[0]>=0.2 && x[0]<=0.8) && ( x[1] >=0.2 && x[1]<=0.8))
       gx += 2.0*(-0.9*si*si + pow(fabs(si), 0.6)); 
     else 
       gx += pow(fabs(ti), 0.6);
  }
  f[0] = (1.0+gx)*x[0]*x[1];
  f[1] = (1.0+gx)*x[0]*(1.0-x[1]);
  f[2] = (1.0+gx)*(1.0-x[0]);
}
#endif

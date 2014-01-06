#pragma once
#include "general.h"
typedef struct {
  double *freq;
  double *like0;
  double *likeF;
  double *F;
}funkyHWE;


class hwe:public general{
public:
  int doHWE;
  //none optional stuff
  gzFile outfileZ;
  hwe(const char *outfiles,argStruct *arguments,int inputtype);
  ~hwe();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

  void estHWE(double *x,double *loglike,int nInd);
  double HWE_like(double *x,double *loglike,int nInd);
  void HWE_EM(double *x,double *loglike,int nInd);

};


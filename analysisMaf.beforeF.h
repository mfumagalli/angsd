#pragma once
#include "general.h"

class frequency:public general{
private:
  char *refName;
  char *ancName;
  float EM_START;
  float minLRT;
  int minKeepInd;
  float defaultMinLRT;
  float defualtMinMaf;

  void prepPrint(funkyPars *pars);

  gzFile outfileZ;
  FILE *outfile;
  //none optional stuff
  //FILE *outfile;
  int doMaf;
  int emIter;
  int doPost;
  int doSNP;
  int doMajorMinor;
  int GL;
  int minInd;
  float minMaf;

  //do gz compression on the fly
  int doZ;


public:

  void run(funkyPars  *pars);
  void likeFreq(funkyPars *pars);
  void postFreq(funkyPars *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  frequency(const char *outfiles,argStruct *arguments,int inputtype);
  ~frequency();
  static double likeFixedMinor_bfgs(double *loglikes,int numInds);
  //other stuff
    static double likeFixedMinor_wrapper(const double *para,const void *dats);
  double emFrequencyNoFixed(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int minor);
  static float likeNoFixedMinor(float p,double *logLikes,int numInds,int major);
  static double likeNoFixedMinor_wrapper(const double *para,const void *dats);
  static double likeNoFixedMinor_bfgs(double *loglikes,int numInds,int major);

  static double emFrequency(double *loglike,int numInds, int iter,double start,int *keep,int keepInd);
  static double likeFixedMinor(double p,double *logLikes,int numInds);
};



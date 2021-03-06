#pragma once
#include "shared.h"
#include "general.h"
#include <list>
typedef struct{
  int major;
  int minor;
  suint *aLine;
}point;

typedef std::list<point> myList;


class error:public general{
private:
   float EM_START;
   int emIter;

  int doCounts;
  FILE *outfile1;
  FILE *outfile2;
  
  int minSites;
  myList aList;
  int doError;
  void consolidate(funkyPars *p);
  char *errorFname;
  double **errors;
  double minPhat;
  double eps;
public:
  //none optional stuff

  
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void addDefault(funkyPars *pars);
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);

  error(const char *outfiles,argStruct *arguments,int inputtype);
  ~error();
  //other stuff 
  
  void setMajorMinor(suint **counts,char *major,char *minor,int nInd,int nSites);
  void likeFixedMinorError_bfgs_tsk2(double **errors,int numSites,int nInd,int *major,int *minor,suint **counts,int emIter,double EM_START);// same as above but only uses the parameters in function.
  static double **** generateErrorPointers(double **error,int nG,int nA);
  static float *logfact(int len ) ;
  static void killGlobalErrorProbs(double ****errorProbs);
  double likeFixedMinorError_wrapper(const double *para,const void *dats);

};

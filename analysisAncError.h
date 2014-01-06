#include "general.h"

class ancErr:public general{
private:
  char *refName;
  char *ancName;
  int doAncError;
  int nInd;
  int sample;
  FILE *outfile;
  FILE *outfile2;
  int currentChr;

  size_t **alleleCounts; //[ind][125]; 
  size_t **alleleCountsChr; //[ind][125]; 

  void model1(funkyPars *pars);
  void model2(funkyPars *pars);

public:
  ancErr(const char *outfiles,argStruct *arguments,int inputtype);
  ~ancErr();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};

#pragma once
#include "general.h"


typedef struct{
  int **dat;
}genoCalls;


class callGenotypes:public general{
private:
  int doGeno;
  float postCutoff;
  gzFile outfileZ;
  int geno_minDepth;
public:
  callGenotypes(const char *outfiles,argStruct *arguments,int inputtype);
  ~callGenotypes();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void getGeno(funkyPars *pars);
  void printGeno(funkyPars *pars);
  
};

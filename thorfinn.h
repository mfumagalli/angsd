#pragma once
#include "general.h"

class thorfinn:public general{
private:
  int doThorfinn;
  char *refName;
  FILE *outfile;
public:
  static double *indF;
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  thorfinn(const char *outfiles,argStruct *arguments,int inputtype);
  ~thorfinn();
};



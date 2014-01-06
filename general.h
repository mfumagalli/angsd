#pragma once
//DEFAULT TRIMMING
#define MINQ 13


#include "shared.h"
#include "bams.h"
#include "argStruct.h"
class general{
public:
  static aHead *header;//contains the header of a single bam;
  static std::map<char *,int,ltstr> *revMap;
  int index;
  static int tot_index;
  //  virtuel general()
  virtual void run(funkyPars *f)=0;
  virtual void print( funkyPars *f)=0;
  virtual void printArg(FILE *fp)=0;
  virtual void clean(funkyPars *f)=0;
  general(){index=tot_index++;};
  virtual ~general(){};
};


general **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments);

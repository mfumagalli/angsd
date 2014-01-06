#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "mrStruct.h"
#include "analysisFunction.h"

class reader_sim{
private:
  int nInd;
  int nInd2;//<-used when slicing out samples, totally addhoc, should be generealized
  int from;
  int to;
  gzFile gz;
  void getOptions(argStruct *arguments);
  char *fname;
  void printArg( FILE *fp);
public:
  void init(argStruct *as);
  mr::funkyPars *fetch(int chunkSize);
  void close() {gzclose(gz);}
  reader_sim(){fname=NULL;nInd=0;};
  ~reader_sim(){free(fname);}
};

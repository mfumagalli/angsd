#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "mrStruct.h"
#include "analysisFunction.h"

class reader_sim{
private:
  int nInd;
  gzFile gz;
  void getOptions(argStruct *arguments);
  char *fname;
  void printArg( FILE *fp);
public:
  void init(argStruct *as);
  mr::funkyPars *fetch(int chunkSize);
  void close() {gzclose(gz);}
  reader_sim(){fname=strdup("dummy");nInd=0;};
};

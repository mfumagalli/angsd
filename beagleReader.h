#pragma once

#include <zlib.h>

#include "mrStruct.h"
class beagle_reader{
private:
  char *fname;
  int fexists(const char* str);
  gzFile openBeagleFile;

  int intName; // intrepret SNP name as chr_pos
  void getOptions(argStruct *arguments);
  void printArg(FILE *fp);
public:
  beagle_reader(argStruct *arguments);
  int nInd;
  mr::funkyPars *fetch(int chunkSize);
};


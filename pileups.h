#pragma once
#include <map>
#include <cstring>
#include <vector>

#include "pileup.h"
#include "mrStruct.h"
#include "argStruct.h"

#include "analysisFunction.h"


class pileups{
private:
  int nInd; //<- used for selecting a subset. 
  char *faifile;
  char *lStart;
  char *lStop;
  char *fnames;
  std::vector<char *> filenames;

  binInput *bobj;
  txtInput *tobj;
  glfClass **obj;
  int *keepGoing;
  int filesLeft;

  lociPileup ultraStart;
  lociPileup ultraStop;

  char**chrnames;
  void print_lociPileup(const lociPileup &l);
  void parseLoci(lociPileup &l,const char *str,const cMap &cmap);
  cMap faiIndex;
  //used for killing the program and flushing the buffers. program will still segfault.
  //initialize signal handler
  cMap buildMap(const char * fname);
  mr::funkyPars *collapse(pileworld::aMap &asso);
  int nFiles;
  int type;
  void getOptions(argStruct *arguments);
  void printArg( FILE *fp);
public:
  //  pileups(char *faifile,char *lStart,char*lStop,std::vector<char *> &filenames,int type);
  pileups(argStruct *arguments);
  ~pileups();
  mr::funkyPars *fetch(int nLines,int chunkSize);
};


#pragma once
#include <map>
#include "pileups.h"//used for old style pileups
#include "soapMaster.h" //used for soapfiles
#include "reader_sim.h"
#include "tglfs.h" //used for flat glf files
#include "beagleReader.h"
#include  "mrStruct.h"
#include "argStruct.h"


class multiReader{
private:
  int nLines;
  int chunkSize;
  reader_sim *mysims;
  soapMaster *sm;
  pileups *pl;
  tglfs *mytglfs;
  beagle_reader *bglObj;
  
  int type;//0=soap;1=glf;2=glfclean,3=tglf,4=simfiles
  void getOptions(argStruct *arguments);
  void printArg( FILE *fp);
  std::map<char*,int> mMap;
  char *fainame;
public:
  multiReader(int type_a,argStruct *arguments);
  mr::funkyPars *fetch();
  ~multiReader();
};

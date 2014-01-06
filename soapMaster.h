#pragma once

#include "soapReader.h"

class soapMaster{
private:
  soapReader *obj;
  int *keepGoing;
  int filesLeft;
  FILE *posiFilepointer;

  tmpList *suYeonsList;
    void print_loci(const loci &l){
  fprintf(stderr,"\r\t-> loci: (%s:%d)\t",l.chromo,l.position);
  }
  int nFiles;
  //used for soapreading
  int qs;
  int maxHits;
  float eps;
  int strand;
  int minDepth;
  int maxDepth;
  float cutOff;
  int downsample;
public:
  soapMaster(argStruct *arguments);
  ~soapMaster();
  void printArg(FILE *argFile);
  mr::funkyPars *fetch(int nLines,int chunkSize);
};

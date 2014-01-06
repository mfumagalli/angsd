#pragma once

#include "bams.h"


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};



typedef struct{
  int nInd;//number of inds inferred from filelists
  int argc;
  char **argv;
  int inputtype;//
  int *usedArgs; //array of ints telling if args been used
  FILE *argumentFile; //logfile
  const char *outfiles; //prefix output
  aHead *hd;
  std::map<char *,int,ltstr> *revMap;
}argStruct;

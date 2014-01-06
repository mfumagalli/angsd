#pragma once
#include "shared.h"



typedef struct{
  int jobtype;// 0=nothing,1=view,2=uppile
  char *inputfile;//used for single sample jobtype=1 or jobtype-=2
  char *inputfiles;//when using -b argument;
  char *region;//when using -r argument;
  char *regionfile;//when using -rf argument;
  int nInd;//used for selecting the first nInd samples from -b;
  int type;
  int nLines; //number of lines to be read from a bam file
  int show;
  void (*callback)(void *);
}args;

args *getArgsBambi(argStruct *arguments);
//void print(args *a,FILE *fp);
void bamInfo(FILE *fp);

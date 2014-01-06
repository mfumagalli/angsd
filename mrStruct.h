#pragma once

typedef struct {
  char * chromo;
  int position;
}loci;
#define LENS 10000


typedef unsigned int suint;

namespace mr{

  typedef struct {
    //dynamic
    loci *sites;//<- this will become obsolete
    int refId;//<- this will be the future 4486
    int *posi;//<- this will be the future 4486
    
    suint **counts;//[nsites][5xiNind] !!!! nope! only 4
    double **likes;
    double **post;
    double *phat;
    
    char *major;
    char *minor;
    char *ref;
    char *anc;
    
    //primitives
    int numSites;
    int nInd;
    int **depth;
    //extra stuff associated with each analysis module
  }funkyPars;

  funkyPars *allocFunkyPars();
  void deallocFunkyPars(funkyPars *p);
}


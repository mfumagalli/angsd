#pragma once
#include "mrStruct.h"
#include "argStruct.h"
//posfiel = chr, pos, ref, anc
/*
  per glf object
 */
class tglf{
  int id;
  gzFile inGlf;
public:
  const char* filename;
  void fetch(double **liks,int nlines);
  
  void init(int num,const char *fname);
  void close() {
    //printf("closing:%s\n",filename);
    gzclose(inGlf);
  }
};

/*
  container class for all glf object and the position file
 */

class tglfs{
  FILE *posFP;
  tglf *tglfObj;
  mr::funkyPars *readposition(int nlines);
  int nInd;
  char *fnames;
  std::vector<char *> filenames;
  
  char *posfname;
  void getOptions(argStruct *arguments);
  void printArg( FILE *fp);
public:
~tglfs(){ for(int i=0;i<filenames.size();i++)
      free(filenames[i]);
  };
  //constructor
  void init(argStruct *arguments);
  //will call readposition, and fetch on all files
  mr::funkyPars *fetch(int nlines);
  void close(){
    for(int i=0;i<nInd;i++)
      tglfObj[i].close();
    delete [] tglfObj;
    fclose(posFP);
  }
};

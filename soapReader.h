#pragma once

#include <cstring>
#include <map>
#include <list>
#include <zlib.h>

//typedef short unsigned int suint;
typedef unsigned int suint;
struct cmp_loci {
  bool operator()(const loci& first,const  loci& second) const{
    int tmp = std::strcmp(first.chromo, second.chromo);
    if(tmp!=0)
      return tmp<0;
    else
      return first.position<second.position;
  }
};



typedef std::map<int,int> iMap;// inner map (witin IN chromosome map)
typedef std::map<char *,iMap*,ltstr> oMap; // map from chromo -> correct INNER map



#ifdef ZLIB_BUF
#include "kseq.h" 
#define BUF_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUF_SIZE)
#endif

typedef struct{
  loci l;
  suint *inf;
  int major;
  int minor;
  double phat;
}aSite;


typedef struct{
  loci from;
  suint *to;
}datum;




//these are just typedefs to make life easier
//typedef std::map<loci,int * ,cmp_loci,myallocator<loci> > aMap;
namespace soapworld{
  typedef std::map<loci,suint * ,cmp_loci > aMap;
}
typedef std::list<datum> aList;
typedef std::list<aSite> tmpList;

#define NUMBASES 4



class soapReader{
private:
  const char* fname;
  gzFile file;
  int id;
  int returnVal;
#ifdef ZLIB_BUF
  kstream_t *ks;
  kstring_t *s;
  int dret;
#endif


  aList listBuffer;
  void insertBuffered(soapworld::aMap &asso);
  int validate(const char * str);//used to validate empty buffer arrays

  //not using 
  int doStuff_getAsso_noTarget(soapworld::aMap &asso,int nlines,int qs, int strand,int maxHits);
  int doStuff_useStop_noTarget(soapworld::aMap &asso,loci &stop,int qs, int strand,int maxHits);

  //new simpler functions
  int doStuff_noTarget(soapworld::aMap &asso,int nlines,int qs, int strand,int maxHits);

public:
  static int numFiles;
  void flush(soapworld::aMap &asso) {insertBuffered(asso);}
  void init(const char *filename);
  
  int doStuff_getAsso(soapworld::aMap &asso,int nlines,int qs, int strand,int maxHits);
  int doStuff_useStop(soapworld::aMap &asso,loci &stop,int qs, int strand,int maxHits);
  int notEof() {return returnVal;}
  void close() {
#ifdef ZLIB_BUF
    ks_destroy(ks);
    free(s->s);
    free(s);
#endif
    gzclose(file);
}
  soapReader(){}//target=NULL; }
  //~soapReader() {//free(s);}
};


bool operator <= (const loci& first,const  loci& second);
bool operator > (const loci& first,const  loci& second);
bool operator == (const loci& first,const  loci& second) ;

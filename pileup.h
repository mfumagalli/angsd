
#pragma once

#include <fstream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include <vector>
#include <stdint.h>




struct cmp_char {
  bool operator()(const char* first,const  char* second) const {
    int tmp = std::strcmp(first, second);
    return tmp<0;
  }
};








typedef std::map<char *,int,cmp_char> cMap; // mapping of the chromoname -> int id

typedef struct {
  int chromo;
  int position;
}lociPileup;





typedef struct {
  //  int* locus; //4 times number of individuals {A,C,T,G} +4 the last 4 will be the sum
  double  *lk; //10 times the number of individuals {AA,Aa,aa}
  short int major; //either ACGT {0...3} a
  short int minor; //either ACGT {0...3} A
  short int ancestral;
  short int derived;
  int ref;
  double phat;
  double emPhat;
  int *datumDepth;
}datumP;


//comparison operator used when comparing lociPileup
struct cmp_lociPileup {
  bool operator()(const lociPileup& first,const  lociPileup& second) const {
      if(first.chromo==second.chromo)
      return first.position<second.position;
      else
       return first.chromo<second.chromo;
  }
};



bool operator <= (const lociPileup& first,const  lociPileup& second);
bool operator > (const lociPileup& first,const  lociPileup& second) ;
bool operator == (const lociPileup& first,const  lociPileup& second) ;


namespace pileworld{
  typedef std::map<lociPileup,datumP,cmp_lociPileup> aMap; // mapping of the reads per locus
}




class glfClass{
protected:
  int id; // representing which file we are using id=0 first file etc
  int buffered; // if we have a buffered element from the last read
  lociPileup last; //the last position
  double likeratios[10]; //all 10 llh, we kept all 10 because we didn't know the min/maj before //from 0.981 we are keeping all likeratios...
  lociPileup ultrastart; //this is the start postion, if we are slicing out a region
  lociPileup ultrastop; //...
  std::ifstream inGlf;
  gzFile gz;
  cMap cmap;
  int depthInd;
  const char* filename; //the name of the file associated with id
  static int total;
public:
  int linenum;
  virtual  int readlines(pileworld::aMap& ,lociPileup ,int )=0;

  virtual void init(int ,const char *,lociPileup ,lociPileup,const cMap &)=0;

  int close() {
    inGlf.close();
    if(gz!=Z_NULL)
      gzclose(gz);
    return linenum;
  }
  int getTotal() {return total;};
  void flush(pileworld::aMap &asso);
};

class txtInput : public glfClass{
public:
  int readlines(pileworld::aMap &asso,lociPileup stop,int nLines);
  void init(int num,const char *fname,lociPileup start,lociPileup stop,const cMap &c);
};

class binInput : public glfClass{
private:
  char chrName[100];
  int curPos;
  int refLen;
  int readHeader;
  int myChrId;
public:
  int readlines(pileworld::aMap &asso,lociPileup stop,int nLines);
  void init(int num,const char *fname,lociPileup start,lociPileup stop,const cMap &c);
};


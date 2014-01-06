
#include "shared.h"
#include "general.h"
#include "analysisFunction.h" 


//these are the major builtin analysis that angsd can perform
#include "analysisKeepList.h"
#include "analysisMajorMinor.h"
#include "analysisMaf.h"
#include "analysisEstError.h"
#include "analysisEstLikes.h"
#include "analysisAsso.h"
#include "analysisHWE.h"
#include "analysisAncError.h"
#include "analysisAbbababa.h"
#include "analysisFasta.h"
#include "analysisCallGenotypes.h"
#include "getFasta.h"//for reading fasta; ancestral and refernce
#include "analysisCount.h" //generate counts from reads
//extrastuff
#include "angsd_realSFS.cpp"
#include "analysisCovar.cpp" 

#include "thorfinn.h"
#include "snpStat.h"

#include "snptools.cpp" //<-implementation of some stuff from snptools. 
#include "hetplas.cpp" //<-implementation of hetero plasmic
#include "writePlink.cpp" //<- dump plink files.

int general::tot_index =0;
aHead *general::header = NULL;
std::map<char *,int,ltstr> *general::revMap = NULL;
general **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments){
  int nit=0;
  //  printHd(hd,stderr);
  //change the number of method when adding a new one
  general **tskStuff =new general*[20];
  tskStuff[nit++] = new filter(arguments);//0
  tskStuff[nit++] = new getFasta(arguments);
  tskStuff[nit++] = new countCls(outfiles,arguments,inputtype);
  tskStuff[nit++] = new error(outfiles,arguments,inputtype);
  tskStuff[nit++] = new likeClass(outfiles,arguments,inputtype);
  tskStuff[nit++] = new majorminor(outfiles,arguments,inputtype);
  tskStuff[nit++] = new frequency(outfiles,arguments,inputtype);
  tskStuff[nit++] = new asso(outfiles,arguments,inputtype);
  tskStuff[nit++] = new hwe(outfiles,arguments,inputtype);
  tskStuff[nit++] = new ancErr(outfiles,arguments,inputtype);//9
  tskStuff[nit++] = new callGenotypes(outfiles,arguments,inputtype);//10
  tskStuff[nit++] = new realSFS(outfiles,arguments,inputtype);
  tskStuff[nit++] = new covar(outfiles,arguments);
  tskStuff[nit++] = new thorfinn(outfiles,arguments,inputtype);
  tskStuff[nit++] = new snpStat(outfiles,arguments,inputtype);
  tskStuff[nit++] = new snptools(outfiles,arguments,inputtype);
  tskStuff[nit++] = new hetplas(outfiles,arguments,inputtype);//16
  tskStuff[nit++] = new plink(outfiles,arguments,inputtype);//17
  tskStuff[nit++] = new abbababa(outfiles,arguments,inputtype);//18
  tskStuff[nit++] = new fasta(outfiles,arguments,inputtype);//19
  //add yours here:


  //don't touch below
  nItem = nit;
  return tskStuff;
}

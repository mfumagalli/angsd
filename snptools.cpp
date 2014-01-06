
#include "analysisFunction.h"
#include "shared.h"


class snptools:public general{
private:
  uint16_t *ebd;
  int curChr;
public:
  int doSnpTools;
  
  snptools(const char *outfiles,argStruct *arguments,int inputtype);
  ~snptools();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
void snptools::printArg(FILE *fp){
  fprintf(fp,"doSnpTools=%d\n",doSnpTools);
  
}
void save_chr(uint16_t *p){

}
void snptools::run(funkyPars *pars){
  if(!doSnpTools)
    return;
  if(curChr!=pars->refId){
    save_chr(ebd);

  }
  
    
}

void snptools::clean(funkyPars *fp){
  if(!doSnpTools)
    return;

  
}

void snptools::print(funkyPars *fp){
  if(!doSnpTools)
    return;
      
}


void snptools::getOptions(argStruct *arguments){
  //default
  doSnpTools=0;

  //from command line
  doSnpTools=angsd::getArg("-doSnpTools",doSnpTools,arguments);
  if(doSnpTools==-999){
    doSnpTools=0;
    printArg(stderr);
    exit(0);
  }
  if(doSnpTools==0)
    return;

  printArg(arguments->argumentFile);

}


snptools::snptools(const char *outfiles,argStruct *arguments,int inputtype){
  curChr =-1;
  getOptions(arguments);
  if(doSnpTools)
    fprintf(stderr,"running doSnpTools=%d\n",doSnpTools);
}

snptools::~snptools(){


}

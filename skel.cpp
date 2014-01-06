/*
  This is a skeleton class for making new object that uses the ANGSD framework

 */


#include "analysisFunction.h"
#include "shared.h"
#include "general.h"

class skel:public general{
private:
  int putSomeStuffHere;
public:
  int doSkel;
  skel(const char *outfiles,argStruct *arguments,int inputtype);
  ~skel();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};

void skel::printArg(FILE *fp){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(fp,"\t-doSkel=%d\n",doSkel);
}

void skel::run(funkyPars *pars){

  if(doSkel==0)
    return ;
   
}

void skel::clean(funkyPars *pars){
  if(doSkel==0)
    return;
    

}

void skel::print(funkyPars *pars){
  if(doSkel==0)
    return;

}

void skel::getOptions(argStruct *arguments){

  doSkel=angsd::getArg("-doSkel",doSkel,arguments);
  printArg(arguments->argumentFile);

  if(doSkel==0)
    return;
 
}


skel::skel(const char *outfiles,argStruct *arguments,int inputtype){
  doSkel =0;
  
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doSkel")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
 
}


skel::~skel(){
  if(doSkel==0)
    return;
  
}



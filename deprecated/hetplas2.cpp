
#include "analysisFunction.h"
#include "shared.h"

double **gen_probs(){
  double **ret =new double*[256];
  for(int i=0;i<256;i++)
    ret[i] = new double[16];

  

  for(int i=0;i<256;i++){
    double p1 = log(1-pow(10,-i/10.0));
    double p3 = log((1-exp(p1))/3.0);
    double val= std::max(p1,p3);
    val=std::min(p1-val,p3-val);
    for(int ii=0;ii<16;ii++){
      ret[i][ii] = val;
    }
    ret[i][0]=ret[i][5]=ret[i][10]=ret[i][15] =0;
  }
  
  return ret;
}


class hetplas:public general{
private:
  double **probs;
  FILE *fp;
  int minQ;
public:
  int doHetPlas;
  
  hetplas(const char *outfiles,argStruct *arguments,int inputtype);
  ~hetplas();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
void hetplas::printArg(FILE *fp){
  fprintf(fp,"--------------\n%s:\n",__FILE__);
  fprintf(fp,"\t-doHetPlas=%d\n",doHetPlas);
  fprintf(fp,"\t-minQ=%d\n",minQ);
}
void hetplas::run(funkyPars *pars){
 
  if(!doHetPlas)
    return;
  assert(pars->chk!=NULL);
  double **res = new double*[pars->numSites];
  for(int s=0;s<pars->numSites;s++)
    if(pars->keepSites[s]){
      res[s] = new double[4*pars->nInd];
      for(int i=0;i<4*pars->nInd;i++)
	res[s][i]=0;
      for(int i=0;i<pars->nInd;i++){
	tNode *tn = &pars->chk->nd[s][i];
	for(int l=0;l<tn->l;l++)
	  for(int ll=0;ll<4;ll++){
	    if(refToInt[tn->seq[l]]!=4&&tn->qs[l]>=minQ)
	      //fprintf(stderr,"i*4+ll=%d qs=%d ofs=%d reftoInt=%d \n",i*4+ll,tn->qs[l],refToInt[tn->seq[l]]*4+ll,refToInt[tn->seq[l]]*4);
	    res[s][i*4+ll] += probs[tn->qs[l]][refToInt[tn->seq[l]]*4+ll];
	  }
      }
	
    }
 
  pars->extras[index] = res;
   
}

void hetplas::clean(funkyPars *pars){
  if(!doHetPlas)
    return;
  double **tmp =(double **) pars->extras[index];
  for(int i=0;i<pars->numSites;i++)
    if(pars->keepSites[i])
      delete [] tmp[i];
      
  delete [] tmp;
  
}

void hetplas::print(funkyPars *pars){
  if(!doHetPlas)
    return;
  double **tmp =(double **) pars->extras[index];
  for(int i=0;i<pars->numSites;i++)
    if(pars->keepSites[i]){
      fprintf(fp,"%s\t%d",header->name[pars->refId],pars->posi[i]);
      double *lik = tmp[i]; 
      for(int ii=0;ii<4*pars->nInd;ii++)
	fprintf(fp,"\t%f",lik[ii]);
      fprintf(fp,"\n");
    }
}


void hetplas::getOptions(argStruct *arguments){
  //from command line
  

  doHetPlas=angsd::getArg("-doHetPlas",doHetPlas,arguments);
 
  if(doHetPlas==0)
    return;

  minQ = angsd::getArg("-minQ",minQ,arguments);
  
  printArg(arguments->argumentFile);

}


hetplas::hetplas(const char *outfiles,argStruct *arguments,int inputtype){
  doHetPlas =0;
  minQ = 13;
  probs=NULL;
 
  
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doHetPlas")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
  if(doHetPlas)
    fprintf(stderr,"running doHetPlas=%d\n",doHetPlas);
  if(doHetPlas){
    probs=gen_probs();
    fp=openFile(outfiles,".hetGL");

  } 
}

hetplas::~hetplas(){
  if(doHetPlas)
    fclose(fp);

}

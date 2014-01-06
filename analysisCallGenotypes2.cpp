/*
  thorfinn@binf.ku.dk april 15 2012
  class to call genotypes
 */

#include "shared.h"
#include <cmath>
#include "analysisFunction.h"

//filename of dumped file
#define GENO ".geno"

//stuff related to genotypecalling
#define GENO_MAJOR_MINOR 0x1 // write the major and the minor
#define GENO_PRINT 2 //print the called genotype 0,1,2
#define GENO_ALLELE 4 //print the called genotype AA, AC AG ...
#define GENO_ALL_POST 8 //print the all posterior
#define GENO_WRITE_POST 16 //write the post of the called genotype
#define GENO_FOR_COVAR 32 //binary dump of the posteriors.

typedef struct{
  int **dat;
}genoCalls;


class callGenotypes:public general{
private:
  int doGeno;
  int doCounts;
  int geno_minDepth;
  float postCutoff;
  FILE *outfile;
  
public:
  callGenotypes(const char *outfiles,argStruct *arguments,int inputtype);
  ~callGenotypes();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void getGeno(funkyPars *pars);
  void printGeno(funkyPars *pars);
  
};
void callGenotypes::printArg(FILE *argFile){
  fprintf(argFile,"-----------------\n%s:\n\n",__FILE__);
  fprintf(argFile,"-doGeno\t%d\n",doGeno);
  fprintf(argFile,"\t1: write major and minor\n");
  fprintf(argFile,"\t2: write the called genotype encoded as -1,0,1,2, -1=not called otherwise counts of derived\n");
  fprintf(argFile,"\t4: write the called genotype directly: eg AA,AC etc \n");
  fprintf(argFile,"\t8: write the posterior probability of all possible genotypes\n");
  fprintf(argFile,"\t16: write the posterior probability of called gentype\n");
  fprintf(argFile,"\t32: write the posterior probability of called gentype as binary\n");
  //  fprintf(argFile,"\t64: write the three posterior probability (Beagle style)\n");
  fprintf(argFile,"\t-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3\n");
  fprintf(argFile,"\t-postCutoff=%f\t(-1 indicates no cutof)\n",postCutoff);
  fprintf(argFile,"\n\tNB When writing the posterior the -postCutoff is not used\n\n");
  fprintf(argFile,"\t-doCounts\t%d\t(Count the number A,C,G,T. All sites, All samples)\n",doCounts);
  fprintf(argFile,"\t-geno_minDepth=%d\t(-1 indicates no cutof)\n",geno_minDepth);
}


void callGenotypes::getOptions(argStruct *arguments){
  //default
  //from command line
  doGeno=angsd::getArg("-doGeno",doGeno,arguments);
  if(doGeno==0)
    return;
  int doMaf,doPost;
  doMaf=doPost= 0;

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  postCutoff=angsd::getArg("-postCutoff",postCutoff,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  geno_minDepth=angsd::getArg("-geno_minDepth",geno_minDepth,arguments);

  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
  
  if(arguments->inputtype!=5&&doPost==0){
    fprintf(stderr,"\n\t-> You need -doPost to call genotypes \n");
    exit(0);
  }

  if(doPost!=0&&doMaf==0){
    fprintf(stderr,"\n\t-> You need -doMaf in order to get posterior probabilities\n");
    exit(0);
  }

  if(geno_minDepth>=0&&doCounts==0){
    fprintf(stderr,"\n\t-> You need -doCounts in order to filter based on per-genotype depth\n");
    exit(0);
  }
}

callGenotypes::callGenotypes(const char *outfiles,argStruct *arguments,int inputtype){
  doGeno=0;
  doCounts=0;
  postCutoff=-1;
  geno_minDepth=-1;

  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doGeno")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doGeno==0)
    return;



  //make output files
  const char* postfix;
  postfix=GENO;
  outfile = openFile(outfiles,postfix);
  FILE *oFileGenos=outfile;

  
}


callGenotypes::~callGenotypes(){
  if(doGeno==0)
    return;
}


void callGenotypes::getGeno(funkyPars *pars){
  //pp is eiter the genotype likelihoods or the posterior probablity
  genoCalls *geno =(genoCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    for( int i =0;i<pars->nInd;i++){
      double maxPP=pars->post[s][i*3+0];
      int maxGeno=0;
      if(pars->post[s][i*3+1]>maxPP){
	maxPP=pars->post[s][i*3+1];
	maxGeno=1;
      }
      if(pars->post[s][i*3+2]>maxPP){
	maxPP=pars->post[s][i*3+2];
	maxGeno=2;
      }
      if(maxPP<postCutoff)
	maxGeno=-1;
      if(doCounts!=0) {
	int geno_Depth = pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3];
	if(geno_Depth<geno_minDepth)
	  maxGeno=-1;
      }
      geno->dat[s][i]=maxGeno;
    }
  }
}



void callGenotypes::printGeno(funkyPars *pars){
  genoCalls *geno =(genoCalls *) pars->extras[index];

  if(pars->keepSites==NULL){
    fprintf(stderr,"\t->[%s] We would expect to have an array of keepsites\n",__FUNCTION__);
    exit(0);
  }

  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    if(!(doGeno&GENO_FOR_COVAR))
      fprintf(outfile,"%s\t%d\t",header->name[pars->refId],pars->posi[s]+1);
    if(doGeno&GENO_MAJOR_MINOR&&!(doGeno&GENO_FOR_COVAR))
      fprintf(outfile,"%c\t%c\t",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
    if(doGeno&GENO_FOR_COVAR){
      fwrite(pars->post[s],sizeof(double),3*pars->nInd,outfile);
      continue;
    }

    for(int i=0;i<pars->nInd;i++){
      if(doGeno&GENO_PRINT)
	fprintf(outfile,"%d\t",geno->dat[s][i]);
      if(doGeno&GENO_ALLELE){
	if(geno->dat[s][i]==0)
	  fprintf(outfile,"%c%c\t",intToRef[pars->major[s]],intToRef[pars->major[s]]);
	else if(geno->dat[s][i]==1)
	  fprintf(outfile,"%c%c\t",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==2)
	  fprintf(outfile,"%c%c\t",intToRef[pars->minor[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==-1)
	  fprintf(outfile,"NN\t");

      }
      if(doGeno&GENO_WRITE_POST)
	fprintf(outfile,"%f\t",pars->post[s][i*3+angsd::whichMax(&pars->post[s][i*3],3)]);
      if(doGeno&GENO_ALL_POST)
	for(int n=0;n<3;n++)
	  fprintf(outfile,"%f\t",pars->post[s][i*3+n]);
    }
    fprintf(outfile,"\n");

  }

}




void callGenotypes::clean(funkyPars *pars){
  if(doGeno==0)
    return;
  genoCalls *geno =(genoCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++)
    delete[] geno->dat[s];
  delete[] geno->dat;

  delete geno;
}


void callGenotypes::print(funkyPars *pars){
  if(doGeno==0)
    return;
  
  genoCalls *geno =(genoCalls *) pars->extras[index];

  printGeno(pars);
  
}

void callGenotypes::run(funkyPars *pars){
  
  if(doGeno==0)
    return;
  
  //allocate genoCall struct
  genoCalls *geno = new genoCalls;
  geno->dat = new int*[pars->numSites];   
  for(int s=0;s<pars->numSites;s++)
    geno->dat[s]=new int[pars->nInd];
  //genoCall struct to pars
  pars->extras[index] = geno;
  

  getGeno(pars);
}


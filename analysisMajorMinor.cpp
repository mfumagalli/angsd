/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd

  
  This class will populate the major/minor entries of the funkypars

  1) set majorminor according to the GL.
  2) set majorminor accordning to the counts of alleles
  3) wont populate the majorminor but will use the information from -filter
  4) set the major to the reference allele
  5) set the major to the ancestral allele


  models used are:
  1)
  Line Skotte, Thorfinn Sand Korneliussen, Anders Albrechtsen. Association testing for next-generation sequencing data using score statistics. Genet Epidemiol. 2012 Jul;36(5):430-7. 

  2)
  Li Y, Vinckenbosch N, Tian G, Huerta-Sanchez E, Jiang T, Jiang H, Albrechtsen A, Andersen G, Cao H, Korneliussen T, et al., 2010. Resequencing of 200 human exomes identifies an excess of low-frequency non-synonymous coding variants. Nat Genet 42:969–972. 


*/


#include <cmath> //<- for log,exp
#include <cassert>
#include <cfloat>
#include "shared.h"
#include "analysisFunction.h"
#include "general.h"
#include "analysisMajorMinor.h"

void majorminor::printArg(FILE *argFile){
  fprintf(argFile,"-------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doMajorMinor\t%d\n",doMajorMinor);
  fprintf(argFile,"\t1: Infer major and minor from GL\n");
  fprintf(argFile,"\t2: Infer major and minor from allele counts\n");
  fprintf(argFile,"\t3: use major and minor from bim file (requires -filter afile.bim)\n");
  fprintf(argFile,"\t4: Use reference allele as major (requires -ref)\n");
  fprintf(argFile,"\t5: Use ancestral allele as major (requires -anc)\n");
}

void majorminor::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;
  //default
  doMajorMinor=0;

  //below is used only validating if doMajorMinor has the data it needs
  int GL=0;
  int doCounts=0;

  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
   
  char *ref = NULL;
  char *anc = NULL;
  ref=angsd::getArg("-ref",ref,arguments);
  anc=angsd::getArg("-anc",anc,arguments);

  if(doMajorMinor==4&&ref==NULL){
    fprintf(stderr,"Must supply reference (-ref) when -doMajorMinor 4");
    exit(0);
  }
  if(doMajorMinor==5&&anc==NULL){
    fprintf(stderr,"Must supply ancestral (-anc) when -doMajorMinor 5");
    exit(0);
  }
  free(ref);
  free(anc);
  GL=angsd::getArg("-GL",GL,arguments);


  /*
    inputtype 
    0) soap
    1) samglf 
    2) samglfClean 
    3) tglf 
    4) simulation type 
    5) beagle/post
    6) mpileup
    7) bam
  */
 if(inputtype==5&&doMajorMinor){
   fprintf(stdout,"Cannot estimate the major and minor based on posterior probabilities\n");
   exit(0);
 }
 if((inputtype==0 || inputtype==6 || inputtype==7)&&GL==0&&doMajorMinor==1){
   fprintf(stderr,"-doMajorMinor 1 is based on genotype likelihoods, you must specify a genotype likelihood model -GL \n");
   exit(0);
 }
 if(doMajorMinor==2&&doCounts==0){
   fprintf(stderr,"-doMajorMinor 2 is based on allele counts, you must specify -doCounts 1\n");
   exit(0);
   
 }
 
}

majorminor::majorminor(const char *outfiles,argStruct *arguments,int inputtype){
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doMajorMinor")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  printArg(arguments->argumentFile);
}

majorminor::~majorminor(){
}

void majorminor::clean(funkyPars *pars){
  if(doMajorMinor){
    delete [] pars->major;
    delete [] pars->minor;

  }
}

void majorminor::print(funkyPars *pars){
}

// Function to swap major/minor reference/ancestral
void modMajorMinor(funkyPars *pars,int doMajorMinor){
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    
    if((doMajorMinor==4 &&pars->ref[s]==4)||(doMajorMinor==5 &&pars->anc[s]==4)){
      pars->keepSites[s]=0;
      continue;
    }
    int maj = pars->major[s];
    int mmin = pars->minor[s];
    

    if(doMajorMinor==4){
      if(pars->ref[s]!=maj&&pars->ref[s]!=mmin){
	//inferred major and minor is not part of the ref;
	pars->keepSites[s] =0;
	continue;
      }
      if(pars->ref[s]==mmin){
	//if reference is the inferre3d minor, then swap
	pars->major[s] = mmin;
	pars->minor[s] = maj;
      }
      //otherwise don't do anything
    }
    if(doMajorMinor==5){
      if(pars->anc[s]!=maj&&pars->anc[s]!=mmin){
	//inferred major and minor is not part of the anc;
	pars->keepSites[s] =0;
	continue;
      }
      if(pars->anc[s]==mmin){
	//if reference is the inferre3d minor, then swap
	pars->major[s] = mmin;
	pars->minor[s] = maj;
      }
      //otherwise don't do anything
    }    
  }
 
}

void majorMinorGL(funkyPars *pars){
  
  float lmax;
  float totalLike;
  int choiceMajor;
  int choiceMinor;
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0){
      pars->major[s] = 4;
      pars->major[s] = 4;
      continue;
    }
    
    //if we have data then estimate the major/minor
    lmax=-FLT_MAX;
    for(int Imajor=0;Imajor<3;Imajor++) {
      for(int Iminor=(Imajor+1);Iminor<4;Iminor++){
	totalLike=0;
	for(int i=0;i<pars->nInd;i++)
	  totalLike+=angsd::addProtect3(pars->likes[s][i*10+angsd::majorminor[Imajor][Imajor]]+log(0.25),
					pars->likes[s][i*10+angsd::majorminor[Imajor][Iminor]]+log(0.5),
					pars->likes[s][i*10+angsd::majorminor[Iminor][Iminor]]+log(0.25)
					);
	if(totalLike>lmax){
	  lmax=totalLike;
	  choiceMajor=Imajor;
	  choiceMinor=Iminor;
	}
      }
    }
    float W0;
    float W1;
    float W2;
    float sum=0;
    for(int i=0;i<pars->nInd;i++){
      W0=exp(pars->likes[s][i*10+angsd::majorminor[choiceMajor][choiceMajor]])*0.25;
      W1=exp(pars->likes[s][i*10+angsd::majorminor[choiceMajor][choiceMinor]])*0.5;
      W2=exp(pars->likes[s][i*10+angsd::majorminor[choiceMinor][choiceMinor]])*0.25;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
    }
    if(sum/pars->nInd<0.5){
      pars->major[s]=choiceMajor;
      pars->minor[s]=choiceMinor;
    }
    else{
      pars->major[s]=choiceMinor;
      pars->minor[s]=choiceMajor;
    }
  }
}

void majorMinorCounts(suint **counts,int nFiles,int nSites,char *major,char *minor,int *keepSites) {
  assert(counts!=NULL);

  for(int s=0;s<nSites;s++){
    if(keepSites==0)
      continue;

    //part one
    //first lets get the sum of each nucleotide
    int bases[4] = {0,0,0,0};
    for(int i=0;i<nFiles;i++)
      for(size_t j=0;j<4;j++){
	bases[j] += counts[s][i*4+j];
	//	fprintf(oFile,"%d\t",bases[j]);
      }


    //now get the major/minor
    int majorA = 0;
    for(int i=1;i<4;i++)
      if (bases[i]>bases[majorA])
        majorA = i;
    //maj is now the major allele (most frequent)
    int temp=0;
    int minorA= majorA;
    for(int i=0;i<4;i++){
      if(majorA==i) //we should check the against the major allele      
        continue;
      else if (bases[i]>temp){
        minorA = i;
        temp=bases[i];
      }
    }
    major[s] = majorA;
    minor[s] = minorA;
  }
}



void majorminor::run(funkyPars *pars){
  if(doMajorMinor==0 || doMajorMinor ==3 )
    return;
  if(doMajorMinor==1 && pars->likes==NULL){
    fprintf(stderr,"[%s.%s():%d] Problem here:\n",__FILE__,__FUNCTION__,__LINE__);
    exit(0);
  }

  
  //allocate and initialize
  pars->major = new char [pars->numSites];
  pars->minor = new char [pars->numSites];
  memset(pars->major,4,pars->numSites);
  memset(pars->minor,4,pars->numSites);
  
  if(doMajorMinor!=2)
    majorMinorGL(pars);
  else if(doMajorMinor==2)
    majorMinorCounts(pars->counts,pars->nInd,pars->numSites,pars->major,pars->minor,pars->keepSites);
  else
    fprintf(stderr,"[%s.%s()%d] Should never happen\n",__FILE__,__FUNCTION__,__LINE__);

  //if user has requested reference/ancestral as major then swap if needed
  if(doMajorMinor==4||doMajorMinor==5)
    modMajorMinor(pars,doMajorMinor);
}

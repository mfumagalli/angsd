/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
 
    
  anders albrecht@binf.ku.dk made this.

  part of angsd
  ans -> anc dec 7 2013, added -ref -ancb
*/

#include <cmath>
#include <cstdlib>

#include "analysisFunction.h"
#include "analysisAncError.h"

void ancErr::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAncError\t%d\n",doAncError);
  fprintf(argFile,"\t1: EM v1 \n");
  fprintf(argFile,"\t2: EM v2 (Under development, don't use)\n");
  fprintf(argFile,"\t-sample\t%d\t(Sampling strategies, for EM v1)\n",sample);
  fprintf(argFile,"\t 0:\t (Use all bases)\n");
  fprintf(argFile,"\t 1:\t (Sample single base)\n");
  fprintf(argFile,"\t 2:\t (Sample first base)\n");
  fprintf(argFile,"Required:\n\t-ref\t%s\t(fastafile containg \'perfect\' sample)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(fastafile containg outgroup)\n",ancName);
  fprintf(argFile,"\nNB: the -ref should be a fasta for a sample where you assume no errors.\nWe measure the difference between the outgroup and your -ref sample.\nThe statistic is then the excess of substitutions between your BAM file and outgroup, compared to the perfect sample.\n");
}

void ancErr::getOptions(argStruct *arguments){

  //from command line
  doAncError=angsd::getArg("-doAncError",doAncError,arguments);
  sample=angsd::getArg("-sample",sample,arguments);
  nInd=arguments->nInd;
  refName=angsd::getArg("-ref",refName,arguments);
  ancName=angsd::getArg("-anc",ancName,arguments);
  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
  if(doAncError){
    if(arguments->inputtype!=0&&arguments->inputtype!=7){
      fprintf(stderr,"Error: bam or soap input needed for -doAncError \n");
      exit(0);
    }
    if(refName==NULL||ancName==NULL){
      fprintf(stderr,"Must supply -ref and -anc for analysis\n");
      exit(0);
    }
  }

}

ancErr::ancErr(const char *outfiles,argStruct *arguments,int inputtype){
  doAncError=0;
  sample=1;
  currentChr=-1;
  refName=NULL;
  ancName=NULL;
  outfile = NULL;
  outfile2= NULL;
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doAncError")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAncError==0)
    return;

  //make output files
  const char* postfix;
  postfix=".ancError";
  outfile = aio::openFile(outfiles,postfix);
  const char* postfix2;
  postfix2=".ancErrorChr";
  outfile2 = aio::openFile(outfiles,postfix2);

  //allocate allele counts
  alleleCounts = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCounts[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCounts[i][j]=0;

  alleleCountsChr = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCountsChr[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCountsChr[i][j]=0;
}


ancErr::~ancErr(){

  if(doAncError==0)
    return;

  if(doAncError==1){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	fprintf(outfile,"%lu\t",alleleCounts[i][j]);
      fprintf(outfile,"\n");
    }
  }

  for(int i=0;i<nInd;i++)
    delete[]  alleleCounts[i];
  delete [] alleleCounts; 

  for(int i=0;i<nInd;i++)
    delete[]  alleleCountsChr[i];
  delete [] alleleCountsChr; 

  if(outfile) fclose(outfile);
  if(outfile2) fclose(outfile2);
}


void ancErr::clean(funkyPars *pars){

}



void ancErr::model2(funkyPars *pars){

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    
    fprintf(outfile,"%d\t%c",pars->posi[s]+1,intToRef[pars->anc[s]]);
    fprintf(outfile,"\t%c",intToRef[pars->ref[s]]);
    
    for(int i=0;i<pars->nInd;i++){
      int allele=4;
      if(pars->chk->nd[s][i].l!=0)
	refToInt[pars->chk->nd[s][i].seq[0]];
      fprintf(outfile,"\t%c",intToRef[allele]);
    }
      fprintf(outfile,"\n");
    }

}


void ancErr::model1(funkyPars *pars){
    if(sample==0){//use all bases
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  for(int j=0;j<pars->chk->nd[s][i].l;j++){
	      alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	      alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	  }
	}
      }
    }
    if(sample==1){//random base
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  if(pars->chk->nd[s][i].l==0)
	    continue;
	  int j = std::rand() % pars->chk->nd[s][i].l;
	  alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	  alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	}
      }
    }
    if(sample==2){//first base
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  if(pars->chk->nd[s][i].l!=0){
	    alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[0]]]++;
	    alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[0]]]++;
	  }
	}
      }
    }


}



void ancErr::print(funkyPars *pars){

  if(doAncError==0)
    return;

  if(doAncError==2){
    model2(pars);
  } else if(doAncError==1) {
    model1(pars);
    if(currentChr==-1)
      currentChr=pars->refId;
    if(currentChr!=pars->refId){
      fprintf(outfile2,"Chr: \t %s\n",header->name[currentChr]);
      for(int i=0;i<nInd;i++){
	for(int j=0;j<125;j++)
	  fprintf(outfile2,"%lu\t",alleleCountsChr[i][j]);
	fprintf(outfile2,"\n");
      }
      for(int i=0;i<nInd;i++)
	for(int j=0;j<256;j++)
	  alleleCountsChr[i][j]=0;
      currentChr=pars->refId;
    }
  }
   
}


void ancErr::run(funkyPars *pars){

  if(doAncError==0)
    return;

}






/*
s<-scan("tmp.ancError")
#pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]


fun<-function(a,p,samp)
  s[25*a+5*p+samp+1]

##perfect vs. ancestry
 sum(fun(1,1,0:3));sum(fun(1,0,0:3))

##sample vs. ancestry
 sum(fun(1,0:3,1));sum(fun(1,0:3,0))

##perfect vs. sample
 sum(fun(0:3,1,1));sum(fun(0:3,0,1))

s<-scan("tmp2.ancError")

fun<-function(a,p,samp)
  s[25*a+5*p+samp+1]

##perfect
 sum(fun(1,1,0:3));sum(fun(1,0,0:3))

##sample
 sum(fun(1,0:3,1));sum(fun(1,0:3,0))


s<-scan("tmp3.ancError")
#pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]


fun<-function(a,p,samp)
  s[25*a+5*p+samp+1]

##perfect vs. ancestry
 sum(fun(1,1,0:3));sum(fun(1,0,0:3))

##sample vs. ancestry
 sum(fun(1,0:3,1));sum(fun(1,0:3,0))

##perfect vs. sample
 sum(fun(0:3,1,1));sum(fun(0:3,0,1))

*/

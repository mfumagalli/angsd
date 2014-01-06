/*
  fix parsing of arugment eg -doMAf ==16 -> doMAf &16 etc
  DRAGON

 */

#include <cassert>
#include <cmath>
#include "bfgs.h"
#include "shared.h"
#include "analysisFunction.h"
#include "kstring.h"//<-used for buffered output
#include "analysisMaf.h"

typedef struct{
  int major;
  double *loglikes;
  int numInds;
}bfgs_vars;


//simple phat estimator from the 200danes article, one site function
double phatFun(suint *counts,int nFiles,double eps,char major,char minor) {

  double pis[nFiles];
  double wis[nFiles];

  if(major==minor)
    return 0.0;

  //if the site is variable
  for(int i=0;i<nFiles;i++){
    int ni = counts[i*4+minor];
    int nt= counts[i*4+minor] + counts[i*4+major]; 
    if(nt==0){//if we dont have any reads for individual 'i'
      pis[i] = 0;
      wis[i] = 0;
      continue;
    }
    pis[i] = (ni-eps*nt)/(nt*(1-2*eps));
    wis[i] = 2.0*nt/(nt+1.0);
  }
  
  double tmp=0;
  for(int i=0;i<nFiles;i++)
    tmp+= (pis[i]*wis[i]);
  
  double phatV = std::max(0.0,tmp/angsd::sum<double>(wis,nFiles));
  return phatV;

}

//simple phat estimator will loop over all sites in the pars
void phatLoop(funkyPars *pars,double eps,double nInd){
  pars->phat = new double[pars->numSites];
  for(int s=0;s<pars->numSites;s++)
    if(pars->keepSites)
      pars->phat[s] = phatFun(pars->counts[s],nInd,eps,pars->major[s],pars->minor[s]);
}

void frequency::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doMaf\t%d\n",doMaf);
  fprintf(argFile,"\t1: BFGS frequency (known major minor)\n");
  fprintf(argFile,"\t2: EM frequency (known major minor)\n");
  fprintf(argFile,"\t4: BFGS frequency (unknown major minor)\n");
  fprintf(argFile,"\t8: EM frequency (unknown major minor)\n");
  fprintf(argFile,"\t16: Frequency from genotype probabilities\n");
  fprintf(argFile,"\t32: AlleleCounts based method (known major minor)\n");
  fprintf(argFile,"\t-doSNP\t%d\n",doSNP);
  fprintf(argFile,"\t-minMaf\t%f %d\n",minMaf,filtMaf);
  fprintf(argFile,"\t-minLRT\t%f %d\n",minLRT,filtLrt);
  fprintf(argFile,"\t-ref\t%s\n",refName);
  fprintf(argFile,"\t-anc\t%s\n",ancName);
  fprintf(argFile,"\t-eps\t%f [Only used for -doMaf &32]\n",eps);
  fprintf(argFile,"-doPost\t%d\t(Calculate posterior prob 3xgprob)\n",doPost);
  fprintf(argFile,"\t1: Using frequency as prior\n");
  fprintf(argFile,"\t2: Using uniform prior\n");
  fprintf(argFile,"-beagleProb\t%d (Dump beagle style postprobs)\n",beagleProb);
  fprintf(argFile,"NB these frequency estimators requires major/minor -doMajorMinor\n");
  fprintf(argFile,"\n");

}

//fancy little function
int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

 
void frequency::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;
  if(inputtype==5)
    inputIsBeagle =1;

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  doSNP=angsd::getArg("-doSNP",doSNP,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);

  if(doMaf==0)
    return;
  double tmp=-1;
  tmp=angsd::getArg("-minMaf",tmp,arguments);
  if(tmp!=-1){
    filtMaf=1;
    assert(tmp<=1||tmp>=0);
    minMaf = tmp; 
  }
  tmp=-1;
  tmp=angsd::getArg("-minLRT",tmp,arguments);
  if(tmp!=-1){
    filtLrt=1;
    minLRT = tmp;
  }

  refName = angsd::getArg("-ref",refName,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);

  if(doMaf&& !isPowerOfTwo((unsigned int) doMaf)){
    fprintf(stderr,"\n[%s] You have selected filters for maf/lrt\n",__FILE__);
    fprintf(stderr,"If you have selected more than one MAF estimator we will choose in following order\n");
    fprintf(stderr,"\t1. knownminor bfgs\n");
    fprintf(stderr,"\t2. knownminor EM\n");
    fprintf(stderr,"\t3. unknownminor bfgs\n");
    fprintf(stderr,"\t4. unknownminor EM\n");
    fprintf(stderr,"\t5. Posterior maf\n\n");
  }
  if(doSNP&&doMaf==0){
    fprintf(stderr,"You've selected snpcalling but no maf estimator please select -doMaf INT\n");
    exit(0);
  }
  if(filtLrt&&(doSNP==0||doMaf==0)){
    fprintf(stderr,"You've selected minLRT threshold please also select -doMaf INT and doSNP INT\n");
    exit(0);
  }
  if(filtMaf &&(doMaf==0)){
    fprintf(stderr,"\nYou've selected minmaf/minlrt and no MAF estimator, choose -doMaf\n\n");
    exit(0);
  }

  int doCounts =0;
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);

  beagleProb=angsd::getArg("-beagleProb",beagleProb,arguments);

  if(doMaf==0 &&doPost==0)
    return;
  if(inputtype!=5&&inputtype!=0&&doMajorMinor==0){
    fprintf(stderr,"You must specify \'-doMajorMinor\' to infer major/minor \n");
    exit(0);
  }

  if(inputtype==5&&doMaf!=16){
    fprintf(stderr,"Only \'-doMaf 16\' can be performed on posterior input\n");
    exit(0);
  }
  if((inputtype==0 || inputtype==6 || inputtype==7)&&GL==0){
    fprintf(stderr,"Error: For sequence data (bam,SOAP) likehoods (-GL) must be specified for frequency estimation\n");
    exit(0);
  }
  if(inputtype!=5&&doMaf==16){
    fprintf(stderr,"Error: doMaf=16 can be performed on posterior input\n");
    exit(0);
  }
  if(doSNP&&doMaf==16){
    fprintf(stderr,"Error: doMaf=16 cannot be used for a likelihood ratio test (doSNP) \n");
    exit(0);
  }
  if(inputtype==5&&doPost){
    fprintf(stderr,"Error: Cannot estimate post (doPost) based on posterior probabilites\n");
    exit(0);
  }
  if(doMaf&32&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts for MAF estimator based on counts\n");
    exit(0);
  }
  if(beagleProb && doPost==0){
    fprintf(stderr,"Must supply -doPost 1 to write beaglestyle postprobs\n");
    exit(0);
  }

  if(doPost!=0 && doMajorMinor==0){
    fprintf(stderr,"Do post requires major and minor: supply -doMajorMinor \n");
    exit(0);
  }
}

//constructor
frequency::frequency(const char *outfiles,argStruct *arguments,int inputtype){
  inputIsBeagle =0;
  beagleProb = 0; //<-output for beagleprobs
  filtLrt=filtMaf =0;
  minMaf =0.01;
  minLRT =24;
  nInd = arguments->nInd;
  eps = 0.001;
  outfileZ2 = Z_NULL;
  outfileZ = Z_NULL;
    
  doMaf=0;
  GL=0;
  doSNP=0;
  doPost=0;

  emIter=100;
  doMajorMinor=0;
  refName = NULL;
  ancName = NULL;
  EM_START = 0.001;

  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doMaf")||!strcmp(arguments->argv[1],"-doPost")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }



  
  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doMaf==0)
    return;
  //make output files
  const char* postfix;
  postfix=".mafs.gz";
  outfileZ = aio::openFileGz(outfiles,postfix,GZOPT);
  if(beagleProb){
    postfix=".beagle.gprobs.gz";
    outfileZ2 = aio::openFileGz(outfiles,postfix,GZOPT);
  }
  //print header
  kstring_t bufstr;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  kputs("chromo\tposition\tmajor\tminor\t",&bufstr);
  if(refName!=NULL||inputtype==3)// inputtyp=3 is tglf, if tglf then -posi has been supplied
    kputs("ref\t",&bufstr);
  if(ancName||inputtype==3||inputtype==4)// inputtyp=3 is tglf, if tglf then -posi has been supplied
    kputs("anc\t",&bufstr);
  
  if(doMaf &1)
    kputs("knownBFGS\t",&bufstr);
  if(doMaf &2)
    kputs("knownEM\t",&bufstr);
  if(doMaf &4)
    kputs("unknownBFGS\t",&bufstr);
  if(doMaf &8)
    kputs("unknownEM\t",&bufstr);
  if(doMaf &16)
    kputs("PPmaf\t",&bufstr);
  if(doMaf &32)
    kputs("phat\t",&bufstr);
  
  if(doSNP){
    if(doMaf &1)
      kputs("pK-BFGS\t",&bufstr);
    if(doMaf &2)
      kputs("pK-EM\t",&bufstr);
    if(doMaf &4)
      kputs("pU-BFGS\t",&bufstr);
    if(doMaf &8)
      kputs("pu-EM\t",&bufstr);
  }
  kputs("nInd\n",&bufstr);
  gzwrite(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;
  if(beagleProb){
    kputs("marker\tallele1\tallele2",&bufstr);
    for(int i=0;i<arguments->nInd;i++){
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
    }
    kputc('\n',&bufstr);
    gzwrite(outfileZ2,bufstr.s,bufstr.l);
  }

  free(bufstr.s);
}


frequency::~frequency(){
  if(outfileZ!=Z_NULL)     gzclose(outfileZ);
  if(outfileZ2!=Z_NULL)    gzclose(outfileZ2);
  free(refName);
  free(ancName);
}



void frequency::print(funkyPars *pars) {
  if(doMaf==0)
    return;
  kstring_t bufstr;
  bufstr.s=NULL; bufstr.l=bufstr.m=0;

  freqStruct *freq =(freqStruct *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    //plugin chr,pos,major,minor
    kputs(header->name[pars->refId],&bufstr);kputc('\t',&bufstr);
    kputw(pars->posi[s]+1,&bufstr);kputc('\t',&bufstr);
    kputc(intToRef[pars->major[s]],&bufstr);kputc('\t',&bufstr);
    kputc(intToRef[pars->minor[s]],&bufstr);kputc('\t',&bufstr);


    //plugin ref, anc if exists
    if(pars->ref!=NULL)
      {kputc(intToRef[pars->ref[s]],&bufstr);kputc('\t',&bufstr);}
    if(pars->anc!=NULL)
      {kputc(intToRef[pars->anc[s]],&bufstr);kputc('\t',&bufstr);}

    
    
    if(doMaf &1)
      ksprintf(&bufstr,"%f\t",freq->pml[s]);
    if(doMaf &2)
      ksprintf(&bufstr,"%f\t",freq->pEM[s]);
    if(doMaf &4)
      ksprintf(&bufstr,"%f\t",freq->pmlun[s]);
    if(doMaf &8)
      ksprintf(&bufstr,"%f\t",freq->pEMun[s]);
    if(doMaf &16)
      ksprintf(&bufstr,"%f\t",freq->freq[s]);
    if(doMaf &32)
      ksprintf(&bufstr,"%f\t",pars->phat[s]);


    if(doSNP){
      if(doMaf &1)
	ksprintf(&bufstr,"%f\t",freq->pmlSNP[s]);
      if(doMaf &2)
	ksprintf(&bufstr,"%f\t",freq->pEMSNP[s]);
      if(doMaf &4)
	ksprintf(&bufstr,"%f\t",freq->pmlunSNP[s]);
      if(doMaf &8)
	ksprintf(&bufstr,"%f\t",freq->pEMunSNP[s]);
    }

    kputw(pars->keepSites[s],&bufstr);kputc('\n',&bufstr);
  }

  gzwrite(outfileZ,bufstr.s,bufstr.l);  
  bufstr.l=0;

  if(beagleProb){
    //beagle format
    for(int s=0;s<pars->numSites;s++) {
      
      if(pars->keepSites[s]==0)
	continue;
      //	fprintf(stderr,"keepsites=%d\n",pars->keepSites[s]);
      kputs(header->name[pars->refId],&bufstr);
      kputc('_',&bufstr);
      kputw(pars->posi[s]+1,&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->major[s],&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->minor[s],&bufstr);

      int major = pars->major[s];
      int minor = pars->minor[s];
      assert(major!=4&&minor!=4);
	
      for(int i=0;i<3*pars->nInd;i++) {
	ksprintf(&bufstr, "\t%f",pars->post[s][i]);
      }
      
      kputc('\n',&bufstr);
      gzwrite(outfileZ2,bufstr.s,bufstr.l);
    }
    
  }
  free(bufstr.s);
}



void frequency::clean(funkyPars *pars) {
  if(doMaf==0)
    return;

  freqStruct *freq =(freqStruct *) pars->extras[index];

  //cleaning
  delete [] freq->freq;
  delete [] freq->lrt_snp;//maybe in doSNP
  delete [] freq->pml;
  delete [] freq->pmlun;
  delete [] freq->pEM;
  delete [] freq->pEMun;
  delete [] freq->pmlSNP;
  delete [] freq->pmlunSNP;
  delete [] freq->pEMSNP;
  delete [] freq->pEMunSNP;
  delete freq;

  
  if(pars->post!=NULL){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->post[i];
    delete [] pars->post;
  }


}

freqStruct *allocFreqStruct(){
  freqStruct *freq = new freqStruct;
  
 freq->pml=NULL;
 freq->pmlun= NULL;
 freq->pEM=NULL;
 freq->pEMun=NULL;
 freq->pmlSNP = NULL;
 freq->pmlunSNP = NULL;
 freq->pEMSNP = NULL;
 freq->pEMunSNP = NULL;
 freq->freq = NULL;
 freq->lrt_snp = NULL;
 return freq;
}

void frequency::run(funkyPars *pars) {
 
  if(doMaf==0&&doPost==0)
    return;
  freqStruct *freq = allocFreqStruct();
  pars->extras[index] = freq;

  if(doMaf!=0) {

    if(doMaf&32)
      phatLoop(pars,eps,nInd);
    
    if(doMaf&16)
      postFreq(pars,freq);
    if(doMaf % 16){
      likeFreq(pars,freq);
    }
    
    if(freq->freq==NULL){
      fprintf(stderr,"%s:Frequency not initilized\n",__FUNCTION__);
      exit(1);
    }
    
    for(int s=0;s<pars->numSites;s++){
      //   fprintf(stderr,"keepSites[%d]=%d\n",s,pars->keepSites[s]);
      if(pars->keepSites[s]==0)
	continue;

      if(filtMaf==1){
	if(freq->freq[s] < minMaf)
	  pars->keepSites[s]=0;
	else if(freq->freq[s] > 1 - minMaf)
	  pars->keepSites[s]=0;
      }
      if(filtLrt==1 && (doSNP&&(freq->lrt_snp[s] < minLRT)))
      	pars->keepSites[s]=0;

    }
  }
  if(doPost){
    if(pars->likes==NULL){
      fprintf(stderr,"[%s.%s()] likelihoods missing\n",__FILE__,__FUNCTION__);
      exit(0);
    }
    double **post = new double*[pars->numSites];
    double **like=angsd::get3likes(pars);
    for(int s=0;s<pars->numSites;s++){
      post[s] = new double [3*pars->nInd];

      if(pars->keepSites[s]==0){
	delete [] like[s];
	continue;
      }if(doPost==1) {//maf prior
	double freqEst=freq->freq[s];
	for(int i=0;i<pars->nInd;i++){
	  //	  fprintf(stderr,"[%d]\nlik= %f %f %f\n",i,like[0][i*3+0],like[0][i*3+1],like[0][i*3+2]);
	  post[s][i*3+2]=like[s][i*3+2]+2*log(freqEst);
	  post[s][i*3+1]=like[s][i*3+1]+log(2)+log(1-freqEst)+log(freqEst);
	  post[s][i*3+0]=like[s][i*3+0]+2*log(1-freqEst);
	  //	  fprintf(stderr,"likmod= %f %f %f\n",post[0][i*3+0],post[0][i*3+1],post[0][i*3+2]);
	  double norm = angsd::addProtect3(post[s][i*3+0],post[s][i*3+1],post[s][i*3+2]);
	  post[s][i*3+0]=exp(post[s][i*3+0]-norm);
	  post[s][i*3+1]=exp(post[s][i*3+1]-norm);
	  post[s][i*3+2]=exp(post[s][i*3+2]-norm);
	  //	  fprintf(stderr,"%f %f %f\n",post[0][i*3+0],post[0][i*3+1],post[0][i*3+2]);
	}
	  
	//	exit(0);
 

      }
      else if(doPost==2){//uniform prior
	for(int i=0;i<pars->nInd;i++){
	  double norm = angsd::addProtect3(like[s][i*3+0],like[s][i*3+1],like[s][i*3+2]);
	  for(int g=0;g<3;g++)
	    post[s][i*3+g]=exp(like[s][i*3+g]-norm);
	}
      }
      else{
	fprintf(stderr,"[%s] doPost must be 1 or 2 \n",__FUNCTION__);
	exit(0);
      
      }
      delete [] like[s];
    }
    pars->post = post; 
    delete[] like;
   
  }
  
}


/*
Estimate the allele frequency from the posterior probabilities and add them to funkyPars
*/
void frequency::postFreq(funkyPars  *pars,freqStruct *freq){
  double *returnFreq = new double[pars->numSites]; 

  for(int s=0;s<pars->numSites;s++){
    returnFreq[s]=0;
    for(int i=0;i<pars->nInd;i++){
      returnFreq[s]+=pars->post[s][i*3+1]+2*pars->post[s][i*3+2];
    }
    returnFreq[s] = returnFreq[s]/(pars->nInd*2);
    //fprintf(stderr,"returnfreq %f\n",returnFreq[s]);
  }
  freq->freq=returnFreq;
}


void frequency::likeFreq(funkyPars *pars,freqStruct *freq){//method=1: bfgs_known ;method=2 em;method=4 bfgs_unknown

  //here only the likelihoods for the three genotypes are used. 
  double **loglike = NULL;
  if(inputIsBeagle==1)
    loglike= pars->likes;
  else
    loglike=angsd::get3likesRescale(pars);
  assert(loglike!=NULL);
  //the pml frequencies 
  double *pml=NULL;
  double *pEM =NULL;
  double *pmlun =NULL;
  double *pEMun =NULL;
  double *pmlSNP=NULL;
  double *pEMSNP =NULL;
  double *pmlunSNP =NULL;
  double *pEMunSNP =NULL;
  double *returnFreq = new double[pars->numSites]; //return value not to be deleted
  double *lrt_snp = new double[pars->numSites]; //return value not to be deleted
  //  fprintf(stderr,"method:=%d\n",method);
  if(doMaf &1){ 
    pml = new double[pars->numSites]; // numeric optimisation
    if(doSNP)
      pmlSNP = new double[pars->numSites]; // numeric optimisation
  }if(doMaf &2){
    pEM =new double[pars->numSites]; //em algorithm
    for(int i=0;i<pars->numSites;i++)
      pEM[i] = 0.0;
    if(doSNP)
      pEMSNP =new double[pars->numSites]; //em algorithm
  }if(doMaf &4){
    pmlun = new double[pars->numSites]; // numeric optimisation
    if(doSNP)
      pmlunSNP = new double[pars->numSites]; // numeric optimisation
  }if(doMaf &8){
    pEMun = new double[pars->numSites]; // EM optimisation
    if(doSNP)
      pEMunSNP = new double[pars->numSites]; // EM optimisation
  }
  // number of individuals with data
  int *keepInd = pars->keepSites;
  int keepList[pars->nInd];  

  //loop though all sites and check if we have data.
  //fprintf(stderr,"keepSites[0] %d\n",pars->keepSites[0]);
  for(int s=0;s<pars->numSites;s++) {
    if(keepInd[s]==0)//if we dont have any information
      continue;
    keepInd[s]=0;//
    for(int i=0 ; i<pars->nInd ;i++) {//DRAGON CHECK THIS
      //fprintf(stderr,"size %d\nind %d\t loglike:%f\t%f\t%f\n",s,i,loglike[s][i*3+0],loglike[s][i*3+1],loglike[s][i*3+2]);
      keepList[i]=1;
      if(loglike[s][i*3+0]+loglike[s][i*3+1]+loglike[s][i*3+2]>-0.0001){
	//	fprintf(stderr,"size %d\nind %d\t loglike:%f\t%f\t%f\n",s,i,loglike[s][i*3+0],loglike[s][i*3+1],loglike[s][i*3+2]);
	keepList[i]=0;
      }
      else{
	keepInd[s]++;
      }

    }

    if(doMaf &1) {

      pml[s]= likeFixedMinor_bfgs(loglike[s],pars->nInd);

      if(doSNP &1)
    	pmlSNP[s] = 2*likeFixedMinor(0.0,loglike[s],pars->nInd)-2*likeFixedMinor(pml[s],loglike[s],pars->nInd);
    }

    if( doMaf &2) {

      if(pars->phat!=NULL)
	pEM[s]=emFrequency(loglike[s],pars->nInd,emIter,pars->phat[s],keepList,keepInd[s]);
      else
	pEM[s]=emFrequency(loglike[s],pars->nInd,emIter,EM_START,keepList,keepInd[s]);
      
      if(doSNP)
      	pEMSNP[s] = 2*likeFixedMinor(0.0,loglike[s],pars->nInd)-2*likeFixedMinor(pEM[s],loglike[s],pars->nInd);
    }
    if( doMaf &4) {
      pmlun[s]= likeNoFixedMinor_bfgs(pars->likes[s],pars->nInd,pars->major[s]);
      if(doSNP)
	pmlunSNP[s]= 2*likeNoFixedMinor(0.0,pars->likes[s],pars->nInd,pars->major[s])-2*likeNoFixedMinor(pmlun[s],pars->likes[s],pars->nInd,pars->major[s]);
      
    }
    if( doMaf &8 ){
      if(pars->phat!=NULL)
	pEMun[s]=emFrequencyNoFixed(pars->likes[s],pars->nInd,emIter,pars->phat[s],keepList,keepInd[s],pars->major[s],pars->minor[s]);
      else
      	pEMun[s]=emFrequencyNoFixed(pars->likes[s],pars->nInd,emIter,EM_START,keepList,keepInd[s],pars->major[s],pars->minor[s]);
      if(doSNP)
	pEMunSNP[s] = 2*likeNoFixedMinor(0.0,pars->likes[s],pars->nInd,pars->major[s])-2*likeNoFixedMinor(pEMun[s],pars->likes[s],pars->nInd,pars->major[s]);

    }

  }


  freq->pml=pml;
  freq->pEM=pEM;
  freq->pmlun=pmlun;
  freq->pEMun=pEMun;
  freq->pmlSNP=pmlSNP;
  freq->pEMSNP=pEMSNP;
  freq->pmlunSNP=pmlunSNP;
  freq->pEMunSNP=pEMunSNP;

  if(inputIsBeagle!=1){
    for(int i=0;i<pars->numSites;i++)
      delete [] loglike[i];
    delete [] loglike;
  }
 
  //anders: this is less stupid but slower
  //this is abit stupid, but in the case where multiple MAFS has been selected, we use the max(method)
  //thorfinn: this is very stupid
  for(int s=0;s<pars->numSites;s++){
    if(doMaf &8)
      returnFreq[s]=pEMun[s];
    else if(doMaf &4)
      returnFreq[s]=pmlun[s];
    else if(doMaf &2 )
      returnFreq[s]=pEM[s];
    else if(doMaf &1 )
      returnFreq[s]=pml[s];
  }
  //thorfinn april 16 2012
  for(int s=0;s<pars->numSites;s++){
    if((doMaf &8) && doSNP)
      lrt_snp[s]=pEMunSNP[s];
    else if((doMaf &4) && doSNP)
      lrt_snp[s]=pmlunSNP[s];
    else if((doMaf &2) && doSNP )
      lrt_snp[s]=pEMSNP[s];
    else if((doMaf &1) && doSNP)
      lrt_snp[s]=pmlSNP[s];
  }

  freq->freq=returnFreq;
  freq->lrt_snp = lrt_snp;//thorfinn add 16april 2012

}

double frequency::likeFixedMinor_bfgs(double *loglikes,int numInds){

  bfgs_vars *bv = new bfgs_vars;
  bv->loglikes = loglikes;
  bv->numInds = numInds;
  double start[1] = {0.01};
  double lbound[1] = {0};
  double ubound[1] = {1};
  int lims[1] = {2};
  double res=findmax_bfgs(1,start,bv,&likeFixedMinor_wrapper,NULL,lbound,ubound,lims,-1);

  delete bv;
  return start[0];

}
///* unstaticing when possible
double frequency::likeFixedMinor_wrapper(const double *para,const void *dats){
  bfgs_vars *bf =(bfgs_vars *) dats;
  double *loglikes = bf->loglikes;
  int numInds = bf->numInds;
  double res = likeFixedMinor(para[0],loglikes,numInds);
  return res;
}



double frequency::emFrequencyNoFixed(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int minor){

  int iMinor[3];
  int n=0;
  for(int i=0;i<4;i++){
    if(i!=major){
      iMinor[n]=i;
      n++;
    }
  }

  double **loglikeGeno;
  loglikeGeno =new double*[3];

  for(int j=0;j<3;j++){
    loglikeGeno[j] = new double[numInds*3]; 
    for(int i=0;i<numInds;i++){
      loglikeGeno[j][i*3+0]=loglike[i*10+angsd::majorminor[major][major]];
      loglikeGeno[j][i*3+1]=loglike[i*10+angsd::majorminor[major][iMinor[j]]];
      loglikeGeno[j][i*3+2]=loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]];
    }
  }

  float W0[3];
  float W1[3];
  float W2[3];
  start=0.4;
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;
  
  int it=0;
  float weight[3];
  float normWeight[3];
  float norm;
  int correctInd;

  for(it=0;it<iter;it++){
    for(int j=0;j<3;j++)
      weight[j]=-likeFixedMinor(p,loglikeGeno[j],keepInd);
    norm=angsd::addProtect3(weight[0],weight[1],weight[2]);
    for(int j=0;j<3;j++)
      normWeight[j]=exp(weight[j]-norm);
    sum=0;
    correctInd=keepInd;
    for(int i=0;i<numInds;i++){
      if(keep[i]==0)
        continue;
      for(int j=0;j<3;j++){
	W0[j]=exp(loglike[i*10+angsd::majorminor[major][major]])*pow(1-p,2);
	W1[j]=exp(loglike[i*10+angsd::majorminor[major][iMinor[j]]])*2*p*(1-p);
	W2[j]=exp(loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]])*(pow(p,2));
	if(W0[j]+W1[j]+W2[j]<0.0000001)
	  continue;
	sum+=(W1[j]+2*W2[j])/(2*(W0[j]+W1[j]+W2[j]))*normWeight[j];
      }
    }

    p=sum/correctInd;
       
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
       
    temp_p=p;
    if(std::isnan(p)){
      fprintf(stderr,"[%s] caught nan will exit\n",__FUNCTION__);
      for(int j=0;j<3;j++)
	fprintf(stderr,"w%d %f ",j,weight[j]);
      fflush(stderr);

      for(int i=0;i<numInds;i++){
	for(int j=0;j<10;j++){
	  if(keep[i]==0)
	    continue;
	  fprintf(stderr,"%f\t",loglike[i*10+j]);
	}
	fprintf(stderr,"\n");
      }
      fprintf(stderr,"%f | %f %f %f |%f | %d %d | \n",p, normWeight[0], normWeight[1], normWeight[2],angsd::addProtect3(weight[0],weight[1],weight[2]),major,minor);

    }
  }
  
  for(int j=0;j<3;j++)
    delete[] loglikeGeno[j]; 
  delete[] loglikeGeno;
  
  return(p);
}



double frequency::likeNoFixedMinor(double p,double *logLikes,int numInds,int major){
  //logLikes contains the 10 likelihoods for one site
  //bfgs must be used
  float partialLike[3];
  for(int j=0;j<3;j++)
    partialLike[j]=0;
  float totalLike=0;
  int iMinor[3];
  int n=0;
  for(int i=0;i<4;i++){
    if(i!=major){
      iMinor[n]=i;
      n++;
    }
  }

  for(int i=0;i<numInds;i++){
    for(int j=0;j<3;j++) 
      partialLike[j]+=
	angsd::addProtect3(
			   logLikes[i*10+angsd::majorminor[major][major]]+2*log(1-p),
			   logLikes[i*10+angsd::majorminor[major][iMinor[j]]]+log(2)+log(p)+log(1-p),
			   logLikes[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]]+2*log(p));
  }
  totalLike=angsd::addProtect3(partialLike[0],partialLike[1],partialLike[2]);
  return -totalLike;

}



double frequency::likeNoFixedMinor_wrapper(const double *para,const void *dats){
  bfgs_vars *bf = (bfgs_vars*) dats;

  return likeNoFixedMinor(para[0],bf->loglikes,bf->numInds,bf->major);
}

double frequency::likeNoFixedMinor_bfgs(double *loglikes,int numInds,int major){
  bfgs_vars bv;
  bv.loglikes = loglikes;
  bv.numInds = numInds;
  bv.major=major;
  double start[1] = {0.01};
  double lbound[1] = {0};
  double ubound[1] = {1};
  int lims[1] = {2};
  double res=findmax_bfgs(1,start,(void*) &bv,likeNoFixedMinor_wrapper,NULL,lbound,ubound,lims,-1);
  return start[0];
}

double frequency::likeFixedMinor(double p,double *logLikes,int numInds){
  double totalLike=0;
  for(int i=0;i<numInds;i++){
    totalLike+=angsd::addProtect3(logLikes[i*3+0]+log(pow(1-p,2)),
				  logLikes[i*3+1]+log(2)+log(p)+log(1-p),
				  logLikes[i*3+2]+log(pow(p,2)));
  }
    return -totalLike;
}

double frequency::emFrequency(double *loglike,int numInds, int iter,double start,int *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;


  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }

  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);

    for(int ii=0;1&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	    fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    exit(0);
  }
  
  return(p);
}


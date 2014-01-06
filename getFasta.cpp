/*
  thorfinn thorfinn@binf.ku.dk 18dec 2012
  assumes zero index positions

  should remove mutexs, and make function thread safe. dragon daemon

  part of angsd using faidx
 
  Actually no need to mutex, except for the first chunk. Lets fix it at some point.
*/

#include <cassert>

#include "analysisFunction.h"
#include "shared.h"

#include "pthread.h"
#include "general.h"
#include "getFasta.h"


//this will initialize our data
perFasta *init(const char *fname){
  //  fprintf(stderr,"[%s.%s:%d] %s\n",__FILE__,__FUNCTION__,__LINE__,fname);
  fprintf(stderr,"\t-> Reading fasta: %s\n",fname);
  perFasta *r= new perFasta;
  r->fai = NULL;
  r->seqs = NULL;
  r->curChr = -1;
  if(pthread_mutex_init(&r->aMut,NULL)){
    fprintf(stderr,"[%s:%s] error initializing mutex \n",__FILE__,__FUNCTION__);
    exit(0);
  }

  if(NULL==(r->fai = fai_load(fname))){
    fprintf(stderr,"[%s:%s] error reading fai file:%s\n",__FILE__,__FUNCTION__,fname);
    exit(0);
  }
  
  return r;
}



void getFasta::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n\n",__FILE__);
  fprintf(argFile,"\t-ref\t%s\t(afile.fasta)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(afile.fasta)\n",ancName);
  fprintf(argFile,"\tNB these fasta files should be indexed 'samtools faidx'\n");
  fprintf(argFile,"\n");

}


//this will destroy a perfasta structure
void destroy(perFasta *f){
  fai_destroy(f->fai);
  free(f->seqs);
  delete f;
  f=NULL;
}

getFasta::~getFasta(){
  if(ref)
    destroy(ref);
  if(anc)
    destroy(anc);
  free(refName);
  free(ancName);
}


void getFasta::getOptions(argStruct *arguments){

  refName = NULL;
  ancName = NULL;
  //from command line
  ancName = angsd::getArg("-anc",ancName,arguments);
  refName = angsd::getArg("-ref",refName,arguments);
  //  fprintf(stderr,"refName=%s\tancName=%s\n",refName,ancName);

  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-ref")|| !strcmp(arguments->argv[1],"-anc") ){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
}

char *getFasta::loadChr(perFasta *f, char*chrName,int chrId){
  //  fprintf(stderr,"[%s] \t->loading chr:%s from faidx=%p curchr=%d\n",__FUNCTION__,chrName,f,f->curChr);
  free(f->seqs);
  f->seqs=NULL;
  //fprintf(stderr,"f->curChr=%d chrId=%d\n",f->curChr,chrId);
  f->curChr = chrId;
  f->seqs = faidx_fetch_seq(f->fai, chrName, 0, 0x7fffffff, &f->chrLen);
  if(f->seqs==NULL){
    fprintf(stderr,"[%s] Error loading fasta info from chr:%s alleles \n",__FUNCTION__,chrName);  
  }
  //  fprintf(stderr,"[%s] done\n",__FUNCTION__);
  return f->seqs;
}



char *getFasta::magic(int refId,int *posi,int numSites,perFasta *f){
  pthread_mutex_lock(&f->aMut);
  assert(refId!=-1);
  //load new chr if different from last or if nothing has been loaded before
  if(f->curChr==-1||refId!=f->curChr){
    //    fprintf(stderr,"[%s] chaning to chr: %d\n",__FUNCTION__,refId);
    //fflush(stderr);
    loadChr(f,header->name[refId],refId);
  }
  //first check that last position is not to long after.
  if(posi[numSites-1] > f->chrLen+200){
    fprintf(stderr,"Trying to access fasta efter end of chromsome+200:%s/%s pos=%d ref_len=%d\n",header->name[f->curChr],header->name[refId],posi[numSites-1],f->chrLen);
  }
  //now loop over all positions. We allow +200 since a read might match to the end of a chr
  char *ret=new char[numSites];
    
  for(int i=0;i<numSites;i++){
    if(f->seqs&&posi[i]<f->chrLen)
      ret[i] = refToInt[f->seqs[posi[i]]];//offseting by one, to make in work like samtools
    else{
      ret[i] = 4;
    }
  }  

  pthread_mutex_unlock(&f->aMut);
  return ret;
}


 getFasta::getFasta(argStruct *arguments){
  refName = NULL;
  ancName = NULL;
  ref=NULL;
  anc=NULL;
   if(arguments->argc==2){
     if(!strcmp(arguments->argv[1],"-ref")||!strcmp(arguments->argv[1],"-anc")){
       printArg(stdout);
       exit(0);
    }else
       return;
   }

   getOptions(arguments);
   printArg(arguments->argumentFile);


   if(ancName!=NULL)
     anc = init(ancName);
  
   if(refName!=NULL)
    ref = init(refName);
   
}

void getFasta::run(funkyPars *pars){
  if(pars->numSites==0)
    return;
  if(ref)
    pars->ref = magic(pars->refId,pars->posi,pars->numSites,ref);
  if(anc)
    pars->anc = magic(pars->refId,pars->posi,pars->numSites,anc);
  
}


void getFasta::print(funkyPars *f){}
void getFasta::clean(funkyPars *f){}

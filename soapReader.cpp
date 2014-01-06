#include <cstdio>
#include <cstdlib>
#include "argStruct.h"
#include "mrStruct.h"
#include "soapReader.h"
#include "analysisFunction.h"

using namespace soapworld;

bool operator <= (const loci& first,const  loci& second) {

  if(strcmp(first.chromo,second.chromo)==0)
    return first.position<=second.position;
  else
    return strcmp(first.chromo,second.chromo)<=0;
}


bool operator > (const loci& first,const  loci& second) {
  return (!(first <= second));
}

bool operator == (const loci& first,const  loci& second) {

    if(first.chromo==second.chromo)
      return first.position==second.position;
    else
      return first.chromo==second.chromo;
}



int soapReader::numFiles =0;

enum {
  A = 0,
  C = 1,
  G = 2,
  T = 3
};



int getOffset(char c){
  if(c=='A')
    return A;
  else if(c=='C')
    return C;
  else if(c=='T')
    return T;
  else if(c=='G')
    return G;
  else 
    fprintf(stderr,"This error should never occur, a nonseq read char occured: %c\n",c);
  exit(0);
  return -1;
}



int soapReader::doStuff_useStop_noTarget(aMap &asso,loci &stop,int qs,int strand,int maxHits) {
  //  validate("usestop");
  if(asso.size()==0){
    fprintf(stderr,"\t-> Huge error asso is zero big\n");
    exit(0);
  }
  //fprintf(stderr,"\t-> soapReader[%d] talking\n",id);
  //make sure if we have something in the buffer, that we input those.

  aMap::iterator it;

  //char *buffer=new char[LENS];
  char *qc; //quality-score
  char *seq;
  int len; //length of read
  int pos; //startposition;
  int keepGoing = 1;

  char *chr=NULL;
#ifndef ZLIB_BUF
  //  fprintf(stderr,"not using buffered\n");
  char buf[LENS];
#endif

  int linesRead =0;
  while(keepGoing ) {
#ifdef ZLIB_BUF
    if(ks_getuntil(ks,'\n',s,&dret)<0){
#else
    if(0==gzgets(file,buf,LENS)){
#endif
      returnVal =0;
      keepGoing =0;
      //      fprintf(stderr,"\neof in file: %s\n",fname);
      fflush(stderr);
      //      fprintf(stderr,"isNULL\n");
      break;
    }
    linesRead++;
    
#ifdef ZLIB_BUF
    chr=strtok(s->s," \t");//skip first column this is the readID
#else
    chr=strtok(buf," \t\n");
#endif
    if(chr==NULL)
      break;
    seq=strtok(NULL," \t");//sequence read;
    qc = strtok(NULL," \t");//this is the qualityscore
    int numHits = atoi(strtok(NULL," \t"));
    if(numHits>maxHits){//skip line if it is a multihit
      //      fprintf(stderr,"skipping\n");
      continue;
    }
    strtok(NULL," \t");//skip pair end alignment;
    len = atoi(strtok(NULL," \t"));
    if(len!=(int)strlen(seq)||len!=(int)strlen(qc)){
      fprintf(stderr,"Error in lengths in file: %s with fileid[%d]\n",fname,id);
      //      fprintf(stderr,"len=%d\tseq=%s\tqc=%s\nbuffer=%s",len,seq,qc,duped);

      exit(0);
    }


    if(strand==-1)//skip plus minus strand;
      strtok(NULL," \t");
    else{
      chr = strtok(NULL," \t");
      if(strand!=chr[0])
      continue;
    }  
      chr = strtok(NULL," \t");
    
    pos = atoi(strtok(NULL," \t"));
    int Start=0;
    if(pos<Start&&Start)
      continue;

    it = --asso.end();
    if((loci){chr,pos} > it->first)
      keepGoing =0;
    
    //lets loop through the quality score 
    for(int i=0; i<len; i++){
      if(qc[i]>qs){
	loci tmp_loci = {chr,pos+i};

	if(tmp_loci>stop){
	  loci l = (loci) {strdup(chr),pos+i};
	  suint  *aa= new suint[4];
	  aa[0]=0;aa[1]=0;aa[2]=0;aa[3]=0;
	  aa[getOffset(seq[i])]++;
	  datum tmp =(datum) {l,aa};
	  listBuffer.push_back(tmp);
	  continue;
	}

	it = asso.find(tmp_loci);
	suint *ary; //the array containing the basecount across all individuals
	
	if(it==asso.end()){
	  ary = new suint[NUMBASES*numFiles];
	  memset(ary,0,NUMBASES*numFiles*sizeof(suint));
	  asso.insert(aMap::value_type((loci){strdup(chr),pos+i}, ary));
	  //	  fprintf(stderr,"inserting pos:%d\n",pos+i);
	}else
	  ary=it->second;
	//fprintf(stderr,"%d\n",NUMBASES*id);
	ary[NUMBASES*id+getOffset(seq[i])]++;
	
      }
    }

  }
  
  //delete [] buffer;
  return linesRead;
}



void soapReader::insertBuffered(aMap &asso){
  aList::iterator lt;
  aMap::iterator mt;
  if(listBuffer.empty())
    return;
  for(lt=listBuffer.begin();lt!=listBuffer.end();lt++){
    //print_loci(lt->from);
    
    mt = asso.find(lt->from); //mt is a map iterator
    suint *ary; //the array containing the basecount across all individuals
    
    if(mt==asso.end()){
      ary = new suint[NUMBASES*numFiles];
      memset(ary,0,NUMBASES*numFiles*sizeof(suint));
      asso.insert(aMap::value_type((lt->from), ary));
    }else{
      free(lt->from.chromo);
      ary=mt->second;
    }
    for(int j=0;j<4;j++)
      ary[id*NUMBASES+j] += lt->to[j];
    
    delete [] lt->to;
  }
  listBuffer.clear();
}

void soapReader::init(const char *filename){
  id=numFiles;
  numFiles++;
  fname = filename;
  returnVal =1;
  

  if(0)
    fprintf(stderr,"trying to open file: %s with id:%d\n",filename,numFiles);
  //  int fexists(const char* str);
  if(aio::fexists(filename)) {
    if(Z_NULL==(file=gzopen(filename,"rb"))){
      fprintf(stderr,"\t-> Error opening file:%s\t with id:%d\n",filename,id);
      fprintf(stderr,"\t-> Have you exceeded the maximum number of allowed open files per process\n");
      exit(0);
    }
#ifdef ZLIB_BUF
    s =(kstring_t *) calloc(1, sizeof(kstring_t));
    ks = ks_init(file);
#endif
    //	    fprintf(stderr,"lengts to allocate:%d\n",len-i);
        //    validate("from init");
    
  }else{
    fprintf(stderr,"\t-> Error opening filepointer for file:%s\t with id:%d\n",filename,id);
    exit(0);
  }

}

int soapReader::doStuff_getAsso_noTarget(aMap &asso,int nLines,int qs,int strand,int maxHits){
  //  validate("from getAsso ");
  //  fprintf(stderr,"size of target:%p\n",target);
  size_t pre_size = asso.size();
  aMap::iterator end_of_asso;
  if(pre_size)
    end_of_asso=--asso.end();
  aMap::iterator it;

  //char *buffer=new char[LENS];
  char *qc; //quality-score
  char *seq;
  int len; //length of read
  int pos; //startposition;
  char *chr=NULL;
#ifndef ZLIB_BUF
  char buf[LENS];
#endif

  int lines_done = 0;
  int linesread =0;
  while(pre_size || (lines_done++<nLines) ) {
    
    //    fprintf(stderr,"lines_done:%d\n",lines_done);
#ifdef ZLIB_BUF
    if(ks_getuntil(ks,'\n',s,&dret)<0){
#else
    if(0==gzgets(file,buf,LENS)){
#endif
      returnVal =0;
      //      fprintf(stderr,"\neof in file: %s\n",fname);
      fflush(stderr);
      //      fprintf(stderr,"\t-> gzGets Z_NULL\n");
      break;
    }
    linesread++;
    
    //    returnVal = (int) strlen(gzgets(file,buffer,LENS));
#ifdef ZLIB_BUF
    chr=strtok(s->s," \t");//skip first column this is the readID
#else
    chr=strtok(buf," \t\n");
#endif

    if(chr==NULL)
      break;
    seq=strtok(NULL," \t");//sequence read;
    qc = strtok(NULL," \t");//this is the qualityscore
    int numHits = atoi(strtok(NULL," \t"));
    if(numHits >maxHits){//skip line if it is a multihit
      //      fprintf(stderr,"skipping\n");
      continue;
    }
    strtok(NULL," \t");//skip pair end alignment;
    len = atol(strtok(NULL," \t"));
    if(len!=(int)strlen(seq)||len!=(int)strlen(qc)){//this shouldnt be needed anymore
      fprintf(stderr,"Error in lengths\n");
      exit(0);
    }
    if(strand==-1)//skip plus minus strand;
      strtok(NULL," \t");
    else{
      chr = strtok(NULL," \t");
      if(strand!=chr[0])
      continue;
    }  
      chr = strtok(NULL," \t");
    
    pos = atoi(strtok(NULL," \t"));
    //    fprintf(stderr,"%d\n",pos);
    //lets loop through the quality score 
    for(int i=0; i<len; i++){
      if(qc[i]>qs){
	loci tmp_loci = {chr,pos+i};

	it = asso.find(tmp_loci);
	suint *ary; //the array containing the basecount across all individuals
	
	if(it==asso.end()){
	  ary = new suint[NUMBASES*numFiles];
	  memset(ary,0,NUMBASES*numFiles*sizeof(suint));
	  asso.insert(aMap::value_type((loci){strdup(chr),pos+i}, ary));
	  //	  fprintf(stderr,"inserting pos:%d\n",pos+i);
	}else
	  ary=it->second;
	ary[NUMBASES*id+getOffset(seq[i])]++;

      }
    }

    if(pre_size){
      it = --asso.end();
      if(it->first > end_of_asso->first)
	pre_size  = 0;

    }
    
  }
  
  //delete [] buffer;

  return linesread;
}


int soapReader::doStuff_noTarget(aMap &asso,int nLines,int qs,int strand,int maxHits){
  loci stop;
  int usingStop =0;
  if(asso.size()!=0){
    usingStop = 1;
    stop = (--asso.end())->first;
  }
  //if usingStop = 1; then we are using the last element in our maps
  //if usingStop = 0; then we will continue to read nLines from the file

  char *qc; //quality-score
  char *seq;
  int len; //length of read
  int pos; //startposition;
  char *chr=NULL;
#ifndef ZLIB_BUF
  char buf[LENS];
#endif

  int linesread =0;
  int keepgoing =1;
  while(keepgoing) {
    linesread++;
#ifdef ZLIB_BUF
    if(ks_getuntil(ks,'\n',s,&dret)<0){
#else
    if(0==gzgets(file,buf,LENS)){
#endif
      returnVal =0;
      fflush(stderr);
      break;
    }
    
#ifdef ZLIB_BUF
    chr=strtok(s->s," \t");//skip first column this is the readID
#else
    chr=strtok(buf," \t\n");
#endif

    if(chr==NULL)
      break;
    seq=strtok(NULL," \t");//sequence read;
    qc = strtok(NULL," \t");//this is the qualityscore
    int numHits = atoi(strtok(NULL," \t"));
    if(numHits >maxHits){//skip line if it is a multihit
      linesread--;//decrement counter;
      continue;
    }
    strtok(NULL," \t");//skip pair end alignment;
    len = atol(strtok(NULL," \t"));
    if(len!=(int)strlen(seq)||len!=(int)strlen(qc)){//this shouldnt be needed anymore
      fprintf(stderr,"Error in lengths\n");
      exit(0);
    }
    if(strand==-1)//skip plus minus strand;
      strtok(NULL," \t");
    else{
      chr = strtok(NULL," \t");
      if(strand!=chr[0])
      continue;
    }  
    chr = strtok(NULL," \t");
    pos = atoi(strtok(NULL," \t"));

    for(int i=0; i<len; i++){
      if(qc[i]>qs) { 
	loci tmp_loci = {chr,pos+i};

	if(usingStop){//if we are after the stop
	  if(tmp_loci>stop){
	    keepgoing = 0;
	    loci l = (loci) {strdup(chr),pos+i};
	    suint  *aa= new suint[4];
	    aa[0]=0;aa[1]=0;aa[2]=0;aa[3]=0;
	    aa[getOffset(seq[i])]++;
	    datum tmp =(datum) {l,aa};
	    listBuffer.push_back(tmp);
	    continue;
	  }
	}

	aMap::iterator it = asso.find(tmp_loci);
	suint *ary; //the array containing the basecount across all individuals
	
	if(it==asso.end()){
	  ary = new suint[NUMBASES*numFiles];
	  memset(ary,0,NUMBASES*numFiles*sizeof(suint));
	  asso.insert(aMap::value_type((loci){strdup(chr),pos+i}, ary));
	  //	  fprintf(stderr,"inserting pos:%d\n",pos+i);
	}else
	  ary=it->second;
	ary[NUMBASES*id+getOffset(seq[i])]++;

      }
    }    
    
    if(linesread>nLines)
      break;
    }
    if(usingStop && asso.size()==0){
      fprintf(stderr,"possible error in indexing\n");
    }
    return linesread;
}



int soapReader::doStuff_useStop(aMap &asso,loci &stop,int qs,int strand,int maxHits) {
  return doStuff_useStop_noTarget(asso,stop,qs,strand,maxHits);

}


int soapReader::doStuff_getAsso(aMap &asso,int nLines,int qs,int strand,int maxHits){
  return doStuff_getAsso_noTarget(asso,nLines,qs,strand,maxHits);

}

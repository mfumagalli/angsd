#include "reader_sim.h"


void reader_sim::getOptions(argStruct *arguments){
  //  fprintf(stderr,"asdf:%p\n",arguments->argumentFile);
  
  fname=angsd::getArg("-sim1",fname,arguments);
  
  if(strcmp(fname,"-999")==0){
    printArg(stdout);
    exit(0);
  }

  nInd = angsd::getArg("-nInd",nInd,arguments);

  from = angsd::getArg("-from",from,arguments);
  to = angsd::getArg("-to",to,arguments);
  if(from!=-1&&to!=-1){
    nInd2=to-from+1;

  }

  if(nInd==0){
    fprintf(stderr,"Must supply -nInd when using simfiles\n");
    printArg(arguments->argumentFile);
  }

  if(nInd2!=-1)
    arguments->nInd = nInd2;
  else
    arguments->nInd = nInd;  
  //  fprintf(stderr,"arguments->nind==%d\n",arguments->nInd);
  //fprintf(stderr,"nind2=%d\n",nInd2);
  //make dummy header and reverse lookup table
  arguments->hd = new aHead;
  arguments->hd->text=NULL;arguments->hd->l_text=0;
  arguments->hd->n_ref =1;
  arguments->hd->name = new char*[1];
  arguments->hd->name[0]= strdup("dummy1");
  arguments->hd->l_name = new int[1];arguments->hd->l_name[0]=strlen(arguments->hd->name[0]);
  arguments->hd->l_ref = new int[1];arguments->hd->l_ref[0] = 300e6;//some random large value, shouldn't be used; DRAGON
  std::map<char*,int,ltstr> *revMap  =  new   std::map<char*,int,ltstr>;
  revMap->insert(std::pair<char*,int>(strdup("dummy1"),0));
  arguments->revMap = revMap;

  
  printArg(arguments->argumentFile);
}


void reader_sim::init(argStruct *as){
  nInd =0;
  nInd2=-1;
  from =to=-1;
  getOptions(as);
  gz = aio::getGz(fname,"r");
}

void reader_sim::printArg(FILE *argFile){
  fprintf(argFile,"reader_sim:\n\n");
  fprintf(argFile,"-sim1\t%s\n",fname);
  fprintf(argFile,"-nInd\t%d\n",nInd);
  fprintf(argFile,"-from\t%d\n",from);
  fprintf(argFile,"-to\t%d (inclusive)\n",to);
  fprintf(argFile,"-nInd2=%d (inferred :to-from+1)\n",nInd2);
  fprintf(argFile,"\n");

}


mr::funkyPars *reader_sim::fetch(int chunkSize){

  mr::funkyPars *r = mr::allocFunkyPars();  
  r->likes=new double*[chunkSize];
  r->posi=new int[chunkSize];
  r->anc = new char[chunkSize];
  memset(r->anc,0,chunkSize);
  r->refId = 0;
  static int pos = 0;
  int i;

  for(i=0;i<chunkSize;i++){
    //    fprintf(stderr,"fechgingigngng:nind=%d nind2=%d\n",nInd,nInd2);
    double *lk = new double [10*nInd];
    size_t bytesRead = gzread(gz,lk,sizeof(double)*10*nInd);
    if(bytesRead==0){
      delete [] lk;
      break;
    }else if(bytesRead!=sizeof(double)*10*nInd){
      fprintf(stderr,"\t-> Error reading full chunk\n");
      exit(0);
    }else{ 
      if(from==-1)
	r->likes[i] = lk;
      else{
	double *lk2=new double[10*nInd2];
	memcpy(lk2,lk+10*from,sizeof(double)*10*nInd2);
	r->likes[i] = lk2;
	delete [] lk;
      }
      r->posi[i] = i+pos;
    }
  }
  if(i==0){
    delete [] r->likes;
    delete [] r->posi;
    delete [] r->anc;
    delete r;
    return NULL;
  }
  if(nInd2==-1)
    r->nInd=nInd;
  else
    r->nInd=nInd2;
  r->numSites=i;
  pos += i;
  
  return r;
}


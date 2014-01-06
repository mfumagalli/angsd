#include <cstdio>
#include <cstdlib>
#include <zlib.h>
#include <vector>
#include <cmath>
#include <fstream>
#include "tglfs.h"
#include "analysisFunction.h"


int bgiToSam[10] = {0,4,5,6,1,7,8,2,9,3};//convert to samtools ordering AA,AC,...


void tglfs::printArg(FILE *argFile){
  fprintf(argFile,"%s:\n\n",__FILE__);
  fprintf(argFile,"-tglf=%s AND -posfile\t%s\n",fnames,posfname);
  fprintf(argFile,"-nInd=%d \n",nInd);
  fprintf(argFile,"\n");
}
void tglfs::getOptions(argStruct *arguments){
  nInd =-1;
  fnames=angsd::getArg("-tglf",fnames,arguments);
  posfname=angsd::getArg("-posfile",posfname,arguments);
  nInd = angsd::getArg("-nInd",nInd,arguments);
  fprintf(stderr,"fnames=%s\n",fnames);
  if(strcmp(fnames,"-999")==0||strcmp(posfname,"-999")==0){
    printArg(stdout);
    exit(0);
  }


  filenames= angsd::getFilenames(fnames,nInd);
  arguments->nInd=filenames.size();
  printArg(stdout);
}



void tglf::init(int num,const char *fname){
  filename = fname;
  id=num;
  //  fprintf(stderr,"glf[%d]=initing with name:%s\n",id,fname);
  inGlf =NULL;
  inGlf = gzopen(fname,"r");
  if(Z_NULL==inGlf){
    fprintf(stderr,"\t-> Your glffilereader object is notworking for fileid[%d]: %s\n",id,filename);
    fprintf(stderr,"\t-> You might need to increase the max number of open files allowed per process\n");
    exit(1);
  }
}

void simplefy(mr::funkyPars *myfunky,int nlines){
  myfunky->likes = new double*[nlines];
  myfunky->sites = new loci[nlines];
  myfunky->ref = new char[nlines];
  myfunky->anc = new char[nlines];
}

void unsimplefy(mr::funkyPars *r,int nlines){
  delete [] r->likes;
  delete [] r->sites;
  delete [] r->ref;
  delete [] r->anc;
  
  delete r;
}

extern int refToInt[256];
mr::funkyPars *tglfs::readposition(int nlines){
  mr::funkyPars * myfunky = mr::allocFunkyPars();
  simplefy(myfunky,nlines);
  char buf[LENS];
  int i=0;

  while(fgets(buf,LENS,posFP)){
    myfunky->likes[i] = new double[10*nInd];
    myfunky->sites[i].chromo = strdup(strtok(buf," \t\n"));
    myfunky->sites[i].position = atoi(strtok(NULL,"\t \n"));
    myfunky->ref[i] = refToInt[strtok(NULL,"\t \n")[0]];
    myfunky->anc[i] = refToInt[strtok(NULL,"\t \n")[0]];
    i++;
    if(i==nlines)
      break;
  }
  myfunky->numSites= i;
  //  fprintf(stderr,"[%s] done reading: %lu\n",__FUNCTION__,myfunky->numSites);

  return myfunky;
}

mr::funkyPars *tglfs::fetch(int nlines){
  mr::funkyPars *fp = readposition(nlines);
  if(fp->numSites==0){
    delete [] fp->likes;
    mr::deallocFunkyPars(fp);
    return NULL;
  }
  for(int i=0;i<nInd;i++)
    tglfObj[i].fetch(fp->likes,fp->numSites);
  return fp;
}


void tglf::fetch(double **liks,int nlines){
  //  printf("glf[%d]=starting readlines size of asso: %lu\n",id,asso.size());

  for(int i=0;i<nlines;i++){
    double ten[10];
    int bytesRead = gzread(inGlf,ten ,sizeof(double)*10);
    if(bytesRead==0){
      fprintf(stderr," EOF to soon:site=%d bytesRead=%d \n",i,bytesRead);
      break;
    }else if(sizeof(double)*10!=bytesRead){
      fprintf(stderr," error reading chunk from glf:%d\n",bytesRead);
      exit(0);
    }
    //now rescale the data to the natural log
    //we are assuming its in log10 the inputfiles
    double *tmp_lk = liks[i] + 10*id;
    for(int ii=0;ii<10;ii++){
      tmp_lk[ii] = ten[bgiToSam[ii]]*log(10);
      //tmp_lk[ii] = ten[ii]*log(10);
      //    fprintf(stderr,"%f,",tmp_lk[ii]);
    }
    //now rescale to the most likely
    //continue;
    double max = tmp_lk[0];
    for(int index=1;index<10;index++)
      if(tmp_lk[index]>max)
	max = tmp_lk[index];
    for(int index=0;index<10;index++){
      tmp_lk[index] = tmp_lk[index]-max;
      //      fprintf(stderr,"%f,",tmp_lk[index]);
    }
  }
}




void tglfs::init(argStruct *arguments){
  getOptions(arguments);
  nInd =(int) filenames.size();
  posFP = aio::getFILE(posfname,"r");
  tglfObj = new tglf[nInd];
  for(int i=0;i<nInd;i++)
    tglfObj[i].init(i,filenames[i]);
}

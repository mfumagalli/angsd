#include <sys/stat.h>
#include <zlib.h>
#include <cassert>
#include "argStruct.h"
#include "beagleReader.h"
#include "analysisFunction.h"

void beagle_reader::printArg(FILE *argFile){
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
  fprintf(argFile,"\t-beagle\t%s\t(Beagle Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-intName=%d\t(Assume First column is chr_position)\n",intName);
  fprintf(argFile,"\tUse -chunkSize for defining how many sites to use at a time\n\tUse -fai to supply a fai file\n");
}
void beagle_reader::getOptions(argStruct *arguments){
  fname=NULL;
  intName=1;
  fname=angsd::getArg("-beagle",fname,arguments);
  intName=angsd::getArg("-intName",intName,arguments);
  printArg(arguments->argumentFile);
}





beagle_reader::beagle_reader(argStruct *arguments){
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);
  fflush(stderr);
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-beagle")){
      printArg(stdout);
      exit(0);
    }else if(!strcmp(arguments->argv[1],"-intName")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);


  if(!fexists(fname)) {
    fprintf(stderr,"file %s does not exits\n",fname);
    exit(0);
  }
  openBeagleFile = gzopen(fname,"rb");
  //  fprintf(stderr,"beagle file is open\n");
  //fflush(stderr);
  size_t lens = 1000000;
  char *buffer = new char[lens];
  const char *delims = "\t \n";

  int nCol=1;

  gzgets(openBeagleFile,buffer,lens); 
  strtok(buffer,delims);
  while(strtok(NULL,delims))
    nCol++;
  if(nCol % 3 ){
    fprintf(stderr,"\t-> Number of columns should be a multiple of 3, nCol=%d\n",nCol);
    exit(0);
  } 
  nInd=nCol/3-1;
  arguments->nInd=nInd;
  delete[] buffer;

  char *faiName = NULL;
  faiName=angsd::getArg("-fai",faiName,arguments);
  if(faiName==NULL){
    fprintf(stderr,"Must supply -fai for beagle input\n");
  }else
    free(faiName);
  //fprintf(stderr,"nCol=%d\tnSamples=%d\n",nCol,nInd);

}

mr::funkyPars *beagle_reader::fetch(int chunksize){

  //  fprintf(stdout,"nInd %d\tchunksize %d\n",nInd,chunksize);
  char refToChar[256] = {
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
  };
  
  size_t lens = 1000000;
  char *buffer = new char[lens];
  const char *delims = "\t \n";
  const char *delims2 = "_\t \n";
  double **post = new double*[chunksize];
  
  mr::funkyPars * myfunky = mr::allocFunkyPars();
  myfunky->sites = new loci[chunksize];
  myfunky->major = new char[chunksize];
  myfunky->minor = new char[chunksize];
  
  
  for(int s=0;s<chunksize;s++)
    post[s] = new double[nInd*3];
  
  int nSites=0;
  static int positions =0;//every site is a new position across different chunks
  
  while(gzgets(openBeagleFile,buffer,lens)){
    
    assert(buffer!=NULL);
    if(intName){
      myfunky->sites[nSites].chromo = strdup(strtok(buffer,delims2));
      myfunky->sites[nSites].position = atoi(strtok(NULL,delims2))-1;
      //  fprintf(stderr,"poisi=%d\n",myfunky->sites[nSites].position);
    }
    else{
      myfunky->sites[nSites].chromo = strdup(strtok(buffer,delims));
      myfunky->sites[nSites].position = positions++;
    }
    
    // char *SNPname=strtok(buffer,delims);


    myfunky->major[nSites] = refToChar[strtok(NULL,delims)[0]];
    myfunky->minor[nSites] = refToChar[strtok(NULL,delims)[0]];
    
    //fprintf(stdout,"chr %s\t pos %d\t major %d\tminor %d\n",myfunky->sites[nSites].chromo,myfunky->sites[nSites].position,myfunky->major[nSites],myfunky->minor[nSites]);
    for(int i=0;i<nInd*3;i++){
      char *tsk = strtok(NULL,delims);
      assert(tsk!=NULL);
      post[nSites][i] = atof(tsk);
    
    }
    
    nSites++;
    if(nSites>=chunksize)
      break;
  }

  if(nSites<chunksize){
    for(int s=nSites;s<chunksize;s++)
      delete[] post[s];
  }
  myfunky->nInd=nInd;
  myfunky->post=post;
  myfunky->numSites = nSites;
  delete[] buffer;
  if(nSites==0)
   return(NULL);
 return(myfunky);
}

int beagle_reader::fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

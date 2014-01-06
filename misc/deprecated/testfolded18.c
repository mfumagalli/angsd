/*


  21 may 2011 testfold9
  ,added whole region optim as default
  added 3 kinds of startsites 0=defualt,1=flat,2=random
  added log output of result, program also prints out the startsite
  
  24may  testfold 10
  added possiblity to pass a point used for likecalculation.
  -calcLike fname
  -doUp will uptransform if the input is in the shifted space

  26may testfold11
  fixed factor in gradient such that all startpoints gives the same
  

  
  29 may testfold12
  added threading for faster computation (should have done much before...) -threads

  30 may testfold 14 (we dont want bad luck)
  changed multiprogramming pattern to static worker pool instead

  9 oct 2011 testfold 15, change screen output, if -nThreads = 1 not so much output  

  10 nov 2011 When given a -nsites argument to high the program hanged

  5 dec 2011 int overflow in fsize printout 

  15 may 2012 changed -nThread to -nThreads and -ncat to -nChr when passing args;
*/

//tst
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "bfgs.h"
#include <sys/stat.h>
#include <pthread.h>


pthread_mutex_t mut1;
pthread_mutex_t mut2;

#define LENS 10000 //maximum bytelength per line,

int printVal = 0;
int nThread =1;
int *xOffsets = NULL;

typedef struct dMatrix_t {
  size_t x; //size_t is the largest primitive datatype.
  size_t y; //we might need that later
  double **data;
} dMatrix;


dMatrix *likedat;//the matrix of likevalues 
int dim = -1; //this is being set according to useFolded and the dimension of the file beeing read.



int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}


FILE *getFILE(const char*fname,const char *mode){
  int writeFile = 0;
  for(int i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  if(writeFile&&fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    fprintf(stderr,"\t-> Try settting prefix for output with -outnames\n");
    exit(0);
  }
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening File handle for file:%s\n",fname);
    exit(0);
  }
  return fp;
}

dMatrix *allocMatrix(size_t x,size_t y){


  double **mat=malloc(sizeof(double*)*x);

  for(int i=0;i<x;i++)
    mat[i]= malloc(sizeof(double)*y);
  
  dMatrix *ret = malloc(sizeof(dMatrix));
  ret->data = mat;
  ret->x=x;
  ret->y=y;
  return ret;
}


void normalize(double ary[],int len){
  double tsum=0;			       
  for(int i=0;i<len;i++)
    tsum += ary[i];
  
  for(int i=0;i<len;i++)
    ary[i] = ary[i]/tsum;
}

void printAry (double ary[],FILE *fp,int len){
  //  fprintf(stderr,"%s dim=%d\n",__FUNCTION__,len);
  for(int i=0;i<len;i++)
    fprintf(fp,"[%d]=%f\t",i,ary[i]);
  fprintf(fp,"\n");
  fflush(fp);
}



char *append(const char* a,const char *b){
  char *c = malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

int sit =0;
void readdat_bin(FILE *fp,int ncat,size_t nSites){
 
  
  size_t i;
  for( i=0;i<nSites;i++){
    //fprintf(stderr,"i=%d\n",i);
    if(ncat+1 != fread(likedat->data[i],sizeof(double),ncat+1,fp)){
      //      fprintf(stderr,"\t Done reading:nCols=%d from \"line\"=%d assuming eof\n",ncat+1,i+sit);
      likedat->x=i;
      break;
    }
  }
  sit += i;
}

//fold the matrix
void foldMatrix(dMatrix *mat){
  int newDim = (mat->y-1)/2+1;
  for (int i=0;i<mat->x;i++)
    for(int j=0;j<newDim-1;j++) {// we shouldn't touch the last element
      mat->data[i][j] = log(exp(mat->data[i][j]-log(2.0)) + exp(mat->data[i][mat->y-1-j]-log(2.0)));//THORFINN NEW
      //  mat->data[i][j] = log(exp(mat->data[i][j])/2.0 + exp(mat->data[i][mat->y-1-j])/2.0);//RASMUS OLD
}
  mat->y=newDim;
}



double likefn(const double locpar[],int printVal,int from,int to){
  //  fprintf(stderr,"%s using (%d,%d)\n",__FUNCTION__,from,to);
  fflush(stderr);
  int i, j;
  double sum, totsum = 0.0;
  for (i=from; i<to; i++){
    sum = 0.0;
    for (j=0; j<likedat->y; j++) //we have already taken care of the fold, so we dont need to split it up.
      sum += locpar[j]*(exp(likedat->data[i][j]));
    if(!isinf(sum))
      totsum = totsum + log(sum);
  }

  if(printVal){
    printf("Likefn: %f: ",totsum);
    for (i=0; i<likedat->y; i++)
      printf("%lf ",locpar[i]);
    printf("\n");
  }    
  return totsum;
}

void *likeDifffn(double in[],double out[],int printVal,int from,int to){
  if(printVal){
    fprintf(stderr,"%s:\n",__FUNCTION__);
    for (int i=0; i<likedat->y; i++)
      fprintf(stderr,"%lf ",in[i]);
    fprintf(stderr,"\n");
  }  
  
  
  {
    //    double pik=0;
    for (int s=from; s<to; s++){
      double bottom =0;
      for(int j=0;j<likedat->y;j++){
	bottom += in[j]*(exp(likedat->data[s][j]));
	
      }
      //	fprintf(stderr,"%f\n",bottom);
      for (int k=0; k<likedat->y; k++) //we have already taken care of the fold, so we dont need to split it up.	
	out[k] += exp(likedat->data[s][k]-log(bottom));
    }
   

    
  }
  if(printVal){
    fprintf(stderr,"likefifffnRes:\n");
    for (int i=0; i<likedat->y; i++)
      fprintf(stderr,"%.2f ",out[i]);
    fprintf(stderr,"\n");
  }  
  return NULL;
}


typedef struct like_diff_pars_t {
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
} like_diff;

size_t nSlaves =0;
void *likeDifffn_slave(void *p){
  pthread_mutex_lock(&mut1);
  nSlaves++;
  pthread_mutex_unlock(&mut1);
  like_diff *pars =(like_diff *) p;
  // fprintf(stderr,"%s[%d]\n",__FUNCTION__,pars->threadId);
  fflush(stderr);
  likeDifffn(pars->inparameters,pars->outparameters, 0,xOffsets[pars->threadId],xOffsets[pars->threadId+1]);
  //pars->likeRes = 7777;
  return NULL;
}


void likeDifffn_master(double in[],double out[]){
  // fprintf(stderr,"[%s]\t->doThread with nthread=%d\n",__FUNCTION__,nThread);
  pthread_t likeDiff_slave_threads[nThread];
  like_diff pars[nThread];
  double outs[nThread][dim];
  int rc;
  for(int i=0;i<nThread;i++){
    pars[i].threadId=i;
    pars[i].inparameters = in;
    for(int j=0;j<dim;j++)
      outs[i][j] =0;
    pars[i].outparameters = outs[i];
    //    printAry(outs[i],stderr,dim);
    rc = pthread_create(&likeDiff_slave_threads[i],NULL,likeDifffn_slave,&pars[i]);
    
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }

  for(int i=0;i<nThread;i++){
    pthread_join(likeDiff_slave_threads[i], NULL);
    //    printAry(outs[i],stderr,dim);
    for (int j=0;j<dim;j++)
      out[j] += outs[i][j];
    //       fprintf(stderr,"thread done running:%f\n",pars[i].likeRes);
    //loglike += -pars[i].likeRes;
  }
     
}
void likeDifffnShift(double in[],double out[]){
  if(printVal){
    fprintf(stderr,"down Input: dim=%d \n",dim);
    for(int i=0;i<dim-1;i++)
      fprintf(stderr,"%f\t",in[i]);
    fprintf(stderr,"\n");
    fflush(stderr);
  }

  double sum = 1.0;
  double locin[dim];
  double locout[dim];
  for(int i=0;i<dim;i++)
    locout[i]=0;

  for (int i=0; i<dim-1; i++)
    sum = sum + in[i];
  locin[0]=1.0/sum;

  for (int i=0; i<dim-1; i++)
    locin[i+1]=in[i]/sum;
  if(nThread==1)
    likeDifffn(locin,locout,0,0,likedat->x);
  else
    likeDifffn_master(locin,locout);

  if(printVal){

    fprintf(stderr,"diff likelihood: \n");
    for(int i=0;i<dim;i++)
      fprintf(stderr,"%f\t",locout[i]);
    fprintf(stderr,"\n");

  }


  double myfac=1;
  for(int i=0;i<dim-1;i++)
    myfac = myfac+in[i];
  myfac=myfac*myfac;
  //myfac=1;
  //  fprintf(stderr,"\nMYFAC=%f\n",myfac);
  for(int i=0;i<dim-1;i++){
    out[i] = locout[0]/myfac;
    for(int j=1;j<dim;j++)
      if((j-1)!=i)
	out[i] += in[j-1]*locout[j]/myfac;
      else
	out[i] += -(sum - in[j-1])*locout[j]/myfac;
  }
  if(printVal){
    fprintf(stderr,"gradient of composite: \n");
    for(int i=0;i<dim-1;i++)
      fprintf(stderr,"%f\t",out[i]);
    fprintf(stderr,"\n");
  }
}

 typedef struct neg_fn_pars_t {
   int threadId; //size_t is the largest primitive datatype.
   double *parameters;
   double likeRes; //we might need that later
 } neg_fn_pars;



size_t nNegFnSlave =0;
 void *neg_fn_slave(void *p){
   pthread_mutex_lock(&mut2);
   nNegFnSlave ++;
   pthread_mutex_unlock(&mut2);
   neg_fn_pars *pars =(neg_fn_pars *) p;
   //   fprintf(stderr,"\nin slave: threadid=%d\n",pars->threadId);
   // fflush(stderr);
   pars->likeRes = likefn(pars->parameters,0,xOffsets[pars->threadId],xOffsets[pars->threadId+1]);
   //pars->likeRes = 7777;
   return NULL;
 }

double neg_fn (const double para[]){

   double fn,  sum = 1.0;
   double locpar[dim];
   for (int i=0; i<dim-1; i++)
       sum = sum + para[i];
   locpar[0]=1.0/sum;
   for (int i=0; i<dim-1; i++)
     locpar[i+1]=para[i]/sum;

   if(nThread==1){//nothreading
     //     fprintf(stderr,"in normal\n");
     fn = likefn(locpar,0,0,likedat->x);
     return -1.0*fn;
   }else{//threading
     // fprintf(stderr,"\t->doThread with nthread=%d\n",nThread);
     pthread_t neg_fn_threads[nThread];
     neg_fn_pars pars[nThread];
     int rc;
     //     fprintf(stderr,"here\n");
     for(int i=0;i<nThread;i++){
       pars[i].threadId=i;
       pars[i].parameters = locpar;
       rc = pthread_create(&neg_fn_threads[i],NULL,neg_fn_slave,&pars[i]);

       if(rc)
	 fprintf(stderr,"error creating thread\n");

     }
     //  fprintf(stderr,"after\n");
     double loglike =0;
     for(int i=0;i<nThread;i++){
       pthread_join(neg_fn_threads[i], NULL);
       //       fprintf(stderr,"thread done running:%f\n",pars[i].likeRes);
       loglike += -pars[i].likeRes;
     }
     return loglike;
  }
}

void writematrix(FILE *out,dMatrix *mat){
  fprintf(stderr,"\t ->writing matrix with dim(%lu,%lu)\n",mat->x,mat->y);
  for (int i=0;i<mat->x;i++){
    for(int j=0;j<mat->y;j++)
      fprintf(out,"%.15f\t",mat->data[i][j]);
    fprintf(out,"\n");
  }

}

int nOptim =1;


void getRand(double ary[],int len){
  srand((unsigned)time(NULL));
  double tmp[len+1];
  for(int i=0;i<len+1;i++)
    tmp[i]=((double) rand() / (RAND_MAX)) ;
  normalize(tmp,len+1);
  for(int i=0;i<len;i++)
    ary[i]=tmp[i+1]/tmp[0];
}


void runOptim(FILE *of,FILE *of_start,double *res,int startType){
  fprintf(stderr,"\t-> Dimension of optimization:%d on nSites=%lu\n",dim-1,likedat->x);
  fflush(stderr);

  //  exit(0);
  int i;  /* number of parameters */
  
  double para[dim-1];  /* parameters */
  double min[dim-1]; /* lower bound */
  double max[dim-1]; /* upper bound */
  int nbd[dim-1]; /* boundary constraint option; nbd[i] = 0 or 1 or 2 or 3; see bfgs.h for details */
  int noisy=-1; /* bigger value produces more detailed output */ 
  
  getRand(para,dim-1);
  

  for (i=0; i<dim-1; i++){
    if(startType==0)
      para[i] = 0.01/(1.0+(double)i);
    else if(startType==1)//flat prob
      para[i] = 1;
    
    min[i] = 0.0000001;
    max[i] = 10.0;
    nbd[i] = 2;
  }
  printAry(para,of_start,dim-1);
  
  
  /* derivative function can replace 'NULL' below; see bfgs.h */
  /* need to provide negative value of the function you would like to maximize (neg_fn) */ 
  double in[dim];
  double out[dim];

  double fnmax = findmax_bfgs( dim-1, para,NULL, neg_fn, likeDifffnShift, min, max, nbd, noisy );
  //double fnmax = findmax_bfgs( dim-1, para, neg_fn, NULL, min, max, nbd, noisy );  
  printf("\t-> Likelihood function value(%d) =%f \n",nOptim++, fnmax); 
  double sum = 1.0;
  for (i=0; i<dim-1; i++)
    sum += para[i];



  double opt[dim];
  opt[0] = 1.0/sum;
  for(int i=0;i<dim-1;i++)
    opt[i+1] = para[i]/sum;
  
  
  fprintf(of,"%f\t",1.0/sum);
  for(int i=0;i<dim-1;i++)
    fprintf(of,"%f\t",para[i]/sum);
  fprintf(of,"\n");
  fflush(of);


  //update the global estimate
  res[0] += 1.0/sum;
  for(int i=0;i<dim-1;i++)
    res[i+1] += para[i]/sum;
}

void info(){
   fprintf(stderr,"-> supply -binput JOINTFILE -nChr NUMBER_OF_CATS (-doFold) (binary jointfile) -outnames PREFIX\n");
}

void print_time(){
  struct tm *local;
  time_t t;
  t=time(NULL);
  local=localtime(&t);
  printf("\t-> %s",asctime(local));
  
}
//ary length is newDim-1
double *up(double *ary,int newDim){
  double *tmp=malloc(newDim *sizeof(double*));
  double ts=1;
  for(int i=0;i<newDim-1;i++)
    ts += ary[i];
  tmp[0] = 1/ts;
  for(int i=0;i<newDim-1;i++)
    tmp[i+1] = ary[i]/ts;
  return tmp;
}


double *getTrue(const char *fname){
  fprintf(stderr,"opening: %s\n",fname);
  double *ret = malloc(LENS*sizeof(double));
  FILE *of=NULL;
  if(NULL==(of=fopen(fname,"r"))){
    printf("error opening:%s\n",fname);
    exit(-1);
  }
  char buffer[LENS];
  int index=0;
  
  while(fgets(buffer,LENS,of)!=NULL){;//loop through all lines
    //printf("buffer:%s\n",buffer);
    ret[index++]=atof(strtok(buffer," \tn"));
    char *tok=NULL;
    while(NULL!=(tok=strtok(NULL," \tn")))
      ret[index++] = atof(tok);
    //fprintf(stderr,"Len of true frea:%d\n",index);
    
  }
  //  fprintf(stderr,"Length of parfile: %d\n",index);
  fclose(of);
  return ret;
}

int main(int argc, char* argv[]){
  pthread_mutex_init(&mut1, NULL);
  pthread_mutex_init(&mut2, NULL);

  if(argc==1){
    info();
    return 0;
  }
  print_time();
  fprintf(stderr,"\t-> (build: %s,%s)\n",__DATE__,__TIME__);

  int argPos=1;//the first usable argument is the third one.

  //loop through program options
  double *true =NULL;
  char *binput=NULL;
  int ncat = 0;
  size_t nSites = 0;
  int doFold = 0;
  int startType =0;
  char *outnames = NULL;
  char *parFile = NULL;
  char *parFile2 = NULL;
  double doUp =0;
  while(argPos <argc){
    if(strcmp(argv[argPos],"-useFolded")==0)
      doFold  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-true")==0)
      true  = getTrue(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-binput")==0)
      binput  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nChr")==0)
      ncat  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nSites")==0)
      nSites  = atol(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outnames")==0)
      outnames  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-startType")==0)
      startType  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-calcLike")==0)
      parFile  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-doUp")==0)
      doUp  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-calcGrad")==0)
      parFile2  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nThread")==0){
      fprintf(stderr,"supply -nThreads instead");
      return 0;
    }else if(strcmp(argv[argPos],"-nThreads")==0)
      nThread  = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0;
    }
    argPos+=2;
  }


  //check which input

  if(binput==NULL ){
    info();
    return 0;
  }
 

  FILE *of=NULL;
  FILE *of_log=NULL;
  FILE *of_local=NULL;
  FILE *of_local_start=NULL;
  char *of_filename = NULL;
  char *of_filename_log = NULL;
  char *of_filename_local =NULL;
  char *of_filename_local_start=NULL;
  FILE *fp=NULL;
 
  if(ncat==0){
    fprintf(stderr,"\t-> Must supply -ncat argument\n");
    return 0;
  }
  //do simple sanity check to check if the size of file matches
 
  size_t filesize =fsize(binput);
  size_t sitesInFile = filesize/(sizeof(double)*(ncat+1));
  fprintf(stderr,"\t-> Filesize: %lu",filesize);
  
  if((filesize %(sizeof(double)*(ncat+1)) )){
    fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
    exit(-1);
  }
  fprintf(stderr,"\t nSites in file: %lu\n",sitesInFile);
  if(nSites==0 ||nSites>sitesInFile){
    nSites = sitesInFile;
  }

  fprintf(stderr,"\t Number of sites within each optimization to: %lu\n",nSites);
  fp = getFILE(binput,"r");

  if(outnames==NULL)
    of_filename =append(binput,".ml");
  else
    of_filename = append(outnames,".ml");
  of_filename_log = append(of_filename,".log");

  of_filename_local = append(of_filename,".local");
  of_filename_local_start = append(of_filename_local,".start");
  fprintf(stderr,"\t-> Dumping results in: %s\t%s (useFold=%d) regionOptim:%lu\n",of_filename,of_filename_local,doFold,nSites);
  fprintf(stderr,"\t-> startoptim sites in: %s log results in:%s\n",of_filename_local_start,of_filename_log);
  
  if(parFile==NULL && parFile2==NULL){
    of=getFILE(of_filename,"w");
    of_log=getFILE(of_filename_log,"w");
    of_local=getFILE(of_filename_local,"w");
    of_local_start =getFILE(of_filename_local_start,"w");
  }

  likedat = allocMatrix(nSites,ncat+1);
  if(nThread!=1){
    xOffsets = malloc(sizeof(int)*(nThread+1));
    xOffsets[0] =0;
    for(int i=1;i<nThread;i++)
      xOffsets[i] = likedat->x/nThread + xOffsets[i-1];
    xOffsets[nThread] = likedat->x;
    for(int i=0;0&i<nThread+1;i++)
      fprintf(stderr,"%d\t",xOffsets[i]);
    //    fprintf(stderr,"\n");
  }
  if(parFile!=NULL){
    fprintf(stderr,"in parfile\n");
    printVal =1 ;
    likedat->y = ncat+1;
    readdat_bin(fp,ncat,nSites);
    dim = likedat->y;
    double *myTmp = getTrue(parFile);
    if(1){
      if(doUp==0){
	//fprintf(stdout,"%s\t%f\n",__FUNCTION__,-likefn(myTmp,printVal,0,likedat->x));
	fprintf(stdout,"%s\t%f\n",__FUNCTION__,neg_fn(myTmp));
      }else
	fprintf(stdout,"%f\n",-likefn(up(myTmp,likedat->y),printVal,0,likedat->x));
      return 0;
    }
  }
  if(parFile2!=NULL){
    fprintf(stderr,"in parfile2\n");
    likedat->y = ncat+1;
    dim = likedat->y;
    readdat_bin(fp,ncat,nSites);
    double *myTmp = getTrue(parFile2);
    
    printVal =1;
    
    if(doUp==0){
      fprintf(stderr,"in par2 doUP=0\n");
      double myRes[likedat->y-1];
      for(int i=0;i<likedat->y-1;i++)
	myRes[i] = 0;
      likeDifffnShift(myTmp,myRes);//smaller space
    }else{//bigger space /prob space 2k+1
      double myRes[likedat->y];
      for(int i=0;i<likedat->y;i++)
	myRes[i] = 0;
      likeDifffn(myTmp,myRes,1,0,likedat->x);
    }
    return 0;
  }



  double *res =(double *) malloc((ncat+1)*sizeof(double));
  for(int i=0;i<(ncat+1);i++)
    res[i] = 0.0;
  while(likedat->x!=0){
    likedat->y = ncat+1;
    readdat_bin(fp,ncat,nSites);
    
    if(doFold) //if we are only interested in the folded sfs, then fold the matrix once and for all!
      foldMatrix(likedat);
    
    dim = likedat->y; //the dimension of the matrix reflects the  dimension of the optimization
    if(likedat->x==0)
      continue;
    
    //writematrix(stdout,likes);
    if(NULL!=true)
      likefn(true,1,0,likedat->y);
    
  
    runOptim(of_local,of_local_start,res,startType);
  }
  

  //global estimate
  double ss =0;
  for(int i=0;i<dim;i++)
    ss+=res[i];
  //the res in this step in the sum of all the .ml.locals, so we just normalize it
  for(int i=0;i<dim;i++){
    fprintf(of,"%f\t",res[i]/ss);
    fprintf(of_log,"%f\t",log(res[i]/ss));
  }
  fprintf(of,"\n");
  fprintf(of_log,"\n");

  //cleanup
  fclose(of_local);
  fclose(of);
  


  for (size_t i=0; i<nSites; i++)
    free(likedat->data[i]);
  free(likedat->data);
  free(likedat);
  print_time();
  if(nThread!=1)
    fprintf(stderr,"nSlaves=%lu\tnNegFnSlaves=%lu\n",nSlaves,nNegFnSlave);
  return 0;
}

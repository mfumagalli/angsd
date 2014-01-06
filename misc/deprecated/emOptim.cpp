#include <iostream>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>

#define LENS 10000

//generic container for a matrix style datastructure
template<typename T>
struct Matrix {
  size_t x;
  size_t y;
  T **matrix;
};



double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

typedef struct emPars_t{
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
  Matrix<double> *likedat;
  int from;
  int to;
  double *lik;
} emPars;


//put 2 cstrings together;
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}



//does a file exists
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}



//get the filesize of file
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


//print the time very beautifull
void print_time(FILE *fp){
  struct tm *local;
  time_t t;
  t=time(NULL);
  local=localtime(&t);
  fprintf(fp,"\t-> %s",asctime(local));
  
}

//utilty functio nfor opening a FILE
FILE *getFILE(const char*fname,const char *mode){
  int writeFile = 0;
  for(int i=0;i<(int)strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  if(0&&writeFile&&fexists(fname)){
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



//lets read the huge file generated from the recursion
Matrix<double> getHugeFile(FILE *fp, size_t ncat,size_t nSites){
  double **mat = new double*[nSites];

  for(size_t i=0;i<nSites;i++){
    mat[i] = new double[ncat+1];
    if(ncat+1 != fread(mat[i],sizeof(double),ncat+1,fp)){
      fprintf(stderr,"\t Error reading:nCols=%lu from \"line\"=%lu assuming eof\n",ncat+1,i);
      exit(0);
    }
  }
  Matrix<double> ret;
  ret.x=nSites;
  ret.y=ncat+1;
  ret.matrix = mat;
  return ret;
}


//print an array to a file 
void printAry(double *ary,int len,FILE *fp,int doLog){
  for(int i=0;i<len-1;i++)
    if(doLog)
      fprintf(fp,"%f\t",log(ary[i]));
    else
      fprintf(fp,"%f\t",(ary[i]));
  if(doLog)
    fprintf(fp,"%f\n",log(ary[len-1]));
  else
    fprintf(fp,"%f\n",(ary[len-1]));
}




void optim_em(Matrix<double> &likedat,int iter,double *start){
  double colRes[likedat.y];
  double likeRes =log(0);
  int ite;
  for(ite=0;ite<iter;ite++){
    double likes =0;
    for(size_t i=0; i<likedat.y; i++)
      colRes[i] = 0;
    for(size_t i=0;i<likedat.x;i++){
      double tmpCol[likedat.y];
      double tmpSum = 0;
      for(size_t j=0;j<likedat.y;j++){
	tmpSum += exp(likedat.matrix[i][j]) * start[j];
	tmpCol[j] = exp(likedat.matrix[i][j]) * start[j];
      }
      likes += log(tmpSum);
      for(size_t j=0;j<likedat.y;j++)
	colRes[j] += tmpCol[j]/tmpSum;
    }
    for(size_t j=0;j<likedat.y;j++)
      start[j] = colRes[j]/likedat.x;
    if((likes-likeRes)<0.01)
      break;
    likeRes = likes;
    fprintf(stderr,"likes[%d]: %f\n",ite,likeRes);
  }
  fprintf(stderr,"likes[%d]: %f\n",ite,likeRes);
  printAry(start,likedat.y,stderr,0);
}


void *doEmStep(double *in,double *out,Matrix<double> *likedat,int from, int to,double *lik){
  *lik =0;
  for(size_t i=0;i<likedat->y;i++)
    out[i] = 0;
  for(int i=from;i<to;i++){
    double tmpCol[likedat->y];
    double tmpSum = 0;
    for(size_t j=0;j<likedat->y;j++){
      tmpSum += exp(likedat->matrix[i][j]) * in[j];
      tmpCol[j] = exp(likedat->matrix[i][j]) * in[j];
    }
    *lik += log(tmpSum);
    for(size_t j=0;j<likedat->y;j++)
      out[j] += tmpCol[j]/tmpSum;
  }
  return NULL;
}



void *doEmStep_slave(void *p){
   emPars *ep = (emPars *) p;
   doEmStep(ep->inparameters,ep->outparameters,ep->likedat,ep->from,ep->to,ep->lik);
   return NULL;
 }


double factln(int n)
{
  static double a[101]; //A static array is automatically initialized to zero.
  if (n < 0) printf("Negative factorial in routine factln: %d \n", n);
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); 
  else return gammln(n+1.0);
}


double bico(int n, int k) {
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

void run_optim(Matrix<double> &likedat,double *start,int iter,int nThread,FILE *em,FILE *em_log,FILE *finRes, float tole) {
  double likeRes = log(0);
  int *xOffsets =(int *) malloc(sizeof(int)*(nThread+1));
  xOffsets[0] =0;
  for(int i=1;i<nThread;i++)
    xOffsets[i] = likedat.x/nThread + xOffsets[i-1];
  xOffsets[nThread] = likedat.x;
  for(int i=0;1&&i<nThread+1;i++)
    fprintf(stderr,"%d\t",xOffsets[i]);
  fprintf(stderr,"\n");
  

  pthread_t thd[nThread];
  emPars pars[nThread];
  double outs[nThread][likedat.y];
  double likes[likedat.y ];
  double lik[likedat.y];
  
  double *y = new double[likedat.y];
  double *tmp;
  //int ite =0;
  //  while(1){
  //   ite++;
  int ite;
  float diff=0;
  for(ite=0;ite<iter;ite++){
    for(int i=0;i<likedat.y;i++){
      y[i] =0 ;
      lik[i] =0;
    }
    double likTotal=0;
    for(int i=0;i<nThread;i++){
      pars[i].threadId=i;
      pars[i].inparameters = start;
      pars[i].outparameters = outs[i];
      pars[i].lik = lik+i;
      pars[i].from = xOffsets[i];
      pars[i].to = xOffsets[i+1];
      pars[i].likedat = &likedat;
      //    printAry(outs[i],stderr,dim);
      int rc = pthread_create(&thd[i],NULL,doEmStep_slave,&pars[i]);
      if(rc)
	fprintf(stderr,"error creating thread\n");
      
    }
    
    for(int i=0;i<nThread;i++){
      pthread_join(thd[i], NULL);
      //    printAry(outs[i],stderr,dim);
      for (int j=0;j<likedat.y;j++)
	y[j] += outs[i][j];
      //       fprintf(stderr,"thread done running:%f\n",pars[i].likeRes);
      //loglike += -pars[i].likeRes;
    }
    for (int j=0;j<likedat.y;j++)
      y[j] = y[j]/likedat.x; 
      
    for(int i=0;i<likedat.y;i++)
      likTotal += lik[i];
    
    fprintf(em,"%f\t",likTotal);
    fprintf(em_log,"%f\t",likTotal);
    printAry(start,likedat.y,em,0);
    printAry(start,likedat.y,em_log,1);
    //    fprintf(stderr,"liktotal[%d]:%f\n",ite,likTotal);
    tmp=start;
    start=y;
    y=tmp;
    //if(likeRes>stopLik)
    diff=fabs(likTotal-likeRes);
    if(diff <tole)
      break;
    likeRes = likTotal;
  }
  fprintf(stderr,"nIter run:%d diff in likes=%f\n",ite,diff);
  for(int i=0;i<likedat.y-1;i++)
    fprintf(finRes,"%f\t",start[i]);
  fprintf(finRes,"%f\n",start[likedat.y-1]);
  
  //  delete [] start; //minor leak
  //  delete [] y;//minor leak
  
  double *scalings =new double[likedat.y];
  for(int i=0;i<likedat.y;i++){
    scalings[i] = (i)*(likedat.y-i-1)/bico(likedat.y-1,2);
    //    fprintf(stderr,"[%d] %f\n",i,scalings[i]);
  }

  double nSeg=1-start[0]-start[likedat.y-1];
  double p =0;
  for(int i=0;i<likedat.y;i++)
    p+=scalings[i]*start[i];


  //  fprintf(stderr,"nSeg=%f p=%f D=%f thetah=%f fayH=%f\n",nSeg,p,D,theh,fayH);

  free(xOffsets);
}


void setStart(double *ary,char *fname,int len){
  FILE *fp = getFILE(fname,"r");
  char buf[LENS];
  fgets(buf,LENS,fp);
  ary[0] = atof(strtok(buf," \t"));
  for(int i=1;i<len;i++)
    ary[i] = atof(strtok(NULL," \t"));
}




void info(){
  fprintf(stderr,"\t-> -binput -nChr -maxIter -nThread -outnames\n");
}


int main(int argc,char **argv){
  if(argc==1){
    info();
    return 0;
  }
  print_time(stderr);
  

  int argPos=1;//the first usable argument is the third one.

  //loop through program options
  char *binput=NULL;
  int ncat = 0;
  int nThread = 4;
  int maxIter =100;
  char *outnames=NULL;
  float tole = 0.00001;
  while(argPos <argc){
    if(strcmp(argv[argPos],"-binput")==0)
      binput  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nChr")==0)
      ncat  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-maxIter")==0)
      maxIter  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-tole")==0)
      tole  = atof(argv[argPos+1]);
    

    else if(strcmp(argv[argPos],"-nThread")==0)
      nThread  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outnames")==0)
      outnames  = (argv[argPos+1]);
 
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0;
    }
    argPos+=2;
  }
  if(binput==NULL ||ncat==0){
    info();
    return 0;
  }
  if(outnames==NULL)
    outnames=binput;

  size_t filesize =fsize(binput);
  fprintf(stderr,"\t-> Filesize: %lu",filesize);
  

  if((filesize %(sizeof(double)*(ncat+1)) )){
    fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
    exit(-1);
  }
  fprintf(stderr,"\t nSites in file: %lu\n",filesize/(sizeof(double)*(ncat+1)));
  size_t nSites = filesize/(sizeof(double)*(ncat+1));
  FILE *fp = getFILE(binput,"r");
  
  Matrix<double> likedat = getHugeFile(fp,ncat,nSites);
  // likedat.x=10000;
  fprintf(stderr,"\t-> Done reading huge files\n");
  
  fprintf(stderr,"\t-> begin em on datadim = (%lu,%lu)\n",likedat.x,likedat.y);
  double start[likedat.y];
  for(int i=0;i<likedat.y;i++)
    start[i] = 1.0/likedat.y;
  
  
  char *fname_em = append(outnames,".em");
  FILE *em =getFILE(fname_em,"w");
  char *fname_em_log = append(fname_em,".log");
  FILE *em_log = getFILE(fname_em_log,"w");
  char *fname_em_ml = append(outnames,".em.ml");
  FILE *em_ml =getFILE(fname_em_ml,"w");
  run_optim(likedat,start,maxIter,nThread,em,em_log,em_ml,tole);

  for(int i=0;i<likedat.x;i++)
    delete [] likedat.matrix[i];
  delete [] likedat.matrix;

  
  print_time(stderr);
  fprintf(stderr,"dumped:%s %s %s \n",fname_em,fname_em_ml,fname_em_log);
  
  return 0;
}

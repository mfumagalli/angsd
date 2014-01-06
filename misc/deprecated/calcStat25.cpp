/*
  make constants static. This will speed up alot

 */

#include <iostream>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <vector>
#include <zlib.h>
#define LENS 10000

typedef struct{
  int pos;
  double tW;//watterson
  double tP;//pairwise
  double tFL;//FU and LI, singleton category
  double tH;// fay thetaH
  double tL;// zeng thetaL
}thetas ;


double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1, double c1) {
  double e1;
  e1 = c1/a1;
  return (e1);
}

double e2f(double a1, double a2, double c2){ 
  double e2;
  e2 = c2/((a1*a1)+a2);
  return (e2);
}


double c1f(double a1, double b1) {
  double c1;
  c1 = b1 - (1/a1);
  return (c1);
}


double c2f(int nsam, double a1, double a2, double b2) {
  double c2;
  c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
  return (c2);
}


double cn(int n){
  double top =  2*n*a1f(n)-4*(n-1);
  double bot =  (n-1)*(n-2);
  return(top/bot);
}


double vd(int n){
  double led1 = a1f(n)*a1f(n)/(a2f(n)+a1f(n)*a1f(n));
  double led2 = cn(n)-(n+1)/(1.0*(n-1));
  return(1+led1*led2);

}


double ud(int n){
  return(a1f(n)-1-vd(n));
  
}


double vf(int n){
  double top = cn(n)+2*(n*n+n+3)/(1.0*(9*n*(n-1)))-2/(1.0*(n-1));
  double bot = a1f(n)*a1f(n)+a2f(n);
  
  return(top/bot);
}

double uf(int n){
  double top = 1+(n+1)/(3.0*(n-1))-4*((n+1)/(1.0*(n-1)*(n-1)))*(a1f(n+1)-(2*n)/(1.0*(n+1)));
  double bot = a1f(n);
  return(top/(bot)-vf(n)  );
}

double H_var(int n){

  double top = 18*n*n*(3*n+2)*a2f(n+1)-(88*n*n*n+9*n*n-13*n+6);
  double bot = 9*n*(n-1)*(n-1);
  return (top/bot);
}

double fayh (int n, double thetaW, double thetaH,double thetaPi){
  double en = thetaH;
  double to = thetaPi;
  double S  =  thetaW*a1f(n);

  double led1 = (n-2)/(6.0*(n-1))*S;
  double led2 = H_var(n)*S*S;
  return ((to-en)/sqrt(led1+led2));

}


double E_var1(int n){
  return(n/(2.0*(n-1))-1/a1f(n));
}

double E_var2(int n){
  double led1 = a2f(n)/(a1f(n)*a1f(n));
  double tmp = n/(1.0*(n-1));
  double led2 = 2*tmp*tmp*a2f(n);
  double led3 = (2*(n*a2f(n)-n+1)/(1.0*(n-1)*a1f(n)));
  double led4 = (3*n+1)/(1.0*(n-1));
  return (led1+led2-led3-led4);
}
  


double zenge(int n, double thetaW, double thetaL){

  double S = thetaW*a1f(n);

  double top =thetaL-thetaW;
  double bot = E_var1(n)*S+E_var2(n)*S*S;
  return(top/sqrt(bot));
}



double tajd(int nsam, double thetaW, double sumk){
  double  a1, a2, b1, b2, c1, c2, e1, e2; 
  
  
  
  a1 = a1f(nsam);
  double segsites  = thetaW*a1;
  if( segsites == 0 ) return( 0.0) ;
  a2 = a2f(nsam);
  b1 = b1f(nsam);
  b2 = b2f(nsam);
  c1 = c1f(a1, b1);
  c2 = c2f(nsam, a1, a2, b2);
  e1 = e1f(a1, c1);
  e2 = e2f(a1, a2, c2);

  return( (sumk - (thetaW))/sqrt((e1*segsites) + ((e2*segsites)*(segsites-1))) ) ;
}


double fulid(int n, double thetaW, double thetaFL){


  double S = thetaW*a1f(n);
  double L = thetaFL*a1f(n);

  double top = S-L;
  double bot = ud(n)*S+vd(n)*S*S;
  //  fprintf(stderr,"S=%f L=%f top=%f bot=%f\n",S,L,top,bot);
  //  fflush(stderr);
  return(top/sqrt(bot));
}

double fulif(int n, double thetaW, double thetaFL,double thetaPi){

  double S = thetaW*a1f(n);
  double top = thetaPi-thetaFL;
  double bot = uf(n)*S+vf(n)*S*S;
  return (top/sqrt(bot));
    
}




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



void info(){
  fprintf(stderr,"\t-> -tinput -nChr -outnames -win -step\n");
}

thetas subtheta(const std::vector<thetas> &the,int begI,int endI){
  thetas theRes;
  theRes.tW=0;
  theRes.tP=0;
  theRes.tFL =0;  
  theRes.tH = 0;
  theRes.tL = 0; 
  for(int i=begI;i<endI;i++){
    theRes.tW += the[i].tW;
    theRes.tP += the[i].tP;
    theRes.tFL += the[i].tFL;
    theRes.tH += the[i].tH;
    theRes.tL += the[i].tL;
  }

  return theRes;
}

char *chrCur =0;


size_t readNew(gzFile fpgz,std::vector<thetas> &res){
  const char *delims = "\t\n ;";
  res.clear();
  static char *chrLast = NULL;
  static thetas t;
  if(chrLast!=NULL)//first run
    res.push_back(t);
  char buf[LENS];

  char *tok;
  while(gzgets(fpgz,buf,LENS)) {
    tok = strtok(buf,delims);
    t.pos = atoi(strtok(NULL,delims));
    t.tW = exp(atof(strtok(NULL,"\t\n ")));
    t.tP = exp(atof(strtok(NULL,"\t\n ")));
    t.tFL = exp(atof(strtok(NULL,"\t\n ")));
    t.tH = exp(atof(strtok(NULL,"\t\n ")));;
    t.tL = exp(atof(strtok(NULL,"\t\n ")));

    if(chrLast!=NULL && strcmp(chrLast,tok)==0)//true for most of the cases
      res.push_back(t);
    else if(chrLast!=NULL) {//true if we are shifting a chr
      chrCur = strdup(chrLast);
      chrLast = strdup(tok);
      break;
    }else{//only happens first time
      //      fprintf(stderr,"plugging back\n");
      chrCur = strdup(tok);
      chrLast = strdup(tok);
      res.push_back(t);
    }
  }
  if(gzeof(fpgz)==1){
    chrCur = strdup(chrLast);
    //    chrLast = strdup(tok);
  }
  if(res.size()==1)
    return 0;
  return res.size();
}

void doStat (const thetas& theSlice,FILE *estpgF,int nChr){
  double D=tajd(nChr,theSlice.tW,theSlice.tP);
  double fuliD=fulid(nChr,theSlice.tW,theSlice.tFL);
  double fuliF=fulif(nChr,theSlice.tW,theSlice.tFL,theSlice.tP);
  double fayH=fayh(nChr,theSlice.tW,theSlice.tH,theSlice.tP);
  double zengE=zenge(nChr,theSlice.tW,theSlice.tL);
  fprintf(estpgF,"%f\t%f\t%f\t%f\t%f\t%f\t%f",D,fuliD,fuliF,fayH,zengE,theSlice.tW,theSlice.tP);
}

void subs(thetas &t, thetas &t2){
  t.tW -= t2.tW;
  t.tP -= t2.tP;
  t.tFL -= t2.tFL;
  t.tH -= t2.tH;
  t.tL -= t2.tL;
}

void subs2(thetas &t, thetas &t2){
  t.tW += t2.tW;
  t.tP += t2.tP;
  t.tFL += t2.tFL;
  t.tH += t2.tH;
  t.tL += t2.tL;
}


int main(int argc,char **argv){

  if(argc==1){
    info();
    return 0;
  }
  print_time(stderr);
  

  int argPos=1;//the first usable argument is the third one.

  //loop through program options
  char *tinput=NULL;
  int nChr = 0;

  char *outnames=NULL;
  int doBlocks =0;
  int win,step;
  win=step=0;
  while(argPos <argc){
    if(strcmp(argv[argPos],"-tinput")==0)
      tinput  = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nChr")==0)
      nChr  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-win")==0)
      win  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-step")==0)
      step  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outnames")==0)
      outnames  = (argv[argPos+1]);
 
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0;
    }
    argPos+=2;
  }
  if(tinput==NULL||nChr==0){
    info();
    return 0;
  }
    
  if(win==0||step==0){
    fprintf(stderr,"No win/step size specified will do scan with window=genome (not done properly yet)\n");
  }
  if(win==step){
    fprintf(stderr,"win=step we will do blockwise estimation\n");
    doBlocks =1;
  }
  if(step>win){
    fprintf(stderr,"stepsize should be smaller than winsize\n");
    return 0;

  }
  if(outnames==NULL)
    outnames=tinput;

  gzFile fpgz = Z_NULL;
  fpgz = gzopen(tinput,"r");
  if(fpgz==Z_NULL){
    fprintf(stderr,"Problems opening file: %s\n",tinput);
  }
  char *estpgS = append(outnames,".pestPG");
  FILE *estpgF =getFILE(estpgS,"w");


  std::vector<thetas> the;
  while(readNew(fpgz,the)) {
    fprintf(stderr,"nSites in new chr: %s %lu\n",chrCur,the.size());
    fflush(stderr);
    if(win==0&&step==0){//if now step/size winsize then use entire chr
      win=the.size();
      step=win;
    }
    //superfunky window approach
    int begPos = (the[0].pos/step)*step;
    int endPos = begPos+win;
    fprintf(stderr,"begpos=%d endpos=%d\n",begPos,endPos);
    int begI = 0;
    int endI = 0;
    
    while(endI<the.size()&&the[endI].pos<=endPos)
      endI++;
    endI--;
    thetas theSlice = subtheta(the,begI,endI);
    while(endI<the.size()-1){
      //   fprintf(stderr,"start=%d stop=%d begi=%d endi=%d the.size=%lu\n",begPos,endPos,begI,endI,the.size());

      fprintf(estpgF,"(%d,%d)(%d,%d)(%d,%d)\t%s\t%d\t",begI,endI,the[begI].pos,the[endI].pos,begPos,endPos,chrCur,the[begI].pos+(the[endI].pos-the[begI].pos)/2);
      doStat(theSlice,estpgF,nChr);
      fprintf(estpgF,"\t%d\n",endI-begI);
      if(doBlocks==0) {//clever fast approach for overlapping windows
	while(1){
	  if(the[begI].pos<begPos+step){
	    subs(theSlice,the[begI]);
	    begI++;
	  }else
	    break;
	}
	while(endI<the.size()){
	  if(the[endI].pos<endPos+step){
	    subs2(theSlice,the[endI]);
	    endI++;
	  }else
	    break;
	}
	endI--;
      }else{//blokwise estimate, don't reuse old estimate
	begI=endI;//new start pos is old end pos
	while(endI<the.size()&&the[endI].pos<=endPos+step)
	  endI++;
	endI--;
	theSlice = subtheta(the,begI,endI);
      }
      begPos+=step;
      endPos+=step;
    }
    fflush(estpgF);
  }
  
  print_time(stderr);
  fprintf(stderr,"dumped:%s\n",estpgS);
  
  return 0;
}

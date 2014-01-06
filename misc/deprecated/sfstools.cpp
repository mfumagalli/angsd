/*
  thorfinn

  Jan 25 2012

  program takes output from realSFS 1, and adds a prior
  
  program can dump both binary and text (-dumpBinary INT)
  
  program can estimate tajima

  example runs: 
1)
./misc/sfstools.g++ -nChr 40 -priorFile pops1.sfs.ml -sfsFile pops1.sfs >pops1.sfstest
This adds the prior to the .sfs file, dumps the output in pop1.sfstest as binary,
2)
./misc/sfstools.g++ -nChr 40 -priorFile pops1.sfs.ml -sfsFile pops1.sfs -tajima tajima.txt -dumpBinary 0 >pops1.sfstest
This adds the prior to the .sfs file, dumps the output in pop1.sfstest as text, and calculate tajima stuff and puts this in tajima.txt

 */

#include <sys/stat.h>
#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

//struct to keep all arguments supplied
typedef struct{
  const char *priorfile;
  int nChr;
  const char *sfsfile;
  int dumpBinary;
  const char *tajima;
}args;




// Checking for file existence, using stat.h.
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

//return filesize of file
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


FILE *getFile(const char*fname,const char *mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
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



//return a filehandle
gzFile getGz(const char*fname,const char* mode){
  
  if(0&&strchr(mode,'w') && fexists(fname)){
    fprintf(stderr,"\t->File exists: %s exiting...\n",fname);
    exit(0);
  }
  gzFile fp;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}
//read program parameters
args *getArgs(char **argv){
  args *ret = new args;
  ret->sfsfile=ret->priorfile=NULL;
  ret->nChr = 0;
  ret->dumpBinary = 1;
  while(*argv){
    if(!strcmp("-nChr",*argv)) ret->nChr = atoi(*++argv);
    else if(!strcmp("-priorFile",*argv)) ret->priorfile = *++argv;
    else if(!strcmp("-sfsFile",*argv)) ret->sfsfile = *++argv;
    else if(!strcmp("-dumpBinary",*argv)) ret->dumpBinary = atoi(*++argv);
    else if(!strcmp("-tajima",*argv)) ret->tajima = *++argv;
    else {fprintf(stderr,"unknown arg \"%s\"\n",*argv);return NULL;}
    ++argv;
  }
  return ret;
}

//print program parameters
void print(const args *a,FILE *fp){
  fprintf(fp,"nChr=%d\n",a->nChr);
  fprintf(fp,"priorfile=%s\n",a->priorfile);
  fprintf(fp,"sfsfile=%s\n",a->sfsfile);
  fprintf(fp,"dumpBinary=%d\n",a->dumpBinary);
  fprintf(fp,"tajima=%s\n",a->tajima);
  if(a->priorfile==a->sfsfile || a->sfsfile==NULL || (a->nChr % 2) ){
    fprintf(stderr,"probs with pars\n");
    exit(0);
  }
}

int isGz(const char* fname){

  FILE *fp=fopen(fname,"r");
  char magic[2];
  if(2!=fread(magic,1,2,fp)) fprintf(stderr,"Problems reading from file\n");
  fclose(fp);
  return !strncmp("magic","\037\213",2);
}


void normalize_array(double *d, int len){
  double s =0;
  for(int i=0;i<len;i++)
    s+=exp(d[i]);
  //  fprintf(stderr,"s=%f\n",s);
  for(int i=0;i<len;i++)
    d[i]=exp(d[i])/s;
}


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



double factln(int n)
{
  static double a[101]; //A static array is automatically initialized to zero.
  if (n < 0) printf("Negative factorial in routine factln: %d \n", n);
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); 
  else return gammln(n+1.0);
}


//Func name: bico
//Purpose  : To calculate the nCk ( n combination k)
double bico(int n, int k) {
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}



int main(int argc, char **argv){

  if(argc==1){
    fprintf(stderr,"-nChr -priorFile -sfsFile -dumpBinary -tajima\n");
    return 0;
  }
  args *pars = NULL;
  //parse arguments
  if(!(pars=getArgs(++argv))) return 0;
  print(pars,stderr);

  //if file isn't a gz then validate the filesize.
  #if 0
  if(!isGz(pars->sfsfile)){
    size_t fs =fsize(pars->sfsfile);
    if(fs % ((pars->nChr+1)*sizeof(double))  ){
      fprintf(stderr,"problem with size of dimension truncated or wrong nChr ?\n");
      return 0;
    }else
      fprintf(stderr,"sfsfile contains: %lu\n",fs/(sizeof(double)*(pars->nChr+1)));
  }
#endif
  //readprior
  double *prior = NULL;
  if(pars->priorfile){
    prior = new double [pars->nChr+1];
    FILE *fp = getFile(pars->priorfile,"r");
    char buf[fsize(pars->priorfile)];
    if(fsize(pars->priorfile)!=fread(buf,1,fsize(pars->priorfile),fp))
      fprintf(stderr,"problems reading from priorfile\n");

    //now tokenize through the cstring
    char *tok = strtok(buf,"\t\r \n");
    int nTok =0;
    if(tok) prior[nTok] =atof(tok);
    
    for(nTok=1;nTok<pars->nChr+1;nTok++){
      tok = strtok(NULL,"\t\r \n");
      prior[nTok] = atof(tok);
    }
    fclose(fp);
  }
  for(int i=0;0&&i<pars->nChr+1;i++)
    fprintf(stderr,"prior[%d]=%f\n",i,prior[i]);
  
  //prepare sfsread
  double tmp[pars->nChr+1];
  gzFile fp = getGz(pars->sfsfile,"r");
  
  int nRead =0;
  int toRead = sizeof(double)*(pars->nChr+1);
  int nSites=0;

  //prepare tajima output
  double *scalings = NULL;
  double aConst=0;
  for(int i=1;i<(pars->nChr);i++){
    aConst += 1.0/i;
  }
  aConst = log(aConst);
  fprintf(stderr,"aConst=%f\n",aConst);
  FILE *tajimaFile = NULL;
  if(pars->tajima!=NULL){
    tajimaFile = getFile(pars->tajima,"w");
    scalings = new double[pars->nChr-1];
    for(int i=0;i<pars->nChr-1;i++){
      double ii=i+1;
      scalings[i] = (ii)*(pars->nChr-ii)/bico(pars->nChr,2);
    }
  }

  //readfile
  while((nRead=gzread(fp,tmp,sizeof(double)*(pars->nChr+1)))){
    if(nRead!=toRead)
      fprintf(stderr,"file looks truncated\n");

    if(prior){//only use prior if we have read a file
      for(int i=0;i<pars->nChr+1;i++)
	tmp[i] += log(prior[i]); //take loglike from sfs, and add the logprior
    }
    normalize_array(tmp,pars->nChr+1);
    if(pars->dumpBinary==0){
      for(int i=0;i<pars->nChr;i++)
	fprintf(stdout,"%f\t",tmp[i]);      
      fprintf(stdout,"%f\n",tmp[pars->nChr]);
    }else{
      fwrite(tmp,sizeof(double),pars->nChr+1,stdout);
    }
    
    if(tajimaFile!=NULL){
      double *workarray = tmp;
      double seq =  1-workarray[0]-workarray[pars->nChr];
      double pairwise=0;
      for(size_t i=1;i<pars->nChr;i++)//skip last element
	pairwise += workarray[i] *scalings[i-1];
      double t1 = exp(log(seq)-aConst);
      double t2 = (pairwise);
      
      
      double tajima = 2.0*(t2-t1)/(t1+t2);
      fprintf(tajimaFile,"%f\t%f\t%f\n",log(seq),log(pairwise),tajima);
    }

    nSites++;
    //return 0;
  }
  if(tajimaFile!=NULL){
    fclose(tajimaFile);
    delete [] scalings;
  }
  if(prior!=NULL)
    delete [] prior;
  fprintf(stderr,"nSites: %d parsed\n",nSites);
  delete pars;
  gzclose(fp);
  return 0;
}

/*
 *  simnextgen.h
 *  
 *
 *  Created by Rasmus Nielsen on 1/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
copied from simnextgen2.c
changed code such that ancestral is always A


moved minfreq after argparsing
10 changed .geno output such that it reflects the "real" genotypes, an doesnt depend on what I've sampled. Added a much faster uniform sampling function. The old one can still be used by -simpleRand 0
9 added inbreeding
8 added a variable errate, mean value is as supplied by user, but range=[0,2*errate]
7 added inbreeding
6 added true genotype output

 */




/*Note 0=A, 1=C, 2=G, 3=T*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>
#define PI 3.141592654


int static z_rndu=137;
double minfreq, myConst;

int dumpBinary = 1;
int simpleRand = 0;


void SetSeed(int seed)
{
  z_rndu = 170*(seed%178) + 137;
}

/*U(0,1): AS 183: Appl. Stat. 31:188-190
Wichmann BA & Hill ID.  1982.  An efficient and portable
pseudo-random number generator.  Appl. Stat. 31:188-190
x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.

Suggested to me by Z. Yang who also provided me with the source code used here. */
double uniform_rasmus()
{
  static int x_rndu=11, y_rndu=23;
  double r;

  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  if (x_rndu<0) x_rndu+=30269;
  if (y_rndu<0) y_rndu+=30307;
  if (z_rndu<0) z_rndu+=30323;
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}

double uniform_thorfinn(){
  return((double)rand() / (double)RAND_MAX);
}

double uniform(){
  double sampled;
  if(simpleRand)
    sampled = uniform_thorfinn();    
  else
    sampled = uniform_rasmus();
  //  fprintf(stdout,"%f\n",sampled);
  return sampled;
}


// Checking for file existence, using stat.h.
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}



gzFile getGz(const char*fname,const char* mode){
  if(fexists(fname)){
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

FILE *getFile(const char*fname,const char *mode){
  if(fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening File handle for file:%s\n",fname);
    exit(0);
  }
  return fp;
}


double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947,-86.50532032942, 24.01409824083,-1.231739572450, 0.1208650973866e-2,-0.5395239385e-5};

  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -=(x+0.5)*log(tmp);
  ser=1.00000000019015;
  for (j=0; j<5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


double Poisson(double xm)
{
  double gammln(double xx);
  static double sq,alxm,g,oldm=(-1.0);
  double em, t, y;
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      ++em;
      t *=uniform();
    } while (t>g);
  } 
  else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(PI*uniform());
	em=sq*y+xm;
      } while (em< 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (uniform()>t);
  }
  return em;
}



int basepick_with_errors(double errate, int inbase){
  int outbase;
  
  if (uniform()<errate){ 
    while ((outbase=(floor(4*uniform()))) == inbase);
    return outbase;
  }
  else return inbase;
}

int pick_a_base(double errate, int genotype[2])	{
  if (uniform()<0.5) 
    return basepick_with_errors(errate, genotype[0]);
  else 
    return basepick_with_errors(errate, genotype[1]);
}

void calclike(int base, double errate, double *like)	{
  /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
  double  prob;
  int i, j, k=0;
  
  for (i=0; i<4; i++){
    for (j=i; j<4; j++){
      if (base==i) 
	prob = (1.0-errate)/2.0;
      else 
	prob = errate/6.0;
      if (base==j) 
	prob = prob +(1.0-errate)/2.0;
      else 
	prob = prob + errate/6.0;
       if (prob <= 0.0) 
	like[k] = -1000000.0;
      else 
	like[k] = like[k] + log(prob);
      k++;
    }
  }
}

char int_to_base(int b)
{
  if (b==0) return 'A';
  else if (b==1) return 'C';
  else if (b==2) return 'G';
  else return 'T';
}



int print_ind_site(double errate, double meandepth, int genotype[2],gzFile resultfile,gzFile glffile){
  int i, b, numreads;
  double like[10];

  numreads = Poisson(meandepth);
  char res[numreads];

  for (i=0; i<10; i++)
    like[i] = 0.0;
  
  
  for (i=0; i<numreads; i++){
    b = pick_a_base(errate,genotype);
    res[i] = int_to_base(b);
    //    fprintf(resultfile,"%c",int_to_base(b));
    calclike(b, errate, like);
  }
  
  if(dumpBinary){
    gzwrite(resultfile,res,numreads);
    gzwrite(glffile,like,sizeof(double)*10);
  }else{
    fprintf(stderr,"textout not implemented\n");
    exit(0);
    for (i=0; i<10; i++)
      fprintf(glffile," %f",like[i]);
  }
  return numreads;
}


int pick_base_from_prior(double basefreq[4])	{
  int i = 0;
  double U, p;
  
  U = uniform();
  p = basefreq[0];
  while (U>p)
    {
      i++;
      p = p + basefreq[i];
    }
  return i;
}

/*
double simsampfreq(int n)
{
int i;
double U, p, c = 0.0;

for (i=1; i<n; i++)
	c = c + 1/(double)i;
i=1;
p=1.0/c;
U=uniform();
while (U>p)
	{
	i++;
	p=p+(1/(double)i)/c;
	}
return i;
}*/



double simfreq()
{
return exp(uniform()*myConst-myConst);
}


char *append(const char* a,const char *b){
  char *c = malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}



/* 
nsites: number of sites
nind: number of diploid individuals
errate: the sequencing error rate 
meandepth; mean sequencing depth
pvar: probability that a site is variable in the population
minfreq: minimum populaiton frequency

please notice that the chance that a site is variable in the sample is determined by both minfreq and pvar
*/

void info(){
  fprintf(stderr,"\t -> Required arg:\n\t\t-outfiles PREFIX\t PREFIX.seq PREFIX.glf PREFIX.frq PREFIX.arg\n");
  fprintf(stderr,"\t -> Optional arg:\n\t\t-nind\tNumber of diploid individuals\n");
  fprintf(stderr,"\t\t-errate\tThe sequencing error rate\n");
  fprintf(stderr,"\t\t-depth\tMean sequencing depth\n");
  fprintf(stderr,"\t\t-pvar\tProb that a site is variable in the population\n");
  fprintf(stderr,"\t\t-mfreq\tMinimum population frequency \n");
  fprintf(stderr,"\t\t-nsites\tNumber of sites\n");
  fprintf(stderr,"\t\t-F\tinbreeding coefficient\n");
  fprintf(stderr,"\t\t-model\t0=fixed errate 1=variable errate\n");
  fprintf(stderr,"\t\t-simpleRand\n");
}


int main(int argc,char *argv[]) {
  if(argc==1){
    info();
    return 0;
  }

  int i, j, k, b1, b2,  var, nsites = 500000, nind = 10;
  static int genotype[2];
  double pfreq, pvar= 0.015, meandepth = 5, errate = 0.0075;
  double basefreq[4] = {0.25, 0.25, 0.25, 0.25};
  double F=0.0;
  int model =1;
  int *freqspec =NULL;
  /*debug code*/

  static int basecheck[4];//initialize to zero
  
  //filehandles and their associated names
  gzFile glffile, resultfile;
  FILE  *freqfile, *argfile, *genofile;
  char *fGlf=NULL,*fFreq=NULL,*fSeq=NULL,*fArg=NULL,*fGeno=NULL;
  char *outfiles=NULL;

  
  int argPos=1;
  while (argPos<argc){
    if(strcmp(argv[argPos],"-nind")==0)
      nind  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-errate")==0)
      errate  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-depth")==0)
      meandepth  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-pvar")==0)
      pvar  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-mfreq")==0)
      minfreq  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outfiles")==0)
      outfiles  = (argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nsites")==0)
      nsites  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-F")==0)
      F  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-model")==0)
      model  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-simpleRand")==0)
      simpleRand  = atoi(argv[argPos+1]);
  
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0;
    }
    argPos+=2;
  }
  if(outfiles==NULL){
    fprintf(stderr,"Must supply -outfiles\n");
    return 0;
  }
  minfreq=0.0001;
  myConst = -log(minfreq);
  

  fGlf = append(outfiles,".glf.gz");
  fFreq = append(outfiles,".frq");
  fSeq = append(outfiles,".seq.gz");
  fArg = append(outfiles,".args");

  fGeno = append(outfiles,".geno");

  fprintf(stderr,"\t->Using args: -nind %d -errate %f -depth %f -pvar %f -nsites %d -F %f -model %d\n",nind,errate,meandepth,pvar,nsites,F,model);
  fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile: %s\ttruefreq: %s args:%s geno:%s\n",fSeq,fGlf,fFreq,fArg,fGeno);
  
  resultfile=getGz(fSeq,"w");
  glffile = getGz(fGlf,"w");
  freqfile= getFile(fFreq,"w"); 
  argfile=getFile(fArg,"w");
  genofile =getFile(fGeno,"w");

  freqspec = malloc((nind+1)*sizeof(int));  
  fprintf(argfile,"\t->Using args: -nind %d -errate %f -depth %f -pvar %f -nsites %d -F %f -model -%d\n",nind,errate,meandepth,pvar,nsites,F,model);
  //  fprintf(resultfile,"%i %i\n",nind, nsites);
  

  /*debug code*/
  for (i=0; i<nind+1; i++)
    freqspec[i]=0;
  /*	 
	 for (i=0; i<100000; i++){
	 pfreq = simfreq();
	 if (pfreq<0.0001) printf("ERROR\n");
	 for (j=0; j<nind; j++){
	 if (pfreq<0.0001+pow(1.5,j)*0.0001)
	 {
	 freqspec[j]++;
	 }
	 }
	 }
	 for (j=0; j<nind; j++){
	 printf("%lf %lf\n",freqspec[j]/(double)100000, (log(0.0001+pow(1.5,j)*0.0001) -log(0.0001))/(-log(0.0001)));
	 }
	 exit(-1);	
  */	
  

  for (i=0; i<nsites; i++) {
    if (uniform() > pvar) { 
      genotype[0]=genotype[1]=0;//pick_base_from_prior(basefreq);
      var=0;
      /*debug code*/
      basecheck[genotype[0]] = 2*nind;
      
    }
    else {
      
      /*debug code*/
      for (k=0; k<4; k++)
	basecheck[k]=0;
      var = 1;
      pfreq=simfreq();
      //      b1 = pick_base_from_prior(basefreq);
      b2 = 0;//pick_base_from_prior(basefreq);//changed such that reference is always A
      //      while ((b2=pick_base_from_prior(basefreq))==b1);
      while ((b1=pick_base_from_prior(basefreq))==b2);
    }
    for (j=0; j<nind; j++) {
      if (var==1) {
	if(uniform()>=F){//no inbreeding case
	  for (k=0; k<2; k++){				
	    if (uniform()<pfreq) 
	      genotype[k] = b1;
	    else genotype[k] = b2; 
	  }
	}else{//inbreding case
	  if (uniform()<pfreq) {
	    genotype[0] = b1;
	    genotype[1] = b1;
	  }
	  else {
	    genotype[0] = b2; 
	    genotype[1] = b2;
	  }
	}
	
	/*debug code*/ 
	basecheck[genotype[0]]++;
	basecheck[genotype[1]]++;
	
      }
      fprintf(genofile,"%d %d\t",genotype[0],genotype[1]);      
      int has_reads =0;
      if(model==0)
	has_reads=print_ind_site(errate, meandepth, genotype,resultfile,glffile);
      else
	has_reads=print_ind_site(2*errate*uniform(), meandepth, genotype,resultfile,glffile);

      if (j<nind-1){
	if(dumpBinary){
	  char sep[1]={'\t'};
	  gzwrite(resultfile,sep,1);
	}else{
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	  fprintf(glffile,",");
	  fprintf(resultfile,", ");
	}
	}
    }
    fprintf(genofile,"\n");
    //    printf("basecount (%d %d %d %d)\n",basecheck[0],basecheck[1],basecheck[2],basecheck[3]); 
    /*debug code*/
    if (1) {
      j=0;

      for (k=0; k<4; k++){
	if (basecheck[k]>j)
	  j=basecheck[k];
      }

      if (j<nind) {
	printf("error in freqspec calculation (%i %i %i %i)\n",basecheck[0],basecheck[1],basecheck[2],basecheck[3]); 
	exit(-1);
      }{

	freqspec[2*nind-j]++;
      }
    }
    
    if(dumpBinary){
      char sep[1] = {'\n'};
      gzwrite(resultfile,sep,1);
    }else{
      fprintf(stderr,"non binary output disabled\n");
      exit(0);
      fprintf(glffile,"\n");
      fprintf(resultfile,"\n");
    }
  }
 

  k=0;
  for (i=0; i<nind+1; i++)
    k=k+freqspec[i];
  for (i=0; i<nind+1; i++)
    fprintf(freqfile,"%f\t",(double)freqspec[i]/(double)k);	
  fprintf(freqfile,"\n");
  
  free(freqspec);
  free(fGlf);
  free(fFreq);
  free(fSeq);
  free(fArg);
  gzclose(resultfile);//fclose flushed automaticly
  gzclose(glffile);
  fclose(argfile);
  fclose(freqfile);
  fclose(genofile);
  return 0;
}




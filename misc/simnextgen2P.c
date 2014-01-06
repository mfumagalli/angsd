/*
 *  simnextgen.h
 *  
 *
 *  Created by Rasmus Nielsen on 1/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *

// Matteo:
  // 2P version: 1 or 2 populations allowed:
  // main changes:
  //	function for generating random values from Beta distribution
  //	independent drawns from Balding-Nichols distribution for estimating pop alle freq for 2 subpopulations given an ancestral pop alle freq and a FST

/ to comiple: cc -lm -lz -O simnextgen2P.c -o simnextgen2P

// Thorfinn:

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
  if(0&&fexists(fname)){
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
  if(0&&fexists(fname)){
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

// gamma distribution
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

// beta distribution using gamma distribution (Matteo)
double beta(double a, double b) {
  double x = uniform();
  return((tgamma(a+b)/(tgamma(a)*tgamma(b)))*pow(x, a-1)*pow(1-x, b-1));
}

// Poisson distribution
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

// generate sequencing errors by picking a different base
int basepick_with_errors(double errate, int inbase){
  int outbase;
  
  if (uniform()<errate){ // if error occurs
    while ((outbase=(floor(4*uniform()))) == inbase); // then take a different base
    return outbase;
  }
  else return inbase;
}

// pick a genotype including errors
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
  
  for (i=0; i<4; i++){  // all 4 possible alleles for 1st base of genotype
    for (j=i; j<4; j++){ // all 4 possible alleles for 2nd base of genotype 
      if (base==i) 
	prob = (1.0-errate)/2.0; // if right base
      else 
	prob = errate/6.0; // if wrong base
      if (base==j) 
	prob = prob +(1.0-errate)/2.0;
      else 
	prob = prob + errate/6.0;
       if (prob <= 0.0) 
	like[k] = -1000000.0;
      else 
	like[k] = like[k] + log(prob); // add all log(prob) as likelihood
      k++;
    }
  }
}

// translate from int code to letters	
char int_to_base(int b)
{
  if (b==0) return 'A';
  else if (b==1) return 'C';
  else if (b==2) return 'G';
  else return 'T';
}

// compute and print results into files
int print_ind_site(double errate, double meandepth, int genotype[2],gzFile resultfile,gzFile glffile){
  int i, b, numreads;
  double like[10];

  numreads = Poisson(meandepth);  // mumber of reads
  char res[numreads];

  for (i=0; i<10; i++) // initialize GL values for all 10 possible genotypes
    like[i] = 0.0;
  
  
  for (i=0; i<numreads; i++){
    b = pick_a_base(errate,genotype); // pick a base including errors
    res[i] = int_to_base(b); 
    //    fprintf(resultfile,"%c",int_to_base(b)); // append
    calclike(b, errate, like); // compute likelihood
  }
  
  if(dumpBinary){
    gzwrite(resultfile,res,numreads);
    gzwrite(glffile,like,sizeof(double)*10); // print likelihoods values
  }else{
    fprintf(stderr,"textout not implemented\n");
    exit(0);
    for (i=0; i<10; i++)
      fprintf(glffile," %f",like[i]);
  }
  return numreads;
}

// randomly pick a base from prior basefreq
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


// simulate sfs for 1 population (or the ancestral one)
double simfreq()
{
return exp(uniform()*myConst-myConst);
}

// simulate sfs for 2 subpopulations using Balding-Nichols distribution (Matteo)
double simfreq2P(double FST, double p_anc) {
  // FST values for the 2 subpops
  // p_anc: ancestral population allele frequency, drawn from a truncated exponential distribution
  return beta( ((1-FST)/FST)*p_anc, ((1-FST)/FST)*(1-p_anc)) ;
}

char *append(const char* a,const char *b){
  char *c = malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// MAIN

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
  fprintf(stderr,"\t\t-nind2\tNumber of diploid individuals for the eventual second population\n");
  fprintf(stderr,"\t\t-errate\tThe sequencing error rate\n");
  fprintf(stderr,"\t\t-depth\tMean sequencing depth\n");
  fprintf(stderr,"\t\t-pvar\tProb that a site is variable in the population\n");
  fprintf(stderr,"\t\t-mfreq\tMinimum population frequency \n");
  fprintf(stderr,"\t\t-nsites\tNumber of sites\n");
  fprintf(stderr,"\t\t-F\tinbreeding coefficient\n");
  fprintf(stderr,"\t\t-FST\tfst estimate for 2 populations\n");
  fprintf(stderr,"\t\t-model\t0=fixed errate 1=variable errate\n");
  fprintf(stderr,"\t\t-simpleRand\n");
}

// Matteo to Thorfinn: should we put base_freq as a parameter too? 

int main(int argc,char *argv[]) {
  if(argc==1){
    info();
    return 0;
  }

  int i, j, k, b1, b2,  var, nsites = 500000, nind = 10, nind2=0;
  static int genotype[2], genotype2[2];
  double pfreq, pfreq_pop1, pfreq_pop2, pvar= 0.015, meandepth = 5, errate = 0.0075, FST=0.1;
  double basefreq[4] = {0.25, 0.25, 0.25, 0.25};
  double F=0.0;
  int model =1;
  int *freqspec =NULL;
  int *freqspec2 =NULL;
  /*debug code*/
  
  static int basecheck[4];//initialize to zero
  static int basecheck2[4];//initialize to zero

//filehandles and their associated names
  gzFile glffile, resultfile;
  gzFile glffile2, resultfile2;
  FILE  *freqfile, *argfile, *genofile;
  FILE  *freqfile2, *genofile2;
  char *fGlf=NULL,*fFreq=NULL,*fSeq=NULL,*fArg=NULL,*fGeno=NULL;
  char *fGlf2=NULL,*fFreq2=NULL,*fSeq2=NULL,*fArg2=NULL,*fGeno2=NULL;
  char *outfiles=NULL;
  
  int argPos=1;
  while (argPos<argc){

    if(strcmp(argv[argPos],"-nind")==0) // (Matteo)
      nind  = atoi(argv[argPos+1]); // (Matteo)
    else if(strcmp(argv[argPos],"-nind2")==0) // (Matteo)
      nind2  = atoi(argv[argPos+1]); // (Matteo)
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
    else if(strcmp(argv[argPos],"-FST")==0) 
      FST = atof(argv[argPos+1]);
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


//  minfreq=0.0001; // Matteo to Thorfinn: why assigning a value to minfreq? Wasn't it an input parameter?
//  myConst = -log(minfreq);
  minfreq=0.0001;
  myConst = -log(minfreq);


  fGlf = append(outfiles,".glf.gz");
  fFreq = append(outfiles,".frq");
  fSeq = append(outfiles,".seq.gz");
  fArg = append(outfiles,".args");

  fGeno = append(outfiles,".geno");

  // eventual second pop
  fGlf2 = append(outfiles,"2.glf.gz");
  fFreq2 = append(outfiles,"2.frq");
  fSeq2 = append(outfiles,"2.seq.gz");
  fGeno2 = append(outfiles,"2.geno");
  
  fprintf(stderr,"\t->Using args: -nind %d -nind2 %d -errate %f -depth %f -pvar %f -nsites %d -F %f -FST %f -model %d\n",nind,nind2,errate,meandepth,pvar,nsites,F,FST,model); //(Matteo)

  fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile: %s\ttruefreq: %s args:%s geno:%s\n",fSeq,fGlf,fFreq,fArg,fGeno);
  fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile2: %s\ttruefreq2: %s args:%s geno2:%s\n",fSeq2,fGlf2,fFreq2,fArg,fGeno2);
  
  resultfile=getGz(fSeq,"w");
  glffile = getGz(fGlf,"w");
  freqfile= getFile(fFreq,"w"); 
  argfile=getFile(fArg,"w");
  genofile =getFile(fGeno,"w");

  resultfile2=getGz(fSeq2,"w");
  glffile2 = getGz(fGlf2,"w");
  freqfile2= getFile(fFreq2,"w"); 
  genofile2 =getFile(fGeno2,"w");
  
  freqspec = malloc((nind+1)*sizeof(int));  
  freqspec2 = malloc((nind2+1)*sizeof(int));
  
  fprintf(argfile,"\t->Using args: -nind %d -nind2 %d -errate %f -depth %f -pvar %f -nsites %d -F %f -FST %f -model -%d\n",nind,nind2,errate,meandepth,pvar,nsites,F,FST,model);
  //  fprintf(resultfile,"%i %i\n",nind, nsites);
  
  /*debug code*/
  for (i=0; i<nind+1; i++) // initialize folded sfs (from 0 to nind+1)
    freqspec[i]=0;
  for (i=0; i<nind2+1; i++) // initialize folded sfs (from 0 to nind+1)
    freqspec2[i]=0;
  
  // cycle on each site
  for (i=0; i<nsites; i++) {
    
    if (uniform() > pvar) { // if it is NOT variable... 
      
//      printf("%s\n", "Not variable");
      
      genotype[0]=genotype[1]=0; //pick_base_from_prior(basefreq);
      genotype2[0]=genotype2[1]=0; //pick_base_from_prior(basefreq);
      var=0;
      /*debug code*/
      basecheck[genotype[0]] = 2*nind; // all individuals are monorphic
      basecheck2[genotype2[0]] = 2*nind2; // all individuals are monorphic

    }
    else { // if it IS variable: a SNP!
      
      /*debug code*/
      for (k=0; k<4; k++) // initialize base checks
	basecheck[k]=0;
	basecheck2[k]=0;
      var = 1;
      //      b1 = pick_base_from_prior(basefreq);
      b2 = 0;//pick_base_from_prior(basefreq);//changed such that reference is always A
      //      while ((b2=pick_base_from_prior(basefreq))==b1);
      while ((b1=pick_base_from_prior(basefreq))==b2); // // take second allele, different from first

      // simulate population allele frequency (or ancestral if 2 subpops)
      pfreq=simfreq(); // if site is not variable, you don't need to compute pfreq

//      printf("%s\n", "Variable");
//      fprintf(stderr,"\t->PFREQ: %f\n",pfreq); //(Matteo)
	    
    } // end if it is variable

// just divide the case where una have 1 pop or 2 pops

if (nind2==0) { // only 1 pop
  
//   fprintf(stderr,"??? ONLY 1 POP\n"); 
  
// now for each individual, // assigned genotypes at each individual as couples of alleles/bases
for (j=0; j<nind; j++) {
      if (var==1) {
	if(uniform()>=F){ //no inbreeding case
	  for (k=0; k<2; k++){				
	    if (uniform()<pfreq) 
	      genotype[k] = b1;
	    else genotype[k] = b2; 
	  }
	}else{//inbreeding case
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
    } // end for j in nind
} else { // end if 1 pop, start if 2 pops

// fprintf(stderr,"!!! 2 POPULATIONS\n"); 
  
  // not considering Inbreeding within subpops: think about it?
  // then assign subpop freq based on time-split? theory of pop genet
  
  // use now balding-nichols and write one set of files for each population
 // simulate distinct pops allele freq from Balding-Nickols distribution, given an ancestral population allele frequency and an FST value
  // 2 independent draws
  pfreq_pop1=beta(FST, pfreq);
  pfreq_pop2=beta(FST, pfreq);

//   fprintf(stderr,"\t->PFREQ: %f\n",pfreq); //(Matteo)
//   fprintf(stderr,"\t->P1: %f\n",pfreq_pop1); //(Matteo)
//   fprintf(stderr,"\t->P2: %f\n",pfreq_pop2); //(Matteo)
  
  // POP 1
  for (j=0; j<nind; j++) {

//    printf("%s\n","FIRST"); 
    
    if (var==1) { // if polymorphic
	  for (k=0; k<2; k++){ // for the 2 'positions' of genotype		
	    if (uniform()<pfreq_pop1) 
	      genotype[k] = b1;
	    else genotype[k] = b2; 
	  } // end for in k
	
	/*debug code*/ 
	basecheck[genotype[0]]++;
	basecheck[genotype[1]]++;
	
      } // end if polymorphic

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
  } // end for j in nind

// POP 2
  for (j=0; j<nind2; j++) {
    
//    printf("%s\n","SECOND"); 

    if (var==1) { // if polymorphic
	  for (k=0; k<2; k++){ // for the 2 'positions' of genotype		
	    if (uniform()<pfreq_pop2) 
	      genotype2[k] = b1;
	    else genotype2[k] = b2; 
	  } // end for in k
	
	/*debug code*/ 
	basecheck2[genotype2[0]]++;
	basecheck2[genotype2[1]]++;
	
      } // end if polymorphic

      fprintf(genofile2,"%d %d\t",genotype2[0],genotype2[1]);      
      int has_reads =0;
      if(model==0)
	has_reads=print_ind_site(errate, meandepth, genotype2,resultfile2,glffile2);
      else
	has_reads=print_ind_site(2*errate*uniform(), meandepth, genotype2,resultfile2,glffile2);

      if (j<nind2-1){
	if(dumpBinary){
	  char sep[1]={'\t'};
	  gzwrite(resultfile2,sep,1);
	}else{
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	  fprintf(glffile2,",");
	  fprintf(resultfile2,", ");
	}
	}
	
  } // end for j in nind2
  
} // end if 2 pops
    
    fprintf(genofile,"\n");
    fprintf(genofile2,"\n");
    //    printf("basecount (%d %d %d %d)\n",basecheck[0],basecheck[1],basecheck[2],basecheck[3]); 
/*    printf("basecount (%d %d %d %d)\n",basecheck[0],basecheck[1],basecheck[2],basecheck[3]); 
    printf("basecount (%d %d %d %d)\n",basecheck2[0],basecheck2[1],basecheck2[2],basecheck2[3]); */
    
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

   /*debug code*/
    if (1 & nind2>0) {
      j=0;

      for (k=0; k<4; k++){
	if (basecheck2[k]>j)
	  j=basecheck2[k];
      }

      if (j<nind2) {
	printf("error in freqspec calculation (%i %i %i %i)\n",basecheck2[0],basecheck2[1],basecheck2[2],basecheck2[3]); 
	exit(-1);
      }{
	fprintf(stderr,"nind=%d\t2*nind=%d\tj=%d:\t 2*nind2-j=%d\n",nind,2*nind,j,2*nind2-j);
	freqspec2[2*nind2-j]++;
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

    if (nind2>0) {
    if(dumpBinary){
      char sep[1] = {'\n'};
      gzwrite(resultfile2,sep,1);
    }else{
      fprintf(stderr,"non binary output disabled\n");
      exit(0);
      fprintf(glffile2,"\n");
      fprintf(resultfile2,"\n");
    }
    }
  
  
  } // end cycle on each site
 

  k=0;
  for (i=0; i<nind+1; i++)
    k=k+freqspec[i];
  for (i=0; i<nind+1; i++)
    fprintf(freqfile,"%f\t",(double)freqspec[i]/(double)k);	
    fprintf(freqfile,"\n");

  k=0;
  for (i=0; i<nind2+1; i++)
    k=k+freqspec2[i];
  for (i=0; i<nind2+1; i++)
    fprintf(freqfile2,"%f\t",(double)freqspec2[i]/(double)k);	
    fprintf(freqfile2,"\n");

  free(freqspec);
  free(freqspec2);

  free(fGlf);
  free(fFreq);
  free(fSeq);

  free(fGlf2);
  free(fFreq2);
  free(fSeq2);
  free(fArg);

  gzclose(resultfile);//fclose flushed automaticly
  gzclose(glffile);
  fclose(argfile);
  fclose(freqfile);
  fclose(genofile);

  gzclose(resultfile2);//fclose flushed automaticly
  gzclose(glffile2);
  fclose(freqfile2);
  fclose(genofile2);

  return 0;
}





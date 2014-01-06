//must set ntrheads to 1
#include <cstdio>

#include "shared.h"
#include "analysisFunction.h"

#include "general.h"
#define MINLIKE -1000.0 //this is for setting genotypelikelhoods to missing (EXPLAINED BELOW)


typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  double **pLikes;
}realRes;



/*
  From 0.501 the realSFS method can take into account individual inbreeding coefficients
  Filipe guera
*/
namespace filipe{
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,char *keepSites,realRes *r,int noTrans,int doRealSFS,char *major,char *minor,double *freq,double *indF);
}


class realSFS : public general{
  int doRealSFS;
  FILE *outfileSFS;
  FILE *outfileSFSPOS;
  gzFile fpgz;
  int underFlowProtect;
  int fold;
  int isSim;
  int noTrans;
  char *anc;
  char *pest;
  double *prior; //<- outputfile form pest;
  int doThetas;
  void fin(funkyPars *p,int index,double *prior,gzFile fpgz);

  double aConst;
  double aConst2;
  double aConst3;
  double *scalings;


  double *filipeIndF;

public:
  //none optional stuff
  FILE *outfile;
  realSFS(const char *outfiles,argStruct *arguments,int inputtype);
  ~realSFS();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

  
};

void realSFS::printArg(FILE *argFile){
  fprintf(argFile,"%s:\n\n",__FILE__);
  fprintf(argFile,"-realSFS\t%d\n",doRealSFS);
  fprintf(argFile,"-underFlowProtect\t%d\n",underFlowProtect);
  fprintf(argFile,"-fold\t%d (deprecated)\n",fold);
  fprintf(argFile,"-anc=%s\n",anc);
  fprintf(argFile,"-pest=%s\n",pest);
  fprintf(argFile,"-noTrans=%d\n",noTrans);
  fprintf(argFile,"\n");

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


//log version of binomial coef
double lbico(int n, int k){
  //  return log(bico(n,k));
  return factln(n)-factln(k)-factln(n-k);
}

double bico(int n, int k){
   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

double lbico(double n, double k){
  return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);
}




void realSFS::getOptions(argStruct *arguments){
  underFlowProtect = 0;
  fold =0;
  isSim =0;
  //from command line
  anc=NULL;
  pest=NULL;
  noTrans = 0;
  prior = NULL;
  doRealSFS=angsd::getArg("-realSFS",doRealSFS,arguments);
  noTrans = angsd::getArg("-noTrans",noTrans,arguments);
  pest = angsd::getArg("-pest",pest,arguments);
  int GL = 0;
  GL = angsd::getArg("-GL",GL,arguments);
  doThetas= angsd::getArg("-doThetas",doThetas,arguments);
  if(doRealSFS==0&&doThetas==0)
    return;

  underFlowProtect=angsd::getArg("-underFlowProtect",underFlowProtect,arguments);
  fold=angsd::getArg("-fold",fold,arguments);
  char *sim1 = NULL;
  sim1=angsd::getArg("-sim1",sim1,arguments);
  if(sim1!=NULL)
    isSim =1;
  //  fprintf(stderr,"sim1=%p %s\n",sim1,sim1);
  if(doRealSFS==-999){
    doRealSFS=0;
    printArg(stdout);
    exit(0);
  }
  if(doRealSFS==0)
    return;
  anc = angsd::getArg("-anc",anc,arguments);
  if(doRealSFS && (anc==NULL&&isSim==0)){
    fprintf(stderr,"Must supply -anc for polarizing the spectrum\n");
    exit(0);
  }
  if(pest!=NULL)
    prior=angsd::readDouble(pest,arguments->nInd*2+1);

  if(GL==0 &&arguments->inputtype==7){//DRAGON this is not enough, but this is most likely what everyone done...
    fprintf(stderr,"Must supply genotype likelihoods (-GL [INT])\n");
    printArg(arguments->argumentFile);
    exit(0);
  }
  if(doRealSFS==2){
    fprintf(stderr,"\t->(Using Filipe modification of: %s)\n",__FILE__);
    int doMajorMinor =0;
    doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
    int doMaf =0;
    doMaf = angsd::getArg("-doMaf",doMaf,arguments);
    if(doMajorMinor==0||doMaf==0){
      fprintf(stderr,"Must have major/minor and MAF for using the inbreeding version of realSFS\n");
      exit(0);
    }
    char *indF_name = NULL;
    indF_name =  angsd::getArg("-indF",indF_name,arguments);
    if(indF_name==NULL){
      filipeIndF = new double[arguments->nInd];
      memset(filipeIndF, 0, arguments->nInd*sizeof(double));
      //fprintf(stderr,"Need to set -indF\n");
      //exit(0);
    }else
      filipeIndF = angsd::readDouble(indF_name,arguments->nInd);
  }


  printArg(arguments->argumentFile);
}

realSFS::realSFS(const char *outfiles,argStruct *arguments,int inputtype){

  const char *SFS = ".sfs";
  const char *SFSPOS =".sfs.pos";
  const char *THETAS =".thetas.gz";

  //default
  doRealSFS=0;
  doThetas = 0;
  outfileSFS = NULL;
  outfileSFSPOS = NULL;
  fpgz = Z_NULL;
  getOptions(arguments);
  
  if(doRealSFS==0)
    return;
  if(doRealSFS!=0&&doThetas==0){
    outfileSFS =  openFile(outfiles,SFS);
    outfileSFSPOS =  openFile(outfiles,SFSPOS);
  }else{
    fpgz = openFileGz(outfiles,THETAS,"w6h");;

    aConst=0;
    int nChr = 2*arguments->nInd;
    for(int i=1;i<nChr;i++)
    aConst += 1.0/i;
    aConst = log(aConst);//this is a1

  
   aConst2 = log((nChr*(nChr-1))/2.0);//choose(nChr,2)
   aConst3 = log((1.0*nChr-1.0));

   scalings = new double [nChr+1];
  for(int i=0;i<nChr+1;i++)
    scalings[i] = log(i)+log(nChr-i);



  }
}


realSFS::~realSFS(){
  if(doRealSFS==0)
    return;

  if(doRealSFS!=0&&doThetas==0){
    fclose(outfileSFS);
    fclose(outfileSFSPOS);
  }else if(doRealSFS!=0&&doThetas==1)
    gzclose(fpgz);
}


int whichMax(double *d,int len){
  int r=0;
  for(int i=1;i<len;i++)
    if(d[i]>d[r])
      r=i;
  //now check if site doesnt have data.
  
  if(r==0){//only check if nothing is higher than the first
    for(int i=1;i<len;i++)
      if(d[i]!=d[0])//we see a diffrence so we have information
	return r;
    return -1;//we didnt have information 
  }else
    return r;
}


int isSame(double a,double b,double tolerance){
  return (fabs(a-b)<tolerance);
}


void filipe::algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,char *keepSites,realRes *r,int noTrans,int doRealSFS,char *major,char *minor,double *freq,double *indF) {
  //  fprintf(stderr,"liks=%p anc=%p nsites=%d nInd=%d underflowprotect=%d fold=%d keepSites=%p r=%p\n",liks,anc,nsites,numInds,underFlowProtect,fold,keepSites,r);
  int myCounter =0;
  fflush(stderr);
  if(anc==NULL||liks==NULL){
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }
  double sumMinors[2*numInds+1];  //the sum of the 3 different minors
  double m[numInds][3];  //HWE priors

  for(int it=0; it<nsites; it++) {//loop over sites
    int major_offset = anc[it];
    if(major_offset==4 || keepSites[it]==0){//skip if no ancestral information
      //      r->oklist is zero no need to update
      continue;
    }
    //set the resultarray to zeros
    for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
      sumMinors[sm] = 0;
    
    //loop through the 3 different minors
    for(int minor_offset=0;minor_offset<4;minor_offset++) {
      if(minor_offset == major_offset)
	continue;

      if(noTrans){
	if(major_offset==2&&minor_offset==0||major_offset==0&&minor_offset==2)
	  continue;
	if(major_offset==1&&minor_offset==3||major_offset==3&&minor_offset==1)
	  continue;
      }

      if(doRealSFS == 2 && minor_offset != major[it] && minor_offset != minor[it]) continue;
      
      double totmax = 0.0;
      //hook for only calculating one minor
      int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
      int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
      int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
      //      fprintf(stderr,"%d:%d\t%d\t%d\n",major_offset,Aa_offset,AA_offset,aa_offset);
      //part two
      double hj[2*numInds+1];
      for(int index=0;index<(2*numInds+1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,PAa,Paa;

      // Assign freqs to ancestral and non-ancestral alleles
      double anc_freq, non_anc_freq;
      if(doRealSFS == 2){
	anc_freq = non_anc_freq = freq[it];
	if(major_offset == major[it])
	  anc_freq = 1-freq[it];
	if(minor_offset == major[it])
	  non_anc_freq = 1-freq[it];
      }

      for(int i=0 ; i<numInds ;i++) {
	if(doRealSFS == 2) { //calculate priors based on HWE + indF
	  m[i][0] = pow(anc_freq,2) + anc_freq*non_anc_freq*indF[i];
	  m[i][1] = 2*anc_freq*non_anc_freq - 2*anc_freq*non_anc_freq*indF[i];
	  m[i][2] = pow(non_anc_freq,2) + anc_freq*non_anc_freq*indF[i];
	} else { // If not realSFS2, assume uniform prior
	  m[i][0] = m[i][1] = m[i][2] = 1;
	}

	//	printf("pre scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	double GAA,GAa,Gaa;
#ifdef RESCALE
	if(0){//The rescaling is now done in 'getMajorMinor()'
	  double max = liks[it][i*10+0];
	  for(int index=1;index<10;index++)
	    if(liks[it][i*10+index]>max)
	      max = liks[it][i*10+index];
	  for(int index=0;index<10;index++)
	    liks[it][i*10+index] = liks[it][i*10+index]-max;
	}
#endif

	GAA = ((m[i][2] != 0) ? (log(m[i][2])+liks[it][i*10+AA_offset]) : (-999999999));
	GAa = ((m[i][1] != 0) ? (log(m[i][1])+liks[it][i*10+Aa_offset]) : (-999999999));
	Gaa = ((m[i][0] != 0) ? (log(m[i][0])+liks[it][i*10+aa_offset]) : (-999999999));
	if(doRealSFS==1) GAa += log(2);
	//printf("post scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	//printf("[GAA] GAA=%f\tGAa=%f\tGaa=%f\n",GAA,GAa,Gaa);

	//do underlfow protection (we are in logspace here) (rasmus style)
	if(1){
	  double mymax;
	  if (Gaa > GAa && Gaa > GAA) mymax = Gaa;
	  else if (GAa > GAA) mymax = GAa;
	  else mymax = GAA;
	  // fprintf(stdout,"mymax[%d]=%f\t",i,mymax);
	  
	  if(mymax<MINLIKE){
	    //	    fprintf(stderr,"\n%f %f %f\n",GAA, GAa,Gaa);
	    Gaa = 0;
	    GAa = 0;
	    GAA = 0;
	    totmax = totmax + mymax;
	  }else{
	    Gaa=Gaa-mymax;
	    GAa=GAa-mymax;
	    GAA=GAA-mymax;
	    totmax = totmax + mymax;
	  }
	//	fprintf(stderr,"totmax=%f\n",totmax);
	//END underlfow protection (we are in logspace here) (rasmus style)
	}

	if(underFlowProtect==0){
	  PAA=exp(GAA);
	  PAa=exp(GAa);
	  Paa=exp(Gaa);
	}else{
	  PAA =(GAA);///(MAA+MAa+Maa);
	  PAa =(GAa);///(MAA+MAa+Maa);
	  Paa =(Gaa);///(MAA+MAa+Maa);
	}

	//check for underflow error, this should only occur once in a blue moon
	if(isnan(Paa)||isnan(PAa)||isnan(Paa)){
	  //fprintf(stderr,"Possible underflow at: \t%d\t%d\n",it->first.chromo,it->first.position);
	  fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
	}
	//	fprintf(stdout,"it=%d PAA=%f\tPAa=%f\tPaa=%f\n",it,PAA,PAa,Paa);
	if(i==0){
	  hj[0] =Paa;
	  hj[1] =PAa;
	  hj[2] =PAA;
	}else{
	  //fprintf(stderr,"asdf\n");
	  for(int j=2*(i+1); j>1;j--){
	    //  print_array(stdout,hj,2*numInds+1,0);
	    //print_array(hj,2*numInds+1);
	    double tmp;
	    if(underFlowProtect==1)
	      tmp = angsd::addProtect3(log(m[i][2])+PAA+hj[j-2],log(m[i][1])+PAa+hj[j-1],log(m[i][0])+Paa+hj[j]);
	    else
	      tmp = m[i][2]*PAA*hj[j-2] + m[i][1]*PAa*hj[j-1] + m[i][0]*Paa*hj[j];
	    
	    if(isnan(tmp)){
	      fprintf(stderr,"is nan:%d\n",j );
	      
	      hj[j] = 0;
	      break;
	    }else
	      hj[j] = tmp;
	  }
	  if(underFlowProtect==1){
	    hj[1] = angsd::addProtect2(log(m[i][0])+Paa+hj[1],log(m[i][1])+PAa+hj[0]);
	    hj[0] = log(m[i][0]) + Paa + hj[0];
	  }
	  else{
	    hj[1] = m[i][0]*Paa*hj[1] + m[i][1]*PAa*hj[0];
	    hj[0] = m[i][0]*Paa*hj[0];
	  }
	}
	//ifunderflowprotect then hj is in logspace
	
      }
      //      fprintf(stderr,"%f %f %f\n",hj[0],hj[1],hj[2]);
      //fprintf(stdout,"\nscaledLikes=");
      for(int ii=0;0&&ii<10*numInds;ii++)
	fprintf(stdout,"%f\t",liks[it][ii]);

      //      totmax=0;
      for(int i=0;i<(2*numInds+1);i++)
	if(underFlowProtect==0){
	  if(doRealSFS == 2)
	    sumMinors[i] += hj[i];
	  //sumMinors[i] += hj[i]/(1-M0-M2); //As in the PLoS ONE paper
	  else
	    sumMinors[i] += exp(log(hj[i])-lbico(2*numInds,i)+totmax);
	}else{
	  if(doRealSFS == 2)
	    sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]));
	  //sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-log(1-M0-M2))); //As in the PLoS ONE paper
	  else
	    sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-lbico(2*numInds,i)+totmax));
	}
    }

    //sumMinors is in normal space, not log
    /*
      we do 3 things.
      1. log scaling everyting
      2. rescaling to the most likely in order to avoid underflows in the optimization
      3. we might do a fold also.
     */    

    if(fold) {
      int newDim = numInds+1;
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	sumMinors[i] = log(sumMinors[i] + sumMinors[2*numInds-i]);//THORFINN NEW
      sumMinors[newDim-1] = log(sumMinors[newDim-1])+log(2.0);
      angsd::logrescale(sumMinors,newDim);
      if(isnan(sumMinors[0]))
	r->oklist[it] = 2;
      //fprintf(stderr,"TAYLOR THIS IS AN ERROR: (%s,%d)\n",locs[it].chromo,locs[it].posit`ion+1);
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(numInds+1));
	myCounter++;
      }
    }else{
      for(int i=0;i<2*numInds+1;i++)
	sumMinors[i] = log(sumMinors[i]);
      angsd::logrescale(sumMinors,2*numInds+1);
      if(isnan(sumMinors[0]))
	r->oklist[it] = 2;
	//fprintf(stderr,"TAYLOR THIS IS AN ERROR: (%s,%d)\n",locs[it].chromo,locs[it].position+1);
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[2*numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
	myCounter++;
      }
    }
    
    //    fprintf(testFP,"\n");
  }
  
}



void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,char *keepSites,realRes *r,int noTrans) {
  //  fprintf(stderr,"liks=%p anc=%p nsites=%d nInd=%d underflowprotect=%d fold=%d keepSites=%p r=%p\n",liks,anc,nsites,numInds,underFlowProtect,fold,keepSites,r);
  int myCounter =0;
  if(anc==NULL||liks==NULL){
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }
  double sumMinors[2*numInds+1];  //the sum of the 3 different minors

  for(int it=0; it<nsites; it++) {//loop over sites
    int major_offset = anc[it];
    if(major_offset==4||(keepSites[it]==0)){//skip of no ancestral information
      //      r->oklist is zero no need to update
      continue;
    }
    //set the resultarray to zeros
    for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
      sumMinors[sm] = 0;
    
    //loop through the 3 different minors
    for(int minor_offset=0;minor_offset<4;minor_offset++) {
      if(minor_offset == major_offset)
	continue;
      if(noTrans){
	if(major_offset==2&&minor_offset==0||major_offset==0&&minor_offset==2)
	  continue;
	if(major_offset==1&&minor_offset==3||major_offset==3&&minor_offset==1)
	  continue;
      }
      double totmax = 0.0;
      //hook for only calculating one minor
      int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
      int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
      int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
      //      fprintf(stderr,"%d:%d\t%d\t%d\n",major_offset,Aa_offset,AA_offset,aa_offset);
      //part two
      double hj[2*numInds+1];
      for(int index=0;index<(2*numInds+1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,PAa,Paa;

      for(int i=0 ; i<numInds ;i++) {
	//	printf("pre scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	double GAA,GAa,Gaa;
#ifdef RESCALE
	if(0){//The rescaling is now done in 'getMajorMinor()'
	double max = liks[it][i*10+0];
	for(int index=1;index<10;index++)
	  if(liks[it][i*10+index]>max)
	    max = liks[it][i*10+index];
	for(int index=0;index<10;index++)
	  liks[it][i*10+index] = liks[it][i*10+index]-max;
	}
#endif
	//printf("post scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	GAA = liks[it][i*10+AA_offset];

	GAa = log(2.0)+liks[it][i*10+Aa_offset];
	Gaa = liks[it][i*10+aa_offset];
	//printf("[GAA] GAA=%f\tGAa=%f\tGaa=%f\n",GAA,GAa,Gaa);
	//do underlfow protection (we are in logspace here) (rasmus style)
	if(1){
	  double mymax;
	  if (Gaa > GAa && Gaa > GAA) mymax = Gaa;
	  else if (GAa > GAA) mymax = GAa;
	  else mymax = GAA;
	  // fprintf(stdout,"mymax[%d]=%f\t",i,mymax);
	  
	  if(mymax<MINLIKE){
	    //	    fprintf(stderr,"\n%f %f %f\n",GAA, GAa,Gaa);
	    Gaa = 0;
	    GAa = 0;
	    GAA = 0;
	    totmax = totmax + mymax;
	  }else{
	    Gaa=Gaa-mymax;
	    GAa=GAa-mymax;
	    GAA=GAA-mymax;
	    totmax = totmax + mymax;
	  }
	//	fprintf(stderr,"totmax=%f\n",totmax);
	//END underlfow protection (we are in logspace here) (rasmus style)
	}
	if(underFlowProtect==0){
	  PAA=exp(GAA);
	  PAa=exp(GAa);
	  Paa=exp(Gaa);
	}else{
	  PAA =(GAA);///(MAA+MAa+Maa);
	  PAa =(GAa);///(MAA+MAa+Maa);
	  Paa =(Gaa);///(MAA+MAa+Maa);
	}

	//check for underflow error, this should only occur once in a blue moon
	if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	  fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
	}
	//	fprintf(stdout,"it=%d PAA=%f\tPAa=%f\tPaa=%f\n",it,PAA,PAa,Paa);
	if(i==0){
	  hj[0] =Paa;
	  hj[1] =PAa;
	  hj[2] =PAA;
	}else{
	  //fprintf(stderr,"asdf\n");
	  for(int j=2*(i+1); j>1;j--){
	    //  print_array(stdout,hj,2*numInds+1,0);
	    //print_array(hj,2*numInds+1);
	    double tmp;
	    if(underFlowProtect==1)
	      tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	    else
	      tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	    
	    if(std::isnan(tmp)){
	      fprintf(stderr,"is nan:%d\n",j );
	      
	      hj[j] = 0;
	      break;
	    }else
	      hj[j]  =tmp;
	  }
	  if(underFlowProtect==1){
	    hj[1] = angsd::addProtect2(Paa+hj[1],PAa+hj[0]);
	    hj[0] = Paa+hj[0];
	  }
	  else{
	    hj[1] = Paa*hj[1] + PAa*hj[0];
	    hj[0] = Paa*hj[0];
	  }
	}
	//ifunderflowprotect then hj is in logspace
	
      }
      //      fprintf(stderr,"%f %f %f\n",hj[0],hj[1],hj[2]);
      //fprintf(stdout,"\nscaledLikes=");
      for(int ii=0;0&&ii<10*numInds;ii++)
	fprintf(stdout,"%f\t",liks[it][ii]);

      //      totmax=0;
      for(int i=0;i<(2*numInds+1);i++)
	if(underFlowProtect==0)
	  sumMinors[i] +=  exp(log(hj[i])-lbico(2*numInds,i)+totmax);
	else
	  sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-lbico(2*numInds,i)+totmax));
    }
    //sumMinors is in normal space, not log
    /*
      we do 3 things.
      1. log scaling everyting
      2. rescaling to the most likely in order to avoid underflows in the optimization
      3. we might do a fold also.
      
     */    

    if(fold) {
      fprintf(stderr,"fold is deprecated: will exit\n");
      exit(0);

#ifdef strip
      int newDim = numInds+1;
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	  sumMinors[i] = log(sumMinors[i] + sumMinors[2*numInds-i]);//THORFINN NEW
	sumMinors[newDim-1] = log(sumMinors[newDim-1])+log(2.0);
	angsd::logrescale(sumMinors,newDim);
	fwrite(sumMinors,sizeof(double),newDim,sfsfile);
#endif
    }else{
      for(int i=0;i<2*numInds+1;i++)
	sumMinors[i] = log(sumMinors[i]);
      angsd::logrescale(sumMinors,2*numInds+1);
      if(std::isnan(sumMinors[0]))
	r->oklist[it] = 2;
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[2*numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
	myCounter++;
      }
    }
    
    //    fprintf(testFP,"\n");
  }
  
}


void normalize_array(double *d, int len){
  double s =0;
  for(int i=0;i<len;i++)
    s+=d[i];

  for(int i=0;i<len;i++)
    d[i]=d[i]/s;
}


void normalize_array2(double *d, int len){
  double s =0;
  for(int i=0;i<len;i++)
    s+=exp(d[i]);
  s=log(s);

  for(int i=0;i<len;i++)
    d[i]=d[i]-s;
}



//Basicly the same as algoBayAll but only looping through the 3 derived given that an ancestral exists;




void print_array(FILE *fp,double *ary,int len,int doLogTransform){
  if(doLogTransform){
    for (int i=0;i<len-1;i++)
      fprintf(fp,"%f\t",log(ary[i]));
    fprintf(fp,"%f\n",log(ary[len-1]));
  }else{
    for (int i=0;i<len-1;i++)
      fprintf(fp,"%f\t",(ary[i]));
    fprintf(fp,"%f\n",(ary[len-1]));
  }
}



void realSFS::run(funkyPars  *p){
  if(doRealSFS==0||p->numSites==0)
    return;
  else if(doRealSFS==1||doRealSFS==2){
    realRes *r = new realRes;
    r->oklist=new char[p->numSites];
    memset(r->oklist,0,p->numSites);
    r->pLikes=new double*[p->numSites];
    if(doRealSFS==1)
      algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans);
    else
      filipe::algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans,doRealSFS,p->major,p->minor,p->results->asso->freq,filipeIndF);
    p->extras[index] = r;
  }else{
    fprintf(stderr,"unsupported realSFS arg=%d\n",doRealSFS);
  }

}

void realSFS::clean(funkyPars *p){
 if(doRealSFS==0||p->numSites==0)
    return;
 
  realRes *r=(realRes *) p->extras[index];
  
  //  realRes *r=(realRes *) p->extras[index];
  int id=0;
  for(int i=0;i<p->numSites;i++)
    if(r->oklist[i]==1)
      delete [] r->pLikes[id++];
  delete [] r->pLikes;
  delete [] r->oklist;
  delete r;

}

void printFull(funkyPars *p,int index,FILE *outfileSFS,FILE *outfileSFSPOS,char *chr,int folded){

 realRes *r=(realRes *) p->extras[index];
 int id=0;
 
 for(int i=0; (i<p->numSites);i++)
   if(r->oklist[i]==1)
     if(folded==0)
       fwrite(r->pLikes[id++],sizeof(double),2*p->nInd+1,outfileSFS);
     else
       fwrite(r->pLikes[id++],sizeof(double),p->nInd+1,outfileSFS);
 for(int i=0;i<p->numSites;i++)
   if(r->oklist[i]==1)
     fprintf(outfileSFSPOS,"%s\t%d\n",chr,p->posi[i]+1);
   else if (r->oklist[i]==2)
     fprintf(stderr,"PROBS at: %s\t%d\n",chr,p->posi[i]+1);

}


void realSFS::fin(funkyPars *pars,int index,double *prior,gzFile fpgz){
 realRes *r=(realRes *) pars->extras[index];
 int id=0;
 
 for(int i=0; (i<pars->numSites);i++){
   
   if(r->oklist[i]==1){
     double *tmp = r->pLikes[id++];
     for(int ii=0;i<2*pars->nInd+1;ii++)
       tmp[ii] += log(prior[ii]); //take loglike from sfs, and add the logprior
   
     normalize_array2(tmp,2*pars->nInd+1);

     //now tmp contains our posterior expectation of the different classes of frequencies
     double *workarray = tmp;
    
     //First find thetaW: nSeg/a1
     double pv = 1-exp(workarray[0])-exp(workarray[2*pars->nInd]);
     double seq;
     if(pv<0)
       seq=log(0.0);
     else
       seq =log(1-exp(workarray[0])-exp(workarray[2*pars->nInd]))-aConst;
     //     fprintf(stderr,"seq=%f exp0=%f exp2=%f \n",seq,exp(workarray[0]),exp(workarray[2*pars->nInd]));

    double pairwise=0;    //Find theta_pi the pairwise
    double thL=0;    //Find thetaL sfs[i]*i;
    double thH=0;//thetaH sfs[i]*i^2
    for(size_t ii=1;ii<2*pars->nInd;ii++){

      pairwise += exp(workarray[ii]+scalings[ii] );
      double li=log(ii);
      
      thL += exp(workarray[ii])*ii;
      thH += exp(2*li+workarray[ii]);
    }
    gzprintf(fpgz,"%s\t%d\t%f\t%f\t%f\t%f\t%f\n",header->name[pars->refId],pars->posi[i]+1,seq,log(pairwise)-aConst2,workarray[1],log(thH)-aConst2,log(thL)-aConst3);
   }else if(r->oklist[i]==2)
     fprintf(stderr,"PROBS at: %s\t%d\n",header->name[pars->refId],pars->posi[i]+1);
 }
}




void realSFS::print(funkyPars *p){
 if(doRealSFS==0||p->numSites==0)
    return;
 
 if(doRealSFS>0&&prior==NULL)
   printFull(p,index,outfileSFS,outfileSFSPOS,header->name[p->refId],fold);
 else if(doThetas==1&&doRealSFS==1&&prior!=NULL)
   fin(p,index,prior,fpgz);
}




/*
//Functions below  should be added 
*/
/*


double myComb2(int k,int r, int j){
  //return 1.0;
  //  return 1.0/bico(2*k,2);
  if(j>r)
    fprintf(stderr,"%s error in k=%d r=%d j=%d\n",__FUNCTION__,k,r,j);
      //    return 0;
  double fac1= bico(r,j)*bico(2*k-r,2-j);
  double fac2=bico(2*k,2);
  //  fprintf(stderr,"[%s]=%f\t=%f\tres=%f\n",__FUNCTION__,fac1,fac2,fac1/fac2);
  return fac1/fac2;
  //return fac1;
}





void algoGeno(double **liks,int *major,int *minor,int nsites,int numInds,FILE *sfsfile,int underFlowProtect, int fold,loci *locs,FILE *sfssites,int *keepSites,double *pest) {

  //void algoGeno(aMap &asso,int numInds,FILE *sfsfile,int underFlowProtect, double *pest) {
  fprintf(stderr,"UNDERFLOWPROTECT: %d\n",underFlowProtect);
  
  for(int r=0;0&&r<(2*numInds-1);r++)
    for(int j=0;j<3;j++){
      double res = myComb2(numInds,r,j);
      fprintf(stderr,"(%d,%d,%d) =%f\n",numInds,r,j,res);
    }

  
  //algorithm goes on by a site on site basis //pretty much the same as 'algo' without the prior phat
  
  //  underFlowProtect = 0;
  for(int it=0; it<nsites; it++) {//loop over sites
    int minor_offset=minor[it];
    int major_offset=major[it];
    if(keepSites[it]==0)
      continue;

    //    fprintf(sfsfile,"%s %d\t",chrnames[it->first.chromo],it->first.position);
    //1. first loop through all possible major/minor

    //    fprintf(sfsfile,"%d %d\t0 0 0 0\t",major_offset,minor_offset);
    int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
    int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
    int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
    //fprintf(stderr,"OFFSET_INFO: min=%d\tmaj=%d aa=%d Aa=%d AA=%d\n",minor_offset,major_offset,aa_offset,Aa_offset,AA_offset);
    //part two
    double hj[2*numInds+1];
    for(int index=0;index<(2*numInds+1);index++)
      if(underFlowProtect==0)
	hj[index]=0;
      else
	hj[index]=log(0);
    double PAA,PAa,Paa;
    //    numInds =5;
    for(int i=0 ; i<numInds ;i++) {
      //	printf("AA=%f\tAa=%f\taa=%f\n",p.lk[i*3+AA_offset],p.lk[i*3+Aa_offset],p.lk[i*3+aa_offset]);
      double GAA,GAa,Gaa;
      GAA = liks[it][i*10+AA_offset];
      GAa = log(2.0)+liks[it][i*10+Aa_offset];
      Gaa = liks[it][i*10+aa_offset];

      //fprintf(stderr,"GAA=%f\tGAa=%f\tGaa=%f\n",GAA,GAa,Gaa);
      if(underFlowProtect==0){
	GAA=exp(GAA);
	GAa=exp(GAa);
	Gaa=exp(Gaa);
      }
	  
      
      PAA =(GAA);///(MAA+MAa+Maa);
      PAa =(GAa);///(MAA+MAa+Maa);
      Paa =(Gaa);///(MAA+MAa+Maa);
      
      
      //check for underflow error, this should only occur once in a blue moon
      if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	printf("PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
      }
      
      if(i==0){
	hj[0] =Paa;
	hj[1] =PAa;
	hj[2] =PAA;
      }else{
	
	for(int j=2*(i+1); j>1;j--){

	  double tmp;
	  if(underFlowProtect==1)
	    tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	  else
	    tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	      
	  if(std::isnan(tmp)){
	    fprintf(stderr,"jis nan:%d\n",j );
	    hj[j] = 0;
	    break;
	  }else
	    hj[j]  =tmp;

	}
	if(underFlowProtect==1){
	  hj[1] = angsd::addProtect2(Paa+hj[1],PAa+hj[0]);
	  hj[0] = Paa+hj[0];
	}
	else{
	  hj[1] = Paa*hj[1] + PAa*hj[0];
	  hj[0] = Paa*hj[0];
	}
	//	print_array(stdout,hj,2*numInds+1,!underFlowProtect);
      }
      
      
    }//after recursion
    //if we are underflowprotecting ht hj is in logspace
    // print_array(stdout,hj,2*numInds+1,!underFlowProtect);
    for(int i=0;i<(2*numInds+1);i++){
      //fprintf(stdout,"BICO: %f\n",log(bico(2*numInds,i)));
      if(underFlowProtect)
	hj[i] =  (hj[i]-log(bico(2*numInds,i)));
      else
	hj[i] =  exp(log(hj[i])-log(bico(2*numInds,i)));

    }
    //fprintf(stdout,"the full after update hj\n");
    //print_array(stdout,hj,2*numInds+1,!underFlowProtect);
    double denominator = 0;
    if(underFlowProtect)
      denominator = log(denominator);
    for(int i=0;i<=2*numInds;i++)
      if(underFlowProtect)
	denominator = angsd::addProtect2(denominator,log(pest[i])+angsd::addProtect2(hj[i],hj[2*numInds-i]));
      else
	denominator += pest[i] * (hj[i] +hj[2*numInds-i]);

    int whichGeno[numInds];
    double whichProb[numInds];

    for(int select=0;select<numInds;select++) {
      double *hj = new double[2*numInds-1];
      for(int index=0;index<(2*numInds-1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,PAa,Paa;
      
      int ishadow =-1;
      for(int i=0 ; i<numInds ;i++) {
	//	printf("AA=%f\tAa=%f\taa=%f\n",p.lk[i*3+AA_offset],p.lk[i*3+Aa_offset],p.lk[i*3+aa_offset]);
	if(i!=select){
	  ishadow++;
	}else
	  continue;
	double GAA,GAa,Gaa;
	GAA = liks[it][i*10+AA_offset];
	GAa = log(2.0)+liks[it][i*10+Aa_offset];
	Gaa = liks[it][i*10+aa_offset];
	
	if(underFlowProtect==0){
	  GAA=exp(GAA);
	  GAa=exp(GAa);
	  Gaa=exp(Gaa);
	}
	  
	
	PAA =(GAA);///(MAA+MAa+Maa);
	PAa =(GAa);///(MAA+MAa+Maa);
	Paa =(Gaa);///(MAA+MAa+Maa);


	//check for underflow error, this should only occur once in a blue moon
	if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	  fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
	}
	
	if(ishadow==0){
	  hj[0] =Paa;
	  hj[1] =PAa;
	  hj[2] =PAA;
	}else{
	  
	  for(int j=2*(ishadow+1); j>1;j--){
	    //	    fprintf(stderr,"j=%d\n",j);
	    //print_array(hj,2*numInds+1);
	    double tmp;
	    if(underFlowProtect==1)
	      tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	    else
	      tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	    
	    if(std::isnan(tmp)){
	      fprintf(stderr,"jis nan:%d\n",j );

	      hj[j] = 0;
	      break;
	    }else
	      hj[j]  =tmp;
	  }
	  if(underFlowProtect==1){
	    hj[1] = angsd::addProtect2(Paa+hj[1],PAa+hj[0]);
	    hj[0] = Paa+hj[0];
	  }
	  else{
	    hj[1] = Paa*hj[1] + PAa*hj[0];
	    hj[0] = Paa*hj[0];
	  }
	}
	//ifunderflowprotect then hj is in logspace
	
      }//after recursion
      for(int i=0;i<(2*(numInds-1)+1);i++){
	//fprintf(stdout,"BICO: %f\n",log(bico(2*numInds,i)));
	if(underFlowProtect)
	  hj[i] =  (hj[i]-log(bico(2*(numInds-1),i)));
	else
	  hj[i] =  exp(log(hj[i])-log(bico(2*(numInds-1),i)));
	
      }
       
      //now do all the genocalling for individual =select
      // fprintf(stderr,"seelct=%d\tmyMaj=%d\tmyMin=%d\tAA_offset=%d Aa_offset=%d aa_offset=%d\n",select,p.major,p.minor,AA_offset,Aa_offset,aa_offset);
      double *asdf = hj;

      //print_array(stdout,asdf,2*numInds-1,!underFlowProtect);
      double res[3]; //is always logged


      for (int j=0;j<3;j++){//loop through 3 genotypes
	//	fprintf(stderr,"jis: %d\n",j);
	double g;
	if(j==0)
	  g=liks[it][10*select+angsd::majorminor[major[it]][major[it]]];
	else if(j==1)
	  g=liks[it][10*select+angsd::majorminor[major[it]][minor[it]]];
	else
	  g=liks[it][10*select+angsd::majorminor[minor[it]][minor[it]]];
	
	double tmp=0;
	if(underFlowProtect)
	  tmp = log(tmp);
	for(int r=0;r<2*numInds-1;r++){
	  //	  fprintf(stderr,"mycomb2:%f\tpes1: %f\tpes2: %f\n",myComb2(numInds,r+j,j),pest[r+j],pest[2*numInds-r-j]);
	  if(underFlowProtect)
	    tmp = angsd::addProtect2(tmp, log(myComb2(numInds,r+j,j)) + (asdf[r])+log(pest[r+j]+pest[2*numInds-r-j]));
	  else
	    tmp += myComb2(numInds,r+j,j)*(asdf[r])*(pest[r+j]+pest[2*numInds-r-j]);
	}
	//g is directly from glf file, tmp is in log depending on underflowprotect

	if(underFlowProtect)
	  res[j]=g+(tmp);
	else
	  res[j]=g+log(tmp);
	//	fprintf(stderr,"j:%d log(g)=%f tmp=%f g*tmp=%f res[%d]=%f logres[%d]=%f \n",j,g,log(tmp),(g*tmp),j,res[j],j,log(res[j]));

	//res is always in log
	
	

	//CODE BELOW IS A CHECK
	/*
	double tmp1=0;
	double tmp2=0;
	if(underFlowProtect){
	  tmp1 = log(tmp1);
	  tmp2 = log(tmp2);
	}
	  
	for(int r=j;r<=2*numInds-2+j;r++){{
	    // fprintf(stderr,"tmp1: r=%d index: %d\n",r,r-j);
	    //fprintf(stderr,"logcomb: %f\n",myComb2(numInds,r,j));
	    //    fprintf(stderr,"tmp1 in preloop: %f\n",tmp1);
	    if(underFlowProtect)
	      tmp1 = angsd::addProtect2(tmp1, log(myComb2(numInds,r,j))+ asdf[r-j]+log(pest[r]));
	    else
	      tmp1 += myComb2(numInds,r,j)*(asdf[r-j])*(pest[r]);
	    //fprintf(stderr,"tmp1 in reloop: %f\n",tmp1);
	  }
	  
	}

	for(int r=2-j;r<=2*numInds-j;r++){
	  //fprintf(stderr,"tmp2: r=%d index=%d\n",r,r-2+j);
	  //fprintf(stderr,"logcomb: %f\n",myComb2(numInds,r,2-j));
	  if(underFlowProtect)
	    tmp2 = angsd::addProtect2(tmp2,log(myComb2(numInds,r,2-j))+ asdf[(2*(numInds-1))-(r-2+j)]+ log(pest[r]));
	  else
	    tmp2 += myComb2(numInds,r,2-j)*(asdf[2*(numInds-1)-(r-2+j)])*(pest[r]);
	}
	//	fprintf(stderr,"tmp1: %f tmp2: %f\n",(tmp1),(tmp2));	
	double res2;//always log

	if(underFlowProtect)
	  res2 = g+(angsd::addProtect2(tmp1,tmp2));
	else
	  res2 = g+log(tmp1+tmp2);
	//fprintf(stdout,"%f vs %f\n",res[j],res2);
	if(1&&!isSame(res[j],res2,0.000001)){
	  fprintf(stdout,"%f vs %f\n",res[j],res2);
	}
	//check that is sums to one according the model

      
	  /*
	if(underFlowProtect)
	  fprintf(stdout,"\tASDFASDFASDF: %f\n",log((res2-denominator)));
	else
	  fprintf(stdout,"\tASDFASDFASDF: %f\n",log(exp(res2)/denominator));	

	//CODE ABOVE IS CHECK



      }//after the loop of the tree genoypes
      
      
      double shouldBeOne =0;
      if(underFlowProtect)
	shouldBeOne = exp(res[0]-(denominator))+exp(res[1]-(denominator))+exp(res[2]-(denominator));
      else
	shouldBeOne = exp(res[0]-log(denominator))+exp(res[1]-log(denominator))+exp(res[2]-log(denominator));
      if((fabs(shouldBeOne-1)>0.000001)){
	fprintf(stderr,"this should be one: %f\n",shouldBeOne);
	exit(0);
      }
      


      //      print_array(stdout,res,3,0);      
      
      //logrescale(res,3);
      double mySum=exp(res[0])+exp(res[1])+exp(res[2]);
      for(int i=0;i<3;i++)
	res[i] =exp(res[i])/mySum;
      int best = whichMax(res,3);
      if(best==0)
	whichGeno[select] = angsd::majorminor[major[it]][major[it]];
      if(best==1)
	whichGeno[select] = angsd::majorminor[major[it]][minor[it]];
      if(best==2)
	whichGeno[select] = angsd::majorminor[minor[it]][minor[it]];
      whichProb[select] = res[best];
      //print_array(sfsfile,res,3,0);      

    }//after select loop  
    /*
    print_array(sfsfile,whichGeno,numInds);
    print_array(sfsfile,whichProb,numInds,0);

    for(int i=0;i<numInds-1;i++)
      fprintf(sfsfile,"%d ",whichGeno[i]);
    fprintf(sfsfile,"%d\t",whichGeno[numInds-1]);
    for(int i=0;i<numInds-1;i++)
      fprintf(sfsfile,"%f ",whichProb[i]);
    //    fprintf(sfsfile,"%f\t%f\n",whichProb[numInds-1],it->second.emPhat);
  }//after all sites
  
}
	*/

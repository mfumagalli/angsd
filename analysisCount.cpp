 /*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd
  Class that works with counts
*/

#include <zlib.h>
#include <cassert>
#include "analysisFunction.h"
#include "general.h"
#include "analysisCount.h"
#include "kstring.h"




void countCls::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doCounts\t%d\t(Count the number A,C,G,T. All sites, All samples)\n",doCounts);
  
  
  fprintf(argFile,"\t-minQfile\t%s\t file with individual quality score thresholds)\n",minQfile);
  fprintf(argFile,"\t-setMaxDepth\t%d\t(If total depth is larger then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMaxDepth);
  fprintf(argFile,"\t-setMinDepth\t%d\t(If total depth is smaller then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMinDepth);
  fprintf(argFile,"\t-trim\t\t%d\t(trim ends of reads)\n",trim);

  fprintf(argFile,"\t-minInd\t\t%d\t(Discard site if effective sample size below value.\n\t\t\t\t 0 indicates no filtering)\n",minInd);
  fprintf(argFile,"Filedumping:\n");
  fprintf(argFile,"\t-doDepth\t%d\t(dump distribution of seqdepth)\t%s,%s\n",doDepth,postfix4,postfix5);
  fprintf(argFile,"\t  -maxDepth\t%d\t(bin together high depths)\n",maxDepth);
  
  fprintf(argFile,"\t-doQsDist\t%d\t(dump distribution of qscores)\t%s\n",doQsDist,postfix3);  
  fprintf(argFile,"\t-dumpCounts\t%d\n",dumpCounts);
  fprintf(argFile,"\t  1: total seqdepth for site\t%s\n",postfix1);
  fprintf(argFile,"\t  2: seqdepth persample\t\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  3: A,C,G,T sum over samples\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  4: A,C,G,T sum every sample\t%s,%s\n",postfix1,postfix2);
}

int calcSum(suint *counts,int len){
  int tmp=0;
  for(int i=0;i<len;i++)
    tmp += counts[i];
  return tmp;
}

void printCounts(char *chr,int *posi,suint **counts,int nSites,size_t nInd,kstring_t &bpos,kstring_t &bbin,int dumpType,int *keepSites){
  bpos.l=bbin.l=0;

  for(int s=0;s<nSites;s++){
    if(keepSites[s]==0)
      continue;
    ksprintf(&bpos, "%s\t%d\t%d\n",chr,posi[s]+1,calcSum(counts[s],4*nInd));
    
    //if we need per sample info
    if(dumpType>1) {
      if(dumpType==4)//count A,C,G,T
	for(int i=0;i<4*nInd;i++)
	  ksprintf(&bbin,"%u\t",counts[s][i]);
      else if(dumpType==2){//print A+C+G+T
	for(int n=0;n<nInd;n++)
	  ksprintf(&bbin,"%u\t",counts[s][n*4]+counts[s][n*4+1]+counts[s][n*4+2]+counts[s][n*4+3]);
      }else{//overall sum of A,C,G,T
	size_t tsum[4]={0,0,0,0};
	for(int i=0;i<4*nInd;i++)
	  tsum[i%4] +=counts[s][i];
	ksprintf(&bbin,"%zu\t%zu\t%zu\t%zu",tsum[0],tsum[1],tsum[2],tsum[3]);
      }
      kputc('\n',&bbin);	
    }
  }

}


void countCls::getOptions(argStruct *arguments){

  //from command line
  minQfile=angsd::getArg("-minQfile",minQfile,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  dumpCounts=angsd::getArg("-dumpCounts",dumpCounts,arguments);
    trim=angsd::getArg("-trim",trim,arguments);
  doQsDist=angsd::getArg("-doQsDist",doQsDist,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);
  setMaxDepth = angsd::getArg("-setMaxDepth",setMaxDepth,arguments);
  setMinDepth=angsd::getArg("-minDepth",setMinDepth,arguments);
  doDepth=angsd::getArg("-doDepth",doDepth,arguments);
  maxDepth=angsd::getArg("-maxDepth",maxDepth,arguments);
  

  if(dumpCounts&&doCounts==0){
    fprintf(stderr,"You must supply -doCounts if you want to dumpcounts\n");
    exit(0);
  }

  if(doDepth!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want depth distribution");
    exit(0);
  }

 if(doQsDist!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want qscore distribution");
    exit(0);
  }



}
//constructor
countCls::countCls(const char *outfiles,argStruct *arguments,int inputtype){
  
  nInd=arguments->nInd;
  minInd = 0;
  setMinDepth =-1;
  trim =0;
  dumpCounts =0;
  doCounts = 0;
  doQsDist = 0;
  
  doDepth = 0;
  maxDepth = 100;
  setMaxDepth = -1;
  minQfile=NULL;

  //make output files
  postfix1=".pos.gz";
  postfix2=".counts.gz";
  postfix3=".qs";
  postfix4=".depthSample";
  postfix5=".depthGlobal";
  bpos.s=NULL;bpos.l=bpos.m=0;
  bbin.s=NULL;bbin.l=bbin.m=0;

  //from command line
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doCounts")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  oFiles=strdup(outfiles);
  printArg(arguments->argumentFile);
  
  //  oFileCountsPos = oFileCountsBin = oFileQs = NULL;
  oFileCountsPos = oFileCountsBin =  Z_NULL;

  if(dumpCounts){
    oFileCountsPos = aio::openFileGz(outfiles,postfix1,GZOPT);
    gzprintf(oFileCountsPos,"chr\tpos\ttotDepth\n");
    if(dumpCounts>1)
      oFileCountsBin = aio::openFileGz(outfiles,postfix2,GZOPT);
    if(dumpCounts==2)
      for(int i=0;i<arguments->nInd;i++)
	gzprintf(oFileCountsBin,"ind%dTotDepth\t",i);
    if(dumpCounts==3)
      gzprintf(oFileCountsBin,"totA\ttotC\ttotG\ttotT");
    if(dumpCounts==4)
      for(int i=0;i<arguments->nInd;i++)
	gzprintf(oFileCountsBin,"ind%d_A\tind%d_C\tind%d_G\tind%d_T\t",i,i,i,i);
    if(dumpCounts>1)
      gzprintf(oFileCountsBin,"\n");
  }

  if(doQsDist){
    //datastructures needed
    qsDist = new size_t[256];
    memset(qsDist,0,256*sizeof(size_t));
    //prepare outputfile
    
  }
  if(doDepth){
    depthCount=new size_t *[arguments->nInd];
    for(int i=0;i<nInd;i++)
      depthCount[i]=new size_t[maxDepth+1];
    for(int i=0;i<nInd;i++)
      for(int j=0;j<maxDepth+1;j++)
	depthCount[i][j]=0;
    
    globCount = new size_t[maxDepth+1];
    memset(globCount,0,sizeof(size_t)*(maxDepth+1));

    
  }
  if(minQfile!=NULL){
    minQmat = angsd::getMatrix(minQfile,0,100000);
    if(minQmat.x!=nInd){
      fprintf(stderr,"Number of lines in the minQfile does not match the number of individuals \n");
      exit(0);
    }
    if(!(minQmat.y==1||minQmat.y==4)){
      fprintf(stderr,"Number of colums in the minQfile has to be 1 or 4 \n");
      exit(0);
    }
  }
  
  
}


void printQs(FILE *fp,size_t *ary){
  int firstidx=0;
  for(int i=0;i<256;i++)
     if(ary[i]!=0){
      firstidx=i;
      break;
    }
  int lastidx=255;
  for(int i=255;i>=0;i--)
    if(ary[i]!=0){
      lastidx=i;
      break;
    }
  for(int i=firstidx;i<=lastidx;i++)
    fprintf(fp,"%d\t%lu\n",i,ary[i]);
  
}


countCls::~countCls(){
  if(oFileCountsBin!=Z_NULL)    gzclose(oFileCountsBin);
  if(oFileCountsPos!=Z_NULL)    gzclose(oFileCountsPos);
  if(doQsDist){
    FILE *oFileQs = NULL;
    oFileQs = aio::openFile(oFiles,postfix3);
    fprintf(oFileQs,"qscore\tcounts\n");
    printQs(oFileQs,qsDist);
    if(oFileQs) fclose(oFileQs);
    delete[] qsDist;

  }

  if(doDepth){
    FILE *oFileSamplDepth = aio::openFile(oFiles,postfix4);
    FILE *oFileGlobDepth = aio::openFile(oFiles,postfix5);
    for(int i=0;i<nInd;i++){
      for(int j=0;j<maxDepth+1;j++){
	fprintf(oFileSamplDepth,"%lu\t",depthCount[i][j]);
      }
      fprintf(oFileSamplDepth,"\n");
    }
    //thorfinn
    for(int j=0;j<maxDepth+1;j++)
      fprintf(oFileGlobDepth,"%lu\t",globCount[j]);
    fprintf(oFileGlobDepth,"\n");
  

    //clean depthCount
    for(int i=0;i<nInd;i++)
      delete[]  depthCount[i];
    delete[] depthCount; 
    
    if(oFileSamplDepth) fclose(oFileSamplDepth);
    if(oFileSamplDepth) fclose(oFileGlobDepth);
  }
  
  if(minQfile!=NULL){
    //  angsd::printMatrix(minQmat,stderr);
    angsd::deleteMatrix(minQmat);
  }

  free(oFiles);
  free(bpos.s);
  free(bbin.s);
}

void countQs(const chunkyT *chk,size_t *ret,int trim,int *keepSites){
  
  suint **cnts = new suint*[chk->nSites];
  for(int s=0;s<chk->nSites;s++){
    if(keepSites[s]==0)
      continue;
    //loop over sites
    for(int n=0;n<chk->nSamples;n++){
      //loop over samples
      for(int l=0;l<chk->nd[s][n].l;l++){
	//loop over persample reads for this position/sample
	if(chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim)
	  continue;
	ret[chk->nd[s][n].qs[l]]++;
	
      }
    }
  }
}



void countCls::print(funkyPars *pars){

  if(dumpCounts)
    printCounts(header->name[pars->refId],pars->posi,pars->counts,pars->numSites,pars->nInd,bpos,bbin,dumpCounts,pars->keepSites);
  gzwrite(oFileCountsBin,bbin.s,bbin.l);
  gzwrite(oFileCountsPos,bpos.s,bpos.l);

  if(doQsDist)
    countQs(pars->chk,qsDist,trim,pars->keepSites);
  
  if(doDepth!=0){
    assert(pars->counts!=NULL);
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      for(int i=0;i<pars->nInd;i++){
	int sum=0;
	for(int a=0;a<4;a++)
	  sum+=pars->counts[s][i*4+a];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
	depthCount[i][sum]++;	
      }
    }
    //thorfinn below
    
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      int sum=0;
      for(int i=0;i<4*pars->nInd;i++){
	sum+=pars->counts[s][i];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
      }
      globCount[sum]++;
    } 
  }

}


void countCls::clean(funkyPars *pars){
  if(doCounts||dumpCounts){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->counts[i];
    delete [] pars->counts;
  }
}


//dragon update with keeplist so we only count necessary sites
suint **countCls::countNucs(const chunkyT *chk,int trim,int *keepSites){
  suint **cnts = new suint*[chk->nSites];
  if(minQfile==NULL){
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;l<chk->nd[s][n].l;l++){
	  int allele = refToInt[chk->nd[s][n].seq[l]];
	  if(chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim||allele==4){
	    //fprintf(stderr,"ind %d,allele %d\tminQ %d\tQ %d\tMapQ %d\n",n,allele,minQ,chk->nd[s][n].mapQ[l]);
	    continue;
	  }
	  cnts[s][4*n+allele]++;
	}
      }
    }
  }
  else if(minQmat.y>1){
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;l<chk->nd[s][n].l;l++){
	  int allele = refToInt[chk->nd[s][n].seq[l]];	  
	  if(allele==4)
	    continue;
	  if(chk->nd[s][n].qs[l] < minQmat.matrix[n][allele]||chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim){
	    //fprintf(stderr,"ind %d,allele %d\tminQ %d\tQ %d\tMapQ %d\n",n,allele,minQ,chk->nd[s][n].mapQ[l]);
	  continue;
	  }
	  cnts[s][4*n+allele]++;
	}
      }
    }

  }
 else{
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;l<chk->nd[s][n].l;l++){
	  int allele = refToInt[chk->nd[s][n].seq[l]];
	  if(chk->nd[s][n].qs[l] < minQmat.matrix[n][0]||chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim||allele==4){
	    //fprintf(stderr,"ind %d,allele %d\tminQ %d\tQ %d\tMapQ %d\n",n,allele,minQ,chk->nd[s][n].mapQ[l]);
	  continue;
	  }
	  cnts[s][4*n+allele]++;
	}
      }
    }

  }
  return cnts;
}




void countCls::run(funkyPars *pars){
  if(doCounts==0)
    return;
  assert(pars->chk!=NULL&&pars->counts==NULL);
  pars->counts = countNucs(pars->chk,trim,pars->keepSites);
  // fprintf(stderr,"%d\n",pars->keepSites[0]);
  //modify keepsites;
  if(minInd!=0) {
    for(int i=0;i<pars->numSites;i++){
      if(pars->keepSites[i]==0)
	continue;
      //THIS IS REALLY STUPID but lets count number of samples wiht info
      int nDep =0;
      for(int s=0;s<pars->nInd;s++){
	int dep=0;
	for(int j=0;j<4;j++)
	  dep += pars->counts[i][s*4+j];
	if(dep)
	  nDep++;
      }
      //nDep is now the number of sapmles wiht info
      if(nDep<minInd)
	pars->keepSites[i] = 0;
      else
	pars->keepSites[i] = nDep;
      
    }
  }
  if(setMaxDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      size_t totSum = calcSum(pars->counts[s],4*nInd);
      if(totSum>setMaxDepth)
	pars->keepSites[s]=0;
      else{
	suint *ps = pars->counts[s];
	for(int i=0;i<pars->nInd;i++){
	  int iSum = ps[i*4]+ps[i*4+1]+ps[i*4+2]+ps[i*4+2]+ps[i*4+3];
	  if(totSum>iSum){
	    pars->keepSites[s]=0;
	    break;
	  }
	}
      }
    }
  }
  //fprintf(stderr,"minDepth=%d nInd=%d pars->numsites=%d\n",minDepth,nInd,pars->numSites);
  if(setMinDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      size_t totSum = calcSum(pars->counts[s],4*pars->nInd);
      if(totSum<setMinDepth)
	pars->keepSites[s]=0;
      else{
	int nDep =0;
	for(int i=0;i<pars->nInd;i++){
	  int dep=0;
	  for(int j=0;j<4;j++)
	    dep += pars->counts[s][i*4+j];
	  if(dep)
	    nDep++;
	}
	pars->keepSites[s]= nDep;

      }
    }
  }
  
}

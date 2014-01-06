/*
  small class to filter positions and assign major minor.
  This class can be improved very much.

  Thorfinn 7march 2013

  On the most basic level then this class will remove those sites with an effective sample size<minInd

  It can also be used for filtering away those sites not included in the -filter file.keep

  Or you can supply a bimfile, and then the major/minor from the bim file will be used as major/minor throughout the program


 */
#include <cassert>
#include "pthread.h"
#include "shared.h"
#include "general.h"
#include "analysisFunction.h"
#include "analysisKeepList.h"


//this function is much to slow on a genome scale should be improved
fMap getMap(char *fname,std::map<char*,int,ltstr> *revMap){
  const char *delims = "\t\n ";
  FILE *fp=getFILE(fname,"r");
  
  char buf[LENS];
  int nsites=0;
  fMap ret;
  std::map<char *,int>::iterator rit;
  while(fgets(buf,LENS,fp)){
    char *chr = strtok(buf,delims);
    strtok(NULL,delims);//rsnumber
    strtok(NULL,delims);//centimorgan
    char *tok = strtok(NULL,delims);
    if(tok==NULL){
      fprintf(stderr,"Problem with fileformat in .bim file\n");
      exit(0);
    }
    int pos = atoi(tok)-1;//genomic position in bp
    mm value;
    value.major = refToInt[strtok(NULL,delims)[0]];
    value.minor = refToInt[strtok(NULL,delims)[0]];

    //check for N if this exists;
    if(value.major==4||value.minor==4){
      fprintf(stderr,"N extists in major minor defintion\n");
      break;
    }
    //    fprintf(stderr,"chr=%s pos=%d major=%d minor=%d\n",chr,pos,mymm.major,mymm.minor);
    rit=revMap->find(chr);
    if(rit==revMap->end()){
      fprintf(stderr,"Problem finding chromosome: %s in lookuptable\n",chr);
      exit(0);
    }

    mm key;
    key.major = rit->second;
    key.minor = pos;
    fMap::iterator it = ret.find(key);
    if(it!=ret.end()){
      fprintf(stderr,"duplicate entry in filterlist:%s : will exit offending position below\n",fname);
      fprintf(stderr,"chr=%s pos=%d major=%d minor=%d\n",chr,pos,value.major,value.minor);
      exit(0);
    }else
      ret.insert(fMap::value_type(key, value));
  }
  fclose(fp);
  return ret;
}


filter::~filter(){
  if(doFilter==0)
    return;
  //delete [] keepsChr;

}


filter::filter(argStruct *arguments){
  //below if shared for all analysis classes
  header = arguments->hd;
  revMap = arguments->revMap;

  //his is used by this class
  keepsChr = NULL;
  curChr = -1;
  fp = NULL;
  minInd = -1;
  fname = NULL;
  doMajorMinor =0;
  doFilter =0;

  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-filter")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  //get options and print them
  getOptions(arguments);
  printArg(arguments->argumentFile);

}


void filter::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-filter\t\t%s dofilter=%d\n",fname,doFilter);
  fprintf(argFile,"\t-doMajorMinor\t%d\t\n",doMajorMinor);
  fprintf(argFile,"\t1: Infer major and minor from GL\n");
  fprintf(argFile,"\t2: Infer major and minor from allele counts\n");
  fprintf(argFile,"\t3: use major and minor from bim file (requires -filter afile.bim)\n");
  fprintf(argFile,"\t4: Use reference allele as major (requires -ref)\n");
  fprintf(argFile,"\t5: Use ancestral allele as major (requires -anc)\n");
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if atleast minInd of samples has data\n",minInd);  fprintf(argFile,"\n");
}

//1=bim file 2=keep file
int findType( char *fname){
  char *dotStart = strrchr(fname,'.');
  //  fprintf(stderr,"%s",dotStart);
  if(0==strcmp(dotStart,".bim"))
    return 1;
  if(0==strcmp(dotStart,".keep"))
    return 2;
  fprintf(stderr,"Unknown filterfile, should be either .bim or .keep\n");
  exit(0);
}


void filter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-filter",fname,arguments);
  
  if(fname!=NULL)
    doFilter = findType(fname);
  //1=bim 2=keep

  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && doFilter!=1){
    fprintf(stderr,"Must supply -filter with .bim file if -doMajorMinor 3\n");
    exit(0);
  }
  

  if(doFilter==1){
    fm = getMap(fname,revMap);
    fprintf(stderr,"\t-> number of sites in filter: %lu\n",fm.size());
  }else if(doFilter==2){
    fp = getFILE(fname,"r");
    //  readSites();
    fprintf(stderr,"Filtering with .keep is still beta\n");
  }
  minInd = angsd::getArg("-minInd",minInd,arguments);
}


void filter::run(funkyPars *p){
  //  fprintf(stderr,"nsites=%d\n",p->numSites);
  p->keepSites=new int[p->numSites];
  //  p->results->freq->keepInd = new int[p->numSites];

  for(int s=0;s<p->numSites;s++){
    p->keepSites[s]=p->nInd;
    //    p->results->freq->keepInd[s]=nInd;  
  }
  if(doFilter==1 && doMajorMinor==3){
    p->major= new char [p->numSites];
    p->minor= new char [p->numSites];
    for(int i=0;i<p->numSites;i++){
      p->major[i] = 4;
      p->minor[i] = 4;
    }
  }

  if(doFilter==1) {
    for(int s=0;s<p->numSites;s++){
      mm key;
      key.major=p->refId;key.minor = p->posi[s];

      fMap::iterator it = fm.find(key);
      if(it==fm.end()){
	//fprintf(stderr,"not here\n");
	p->keepSites[s] =0;
	continue;
      }
      else if(doMajorMinor==3) {
	p->major[s] = it->second.major;
	p->minor[s] = it->second.minor;
	if(p->major[s]<0||p->major[s]>3)
	  p->keepSites[s] =0;
	if(p->minor[s]<0||p->minor[s]>3||p->minor==p->major)
	  p->keepSites[s] =0;
	

      }
    }
  }else if(doFilter==2){
    if(p->refId!=curChr){
      fprintf(stderr,"Problem with mismatch of filtering dataChr=%d filterChr=%d check ordering of chrs in filter file\n",p->refId,curChr);
      fflush(stderr);
      exit(0);
    }

    for(int s=0;s<p->numSites;s++){
      //      fprintf(stderr,"posi=%d reflength=%d\n",p->posi[s],header->l_ref[p->refId]);
      if(keepsChr==NULL||keepsChr[p->posi[s]]==0)
	p->keepSites[s] = 0;
      
    }
  }

  //how set the keepsites according the effective sample size persite
  if(-1!=minInd){
    if(p->chk!=NULL){
      //loop over sites;
      for(int s=0;s<p->numSites;s++){
	if(p->keepSites[s]==0)
	  continue;
	int nInfo =0;
	tNode *tn = p->chk->nd[s];
	//loop over samples;
	for(int i=0;i<p->nInd;i++)
	  if(tn[i].l!=0)
	    nInfo++;
	//fprintf(stdout,"%d %d\n",nInfo,minInd);
	if(minInd<=nInfo)
	  p->keepSites[s] =nInfo;
	else
	  p->keepSites[s] =0;
      }
    
    }
  }
}
void filter::print(funkyPars *p){
}

void filter::clean(funkyPars *p){
  
}


void filter::readSites() {
  if(doFilter!=2)
    return;
  fprintf(stderr,"[%s]\n",__FUNCTION__);
  assert(doFilter==2&&fp!=NULL);
  static int bufPos = -1;//position
  static int bufChr = -1;//chromosome id refID
  static char *chrName = NULL;//chromosome name
  delete [] keepsChr;   keepsChr = NULL;
  
  char buf[1024];
  std::map<char *,int,ltstr>::iterator it;

  //this conditional only happens at beginning of file.
  if(bufPos==-1&&bufChr==-1){
    fgets(buf,1024,fp);
    char *tok = strtok(buf,"\n \t");//chrname
    it = revMap->find(tok);
    int posi = atoi(strtok(NULL,"\n \t")) -1; //offset by one
    if(it==revMap->end()){
      fprintf(stderr,"Problem finding chr=%s\n",tok);
      exit(0);
    }
    if(posi>header->l_ref[it->second]){
      fprintf(stderr,"Position in keep list exeeds reference\n");
      exit(0);
    }
    bufPos =posi;
    bufChr = it->second;
  }
  
  if(curChr==bufChr){
    fprintf(stderr,"Couldn't read more positions for filtering, won't use more sites for analysis\n");
    extern int SIG_COND;
    SIG_COND=0;
    return;
  }
  
  fprintf(stderr,"reading a chr from posi file bufPos=%d bufCHr=%d revmapSize=%lu chrName=%s\n",bufPos,bufChr,revMap->size(),bufChr!=-1?header->name[bufChr]:NULL);
  fflush(stderr);

  keepsChr = new char[header->l_ref[bufChr]];
  memset(keepsChr,0,header->l_ref[bufChr]);
  keepsChr[bufPos] = 1;
  curChr = bufChr;
  chrName=header->name[bufChr];
  
  while(fgets(buf,1024,fp)) {
    char *tok = strtok(buf,"\n \t");
    int posi = atoi(strtok(NULL,"\n \t")) -1; //offset by one

    //most general case, the next position is on same chromosome as the last line.
    if(0==strcmp(chrName,tok)){
      if(posi>header->l_ref[bufChr]){
	fprintf(stderr,"Position in keep list exeeds reference\n");
	exit(0);
      }
      keepsChr[posi] = 1;
      continue;
    }

    //check if the chromosomeID matches the header of our data
    it = revMap->find(tok);
    if(it==revMap->end()){
      fprintf(stderr,"Problem finding chr=%s\n",tok);
      exit(0);
    }
    //check if length does no exceed reflength
    if(posi>header->l_ref[it->second]){
      fprintf(stderr,"Position in keep list exeeds reference\n");
      exit(0);
    }
    //if we are here then we should buffer whatthe site and break;
    bufChr=it->second;
    bufPos=posi;
    break;
  }
  int tsum= 0;
  for(int i =0;i<header->l_ref[bufChr];i++)
    if(keepsChr[i])tsum++;
  fprintf(stderr,"Done reading a chrID=%d with name=%s nSites=%d\n",curChr,chrName,tsum);

}

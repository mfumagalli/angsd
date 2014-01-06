
extern int SIG_COND;

#include "mrStruct.h"
#include "pileups.h"
#include "analysisFunction.h"

using namespace pileworld;
//to check if our std::map works
void print_cMap(const cMap& m){
  for(cMap::const_iterator it=m.begin(); it!=m.end(); it++)
    printf("%s\t%d\n",it->first,it->second);

}



cMap pileups::buildMap(const char* fname){
  fprintf(stderr,"\t Opening chromosome index file: :%s\n",fname);
  std::ifstream pfile;
  pfile.open(fname,std::ios::in);
  char buffer[LENS];
  int id=0;
  cMap ret;
  fprintf(stderr,"\t-> Building chromosome index\n");
  while(pfile.getline(buffer,LENS) ) {
    char* chr = strdup(strtok(buffer,"\t"));
    ret.insert(std::make_pair(chr,id++));
    fprintf(stderr,"\t-> %s->%d\n",chr,id-1);
  }
  fprintf(stderr,"\t-> %s->%d\n","Last",id);
  fprintf(stderr,"\t-> Done building chromosome index: nchromo=%lu\n",ret.size());
  const char *Last = "Last";
  if(ret.size()!=0){//nescearry hook for rasmus input
    chrnames = new char*[ret.size()+1];
    for(cMap::const_iterator cit=ret.begin();cit!=ret.end();++cit){
      chrnames[cit->second] = strdup(cit->first);
    }
    chrnames[ret.size()] = strdup(Last);
  }else{//htis is rasmus input
    fprintf(stderr,"\t-> supply -fai\n");
    exit(0);
  }



  return ret;
}

mr::funkyPars *pileups::collapse(aMap &asso){

  if(asso.size()==0)
    return NULL;
  mr::funkyPars *r = mr::allocFunkyPars();
  r->numSites = asso.size();
  r->sites = new loci[r->numSites];
  r->likes = new double*[r->numSites];

  int **mynewDepth = new int*[r->numSites];
  size_t pos =0;
  
  for(aMap::iterator it=asso.begin();it!=asso.end();it++){
    loci l = (loci) {strdup(chrnames[it->first.chromo]),it->first.position};
    r->sites[pos] = l;
    r->likes[pos] = it->second.lk;
    mynewDepth[pos] = it->second.datumDepth;
    pos++;
  }
  r->depth = mynewDepth;
  return r;
}


void pileups::parseLoci(lociPileup &l,const char *str,const cMap &cmap){
  //  fprintf(stderr,"adfasdfadfasfsadfsadfasdfasdf\n");
  const char *delims = ":";
  //free(l.chromo);
  char *tmptmp = strdup(str);//TSK
  char *chr = strtok(tmptmp,delims);//TSK
  int pos = atoi(strtok(NULL,delims));

  cMap::const_iterator cit = cmap.find(chr);

  if(cit==cmap.end()){
      printf("chr: %s doesn't exist in faiIndex will exit\n",chr);
      exit(0);
  }



  l.chromo=cit->second;
  l.position = pos;
  free(tmptmp);//TSK
}


void pileups::printArg(FILE *argFile){
  fprintf(argFile,"%s:\n\n",__FILE__);
  fprintf(argFile,"-samglf OR samglfclean\t%s\n",fnames);
  fprintf(argFile,"-fai\t%s\n",faifile);
  fprintf(argFile,"-lStart %s -lStop\t%s\n",lStart,lStop);
  fprintf(argFile,"-nInd\t%d\n",nInd);
  fprintf(argFile,"\n");
}


void pileups::getOptions(argStruct *arguments){
  nInd = -1;
  type=arguments->inputtype;
  faifile=angsd::getArg("-fai",faifile,arguments);
  lStart=angsd::getArg("-lStart",lStart,arguments);
  lStop=angsd::getArg("-lStop",lStop,arguments);
  nInd = angsd::getArg("-nInd",nInd,arguments);
  if(type==1)
    fnames=angsd::getArg("-samglf",fnames,arguments);
  else
    fnames=angsd::getArg("-samglfclean",fnames,arguments);
  
  fprintf(stderr,"fnames=%s\n",fnames);
  if(strcmp(fnames,"-999")==0){
    printArg(stdout);
    exit(0);
  }


  filenames= angsd::getFilenames(fnames,nInd);
  faiIndex = buildMap(faifile);
  arguments->nInd=filenames.size();
  printArg(stdout);
}


//pileups::pileups(progArgs *pg){
//pileups::pileups(char *faifile,char *lStart,char*lStop,std::vector<char *> &filenames,int type_a){
pileups::pileups(argStruct *arguments){
  nFiles =0;
  getOptions(arguments);
  ultraStart=(lociPileup){0,0};
  ultraStop=(lociPileup){0,0};


  if(lStart!=NULL)
    parseLoci(ultraStart,lStart,faiIndex);
  if(lStop!=NULL)
    parseLoci(ultraStop,lStop,faiIndex);
  else
    ultraStop.chromo = faiIndex.size();
  
  bobj = NULL;
  tobj = NULL;
  //  fprintf(stderr,"[%s] allocing:%lu\n",__FUNCTION__,filenames.size());
  obj = new glfClass*[filenames.size()];
  
  if(type==2){
    fprintf(stderr,"\t Assuming texttype glf\n");
    tobj = new txtInput[filenames.size()];
  }else{
    fprintf(stderr,"\t Assuming binary glf\n");
    bobj = new binInput[filenames.size()];
  }
  for(int i=0;i<(int)filenames.size();i++)
    if(type==2){
      //  fprintf(stderr,"\t allocing tobj");
      obj[i] = &tobj[i];
    }else{
      // fprintf(stderr,"\t allocing bobj");
      obj[i] = &bobj[i];
    }
  for(int i=0;1&&i<(int)filenames.size();i++){
    obj[i]->init(i,filenames[i],ultraStart,ultraStop,faiIndex);
  }
  //  return 0;
  keepGoing = new int[filenames.size()];
  for(int i=0;i<filenames.size();i++)
    keepGoing[i] =1;
  filesLeft = (int) filenames.size() +1;//we need one final iteration for cleaning up the buffers
  
  nFiles = filenames.size();

}



mr::funkyPars *pileups::fetch(int nLines,int chunkSize){
  //fprintf(stderr,"nfiles:%d\n",nFiles);
  //used for killing the program and flushing the buffers. program will still segfault.
  //initialize signal handler
  
  aMap asso;
  while(SIG_COND&&(filesLeft)) {
    for(int i=0;i<nFiles;i++){
      //      fprintf(stderr,"[%lu]is total:%d\n",i,obj[i]->getTotal());
      obj[i]->flush(asso);//flush buffered data from previous seq line.
    }
    int indexed =0;
    
    //    fprintf(stderr,"\t-> filesLeft: %d\n",filesLeft);
    for(size_t i=0;i<nFiles;i++){//loop through all files
      
      if(keepGoing[i]==0)//if where are done with file i skip it.
	continue;
      else if(indexed==0) {//if we havent read nLines from any file yet
	keepGoing[i] = obj[i]->readlines(asso,ultraStop,nLines);
	if(keepGoing[i]!=0)
	  indexed=1;
	else
	  filesLeft--;
      }
      
      if(asso.size()==0){//this should never happen
	fprintf(stderr,"\t-> This should never happen\n");
	continue;
      }
  
      aMap::iterator it=--asso.end();
      lociPileup stop = it->first;
      if(0&&i==0)//print stop
	print_lociPileup(stop);
      if(keepGoing[i]==0){
	//	fprintf(stderr,"\t -> file[%lu] is done so skipping\n",i);
	continue;
      }
      else{
	//	fprintf(stderr,"\t-> File[%lu] will try to use stop\n",i);
	keepGoing[i] = obj[i]->readlines(asso,stop,0);
	//fprintf(stderr,"file[%lu] = lines actually read:%d\n",i,keepGoing[i]);
	if(0==keepGoing[i])
	    filesLeft--;
      }
      
    }
      //printMap(asso,obj[0].numFiles);
    if(asso.size()==0){
      filesLeft--;
    } 

    if(asso.size()==0)
      return NULL;
    if((int)asso.size()> chunkSize){
      mr::funkyPars *tmp = collapse(asso);
      asso.clear();
      return tmp;
    }
     //    exit(0);
    
    if(filesLeft==1){//LAST segment before EOF
      
      return (collapse(asso));

    }
   
  }
  return NULL;
}

pileups::~pileups(){
  if(nFiles==0)
    return;
  for(int i=0;i<filenames.size();i++)
    free(filenames[i]);

  delete [] keepGoing;
  for(int i =0;i<=(int) faiIndex.size();i++)
    free(chrnames[i]);
  delete[] chrnames;

  //cleanup faiIndex
  for(cMap::iterator it=faiIndex.begin();it!=faiIndex.end();++it)
    free(it->first);
  for(int i=0;i<nFiles;i++)
    obj[i]->close();

  if(type==2)
    delete [] tobj;
  else
    delete [] bobj;
  delete [] obj;
  
  
}



void pileups::print_lociPileup(const lociPileup &l){
  fprintf(stderr,"\r\t-> loci: (%d:%d)\t",l.chromo,l.position);
}

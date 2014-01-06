#include <vector>
#include <zlib.h>
#include <cstring>
#include <list>
#include <map>


#include "mrStruct.h"
#include "argStruct.h"

#include "soapReader.h"
#include "soapMaster.h"
#include "analysisFunction.h"

extern int SIG_COND;

#define POSITIONS ".pos"

using namespace soapworld;
/*
double sum(const double *ary,int len){
  double s =0;
  for(int i=0;i<len ; i++)
    s+=ary[i];
  //  printf("sum:%f\n",s);
  return s;
}
*/


void swap (int& first, int& second)
{
        int temp = first;
        first = second;
        second = temp;
}



void randomDownSample(suint *count,int depth,int minDepth,int nFiles){


  int totalArray[depth];
  int counter=0;
  suint newCount[nFiles*4];
  for(int i=0;i<nFiles;i++)
    for(int b=0;b<4;b++){
      newCount[i*4+b]=0;
      for(int j=0;j<count[i*4+b];j++){
	totalArray[counter]=i*4+b; 
	//	fprint(stderr,"%d\t",totalArray[counter]);
	counter++;
      }
    }
  if(counter!=depth){
    fprintf(stderr,"Depth to euqal to counter %d\t%d\t",depth,counter);   
    exit(0);

  }
  //permute totalArray
  for(int i=0;i<minDepth;i++) {
    int j = rand() % (depth);
    if(j>=depth||j<0){
      fprintf(stderr,"wrong random \n");
      exit(0);
    }  
    swap(totalArray[i],totalArray[j]);
  }


 for(int i=0;i<minDepth;i++)
   newCount[totalArray[i]]++;

 for(int i=0;i<nFiles;i++)
   for(int b=0;b<4;b++){
     count[i*4+b]=newCount[i*4+b];
     
   }
  //count=newCount;
  //fprintf(stderr,"\n");
  


}





void phat1(aMap asso,int nFiles,double eps,int minDepth,int maxDepth,FILE *oFile,double cutOff,int downsample,tmpList *suYeonsList){
  
  //algorithm goes on by a site on site basis
  
  double *pis = new double[nFiles];
  double *wis = new double[nFiles];
  
  //this is tricky, we shouldnt loop through the entire matrix, but just the number of lines that were actually read into the matrix
  for (aMap::iterator it=asso.begin();it!=asso.end();it++){
    //    fprintf(oFile,"%s\t%d\t",it->first.chromo,it->first.position);
    //part one
    //first lets get the sum of each nucleotide
    int bases[NUMBASES] = {0,0,0,0};
    for(int i=0;i<nFiles;i++)
      for(size_t j=0;j<NUMBASES;j++){
	bases[j] += it->second[i*NUMBASES+j];
	//	fprintf(oFile,"%d\t",bases[j]);
      }
    int depth = 0;
    for(size_t j=0;j<NUMBASES;j++)
      depth += bases[j];
    if(depth<minDepth||depth>maxDepth){
      //      fprintf(oFile,"-1\n");
      free(it->first.chromo);
      //fprintf(stderr,"herer\n");
      delete [] it->second; 
      continue;
    }
    if(downsample){
      srand(it->first.position);
      randomDownSample(it->second,depth,minDepth,nFiles);
      for(size_t j=0;j<NUMBASES;j++)
	bases[j]=0;
      for(int i=0;i<nFiles;i++)
	for(size_t j=0;j<NUMBASES;j++){
	  bases[j] += it->second[i*NUMBASES+j];
	//	fprintf(oFile,"%d\t",bases[j]);
	}
      depth = 0;
      for(size_t j=0;j<NUMBASES;j++)
	depth += bases[j];
      if(depth!=minDepth){
	fprintf(stderr,"depht is not equual to minDepth\n");
	continue;
      }

    }  
    //    fprintf(oFile,"depth=%d\t",depth);
    //now get the major/minor
    
    int major = 0;
    for(int i=1;i<4;i++)
      if (bases[i]>bases[major])
        major = i;
    //maj is now the major allele (most frequent)
    
    int temp=0;
    int minor= major;
    for(int i=0;i<4;i++){
      if(major==i) //we should check the against the major allele      
        continue;
      else if (bases[i]>temp){
        minor = i;
        temp=bases[i];
      }
    }
    if (minor==major){
      //      fprintf(oFile,"%d\t%d\t",major,minor);
      for(int i=0;0&&i<4;i++)
        fprintf(oFile,"%d\t",bases[i]);
      //      fprintf(oFile,"0\n");
      free(it->first.chromo);
      delete [] it->second;
      continue;
    }

    //    fprintf(oFile,"%d\t%d\t",major,minor);

    //if the site is variable
    
    for(int i=0;i<nFiles;i++){
      int ni = it->second[i*4+minor];
      int nt= it->second[i*4+minor] + it->second[i*4+major]; 
      if(nt==0){//if we dont have any reads for individual 'i'
        pis[i] = 0;
        wis[i] = 0;

        continue;
      }
      pis[i] = (ni-eps*nt)/(nt*(1-2*eps));
      wis[i] = 2.0*nt/(nt+1.0);
    }
    
    double tmp=0;
    for(int i=0;i<nFiles;i++)
      tmp+= (pis[i]*wis[i]);
    
    double phat = std::max(0.0,tmp/::angsd::sum<double>(wis,nFiles));
    if(phat>cutOff){
      aSite s;
      s.l = it->first;
      s.inf = it->second;
      s.phat = phat;
      s.minor = minor;
      s.major = major;
      suYeonsList->push_back(s);
      fprintf(oFile,"%s\t%d\t",it->first.chromo,it->first.position);//print pos
      fprintf(oFile,"%d\t%d\t",major,minor);
      for(int i=0;i<4;i++)
	fprintf(oFile,"%d\t",bases[i]);
      fprintf(oFile,"%f\n",phat);
    }else{
      free(it->first.chromo);
      delete [] it->second;
    }
  }
  
  
  delete [] pis;
  delete [] wis;
}







//funkyPars *collapse(tmpList *lis,size_t nInd,Matrix<double> ymat,int minHigh,int minCount,int emIter,int doMaf,int doAsso,double assoCutoff,double **errors,Matrix<double> covmat,int getCovar){
mr::funkyPars *collapse(tmpList *lis){
  if(lis->size()==0)
    return NULL;
  
  mr::funkyPars *r = mr::allocFunkyPars();

  r->numSites = lis->size();
  r->sites = new loci[r->numSites];
  r->phat = new double[r->numSites];
  r->minor = new char[r->numSites];
  r->major = new char[r->numSites];

  //  r->likes = NULL;
  r->counts  = new suint*[r->numSites];

  size_t pos =0;
  for(tmpList::iterator it=lis->begin();it!=lis->end();it++){
    r->sites[pos] = it->l;
    r->phat[pos] = it->phat;
    r->minor[pos] = it->minor;
    r->major[pos] = it->major;
    r->counts[pos++] = it->inf;
  }

  //printFunkyPars1(r);
  lis->clear();
  return r;
}




void soapMaster::printArg(FILE *argFile){
  fprintf(argFile,"%s:\n",__FILE__);
  fprintf(argFile,"-soap\t\n");
  fprintf(argFile,"-qs\t%d\n",qs);
  fprintf(argFile,"-maxHits\t%d\n",maxHits);
  fprintf(argFile,"-eps\t%f\n",eps);
  fprintf(argFile,"-strand\t%d\n",strand);
  fprintf(argFile,"-minDepth\t%d\n",minDepth);
  fprintf(argFile,"-maxDepth\t%d\n",maxDepth);
  fprintf(argFile,"-cutOff\t%f\n",cutOff);
  fprintf(argFile,"-downsample\t%d\n",downsample);
  fprintf(argFile,"\n");
}


soapMaster::soapMaster(argStruct *arguments){
  const char *outfiles = arguments->outfiles;
  char *soapFileName=0;
  soapFileName = angsd::getArg("-soap",soapFileName,arguments);
  if(strcmp(soapFileName,"-999")==0){
    printArg(stdout);
    exit(0);
  }
  std::vector<char *> filenames;
  filenames = angsd::getFilenames(soapFileName,-1);
  
  posiFilepointer = aio::openFile(outfiles,POSITIONS);
  
  //initialize all soapReader objects;
  obj = new soapReader[filenames.size()];
  for(size_t i=0;i<filenames.size();i++)
    obj[i].init(filenames[i]);
  keepGoing = new int[filenames.size()];
  for(int i=0;i<filenames.size();i++)
    keepGoing[i] =1;
  //memset(keepGoing,1,filenames.size());
  filesLeft = (int) filenames.size() +1;//we need one final iteration for cleaning up the buffers

  suYeonsList = new tmpList;
  //fprintf(stderr,"[%s] done\n",__FUNCTION__);
  
  nFiles = filenames.size();
  //
  qs = angsd::getArg("-qs",qs,arguments);
  maxHits = angsd::getArg("-maxHits",maxHits,arguments);
  eps = angsd::getArg("-eps",eps,arguments);
  strand = angsd::getArg("-strand",strand,arguments);
  minDepth = angsd::getArg("-minDepth",minDepth,arguments);
  maxDepth = angsd::getArg("-maxDepth",maxDepth,arguments);
  cutOff = angsd::getArg("-cutOff",cutOff,arguments);
  downsample = angsd::getArg("-downsample",downsample,arguments);

  printArg(arguments->argumentFile);
}

mr::funkyPars *soapMaster::fetch(int nLines,int chunkSize){
  

  while(SIG_COND&&(filesLeft)) {
    //    fprintf(stderr,"filesLeft: %d\n",filesLeft);
    if(filesLeft <0){
      fprintf(stderr,"\t corruption in [%s]\n",__FUNCTION__);
      exit(0);
    }
    
    aMap asso ;
    for(size_t i=0;i<nFiles;i++)
      obj[i].flush(asso);//flush buffered data from previous seq line.
    
    int indexed =0;
    for(size_t i=0;i<nFiles;i++){//loop through all files
      //	fprintf(stderr,"FILEID=%lu filesleft:%d\n",i,filesLeft);
      if(keepGoing[i]==0)//if where are done with file i skip it.
	continue;
      else if(indexed==0) {//if we havent read nLines from any file yet
	int lread = obj[i].doStuff_getAsso(asso,nLines,qs,strand,maxHits);
	indexed=1;
	if(0==(keepGoing[i]=obj[i].notEof())){
	  //  fprintf(stderr,"fileid=%lu decreminging \n",i);
	  filesLeft--;
	}
	if(lread==0)
	  continue;
      }
      if(asso.size()==0)
	continue;
      aMap::iterator it=--asso.end();
      loci stop = it->first;
      if(0&&i==0)//print some stuff to the screen
	print_loci(stop);
      if(keepGoing[i]==0){
	fprintf(stderr,"\t -> file[%lu] is done so skipping\n",i);
	continue;
      }
      else{
	int lread = obj[i].doStuff_useStop(asso,stop,qs,strand,maxHits);
	if(0==(keepGoing[i]=obj[i].notEof())){
	  //	    fprintf(stderr,"fileid=%lu decreminging obj[%lu].notEof=%d\n",i,i,obj[i].notEof());
	  filesLeft--;
	}
      }
    }
    //printMap(asso,obj[0].numFiles);
    if(asso.size()==0){
      if(filesLeft ==1){
	//fprintf(stderr,"decreminging due to zero indexed: \n");
	filesLeft--;
      }
    }// else
    //	print_loci((--asso.end())->first);
    
    if(asso.size()!=0){//only do this, we have new reads.
      loci loci1 = asso.begin()->first;
      loci loci2 = (--asso.end())->first;
      fprintf(stderr,"\t[phatinfo] nSites in region (%s,%d)-(%s,%d)=%lu ",loci1.chromo,loci1.position,loci2.chromo,loci2.position,asso.size());
      
      phat1(asso,obj[0].numFiles,eps,minDepth,maxDepth,posiFilepointer,cutOff,downsample,suYeonsList);
      fprintf(stderr,"(%lu sites above cutOff)\n",suYeonsList->size());
      if(suYeonsList->size()>chunkSize)
	return (collapse(suYeonsList));//,nFiles,ymat,minHigh,minCount,emIter,doMaf,doAsso,assoCutoff,errors,covmat   ,getCovar  ));
      asso.clear();
    }
    
  }
  if(suYeonsList->size()!=0)//LAST segment before EOF
    return (collapse(suYeonsList));//,nFiles,ymat,minHigh,minCount,emIter,doMaf,doAsso,assoCutoff,errors,covmat   ,getCovar));
  else{
    fprintf(stderr,"Wauv");
    return NULL;
  }
}

soapMaster::~soapMaster(){
  delete suYeonsList;
  delete [] keepGoing;
  for(int i=0;i<nFiles;i++)
    obj[i].close();
  delete [] obj;
  if(posiFilepointer) fclose(posiFilepointer);
}

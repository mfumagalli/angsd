
#include "multiReader.h"

mr::funkyPars *mr::allocFunkyPars(){
  //  fprintf(stderr,"allocing funkypars\n");
  mr::funkyPars *r = new mr::funkyPars;
  
  //  r->sites = NULL;
  r->counts = NULL;//soap
  r->likes = NULL;//glfs,sim,tgfl
  r->post = NULL;//beagle
  r->phat = NULL;//soap

  r->major = NULL;//tglf
  r->minor = NULL;//tglf
  r->depth = NULL;//pileups

  r->ref = NULL;//tglf
  r->anc= NULL;//tglf

  r->posi = NULL;//this might only be populated by one of them
  r->refId = -1;
  return r;
}


void mr::deallocFunkyPars(mr::funkyPars *p) {
  
  //cleanup begin  
  if(p->major!=NULL)
    delete[] p->major;
  if(p->minor!=NULL)
    delete[] p->minor;


  
  //  fprintf(stderr,"[%s] cleaning funkypars\n",__FUNCTION__);  

  if(p->depth!=NULL){
    for(int i=0;i<p->numSites;i++)
      delete [] p->depth[i];
    delete [] p->depth;
  }

  if(1){
    if(p->anc!=NULL)
      delete []  p->anc;
    if(p->ref!=NULL)
      delete [] p->ref;
  }


  delete [] p->posi;//<-new way
  
  delete p;

}




aHead *getHeadFromFai(const char *fname){
  std::vector<char *> chrs;
  std::vector<int> lengths;
  FILE *fp = aio::getFILE(fname,"r");
  char buf[1024];
  while(fgets(buf,1024,fp)){ 
    chrs.push_back(strdup(strtok(buf,"\t \n")));//<-strdup so don't clean here
    lengths.push_back(atoi(strtok(NULL,"\t \n")));
  }
  
  aHead *ret = new aHead;
  ret->l_text = strlen(fname);
  ret->text = new char[strlen(fname)+1];
  ret->text = strcpy(ret->text,fname);
  ret->n_ref = chrs.size();
  ret->l_name = new int[chrs.size()];
  ret->l_ref = new int [chrs.size()];
  ret->name = new char*[chrs.size()];
  for(size_t i=0;i<chrs.size();i++){
    ret->l_name[i] = strlen(chrs[i]);
    ret->l_ref[i] = lengths[i];
    //    ret->name[i] = chrs[i];
    ret->name[i] = new char[strlen(chrs[i])+1];
    strcpy(ret->name[i],chrs[i]);
  }
  return ret;
}

std::map<char *,int,ltstr> *buildRevTable(aHead *hd){
  std::map<char *,int,ltstr> *ret = new std::map<char *,int,ltstr>;
  for(int i=0;i<hd->n_ref;i++){
    ret->insert(std::pair<char *,int>(strdup(hd->name[i]),i));
  }
  for(std::map<char *,int,ltstr>::iterator it= ret->begin();0&&it!=ret->end();++it)
    fprintf(stderr,"%s %d\n",it->first,it->second);
  return ret;
}






void multiReader::printArg(FILE *argFile){
  fprintf(argFile,"%s:\n\n",__FILE__);
  fprintf(argFile,"-nLines=%d\n",nLines);
  fprintf(argFile,"-chunkSize=%d \n",chunkSize);
  fprintf(argFile,"-fai=%s\n",fainame);
  fprintf(argFile,"\n");
}
void multiReader::getOptions(argStruct *arguments){

  nLines=angsd::getArg("-nLines",nLines,arguments);
  chunkSize=angsd::getArg("-chunkSize",chunkSize,arguments);
  fainame = angsd::getArg("-fai",fainame,arguments);

  printArg(arguments->argumentFile);
}


multiReader::multiReader(int type_a,argStruct *arguments) {
  type = type_a;
  if(type>5||type==-1)
    return;
  fainame=NULL;
  nLines=1000;
  chunkSize = 10000;
  if(arguments->argc!=2)
    getOptions(arguments);
  
  sm=NULL;
  pl=NULL;
  mytglfs=NULL;

  // pg=pg_in;

  switch(type){
  case 0:{
    //    fprintf(stderr,"\t[%s] assuming soapfile\n",__FUNCTION__);
    sm = new soapMaster (arguments);
    break;
  }
  case 1:{
    //fprintf(stderr,"\t[%s] assuming glfv3 files \n",__FUNCTION__);
    //    pl = new pileups (faifile,lStart,lStop,filenames,type);
    pl = new pileups (arguments);
    break;
  }
  case 2:{
    //fprintf(stderr,"\t[%s] assuming glfv3 files \n",__FUNCTION__);
    //    pl = new pileups (faifile,lStart,lStop,filenames,type);
    pl = new pileups (arguments);
    break;
  }
  case 3:{
    //fprintf(stderr,"\t[%s] assuming tglf\n",__FUNCTION__);
    mytglfs = new tglfs;
    mytglfs->init(arguments);
    break;
  }
  case 4:{
    //  fprintf(stderr,"\t[%s] assuming simulation files (single pop)\n",__FUNCTION__);
    mysims = new reader_sim;
    mysims->init(arguments);
    break;
  }
  case 5:{
    //    fprintf(stderr,"\t[%s] assuming beagle gprobs input \n",__FUNCTION__);
    bglObj = new beagle_reader(arguments);
    //    bglObj->init(arguments);
    break;
  }
    
  default:{
    //    fprintf(stderr,"\t[%s] assuming what ?\n",__FUNCTION__);
    //pl = new pileups (faifile,lStart,lStop,filenames,type);
    break;
  }
  }
  if(arguments->inputtype==4||arguments->inputtype==7){
    if(fainame!=NULL){
      fprintf(stderr,"Bams or simfiles does not require -fai argument\n");
      exit(0);
    }

  }else{
    if(fainame==NULL){
      fprintf(stderr,"For non-bams or simfiles you should include -fai arguments\n");
      exit(0);
    }else{
      arguments->hd=getHeadFromFai(fainame);
      arguments->revMap = buildRevTable(arguments->hd);
    }
  }
  
}

multiReader::~multiReader(){
 
  switch(type){
  case 0:{
    delete sm;
    break;
  }
  case 3:{
    mytglfs->close();
    delete mytglfs;
    break;
  }
  case 4:{
    mysims->close();
    delete mysims;
    break;
  }    
  case 5:{
    delete bglObj;
      //    delete mysims;
    break;
  }    
  default:{
    break;
  }
  }
  if(type==2||type==3)
    delete pl;

}
mr::funkyPars *multiReader::fetch(){
  mr::funkyPars *fp = NULL;
  switch(type){
  case 0:{
    fp =  sm->fetch(nLines,chunkSize);
    break;
  }
  case 3:{
    fp =  mytglfs->fetch(nLines);
    break;
  }
  case 4:{
    fp = mysims->fetch(chunkSize); 
    break;
  }
  case 5:{
    fp = bglObj->fetch(chunkSize); 
    break;
  }
  default:{
    fp =  pl->fetch(nLines,chunkSize);
    break;
  }
  }

    
  return fp;

}

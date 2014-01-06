/*
  2 leaks
  1) when choosing region
  2) when doing strdup in indexing

 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <csignal>
#include <vector>
#include <stdint.h>

#include "bams.h"
#include "mUpPile.h"
#include "parseArgs_bambi.h"
#include "indexer.h"
#include "general.h"
#include "getFasta.h"
#include "analysisFunction.h"
#include "knetfile.h"
extern general **allMethods;
getFasta *gf=NULL;


extern int SIG_COND;

static const char *bam_nt16_rev_table2 = "=ACMGRSVTWYHKDBN";


void dalloc_bufReader(bufReader &ret){
  free(ret.bamfname);
  free(ret.it.off);//cleanup offsets if used
  bgzf_close(ret.fp);
  dalloc(ret.hd);
}



bufReader initBufReader2(const char*fname){
  bufReader ret;
  ret.bamfname = strdup(fname);
  ret.fp = openBAM(ret.bamfname);
  ret.hd = getHd(ret.fp);
  ret.isEOF =0;
  ret.it.from_first=1;//iterator, used when supplying regions
  ret.it.finished=0;//iterator, used when supplying regions
  ret.it.off = NULL;
  ret.it.dasIndex = NULL;
  return ret;
}


bufReader *initializeBufReaders2(const std::vector<char *> vec){
  bufReader *ret = new bufReader[vec.size()];
  
  for(size_t i =0;i<vec.size();i++)
    ret[i] = initBufReader2(vec[i]);

  //now all readers are inialized, lets validate the header is the same
  for(size_t i=1;i<vec.size();i++)
    if(compHeader(ret[0].hd,ret[i].hd))
      printErr();

  return ret;
}



void printAuxBuffered(uint8_t *s, uint8_t *sStop,kstring_t &str ) {
  //  fprintf(stderr,"\ncomp:%p vs %p\n",s,sStop);
  
  while (s < sStop) {
    uint8_t type;
    kputc('\t', &str);kputc(s[0], &str);kputc(s[1], &str); kputc(':', &str); 
    //    fprintf(stderr,"\t%c%c:",s[0],s[1]);
    s += 2; type = *s; ++s;
    //    fprintf(stderr,"\ntype=%c\n",type);//,(char)*s);
    //    kputc('\t', &str); kputsn((char*)key, 2, &str); kputc(':', &str);
    if (type == 'A') { kputsn("A:", 2, &str); kputc(*s, &str); ++s; }
    else if (type == 'C') { kputsn("i:", 2, &str); kputw(*s, &str); ++s; }
    else if (type == 'c') { kputsn("i:", 2, &str); kputw(*(int8_t*)s, &str); ++s; }
    else if (type == 'S') { kputsn("i:", 2, &str); kputw(*(uint16_t*)s, &str); s += 2; }
    else if (type == 's') { kputsn("i:", 2, &str); kputw(*(int16_t*)s, &str); s += 2; }
    else if (type == 'I') { kputsn("i:", 2, &str); kputuw(*(uint32_t*)s, &str); s += 4; }
    else if (type == 'i') { kputsn("i:", 2, &str); kputw(*(int32_t*)s, &str); s += 4; }
    else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
    else if (type == 'd') { ksprintf(&str, "d:%lg", *(double*)s); s += 8; }
    else if (type == 'Z' || type == 'H') { kputc(type, &str); kputc(':', &str); while (*s) kputc(*s++, &str); ++s; }
    else if (type == 'B') {
      uint8_t sub_type = *(s++);
      int32_t n;
      memcpy(&n, s, 4);
      s += 4; // no point to the start of the array
      kputc(type, &str); kputc(':', &str); kputc(sub_type, &str); // write the typing
      for (int i = 0; i < n; ++i) {
	kputc(',', &str);
	if ('c' == sub_type || 'c' == sub_type) { kputw(*(int8_t*)s, &str); ++s; }
	else if ('C' == sub_type) { kputw(*(uint8_t*)s, &str); ++s; }
	else if ('s' == sub_type) { kputw(*(int16_t*)s, &str); s += 2; }
	else if ('S' == sub_type) { kputw(*(uint16_t*)s, &str); s += 2; }
	else if ('i' == sub_type) { kputw(*(int32_t*)s, &str); s += 4; }
	else if ('I' == sub_type) { kputuw(*(uint32_t*)s, &str); s += 4; }
	else if ('f' == sub_type) { ksprintf(&str, "%g", *(float*)s); s += 4; }
      }
    }
  }
  //  fprintf(stderr,"done\n");
}





void printReadBuffered(aRead &rd,aHead *hd,kstring_t &str) {
   str.l = 0;
   
   if(bam_validate1(hd,rd)==0){
     fprintf(stderr,"problems validateing\n");
     exit(0);
   }
   kputsn((char *)rd.vDat,rd.l_qname-1,&str);kputc('\t', &str);
   kputw((int)rd.flag_nc>>16, &str); kputc('\t', &str); 
   
   
   if(rd.refID==-1)//unmatched read
     kputc('*', &str);
   else
     kputs(hd->name[rd.refID] , &str);
   kputc('\t', &str); 
   
   kputw(rd.pos+1, &str);   kputc('\t', &str); 
   kputw(rd.mapQ, &str);kputc('\t', &str); 


   int nCigs = rd.nCig;

   if(nCigs==0)
     kputc('*', &str);// if no cigars
   else{
     for (int i = 0; i < nCigs; ++i) {//print cigars
       uint32_t *cigs =getCig(&rd);
       kputw(cigs[i]>>BAM_CIGAR_SHIFT, &str);
       kputc("MIDNSHP"[cigs[i]&BAM_CIGAR_MASK], &str);
     }
   }
   kputc('\t', &str); 
   
   if(rd.next_refID==-1)
     kputc('*', &str);// if no cigars     
   else if(rd.refID==rd.next_refID)
     kputc('=', &str);
   else
     kputs(hd->name[rd.next_refID] , &str);
   kputc('\t', &str); 

   kputw(rd.next_pos+1, &str);   kputc('\t', &str); 
   kputw(rd.tlen, &str);   kputc('\t', &str); 


   //start seq
   char *seq = (char *)getSeq(&rd);
   for(int i=0;i<rd.l_seq;i++)
     kputc(bam_nt16_rev_table2[bam1_seqi(seq, i)], &str); 
   
   kputc('\t', &str); 

   char *quals =(char *)getQuals(&rd);
   for(int i=0;i<rd.l_seq;i++)
     kputc(quals[i]+33, &str); 
   

   //below is taken directly from samtools,(not to steal, to preserve ordering etc, all credits go where credit is due)
   //from aux start to the last memadrs in chunk
   printAuxBuffered(getAuxStart(&rd),rd.vDat+rd.block_size,str);
   kputc('\n', &str); 
}

//this function will print a samfile from the bamfile
int motherView(bufReader *rd,int nFiles,std::vector<regs>regions) {
  aRead b;
  b.vDat=new uint8_t[RLEN];
  kstring_t str;  str.s=NULL; str.l=str.m=0;
  
  if(regions.size()==0) {//print all
    int block_size;
    while(SIG_COND && bgzf_read(rd[0].fp,&block_size,sizeof(int))){
      getAlign(rd[0].fp,block_size,b);
      printReadBuffered(b,rd[0].hd,str);
      fprintf(stdout,"%s",str.s);
    }
  }else {
    for(int i=0;i<(int)regions.size();i++){
      int tmpRef = regions[i].refID;
      int tmpStart = regions[i].start;
      int tmpStop = regions[i].stop;
      
      getOffsets(rd[0].bamfname,rd[0].hd,rd[0].it,tmpRef,tmpStart,tmpStop);
      int ret =0;
      while(SIG_COND){
	ret = bam_iter_read(rd[0].fp, &rd[0].it, b);
	if(ret<0)
	  break;
	printReadBuffered(b,rd[0].hd,str);
	fprintf(stdout,"%s",str.s);
      }
      free(rd[0].it.off);//the offsets
    }
    free(str.s);
    delete [] b.vDat;
  }
  return 0;
}





std::vector<char*> getRegions(const char * name){
  if(!aio::fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t\n\r";
  std::vector<char*> ret;
  FILE *fp =aio::getFILE(name,"r");
  char buffer[aio::fsize(name)+1];
  if(aio::fsize(name)!=fread(buffer,1,aio::fsize(name),fp))
    fprintf(stderr,"[%s] Problems reading %lu from: %s\n",__FUNCTION__,aio::fsize(name),name);
  buffer[aio::fsize(name)]='\0';
  
  char *tok = strtok(buffer,delims);

  while(tok!=NULL){
    if(tok[0]!='#'){
      ret.push_back(strdup(tok));
    }
    tok = strtok(NULL,delims);
  }

  fprintf(stderr,"\t-> From regionsfile: %s we read %lu\n",name,ret.size());
  if(fp) fclose(fp);
  return ret;
}




void printChunky2(const chunky* chk,FILE *fp,char *refStr) {
  //  fprintf(stderr,"[%s] nsites=%d region=(%d,%d) itrReg=(%d,%d)\n",__FUNCTION__,chk->nSites,chk->refPos[0],chk->refPos[chk->nSites-1],chk->regStart,chk->regStop);
  if(chk->refPos[0]>chk->regStop){
    fprintf(stderr,"\t->Problems with stuff\n");
    exit(0);
  }
  int refId = chk->refId;
  for(int s=0;s<chk->nSites;s++) {
    if(chk->refPos[s]<chk->regStart || chk->refPos[s]>chk->regStop-1 ){
      for(int i=0;i<chk->nSamples;i++)
	dalloc_node(chk->nd[s][i]);
      delete [] chk->nd[s];
      continue;
    }
    fprintf(fp,"%s\t%d",refStr,chk->refPos[s]+1);     
    if(gf->ref!=NULL){
      if(refId!=gf->ref->curChr)
	gf->loadChr(gf->ref,refStr,refId);
      if(gf->ref->seqs!=NULL){
	fprintf(fp,"\t%c",gf->ref->seqs[chk->refPos[s]]);
      }
    }
    for(int i=0;i<chk->nSamples;i++) {

      //      fprintf(stderr,"seqlen[%d,%d]=%lu\t",s,i,chk->nd[s][i].seq->l);
      if(chk->nd[s][i].seq.l!=0){
	fprintf(fp,"\t%d\t",chk->nd[s][i].depth);
	for(size_t l=0;l<chk->nd[s][i].seq.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].seq.s[l]);
	fprintf(fp,"\t");
		
	for(size_t l=0;l<chk->nd[s][i].qs.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].qs.s[l]);
	//	fprintf(fp,"\t");

      }else
	fprintf(fp,"\t0\t*\t*");
      dalloc_node(chk->nd[s][i]);
    }
    //    fprintf(stderr,"\n");
    fprintf(fp,"\n");
    delete [] chk->nd[s];
  }
  delete [] chk->nd;
  delete [] chk->refPos;
  delete chk;
} 

void (*func)(void *) = NULL;

void printReg(FILE *fp,std::vector<regs> &regions){
  fprintf(fp,"-------------\n");
  fprintf(fp,"regions.size()=%lu\n",regions.size());
  for(size_t i=0;i<regions.size();i++)
    fprintf(fp,"reg[%zu]= %d %d %d\n",i,regions[i].refID,regions[i].start,regions[i].stop);
  fprintf(fp,"-------------\n");
}

char *download_from_remote(const char *url)
{
  const int buf_size = 1 * 1024 * 1024;
  char *fn;
  FILE *fp=NULL;
  uint8_t *buf;
  knetFile *fp_remote;
  int l;
  l = strlen(url);
  for (fn = (char*)url + l - 1; fn >= url; --fn)
    if (*fn == '/') break;
  ++fn; // fn now points to the file name
  if(aio::fexists(fn))
    return fn;
  fprintf(stderr, "attempting to download the remote index file: %s\n",url);
  extern std::vector <char *> dumpedFiles;
  dumpedFiles.push_back(strdup(fn));
  fp_remote = knet_open(url, "r");
  if (fp_remote == 0) {
    fprintf(stderr, "Fail to open remote file:%s\n",url);
    exit(0);
  }
  if ((fp = fopen(fn, "wb")) == 0) {
    fprintf(stderr, "Fail to save remote file:%s in working dir\n",url);
    exit(0);
    knet_close(fp_remote);
  }
  buf = (uint8_t*)calloc(buf_size, 1);
  while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
    fwrite(buf, 1, l, fp);
  free(buf);
  if(fp) fclose(fp);
  knet_close(fp_remote);
  return fn;
}



void modNames(bufReader *rd,int nFiles){
  for(int i=0;i<nFiles;i++){
    char *fn = rd[i].bamfname;
    if (strstr(fn, "ftp://") == fn || strstr(fn, "http://") == fn){
      //check if bai file exist on server
      char *fnidx =(char*) calloc(strlen(fn) + 5, 1);
      strcat(strcpy(fnidx, fn), ".bai");
      char* tmp =strdup(download_from_remote(fnidx));
      free(rd[i].bamfname);
      free(fnidx);
      rd[i].bamfname=tmp;
    }else{
      char *fnidx =(char*) calloc(strlen(fn) + 5, 1);
      strcat(strcpy(fnidx, fn), ".bai");
      free(rd[i].bamfname);
      rd[i].bamfname=fnidx;
    }


  }

}



int bammer_main(args *pars){
  gf=(getFasta *) allMethods[1];
  func = pars->callback;
  std::vector<char *> nams;

  //read bamfiles
  if(pars->inputfiles)
    nams = angsd::getFilenames(pars->inputfiles,pars->nInd);
  else
    nams.push_back(strdup(pars->inputfile));//if only one file just push bamfile
  
  bufReader *rd = initializeBufReaders2(nams);
 
  //read regions

  std::vector<char *> regionsRaw;
  if(pars->regionfile)
    regionsRaw = getRegions(pars->regionfile);
  else if(pars->region!=NULL){
    regionsRaw.push_back(strdup(pars->region));//if only one file just push bamfile
  }
  std::vector<regs> regions;
  for(size_t i=0;i<regionsRaw.size();i++){
    regs tmpRegs;
    if(parse_region(regionsRaw[i],*rd[0].hd,tmpRegs.refID,tmpRegs.start,tmpRegs.stop)<0||tmpRegs.stop<tmpRegs.start){
      fprintf(stderr,"[%s] problems with indexing: %s\n",__FUNCTION__,regionsRaw[i]);
      exit(0);
    }else
      regions.push_back(tmpRegs);
  }

  //printReg(stderr,regions)
  //each filereader contains a filename for a bam, if we try to read remote files, then we need to download the bai file, and update the filename
  if(regions.size()!=0)
    modNames(rd,nams.size());


  extern int maxThreads;

  if(pars->jobtype==2)
    uppile(pars->show,maxThreads,rd,pars->nLines,nams.size(),regions);
  else if(pars->jobtype==1)    	
    while(motherView(rd,nams.size(),regions));
  else{
    fprintf(stderr,"nothing to do? what should program do?\n");
  }

  //cleanup stuff
  for(int i=0;i<(int)nams.size();i++){
    free(nams[i]);
    dalloc_bufReader(rd[i]);

  }
  for(size_t i=0;i<regions.size();i++)
    free(regionsRaw[i]);
  delete [] rd;
  free(pars->inputfiles);free(pars->inputfile);free(pars->region);
  delete pars;
  return 0;
}

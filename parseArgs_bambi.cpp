#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "general.h"
#include "parseArgs_bambi.h"
#include "shared.h"
#include "analysisFunction.h"


// below is default samtools parameters
int uniqueOnly = 0;
int only_proper_pairs = 1;
int remove_bads = 1;
int minMapQ =0;
int minQ = MINQ;
int adjustMapQ =0;
int baq =0;
void printArg(FILE *argFile,args *ret){
  fprintf(argFile,"---------------\n%s: bam reader:\n",__FILE__);
  fprintf(argFile,"\t-bam\t\t%s\tSupply a file list of BAMfiles\n",ret->inputfile);
  fprintf(argFile,"\t-i\t\t%s\tSupply a single BAMfile\n",ret->inputfile);
  fprintf(argFile,"\t-r\t\t%s\tSupply a single region in commandline (see examples below)\n",ret->region);
  fprintf(argFile,"\t-rf\t\t%s\tSupply multiple regions in a file (see examples below)\n",ret->regionfile);
  fprintf(argFile,"\t-remove_bads\t%d\tDiscard \'bad\' reads, (flag >=255) \n",remove_bads);
  fprintf(argFile,"\t-nInd\t\t%d\tOnly use first nInd from the filelist from the -bam argument\n",ret->nInd);
  fprintf(argFile,"\t-nLines\t\t%d\tRead nLines from files at a time\n",ret->nLines);
  //  fprintf(argFile,"\t-type\t\t%d\n",ret->type);
  fprintf(argFile,"\t-uniqueOnly\t%d\tDiscards reads that doesn't map uniquely\n",uniqueOnly);
  fprintf(argFile,"\t-show\t\t%d\tMimic 'samtools mpileup' also supply -ref fasta for printing reference column\n",ret->show);
  fprintf(argFile,"\t-minMapQ\t%d\tDiscard reads with mapping quality below\n",minMapQ);
  fprintf(argFile,"\t-minQ\t%d\tDiscard reads with mapping quality below\n",minQ);
  fprintf(argFile,"\t-only_proper_pairs\t%d\tOnly use reads where the mate could be mapped\n",only_proper_pairs);
  fprintf(argFile,"\t-C\t\t%d\tadjust mapQ for excessive mismatches (as SAMtools), supply -ref\n",adjustMapQ);
  fprintf(argFile,"\t-baq\t\t%d\tadjust qscores around indels (as SAMtools), supply -ref\n",baq);
  fprintf(argFile,"\n");
  fprintf(argFile,"Examples for region specification:\n");
  fprintf(argFile,"\t\tchr:\t\tUse entire chromosome: chr\n");
  fprintf(argFile,"\t\tchr:start-\tUse region from start to end of chr\n");
  fprintf(argFile,"\t\tchr:-stop\tUse region from beginning of chromosome: chr to stop\n");
  fprintf(argFile,"\t\tchr:start-stop\tUse region from start to stop from chromosome: chr\n");
  fprintf(argFile,"\t\tchr:site\tUse single site on chromosome: chr\n");
}

args *allocArgs (){
  args *ret = new args;
  ret->inputfile = NULL;
  ret->inputfiles = NULL;
  ret->region = NULL;
  ret->regionfile= NULL;
  ret->callback= NULL;
  ret->show = 0;
  ret->nInd = 0;
  ret->type = 1;
  ret->nLines = 50;
  ret->inputfile=ret->inputfiles=ret->region=ret->regionfile=NULL;
  ret->jobtype = 2;
  return ret;
}

//read program parameters
args *getArgsBambi(argStruct *arguments){
  args *ret = allocArgs();
  ret->inputfiles=angsd::getArg("-bam",ret->inputfiles,arguments);
  ret->inputfiles=angsd::getArg("-b",ret->inputfiles,arguments);
  ret->inputfile=angsd::getArg("-i",ret->inputfile,arguments);
  if(ret->inputfiles==NULL && ret->inputfile==NULL)
    return NULL;
  if(ret->inputfiles!=NULL&&strcmp(ret->inputfiles,"-999")==0){//DRAGON fix this
    ret->inputfiles=NULL;
    printArg(stdout,ret);
    exit(0);
  }
  remove_bads = angsd::getArg("-remove_bads",remove_bads,arguments);
  uniqueOnly = angsd::getArg("-uniqueOnly",uniqueOnly,arguments);
  only_proper_pairs =angsd::getArg("-only_proper_pairs",only_proper_pairs,arguments);
  minMapQ = angsd::getArg("-minMapQ",minMapQ,arguments);
  minQ = angsd::getArg("-minQ",minQ,arguments);
  adjustMapQ = angsd::getArg("-C",adjustMapQ,arguments);
  baq = angsd::getArg("-baq",baq,arguments);
  ret->region =angsd::getArg("-r",ret->region,arguments);
  ret->regionfile = angsd::getArg("-rf",ret->regionfile,arguments);
  ret->nInd = angsd::getArg("-nInd",ret->nInd,arguments);
  ret->nLines = angsd::getArg("-nLines",ret->nLines,arguments);
  ret->type = angsd::getArg("-type",ret->type,arguments);
  ret->show = angsd::getArg("-show",ret->show,arguments);
  //  minQ = angsd::getArg("-minQ",minQ,arguments);
  if(ret->inputfile && ret->inputfile[0]=='-'){
    fprintf(stderr,"[Error] Problem understanding inputfilename: %s\n",ret->inputfile);
    return NULL;
  }
  char *tmp = NULL;
  tmp = angsd::getArg("-ref",tmp,arguments);
  if(tmp!=NULL && adjustMapQ!=0){
    fprintf(stderr,"Must also supply -ref for adjusting the mapping quality\n");
    exit(0);
  }
  if(tmp!=NULL&&baq!=0){
    fprintf(stderr,"Must also supply -ref for adjusting base qualities (baq)\n");
    exit(0);
  }
  free(tmp);
  printArg(arguments->argumentFile,ret);
  return ret;
}

void bamInfo(FILE *fp){
  args *ret = allocArgs();
  printArg(fp,ret);
}

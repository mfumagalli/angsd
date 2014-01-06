/*
  copyright (gpl)   19 sep 2012.
  Thorfinn Sand Korneliussen thorfinn@binf.ku.dk
  Anders albrechtsen albrecht@binf.ku.dk 
  0.12 more object oriented version of 0.03
  0.2 implemented bambi the native bamreader and GL estimator
  0.3 first real public, cleancing of program
  0.4 much is now better concerning screen output
  0.4485 Version which are released together with submit of application note.
  0.500 first major fixup will begin to remove loci
  0.501 continuing above
  0.534 many minor bugfixes
  use wiki for updated list
*/

#define VERSION 0.570
#define ARGS ".arg"

#include<iostream>//for printing time
#include<cstring> //for cstring functions
#include<cstdlib> //for exit()
#include<cstdio> //for fprintf
#include <signal.h>//for catching ctrl+c, allow threads to finish
#include <libgen.h>//for checking if output dir exists 'dirname'

#include "shared.h"
#include "multiReader.h"
#include "mrStruct.h"
#include "parseArgs_bambi.h"
#include "bammer_main.h"
#include "analysisFunction.h"

extern std::vector <char *> dumpedFiles;//a vector of outputfilenames from shared.cpp

int SIG_COND =1;//if we catch signal then quit program nicely
int VERBOSE =1;

extern pthread_mutex_t mUpPile_mutex;//mutex for not deleting mUppile data untill end of data

void printTime(FILE *fp){
 time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}


//this function is called from within the bamreader
void callBack_bambi(fcb *fff){

  if(fff==NULL){
    //    fprintf(stderr,"SEnding NULL this is a killswitch");
    selector(NULL);//<-send NULL which acts as a killswitch
  }else{
    funkyPars *fp = allocFunkyPars();
    fp->for_callback=fff;
    fp->refId = fp->for_callback->refId;
    selector(fp);
  }
}



void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)");
  VERBOSE=0;
  SIG_COND=0;
  pthread_mutex_unlock(&mUpPile_mutex);
}

 //we are threading so we want make a nice signal handler for ctrl+c
void catchkill(){
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

}

//print nice info
void printProgInfo(FILE *fp){

  fprintf(fp,"\t-> Please use the website \"http://www.popgen.dk/angsd\" as reference\n");
  fprintf(fp,"\t-> Use -nThreads or -P for number of threads allocated to the program\n\n"); 
 fprintf(fp,"Overview of methods:\n");
  fprintf(fp,"\t-GL\t\tEstimate genotype likelihoods\n");
  fprintf(fp,"\t-doCounts\tCalculate various counts statistics\n");
  fprintf(fp,"\t-doAsso\t\tPerform association study\n");
  fprintf(fp,"\t-doMaf\t\tEstimate allele frequencies\n");
  fprintf(fp,"\t-doError\tEstimate the type specific error rates\n");
  fprintf(fp,"\t-doAnsError\tEstimate the errorrate based on perfect fastas\n");
  fprintf(fp,"\t-doHWE\t\tEst inbreedning per site\n");
  fprintf(fp,"\t-doGeno\t\tCall genotypes\n");
  fprintf(fp,"\t-doFasta\tGenerate a fasta for a BAM file\n");
  fprintf(fp,"\t-doAbbababa\tPerform an ABBA-BABA test\n");
  fprintf(fp,"\t-sites\t\tAnalyse specific sites (can force major/minor)\n");
  fprintf(fp,"\t-doSaf\t\tEstimate the SFS and/or neutrality tests genotype calling\n");
  fprintf(fp,"\n\tBelow are options that can be usefull\n");
  fprintf(fp,"\t-bam\t\tOptions relating to bam reading\n\t-doMajorMinor\tInfer the major/minor using different approaches\n");  
  fprintf(fp,"\t-ref/-anc\tRead reference or ancestral genome\n");
  fprintf(fp,"\tmany others\n\n");
 
  fprintf(fp,"For information of specific options type: \n\t./angsd METHODNAME eg \n\t\t./angsd -GL\n\t\t./angsd -doMaf\n\t\t./angsd -doAsso etc\n");
  fprintf(fp,"Examples:\n\tEstimate MAF for bam files in 'list'\n\t\t\'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1\'\n");
}




void parseArgStruct(argStruct *arguments){
  //validate that there are no duplicate parameters
  for(int i=0;i<arguments->argc-1;i++)
    if(strcmp("-realSFS",arguments->argv[i])==0){
	fprintf(stderr,"-realSFS is now called -doSaf\n");
	exit(0);
    }  
  
  for(int i=0;i<arguments->argc-1;i++)
    if(arguments->argv[i][0]=='-') {
      for(int ii=i+1;ii<arguments->argc;ii++){
	if(arguments->argv[ii][0]=='-'){
	  if(0==strcmp(arguments->argv[i],arguments->argv[ii])){
	    fprintf(stderr,"Duplicate parameter: %s  supplied, will exit\n",arguments->argv[i]);
	    exit(0);
	  }
	}
      }
    }

  for(int i=1;i<arguments->argc;i++){
    if(arguments->usedArgs[i]==0){
      fprintf(stderr,"%d argument \t%s is unknown will exit\n",i,arguments->argv[i]);
      fflush(stderr);
      exit(0);
    }
  }
  if(arguments->hd==NULL){
    fprintf(stderr,"\t-> Error: You must supply a -fai file such that we know the sizes of the genomes (version .4486)\n");
     exit(0);
  }


}




 void getNind(const char *name,argStruct *arguments){
   if(!aio::fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(0);
  }
  int nInd = -1;
  nInd = angsd::getArg("-nInd",nInd,arguments);
  std::vector<char*> delme = angsd::getFilenames(name,nInd);
  arguments->nInd=delme.size();
  
  for(size_t i=0;i<delme.size();i++)
    free(delme[i]);
 }


argStruct *setArgStruct(int argc,char **argv) { 
  argStruct *arguments = new argStruct;
  arguments->outfiles = NULL;
  arguments->hd = NULL;
  arguments->argc=argc;
  arguments->argv=argv;
  arguments->usedArgs= new int[argc+1];//well here we allocate one more than needed, this is only used in the ./angsd -beagle version
  for(int i=0;i<argc;i++)
    arguments->usedArgs[i]=0;
  
  arguments->inputtype=-1;
  
  if(0&&argc==2)//if help should be printed
    return arguments;

  int argPos = 1;
  int extra= argc==2?1:0;
  //  fprintf(stderr,"argpos=%d cond=%d\n",argPos,argc+extra);
  while(argPos <argc+extra) { 
    //fprintf(stderr,"argv[%d]=%s\n",argPos,argv[argPos]);
    if (strcmp(argv[argPos],"-samglf")==0)
      arguments->inputtype=1;
    else if (strcmp(argv[argPos],"-samglfclean")==0)
      arguments->inputtype=2;
    else if (strcmp(argv[argPos],"-sim1")==0)
      arguments->inputtype=4;
    else if (strcmp(argv[argPos],"-beagle")==0){
      arguments->inputtype=5;
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
    }
     else if (strcmp(argv[argPos],"-soap")==0){
       arguments->inputtype=0;
       arguments->usedArgs[argPos]=1;
       if(argc==2)
	 break;
     } else if (strcmp(argv[argPos],"-tglf")==0){
       arguments->inputtype=3;
       arguments->usedArgs[argPos]=1;
       arguments->usedArgs[argPos+1]=1;
     }     else if(strcmp(argv[argPos],"-out")==0){
       arguments->outfiles = argv[argPos+1];
       arguments->usedArgs[argPos]=1;
       arguments->usedArgs[argPos+1]=1;
     //check if folder exists, otherwise program crashes
       char *dirNam2 = strdup(arguments->outfiles);
       char *dirNam=dirname(dirNam2);
       //       fprintf(stderr,"dirname: %s strlen(dirNam):%zu\n",dirNam,strlen(dirNam));
       if(strlen(dirNam)>1 &&!aio::fexists(dirNam)){
	 fprintf(stderr,"\t Folder: \'%s\' doesn't exist, please create\n",dirNam);
	 exit(0);
       }
       free(dirNam2);
     }
     else if(strcmp(argv[argPos],"-bam")==0||strcmp(argv[argPos],"-i")==0||strcmp(argv[argPos],"-b")==0 ){
       if(argc==2)
	 break;
       if(strcmp(argv[argPos],"-bam")==0||strcmp(argv[argPos],"-b")==0){
	 getNind(argv[argPos+1],arguments);
	 if(arguments->nInd==0){
	   fprintf(stderr,"BAM list looks empty? will exit\n");
	   exit(0);
	 }
	 char buftmp[1024];//fix this if user gives a bam instead of a bam filelist
	 FILE *fp = fopen(argv[argPos+1],"r");
	 int delme= fread(buftmp,1,1024,fp);
	 buftmp[1023] = '\0';
	 char *tmp = strtok(buftmp,"\t \n");
	 arguments->hd = getHd_andClose(tmp);
	 if(fp) fclose(fp);
	 
       }else{
	 arguments->nInd = 1;
	 arguments->hd = getHd_andClose(argv[argPos+1]);
       }
       std::map<char *,int,ltstr> *buildRevTable(aHead *hd);
       arguments->revMap = buildRevTable(arguments->hd);
       arguments->inputtype=7;
       arguments->usedArgs[argPos]=1;
       arguments->usedArgs[argPos+1]=1;
       
     }
     argPos+=2;

   }
  //  fprintf(stderr,"inputtype=%d\n",arguments->inputtype);
  if(arguments->inputtype==-1 && argc!=2){
    fprintf(stderr,"\t-> Error: supply inputfiles:[-i afile.bam -bam bamlist, -samglf, -beagle, -samglfclean, -tglf]\n");
    exit(0);
  }

   return arguments;
 }

template <typename T>
T* get(T* p,size_t nItems){
  T * r = new T[nItems];
  //fprintf(stderr,"copyinging: nItems=%lu with persize=%lu in totalbytes=%lu\n",nItems,sizeof(T),nItems*sizeof(T));
  memcpy(r,p,sizeof(T)*nItems);
  return r;
}



void splitFunkyPars(mr::funkyPars *fp,std::map<char*,int,ltstr> *revMap){
  std::map<int,int> tmp;//will contain the number of sites on different chrs;
  std::map<char*,int,ltstr>::iterator it;
  std::map<int,int>::iterator cnts;
  for(int i=0;i<fp->numSites;i++){
    it=revMap->find(fp->sites[i].chromo);
    if(it==revMap->end()){
      fprintf(stderr,"Problem finding chromosome: %s in index\n",fp->sites[i].chromo);
      exit(0);
    }
    cnts = tmp.find(it->second);
    if(cnts==tmp.end())
      tmp.insert(std::pair<int, int> (it->second,1));
    else
      (cnts->second)++;
  }
  for(cnts=tmp.begin();cnts!=tmp.end();++cnts)
     fprintf(stderr,"%d %d\n",cnts->first,cnts->second);

  int posi =0;
  for(cnts=tmp.begin();cnts!=tmp.end();++cnts){
    funkyPars *f = allocFunkyPars();
    f->numSites = cnts->second;
    f->posi = new int[f->numSites];
    f->refId = cnts->first;
    f->nInd = fp->nInd;
    if(fp->counts!=NULL)
      f->counts = new suint*[f->numSites];
    if(fp->likes!=NULL)
      f->likes = new double*[f->numSites];
    if(fp->post!=NULL)
      f->post = new double*[f->numSites];
    if(fp->major!=NULL){
      f->major = get<char>(fp->major+posi,f->numSites);
      f->minor = get<char>(fp->minor+posi,f->numSites);
    }
    if(fp->ref!=NULL)
      f->ref = get<char>(fp->ref+posi,f->numSites);
    if(fp->anc!=NULL)
      f->anc = get<char>(fp->anc+posi,f->numSites);
    if(fp->phat!=NULL)
      f->phat = get<double>(fp->phat+posi,f->numSites);
    for(int i=0;i<f->numSites;i++){
      f->posi[i] = fp->sites[posi].position;
      if(fp->counts!=NULL)
	f->counts[i] = fp->counts[posi];
      if(fp->likes!=NULL)
	f->likes[i] = fp->likes[posi];
      if(fp->post!=NULL)
	f->post[i] = fp->post[posi];
      posi++;
    }
    void waiter(int);
    waiter(f->refId);
    selector(f);
  }
  delete fp;

}

void sendToSelector(mr::funkyPars *mfp,std::map<char *,int,ltstr> *revMap,int inputtype){

  if(inputtype==4){
    //simulated data, only ancestral likelihds and posi and refId
    //this is the most simple case, since we know that all data belongs to the same chromosome
    funkyPars *fp = allocFunkyPars();
    fp->nInd = mfp->nInd;
    fp->numSites = mfp->numSites;
    fp->refId = mfp->refId;
    fp->likes = mfp->likes;
    fp->anc = mfp->anc;
    fp->posi = mfp->posi;
    delete mfp;
    selector(fp);
  }else{
    //first check if the chunk belong to the same chromosome
    if(strcmp(mfp->sites[0].chromo,mfp->sites[mfp->numSites-1].chromo)==0){
      funkyPars *fp = allocFunkyPars();
      fp->nInd = mfp->nInd;
      fp->numSites = mfp->numSites;
      fp->likes = mfp->likes;
      fp->anc = mfp->anc;
      fp->ref = mfp->ref;
      
      fp->major = mfp->major;
      fp->minor = mfp->minor;
      fp->phat = mfp->phat;
      fp->post = mfp->post;
      fp->likes = mfp->likes;
      fp->counts = mfp->counts;
      
      std::map<char *,int,ltstr>::iterator it = revMap->find(mfp->sites[0].chromo);
      if(it==revMap->end()){
	fprintf(stderr,"Problem finding chromosome: %s\n",mfp->sites[0].chromo);
	exit(0);
      }else
	fp->refId = it->second;
      
      fp->posi = new int[fp->numSites];
      for(int i=0;i<mfp->numSites;i++)
	fp->posi[i] = mfp->sites[i].position;
      
      selector(fp);
    }else{
      //this is the most ugly part. We have more than one chromosome, so we need to split accordingly
      splitFunkyPars(mfp,revMap);
      
    }
  }
  
}


 int main(int argc, char** argv){
   //   fprintf(stderr,"DEV VERSION might not work thorfinn feb 1 2013 sizeof(funkypars):%zu\n",sizeof(funkyPars));
   //no arguments supplied -> print info
   if(argc==1||(argc==2&&(strcmp(argv[1],"--version")==0||strcmp(argv[1],"--help")==0))){//if haven't been supplied with arguments, load default,print, and exit
     printProgInfo(stderr);
     return 0;
   }
   
   //print time
   clock_t t=clock();
   time_t t2=time(NULL);

   //intialize our signal handler for ctrl+c
   catchkill();

   //print arguments supplied
   fprintf(stderr,"Command:\n");
   for(int i=0;i<argc;i++)
     fprintf(stderr,"%s ",argv[i]);
   fprintf(stderr,"\n\t-> angsd version: %.3f\t build(%s %s)\n",VERSION,__DATE__,__TIME__); 

   

   argStruct *arguments = NULL;
   arguments = setArgStruct(argc,argv);
   //   fprintf(stderr,"here: %d\n",arguments->inputtype);
   arguments->argumentFile=stderr;
   //catch the informative case where we want to view options associated with a command
   if(argc==2){
     //if we dont use bamfiles
     //fprintf(stderr,"here: %d\n",arguments->inputtype);
     multiReader mr(arguments->inputtype,arguments);
     FILE *tempFile;
     fprintf(stderr,"\t-> Analysis helpbox/synopsis information:\n");
     args *bamPars=getArgsBambi(arguments);
     init(arguments);//program dies here after printing info
   }
   
   //check that user supplied an output prefix for files
   if(arguments->outfiles ==NULL){
     fprintf(stderr,"\t-> No \'-out\' argument given, output files will be called \'angsdput\'\n");
     arguments->outfiles = strdup("angsdput");
     //     fprintf(stderr,"\t-> Must supply -out\n");
     //return 0;
   }

   //print arguments into logfile
   FILE *argFp = NULL;
   argFp = aio::openFile(arguments->outfiles,ARGS);
   arguments->argumentFile=argFp;
   multiReader mr(arguments->inputtype,arguments);

   args *bamPars =NULL;
   if(arguments->inputtype==7){
     if(!(bamPars=getArgsBambi(arguments))){
       bamInfo(stderr);
       return 0;
     }
   }
   
   for(int i=0;i<argc;i++)
     fprintf(argFp,"%s ",argv[i]);
   fprintf(argFp,"\n\n");
   fprintf(argFp,"\n\t-> angsd version: %.3f\t build(%s %s)\n",VERSION,__DATE__,__TIME__); 
   fprintf(argFp,"\t->");
   printTime(argFp);  


   init(arguments);
   parseArgStruct(arguments);//validate that all arguments has indeed been parsed in the instantiation 'init' above

   //Below is main loop which will run until nomore data
   if(arguments->inputtype!=7){
     fprintf(stderr,"\t-> multireader\n");
    while(SIG_COND) {

      mr::funkyPars *tmp = mr.fetch(); //YES MISTER FETCH

      if(tmp==NULL)
	break;

      sendToSelector(tmp,arguments->revMap,arguments->inputtype);

    }
  }else 
     bammer_main(bamPars);
     //     bammer_main(bamPars,arguments);
  
   //now we have no more data lets wait and cleanup

  fprintf(stderr,"\t-> Done reading data waiting for calculations to finish\n");

  destroy();

  
  //printout the filenames generated
  fprintf(stderr,"\t-> Output filenames:\n");
  for(int i=0;i<(int)dumpedFiles.size();i++){
    fprintf(stderr,"\t\t->\"%s\"\n",dumpedFiles[i]);
    fprintf(argFp,"\t\t->\"%s\"\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  // fprintf(stderr,"\n");
  fprintf(argFp,"\n");
  // fprintf(stderr,"\t");
  printTime(stderr);
  
  delete []   arguments->usedArgs;
  dalloc(arguments->hd);
  for(  std::map<char *,int,ltstr> ::iterator it=arguments->revMap->begin();it!=arguments->revMap->end();it++){
    free(it->first);
    //    fprintf(stderr,"freeing\n");
  }
  delete arguments->revMap;
  delete arguments;
  //print out nice messages
  
  fprintf(stderr,"\t-> Arguments and parameters for all analysis are located in .arg file\n");
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(argFp, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(argFp, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  if(argFp) fclose(argFp);
  
  return 0;
}

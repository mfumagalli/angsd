
#include "pileup.h"
#include "mrStruct.h"
using namespace pileworld;

int glfClass::total=0;

inline
double calcLikeRatio(int a){
  return (pow(10.0,-a/10.0));
}

//to speed up calculations, we will generate a likelookup table for all values 0-255.
//The index will then be the value from one of the columns in the glf file

double likeLookup[256];


void fileErrorPileup(const char* file){
  fprintf(stderr,"\t-> Problems opening file: %s\t right directory?\n",file);
}



bool operator <= (const lociPileup& first,const  lociPileup& second) {

    if(first.chromo==second.chromo)
      return first.position<=second.position;
    else
      return first.chromo<=second.chromo;
}

bool operator > (const lociPileup& first,const  lociPileup& second) {
    return (!(first <= second));
}

bool operator == (const lociPileup& first,const  lociPileup& second) {

    if(first.chromo==second.chromo)
      return first.position==second.position;
    else
      return first.chromo==second.chromo;
}




datumP allocDatum(int numInds){
  datumP data;

  data.datumDepth= new int[numInds]; //]allocIntArray(numInds);
  for(int i=0;i<numInds;i++)
    data.datumDepth[i] = 0;
  data.lk = new double[10*numInds];
  for(int i=0;i<numInds*10;i++)
    data.lk[i] = 0;
  //  data.lk=allocDoubleArray(10*numInds);//allocUnsignedCharArray(3*numInds);
  data.phat =0.0;
  data.ancestral = 0;
  data.major=-1;
  data.minor=-1;
  //  data.depth =0;
  return data;
}



void generate_likeLookup(){
  //fprintf(stderr,"generating likelookup\n");
  for (short int i=0;i<=255;i++)
    likeLookup[i] =  log(calcLikeRatio(i));
}






void txtInput::init(int num,const char *fname,lociPileup start, lociPileup stop,const cMap& c){

  if(num==0)
    generate_likeLookup();
  cmap =c;
  total++;
  filename = fname;
  id=num;
  buffered=0;
  ultrastart = start;
  ultrastop = stop;
  gz=Z_NULL;//only used for making valgrind perfect :)
  inGlf.open(filename,std::ios::in);
  //  fprintf(stderr,"\t->id:%d of %d \t initializingfile : %s\n",id,total,fname);
  if(start==(lociPileup){0,0} )
    return;

  //  fprintf(stderr,"[%s] is fastforwarding to start:\n",fname);
  //added in 0.8 this is necessary for going through the glf file before the actual readline

  char buffer[LENS];
  char *chr;
  int pos;

  if(!inGlf){
    //if file cant be opened, dont panic, just continue with next file
    fileErrorPileup(filename);
  }


  while(inGlf.getline(buffer,LENS)){

    if(buffer==NULL)
      continue;
    //    printf("%s\n",buffer);
    chr = strtok(buffer," \t");
    pos = atoi(strtok(NULL," \t"));

    char *tmp_ref=strtok(NULL," \t");
    char ref = tmp_ref[0];


    //fixed case of reference beeing unmappable
    if(ref=='N'||ref=='*')//DANIELLE
      continue;

    cMap::const_iterator cit = cmap.find(chr);
    if(cit==cmap.end()){
      printf("chr: %s doesn't exist in faiIndex will assume file is DONE \n",chr);
      exit(0);
    }



    lociPileup tmp_lociPileup = {cit->second,pos}; //this is the key for our associative array

    //if we havent reached the start yet
    if (ultrastart>tmp_lociPileup)
      continue;
    //print("key is larger :",tmp_lociPileup,stop);

      //the key is larger than the stop so put it in buffer and exit loop
    last.chromo = tmp_lociPileup.chromo;
    last.position = tmp_lociPileup.position;
      //      print("glf[BUFFER] inputting this in to the buffer",id,last);

    depthInd = atoi(strtok(NULL," \t"));
	strtok(NULL," \t");
	strtok(NULL," \t");
    for(int i=0;i<10;i++)
      likeratios[i] = likeLookup[atoi(strtok(NULL," \t"))];
    buffered = 1;
    break;
  }

  //return (!inGlf.eof()||buffered);
  
}


void binInput::init(int num,const char *fname,lociPileup start, lociPileup stop,const cMap& c){
  //  fprintf(stderr,"[%s] start\n",__FUNCTION__);
  if(num==0)
    generate_likeLookup();
  cmap =c;
  total++;
  filename = fname;
  id=num;
  buffered=0;
  ultrastart = start;
  ultrastop = stop;
  gz=Z_NULL;
  if(Z_NULL==(gz= gzopen(fname,"r"))){
    fprintf(stderr,"\t Problems opening: %s\n",fname);
    exit(0);
  }

  const char magic[] = "GLF\3";
  char buf[64];
  gzread(gz,buf,4);
  if(memcmp(buf,magic,4)!=0){
    fprintf(stderr,"\t-> Problems reading magic number, file not a glfv3 file:\n%s\nGLF\3\n",buf);
    exit(0);
  }
  
  int32_t hlen;
  gzread(gz,&hlen,sizeof(int32_t));//get length of extra information
  gzread(gz,buf,hlen); //read extra information
  readHeader =1;
  //  fprintf(stderr,"[%s]:%d end readheader=%d\n",__FUNCTION__,id,readHeader);
}




void glfClass::flush(aMap &asso){
  //  fprintf(stderr,"buffered is: %d\n",buffered);
  if(buffered){
    //    printf("a buffered element\n");
    //we got a buffered element but the key doesnt exist yet;

    if(last>ultrastop){
      //      fprintf(stderr,"we got a lockdown");
      return ;

    }

    aMap::iterator it = asso.find(last); //lets find a iterator for it!
    //    print(1"glf[BUFFER] Will insert: \t",id,last);
    double *the_line = NULL;
    int *the_depth = NULL;
    if(it==asso.end()){
      datumP dats= allocDatum(total);
      the_line =dats.lk;
      the_depth=dats.datumDepth;
      asso.insert(std::make_pair(last,dats));
    }else{
      the_line = it->second.lk;
      the_depth=it->second.datumDepth;
    }
    //copy buffered values to the datum
    /*from 0.981 we are keeping all likeratios*/
    for(int i=0;i<10;i++)
      the_line[10*id+i] = likeratios[i];
    //    fprintf(stderr,"%p\n",it->second.datumDepth);
    the_depth[id] = depthInd;
    depthInd =0;
    buffered=0;

  }



}

int txtInput::readlines (aMap &asso,lociPileup stop,int nLines){
  
  //  printf("\t-> id=%d\t%s is reading glfs\n",id,filename);
  //print(stop);
  aMap::iterator it;
  char buffer[LENS];
  char *chr;
  int pos;
  //  double lk[10];//this will contain all 10 likeratios from the glf file AA,..,TT


  int linesRead = 0;
  while(inGlf.getline(buffer,LENS)) {

    chr = strtok(buffer," \t");
    pos = atoi(strtok(NULL," \t"));
    
    char *tmp_ref=strtok(NULL," \t");
    char ref = tmp_ref[0];


    //fixed case of reference beeing unmappable
    if(ref=='N'||ref=='*')//DANIELLE
      continue;

    cMap::const_iterator cit = cmap.find(chr);
    if(cit==cmap.end()){
      printf("chr: %s doesn't exist in faiIndex will assume that file is DONE\n",chr);
      return 0;
    }



    lociPileup tmp_lociPileup = {cit->second,pos}; //this is the key for our associative array
    
    if(tmp_lociPileup>ultrastop){
      //fprintf(stderr,"[%s]we got a lockdown\n",filename);
      //  fprintf(stderr,"[%d]we got a lockdown\n",id);
      return 0;
    }



    int myDepth = atoi(strtok(NULL," \t"));
    strtok(NULL," \t");
    strtok(NULL," \t");


    if (nLines==0){//if this is the first file dont buffer anything
      if(tmp_lociPileup>stop){
	//print("key is larger :",tmp_lociPileup,stop);
	//the key is larger than the stop so put it in buffer and exit loop
	
	last.chromo = tmp_lociPileup.chromo;
	last.position = tmp_lociPileup.position;
      //      print("glf[BUFFER] inputting this in to the buffer",id,last);
	
	for(int i=0;i<10;i++)
	  likeratios[i] = likeLookup[atoi(strtok(NULL," \t"))];//+retVal;
	depthInd = myDepth;
	buffered = 1;
	break;
      }
    }
    it = asso.find(tmp_lociPileup); //lest find a iterator for it!


    double *the_line= NULL;
    int *the_depth = NULL;
    if(it==asso.end()){
      datumP dats= allocDatum(total);
      //      fprintf(stderr,"in alloc part %p\n",dats.datumDepth);
      the_line =dats.lk;
      the_depth = dats.datumDepth;
      asso.insert(std::make_pair(tmp_lociPileup,dats));
    }else{
      the_line = it->second.lk;
      the_depth = it->second.datumDepth;
    }
    //fprintf(stderr,"%p\n",it->second.datumDepth);
    for(int i=0;i<10;i++)
      the_line[id*10+i] = likeLookup[atoi(strtok(NULL," \t"))];//offsetting by id*10
    //fprintf(stderr,"%d\n",id);
    the_depth[id] = myDepth;
    linesRead++;
    if(nLines!=0){
      if(linesRead>nLines)
	break;
    }
  }

  return (buffered+linesRead);
}


void printLik(FILE *fp,uint8_t lik[10]){
  for(int i=0;i<9;i++)
    fprintf(fp,"%u\t",lik[i]);
  fprintf(fp,"%u\n",lik[9]);
}




int binInput::readlines (aMap &asso,lociPileup stop,int nLines){
  //  FILE *of=stdout;
  //  printf("\t-> id=%d\t%s is reading glfs readHeader=%d\n",id,filename,readHeader);
  //print(stop);
  
  int linesRead = 0;

  uint32_t refLength;
  int notEOF =1;
  while(notEOF){//read until end of file
    
    if(readHeader){
      //      fprintf(stderr,"[%d]-> READING CHR HEADER\n",id);
      int32_t chrLen;
      gzread(gz,&chrLen,sizeof(int32_t));//read the length of the chromosomename inc \0
      curPos =1;//set startpos
      gzread(gz,chrName,chrLen);//read chromoname
      gzread(gz,&refLength,sizeof(uint32_t));//read reference length
      readHeader =0;//now we read the header
      cMap::const_iterator cit = cmap.find(chrName);
      if(cit==cmap.end()){
	printf("chr: %s doesn't exist in faiIndex, assuming end of file\n",chrName);
	return 0;
      }
      myChrId=cit->second;
    }
    while(readHeader==0) {

      uint8_t rt;
      uint32_t offset,min_depth;
      uint8_t rms;


      if(0==gzread(gz,&rt,sizeof(uint8_t))){//if EOF
	notEOF =0;//this will break 2 loops
	break;
      }
      int rtype= rt>>4;
      if(rtype!=0){
	gzread(gz,&offset,sizeof(uint32_t));
	curPos  += offset;
	gzread(gz,&min_depth,sizeof(uint32_t));
	gzread(gz,&rms,sizeof(uint8_t));
      }if(rtype==1){
	linesRead++; //we have a site 
	uint8_t lk[10];
	gzread(gz,lk,10*sizeof(uint8_t));
	/*
	//fprintf(of,"%d\t%s\t%d\t",nLines,chromo,position);
	fprintf(of,"%s\t%d\t",chrName,curPos);
	fprintf(of,"%c\t","XACMGRSVTWYHKDBN"[(rt)&0xf]);
	fprintf(of,"%u\t",min_depth<<8>>8);//print readdeptm
	fprintf(of,"%u\t",rms);//print rms
	fprintf(of,"%u\t",min_depth>>24); //print mapqual
	printLik(of,lk);
	*/


	//now populate

	lociPileup tmp_lociPileup = {myChrId,curPos}; //this is the key for our associative array
	
	if(tmp_lociPileup>ultrastop){
	  //      fprintf(stderr,"we got a lockdown");
	  return 0;
	}
	if (nLines==0){//if we use "stop" solely to break
	  
	  if(tmp_lociPileup>stop){
	    //we need to buffer this; and exit function
	    last.chromo = tmp_lociPileup.chromo;
	    last.position = tmp_lociPileup.position;
	    //      print("glf[BUFFER] inputting this in to the buffer",id,last);
	    
	    for(int i=0;i<10;i++)
	      likeratios[i] = likeLookup[lk[i]];//+retVal;
	    depthInd = min_depth<<8>>8;
	    buffered = 1;
	    return linesRead;
	  }
	}
	// we didnt need to buffer so lets put it in the dataset


	double *the_line= NULL;
	int *the_depth = NULL;

	int doNewVersion =1;
	if(doNewVersion){
	  aMap::iterator it = asso.lower_bound(tmp_lociPileup);
	  if(it != asso.end() && !(asso.key_comp()(tmp_lociPileup, it->first))){
	    // key already exists
	    the_line = it->second.lk;
	    the_depth = it->second.datumDepth;
	  }else {
	    // the key does not exist in the map
	    datumP dats= allocDatum(total);
	    //      fprintf(stderr,"in alloc part %p\n",dats.datumDepth);
	    the_line =dats.lk;
	    the_depth = dats.datumDepth;
	    asso.insert(it ,aMap::value_type(tmp_lociPileup, dats));
	  }
	}else{
	  
	  aMap::iterator it = asso.find(tmp_lociPileup); //lest find a iterator for it!
	  if(it==asso.end()){
	    datumP dats= allocDatum(total);
	    //      fprintf(stderr,"in alloc part %p\n",dats.datumDepth);
	    the_line =dats.lk;
	    the_depth = dats.datumDepth;
	    asso.insert(std::make_pair(tmp_lociPileup,dats));
	  }else{
	    the_line = it->second.lk;
	    the_depth = it->second.datumDepth;
	  }
	}
	//fprintf(stderr,"%p\n",it->second.datumDepth);
	for(int i=0;i<10;i++)
	  the_line[id*10+i] = likeLookup[lk[i]];//offsetting by id*10
	//fprintf(stderr,"%d\n",id);
	the_depth[id] = min_depth<<8>>8;
	//	fprintf(stderr,"[%s]:%d =linesread=%d IN\n",__FUNCTION__,id,linesRead);
	if(nLines!=0){
	  if(linesRead>nLines)
	    return linesRead;
	}

      }else if(rtype==2){
	uint8_t lkHom1;
	int16_t indelLen1;
	int16_t indelLen2;
	char indel1[abs(indelLen1)];
	char indel2[abs(indelLen1)];

	gzread(gz,&lkHom1,sizeof(uint8_t));
	gzread(gz,&lkHom1,sizeof(uint8_t));
	gzread(gz,&lkHom1,sizeof(uint8_t));
	gzread(gz,&indelLen1,sizeof(int16_t));
	gzread(gz,&indelLen2,sizeof(int16_t));
	gzread(gz,&indel1,abs(indelLen1));
	gzread(gz,&indel2,abs(indelLen2));
      }else if(rtype==0){
	readHeader =1;
      }
    }
  }
  fprintf(stderr,"[%s] =linesread=%d\n",__FUNCTION__,linesRead);
  
  return buffered+linesRead;

}



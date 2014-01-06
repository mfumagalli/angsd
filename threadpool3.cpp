
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <cstring>
#include <cstdlib>

//#define NTHREADS    5
//#define NPROBS 500

#define __printErr(b) (({fprintf(b,"\t-> Problems: in func:%s at line:%d in file:%s \n",__FUNCTION__,__LINE__,__FILE__);fflush(b);  exit(0);}))



typedef struct{
  int                 STOPTHREADS;//variable used for exiting threads
  
  //cond and mutex for signaling to all threads to start
  pthread_cond_t      *cond_startThreads ;
  pthread_mutex_t     *mutex_startThreads;
  
  //cond and mutex for signalling that threads are done used by maincalculation)
  pthread_cond_t      *cond_done;
  pthread_mutex_t     *mutex_done;

  //counter/mutex for selecting which data_for_threads[NPROBS] to work on:{0,...,(NPROBS-1)}
  pthread_mutex_t     *mutex_which_elem;  
  int which_array_elem;

  //mutex lock for making sure we are only running on iter at a time
  pthread_mutex_t     *mutex_iter;  
  int currentlyRunning; //0 or 1

  //
  //pthread_spinlock_t *spinner;

  int nDone;//counter for checking how many threads are done 

  //2below added 22 feb 2012
  int type;
  int nThreads;

  FILE *fp;
  
//utility function for printout line nr etc in case of error
  pthread_t *threadid;
  void (*callback)(void *,int);

  void *funpara;
  int length;
  int nthread;
  void *thread_args;//all of them, in a inner* array
  int *splitpoints;
} threadpool;


typedef struct{
  int threadNum;
  threadpool *tp;
}inner;




void *threadfunc_spin(void *parm) {
  inner *args = (inner *) parm;
  int i = args->threadNum;
  fprintf(args->tp->fp,"i=%d is init\n",i);
  threadpool *tp = args->tp;
  int  rc;
  
  while (0==tp->STOPTHREADS) {//main loop
    if((rc = pthread_mutex_lock(tp->mutex_startThreads)))
      __printErr(tp->fp);
    
    if((rc = pthread_cond_wait(tp->cond_startThreads, tp->mutex_startThreads)))
      __printErr(tp->fp);
    if((rc = pthread_mutex_unlock(tp->mutex_startThreads)))
      __printErr(tp->fp);
    if(tp->STOPTHREADS==1)
      break;
    while(1){//this is the loop that will be performed when _iter is triggered
      if((rc=pthread_mutex_lock(tp->mutex_which_elem)))
	__printErr(tp->fp);
      int which_index = tp->which_array_elem;
      tp->which_array_elem++;
      if((rc=pthread_mutex_unlock(tp->mutex_which_elem)))
	__printErr(tp->fp);

      if(which_index>=tp->length){
	if((rc=pthread_mutex_lock(tp->mutex_done)))
	  __printErr(tp->fp);
	if((++(tp->nDone))==tp->nthread)
	  if((rc=pthread_cond_signal(tp->cond_done)))
	    __printErr(tp->fp);

	if((rc=pthread_mutex_unlock(tp->mutex_done)))
	  __printErr(tp->fp);
	break;
      }
      fprintf(tp->fp,"\tthread\t%d\t wokring on index=%d\n",(int)i,which_index);
      fflush(tp->fp);
      tp->callback(tp->funpara,which_index);
      //      doCalc(which_index);
    }
  }
  //  fprintf(tp->fp,"exiting thread[%d]\n",i);
  //fflush(tp->fp);
  return NULL;
}

void *threadfunc_seq(void *parm) {
  inner *args = (inner *) parm;
  int threadindex = args->threadNum;
  fprintf(stderr,"[%s] i=%d is init\n",__FUNCTION__,threadindex);
  threadpool *tp = args->tp;
  int  rc;
  
  while (0==tp->STOPTHREADS) {//main loop
    if((rc = pthread_mutex_lock(tp->mutex_startThreads)))
      __printErr(tp->fp);
    
    if((rc = pthread_cond_wait(tp->cond_startThreads, tp->mutex_startThreads)))
      __printErr(tp->fp);
    if((rc = pthread_mutex_unlock(tp->mutex_startThreads)))
      __printErr(tp->fp);
    if(tp->STOPTHREADS==1)
      break;
    while(1){//this is the loop that will be performed when _iter is triggered
      for(int i=tp->splitpoints[threadindex];i<tp->splitpoints[threadindex+1];i++)
	tp->callback(tp->funpara,i);

      //now jobs done increment the nDone
      //      if((rc=pthread_spin_lock(tp->spinner)))
      //	__printErr(tp->fp);
      tp->nDone++;
      //     if((rc=pthread_spin_unlock(tp->spinner)))
      //	__printErr(tp->fp);
      
      if(tp->nDone==tp->nthread){
	//	fprintf(stderr,"[%s] signalling (all threads done with calcs)\n",__FUNCTION__);
	if((rc=pthread_mutex_lock(tp->mutex_done)))
	  __printErr(tp->fp);
	if((rc=pthread_cond_signal(tp->cond_done)))
	  __printErr(tp->fp);
	if((rc=pthread_mutex_unlock(tp->mutex_done)))
	  __printErr(tp->fp);
      }
      break;
    }
  
  }
  delete args;
  //  fprintf(tp->fp,"exiting thread[%d]\n",threadindex);
  //fflush(tp->fp);
  return NULL;
}



threadpool *init_threadpool(int nthread,void (*callback)(void*,int),int nprobs,int type){
  threadpool *ret = new threadpool;
  ret->STOPTHREADS =0;
  
  ret->cond_startThreads = new pthread_cond_t;
  ret->cond_done = new pthread_cond_t;

  ret->mutex_startThreads = new pthread_mutex_t;
  ret->mutex_done = new pthread_mutex_t;
  ret->mutex_which_elem = new pthread_mutex_t;
  ret->mutex_iter = new pthread_mutex_t;
  
  //spinner
  //  ret->spinner = new pthread_spinlock_t;

  pthread_cond_init(ret->cond_startThreads,NULL);
  pthread_cond_init(ret->cond_done,NULL);
  
  pthread_mutex_init(ret->mutex_startThreads,NULL);
  pthread_mutex_init(ret->mutex_done,NULL);
  pthread_mutex_init(ret->mutex_which_elem,NULL);
  pthread_mutex_init(ret->mutex_iter,NULL);  
  // pthread_spin_init(ret->spinner,0);
  ret->fp = stderr;

  ret->nDone =0;
  ret->callback=callback;
  ret->which_array_elem =0;
  ret->nthread = nthread;
  //initialize threads
  ret->threadid = new pthread_t[nthread];
  ret->thread_args = new inner[nthread];

  //if do sequential threading calculate splitpoints
  ret->splitpoints=NULL;
  if(type==1){
    ret->splitpoints = new int[nthread+1];
    ret->splitpoints[0] = 0;
    for(int i=1;i<nthread;i++)
      ret->splitpoints[i] = ret->splitpoints[i-1] + nprobs/nthread;
    ret->splitpoints[nthread]=nprobs;
#if 0
      for(int i=0;i<nthread+1;i++)
	fprintf(stderr,"i=%d\t",ret->splitpoints[i]);
      fprintf(stderr,"\n");
#endif
  }
  


  int rc;


  for(int i=0; i<ret->nthread; ++i) {
    inner *parm = new inner;
    parm->threadNum =i;
    parm->tp =ret;
    if(type ==0 ){
      if((rc = pthread_create(&ret->threadid[i], NULL, threadfunc_spin,(void *) parm)))
	__printErr(ret->fp);
    }
    else if(type==1)
      if((rc = pthread_create(&ret->threadid[i], NULL, threadfunc_seq,(void *) parm)))
	__printErr(ret->fp);
  }
  ret->type = type;
  ret->nThreads = nthread;
  sleep(1);//stupid way to serialize threads;
  return ret;
}

void threadpool_iter(threadpool *tp,void *funpara,int length){
  pthread_mutex_lock(tp->mutex_iter);
  int rc;
   
  if((rc=pthread_mutex_lock(tp->mutex_startThreads)))
    __printErr(tp->fp);
  tp->length=length;
  tp->funpara=funpara;
  /* The condition has occured. Set the flag and wake up any waiters */
  //conditionMet = 1;
  //  fprintf(tp->fp,"\nWake up all waiters...\n");
  // fflush(tp->fp);
  if((rc = pthread_cond_broadcast(tp->cond_startThreads)))
    __printErr(tp->fp);
  
  if((rc = pthread_mutex_unlock(tp->mutex_startThreads)))
    __printErr(tp->fp);
  pthread_mutex_unlock(tp->mutex_iter);
}

void wait(threadpool *tp){
  int rc;
  if((rc=pthread_mutex_lock(tp->mutex_done)))
    __printErr(tp->fp);
  while(tp->nDone<tp->nthread){
    if((rc=pthread_cond_wait(tp->cond_done,tp->mutex_done)))
      __printErr(tp->fp);
  }
  if((rc=pthread_mutex_unlock(tp->mutex_done)))
    __printErr(tp->fp);

  tp->nDone = tp->which_array_elem=0;
  //  fprintf(stderr,"done waiting\n");
  //fflush(stderr);
}


void wait_kill(threadpool *tp){
  int rc;
  
  if((rc=pthread_mutex_lock(tp->mutex_done)))
    __printErr(tp->fp);
  fprintf(stderr,"waitkill prewhile: tp->nDone:%d\n\n",tp->nDone);
  fflush(stderr);
 
  while(tp->nDone<tp->nthread){
    fprintf(stderr,"waitkill inwhile: tp->nDone:%d\n\n",tp->nDone);
    fflush(stderr);
    if((rc=pthread_cond_wait(tp->cond_done,tp->mutex_done)))
      __printErr(tp->fp);
  }
  if((rc=pthread_mutex_unlock(tp->mutex_done)))
    __printErr(tp->fp);
  fprintf(stderr,"done waiting will kill now\n");
  fflush(stderr);
  tp->STOPTHREADS = 1;
  
  pthread_mutex_lock(tp->mutex_iter);
   
  if((rc=pthread_mutex_lock(tp->mutex_startThreads)))
    __printErr(tp->fp);

  /* The condition has occured. Set the flag and wake up any waiters */
  //conditionMet = 1;
  // fprintf(tp->fp,"\nSignalling all threads to die...\n");
  fflush(tp->fp);
  if((rc = pthread_cond_broadcast(tp->cond_startThreads)))
    __printErr(tp->fp);
  
  if((rc = pthread_mutex_unlock(tp->mutex_startThreads)))
    __printErr(tp->fp);
  pthread_mutex_unlock(tp->mutex_iter);
  delete [] tp->splitpoints;
  
  for (int i=0; i<tp->nThreads; ++i) {
    if((rc = pthread_join(tp->threadid[i], NULL)))
      __printErr(tp->fp);
  }

  delete tp->cond_startThreads;
  delete tp->mutex_done;
  delete tp->mutex_startThreads;
  //delete tp->spinner;
  delete tp->cond_done;
  delete tp->mutex_which_elem;
  delete tp->mutex_iter;
  delete [] tp->threadid;
  inner *tmptmp =(inner *) tp->thread_args;
  delete [] tmptmp;
  //delete [] tp->splitpoints;
delete tp;
  
}




//code below is an example

#ifdef WITH_MAIN

typedef struct{
  int l;
  int *ary;

}workdata;



int nCal=0;
void doCalc(void *funpara,int index){
  fprintf(stderr,"index=%d\n",index);
  nCal++;
  int nrep=4000;
  workdata *data_for_threads = (workdata *) funpara;
  workdata wd = data_for_threads[index];
  for(int i=0;i<nrep;i++)
    for(int j=0;j<wd.l;j++)
      wd.ary[j] += wd.ary[j]+index;
  
}


void initData(workdata*data_for_threads, int NPROBS){
  int arrayLen =8096;
  for(int i=0;i<NPROBS;i++){
    data_for_threads[i].l=arrayLen;
    data_for_threads[i].ary = new int[arrayLen];
    memset(data_for_threads[i].ary,0,sizeof(int)*arrayLen);
  }
  
}
void collectData(workdata *data_for_threads,int NPROBS){
  for(int i=0;i<NPROBS;i++)
    delete [] data_for_threads[i].ary;
}

FILE *fp = stdout;


int main(){
  int NPROBS = 1000;
  int NTHREADS = 5;
  int type =1;
  threadpool *tp = init_threadpool(NTHREADS,doCalc,NPROBS,type);

  int ntimes =5;
  int rc;


  workdata data_for_threads[NPROBS];
  while(ntimes){
    fprintf(fp,"ntimes=",ntimes);
    initData(data_for_threads,NPROBS);
    threadpool_iter(tp,data_for_threads,NPROBS);
    wait(tp);
    collectData(data_for_threads,NPROBS);
    ntimes--;

  }
  tp->STOPTHREADS=1;
  threadpool_iter(tp,data_for_threads,NPROBS);
  //  wait(tp);
  for (int i=0; i<NTHREADS; ++i) {
    if((rc = pthread_join(tp->threadid[i], NULL)))
      __printErr(tp->fp);
  }

  fprintf(stderr,"ncal is=%d\n",nCal);
}

#endif

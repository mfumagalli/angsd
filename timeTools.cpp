/*
  stupid little function to check creation times of files


 */
#include <sys/stat.h>

//checks that newer is newer than older
int isNewer(const char *newer,const char *older){
  struct stat one;
  struct stat two;
  stat(newer, &one );
  stat(older, &two );
  
  return one.st_mtime>=two.st_mtime;
}


#ifdef __WITH_MAIN__
#include <cstdio>

int main(int argc,char **argv){
  if(argc!=3){
    fprintf(stderr,"supply 2 filenames\n");
    return 0;
  }
  fprintf(stderr,"is %s modified after: %s = %d\n",argv[1],argv[2],isNewer(argv[1],argv[2]));
  

  return 0;
}
#endif

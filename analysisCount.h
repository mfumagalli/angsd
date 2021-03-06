
class countCls : public general{
private:
  const char* postfix1; //.pos
  const char* postfix2;//.counts;
  const char* postfix3;//.qs
  const char* postfix4;//.depthSample
  const char* postfix5;//.depthGlobal
  char *oFiles;
  int dumpCounts;
  int doCounts;
  int doQsDist;//0=nothing 1=overall 

  char *minQfile;
  angsd::Matrix<double> minQmat;


  kstring_t bpos;
  kstring_t bbin;

  gzFile oFileCountsBin;
  gzFile oFileCountsPos;
 
  size_t *qsDist;
  int nInd;
  int minInd;

  int trim;
  int setMaxDepth;
  int setMinDepth;
  //from depth class
  int doDepth;
  int maxDepth;

  size_t **depthCount;
  size_t *globCount;


public:

  
  //none optional stuff
  suint **countNucs(const chunkyT *chk,int trim,int *keepSites);  
  countCls(const char *outfiles,argStruct *arguments,int inputtype);
  ~countCls();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};

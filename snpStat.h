

class snpStat:public general{
  gzFile outfileZ;
public:
  int doSnpStat;
  
  snpStat(const char *outfiles,argStruct *arguments,int inputtype);
  ~snpStat();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};

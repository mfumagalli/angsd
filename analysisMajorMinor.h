

class majorminor:public general{
  int doMajorMinor;
public:
  void printArg(FILE *argFile);
  void run(funkyPars *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  majorminor(const char *outfiles,argStruct *arguments,int inputtype);
  void getOptions(argStruct *arguments);
  ~majorminor();
};

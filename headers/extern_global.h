#ifndef GLOBAL
#define GLOBAL

extern int status;
extern int verb;
extern int rank;
extern int size;
extern int LikeNum;
extern double NU;
extern char workdir[500];

extern char use_tinker;
extern char use_cole05;

#define STD_STR 500
#define LONG_STR 2000
#define SUCCESS 1
#define FAILURE 0
#define IO_FAILURE 10
#define MEM_FAILURE 11
#define NUM_FAILURE 12

#endif

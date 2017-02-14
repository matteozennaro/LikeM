#ifndef RUNTIMEMSG
#define RUNTIMEMSG

extern void frame(char str[]);
extern void mpi_exit(int mystatus);
extern int check_status(int mystatus);
extern void read_err(char paramname[]);
extern void print_help();

#endif

// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
// KSA
// MPAS_ALL_TASKS_PRINT is undefined as default, but if encountering endrun or
// other run-time error, need to define MPAS_ALL_TASKS_PRINT so that correct stack 
// trace is passed to cesm.log. 
// Without MPAS_ALL_TASKS_PRINT, the stack trace comes from the master proc, which
// is not necessarily where the error occurs.
// Defining MPAS_ALL_TASKS_PRINT allows all processes write into standard error
// unit, which seems to be used to show stack trance, and goes to cesm.log
// KSA


#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#define UNDERSCORE
#define MPAS_NO_LOG_REDIRECT
//#define MPAS_ALL_TASKS_PRINT


#ifdef UNDERSCORE
#define open_streams open_streams_
#define close_streams close_streams_
#else
#ifdef DOUBLEUNDERSCORE
#define open_streams open_streams__
#define close_streams close_streams__
#endif
#endif

int fd_out, fd_err;

void open_streams(int * id)
{
   char fname[128];

#ifndef MPAS_NO_LOG_REDIRECT

#ifndef MPAS_ALL_TASKS_PRINT
   if(*id == 0){
	   sprintf(fname, "log.%4.4i.err", *id);
	   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_err, 2) < 0) {
		   printf("Error duplicating STDERR\n");
		   return;
	   }

	   sprintf(fname, "log.%4.4i.out", *id);
	   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_out, 1) < 0) {
		   printf("Error duplicating STDOUT\n");
		   return;
	   }
   } else {
	   sprintf(fname, "/dev/null");
	   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_err, 2) < 0) {
		   printf("Error duplicating STDERR\n");
		   return;
	   }

	   sprintf(fname, "/dev/null");
	   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_out, 1) < 0) {
		   printf("Error duplicating STDOUT\n");
		   return;
	   }
   }
#else // MPAS_ALL_TASKS_PRINT
   sprintf(fname, "log.%4.4i.err", *id);
   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
   if (dup2(fd_err, 2) < 0) {
      printf("Error duplicating STDERR\n");
      return;
   }

   sprintf(fname, "log.%4.4i.out", *id);
   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
   if (dup2(fd_out, 1) < 0) {
      printf("Error duplicating STDOUT\n");
      return;
   }
#endif // MPAS_ALL_TASKS_PRINT

#else // MPAS_NO_LOG_REDIRECT

#ifndef MPAS_ALL_TASKS_PRINT
   if(*id != 0){
	   sprintf(fname, "/dev/null");
	   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_err, 2) < 0) {
		   fprintf(stderr,"Error duplicating STDERR\n");
		   return;
	   }

//	   sprintf(fname, "/dev/null");
//	   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
//	   if (dup2(fd_out, 1) < 0) {
//		   fprintf(stderr, "Error duplicating STDOUT\n");
//		   return;
//	   }
   }
#endif //MPAS_ALL_TASKS_PRINT

#endif //MPAS_NO_LOG_REDIRECT
}

void close_streams()
{

}

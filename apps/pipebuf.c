#/*
# PipeBuf version 1.2
# To compile run:           sh pipebuf.c
# Help is then available:   ./pipebuf -h

# Copyright (C) 1998 Jan Kratochvil <short@ucw.cz>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; you must use exactly version 2.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may download a copy of the GNU General Public License from URL
#   http://www.opensource.org/gpl-license.html
# If not, write to the Free Software Foundation,
# Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

echo "Compiling PipeBuf..."
a="cc"
b="-s -Wall -O6 -fexpensive-optimizations -fomit-frame-pointer -D_GNU_SOURCE=1"
c="-o `basename "$0" .c` $0"
echo "$a $b $c"
if $a $b $c;then echo -n
else echo "$a $c"
  if $a $c;then echo -n
  else echo "Failed - please check the output.";exit
  fi
fi
echo "done."
exit
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <limits.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/msg.h>
#include <getopt.h>
#include <errno.h>
#include <time.h>

#define BUFSZ (12<<10) /* in KB */
#define BUFWARN (80) /* in percents */
#define WARNTIM (5) /* in seconds */
#define MAX_XFER (PIPE_BUF)
#undef DEBUG

#ifndef SHMMAX
#define SHMMAX 0x2000000
#endif

#ifndef __NORETURN
#if __GNUC__ >= 2
#define __NORETURN __attribute__((__noreturn__))
#else
#define __NORETURN
#endif
#endif

#define RDID (0)
#define WRID (1)

#ifndef DEBUG
#define dbg(cmd)
#else
#define dbg(cmd) cmd
#endif
#define PNAME_LEN (16)
#define bufw(n) ((n)==bufsz?0:(n))
#define failf(name) do { fprintf(stderr,"%s: ",pname); perror(name"()"); exit(EXIT_FAILURE); } while (0)
#ifndef max
#define max(a,b) ((a)>=(b)?(a):(b))
#endif
#ifndef min
#define min(a,b) ((a)<=(b)?(a):(b))
#endif

const char version[]="This is PipeBuf, version 1.2\n";

long bufsz=BUFSZ<<10,bufwarn=BUFWARN;
int prefill,quiet,verbose;

char *pname;
int shmid=-1,dis_cleanup=0,myself;
pid_t other=-1;
volatile struct {
	long rp,wp;
	int eof;
	time_t up;
	struct {
		char sent;
		int msqid;
		} p[2];
	} *comm;
struct msgbuf smsgbuf={1,{0}};
struct msgbuf rmsgbuf;

static void cleanup(void)
{
	dbg(fprintf(stderr,"%s: cleanup()\n",pname));
	if (dis_cleanup) return;
	if (other!=-1) kill(other,SIGTERM);
	shmctl(shmid,IPC_RMID,NULL);
	if (comm) {
		dbg(fprintf(stderr,"%s: msgctl(%d,IPC_RMID), myself=%d\n",pname,comm->p[ myself].msqid, myself));
		msgctl(comm->p[ myself].msqid,IPC_RMID,NULL);
		dbg(fprintf(stderr,"%s: msgctl(%d,IPC_RMID),!myself=%d\n",pname,comm->p[!myself].msqid,!myself));
		msgctl(comm->p[!myself].msqid,IPC_RMID,NULL);
		}
	dis_cleanup=1;
	exit(EXIT_FAILURE);
}

static void wake(void)
{
	dbg(fprintf(stderr,"%s: waking (sent=%d)...\n",pname,comm->p[myself].sent));
	if (comm->p[myself].sent) return;
	dbg(fprintf(stderr,"%s: msgsnd(msqid=%d), myself=%d\n",pname,comm->p[!myself].msqid,myself));
/*	dbg(fprintf(stderr,"%s: sleeping 5 sec...\n",pname));
	sleep(5);
	dbg(fprintf(stderr,"%s: sleeping 5 sec done\n",pname));*/
	if (msgsnd(comm->p[!myself].msqid,(struct msgbuf *)&smsgbuf,0,IPC_NOWAIT)) if (errno!=EAGAIN) failf("msgsnd");
	comm->p[myself].sent++;
	dbg(fprintf(stderr,"%s: waked (sent=%d)\n",pname,comm->p[myself].sent));
}

static void shake(void)
{
	wake();
	dbg(fprintf(stderr,"%s: waiting (sent=%d)...\n",pname,comm->p[myself].sent));
	if (comm->eof) dbg(fprintf(stderr,"%s: EOF => breaking out of shake()!\n",pname));
	else {
		dbg(fprintf(stderr,"%s: msgrcv(msqid=%d),myself=%d\n",pname,comm->p[myself].msqid,myself));
/*		dbg(fprintf(stderr,"%s: sleeping 5 sec...\n",pname));
		sleep(5);
		dbg(fprintf(stderr,"%s: sleeping 5 sec done\n",pname));*/
		if (msgrcv(comm->p[myself].msqid,&rmsgbuf,0,0,0)) failf("msgrcv");
		}
	comm->p[!myself].sent--;
	dbg(if (comm->p[!myself].sent<0) {
		fprintf(stderr,"%s: FATAL - .sent=%d (<0)!\n",pname,comm->p[!myself].sent);
		exit(EXIT_FAILURE);
		});
	dbg(fprintf(stderr,"%s: wait returned (sent=%d)\n",pname,comm->p[myself].sent));
}

static void warnbuf(void)
{
long crp,cwp,bufused,cup;
float percused;
char *tim;
	if (!comm->up||quiet||comm->eof||time(NULL)-comm->up<WARNTIM) return;
	comm->up=cup=time(NULL);
	crp=comm->rp; cwp=comm->wp;
	bufused=crp-cwp+bufsz*!(crp>cwp);
	if ((percused=(float)bufused*100/bufsz)>=bufwarn) return;
	tim=ctime(&cup);
	*strchr(tim,'\n')='\0';
	fprintf(stderr,"%s: %s: WARNING - Low buffer fill-up: %8ld of %8ld (%2.1f%%)\n",
		pname,tim,bufused,bufsz,percused);
}

static __NORETURN void usage(void)
{
	fprintf(stderr,"\
%s\
This command offers the pipe buffering:\n\
\n\
Usage: pipebuf [-b|--buffer <size in KB>] [-p|--prefill] [-w|--warning <percent>]\n\
               [-q|--quiet] [-v|--verbose] [-h|--help] [-V|--version]\n\
\n\
  -b, --buffer <size in KB>\tSpecify buffer size (1-%dKB, def=%dKB)\n\
  -p, --prefill\t\t\tFill the buffer before first write\n\
  -w, --warning <percent>\tNo-buffer-data warnings threshold (0-100%%, def=%d%%)\n\
  -q, --quiet\t\t\tDon't print warnings\n\
  -v, --verbose\t\t\tInform about phases of transfer\n\
  -h, --help\t\t\tPrint a summary of the options\n\
  -V, --version\t\t\tPrint the version number\n\
",version,(SHMMAX>>10)-1,BUFSZ,BUFWARN);
	exit(EXIT_FAILURE);
}

const struct option longopts[]={
{"buffer" ,1,0,'b'},
{"prefill",0,0,'p'},
{"warning",1,0,'w'},
{"quiet"  ,0,0,'q'},
{"verbose",0,0,'v'},
{"help"   ,0,0,'h'},
{"version",0,0,'V'}};

int main(int argc,char **argv)
{
long cfp;
int r,optc;
caddr_t buf;
char *s;

	pname=*argv;
	atexit(cleanup);
	signal(SIGTERM,(void (*)(int))cleanup);
	signal(SIGQUIT,(void (*)(int))cleanup);
	signal(SIGINT ,(void (*)(int))cleanup);
	signal(SIGHUP ,(void (*)(int))cleanup);
	while ((optc=getopt_long(argc,argv,"b:pw:qvhV",longopts,NULL))!=EOF) switch (optc) {
		case 'b':
			errno=EINVAL;
			bufsz=strtol(optarg,&s,0);
			if (*s!='\0'||bufsz<1||bufsz<<10>=SHMMAX) { perror(optarg); usage(); }
			bufsz<<=10;
			break;
		case 'p':
			prefill=1;
			break;
		case 'w':
			errno=EINVAL;
			bufwarn=strtol(optarg,&s,0);
			if (*s!='\0'||bufwarn<0||bufwarn>100) { perror(optarg); usage(); }
			break;
		case 'q':
			quiet=1;
			verbose=0;
			break;
		case 'v':
			verbose=1;
			quiet=0;
			break;
		case 'V':
			fprintf(stderr,version);
			exit(EXIT_FAILURE);
		default: /* also 'h' */
			usage();
			break;
		}
	if ((shmid=shmget(IPC_PRIVATE,bufsz+sizeof(*comm),0600|IPC_CREAT|IPC_EXCL))==-1) failf("shmget");
	if ((int)(buf=shmat(shmid,0,0))==-1) failf("shmat");
	comm=(void *)buf+bufsz;
	bzero((void *)comm,sizeof(*comm));
	if ((comm->p[RDID].msqid=msgget(IPC_PRIVATE,0777|IPC_CREAT|IPC_EXCL))==-1) failf("msgget");
	if ((comm->p[WRID].msqid=msgget(IPC_PRIVATE,0777|IPC_CREAT|IPC_EXCL))==-1) failf("msgget");
	if (!prefill) comm->up=time(NULL);
	other=fork();
	if (other) {
	/* Read process */
		dbg(fprintf(stderr,"%s: started rd\n",pname));
		myself=RDID;
		strncat(pname,"-rd",max(PNAME_LEN-strlen(pname)-1,0));
		if (close(STDOUT_FILENO)) failf("close");
		if (verbose) fprintf(stderr,"%s: Using buffer %ldKB%s...\n",pname,bufsz>>10,(prefill?", filling":""));
		dbg(fprintf(stderr,"%s: pname check rd\n",pname));
		for (;;) {
			if (bufw(comm->rp+1)==comm->wp) {
				if (!comm->up) {
					comm->up=time(NULL);
					if (verbose) fprintf(stderr,"%s: Buffer filled-up, starting transfer...\n",pname);
					}
				shake();
				continue;
				}
			warnbuf();
			cfp=comm->wp;
			dbg(fprintf(stderr,"%s: rp=%ld, wp=%ld",pname,comm->rp,cfp));
			cfp=(cfp<=comm->rp?bufsz-comm->rp-!cfp:cfp-1-comm->rp);
			dbg(fprintf(stderr,", read(%d,%08lx,%ld)\n",STDIN_FILENO,(long)&buf[comm->rp],min(cfp,MAX_XFER)));
			if ((r=read(STDIN_FILENO,&buf[comm->rp],min(cfp,MAX_XFER)))==-1) failf("read");
			dbg(fprintf(stderr,"%s: rp=%ld, wp=%ld, read=%d\n",pname,bufw(comm->rp+r),comm->wp,r));
			if (r) comm->rp=bufw(comm->rp+r);
			else {
				comm->eof=1;
				if (!comm->up) {
					comm->up=time(NULL);
					if (verbose) fprintf(stderr,"%s: Reached EOF before buffer fill-up, starting transfer...\n",pname);
					}
				wake();
				break;
				}
			if (comm->up) wake();
			}
		if (verbose) fprintf(stderr,"%s: All input data read, waiting for write completion...\n",pname);
		if (waitpid(other,NULL,0)!=other) failf("waitpid");
		}
	else {
	/* Write process */
		dbg(fprintf(stderr,"%s: started wr\n",pname));
		myself=WRID;
		strncat(pname,"-wr",max(PNAME_LEN-strlen(pname)-1,0));
		other=getppid();
		if (close(STDIN_FILENO)) failf("close");
		dbg(fprintf(stderr,"%s: pname check wr\n",pname));
		for (;;) {
			if (comm->eof&&comm->rp==comm->wp) break;
			while (!comm->eof&&(!comm->up||comm->rp==comm->wp)) shake();
			if (comm->eof&&comm->rp==comm->wp) break;
			cfp=comm->rp;
			dbg(fprintf(stderr,"%s: rp=%ld, wp=%ld",pname,cfp,comm->wp));
			cfp=(cfp<=comm->wp?bufsz-comm->wp:cfp-comm->wp);
			dbg(fprintf(stderr,", write(%d,%08lx,%ld)\n",STDOUT_FILENO,(long)&buf[comm->wp],min(cfp,MAX_XFER)));
			if ((r=write(STDOUT_FILENO,&buf[comm->wp],min(cfp,MAX_XFER)))==-1) failf("write");
			dbg(fprintf(stderr,"%s: rp=%ld, wp=%ld, write=%d\n",pname,comm->rp,bufw(comm->wp+r),r));
			if (!(comm->wp=bufw(comm->wp+r))&&comm->rp) continue;
			wake();
			warnbuf();
			}
		}
	if (verbose) fprintf(stderr,"%s: Ending operation (sent=%d).\n",pname,comm->p[myself].sent);
	return(EXIT_SUCCESS);
}

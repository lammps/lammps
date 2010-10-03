#include "stdio.h"
#include "string.h"
#include "files.h"

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

void replace(char *file, char *header, int n, char **lines)
{
  FILE *fpr = fopen(file,"r");
  FILE *fpw = fopen("tmp.file","w");

  char line[MAXLINE];
  while (fgets(line,MAXLINE,fpr)) {
    if (strstr(line,header)) {
      fprintf(fpw,"%s",line);
      for (int i = 0; i < n; i++) {
	fgets(line,MAXLINE,fpr);
	fprintf(fpw,"%s",lines[i]);
      }
    } else fprintf(fpw,"%s",line);
  }

  fclose(fpr);
  fclose(fpw);
  rename("tmp.file",file);
}

/* ---------------------------------------------------------------------- */

char **extract(char *file, char *header, int n, char **lines)
{
  FILE *fp = fopen(file,"r");

  char line[MAXLINE];
  while (fgets(line,MAXLINE,fp)) {
    if (strstr(line,header)) {
      for (int i = 0; i < n; i++) {
	fgets(line,MAXLINE,fp);
	sprintf(lines[i],"%s",line);
      }
      break;
    }
  }

  fclose(fp);
}

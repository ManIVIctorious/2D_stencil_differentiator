
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2);

int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2){

    int i, rows, comment_flag, control;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE *fd;

    if(inputfile == NULL){
        fprintf(stderr, "No input-file given, reading from stdin\n");
        fd = stdin;
    }else{
        fd = fopen(inputfile, "r");
        if(fd == NULL){
            fprintf(stderr, "\n(-) ERROR opening input-file: \"%s\"", inputfile);
            fprintf(stderr, "\n    Exiting...\n\n");
            exit(1);
        }
    }

    rows = 0;
    while(fgets(buffer, sizeof(buffer), fd) != NULL){

        // check if the first character in buffer is a comment char,
        //  if yes jump to next line
        comment_flag = 0;
        for(i=0; i<(int)strlen(comment); ++i){
            if(buffer[0] == comment[i]){
                comment_flag = 1;
                break;
            }
        }
        if(comment_flag == 1) continue;

        // copy "buffer" with stripped comments to new buffer "line"
        line = strtok(buffer, comment);
        if(line == NULL) continue;

        // remove leading white spaces and tabulators
        //  and skip empty lines
        while(isspace(*line)) line++;
        if(strlen(line) == 0) continue;

        // At this point the requested input line is stripped of
        //  comments and blank lines. From here on the parsing starts:
        //printf("%s\n", line);
//-----------------------------------------------------------------------------------

        if((strncmp(line, "N", 1) == 0) && (isspace((int) line[1]) != 0)){
            control = sscanf(line, "%*s  %d  %d", &(*nq1), &(*nq2));
            if(control != 2){
                fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
                fprintf(stderr, "\n    Line containing point number information\n");
                fprintf(stderr, "\n    does not provide appropriate information\n");
                fprintf(stderr, "\n    Aborting - please check your input...\n\n");
                exit(control);
            }
            continue;
        }

        (*q1)  = realloc((*q1),  (rows + 1) * sizeof(double));
        (*q2)  = realloc((*q2),  (rows + 1) * sizeof(double));
        (*V)   = realloc((*V),   (rows + 1) * sizeof(double));

        control = sscanf(line, "%lf  %lf  %lf", &(*q1)[rows], &(*q2)[rows], &(*V)[rows]);
        if(control != 3){
          fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
          fprintf(stderr, "\n    Aborting - please check your input...\n\n");
          exit(control);
        }

        ++rows;
    }
    fclose(fd); fd = NULL;

    return rows;
}

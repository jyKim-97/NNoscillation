#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
    FILE *fp;
    char *line = NULL;
    char *tokens;
    char *arg;
    size_t len = 0;
    ssize_t read;
    // fp = fopen("./test.txt", "r");
    fp = fopen("./test_cell.csv", "r");
    // if (fp == NULL)  exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
            printf("%s", line);
            tokens = strtok(line, ",");
            while (tokens != NULL){
              // arg = (char) tokens;
              // arg = atof(tokens);
            printf("%s\n", tokens);
            tokens = strtok(NULL, ",");
            // i = i + 1;

            }

    }
    free(line);
    fclose(fp);

}
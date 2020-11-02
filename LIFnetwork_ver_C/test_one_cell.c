#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ArtNetSimulTools.h"
#include <unistd.h>

int main(int argc, char **argv){
    int opt;

    while ((opt=getopt(argc, argv, "a:b:c:")) != -1) {
        switch (opt) {
            case 'a':
                printf("1\n");
                break;
            case 'b':
                printf("2\n");
                break;
            case '?':
                printf("3\n");
                break;
        }
    }
    
    // set parameter

}
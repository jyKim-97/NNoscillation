#!bin/bash

# create object file
gcc -c ArtNetSimulTools.c

# create static libraries
ar rc libArtNet ArtNetSimulTools.o mt64.o

# compile
gcc -o3 -o runNetwork runNetwork.c -L./ -lm -lArtNet

#!/bin/bash

gcc -O3 mcc.c -ligraph -lstdc++ -o mcc -Wall
gcc -O3 mcstats.c -ligraph -lstdc++ -o mcstats -Wall
gcc -O3 mcextract.c -ligraph -lstdc++ -o mcextract -Wall

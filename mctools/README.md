Here you will find source code for each of the command line applications that makes up mctools. These are all written in C and make extensive use of the the igraph library (http://igraph.sf.net). To compile, igraph must be in the appropriate include and library paths and be version 0.6.5 or later. The following commands can then be used for compilation:

	gcc -O3 mcc.c -ligraph -lstdc++ -o mcc -Wall
	gcc -O3 mcstats.c -ligraph -lstdc++ -o mcstats -Wall
	gcc -O3 mcextract.c -ligraph -lstdc++ -o mcextract -Wall

There are a number of compile time flags that can be used to enable non-standard features:
-DDEBUG        : output debugging information.
-DBRENCHMARK   : output timing information for major steps.
-DEXPERIMENTAL : include experimental features e.g., OpenMP support.

For ready-to-use pre-compiled versions of this code see the bin folder in the project root.

If you make use of this software in your work we request that you cite:

T.E. Gorochowski, C.S. Grierson and M. di Bernardo. "MCTools: motif clustering toolkit for complex networks." 2013. (in preparation) (http://chofski.github.com/mctools/)

This software has been developed by Thomas Gorochowski (@chofski). All code is distributed under the OSI recognised [Non-Profit Open Software License version 3.0 (NPOSL-3.0)](http://www.opensource.org/licenses/NOSL3.0).

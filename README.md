### Motif Clustering Tools for Complex Networks

Network motifs have been proposed as possible building blocks for complex networks. 

The `mctools` package provides several command-line applications to calculate the overall clustering of a given motif within a network (`mcc`), break this clustering down into a distribution of the clustering types present (`mcstats`), and extract these clustered subgraphs from the original network (`mcextract`).

Ready-to-run binary versions of the programs are provided for Windows, Linux and Mac OS X. However, for best performance we advise compilation from source. See documentation for further information.

If you make use of this software in your work we request that you cite:

>T.E. Gorochowski, C.S. Grierson and M. do Bernardo. mctools: Motif Clustering Analysis for Complex Networks. 2012. (in preparation) [[doi](http://chofski.github.com/mctools/)]

This software has been developed by Thomas Gorochowski (@chofski) and makes extensive use of the [`igraph`](http://igraph.sf.net) library. All code is distributed under the OSI recognised [Non-Profit Open Software License version 3.0 (NPOSL-3.0)](http://www.opensource.org/licenses/NOSL3.0).
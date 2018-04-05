# Motif Clustering Tools for Complex Networks

The `mctools` package provides several command-line applications to calculate motif clustering statistics for a chosen motif within a complex network (`mcc`), break this clustering down into a distribution of the clustering types present (`mcstats`), and extract these clustered subgraphs from the original network (`mcextract`).

Ready-to-run binary versions of the programs are provided for Windows, Linux and Mac OS X (see `bin` folder). However, for best performance we advise compilation from source (see `mctools` folder).

If you make use of this software in your work we request that you cite:

>T.E. Gorochowski, C.S. Grierson and M. di Bernardo. "Organization of feed-forward loop motifs reveals architectural principles in natural and engineered networks." Science Advances 4: eaap97512015 (2018)

This software has been developed by Thomas Gorochowski (@chofski) and makes extensive use of the [`igraph`](http://igraph.sf.net) library. All code is distributed under the OSI recognised [Non-Profit Open Software License version 3.0 (NPOSL-3.0)](http://www.opensource.org/licenses/NOSL3.0).

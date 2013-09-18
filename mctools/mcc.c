/*================================================================================================
 *
 *  mcc.c - Motif Clustering Coefficient (MCC)
 *
 *  Calculates the motif clustering coefficient for a graph and a z-score based on a random null
 *  model that maintains the same size of graph (in nodes) and number of motifs. The directedness 
 *  of the input graph is important as this is then used when finding the correct motif from an 
 *  isomorphic class ID. A debug mode can be enabled at compile time, using the flag -DDEBUG which 
 *  should help understand the steps being performed. A benchmarking flag -DBRENCHMARK is also 
 *  included to time various parts of the process.
 *
 *  This command outputs two files:
 *     1. PREFIX_samples.txt - motif clustering coefficient values for the random samples.
 *     2. PREFIX_stats.txt   - statistics from the run.
 *
 *  Warning: The implemented method here is within the confines of the total number of motifs
 *           in a graph being in the range of hundreds of thousands. Use the debug mode to check 
 *           the number of motifs that have been found and to estimate the running time if problems 
 *           arise.
 *
 *------------------------------------------------------------------------------------------------
 *
 *  To compile, use the following command:
 *
 *     gcc -I INC_DIR -L LIB_DIR -O3 mcc.c -ligraph -lstdc++ -o mcc -Wall
 *
 *  where INC_DIR is the include directory and LIB_DIR is the library directory. The igraph
 *  library is required to compile this program and can be found at http://igraph.sourceforge.net/
 *  On some systems it is also necessary to link against the GSL Big Number library using -lgmp.
 *
 *------------------------------------------------------------------------------------------------
 *
 *  Usage: mcc FILENAME PREFIX SAMPLE TRIALS MOTIF_SIZE MOTIF_ID
 *
 *         FILENAME    : Graph filename (GML format).
 *         PREFIX      : Prefix to use on output files.
 *         SAMPLE      : Size of the sample to generate z-score with.
 *         TRIALS      : Number of trails to place motifs in random graph (normally 200).
 *         MOTIF_SIZE  : Size of the motif to consider (3 or 4 nodes).
 *         MOTIF_ID    : Isomorphic class ID (from igraph) for the motif.
 *
 *------------------------------------------------------------------------------------------------
 *
 *  Copyright (C) 2013 Thomas E. Gorochowski <tom@chofski.co.uk>
 *
 *  This software released under the Open Source Initiative (OSI) approved Non-Profit Open 
 *  Software License ("Non-Profit OSL") 3.0. This software is distributed in the hope that 
 *  it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *================================================================================================*/


/* Required libraries */
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <igraph.h>


/* Maximum attempts placing a single motif before giving up */
long int MAX_MOTIF_TRIALS;


/* Calculates the motif clustering coefficient. */
int motif_clustering (double *res, igraph_t* graph, igraph_t *motif);

/* Calculates a z-score for a motif clustering coefficient and a set of random samples. */
int z_score (double *res, double mcc, igraph_vector_t *samples);

/* Generates random graphs of a given number of nodes, containing a specified number of different
 * motif types. Uses the function calc_sample to calculate a graph. */
int calc_samples (igraph_vector_t *res, igraph_t *graph, igraph_t *motif, 
						igraph_integer_t count, igraph_integer_t nodes, int samples);

/* Calculates a single random sample, containing a specified number of different motif types. */
int calc_sample (igraph_t *res, igraph_t* graph, igraph_t *motif, 
					  igraph_integer_t count, igraph_integer_t nodes);

/* Count the number of motifs in a graph. */
igraph_integer_t motif_count (igraph_t *graph, igraph_t *motif);

/* Print usage information. */
void print_usage (void);


/*------------------------------------------------------------------------------------------------*/

/* Main function. */
int main (int argc, const char * argv[])
{
	char filename[1000];
	FILE *gFile, *outFile;
	igraph_t G, motif;
	double resMCC, resZScore;
	int suc;
	igraph_integer_t x, count;
	igraph_vector_t samples;
#ifdef BENCHMARK
	clock_t ctime_1, ctime_2;
	ctime_1 = clock();
#endif
	
	/* Check that there are enough arguments */	
	if (argc == 2 && strcmp(argv[1], "-h") == 0) {
		print_usage();
		return 0;
	}
	if ((argc < 2) || (argc < 7) || (argc > 7)) {
		printf("Invalid number of arguments.\n");
		return 1;
	}
	
	/* Seed the random number generator */
	srand(time(NULL));
	
	MAX_MOTIF_TRIALS = (long int)atoi(argv[4]);
	
	/* Load the user specified topology (GML Format) */
	gFile = fopen(argv[1], "r");
	igraph_read_graph_gml(&G, gFile);
	fclose(gFile);
	
	/* Create list of motif graphs */
	igraph_isoclass_create(&motif, atoi(argv[5]), atoi(argv[6]), igraph_is_directed(&G));
	
	suc = motif_clustering(&resMCC, &G, &motif);
	
	count = motif_count (&G, &motif);
	suc = calc_samples(&samples, &G, &motif, count, igraph_vcount(&G), atoi(argv[3]));
	suc = z_score(&resZScore, resMCC, &samples);
	
	printf("Motif clustering coefficient = %.8f, z-score = %.8f\n", resMCC, resZScore);
	fflush(stdout);
	
	/* Output random samples used to calculate z-score */
	sprintf(filename, "%s_samples.txt", argv[2]);
	outFile = fopen(filename, "w");
	for (x=0; x<igraph_vector_size(&samples); x++) {
		fprintf(outFile, "%.8f\n", (double)VECTOR(samples)[(long int)x]);
	}
	fclose(outFile);
	
	/* Output the statistics from the run */
	sprintf(filename, "%s_stats.txt", argv[2]);
	outFile = fopen(filename, "w");
	fprintf(outFile, "Nodes, Edges, MCC, Z-Score\n");
	fprintf(outFile, "%li, %li, %.8f, %.8f\n", (long int)igraph_vcount(&G), 
			  (long int)igraph_ecount(&G), 
			  resMCC, resZScore);
	fclose(outFile);
	
	/* Free used memory and return */
	igraph_destroy(&motif);
	igraph_destroy(&G);
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Motif clustering and z-score calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
#endif
	
	return 0;
}

/*------------------------------------------------------------------------------------------------*/

int motif_clustering (double *res, igraph_t* graph, igraph_t *motif)
{
	igraph_integer_t rotSym, motifEdges;
	long int i, j, motifSize, mapsCount, uniqueMotifs,
	totSharedVerts, actSharedVerts, posSharedVerts, k, p, found, actMaps;
	igraph_t subGraphs;
	igraph_vs_t vs;
	igraph_vector_t *curi, *curj;
	igraph_vector_ptr_t maps;
	igraph_vector_ptr_init(&maps, 0);
	
#ifdef BENCHMARK
	clock_t ctime_1, ctime_2;
	ctime_1 = clock();
#endif
	
	/* 1. Calculate the size and motif symmetries */
	motifSize = (long int)igraph_vcount(motif);
	igraph_count_subisomorphisms_vf2(motif, motif, NULL, NULL, NULL, NULL, &rotSym, NULL, NULL, NULL);
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Motif count calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	/* 2. Find submorphisms between graph and motif */
	igraph_get_subisomorphisms_vf2(graph, motif, NULL, NULL, NULL, NULL, &maps, NULL, NULL, NULL);
	
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("All mappings calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	mapsCount = (long int)igraph_vector_ptr_size(&maps);
	
	/* Clean up maps list (only required for directed graphs) */
	if (igraph_is_directed(graph) != 0) {
		actMaps = 0;
		motifEdges = igraph_ecount(motif);
		for (i=0; i<mapsCount; i++) {
			
			curi = (igraph_vector_t *)VECTOR(maps)[i];
			
			/* Generate a vertex selector from our IDs vector */
			igraph_vs_vector(&vs, curi);
			
			/* Extract the subgraph using this vertex selector */
			igraph_induced_subgraph(graph, &subGraphs, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
			
			if (igraph_ecount(&subGraphs) != motifEdges) {
				/* Remove mapping set first value to -1*/
				VECTOR(*curi)[0] = -1;
			}
			else {
				actMaps++;
			}

			igraph_vs_destroy(&vs);
			igraph_destroy(&subGraphs);
		}
	}
	else {
		actMaps = mapsCount;
	}

#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Mappings cleaned up, calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	/* 3. Calculate unique motifs */
	uniqueMotifs = actMaps / rotSym;
	
	/* 4. Find actual and total possible shared vertices */
	totSharedVerts = 0;
	
	/* OpenMP parallel calculations for large graphs */
#ifdef EXPERIMENTAL
#pragma omp parallel for default(none) private(i,j,k,p,found,curi,curj) shared(motifSize,mapsCount,totSharedVerts,maps)
#endif
	for (i=0; i<mapsCount; i++) {
		curi = (igraph_vector_t *)VECTOR(maps)[i];
		if (VECTOR(*curi)[0] != -1) {
			
			for (j=i+1; j<mapsCount; j++) {
				curj = (igraph_vector_t *)VECTOR(maps)[j];
				
				if (VECTOR(*curj)[0] != -1) {
					found = 0;
					/* Loop through both vectors and see if matching elements are found */
					for (k=0; k<motifSize; k++) {
						for (p=0; p<motifSize; p++) {
							if (VECTOR(*curi)[k] == VECTOR(*curj)[p]) {
								found++;
								break;
							}
						}
					}
					if (found < motifSize) {
						totSharedVerts += found;
					}
				}
			}
		}
	}
	
	/* OpenMP parallel calculations for large graphs */
#ifdef EXPERIMENTAL
#pragma omp parallel for default(none) private(i) shared(posSharedVerts)
#endif
	posSharedVerts = 0;
	for (i=0; i<uniqueMotifs; i++) {
		posSharedVerts += (motifSize-1)*(uniqueMotifs-i-1);
	}
	
	actSharedVerts = totSharedVerts/(rotSym*rotSym);
	
#ifdef DEBUG
	printf(" mapsCount:%i\n actShVerts:%i\n totShVerts:%i\n rotSym:%i\n motifSize:%i\n uniqueMotifs:%i\n posShVerts:%i\n", 
			 (int)mapsCount, (int)actSharedVerts, (int)totSharedVerts,
			 (int)rotSym, (int)motifSize, (int)uniqueMotifs, (int)posSharedVerts);
	fflush(stdout);
#endif
	
	/* 5. Calculate motif clustering coefficient */
	*res = (double)(actSharedVerts)/(double)posSharedVerts;
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Motif clustering calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
#endif
	
	return 0;
}

/*------------------------------------------------------------------------------------------------*/

int z_score (double *res, double mcc, igraph_vector_t *samples)
{
	igraph_integer_t j, sampleSize, actualSamples;
	double sumMCC, sumMCC2, avgMCC, avgMCC2;
	
	/* Find the sample size */
	sampleSize = igraph_vector_size(samples);
	actualSamples = 0;
	
	/* Calculate z-score from the samples */
	sumMCC  = 0;
	sumMCC2 = 0;
	for (j=0; j<sampleSize; j++) {
		if ((double)VECTOR(*samples)[(long int)j] >= 0.0) {
			sumMCC  += (double)VECTOR(*samples)[(long int)j];
			sumMCC2 += pow((double)VECTOR(*samples)[(long int)j], 2.0);
			actualSamples++;
		}
	}
	avgMCC  = sumMCC / (double)actualSamples;
	avgMCC2 = sumMCC2 / (double)actualSamples;
	*res = (mcc - avgMCC) / sqrt(avgMCC2 - (avgMCC*avgMCC));
	
	return 0;
}

/*------------------------------------------------------------------------------------------------*/

int calc_samples (igraph_vector_t *res, igraph_t *graph, igraph_t *motif, 
						igraph_integer_t count, igraph_integer_t nodes, int samples)
{
	int s, suc, flag;
	double *badRes;
	igraph_t Gs;
#ifdef BENCHMARK
	clock_t ctime_1, ctime_2;
#endif
	
	/* Initialise the results vector to the correct size */
	igraph_vector_init(res, samples);
	
	/* Attempt to generate the number of samples required */
	flag = 0;
	
	/* OpenMP parallelisation */
#ifdef EXPERIMENTAL
#pragma omp parallel for default(none) private(s,j,flag,shared) shared(suc,ctime_1,ctime_2,maps,Gs,motifs,samples,nodes,counts,graph,res)
#endif
	for (s=0; s<samples; s++) {
		
#ifdef BENCHMARK
		ctime_1 = clock();
#endif
		
#ifdef DEBUG
		printf("Generating sample %li of %li\n", (long int)s+1, (long int)samples);
		fflush(stdout);
#endif
		
		/* Generate a sample */
		suc = calc_sample(&Gs, graph, motif, count, nodes);
		
#ifdef BENCHMARK
		ctime_2 = clock();
		printf("Sample calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
				 (double)CLOCKS_PER_SEC);
#endif
		
		if (suc == 1) {
			badRes = (double *)igraph_vector_e_ptr(res, (igraph_integer_t)s);
			*badRes = -1.0;
			flag = 1;
		}
		else {
			/* Calculate the stats on the graph */
			motif_clustering(igraph_vector_e_ptr(res, (igraph_integer_t)s), &Gs, motif);
			/* Free used memory for the sample */
			igraph_destroy(&Gs);
		}
	}
	
	if (flag == 1) {
		/* Some of the samples might not have been generated correctly */
		return 1;
	}	
	
	return 0;
}

/*------------------------------------------------------------------------------------------------*/

int calc_sample (igraph_t *res, igraph_t* graph, igraph_t *motif, 
					  igraph_integer_t count, igraph_integer_t nodes)
{
	igraph_integer_t j, k, x, curCount, curAdd, newAdd, motifPlaceTrial, edgePlaceTrial,
	oldCount;
	igraph_t *G, *altG;
	igraph_eit_t eit;
	igraph_vector_t mNodes, newEdges;
#ifdef BENCHMARK
	clock_t ctime_1, ctime_2;
#endif
	
	/* Generate the initially empty graph to contain the final sample */
	G = (igraph_t *)malloc(sizeof(igraph_t));
	igraph_empty(G, nodes, igraph_is_directed(graph));
	
	igraph_vector_init(&mNodes, igraph_vcount(motif));
	
	/* Update the iterator for the edges in the motif */
	igraph_eit_create(motif, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
	
	/* Keep adding motifs until the count for current motif is correct */
	oldCount = 0;
	curCount = 0;
	motifPlaceTrial = 0;
	edgePlaceTrial = 0;
	curAdd = (igraph_integer_t)((long int)count / 5);
	if ((long int)curAdd < 1) curAdd = 1;
	
	while ((long int)motifPlaceTrial < MAX_MOTIF_TRIALS) {
		
		/* Create the alternative graph adding motif using this mapping */
		altG = (igraph_t *)malloc(sizeof(igraph_t));
		igraph_copy(altG, G);
		
#ifdef BENCHMARK
		ctime_1 = clock();
#endif
			
#ifdef DEBUG
			printf("Attempting to add %li motifs\n", (long int)curAdd);
			fflush(stdout);
#endif
			
		/* Attempt to add each motif we require */
		igraph_vector_init(&newEdges, (long int)curAdd*(long int)igraph_ecount(motif)*2);
		x = 0;
		for (j=0; j<curAdd; j++) {
			/* Generate random node IDs to use as mapping for motif */
			for (k=0; k<igraph_vcount(motif); k++) {
				VECTOR(mNodes)[(long int)k] = (igraph_real_t)(rand() % (int)nodes);
			}
			IGRAPH_EIT_RESET(eit);
			do {
				VECTOR(newEdges)[(long int)x] = VECTOR(mNodes)[(long int)IGRAPH_FROM(motif, 
																											IGRAPH_EIT_GET(eit))];
				VECTOR(newEdges)[(long int)x+1] = VECTOR(mNodes)[(long int)IGRAPH_TO(motif, 
																											IGRAPH_EIT_GET(eit))];
				x = x + 2;
				IGRAPH_EIT_NEXT(eit);
			} while (!IGRAPH_EIT_END(eit));
		}
		igraph_add_edges(altG, &newEdges, 0);
		igraph_vector_destroy(&newEdges);
		
#ifdef BENCHMARK
		ctime_2 = clock();
		printf("Added %li motifs/edge in %f seconds\n", (long int)curAdd, (double)(ctime_2 - 
																											ctime_1) / (double)CLOCKS_PER_SEC);
		ctime_1 = clock();
#endif
		
		/* Count the motifs in the new graph */
		curCount = motif_count(altG, motif);
		
#ifdef BENCHMARK
		ctime_2 = clock();
		printf("Calculated new motif counts in %f seconds\n", (double)(ctime_2 - ctime_1) / 
				 (double)CLOCKS_PER_SEC);
#endif
		
		if (curCount < count && curCount != oldCount) {
			/* Accept change, recalculate number of motifs to add and loop */
			newAdd = (igraph_integer_t)(((long int)count - (long int)curCount) / (long int)3);
			if ((long int)newAdd < 1) newAdd = 1;
			if ((long int)newAdd < (long int)curAdd) curAdd = newAdd;
			
			if ((long int)motifPlaceTrial < MAX_MOTIF_TRIALS) {
				motifPlaceTrial = 0;
			}
			else {
				edgePlaceTrial++;
			}
			igraph_destroy(G);
			free(G);
			G = altG;
			
			oldCount = curCount;
			
#ifdef DEBUG
			printf("Accepting change, %li motifs of %li, trial %li\n", 
					 (long int)curCount, 
					 (long int)count,
					 (long int)motifPlaceTrial);
			fflush(stdout);
#endif
			
		}
		else if (curCount == count) {			
			/* Counts match for the motif being added and all others <= required number */
			motifPlaceTrial = 0;
			edgePlaceTrial = 0;
			igraph_destroy(G);
			free(G);
			G = altG;
			
#ifdef DEBUG
			printf("Accepting change, %li motifs of %li, trial %li\n", 
					 (long int)curCount, 
					 (long int)count,
					 (long int)motifPlaceTrial);
			fflush(stdout);
#endif
			
			break;
		}
		else {
			/* Exceeded number of allowed motifs, reject trail, update motifs to add and retry */
			curAdd = (igraph_integer_t)((long int)curAdd / (long int)3);
			if ((long int)curAdd <= 1) {
				curAdd = 1;
				if ((long int)motifPlaceTrial < MAX_MOTIF_TRIALS) motifPlaceTrial++;
				else edgePlaceTrial++;
			}
			igraph_destroy(altG);
			free(altG);
			
#ifdef DEBUG
			printf("Rejecting change, %li motifs instead of %li, trial %li\n", 
					 (long int)curCount, 
					 (long int)count,
					 (long int)motifPlaceTrial);
			fflush(stdout);
#endif
			
		}
	}
				
		/* Free memory used */
		igraph_eit_destroy(&eit);
		
		/* Could not place the motifs so return with error */
		if (curCount > count) {
#ifdef DEBUG
			printf("Exceeded number of motif and single edge trials\n");
			fflush(stdout);
#endif
			igraph_vector_destroy(&mNodes);
			igraph_destroy(G);
			free(G);
			return 1;
		}
	
	
	/* Return the valid graph sample */
	igraph_copy(res, G);
	
	/* Free used memory */
	igraph_vector_destroy(&mNodes);
	igraph_destroy(G);
	free(G);
	
	return 0;
}

/*------------------------------------------------------------------------------------------------*/

igraph_integer_t motif_count (igraph_t *graph, igraph_t *motif)
{
	igraph_integer_t rotSym, motifEdges;
	long int i, motifSize, mapsCount, actMaps;
	igraph_t subGraphs;
	igraph_vs_t vs;
	igraph_vector_t *curi;
	igraph_vector_ptr_t maps;
	igraph_vector_ptr_init(&maps, 0);
	
#ifdef BENCHMARK
	clock_t ctime_1, ctime_2;
	ctime_1 = clock();
#endif
	
	/* 1. Calculate the size and motif symmetries */
	motifSize = (long int)igraph_vcount(motif);
	igraph_count_subisomorphisms_vf2(motif, motif, NULL, NULL, NULL, NULL, &rotSym, NULL, NULL, NULL);
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Motif count calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	/* 2. Find submorphisms between graph and motif */
	igraph_get_subisomorphisms_vf2(graph, motif, NULL, NULL, NULL, NULL, &maps, NULL, NULL, NULL);
	
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("All mappings calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	mapsCount = (long int)igraph_vector_ptr_size(&maps);	
	
	/* Clean up maps list (only required for directed graphs) */
	if (igraph_is_directed(graph) != 0) {
		actMaps = 0;
		motifEdges = igraph_ecount(motif);
		for (i=0; i<mapsCount; i++) {
			
			curi = (igraph_vector_t *)VECTOR(maps)[i];
			
			/* Generate a vertex selector from our IDs vector */
			igraph_vs_vector(&vs, curi);
			
			/* Extract the subgraph using this vertex selector */
			igraph_induced_subgraph(graph, &subGraphs, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
			
			if (igraph_ecount(&subGraphs) != motifEdges) {
				/* Remove mapping set first value to -1*/
				VECTOR(*curi)[0] = -1;
			}
			else {
				actMaps++;
			}
			igraph_vs_destroy(&vs);
			igraph_destroy(&subGraphs);
		}
	}
	else {
		actMaps = mapsCount;
	}
	
#ifdef BENCHMARK
	ctime_2 = clock();
	printf("Mappings cleaned up calculated in %f seconds\n", (double)(ctime_2 - ctime_1) / 
			 (double)CLOCKS_PER_SEC);
	fflush(stdout);
	ctime_1 = clock();
#endif
	
	/* 3. Calculate unique motifs */	
	return (igraph_integer_t)(actMaps / rotSym);	
}

/*------------------------------------------------------------------------------------------------*/

void print_usage (void)
{
	printf("mcc FILENAME PREFIX SAMPLE TRIALS MOTIF_SIZE MOTIF_ID\n");
	printf("    FILENAME   - Graph filename (GML format).\n");
	printf("    PREFIX     - Prefix to use on output files.\n");
	printf("    SAMPLE     - Size of the sample to generate z-score with.\n");
	printf("    TRIALS     - Number of trails when placing motifs in random sample.\n");
	printf("    MOTIF_SIZE - Size of the 1st motif to consider (3 or 4 nodes).\n");
	printf("    MOTIF_ID   - Isomorphic class ID (from igraph) for the 1st motif.\n");
}

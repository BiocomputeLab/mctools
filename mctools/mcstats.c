/*===============================================================================================
 *  mcstats.c
 *
 *  Calculates motif clustering stats for a given motif and all the types of clustering it can
 *  take part in.
 *
 *------------------------------------------------------------------------------------------------
 *
 *  To compile use the following command:
 *
 *     gcc -I INC_DIR -L LIB_DIR -O3 mcstats.c -ligraph -lstdc++ -o mcstats
 *
 *  where INC_DIR is the include directory and LIB_DIR is the library directory. The igraph
 *  library is required to compile this program and can be found at http://igraph.sourceforge.net/
 *
 *------------------------------------------------------------------------------------------------
 *
 *  Usage:
 *
 *     mcstats GRAPH_IN SIZE MOTIF_ID [OUT_PREFIX]
 *
 *     GRAPH_IN   - GML format file of input graph
 *     SIZE       - Size of the motifs to consider
 *     MOTIF_ID   - The isomorphic class of the motif
 *     OUT_PREFIX - Prefix to output all clustering type and node map files (Optional)
 *
 *------------------------------------------------------------------------------------------------
 *
 *  Copyright (C) 2013 Thomas E. Gorochowski <tom@chofski.co.uk>
 *
 *  This software released under the Open Source Initiative (OSI) approved Non-Profit Open 
 *  Software License ("Non-Profit OSL") 3.0. This software is distributed in the hope that 
 *  it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *===============================================================================================*/


#include <igraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define TRUE -1
#define FALSE 0

/* ---------------------------------------------------------------------------------------------- */

/* Function prototypes */
int motif_clustering_stats (igraph_t *G, igraph_t *M, char *prefix);
int clean_subgraph(igraph_t *res, igraph_t *G, igraph_t *M, igraph_vector_t *m1Nodes, igraph_vector_t *m2Nodes);
int add_cluster_type (igraph_vector_ptr_t *cTypes, igraph_t *M, igraph_vector_t *m1, igraph_vector_t *m2);
int merge_motifs (igraph_t *res, igraph_t *M, igraph_vector_t *m1, igraph_vector_t *m2);
void print_usage (void);

/* ---------------------------------------------------------------------------------------------- */

/* Main function */
int main (int argc, const char * argv[])
{
	FILE *gFile;
	igraph_t G, M;
	
	/* Check that there are enough arguments */	
	if (argc == 2 && strcmp(argv[1], "-h") == 0) {
		print_usage();
		return 0;
	}
	if (argc < 4 || argc > 5) {
		printf("Invalid number of arguments.\n");
		return 1;
	}
	
	/* Load the user specified topology (GML Format) */
	gFile = fopen(argv[1], "r");
	igraph_read_graph_gml(&G, gFile);
	fclose(gFile);
	
	/* Generate a graph of the motif we are interested in - use isomorphic class ID */
	igraph_isoclass_create(&M, atoi(argv[2]), atoi(argv[3]), igraph_is_directed(&G));
	
	if (argc == 5) {
		/* We need to output the clustering types in graphs */
		motif_clustering_stats(&G, &M, (char *)argv[4]);
	}
	else {
		/* No graph outputs */
		motif_clustering_stats(&G, &M, NULL);
	}
		
	/* Free used memory and return */
	igraph_destroy(&G);
	igraph_destroy(&M);
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Claculate motif clustering statistics */
int motif_clustering_stats (igraph_t *G, igraph_t *M, char *prefix)
{
	igraph_integer_t i, j, k, s, t, p, i2, j2, k2, overlap, mapsCount, actMapsCount, count, found, mSize;
	int res;
	igraph_vector_t m1, m2, *curMap, *curAllMap, cTypeCounts, *curi, *curj, m1Nodes, m2Nodes, *curM;
	igraph_t subGraph;
	igraph_vs_t vs;
	igraph_bool_t iso;
	char buf[1000];
	FILE *outFile;
	igraph_vector_ptr_t cTypes, maps, actMaps, nMap;
	igraph_vector_ptr_init(&maps, 0);
	igraph_vector_ptr_init(&cTypes, 0);
	mSize = igraph_vcount(M);
	
	/* ----------------------------------------------------- */
	/* PART I: Generate all motif clustering types in a list */
	/* ----------------------------------------------------- */
	
#ifdef DEBUG
	printf("Generate all motif clustering types.\n");
	fflush(stdout);
#endif

	/* Loop through different size overlaps when clustered */
	for (overlap=1; overlap<igraph_vcount(M); overlap++) {

		if (overlap == 1) {
			/* Create new lists the required size to hold overlap mapping */
			igraph_vector_init(&m1, overlap);
			igraph_vector_init(&m2, overlap);
			/* Go through all possible permutations */
			for (i=0; i<igraph_vcount(M); i++) {
				for (i2=0; i2<igraph_vcount(M); i2++) {
					/* Create mapping and add if new */
					VECTOR(m1)[0] = i;
					VECTOR(m2)[0] = i2;
					add_cluster_type(&cTypes, M, &m1, &m2);
				}
			}
			/* Free the overlap lists */
			igraph_vector_destroy(&m1);
			igraph_vector_destroy(&m2);			
		}
		else if (overlap == 2) {
			/* Create new lists the required size to hold overlap mapping */
			igraph_vector_init(&m1, overlap);
			igraph_vector_init(&m2, overlap);
			/* Go through all possible permutations */
			for (i=0; i<igraph_vcount(M); i++) {
				for (j=0; j<igraph_vcount(M); j++) {
					if (i!=j) {
						/* ------------------- */
						for (i2=0; i2<igraph_vcount(M); i2++) {
							for (j2=0; j2<igraph_vcount(M); j2++) {
								if (i2!=j2) {
									/* Create mapping and add if new */
									VECTOR(m1)[0] = i;
									VECTOR(m2)[0] = i2;
									VECTOR(m1)[1] = j;
									VECTOR(m2)[1] = j2;
									add_cluster_type(&cTypes, M, &m1, &m2);
								}
							}
						}
						/* ------------------- */
					}
				}
			}
			/* Free the overlap lists */
			igraph_vector_destroy(&m1);
			igraph_vector_destroy(&m2);
		}
		else if (overlap == 3) {
			/* Create new lists the required size to hold overlap mapping */
			igraph_vector_init(&m1, overlap);
			igraph_vector_init(&m2, overlap);
			/* Go through all possible permutations */
			for (i=0; i<igraph_vcount(M); i++) {
				for (j=0; j<igraph_vcount(M); j++) {
					for (k=0; k<igraph_vcount(M); k++) {
						if (i!=j && i!=k && j!=k) {
							/* ------------------- */
							for (i2=0; i2<igraph_vcount(M); i2++) {
								for (j2=0; j2<igraph_vcount(M); j2++) {
									for (k2=0; k2<igraph_vcount(M); k2++) {
										if (i2!=j2 && i2!=k2 && j2!=k2) {
											/* Create mapping and add if new */
											VECTOR(m1)[0] = i;
											VECTOR(m2)[0] = i2;
											VECTOR(m1)[1] = j;
											VECTOR(m2)[1] = j2;
											VECTOR(m1)[2] = k;
											VECTOR(m2)[2] = k2;
											add_cluster_type(&cTypes, M, &m1, &m2);
										}
									}
								}
							}
							/* ------------------- */
						}
					}
				}
			}
			/* Free the overlap lists */
			igraph_vector_destroy(&m1);
			igraph_vector_destroy(&m2);		
		}
		else {
			printf("Error: there is currently only support for 3 and 4 node motifs\n");
			return 1;
		}			
	}
	
	/* Check to see if we need to output the clustering types in GML format */
	if (prefix != NULL ) {
		for (i=0; i<igraph_vector_ptr_size(&cTypes); i++) {
			sprintf(buf, "%sType%li.gml", prefix, (long int)i+1);
			outFile = fopen(buf, "w");
			igraph_write_graph_gml((igraph_t *)VECTOR(cTypes)[(long int)i], outFile, NULL, NULL);
			fclose(outFile);
		}
	}
	
#ifdef DEBUG
	printf("Found %li types of motif clustering\n", (long int)igraph_vector_ptr_size(&cTypes));
#endif
	
	/* -------------------------------------------------------- */
	/* PART II: Find motifs and do pairwise comparison to types */
	/* -------------------------------------------------------- */
	
#ifdef DEBUG
	printf("Finding motifs in graph.\n");
	fflush(stdout);
#endif
	
	/* Find submorphisms between graph and motif */
	igraph_get_subisomorphisms_vf2(G, M, NULL, NULL, NULL, NULL, &maps, NULL, NULL, NULL);
	
#ifdef DEBUG
	printf("Cleaning up motif mappings.\n");
	fflush(stdout);
#endif
	
	/* Clean up maps list (only required for directed graphs) */
	if (igraph_is_directed(G) != 0) {
		for (i=0; i<igraph_vector_ptr_size(&maps); i++) {
			curMap = (igraph_vector_t *)VECTOR(maps)[(long int)i];
			/* Generate a vertex selector from our IDs vector */
			igraph_vs_vector(&vs, curMap);
			/* Extract the subgraph using this vertex selector */
			igraph_subgraph(G, &subGraph, vs);
			/* Check to see if the motif is missing edges from original graph => not proper motif */
			if (igraph_ecount(&subGraph) != igraph_vcount(M)) {
				/* Remove mapping set first value to -1*/
				VECTOR(*curMap)[0] = -1;
			}
			/* Free used memory */
			igraph_vs_destroy(&vs);
			igraph_destroy(&subGraph);
		}
	}
	
	/* Keep track of the unique motifs found so far */
	mapsCount = igraph_vector_ptr_size(&maps);
	actMapsCount = 0;
	
	/* Somewhere to hold the actual motifs, resize later once the actual number is known */
	igraph_vector_ptr_init(&actMaps, mapsCount);
	
	/* Go through the whole list and add if not found already */
	for (i=0; i<mapsCount; i++) {
		curAllMap = (igraph_vector_t *)VECTOR(maps)[(long int)i];
		if (VECTOR(*curAllMap)[0] != -1) {
			
			found = 0;
			for (j=0; j<actMapsCount; j++) {
				curMap = (igraph_vector_t *)VECTOR(actMaps)[(long int)j];
				/* Loop through both vectors and see if matching elements are found */
				count = 0;
				for (s=0; s<igraph_vcount(M); s++) {
					for (t=0; t<igraph_vcount(M); t++) {
						if (VECTOR(*curAllMap)[(long int)s] == VECTOR(*curMap)[(long int)t]) {
							count++;
							break;
						}
					}
				}
				if (count == igraph_vcount(M)) {
					/* Found the motif, break and do not add */
					found = 1;
					break;
				}
			}
			if (found == 0) {
				/* New motif so add */
				VECTOR(actMaps)[(long int)actMapsCount] = curAllMap;
				actMapsCount++;
			}
			
		}
	}

#ifdef DEBUG
	printf("Found %li actual motif mappings in graph.\n", (long int)actMapsCount);
	fflush(stdout);
#endif
	
	/* Resize to the actual count number */
	igraph_vector_ptr_resize(&actMaps, actMapsCount);
	
	/* At this point, actMaps contains the clean list of motif mappings; we now look at all pairs
	   compare to the clustering types we generated previously */

#ifdef DEBUG
	printf("Finding all pairs of motif and comparing to types.\n");
	fflush(stdout);
#endif
	
	if (prefix != NULL ) {
		/* List to hold node IDs for each type of clustering */
		igraph_vector_ptr_init(&nMap, igraph_vector_ptr_size(&cTypes));
		for (i=0; i<igraph_vector_ptr_size(&cTypes); i++) {
			VECTOR(nMap)[(long int)i] = (igraph_vector_t *)malloc(sizeof(igraph_vector_t));
			igraph_vector_init(VECTOR(nMap)[(long int)i], 0);
		}
	}
	
	/* Create vector to hold the counts for each clustering type */
	igraph_vector_init(&cTypeCounts, igraph_vector_ptr_size(&cTypes)+1);
	igraph_vector_fill(&cTypeCounts, (igraph_integer_t)0);
	
	igraph_vector_init(&m1Nodes, mSize);
	igraph_vector_init(&m2Nodes, mSize);
	for (i=0; i<actMapsCount-1; i++) {
		curi = (igraph_vector_t *)VECTOR(actMaps)[(long int)i];
		for (j=i+(igraph_integer_t)1; j<actMapsCount; j++) {
			curj = (igraph_vector_t *)VECTOR(actMaps)[(long int)j];
			
			/* Generate the subgraph of the motifs */
			for (k=0; k<mSize; k++) {
				VECTOR(m1Nodes)[(long int)k] = VECTOR(*curi)[(long int)k];
				VECTOR(m2Nodes)[(long int)k] = VECTOR(*curj)[(long int)k];
			}
			res = clean_subgraph(&subGraph, G, M, &m1Nodes, &m2Nodes);

			/* Check to see if any clustering has occured */
			if (res != 1) {
				/* Find clustering type and increment count */
				for (k=0; k<igraph_vector_ptr_size(&cTypes); k++) {
					igraph_isomorphic(&subGraph, (igraph_t *)VECTOR(cTypes)[(long int)k], &iso);
					if ((long int)iso != (long int)0) {
						/* Found the type, update the count and break */
						VECTOR(cTypeCounts)[(long int)k] = (VECTOR(cTypeCounts)[(long int)k])+1;

						/* If outputing node maps then check if already exists and output to file */
						if (prefix != NULL ) {
							/* Check to see if the vector already contains element and add for all */
							for (p=0; p<mSize; p++) {
								if (igraph_vector_contains((igraph_vector_t *)VECTOR(nMap)[(long int)k], 
														(igraph_real_t)VECTOR(m1Nodes)[(long int)p]) == FALSE){
									igraph_vector_push_back(VECTOR(nMap)[(long int)k], 
																	(igraph_real_t)VECTOR(m1Nodes)[(long int)p]);
								}
								if (igraph_vector_contains((igraph_vector_t *)VECTOR(nMap)[(long int)k], 
														(igraph_real_t)VECTOR(m2Nodes)[(long int)p]) == FALSE){
									igraph_vector_push_back(VECTOR(nMap)[(long int)k], 
																	(igraph_real_t)VECTOR(m2Nodes)[(long int)p]);

								}
							}
						}
						
						/* Don't continue searching, break instead */
						break;
					}
				}
				igraph_destroy(&subGraph);
			}
			else {
				VECTOR(cTypeCounts)[igraph_vector_ptr_size(&cTypes)] = 
				(VECTOR(cTypeCounts)[igraph_vector_ptr_size(&cTypes)])+1;
			}
		}
	}
	
	/* Check to see if we need to output the node maps */
	if (prefix != NULL ) {
		sprintf(buf, "%sNodeMaps.txt", prefix);
		outFile = fopen(buf, "w");
		for (i=0; i<igraph_vector_ptr_size(&cTypes); i++) {
			curM = (igraph_vector_t *)VECTOR(nMap)[(long int)i];
			for (j=0; j<igraph_vector_size(curM); j++) {
				fprintf(outFile, "%li", (long int)VECTOR(*curM)[(long int)j]);
				if (j < igraph_vector_size(curM)-1) fprintf(outFile, ",");
			}
			fprintf(outFile, "\n");
		}
		fclose(outFile);
	}
	
	/* Print out the results */
	for (i=0; i<igraph_vector_size(&cTypeCounts); i++) {
		printf("%li", (long int)VECTOR(cTypeCounts)[(long int)i]);
		if (i+1<igraph_vector_size(&cTypeCounts)) printf(",");
	}
	printf("\n");
	
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Generate a clean subgraph of the motif mappings - returns 1 if no clustering */
int clean_subgraph(igraph_t *res, igraph_t *G, igraph_t *M, igraph_vector_t *m1Nodes, igraph_vector_t *m2Nodes)
{
	igraph_vector_t map;
	long int pos1, pos2;
	igraph_es_t es;
	igraph_eit_t eit;
	igraph_bool_t s1, s2;
	igraph_integer_t mSize, i, j, mapSize, eid, found;
	mSize = igraph_vector_size(m1Nodes);
	
	/* Check if any overlap - most will not so worth doing this separately */
	found = 0;
	for (i=0; i<mSize; i++) {
		for (j=0; j<mSize; j++) {
			if ((long int)VECTOR(*m2Nodes)[(long int)i] ==
				 (long int)VECTOR(*m1Nodes)[(long int)j]) {
				found++;
				break;
			}
		}
		if (found == 1) break;
	}
	if (found == 0) return 1;
	
	/* Make the full mapping vector */
	igraph_vector_init(&map, mSize*2);
	igraph_vector_fill(&map, -1.0);
	for (i=0; i<mSize; i++) {
		VECTOR(map)[(long int)i] = (igraph_real_t)VECTOR(*m1Nodes)[(long int)i];
	}
	mapSize = mSize;
	
	/* Copy the motifs */
	igraph_copy(res, M);

	for (i=0; i<mSize; i++) {
		found = 0;
		for (j=0; j<mSize; j++) {
			if ((long int)VECTOR(*m2Nodes)[(long int)i] ==
				 (long int)VECTOR(*m1Nodes)[(long int)j]) {
				found = 1;
				break;
			}
		}
		
		if (found == 0) {
			/* Add the node to graph and update mapping list */
			VECTOR(map)[(long int)mapSize] = (igraph_real_t)VECTOR(*m2Nodes)[(long int)i];
			igraph_add_vertices(res, 1, 0);
			mapSize++;
		}
	}
	/* Update the size of the vector */
	igraph_vector_resize(&map, mapSize);

	/* Add the rest of the edges from the motif */
	for (i=mSize; i<mapSize; i++) {

		/* Get all edge (IN and OUT) from the node */
		igraph_es_adj(&es, (igraph_integer_t)VECTOR(map)[(long int)i], IGRAPH_ALL);

		/* For each edge check to see if other end is in mapping list */
		igraph_eit_create(G, es, &eit);
		IGRAPH_EIT_RESET(eit);
		do {
			
			/* Find the mapping */
			eid = IGRAPH_EIT_GET(eit);		
			if(igraph_vector_contains(m2Nodes, (igraph_real_t)IGRAPH_FROM(G, eid)) != 0
				&& igraph_vector_contains(m2Nodes, (igraph_real_t)IGRAPH_TO(G, eid)) != 0) {
			
				s1 = igraph_vector_search(&map, 0, (igraph_real_t)IGRAPH_FROM(G, eid), &pos1);
				s2 = igraph_vector_search(&map, 0, (igraph_real_t)IGRAPH_TO(G, eid), &pos2);

				/* If so then use mapping to add edge */
				if (s1 != 0 && s2 != 0) {
					igraph_add_edge(res, (igraph_integer_t)pos1, (igraph_integer_t)pos2);
				}
			}
			
			IGRAPH_EIT_NEXT(eit);
		} while (!IGRAPH_EIT_END(eit));
		igraph_eit_destroy(&eit);
		igraph_es_destroy(&es);
	}

	/* Remove any duplicate edges */
	igraph_simplify(res, -1, -1, 0);

	/* Free used memory */
	igraph_vector_destroy(&map);
	
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Given a new cluster mapping, checks if exists and adds to types pointer list */
int add_cluster_type (igraph_vector_ptr_t *cTypes, igraph_t *M, igraph_vector_t *m1, igraph_vector_t *m2)
{
	igraph_t *G, *curG, tempG;
	igraph_integer_t i, morphs, found;
	igraph_vs_t vs;
	igraph_vector_t seq;
	
#ifdef DEBUG
	printf("Entering add_cluster_type\n");
	fflush(stdout);
#endif DEBUG
	
	G = (igraph_t *)malloc(sizeof(igraph_t));
	merge_motifs(G, M, m1, m2);
	
	/* Check that the clustering doesn't upset the motif type (double edges) */
	igraph_vector_init(&seq, igraph_vcount(M));
	for (i=0; i<igraph_vcount(M); i++) {
		VECTOR(seq)[(long int)i] = i;
	}
	igraph_vs_vector(&vs, &seq);
	igraph_subgraph(G, &tempG, vs);	
	if ((long int)igraph_ecount(&tempG) != (long int)igraph_ecount(M)) {
		igraph_destroy(&tempG);
		igraph_vs_destroy(&vs);
		igraph_vector_destroy(&seq);
		return 0;
	}
	igraph_destroy(&tempG);
	igraph_vs_destroy(&vs);
	
	/* Do it for the other motifs mapping as well */
	for (i=0; i<igraph_vector_size(m1); i++) {
		VECTOR(seq)[(long int)i] = VECTOR(*m1)[(long int)i];
	}
	for (i=0; i<igraph_vcount(M)-igraph_vector_size(m1); i++) {
		VECTOR(seq)[(long int)(igraph_vector_size(m1) + i)] = igraph_vcount(M) + i;
	}
	igraph_vs_vector(&vs, &seq);
	igraph_subgraph(G, &tempG, vs);	
	if ((long int)igraph_ecount(&tempG) != (long int)igraph_ecount(M)) {
		igraph_destroy(&tempG);
		igraph_vs_destroy(&vs);
		igraph_vector_destroy(&seq);
		return 0;
	}
	igraph_destroy(&tempG);
	igraph_vs_destroy(&vs);
	igraph_vector_destroy(&seq);
	
	/* See if the motif clustering type is already found */
	found = 0;
	for (i=0; i<igraph_vector_ptr_size(cTypes); i++) {
		curG = (igraph_t *)VECTOR(*cTypes)[(long int)i];
		/* Check to see if the clustering type is the same */
		if (igraph_vcount(G) == igraph_vcount(curG) && igraph_ecount(G) == igraph_ecount(curG)) {
			igraph_count_subisomorphisms_vf2(G, curG, NULL, NULL, NULL, NULL, &morphs, NULL, NULL, NULL);
			if (morphs > 0) {
				/* Must already exist */
				found = 1;
				break;
			}
		}		
	}
	
	/* Didn't find it so append to list */
	if (found == 0) {
		igraph_vector_ptr_push_back(cTypes, (void *)G);
	}

#ifdef DEBUG
	printf("Leaving add_cluster_type\n");
	fflush(stdout);
#endif DEBUG
	
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Merge the same motif twice using a set of vertex overlaps */
int merge_motifs (igraph_t *res, igraph_t *M, igraph_vector_t *m1, igraph_vector_t *m2)
{
	igraph_vector_t map, newEdges;
	igraph_integer_t i, j, remaining;
	igraph_eit_t eit;

#ifdef DEBUG
	printf("Entering merge_motifs\n");
	fflush(stdout);
#endif DEBUG
	
	/* Create the merged graph with correct number of nodes */
	igraph_copy(res, M);
	remaining = igraph_vcount(M) - igraph_vector_size(m1);
	igraph_add_vertices(res, remaining, 0);
	
	/* Create the mapping to the new graph, so that edges can be merged correctly */
	igraph_vector_init(&map, igraph_vcount(M));
	igraph_vector_fill(&map, -1.0);
	
	/* Set the known mappings from the overlap */
	for (i=0; i<igraph_vector_size(m1); i++) {
		VECTOR(map)[(long int)VECTOR(*m2)[(long int)i]] = VECTOR(*m1)[(long int)i];
	}
	
	/* Update the remaining mapping entries with new node IDs */
	j = igraph_vcount(M);
	for (i=0; i<igraph_vector_size(&map); i++) {
		/* See if item needs updating */
		if ((long int)VECTOR(map)[(long int)i] == (long int)-1) {
			VECTOR(map)[(long int)i] = j;
			j++;
		}
	}
	
	/* Use the mappings list we have to add the edges */
	igraph_eit_create(M, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
	igraph_vector_init(&newEdges, igraph_ecount(M)*2);
	i = 0;
	IGRAPH_EIT_RESET(eit);
	do {
		VECTOR(newEdges)[(long int)i] = VECTOR(map)[(long int)IGRAPH_FROM(M, IGRAPH_EIT_GET(eit))];
		VECTOR(newEdges)[(long int)i+1] = VECTOR(map)[(long int)IGRAPH_TO(M, IGRAPH_EIT_GET(eit))];
		i += 2;
		IGRAPH_EIT_NEXT(eit);
	} while (!IGRAPH_EIT_END(eit));
	igraph_add_edges(res, &newEdges, 0);
	
	/* Remove duplicate edges (this interferes with isomorphic tests later) */
	igraph_simplify(res, TRUE, FALSE, 0);
	
	/* Free used memory */
	igraph_vector_destroy(&map);
	igraph_vector_destroy(&newEdges);
	igraph_eit_destroy(&eit);
	
#ifdef DEBUG
	printf("Leaving merge_motifs\n");
	fflush(stdout);
#endif DEBUG
	
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Print usage information */
void print_usage (void)
{
	printf("mcstats GRAPH_IN SIZE MOTIF_ID [OUT_PREFIX]\n");
	printf("  GRAPH_IN   - GML format file of input graph\n");
	printf("  SIZE       - Size of the motifs to consider\n");
	printf("  MOTIF_ID   - The isomorphic class of the motif\n");
	printf("  OUT_PREFIX - Prefix to output all clustering type and nodes files (Optional)\n");
}

/* ---------------------------------------------------------------------------------------------- */

/*===============================================================================================
 *  mcextract.c
 *
 *  Extracts from an input graph a subgraph containing the motif of interest. Motifs are given
 *  via their isomorphic class from igraph. Verbose modes can be enabled at compile time 
 *  (-DDEBUG) which will help understand the steps being performed.
 *
 *------------------------------------------------------------------------------------------------
 *
 *  To compile use the following command:
 *
 *     gcc -I INC_DIR -L LIB_DIR -O3 mcextract.c -ligraph -lstdc++ -o mcextract
 *
 *  where INC_DIR is the include directory and LIB_DIR is the library directory. The igraph
 *  library is required to compile this program and can be found at http://igraph.sourceforge.net/
 *
 *------------------------------------------------------------------------------------------------
 *
 *  Usage:
 *
 *     mcextract GRAPH_IN MOTIF_SIZE MOTIF_ID GRAPH_OUT [MAP_OUT]
 *
 *     GRAPH_IN:    GML format file of input graph.
 *     MOTIF_SIZE:  Size of the motifs to consider.
 *     MOTIF_ID:    The isomorphic class of the 1st motif.
 *     GRAPH_OUT:   File to output the subgraph to (GML format).
 *     MAP_OUT:     File containing mappings of in node -> out node (optional)
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
#include <string.h>

#define TRUE -1
#define FALSE 0

/* ---------------------------------------------------------------------------------------------- */

/* Function prototypes */
int  motif_extract (const igraph_t *G, igraph_t *res, igraph_t *M, igraph_vector_t *nMaps);
void print_usage   (void);

/* ---------------------------------------------------------------------------------------------- */

/* Main function */
int main (int argc, const char * argv[])
{
	igraph_integer_t i;
	FILE *gFile, *sFile;
	igraph_t G, subgraphs, M;
	igraph_vector_t motifs, nMaps;
	igraph_vector_init(&motifs, 0);
	
	/* Check that there are enough arguments */	
	if (argc == 2 && strcmp(argv[1], "-h") == 0) {
		print_usage();
		return 0;
	}
	if (argc < 5 || argc > 6) {
		printf("Invalid number of arguments.\n");
		return 1;
	}
	
	/* Load the user specified topology (GML Format) */
	gFile = fopen(argv[1], "r");
	igraph_read_graph_gml(&G, gFile);
	fclose(gFile);
	
	/* Place motifs from command line into vector */
	igraph_vector_push_back(&motifs, atoi(argv[3]));
	
	/* Generate a graph of the motif we are interested in - use isomorphic class ID */
	igraph_isoclass_create(&M, atoi(argv[2]), atoi(argv[3]), igraph_is_directed(&G));
	
	/* Extract the subgraph */
	motif_extract (&G, &subgraphs, &M, &nMaps);

	/* Write extracted subgraph to file */
	sFile = fopen(argv[4], "w");
	igraph_write_graph_gml(&subgraphs, sFile, NULL, NULL);
	fclose(sFile);
	
	/* Output node mappings to original graph */
	if (argc == 6) {
		sFile = fopen(argv[5], "w");
		for (i=0; i<igraph_vector_size(&nMaps); i++) {
			fprintf(sFile, "%li,%li\n", (long int)i, (long int)VECTOR(nMaps)[(long int)i]);
		}
		fclose(sFile);
	}
	
	/* Free used memory and return */
	igraph_vector_destroy(&nMaps);
	igraph_destroy(&G);
	igraph_destroy(&subgraphs);
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Extract the required motifs from the graph  */
int motif_extract (const igraph_t *G, igraph_t *outG, igraph_t *M, igraph_vector_t *nMaps)
{
	igraph_integer_t i, j, s, t, mapsCount, actMapsCount, count, found, mSize;
	igraph_vector_t *curMap, *curAllMap;
	igraph_t subGraph;
	igraph_vs_t vs;
	igraph_vector_ptr_t cTypes, maps, actMaps;
	igraph_vector_ptr_init(&maps, 0);
	igraph_vector_ptr_init(&cTypes, 0);
	mSize = igraph_vcount(M);
	
	igraph_integer_t toAdd, eid;
	long int nID;
	igraph_bool_t sRes;
	igraph_vector_t newMap;
	igraph_es_t es;
	igraph_eit_t eit;
	
	igraph_empty(outG, 0, igraph_is_directed(G));
	
	
#ifdef DEBUG
	printf("Finding motifs in graph.\n");
	fflush(stdout);
#endif
	
	/* Find submorphisms between graph and motif */
	igraph_get_subisomorphisms_vf2(G, M, &maps);
	
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
	
	/* This holds the mappings from our new node ID to the old ones in G */
	igraph_vector_init(nMaps, 0);
	
	/* Some variable for holding and adding nodes and mappings */
	igraph_vector_init(&newMap, mSize);
	es = igraph_ess_all(IGRAPH_EDGEORDER_ID);
	igraph_eit_create(M, es, &eit);
	
	/* To ensure we only include edges of the motifs we grow the graph one motif at a time */
	for (i=0; i<igraph_vector_ptr_size(&actMaps); i++) {
		
		curMap = (igraph_vector_t *)VECTOR(actMaps)[(long int)i];
		toAdd = 0;
		
		/* Calculate mapping for current motif */
		for (j=0; j<mSize; j++) {
			sRes = igraph_vector_search(nMaps, 0, (igraph_real_t)VECTOR(*curMap)[(long int)j], &nID);
			if ((long int)sRes != FALSE) {
				/* Node already exists so take ID from the mapping list */
				VECTOR(newMap)[(long int)j] = (igraph_real_t)nID;
			}
			else {
				/* Need to add a new node to mapping list so use the new ID */
				VECTOR(newMap)[(long int)j] = igraph_vector_size(nMaps);
				igraph_vector_push_back(nMaps, (igraph_real_t)VECTOR(*curMap)[(long int)j]);
				toAdd++;
			}
		}
		
		/* Add missing nodes */
		igraph_add_vertices(outG, toAdd, (void *)0);
		
		/* Iterate through motifs edges using mapping to add edges to growing graph */
		IGRAPH_EIT_RESET(eit);
		do {
			/* Find the mapping */
			eid = IGRAPH_EIT_GET(eit);		
			igraph_add_edge(outG, (igraph_integer_t)VECTOR(newMap)[(long int)IGRAPH_FROM(M, eid)], 
										 (igraph_integer_t)VECTOR(newMap)[(long int)IGRAPH_TO(M, eid)]);

			
			IGRAPH_EIT_NEXT(eit);
		} while (!IGRAPH_EIT_END(eit));
	}
	igraph_eit_destroy(&eit);
	
	/* Remove any duplicate edges */
	igraph_simplify(outG, -1, -1);
	
	/* Free used memory */
	igraph_vs_destroy(&vs);
	
	return 0;
}

/* ---------------------------------------------------------------------------------------------- */

/* Print usage information */
void print_usage (void)
{
	printf("mcextract GRAPH_IN MOTIF_SIZE MOTIF_ID GRAPH_OUT [MAP_OUT]\n");
	printf("  GRAPH_IN   - GML format file of input graph\n");
	printf("  MOTIF_SIZE - Size of the motif to consider\n");
	printf("  MOTIF_ID   - The isomorphic class of the motif to extract\n");
	printf("  GRAPH_OUT  - File to output the subgraph to (GML format)\n");
	printf("  MAP_OUT    - File containing mappings of in node -> out node (optional)\n");
}

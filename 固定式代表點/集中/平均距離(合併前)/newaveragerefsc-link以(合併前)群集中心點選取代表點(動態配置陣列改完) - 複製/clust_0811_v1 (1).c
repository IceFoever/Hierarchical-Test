
// https://github.com/gyaikhom/agglomerative-hierarchical-clustering/blob/master/agglomerate.c

/**
 * Copyright 2014 Gagarine Yaikhom (MIT License)
 *
 * Implements Agglomerative Hierarchical Clustering algorithm.
 *
 * Aug. 4, 2017: modify coord from .x .y to [NUM_ATTRS]. NUM_ATTRS is 2.
 */

#define SMALL_DATA
//#define LARGE_DATA

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */
#define MAX_LABEL_LEN 16
#define NUM_ATTRS 9

#define AVERAGE_LINKAGE  'a' /* choose average distance */
#define CENTROID_LINKAGE 't' /* choose distance between cluster centroids */
#define COMPLETE_LINKAGE 'c' /* choose maximum distance */
#define SINGLE_LINKAGE   's' /* choose minimum distance */

#define alloc_mem(N, T) (T *) calloc(N, sizeof(T))
#define alloc_fail(M) fprintf(stderr,                                   \
                              "Failed to allocate memory for %s.\n", M)
#define read_fail(M) fprintf(stderr, "Failed to read %s from file.\n", M)
#define invalid_node(I) fprintf(stderr,                                 \
                                "Invalid cluster node at index %d.\n", I)

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;
typedef struct item_s item_t;
typedef struct dist_rec_s dist_rec;

// n... : new
typedef struct nnode_s nnode;
struct nnode_s {
        int num_items; /* number of items that was clustered */
        int first_item;
        int last_item;
};


float (*distance_fptr)(float **, const int *, const int *, int, int);

struct cluster_s {
        int num_items; /* number of items that was clustered */
        int num_clusters; /* current number of root clusters */
        int num_nodes; /* number of leaf and merged clusters */
        cluster_node_t *nodes; /* leaf and merged clusters */
        float **distances; /* distance between leaves */
};

struct cluster_node_s {
        int type; /* type of the cluster node */
        int is_root; /* true if cluster hasn't merged with another */
        int height; /* height of node from the bottom */
        float centroid[NUM_ATTRS]; /* centroid of this cluster */
        char *label; /* label of a leaf node */
        int *merged; /* indexes of root clusters merged */
        int num_items; /* number of leaf nodes inside new cluster */
        int *items; /* array of leaf nodes indices inside merged clusters */
        neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct neighbour_s {
        int target; /* the index of cluster node representing neighbour */
        float distance; /* distance between the nodes */
        neighbour_t *next, *prev; /* linked list entries */
};

struct item_s {
        float coord[NUM_ATTRS]; /* coordinate of the input data point */
        char label[MAX_LABEL_LEN]; /* label of the input data point */
};

struct dist_rec_s {
       int index; // nodex index; smallest dist from a certain node i to node (index)
       float dist;
};

float euclidean_distance(const float *a, const float *b)
{
        return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2));
}

void fill_euclidean_distances(float **matrix, int num_items,
                              const item_t items[])
{
   int i, j;
   for (i = 0; i < num_items; ++i) {
     matrix[i][i] = 0.0;
     for (j = i+1; j < num_items; ++j) {
        matrix[i][j] = euclidean_distance(items[i].coord, items[j].coord);
        matrix[j][i] = matrix[i][j];
     } // j loop
  }// i loop
}

float **generate_distance_matrix(int num_items, const item_t items[])
{
      int i;
        float **matrix = alloc_mem(num_items, float *);
        if (matrix) {
                for (i = 0; i < num_items; ++i) {
                        matrix[i] = alloc_mem(num_items, float);
                        if (!matrix[i]) {
                                alloc_fail("distance matrix row");
                                num_items = i;
                                for (i = 0; i < num_items; ++i)
                                        free(matrix[i]);
                                free(matrix);
                                matrix = NULL;
                                break;
                        }
                }
                if (matrix)
                        fill_euclidean_distances(matrix, num_items, items);
        } else
                alloc_fail("distance matrix");
        return matrix;
}

float single_linkage(float **distances, const int a[],
                     const int b[], int m, int n)
{
  int i, j;
        float min = FLT_MAX, d;
        for (i = 0; i < m; ++i)
                for (j = 0; j < n; ++j) {
                        d = distances[a[i]][b[j]];
                        if (d < min)
                                min = d;
                }
        return min;
}

float complete_linkage(float **distances, const int a[],
                       const int b[], int m, int n)
{
  int i, j;
        float d, max = 0.0 /* assuming distances are positive */;
        for (i = 0; i < m; ++i)
                for (j = 0; j < n; ++j) {
                        d = distances[a[i]][b[j]];
                        if (d > max)
                                max = d;
                }
        return max;
}

float average_linkage(float **distances, const int a[],
                      const int b[], int m, int n)
{
  int i, j;
        float total = 0.0;
        for (i = 0; i < m; ++i)
                for (j = 0; j < n; ++j)
                        total += distances[a[i]][b[j]];
        return total / (m * n);
}

float centroid_linkage(float **distances, const int a[],
                       const int b[], int m, int n)
{
        return 0; /* empty function */
}

float get_distance(cluster_t *cluster, int index, int target)
{
        /* if both are leaves, just use the distances matrix */
        if (index < cluster->num_items && target < cluster->num_items)
                return cluster->distances[index][target];
        else {
                cluster_node_t *a = &(cluster->nodes[index]);
                cluster_node_t *b = &(cluster->nodes[target]);
                if (distance_fptr == centroid_linkage)
                        return euclidean_distance(a->centroid, b->centroid);
                else return distance_fptr(cluster->distances,
                                          a->items, b->items,
                                          a->num_items, b->num_items);
        }
}

void free_neighbours(neighbour_t *node)
{
        neighbour_t *t;
        while (node) {
                t = node->next;
                free(node);
                node = t;
        }
}

void free_cluster_nodes(cluster_t *cluster)
{
  int i;
        for (i = 0; i < cluster->num_nodes; ++i) {
                cluster_node_t *node = &(cluster->nodes[i]);
                if (node->label)
                        free(node->label);
                if (node->merged)
                        free(node->merged);
                if (node->items)
                        free(node->items);
                if (node->neighbours)
                        free_neighbours(node->neighbours);
        }
        free(cluster->nodes);
}

void free_cluster(cluster_t * cluster)
{
  int i;
        if (cluster) {
                if (cluster->nodes)
                        free_cluster_nodes(cluster);
                if (cluster->distances) {
                        for (i = 0; i < cluster->num_items; ++i)
                                free(cluster->distances[i]);
                        free(cluster->distances);
                }
                free(cluster);
        }
}

void insert_before(neighbour_t *current, neighbour_t *neighbours,
                   cluster_node_t *node)
{
        neighbours->next = current;
        if (current->prev) {
                current->prev->next = neighbours;
                neighbours->prev = current->prev;
        } else
                node->neighbours = neighbours;
        current->prev = neighbours;
}

void insert_after(neighbour_t *current, neighbour_t *neighbours)
{
        neighbours->prev = current;
        current->next = neighbours;
}

void insert_sorted(cluster_node_t *node, neighbour_t *neighbours)
{
        neighbour_t *temp = node->neighbours;
        while (temp->next) {
                if (temp->distance >= neighbours->distance) {
                        insert_before(temp, neighbours, node);
                        return;
                }
                temp = temp->next;
        }
        if (neighbours->distance < temp->distance)
                insert_before(temp, neighbours, node);
        else
                insert_after(temp, neighbours);
}

#define init_leaf(cluster, node, item, len)             \
        do {                                            \
                strncpy(node->label, item->label, len); \
                { int i; for (i = 0; i < NUM_ATTRS; i++) \
                node->centroid[i] = item->coord[i]; }           \
                node->type = LEAF_NODE;                 \
                node->is_root = 1;                      \
                node->height = 0;                       \
                node->num_items = 1;                    \
                node->items[0] = cluster->num_nodes++;  \
        } while (0)                                     \

#undef init_leaf
void st(nnode **nodes, int *next_item,int num_items,int num_clusters_remaining,int **arr){
               int item ;
               int node_i,i;


               for(node_i =0; node_i < num_clusters_remaining ; node_i++){
                      item = nodes[node_i]->first_item ;
                       for(i = 0 ; i < nodes[node_i]->num_items; i++){
                        arr[node_i][i] = item ;
                        item = next_item[item] ;
                      //  printf("%d ",arr[node_i][i]) ;
                       }
               }
        }

void print_cluster_items(cluster_t *cluster, int index)
{
  int i;
        cluster_node_t *node = &(cluster->nodes[index]);
        fprintf(stdout, "Items: ");
        if (node->num_items > 0) {
                fprintf(stdout, "%s", cluster->nodes[node->items[0]].label);
                for (i = 1; i < node->num_items; ++i)
                        fprintf(stdout, ", %s",
                                cluster->nodes[node->items[i]].label);
        }
        fprintf(stdout, "\n");
}

void print_cluster_node(cluster_t *cluster, int index)
{
        cluster_node_t *node = &(cluster->nodes[index]);
        fprintf(stdout, "Node %d - height: %d, centroid: (%5.3f, %5.3f)\n",
                index, node->height, node->centroid[0], node->centroid[1]);
        if (node->label)
                fprintf(stdout, "\tLeaf: %s\n\t", node->label);
        else
                fprintf(stdout, "\tMerged: %d, %d\n\t",
                        node->merged[0], node->merged[1]);
        print_cluster_items(cluster, index);
        fprintf(stdout, "\tNeighbours: ");
        neighbour_t *t = node->neighbours;
        while (t) {
                fprintf(stdout, "\n\t\t%2d: %5.3f", t->target, t->distance);
                t = t->next;
        }
        fprintf(stdout, "\n");
}


/*
#define merge_to_one(cluster, to_merge, node, node_idx)         \
        do {                                                    \
                node->num_items = to_merge[0]->num_items +      \
                        to_merge[1]->num_items;                 \
                node->items = alloc_mem(node->num_items, int);  \
                if (node->items) {                              \
                        merge_items(cluster, node, to_merge);   \
                        cluster->num_nodes++;                   \
                        cluster->num_clusters--;                \
                        update_neighbours(cluster, node_idx);   \
                } else {                                        \
                        alloc_fail("array of merged items");    \
                        free(node->merged);                     \
                        node = NULL;                            \
                }                                               \
        } while(0)                                              \
*/

cluster_node_t *merge(cluster_t *cluster, int first, int second)
{
/*        int new_idx = cluster->num_nodes;
        cluster_node_t *node = &(cluster->nodes[new_idx]);
        node->merged = alloc_mem(2, int);
        if (node->merged) {
                cluster_node_t *to_merge[2] = {
                        &(cluster->nodes[first]),
                        &(cluster->nodes[second])
                };
                node->merged[0] = first;
                node->merged[1] = second;
                merge_to_one(cluster, to_merge, node, new_idx);
        } else {
                alloc_fail("array of merged nodes");
                node = NULL;
        }
        return node;
*/
}

#undef merge_to_one

void find_best_distance_neighbour(cluster_node_t *nodes,
                                  int node_idx,
                                  neighbour_t *neighbour,
                                  float *best_distance,
                                  int *first, int *second)
{
/*
        while (neighbour) {
                if (nodes[neighbour->target].is_root) {
                        if (*first == -1 ||
                            neighbour->distance < *best_distance) {
                                *first = node_idx;
                                *second = neighbour->target;
                                *best_distance = neighbour->distance;
                        }
                        break;
                }
                neighbour = neighbour->next;
        }
*/
}


cluster_t *merge_clusters(cluster_t *cluster)
{
        int first, second;
        while (cluster->num_clusters > 1) {
//                if (find_clusters_to_merge(cluster, &first, &second) != -1)
//                        merge(cluster, first, second);
        }
        return cluster;
}

#define init_cluster(cluster, num_items, items)                         \
        do {                                                            \
                cluster->distances =                                    \
                        generate_distance_matrix(num_items, items);     \
                if (!cluster->distances)                                \
                        goto cleanup;                                   \
                cluster->num_items = num_items;                         \
                cluster->num_nodes = 0;                                 \
                cluster->num_clusters = 0;                              \
                if (add_leaves(cluster, items))                         \
                        merge_clusters(cluster);                        \
                else                                                    \
                        goto cleanup;                                   \
        } while (0)                                                     \


#undef init_cluster

int print_root_children(cluster_t *cluster, int i, int nodes_to_discard)
{
  int j, t;
        cluster_node_t *node = &(cluster->nodes[i]);
        int roots_found = 0;
        if (node->type == A_MERGER) {
                for (j = 0; j < 2; ++j) {
                        t = node->merged[j];
                        if (t < nodes_to_discard) {
                                print_cluster_items(cluster, t);
                                ++roots_found;
                        }
                }
        }
        return roots_found;
}

void get_k_clusters(cluster_t *cluster, int k)
{
        if (k < 1)
                return;
        if (k > cluster->num_items)
                k = cluster->num_items;

        int i = cluster->num_nodes - 1;
        int roots_found = 0;
        int nodes_to_discard = cluster->num_nodes - k + 1;
        while (k) {
                if (i < nodes_to_discard) {
                        print_cluster_items(cluster, i);
                        roots_found = 1;
                } else
                        roots_found = print_root_children(cluster, i,
                                                          nodes_to_discard);
                k -= roots_found;
                --i;
        }
}

void print_cluster(cluster_t *cluster)
{
  int i;
        for (i = 0; i < cluster->num_nodes; ++i)
                print_cluster_node(cluster, i);
}

void print_clu_dist(float *clu_distances, int num_items, int num_clusters_remaining) {
     int i, j, ind;

     printf("cluster dist:\n    ");
     for (i = 0; i < num_clusters_remaining; i++) printf("%6d ", i);
     printf("\n");
     for (i = 0; i < num_clusters_remaining; i++) {
         ind = i * num_items;
         printf("%6d", i);
         for (j = 0; j < num_clusters_remaining; j++) {
             if (j <= i) printf("   . . ");
             else printf("%6.3f ", clu_distances[ind]);
             ind++;
         }
         printf("\n");
     }
}


void print_best_dist(dist_rec * smallest_dist, int num_clusters_remaining) {
  int i;
        printf("best from each node: (total %d nodes)\n", num_clusters_remaining);
        for (i = 0; i < num_clusters_remaining - 1; i++)
            printf("%d: %d %f\n", i, smallest_dist[i].index, smallest_dist[i].dist);
        printf("\n");
}

// the items in each node (cluster)
void print_nodes(nnode **nodes, int *next_item, int num_clusters_remaining) {
  int i, n;
  int skew,average,allitem;
  int ptr,ptr2;
  float distance ;
  float max = 0.0;
  FILE *fp;
  fp = fopen("skew.txt","a");
  average = 207 / num_clusters_remaining ;
  printf("items in node\n");
  for (i = 0; i < num_clusters_remaining; i++) {
      ptr = nodes[i]->first_item;
      n = nodes[i]->num_items;
      printf("%d \n",n) ;
      printf("node %d: ", i);
      while (n > 0) {
        printf("%d ", ptr);
        ptr = next_item[ptr];
        n--;
      }
      allitem =  nodes[i]->num_items;
      skew += abs(average - allitem) ;
     // printf("allitem %d average %d skew %d \n",allitem,average,skew) ;
     // printf("\n");

  }
  fprintf(fp,"%d \n",skew) ;
  fclose(fp) ;

}

// best_a may be 0 or last => no problem
// [0, best_a], [1, best_a] ... [best_a -1, best_a], [best_a, best_a + 1], ...
// update smallest_dist:
//   [0, best_a], [1, best_a] ... [best_a -1, best_a]:
//        compare the current smallest_dist[i]->dist and clu_dist[i, best_a]
//   [best_a, best_a + 1] ...:
//        choose the min
// clu_dist[best_a, best_a]: don't care
void link_dist(nnode **nodes, int *next_item,
               float *clu_distances, float *item_distances,
               int num_items, int best_a,
               int num_clusters_remaining, dist_rec *smallest_dist,int **arr) {
   int node_i, i, j, item_first_a, item_i, item_j,itemm;
       // item_first_a: first item of nodes[best_a]
   float clu_dist; // can be single_min or comple_max
   float dist;
   float maxdist , mindist ;
   int clu_dist_index;
    int rrr[12000] ;
   int xr ;
   // single link
   //printf("single link\n");
 //  print_nodes(nodes, next_item, num_clusters_remaining);
   for(node_i = 0;node_i < num_clusters_remaining; node_i++){
           rrr[node_i] = nodes[node_i]->num_items ;
           if(rrr[node_i] > 10){
           rrr[node_i] = 10 ;
           }
     }
   // single link
#ifdef SMALL_DATA
   //printf("single link\n");
  // print_nodes(nodes, next_item, num_clusters_remaining);
  /* for(node_i = 0 ; node_i < num_clusters_remaining ; node_i++){
             printf("choose node_i %d" , node_i) ;
        for(i =0 ; i < floor(sqrt(nodes[node_i]->num_items)) ; i++){
             printf(" %d",arr[node_i][i]) ;

        }
        printf("\n") ;
  }*/

#endif
   clu_dist_index = best_a; // incremental: num_items (before row best_a)
 //  item_first_a = nodes[best_a]->first_item ;
   item_first_a = arr[best_a][0] ;
   for (node_i = 0; node_i < best_a; node_i++) {
       // compute dist between node_i and best_a (node_i < best_a)
      // item_j = nodes[node_i]->first_item ;
       item_j = arr[node_i][0];
       if (item_first_a < item_j)
           clu_dist = item_distances[item_first_a * num_items + item_j];
       else clu_dist = item_distances[item_j * num_items + item_first_a];
       // clu_dist: min for single link
       for (j = 0; j < rrr[node_i]; j++) {
          for (i = 0, item_i = item_first_a; i < rrr[best_a]; i++) {
              item_i = arr[best_a][i] ;
              item_j = arr[node_i][j] ;
              if (item_i < item_j)
                 dist = item_distances[item_i * num_items + item_j];
              else dist = item_distances[item_j * num_items + item_i];
              if(i==0){
              maxdist = dist ;
              mindist = dist ;
              }
              if(maxdist < dist) {
               maxdist = dist ;
              }
              if(mindist > dist) {
               mindist = dist ;
              }
             //    item_i = next_item[item_i];
              //printf("   %d %d %f\n", item_i, item_j, clu_dist);
          } // i loop; best_a node
         // item_j = next_item[item_j];
       } // j loop; each node
       //printf("min dist from node %d to merged node %d: %f\n", node_i, best_a, clu_dist);
       // if (smallest_dist[node_i].index == best_a) // for single link, (node_i, best_a) will be even shorter
       // if (clu_distances[clu_dist_index] == best_b): clu_distances[clu_dist_index] != best_a
       clu_dist = 2*(maxdist*mindist)/(maxdist + mindist) ;
       if (clu_dist < smallest_dist[node_i].dist) {
          smallest_dist[node_i].index = best_a;
          smallest_dist[node_i].dist = clu_dist;
       }
       clu_distances[clu_dist_index] = clu_dist;
       clu_dist_index += num_items;
   } // for node_i; before best_a


   // printf("link computing fter best_a\n");
   // to nodes after best_a
   // clu_dist[best_a, best_a +1], [best_a, best_a +2], ...
   // compute dist between best_a and node_i (best_a < node_i)
   clu_dist_index = best_a * num_items + best_a + 1; // incremental: 1
   // first item of best_a un-changed
   for (node_i = best_a + 1; node_i < num_clusters_remaining; node_i++) {
      //  item_j = nodes[node_i]->first_item ;
       item_j = arr[node_i][0];
       if (item_first_a < item_j)
           clu_dist = item_distances[item_first_a * num_items + item_j];
       else clu_dist = item_distances[item_j * num_items + item_first_a];
       for (j = 0; j < rrr[node_i]; j++) {
          for (i = 0, item_i = item_first_a; i < rrr[best_a]; i++) {
              item_i = arr[best_a][i] ;
              item_j = arr[node_i][j] ;
              if (item_i < item_j)
                 dist = item_distances[item_i * num_items + item_j];
              else dist = item_distances[item_j * num_items + item_i];
            //  printf("best_a %d to %d: item %d %d %f %f\n", best_a, node_i, item_i, item_j, dist, clu_dist);
              if(i==0){
              maxdist = dist ;
              mindist = dist ;
              }
              if(maxdist < dist) {
               maxdist = dist ;
              }
              if(mindist > dist) {
               mindist = dist ;
              }
          //    item_i = next_item[item_i];
          } // i loop; items in best_a
       //   item_j = next_item[item_j];
       } // j loop; items in node_i
       clu_dist = 2*(maxdist*mindist)/(maxdist + mindist) ;
       if (clu_dist < smallest_dist[best_a].dist) {
          smallest_dist[best_a].index = node_i;
          smallest_dist[best_a].dist = clu_dist;
       }
       clu_distances[clu_dist_index++] = clu_dist;
      // printf("%f\n",dist);
       // clu_dist_index++;  above
   } // for node_i; after best_a

   //printf("after link computing\n");
   //print_best_dist(smallest_dist, num_clusters_remaining);
}
//////////////

void nmerge(nnode **nodes, int *next_item, int num_clusters_remaining,
                  int best_a, int best_b) {
     nnode *tmp_node_ptr;
        // merge two clusters
//        n_merge(nnode **nodes, int *next_item, best_a, best_b);
        next_item[nodes[best_a]->last_item] = nodes[best_b]->first_item;
        nodes[best_a]->last_item = nodes[best_b]->last_item;
        nodes[best_a]->num_items += nodes[best_b]->num_items;

        // move the last node the best_b
        tmp_node_ptr = nodes[num_clusters_remaining - 1];
      //  printf("nmerge %d \n",num_clusters_remaining - 1) ;
        nodes[num_clusters_remaining - 1] = nodes[best_b];  // already num_clusters_remaining--
        nodes[best_b] = tmp_node_ptr;
      //  printf("nmerge %d \n",best_b) ;
}
//////////
// quality of clustering result: quantitation error
void eval_qe(item_t *cents, nnode **nodes, item_t *items, int *next_items, int num_clusters) {
   int i, j, attr_i, item_i;
   int num_items;
   float clusterdistance[100] ;
   float sum = 0.0, dist; // dist: dist between two items
   float averagesum ;
   float db ;
   double centroiddist ;
   float clusterdist ;
   float max = 0.0 ;
   float all = 0.0 ;
   float sum2 = 0.0 ;
   FILE *fp;
   fp = fopen("end.txt","a");
   FILE *fc;
   fc = fopen("db.txt","a");
   for (i = 0; i < num_clusters; i++) {
       num_items = nodes[i]->num_items;
       item_i = nodes[i]->first_item;
       for (j = 0; j < num_items; j++) {
           dist = 0.0;
           for (attr_i = 0; attr_i < NUM_ATTRS; attr_i++) {
            if(((cents[i].coord[attr_i] - items[item_i].coord[attr_i]) * (cents[i].coord[attr_i] - items[item_i].coord[attr_i])) > 1000){
                dist += 0  ;
              }else{
              dist += ((cents[i].coord[attr_i] - items[item_i].coord[attr_i]) *
                       (cents[i].coord[attr_i] - items[item_i].coord[attr_i]));
              }
             //    printf("%d \n",item_i) ;
             //    printf("%f %f %d %f \n",cents[i].coord[attr_i],items[item_i].coord[attr_i],item_i,dist);
               //        system("pause") ;
           }
           //system("pause") ;
           sum += sqrt(dist);
           sum2 += dist;
           //printf("sum %f\n",sum) ;
           item_i = next_items[item_i];
       } // j loop
       sum2 /= num_items ;
       clusterdistance[i] = sum2 ;
   } // i loop
   //system("pause") ;
   for(i = 0; i < num_clusters; i++){
       for(j = 0; j < num_clusters; j++){
        max = 0.0 ;
        if(i!=j){
        for (attr_i = 0; attr_i < NUM_ATTRS; attr_i++) {
              // printf("%f %f \n",cents[i].coord[attr_i],cents[j].coord[attr_i]) ;
               centroiddist += sqrt((cents[i].coord[attr_i] - cents[j].coord[attr_i]) *
                               (cents[i].coord[attr_i] - cents[j].coord[attr_i]));
              // printf("%f  \n",centroiddist) ;
          }
        clusterdist = clusterdistance[i]+clusterdistance[j] ;
      //  printf("clusterdis %f centroiddist %f",clusterdist,centroiddist) ;
        db = clusterdist/centroiddist ;

        if(db > max){
          max = db ;
         }
        printf("max %f db %f \n",max ,db) ;
       // system("pause") ;
        }
       }
       all += db ;
   }
   printf("qe %f\n", sum);
   averagesum = sum / 207 ;
   fprintf(fc,"%f \n",all/num_clusters);
   fclose(fc) ;
   fprintf(fp,"%f \n",averagesum);
   fclose(fp) ;
}
//////////
void eval_centroid(item_t *cents, nnode **nodes, item_t *items, int *next_items, int num_clusters) {
   int i, j, attr_i, item_i;
   int num_items;
//   float sum;

   //printf("compute the centroids of all clusters\n");
   for (i = 0; i < num_clusters; i++)
       for (attr_i = 0; attr_i < num_clusters; attr_i++)
           cents[i].coord[attr_i] = 0.0;
   for (i = 0; i < num_clusters; i++) {
       num_items = nodes[i]->num_items;
       item_i = nodes[i]->first_item;
       //printf("\nnode %d, num_items %d, items: ", i, num_items);
       for (j = 0; j < num_items; j++) {
           //printf("  %d,", item_i);
           for (attr_i = 0; attr_i < NUM_ATTRS; attr_i++) cents[i].coord[attr_i] += items[item_i].coord[attr_i];
           item_i = next_items[item_i];
       } // j loop
       for (attr_i = 0; attr_i < NUM_ATTRS; attr_i++) cents[i].coord[attr_i] /= (float)num_items;
       //printf("\n   cent %f %f %f\n", cents[i].coord[0], cents[i].coord[1], cents[i].coord[2]);
   } // i loop
   return;
}

// the last node is merged
// may need to find a new smallest neighbor for a node (smallest_dist)
   // num_clusters_remaining: num of clusters after merge
void update_smallest_dist(dist_rec *smallest_dist,
                          float *clu_distances,
                          int num_items,
                          int num_clusters_remaining) {
  int i, j;
  int dist_index;
  int min_index;
  float min;
  for (i = 0; i < num_clusters_remaining -1; i++) {  // -1: don't care the last node
      // if the last node (merged) is the smallest dist to node i
         // if node i is next to the merged node (best_b)=> no problem
      if (smallest_dist[i].index == num_clusters_remaining) {
         dist_index = i * num_items + i + 1;
         min = clu_distances[dist_index] + 1.0;  // any large dist
         for (j = i+1; j < num_clusters_remaining; j++, dist_index++) {
             if (clu_distances[dist_index] < min) {
               min_index = j;
               min = clu_distances[dist_index];
             }
         } // j loop
         smallest_dist[i].index = min_index;
         smallest_dist[i].dist = min;
      }
  }
}

// updating cluster-to-cluster dist & smallest_dist[]
// update smallest_dist[]
// no re-compute dist from each node to the merged node, which possibly affects smallest_dist[]
// nodes[] updated by others
   // num_clusters_remaining: num of clusters after merge
void move_last_node(float *clu_distances,
                    int num_items,
                    int best_b,
                    int num_clusters_remaining,
                    dist_rec *smallest_dist) {
     int i, j;
     int dist_index, dist_index_base;
     int min_index;
     float min;
     int test ;
#ifdef SMALL_DATA
//printf("moving last node to %d ... clu dist before moving\n", best_b);
//print_clu_dist(clu_distances, num_items, num_clusters_remaining + 1);
#endif
        // choose a cluster-to-cluster dist policy
        // single link
        // every node to the last node
        // before best_b
        dist_index_base = 0;
        for (i = 0; i < best_b; i++) {
           // printf("i %d", i) ;
            test = dist_index_base + best_b ;

           // printf(" %d\n",(dist_index_base + best_b)) ;
            clu_distances[dist_index_base + best_b]
                  = clu_distances[dist_index_base + num_clusters_remaining];
        //    printf(" %f",clu_distances[dist_index_base + best_b]) ;
            if (smallest_dist[i].index == num_clusters_remaining) {
               smallest_dist[i].index = best_b; // maybe best_a
            }
            else if (smallest_dist[i].index == best_b) {
              min = clu_distances[dist_index_base + i + 1];
              min_index = i + 1;
              for (j = i + 2; j < num_clusters_remaining; j++) {
                  //printf("j %d", j) ;
                  if (clu_distances[dist_index_base + j] < min) {
                    min = clu_distances[dist_index_base + j];
                    min_index = j;
                  }
              } // j loop
              smallest_dist[i].index = min_index;
              smallest_dist[i].dist = min;
            } // else if
            dist_index_base += num_items;
        }
       //  printf("³Ó¥XÂI %d\n",best_b) ;
        // last node to position best_b (new values at clu_dist(best_b, best_b + 1), ...
        dist_index = num_items * best_b + best_b+1; // [best_b][best_b ...]
        dist_index_base = num_items * (best_b + 1) + num_clusters_remaining;
                          // clu_distances[best_b+1][num_clusters_remaining]
        min = clu_distances[dist_index_base];
        min_index = best_b + 1;  // new value at clu_dist(best_b, best_b+1)
                    // originally (best_b + 1, best_b)
        //printf("p1 (%d %d) %d %f\n", dist_index / num_items, dist_index % num_items, min_index, min);
        for (j = best_b+1; j < num_clusters_remaining; j++) {
        //    printf(" %d\n", j) ;
            clu_distances[dist_index] = clu_distances[dist_index_base];
            //printf("p2 (%d %d) %d %f\n", i, dist_index / num_items, dist_index % num_items, clu_distances[dist_index]);
            if (clu_distances[dist_index] < min) {
              min = clu_distances[dist_index];
              min_index = j;
            }
            //printf("p3 j %d, min_index %d %f\n", j, min_index, min);
            dist_index++;
            dist_index_base += num_items;
        }
        smallest_dist[best_b].index = min_index;
        smallest_dist[best_b].dist = min;
        //printf("p4 best_b %d, min %d\n", best_b, min_index);
      //  printf("³Ó¥XÂI %d\n",best_b) ;
        // rows (nodes) after best_b
        //printf("after best_b; num_clusters_remaining %d\n", num_clusters_remaining);
        for (i = best_b + 1; i < num_clusters_remaining - 1; i++) {
          if (smallest_dist[i].index == num_clusters_remaining) { // find a new smallest for row i
            dist_index = i * num_items + i + 1;
            min = clu_distances[dist_index];
            min_index = i + 1;
            //printf("i %d, %d %d\n", i, smallest_dist[i].index, num_clusters_remaining);
            for (j = i + 2; j < num_clusters_remaining; j++) {
                dist_index++;
                if (clu_distances[dist_index] < min) {
                  min = clu_distances[dist_index];
                  min_index = j;
                }
            }
            smallest_dist[i].index = min_index;
            smallest_dist[i].dist = min;
          } // if
        } // i loop

        // the new merged node, best_a

#ifdef SMALL_DATA
       // print_clu_dist(clu_distances, num_items, num_clusters_remaining);
     //   print_best_dist(smallest_dist, num_clusters_remaining);
#endif

}
///////////

void z_score(item_t *items, int num_items, int num_attrs) {
	int i, i_attr;

	float *sum_sq_or_sd, *sum_or_mean;  // later: sd, mean
	sum_sq_or_sd = (float *)malloc(sizeof(float *) * num_attrs);
	sum_or_mean = (float *)malloc(sizeof(float *) * num_attrs);
	for (i_attr = 0; i_attr < num_attrs; i_attr++) sum_sq_or_sd[i_attr] = sum_or_mean[i_attr] = 0.0;
	for (i = 0; i < num_items; ++i) {
	  item_t *t = &(items[i]);
  	  for (i_attr = 0; i_attr < num_attrs; i_attr++) {
		  sum_or_mean[i_attr] += t->coord[i_attr];
	  	  sum_sq_or_sd[i_attr] += (t->coord[i_attr] * t->coord[i_attr]);
		  // if (i_attr == 0) printf("data %f, sq = %f, sum %f\n", t->coord[i_attr], sum_sq_or_sd[0], sum_or_mean[0]);
	  }
	}
/*
	for (i_attr = 0; i_attr < num_attrs; i_attr++)
		printf("sum sq = %f, sum %f\n", sum_sq_or_sd[0], sum_or_mean[0]);
*/
	// sqrt(mean(sum_sq) - sq of mean)
	for (i_attr = 0; i_attr < num_attrs; i_attr++) {
	  sum_or_mean[i_attr] /= (float)num_items;
	  sum_sq_or_sd[i_attr] = sqrt(sum_sq_or_sd[i_attr] / (float)num_items
	  								- sum_or_mean[i_attr] * sum_or_mean[i_attr]);
	}
	//printf("sd = %f, mean %f, items %d\n", sum_sq_or_sd[0], sum_or_mean[0], num_items);
	// do normalization
	for (i = 0; i < num_items; ++i) {
	  item_t *t = &(items[i]);
  	  for (i_attr = 0; i_attr < num_attrs; i_attr++) {
		  t->coord[i_attr] = (t->coord[i_attr] - sum_or_mean[i_attr]) / sum_sq_or_sd[i_attr];
	  }
	}

/*
	printf("after norm\n");
	for (i = 0; i < num_items; ++i) {
	  item_t *t = &(items[i]);
  	  for (i_attr = 0; i_attr < num_attrs; i_attr++) {
  	  	printf("%f ", t->coord[i_attr]);
  	  }
  	  printf("\n");
    }
*/

	// system("pause");
	free(sum_or_mean);
	free(sum_sq_or_sd);
	return;
}

int read_items(int count, item_t *items, FILE *f)
{
  int i,j;
        for (i = 0; i < count; ++i) {
                 item_t *t = &(items[i]);
                for(j = 0 ; j < NUM_ATTRS ; j++){
                 fscanf(f,"%f",&(t->coord[j]));
                 }
                continue;
                read_fail("item line");
                return i;
        }
        return count;
}

int read_items_from_file(item_t **items, FILE *f)
{
        int count, r;
        r = fscanf(f, "%d\n", &count);
        if (r == 0) {
                read_fail("number of lines");
                return 0;
        }
        if (count) {
                *items = alloc_mem(count, item_t);
                if (*items) {
                        if (read_items(count, *items, f) != count)
                                free(items);
                } else
                        alloc_fail("items array");
        }
        return count;
}


int process_input(item_t **items, const char *fname)
{
        int count = 0;
        FILE *f = fopen(fname, "r");
        if (f) {
                count = read_items_from_file(items, f);
                fclose(f);
        } else
                fprintf(stderr, "Failed to open input file %s.\n", fname);
        return count;
}
double pi(int n) {
  int i ;
  srand(5);
  int count = 0;
  double x, y;
  for (i = 0; i < n; ++i) {
    x = (double) rand() / RAND_MAX;
    y = (double) rand() / RAND_MAX;
    if (x * x + y * y <= 1) ++count;
  }
  return (double) count / n * 4;
}
void choose(nnode **nodes, int *next_item,
               float *clu_distances, float *item_distances,
               int num_items, int best_a,int best_b,
               int num_clusters_remaining,int rep[1000],int **arr){
        //       printf("choose %d %d \n",best_a,best_b);
               int node_i, i, j, item_first_a, item_i, item_j;
               float clu_dist; // can be single_min or comple_max
               float dist;
               int clu_dist_index;
               //////////////////////////
               int x=0 , y=0, a=0 , t=0,o=0,z=0;
               float temp = 0 ;
               int temp2 = 0,temp3 = 0 ;
               int r, k , rr;
               int a1 ;
               int b1 ;
               int sb ;
               int data[12000] ;
               int index = 0 ;
               int tt = 0 ;
               float sumbesta[12000]={0} ;
               float sumbestb[12000]={0} ;
               float sumbestaandbestb[12000] = {0} ;
               float sum1=0,sum2=0 ;
               int sumabitem[12000] = {0} ;
               int order = 0 ;
               int bestar,bestbr ;
              bestar = nodes[best_a]->num_items ;
              if(bestar > 10 ){
              bestar = 10 ;
              }
              bestbr = nodes[best_b]->num_items ;
              if(bestbr > 10){
              bestbr = 10 ;
              }

              // item_first_a = nodes[best_a]->first_item;
              // item_j = nodes[best_b]->first_item;

    for(i =0 ; i < bestar; i++){
        item_j = arr[best_a][i] ;
        for(j = 0 ; j < bestbr ; j++){
        item_i = arr[best_b][j] ;
        if (item_i < item_j){
        dist = item_distances[item_i * num_items + item_j];
      //  printf("%f \n",dist) ;
        }
        else{
        dist = item_distances[item_j * num_items + item_i];
      //  printf("%f \n",dist) ;
        }

        sum1 += dist ;
       }
       // printf(" %d %d %f\n",item_i,item_j,dist) ;
        sum1 = sum1 / bestbr ;
        sumbesta[i] = sum1 ;
         }

    for(i =0 ; i < bestbr; i++){
        item_j = arr[best_b][i] ;
        for(j = 0 ; j < bestar ; j++){
        item_i = arr[best_a][j] ;
        if (item_i < item_j){
        dist = item_distances[item_i * num_items + item_j];
        }
        else{
        dist = item_distances[item_j * num_items + item_i];
        }
        sum2 += dist ;
       }
     //   printf(" %d %d %f\n",item_i,item_j,dist) ;
        sum2 = sum2 / bestar ;
        sumbestb[i] = sum2 ;
       }


       // printf("\n") ;
    for(i = 0; i < bestar + bestbr ;i++){
      if(i < bestar){
        sumabitem[i] = arr[best_a][i] ;
        sumbestaandbestb[i] = sumbesta[i] ;
      }else
      {
        sumabitem[i] = arr[best_b][order] ;
        sumbestaandbestb[i] = sumbestb[order] ;
        order+=1 ;
      }
    //  printf(" %d %f \n", sumabitem[i], sumbestaandbestb[i]) ;

    }



   for(i =0;i < bestar + bestbr ;i++){
      for(j =i; j < bestar + bestbr;j++){
         if(sumbestaandbestb[i] > sumbestaandbestb[j]){
                temp = sumbestaandbestb[j];
                sumbestaandbestb[j] = sumbestaandbestb[i];
                sumbestaandbestb[i] = temp;
                temp2 = sumabitem[j] ;
                sumabitem[j] =sumabitem[i] ;
                sumabitem[i] = temp2;

       }
      }
   }
//printf("choose \n") ;

for(i =0;i < floor(sqrt(nodes[best_b]->num_items))+floor(sqrt(nodes[best_a]->num_items)) ;i++){
      //  printf("arrange %d %f\n",sumabitem[i],sumbestaandbestb[i]);
}


 // printf("re ") ;
  rr = bestar + bestbr ;
 // printf("\n") ;
 // printf("chooose %d \n",rr) ;
  if(rr > 10){
    rr = 10 ;
  }
 // printf("chooose floor rr %d \n",rr) ;
  for(i = 0 ; i < rr ; i++){
   rep[i] = sumabitem[i] ;

  // printf(" %d ",sumabitem[i]) ;
  }
  // printf("end \n");


}
void nmergearr(nnode **nodes,int num_clusters_remaining,int best_a,int best_b,int **arr,int rep[100]){
     int node_i , i ;
     int rr[12000] ;
     int xr ;
     //printf("nmergearr %d %d \n", best_a,best_b) ;
     for(node_i = 0;node_i < num_clusters_remaining; node_i++){
           rr[node_i] = nodes[node_i]->num_items ;
           if(rr[node_i] > 10){
           rr[node_i] = 10 ;
           }
     }
    // printf("nmergearr %d %d \n", best_a,best_b) ;
     for(node_i = 0 ; node_i < num_clusters_remaining ; node_i++){
      //       printf("node_i %d" , node_i) ;
        for(i =0 ; i < rr[node_i] ; i++){

             if(node_i == best_a){
             arr[node_i][i] = rep[i] ;
             }
             if(node_i == best_b){
             arr[node_i][i] = arr[num_clusters_remaining][i] ;
             }
        //     printf(" %d",arr[node_i][i]) ;

        }
       // printf("\n") ;
       }

}
void dunns_index(nnode **nodes, int *next_item, float *item_distances,
               int num_items,int num_clusters_remaining, dist_rec *smallest_dist,float *clu_distances){
        printf("dunn's \n");
        FILE *fp;
        fp = fopen("dunns.txt","a");
        int i,j,z,ind ;
        int item_i,item_j,node_i ;
        float maxdistance[50] ;
        float clusterdistance ;
        float dunn;
        float alldistance[50][50] ;
        float okdistance[50][50] ;
        float dd ;
        float min ;
        float sum,allsum ;
        int ii,jj ;
        float max ;
        for(node_i =0; node_i < num_clusters_remaining ; node_i++){
                item_j = nodes[node_i]->first_item ;
                for(j = 0 ; j < nodes[node_i]->num_items ; j++){
                item_i = nodes[node_i]->first_item ;
              //  printf("%d ",item_j) ;
                for(i = 0 ; i < nodes[node_i]->num_items ; i++){
              //  printf("%d ",item_i) ;
                if(item_i < item_j){
                        max = item_distances[item_i * num_items + item_j];
                }else{
                        max = item_distances[item_j * num_items + item_i];
                }
                if(j == 0 ){
                maxdistance[node_i] = max;
                }else if(maxdistance[node_i] < max)
                {
                maxdistance[node_i] = max ;
                }


                item_i = next_item[item_i] ;
                }
                item_j = next_item[item_j] ;
        }
               // printf("cluster %d %f \n",node_i, maxdistance[node_i]) ;
               // printf("\n") ;
        }
        //system("pause") ;
        for(i =0; i < num_clusters_remaining; i++){
            for(j = 0; j < num_clusters_remaining ; j++){
                   if(i!=j&&i!=num_clusters_remaining-1){
                   item_i = nodes[i]->first_item ;
                   for(ii = 0 ; ii < nodes[i]->num_items ; ii++){
                       item_j = nodes[j]->first_item ;
                       for(jj = 0; jj < nodes[j]->num_items ; jj++){
                        if(item_i < item_j){
                        dd = item_distances[item_i * num_items + item_j] ;
                        }else{
                         dd = item_distances[item_j * num_items + item_i] ;
                         }
                         if(ii == 0){
                         min = dd ;
                         }
                         if(min > dd){
                         min = dd ;
                         }
                        //printf("one %d %d %f\n", item_i,item_j,dd) ;
                        item_j = next_item[item_j] ;
                      }
                      item_i = next_item[item_i] ;
                    }
                    alldistance[i][j] = min ;
                    if(j < i){
                        alldistance[i][j] = alldistance[j][i] ;
                    }
                    if(j == i){
                        alldistance[i][j] = 0.0 ;
                    }
                    // printf("all %d %d %f \n",i,j,alldistance[i][j]) ;
                   // printf(" min %f\n" , min) ;
                   }
             }

        }
        for(i = 0 ; i < num_clusters_remaining ; i++){
           for (j = 0 ; j < num_clusters_remaining ; j++){
                if(i==j){
                        okdistance[i][j] = 0.0;
                }else if(j < i){
                        okdistance[i][j] = alldistance[j][i] ;
                }else{
                okdistance[i][j] = alldistance[i][j] ;
                }
               // printf("%f \n",okdistance[i][j]) ;
           }

        }
        for(i = 0 ; i < num_clusters_remaining; i++){
           //sum = 0.0 ;
           for(j = 0; j < num_clusters_remaining; j++){

                dunn = okdistance[i][j]/maxdistance[i] ;
                if(okdistance[i][j]==0 || maxdistance[i]==0){
                dunn = 0 ;
                }
                sum+=dunn ;

              //  fprintf(fp,"%f %f\n", okdistance[i][j],maxdistance[i]) ;
              //  printf("%d %f\n", i , dunn) ;
           }
            // fprintf(fp,"cluster %d %f \n",i,sum) ;

        }
          fprintf(fp,"%f \n",sum ) ;
          fclose(fp) ;
}
void sp(item_t *cents, nnode **nodes, int num_clusters){
     int i ,j ,attribute;
     float dist= 0.0 ;
     float alldist = 0.0 ;
     float sum= 0.0 ;
     float sp=0.0 ;
     FILE *fp;
     fp = fopen("sp.txt","a");

    /* for(i = 0; i < num_clusters ; i++){
        for(attribute = 0 ; attribute < NUM_ATTRS; attribute++){
               printf("%f ",cents[i].coord[attribute]) ;
               }
        printf("\n") ;
     }
     system("pause") ;*/
     for(i = 0; i < num_clusters ; i++){
        sum = 0.0 ;
        dist = 0.0 ;
        for(j = 0 ; j < num_clusters ; j++){
            if(i != j){
            for(attribute = 0 ; attribute < NUM_ATTRS; attribute++){
            dist = (cents[i].coord[attribute] - cents[j].coord[attribute]) *
                       (cents[i].coord[attribute] - cents[j].coord[attribute]);
            if(dist > 100){
                dist = 0.0 ;
            }
            alldist += sqrt(dist) ;
               }
             //  printf("%f %f %f \n",cents[i].coord[attribute],cents[j].coord[attribute],dist) ;
               sum += alldist ;
             //  system("pause") ;
               }
             }
          sum /= (num_clusters-1) ;
          sp += sum ;
     }
     sp /= num_clusters ;
     fprintf(fp,"%f \n",sp/num_clusters);
     fclose(fp) ;



}
void sc(nnode **nodes, int num_clusters,float *item_distances, int *next_item,int num_items  ){
   int i,j,k,k2 ;
   int item_i,item_j ;
   double samesum ;
   float samedist ;
   float same[20][10000] ;
   double same2dist ;
   double same2sum ;
   float same2[20][10000] ;
   float allsum = 0.0;
   int time = 0 ;
   float si ;
   float all= 0.0 ;
   float end= 0.0 ;
   FILE *fp;
   fp = fopen("sc.txt","a");
   for(i = 0; i < num_clusters; i++){
       item_i = nodes[i]->first_item ;
       for(j = 0; j < nodes[i]->num_items ; j++){
       item_j = nodes[i]->first_item ;
          for(k = 0 ; k < nodes[i]->num_items ; k++){
               if(item_i < item_j){
                         samedist = item_distances[item_i * num_items + item_j] ;
                        }else{
                         samedist = item_distances[item_j * num_items + item_i] ;
                         }

             samesum += samedist ;
             item_j = next_item[item_j] ;
             }
        samesum /= (nodes[i]->num_items - 1 ) ;//////item - 同cluster item distance
        if(nodes[i]->num_items==1){
        samesum = 0.0 ;
        }
        same[i][j] = samesum ;
      //  printf("%f \n",same[i][j]) ;
        item_i = next_item[item_i] ;
       }

   }
   //////////////////////這段迴圈有問題
   for(i = 0 ; i < num_clusters; i++){
        item_i = nodes[i]->first_item ;
        for(j = 0; j < nodes[i]->num_items; j++){
            for(k = 0 ; k < num_clusters; k++){
               item_j = nodes[k]->first_item ;
               if(i!= k){
               for(k2 = 0; k2 < nodes[k]->num_items; k2++ ){
                if(item_i < item_j){
                         same2dist = item_distances[item_i * num_items + item_j] ;
                        }else{
                         same2dist = item_distances[item_j * num_items + item_i] ;
                         }
                 item_j = next_item[item_j] ;
                 same2sum += same2dist ;
                }
               }
               same2sum /=  (nodes[k]->num_items);
               allsum += same2sum ;
            }
           // printf("allsum %f \n") ;
            same2[i][j] = allsum ; ///// object to different cluster object
            item_i = next_item[item_i] ;
           // system("pause") ;
          }

   }
      for(i = 0; i < num_clusters; i++){
        all = 0.0 ;
        for(j = 0; j < nodes[i]->num_items; j++){
      //     printf("same %f same2 %f \n",same[i][j],same2[i][j]) ;
           if(same[i][j] < same2[i][j]){
               si = 1 - (same[i][j]/same2[i][j]) ;
           }
           if(same[i][j] && same2[i][j]==0){
               si = 0 ;
           }
           if(same[i][j] > same2[i][j]){
               si =  (same2[i][j]/same[i][j])-1 ;
           }
         //  printf("si %f\n",si) ;
           all += si ;
           all /= (nodes[i]->num_items) ;
        }
       // printf("all %f \n",all) ;
       // system("pause") ;
        end += all ;
      }

    //   printf("end %f \n",end) ;
    //   system("pause") ;
       end /= num_clusters ;
    //   printf("end %f \n",end) ;
       fprintf(fp,"%f \n",end);
       fclose(fp) ;

}
int main(int argc, char **argv)
{
    int z ;
    for(z = 1 ; z <= 30 ; z++){
    int i, j, k;
    item_t *items = NULL;
    item_t *centroids;  // clustering result
    int num_items;
    float *item_distances, *clu_distances;
    int dist_index, dist_index_base; // mapping 1D to 2D; maybe data type long int
    dist_rec * smallest_dist; // smallest dist from a node i to other nodes j, j > i.
                         // store the node index
    int min_index;
    float min;

    int best_a, best_b;
    int num_clusters_remaining, num_clusters;  // num_clusters: final number of clusters

    float dist_sum; // for computing distance
    double result = pi(1e8);
    clock_t start, end;
    double diff ;
    int **arr, *pData;
    int node_i ;
    int rep[1000] ;
    char filename[200] ;
//        if (argc != 4) {
        if (argc < 1) {
                fprintf(stderr, "Usage: %s <input file> <num clusters> "
                        "<linkage type>\n", argv[0]);
                exit(1);
        } else {
//                int num_items = process_input(&items, argv[1]);
#ifdef SMALL_DATA
                sprintf(filename, "%d.txt" , z) ;
                num_items = process_input(&items,filename);
#else
                num_items = process_input(&items, "large.txt");
#endif
                // z-score normalize
//                z_score(&items, num_items, num_attrs);
                z_score(items, num_items, NUM_ATTRS);
//                set_linkage(argv[3][0]);
        }

/*
        printf("items:\n");
        for (i = 0; i < num_items; i++) {
            for (k = 0; k < NUM_ATTRS; k++)
                printf("%5f ", items[i].coord[k]);
            printf("\n");
        }
        printf("\n");
*/
        num_clusters = 8;
        printf("set num_clusters %d\n", num_clusters);
        num_clusters_remaining = num_items;

        centroids = (item_t *)malloc(num_clusters * sizeof(item_t));
        if (!centroids) {
          printf("fail in malloc cntroid\n");
          system("pause");
          exit(1);
        }

        // item to item distance
        item_distances = (float *)malloc((num_items * num_items)* sizeof(float *));
        clu_distances = (float *)malloc((num_items * num_items)* sizeof(float *));
        dist_index_base = 0;
        for (i = 0; i < num_items; i++) {
            dist_index = dist_index_base + i + 1;
            for (j = i+1; j < num_items; j++) {
                dist_sum = 0.0;
                for (k = 0; k < NUM_ATTRS; k++) {
                    dist_sum += (items[i].coord[k] - items[j].coord[k])
                                * (items[i].coord[k] - items[j].coord[k]);
                }
                // store
                item_distances[dist_index++] = dist_sum;
            } // j loop
            dist_index_base += num_items;
        } // i loop

/*
        printf("item dist:\n");
        dist_index = 0;
        for (i = 0; i < num_items; i++) {
          for (j = 0; j <= i; j++) printf("%5f ", item_distances[dist_index++]);
          for (; j < num_items; j++) printf("%5f ", item_distances[dist_index++]);
          printf("\n");
        }
*/

        // initialize cluster-to-cluster distances
        dist_index_base = 0;
        for (i = 0; i < num_items; i++) {
            dist_index = dist_index_base + i + 1;
            for (j = i+1; j < num_items; j++, dist_index++) {
                clu_distances[dist_index] = item_distances[dist_index];
            } // j loop
            dist_index_base += num_items;
        } // i loop

        // initialize the nodes
        nnode **nodes = (nnode **)malloc(num_items * sizeof(nnode *));
        nnode *nodes_ = (nnode *)malloc(num_items * sizeof(nnode));
        arr = (int **)malloc(num_items*sizeof(int *)+num_items*num_items*sizeof(int));
        for (i = 0, pData = (int *)(arr+num_items); i < num_items; i++, pData += num_items){
        arr[i]=pData;
        }
        //        nnode *tmp_node_ptr; // be used in merging
        memset(nodes,'0',num_items * sizeof(nnode *)) ;
        memset(nodes_ ,'0',num_items * sizeof(nnode)) ;
        int *next_item = (int *)malloc(num_items * sizeof(int));
        for (i = 0; i < num_items; i++) {
            nodes_[i].first_item = nodes_[i].last_item = i;
            nodes_[i].num_items = 1;
            nodes[i] = &nodes_[i];
            next_item[i] = -1;
            }

        // smallest dist from a node i to other nodes j, j > i
        smallest_dist = (dist_rec *)malloc((num_items)* sizeof(dist_rec));
        dist_index_base = 0;
        for (i = 0; i < num_items-1; i++) {
            dist_index = dist_index_base + i + 1;
            min = clu_distances[dist_index];
            min_index = i+1;
            dist_index++;
            for (j = i+2; j < num_items; j++, dist_index++) {
                if (clu_distances[dist_index] < min) {
                   min = clu_distances[dist_index];
                   min_index = j;
                }
            } // j loop
            // store the min dist node index
            smallest_dist[i].index = min_index;
            smallest_dist[i].dist = min;
            dist_index_base += num_items;
        } // i loop

        //printf("best of each: ");
        //for (i = 0; i < num_items-1; i++) printf("%d, ", smallest_dist[i].index);
        //printf("\n");
        //print_best_dist(smallest_dist, num_clusters_remaining);

#ifdef SMALL_DATA
       // print_best_dist(smallest_dist, num_clusters_remaining);
       // print_clu_dist(clu_distances, num_items, num_clusters_remaining);

#endif
       start = clock() ;
       st(nodes,next_item,num_items,num_clusters_remaining,arr) ;
    /*   for(node_i = 0 ; node_i < num_clusters_remaining ; node_i++){
            for(i = 0 ; i < nodes[node_i]->num_items; i++){
             printf("node_i %d %d",node_i,arr[node_i][i]) ;
            }
            printf("\n") ;
        }*/
  while (num_clusters_remaining > num_clusters) {  // loop for a merge
        // find best pair
        best_a = 0;
        best_b = smallest_dist[0].index;
        min = clu_distances[best_b];  // dist_index = best_b
        for (i = 1; i < num_clusters_remaining -1; i++) {
          dist_index = i * num_items + smallest_dist[i].index;
          // printf("(%d %d %5f, min %5f\n", i, smallest_dist[i].index, clu_distances[dist_index], min);
          if (clu_distances[dist_index] < min) {
            min = clu_distances[dist_index];
            best_a = i;
            best_b = smallest_dist[i].index;
          }
        }
#ifdef SMALL_DATA
       // printf("\n") ;
      //  printf("best %d %d\n", best_a, best_b);
#endif
//best_a= 1; best_b = 2;
        choose(nodes,next_item,clu_distances,item_distances,num_items,best_a,best_b,num_clusters_remaining,rep,arr);
        nmerge(nodes, next_item, num_clusters_remaining, best_a, best_b);
        num_clusters_remaining--;
        nmergearr(nodes,num_clusters_remaining,best_a,best_b,arr,rep) ;
#ifdef SMALL_DATA
       // print_nodes(nodes, next_item, num_clusters_remaining);
       // printf("\n");

#endif

        if (best_b == num_clusters_remaining) { // the last node is merged
           // update smallest_dist if the smallest dist to a node is best_b (last node)
          update_smallest_dist(smallest_dist, clu_distances,
                               num_items, num_clusters_remaining);
#ifdef SMALL_DATA
  //        print_best_dist(smallest_dist, num_clusters_remaining);
#endif
        }
        else
          move_last_node(clu_distances, num_items,
                         best_b, num_clusters_remaining, smallest_dist);
          // no re-compute dist to the merged node

        // compute dist from each node to the merged node,
        // which possibly affects smallest_dist[]

        link_dist(nodes, next_item, clu_distances, item_distances, num_items, best_a, num_clusters_remaining, smallest_dist, arr);

/*
        // updating cluster-to-cluster dist & smallest_dist[]
        // choose a cluster-to-cluster dist policy
        // single link
        // every node to the last node
        // before best_b
        dist_index_base = 0;
        for (i = 0; i < best_b; i++) {
            clu_distances[dist_index_base + best_b]
                  = clu_distances[dist_index_base + num_clusters_remaining];
            if (smallest_dist[i].index == num_clusters_remaining)
               smallest_dist[i].index = best_b; // maybe best_a
            dist_index_base += num_items;
        }

//        print_clu_dist(clu_distances, num_items, num_clusters_remaining);
        // after best_b
        dist_index = num_items * best_b + best_b+1; // [best_b][best_b ...]
        dist_index_base = num_items * (best_b + 1) + num_clusters_remaining;
                          // clu_distances[best_b+1][num_clusters_remaining]
        min = clu_distances[dist_index_base];
        min_index = best_b + 1;  // ?
        printf("after:  %d %d %f\n", dist_index, min_index, min);
        for (j = best_b+1; j < num_clusters_remaining; j++) {
            clu_distances[dist_index] = clu_distances[dist_index_base];
            printf(" %d %d %f\n", i, dist_index, clu_distances[dist_index]);
            if (clu_distances[dist_index] < min) {
              min = clu_distances[dist_index];
              min_index = j;
            }
            printf("j %d, min_index %d %f\n", j, min_index, min);
            dist_index++;
            dist_index_base += num_items;
        }
        smallest_dist[best_b].index = min_index;
        smallest_dist[best_b].dist = min;
        printf("after:  _b %d,  min %d\n", best_b, min_index);
        print_clu_dist(clu_distances, num_items, num_clusters_remaining);
        print_best_dist(smallest_dist, num_clusters_remaining);

        // the new merged node, best_a
*/


  } // while loop for a merge

      end = clock() ;

/*
    pointers in (nnode **nodes) will be moving:
       when two clusters merge, the cluster with larger index (in nnode **nodes)
       will be discarded.
       Node can be re-used in next experiment. So, no need to free it until end of program.

       Example: 100 items
       For example, after some iterations, there are clusters 1 ~ 10.
       3 and 7 merge, then 10 will be new 7.

    maintaining items:
      assuming cluster 3: first item is 3, last item is 20, total 2 items.
               cluster 7 has only 1 item 7.
               next_item[20] = first of cluster 7.
               last of cluster 3 is last of cluster 7.
               total items in new cluster 3 is 2+1

    maintaining cluster to cluster distance:
      for each node in the (nnode **nodes), say node 2, cluster-to-cluster distances
      (2, 3), (2, 4) ... until the last cluster
      Initially malloc a cluster-to-cluster distance matrix (100*100)
      (an upper triangle is ok, but too complicated).

    Updating the matrix:
      distance in (i, 10) moves to (i, 7).
      re-compute distance (i, 3), because new items added to cluster 3.
*/

      printf("clustering result:\n");
      print_nodes(nodes, next_item, num_clusters);
      print_clu_dist(clu_distances, num_items, num_clusters_remaining);
      dunns_index(nodes, next_item, item_distances,num_items,num_clusters_remaining,smallest_dist,clu_distances) ;
      eval_centroid( centroids, nodes, items, next_item, num_clusters);
      sp(centroids,nodes,num_clusters) ;
      sc(nodes,num_clusters,item_distances,next_item,num_items) ;
      eval_qe(centroids, nodes, items, next_item, num_clusters);
      free(arr) ;
      free(nodes);
      free(nodes_);
      free(next_item);
      free(clu_distances);
      free(item_distances);
      free(centroids);
      printf("PI = %f\n", result);
      //system("pause");
      printf(" %f  sec", diff / CLOCKS_PER_SEC) ;
 }
      return 0;
}

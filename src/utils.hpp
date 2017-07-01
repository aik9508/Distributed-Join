#include <assert.h>
#include "Relation.hpp"

int mpi_hash(int a);

void print_array(int *v, int sz);
void print_vector(vector<int> v);
void delete_vector(vector<int *> v);
vector<int> flatten_vector(vector<int *> &v0, int arity);
void count_to_displ(int *count, int *displ, int length);
void merge_tuples(int final_dim, int *src, int *tgt, int *&rst_tuple,
                  int *&src_permu, int &src_l, int *&tgt_permu, int &tgt_l,
                  int &nj);
// Flatten the dimensionnal vectors and decide send_counts for different
// processors
void merge_tuples(int arity1, int *tuple1, int *permu1, int arity2, int *tuple2,
                  int *permu2, int &nj);

#include "utils.hpp"

int mpi_hash(int a) {
  a = (a ^ 61) ^ (a >> 16);
  a = a + (a << 3);
  a = a ^ (a >> 4);
  a = a * 0x27d4eb2d;
  a = a ^ (a >> 15);
  return a;
}

void print_array(int *v, int sz) {
  for (int i = 0; i < sz; i++) {
    cout << v[i] << " ";
  }
  cout << endl;
}

void print_vector(vector<int> v) {
  for (int i = 0; i < v.size(); i++) {
    cout << v[i] << " ";
  }
  cout << endl;
}

void delete_vector(vector<int *> v) {
  for (unsigned int i = 0; i < v.size(); i++) {
    delete[] v[i];
  }
}

/*
* Transforms a 2-dimension vectro into 1-dimension.
*/
vector<int> flatten_vector(vector<int *> &v0, int arity) {
  vector<int> v;
  for (vector<int *>::iterator it = v0.begin(); it != v0.end(); it++) {
    for (int i = 0; i < arity; i++) {
      v.push_back(*((*it) + i));
    }
  }
  return v;
}

void count_to_displ(int *count, int *displ, int length) {
  displ[0] = 0;
  for (unsigned int i = 1; i < length; i++) {
    displ[i] = displ[i - 1] + count[i - 1];
  }
}

/*
* Reconstructs the received flattened vector into 2-dimension.
*/

void merge_tuples(int final_dim, int *src, int *tgt, int *&rst_tuple,
                  int *&src_permu, int &src_l, int *&tgt_permu, int &tgt_l,
                  int &nj) {
  vector<int> src_comm_entries, src_only_entries, tgt_comm_entries,
      tgt_only_entries;
  int src_ind = 0, tgt_ind = 0;
  rst_tuple = new int[final_dim];
  for (unsigned int i = 0; i < final_dim; i++) {
    if (src[i] != -1) src_ind++;
    if (tgt[i] != -1) tgt_ind++;
  }
  int new_ind = src_ind;
  for (unsigned int i = 0; i < final_dim; i++) {
    rst_tuple[i] = src[i];
    if (src[i] != -1 && tgt[i] != -1) {
      src_comm_entries.push_back(src[i]);
      tgt_comm_entries.push_back(tgt[i]);
    } else if (src[i] != -1) {
      src_only_entries.push_back(src[i]);
    } else if (tgt[i] != -1) {
      rst_tuple[i] = new_ind++;
      tgt_only_entries.push_back(tgt[i]);
    }
  }
  nj = src_comm_entries.size();
  for (unsigned int i = 0; i < src_only_entries.size(); i++)
    src_comm_entries.push_back(src_only_entries[i]);
  for (unsigned int i = 0; i < tgt_only_entries.size(); i++)
    tgt_comm_entries.push_back(tgt_only_entries[i]);
  src_l = src_ind;
  tgt_l = tgt_ind;
  src_permu = new int[src_comm_entries.size()];
  tgt_permu = new int[tgt_comm_entries.size()];
  for (unsigned int i = 0; i < src_comm_entries.size(); i++) {
    src_permu[i] = src_comm_entries[i];
  }
  for (unsigned int i = 0; i < tgt_comm_entries.size(); i++) {
    tgt_permu[i] = tgt_comm_entries[i];
  }
}

// Flatten the dimensionnal vectors and decide send_counts for different
// processors
void merge_tuples(int arity1, int *tuple1, int *permu1, int arity2, int *tuple2,
                  int *permu2, int &nj) {
  int tuple[arity1 + arity2];
  nj = 0;
  for (unsigned int i = 0; i < arity1 + arity2; i++) {
    tuple[i] = 0;
  }
  for (unsigned int i = 0; i < arity1; i++) {
    tuple[tuple1[i]]++;
  }
  for (unsigned int i = 0; i < arity2; i++) {
    tuple[tuple2[i]]++;
  }
  unsigned int head = 0, tail = arity1 - 1;
  for (unsigned int i = 0; i < arity1; i++) {
    if (tuple[tuple1[i]] == 1) {
      permu1[tail--] = i;
    } else if (tuple[tuple1[i]] == 2) {
      permu1[head++] = i;
      nj++;
    } else {
      // Scilently
    }
  }
  assert(head = tail + 1);

  head = 0;
  tail = arity2 - 1;
  for (unsigned int i = 0; i < arity2; i++) {
    if (tuple[tuple2[i]] == 1) {
      permu2[tail--] = i;
    } else if (tuple[tuple2[i]] == 2) {
      permu2[head++] = i;
    } else {
      // Scilently
    }
  }
  assert(head = tail + 1);
}

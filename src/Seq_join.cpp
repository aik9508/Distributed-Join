#include <ctime>
#include "math.h"
#include "utils.hpp"

int main(int argc, char **argv) {
  int *permu1, *permu2;
  int arity1, arity2, nj;

  clock_t begin;
  begin = clock();
  cout << "[0.0000] Begins to construct inital data" << endl;
  unsigned int i = 1;
  string const file1 = argv[i++];
  Relation r1(file1);
  arity1 = atoi(argv[i++]);
  int *tuple1 = new int[arity1];
  for (unsigned int j = 0; j < arity1; j++) {
    tuple1[j] = atoi(argv[i++]);
  }
  string const file2 = argv[i++];
  Relation r2 = (file1 == file2) ? r1 : Relation(file2);
  arity2 = atoi(argv[i++]);
  int *tuple2 = new int[arity2];
  for (unsigned int j = 0; j < arity2; j++) {
    tuple2[j] = atoi(argv[i++]);
  }
  permu1 = new int[arity1];
  permu2 = new int[arity2];
  merge_tuples(arity1, tuple1, permu1, arity2, tuple2, permu2, nj);
  delete[] tuple1;
  delete[] tuple2;

  printf("[%6.4f]Begins to join...\n",
         (double)(clock() - begin) / CLOCKS_PER_SEC);

  Relation *res = new Relation(r1, r2, permu1, permu2, nj, true, NONE);

  printf("[%6.4f]Finishes joining and finds %d joins...\n",
         (double)(clock() - begin) / CLOCKS_PER_SEC, res->size());
  delete[] permu1;
  delete[] permu2;
  delete res;
}

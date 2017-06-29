#include <ctime>
#include "math.h"
#include "mpi.h"
#include "utils.hpp"

vector<vector<int *> > distribute_vector(Relation &rel, int const numtasks,
                                         int *permu) {
  vector<vector<int *> > distr;
  for (unsigned int i = 0; i < numtasks; i++) {
    vector<int *> tmp;
    distr.push_back(tmp);
  }
  for (vector<int *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++) {
    int dest = mpi_hash(*((*it) + permu[0])) % numtasks;
    distr[dest].push_back(*it);
  }
  return distr;
}

// Flatten the dimensionnal vectors and decide send_counts for different
// processors
void get_scatter_v(Relation &rel, vector<int> &scatter_v, int const numtasks,
                   int *permu, int *displ, int *send_counts) {
  int arity = rel.get_arity();
  if (rel.size() == 0) return;
  // Here produces 'std::out_of_range' error
  vector<vector<int *> > distr = distribute_vector(rel, numtasks, permu);
  for (unsigned int i = 0; i < numtasks; i++) {
    send_counts[i] = distr[i].size() * arity;
    for (unsigned int j = 0; j < distr[i].size(); j++) {
      for (unsigned int k = 0; k < arity; k++) {
        scatter_v.push_back(
            distr[i][j][k]);  // processor i, line j, k-th integer
      }
    }
  }
  count_to_displ(send_counts, displ, numtasks);
}

int main(int argc, char **argv) {
  const int root = 0;
  int numtasks, taskid;
  MPI_Status status;
  int *permu1, *permu2, *send_displ1, *send_displ2, *send_counts1,
      *send_counts2;
  int arity1, arity2, loc_count1, loc_count2, nj;
  vector<vector<int> > v1;
  vector<vector<int> > v2;
  vector<int> scatter_v1;
  vector<int> scatter_v2;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  clock_t begin;
  if (taskid == root) {
    begin = clock();
    cout << "[0.0000] the root task begins to construct inital data" << endl;
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

    send_displ1 = new int[numtasks];
    send_displ2 = new int[numtasks];
    send_counts1 = new int[numtasks];
    send_counts2 = new int[numtasks];
    printf("[%6.4f]construct scatter vector...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    get_scatter_v(r1, scatter_v1, numtasks, permu1, send_displ1, send_counts1);
    get_scatter_v(r2, scatter_v2, numtasks, permu2, send_displ2, send_counts2);
    printf("[%6.4f]broadcast parameters...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  }
  MPI_Bcast(&arity1, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&arity2, 1, MPI_INT, root, MPI_COMM_WORLD);
  if (taskid != root) {
    permu1 = new int[arity1];
    permu2 = new int[arity2];
  }
  MPI_Bcast(permu1, arity1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(permu2, arity2, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&nj, 1, MPI_INT, root, MPI_COMM_WORLD);
  if (taskid == root)
    printf("[%6.4f]scatter data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  MPI_Scatter(taskid == root ? send_counts1 : NULL, 1, MPI_INT, &loc_count1, 1,
              MPI_INT, root, MPI_COMM_WORLD);
  MPI_Scatter(taskid == root ? send_counts2 : NULL, 1, MPI_INT, &loc_count2, 1,
              MPI_INT, root, MPI_COMM_WORLD);
  // cout << "OK" << endl;
  // cout << loc_count1 << " " << loc_count2 << endl;
  int *loc_data1 = new int[loc_count1];
  int *loc_data2 = new int[loc_count2];
  MPI_Scatterv(scatter_v1.data(), send_counts1, send_displ1, MPI_INT, loc_data1,
               loc_count1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Scatterv(scatter_v2.data(), send_counts2, send_displ2, MPI_INT, loc_data2,
               loc_count2, MPI_INT, root, MPI_COMM_WORLD);
  printf("node%d starts to join...\n", taskid);

  Relation *loc_r1 = new Relation(loc_data1, arity1, loc_count1);
  Relation *loc_r2 = new Relation(loc_data2, arity2, loc_count2);
  Relation *res =
      new Relation(*loc_r1, *loc_r2, permu1, permu2, nj, true, NONE);
  delete loc_r1;
  delete loc_r2;
  printf("node%d finished joining and found %d joins...\n", taskid,
         res->size());
  vector<int> v_joined = flatten_vector(res->dataptr, res->get_arity());
  delete res;
  int join_count = v_joined.size();
  int recv_counts[numtasks];
  MPI_Gather(&join_count, 1, MPI_INT, recv_counts, 1, MPI_INT, root,
             MPI_COMM_WORLD);
  int *array_gather;
  int *displs;
  int total_count = 0;
  if (taskid == root) {
    displs = new int[numtasks];
    for (unsigned int i = 0; i < numtasks; i++) {
      displs[i] = total_count;
      total_count += recv_counts[i];
    }
    printf("[%6.4f]recv_counts gather finished...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    array_gather = new int[total_count];
  }
  MPI_Gatherv(v_joined.data(), join_count, MPI_INT, array_gather, recv_counts,
              displs, MPI_INT, root, MPI_COMM_WORLD);
  v_joined.clear();
  if (taskid == root) {
    printf("[%6.4f]vector gather finished...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

    Relation res_gather(array_gather, (arity1 + arity2 - nj), total_count);
    printf("[%6.4f]final result: %d\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, res_gather.size());
    cout << res_gather << endl;
  }
  delete[] permu1;
  delete[] permu2;
  if (taskid == root) {
    delete[] send_displ1;
    delete[] send_displ2;
    delete[] send_counts1;
    delete[] send_counts2;
  }
  MPI_Finalize();
}

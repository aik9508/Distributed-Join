#include <ctime>
#include "Relation.hpp"
#include "mpi.h"

vector<int> flatten_vector(vector<vector<int> *> &v0) {
  vector<int> v;
  for (vector<vector<int> *>::iterator it = v0.begin(); it != v0.end(); it++) {
    for (int i = 0; i < (*it)->size(); i++) {
      v.push_back((*it)->at(i));
    }
  }
  return v;
}

vector<vector<vector<int> *> > distribute_vector(Relation &rel,
                                                 int const numtasks,
                                                 int *permu) {
  vector<vector<vector<int> *> > distr;
  for (unsigned int i = 0; i < numtasks; i++) {
    vector<vector<int> *> tmp;
    distr.push_back(tmp);
  }
  for (vector<vector<int> *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++) {
    int dest = (*it)->at(permu[0]) % numtasks;
    distr[dest].push_back(*it);
  }
  return distr;
}

void get_scatter_v(Relation &rel, vector<int> &scatter_v, int const numtasks,
                   int const arity, int *permu, int *displ, int *send_counts) {
  vector<vector<vector<int> *> > distr =
      distribute_vector(rel, numtasks, permu);
  for (unsigned int i = 0; i < numtasks; i++) {
    send_counts[i] = distr[i].size() * arity;
    for (unsigned int j = 0; j < distr[i].size(); j++) {
      for (unsigned int k = 0; k < arity; k++) {
        scatter_v.push_back(distr[i][j]->at(k));
      }
    }
  }
  displ[0] = 0;
  for (unsigned int i = 1; i < numtasks; i++) {
    displ[i] = displ[i - 1] + send_counts[i - 1];
  }
}

void reconstruct_loc_vector(vector<vector<int> > &v, int *loc_data,
                            int const count, int const arity) {
  for (unsigned int i = 0; i < count; i += arity) {
    vector<int> current_vector;
    for (unsigned int j = 0; j < arity; j++) {
      current_vector.push_back(loc_data[i + j]);
      // cout << loc_data[i * arity + j] << " ";
    }
    // cout << endl;
    v.push_back(current_vector);
  }
  // cout << endl;
}

void print_vector(vector<int> v) {
  for (int i = 0; i < v.size(); i++) {
    cout << v[i] << " ";
  }
  cout << endl;
}

void print_array(int *v, int sz) {
  for (int i = 0; i < sz; i++) {
    cout << v[i] << " ";
  }
  cout << endl;
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
    string const file1 = argv[1];
    string const file2 = argv[2];
    permu1 = new int[2];
    permu2 = new int[2];
    permu1[0] = 1;
    permu1[1] = 0;
    permu2[0] = 0;
    permu2[1] = 1;
    arity1 = 2;
    arity2 = 2;
    nj = 1;
    Relation r1(file1);
    Relation r2 = file1==file2 ? r1 : Relation(file2);
    send_displ1 = new int[numtasks];
    send_displ2 = new int[numtasks];
    send_counts1 = new int[numtasks];
    send_counts2 = new int[numtasks];
    printf("[%6.4f]construct scatter vector...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    get_scatter_v(r1, scatter_v1, numtasks, arity1, permu1, send_displ1,
                  send_counts1);
    get_scatter_v(r2, scatter_v2, numtasks, arity2, permu2, send_displ2,
                  send_counts2);
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
  MPI_Scatter(send_counts1, 1, MPI_INT, &loc_count1, 1, MPI_INT, root,
              MPI_COMM_WORLD);
  MPI_Scatter(send_counts2, 1, MPI_INT, &loc_count2, 1, MPI_INT, root,
              MPI_COMM_WORLD);
  int loc_data1[loc_count1];
  int loc_data2[loc_count2];
  MPI_Scatterv(scatter_v1.data(), send_counts1, send_displ1, MPI_INT, loc_data1,
               loc_count1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Scatterv(scatter_v2.data(), send_counts2, send_displ2, MPI_INT, loc_data2,
               loc_count2, MPI_INT, root, MPI_COMM_WORLD);
  printf("node%d starts to join...\n", taskid);
  reconstruct_loc_vector(v1, loc_data1, loc_count1, arity1);
  reconstruct_loc_vector(v2, loc_data2, loc_count2, arity2);
  Relation loc_r1(v1);
  Relation loc_r2(v2);
  Relation res(loc_r1, loc_r2, permu1, permu2, nj, true, NONE);
  printf("node%d finished joining and found %d joins...\n", taskid, res.size());
  vector<int> v_joined = flatten_vector(res.dataptr);
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
  if (taskid == root) {
    printf("[%6.4f]vector gather finished...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    vector<vector<int> > v_gather;
    reconstruct_loc_vector(v_gather, array_gather, total_count,
                           (arity1 + arity2 - nj));
    Relation res_gather(v_gather);
    printf("[%6.4f]final result:\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    cout << res_gather << endl;
  }
  delete[] permu1;
  delete[] permu2;
  if (taskid == root) {
    delete[] array_gather;
    delete[] displs;
    delete[] send_displ1;
    delete[] send_displ2;
    delete[] send_counts1;
    delete[] send_counts2;
  }
  MPI_Finalize();
}
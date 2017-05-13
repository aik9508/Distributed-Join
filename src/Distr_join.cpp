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

void distribute_vector(Relation &rel, vector<vector<int> > &distr,
                       vector<vector<int> > &v0, int const arity,
                       int const numtasks, int *permu) {
  for (vector<vector<int> *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++) {
    int dest = (*it)->at(permu[0]) % numtasks;
    if (dest == 0) {
      v0.push_back(*(*it));
    } else {
      for (unsigned int i = 0; i < arity; i++) {
        distr[dest].push_back((*it)->at(i));
      }
    }
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
  int *permu1, *permu2;
  int arity1, arity2, loc_count1, loc_count2, nj;
  vector<vector<int> > v1;
  vector<vector<int> > v2;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (taskid == root) {
    cerr << "the root task begins to construct inital data" << endl;
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
    MPI_Bcast(&arity1, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&arity2, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(permu1, arity1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(permu2, arity2, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nj, 1, MPI_INT, root, MPI_COMM_WORLD);
    Relation r1(arity1, file1);
    Relation r2(arity2, file2);
    vector<vector<int> > distr1;
    vector<vector<int> > distr2;
    for (unsigned int i = 0; i < numtasks; i++) {
      vector<int> tmp;
      distr1.push_back(tmp);
      distr2.push_back(tmp);
    }
    distribute_vector(r1, distr1, v1, arity1, numtasks, permu1);
    distribute_vector(r2, distr2, v2, arity2, numtasks, permu2);
    for (unsigned int i = 1; i < numtasks; i++) {
      int count1 = distr1[i].size();
      int count2 = distr2[i].size();
      printf("root sends data (distr1_size:%d  distr2_size:%d) to node%d ...\n",
             count1, count2, i);
      MPI_Send(&count1, 1, MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Send(&count2, 1, MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Send(distr1[i].data(), count1, MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Send(distr2[i].data(), count2, MPI_INT, i, i, MPI_COMM_WORLD);
    }
    printf("data analyzed by root : (v1_size:%lu v2_size:%lu)...\n",
           v1.size() * arity1, v2.size() * arity2);
  } else {
    MPI_Bcast(&arity1, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&arity2, 1, MPI_INT, root, MPI_COMM_WORLD);
    permu1 = new int[arity1];
    permu2 = new int[arity2];
    MPI_Bcast(permu1, arity1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(permu2, arity2, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nj, 1, MPI_INT, root, MPI_COMM_WORLD);
    printf("node%d waits to recieve data...\n", taskid);
    MPI_Recv(&loc_count1, 1, MPI_INT, root, taskid, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Recv(&loc_count2, 1, MPI_INT, root, taskid, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    int loc_data1[loc_count1];
    int loc_data2[loc_count2];
    MPI_Recv(loc_data1, loc_count1, MPI_INT, root, taskid, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Recv(loc_data2, loc_count2, MPI_INT, root, taskid, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("node%d has recieved all data...\n", taskid);
    reconstruct_loc_vector(v1, loc_data1, loc_count1, arity1);
    reconstruct_loc_vector(v2, loc_data2, loc_count2, arity2);
  }
  Relation loc_r1(arity1, v1);
  Relation loc_r2(arity2, v2);
  Relation res = Relation::join(loc_r1, loc_r2, permu1, permu2, nj, true);
  printf("node%d finished joining...\n", taskid);
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
    cout << "recv_counts gather finished..." << endl;
    array_gather = new int[total_count];
  }
  MPI_Gatherv(v_joined.data(), join_count, MPI_INT, array_gather, recv_counts,
              displs, MPI_INT, root, MPI_COMM_WORLD);
  if (taskid == root) {
    cout << "gather finished..." << endl;
    vector<vector<int> > v_gather;
    reconstruct_loc_vector(v_gather, array_gather, total_count,
                           (arity1 + arity2 - nj));
    Relation res_gather(arity1 + arity2 - nj, v_gather);
    cout << "Final result : " << endl;
    cout << res_gather << endl;
  }
  delete[] permu1;
  delete[] permu2;
  if (taskid == root) {
    delete[] array_gather;
    delete[] displs;
  }
  MPI_Finalize();
}
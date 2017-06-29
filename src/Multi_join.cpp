#include <ctime>
#include "mpi.h"
#include "utils.hpp"

/*
* Distributes data lines to different processors according to the module of
* the key column to the number of processors.
*/
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
  int const num_files = atoi(argv[1]);  // Numbers of relations to join
  int const final_dim = atoi(argv[2]);  // Number of columns of final result
  vector<int *> tuples;                 // vector of arrays of size final_dim
  int **permus = new int *[2 * num_files - 2];
  int *counts = new int[2 * num_files - 2];
  int *arities = new int[num_files];
  int *njs = new int[num_files - 1];
  int *local_counts = new int[num_files - 1];
  int **join_order = new int *[num_files - 1];
  for (unsigned int i = 0; i < num_files - 1; i++) {
    join_order[i] = new int[2];
  }
  vector<int *> send_displs;
  vector<int *> send_counts;
  vector<vector<vector<int> > > vs;
  vector<vector<int> > scatter_vs;
  for (unsigned int i = 0; i < num_files; i++) {
    vector<vector<int> > loc_v;
    vs.push_back(loc_v);
  }
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  clock_t begin = clock();
  if (taskid == root) {
    for (unsigned int l = 0; l < num_files; l++) {
      int *tuple = new int[final_dim];
      for (unsigned int k = 0; k < final_dim; k++) {  // Initialize tuples
        tuple[k] = -1;
      }
      tuples.push_back(tuple);
    }
    for (unsigned int l = 0; l < num_files; l++) {
      vector<int> scatter_v;
      int *send_displ = new int[numtasks];
      int *send_count = new int[numtasks];
      scatter_vs.push_back(scatter_v);
      send_displs.push_back(send_displ);
      send_counts.push_back(send_count);
    }
    vector<Relation> relations;
    printf("[%6.4f]the root process begins to construct inital data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    // Create all relations from files and calculate arities and tuples
    ifstream file;
    file.open(argv[3]);
    if (!file.is_open()) {
      cout << "No such file!" << endl;
      exit(2);
    }
    for (unsigned int i = 0; i < num_files - 1; i++) {
      file >> join_order[i][0];
      file >> join_order[i][1];
    }
    string files[num_files];
    for (unsigned int i = 0; i < num_files; i++) {
      string relationfile;
      file >> relationfile;
      files[i] = relationfile;
      bool b = false;
      Relation *r;
      for (unsigned int j = 0; j < i; j++) {
        if (relationfile == files[j]) {
          r = &(relations[j]);
          b = true;
          break;
        }
      }
      if (!b) r = new Relation(relationfile);
      relations.push_back(*r);
      file >> arities[i];
      for (unsigned int k = 0; k < arities[i]; k++) {
        int index;
        file >> index;
        tuples[i][index] = k;
      }
    }
    file.close();

    // calculate all permutations of each step
    for (unsigned int i = 0; i < num_files - 1; i++) {
      int *rst_tuple;
      int src = join_order[i][0];
      int tgt = join_order[i][1];
      merge_tuples(final_dim, tuples[src], tuples[tgt], rst_tuple, permus[src],
                   counts[src], permus[tgt], counts[tgt], njs[i]);
      tuples.push_back(rst_tuple);
    }

    printf("[%6.4f]construct scatter vector...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

    for (unsigned int l = 0; l < num_files; l++) {
      get_scatter_v(relations[l], scatter_vs[l], numtasks, permus[l],
                    send_displs[l], send_counts[l]);
    }
    printf("[%6.4f]broadcast parameters...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  }

  // Distribute tasks to different nodes
  MPI_Bcast(arities, num_files, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(njs, num_files - 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(counts, 2 * num_files - 2, MPI_INT, root, MPI_COMM_WORLD);
  if (taskid != root) {
    for (unsigned int i = 0; i < 2 * num_files - 2; i++) {
      permus[i] = new int[counts[i]];
    }
  }
  if (taskid == root) {
    printf("[%6.4f]broadcast permutations...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  }
  for (unsigned int i = 0; i < num_files - 1; i++) {
    MPI_Bcast(join_order[i], 2, MPI_INT, root, MPI_COMM_WORLD);
  }
  for (unsigned int i = 0; i < 2 * num_files - 2; i++) {
    MPI_Bcast(permus[i], counts[i], MPI_INT, root, MPI_COMM_WORLD);
  }

  if (taskid == root)
    printf("[%6.4f]scatter data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

  for (unsigned int l = 0; l < num_files; l++) {
    MPI_Scatter(taskid == root ? send_counts[l] : NULL, 1, MPI_INT,
                &local_counts[l], 1, MPI_INT, root, MPI_COMM_WORLD);
  }

  vector<int *> local_datas;
  for (unsigned int l = 0; l < num_files; l++) {
    int *local_data = new int[local_counts[l]];
    local_datas.push_back(local_data);
  }
  for (unsigned int l = 0; l < num_files; l++) {
    if (taskid == root) {
      MPI_Scatterv(scatter_vs[l].data(), send_counts[l], send_displs[l],
                   MPI_INT, local_datas[l], local_counts[l], MPI_INT, root,
                   MPI_COMM_WORLD);
    } else {
      MPI_Scatterv(NULL, NULL, NULL, MPI_INT, local_datas[l], local_counts[l],
                   MPI_INT, root, MPI_COMM_WORLD);
    }
  }
  // Unitl now, all processors have received all files' data needed
  printf("[%6.4f]node%d starts to join......\n",
         (double)(clock() - begin) / CLOCKS_PER_SEC, taskid);
  // Local treatments
  Relation *loc_relations[num_files * 2 - 2];
  for (unsigned int l = 0; l < num_files; l++) {
    loc_relations[l] =
        new Relation(local_datas[l], arities[l], local_counts[l]);
  }

  Relation *res = loc_relations[0];
  vector<vector<vector<int> > > v_inters;
  vector<vector<int> > scatter_v_inters;
  for (unsigned int i = 0; i < num_files - 1; i++) {
    vector<int> sv_inter;
    vector<vector<int> > v_inter;
    v_inters.push_back(v_inter);
    scatter_v_inters.push_back(sv_inter);
  }
  int loc_count_inter;
  int recv_count_inter[numtasks], recv_displ_inter[numtasks],
      send_displ_inter[numtasks], send_count_inter[numtasks];
  for (unsigned int i = 0; i < num_files - 1; i++) {
    int src = join_order[i][0];
    int tgt = join_order[i][1];
    printf("[%6.4f]node%d join %d\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i);
    if (i == num_files - 2) {
      res = new Relation(*(loc_relations[src]), *(loc_relations[tgt]),
                         permus[src], permus[tgt], njs[i], true, TRIANGLE);
      delete loc_relations[src];
      delete loc_relations[tgt];
      break;
    } else {
      res = new Relation(*(loc_relations[src]), *(loc_relations[tgt]),
                         permus[src], permus[tgt], njs[i], true, FRIENDS);
    }
    delete loc_relations[src];
    delete loc_relations[tgt];
    printf("[%6.4f]node%d join %d get_scatter_v\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i);
    get_scatter_v(*(res), scatter_v_inters[i], numtasks, permus[i + num_files],
                  send_displ_inter, send_count_inter);
    printf("[%6.4f]node%d join %d All to all communication starts\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i);
    MPI_Alltoall(send_count_inter, 1, MPI_INT, recv_count_inter, 1, MPI_INT,
                 MPI_COMM_WORLD);
    printf("[%6.4f]node%d join %d All to all communication finished\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i);
    count_to_displ(recv_count_inter, recv_displ_inter, numtasks);
    loc_count_inter =
        recv_displ_inter[numtasks - 1] + recv_count_inter[numtasks - 1];
    int *loc_data_inter = new int[loc_count_inter];
    printf("[%6.4f]node%d join %d All to all v communication\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i);
    MPI_Alltoallv(scatter_v_inters[i].data(), send_count_inter,
                  send_displ_inter, MPI_INT, loc_data_inter, recv_count_inter,
                  recv_displ_inter, MPI_INT, MPI_COMM_WORLD);
    res = new Relation(loc_data_inter, counts[i + num_files], loc_count_inter);
    printf("[%6.4f]node%d join %d found %d intermediate results\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, i, res->size());
    loc_relations[num_files + i] = res;
  }
  printf("[%6.4f]node%d finished joining and found %d joins...\n",
         (double)(clock() - begin) / CLOCKS_PER_SEC, taskid, res->size());
  vector<int> v_joined = flatten_vector(res->dataptr, res->get_arity());
  int join_count = v_joined.size();
  int recv_counts[numtasks];

  // Gather results
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
    Relation res_gather(array_gather, final_dim, total_count);
    printf("[%6.4f]final result (%d joins) :\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, res_gather.size());
    cout << res_gather << endl;
  }
  delete[] arities;
  delete[] njs;
  delete[] counts;
  delete[] local_counts;
  for (unsigned int i = 0; i < num_files - 1; i++) {
    delete[] join_order[i];
    delete[] permus[i];
    delete[] permus[i + num_files - 1];
  }
  delete[] permus;
  delete[] join_order;
  if (taskid == root) {
    delete[] displs;
    delete_vector(tuples);
  }
  delete_vector(send_displs);
  delete_vector(send_counts);
  MPI_Finalize();
}

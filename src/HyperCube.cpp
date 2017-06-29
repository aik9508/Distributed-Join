#include <ctime>
#include "mpi.h"
#include "utils.hpp"

void distribute_vector_(vector<vector<int *> > &distr, int *to_insert, int *m,
                        int *weight, vector<int> &var_pos, int dest,
                        int count) {
  if (count == var_pos.size() - 1) {
    for (unsigned int i = 0; i < m[var_pos[count]]; i++) {
      distr[dest].push_back(to_insert);
      dest += weight[var_pos[count]];
    }
  } else {
    for (unsigned int i = 0; i < m[var_pos[count]]; i++) {
      distribute_vector_(distr, to_insert, m, weight, var_pos, dest, count + 1);
      dest += weight[var_pos[count]];
    }
  }
}

vector<vector<int *> > distribute_vector(Relation &rel, int const numtasks,
                                         int *m, int const len, int *pos) {
  // print_array(pos,len);

  vector<vector<int *> > distr;
  int weight[len];
  weight[0] = 1;
  for (unsigned int i = 1; i < len; i++) {
    weight[i] = m[i - 1] * weight[i - 1];
  }
  vector<int> var_pos;
  vector<int> invar_pos;
  for (unsigned int i = 0; i < len; i++) {
    if (pos[i] == -1) {
      var_pos.push_back(i);
    } else {
      invar_pos.push_back(i);
    }
  }
  // print_vector(var_pos);
  // print_vector(invar_pos);
  for (unsigned int i = 0; i < numtasks; i++) {
    vector<int *> tmp;
    distr.push_back(tmp);
  }
  // int count =0;
  for (vector<int *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++) {
    int p[len];
    for (unsigned int i = 0; i < var_pos.size(); i++) {
      p[var_pos[i]] = 0;
    }
    int dest = 0;
    for (unsigned int i = 0; i < invar_pos.size(); i++) {
      p[invar_pos[i]] =
          mpi_hash(*((*it) + pos[invar_pos[i]])) % m[invar_pos[i]];
    }
    for (unsigned int i = 0; i < len; i++) {
      dest += p[i] * weight[i];
    }
    // if(count<20){
    //   print_vector(**it);
    //   print_array(p,len);
    //   cout << dest << endl;
    //   count++;
    // }
    distribute_vector_(distr, *it, m, weight, var_pos, dest, 0);
  }
  return distr;
}

void get_scatter_v(Relation &rel, vector<int> &scatter_v, int const numtasks,
                   int const arity, int *m, int const len, int *pos, int *displ,
                   int *send_counts) {
  vector<vector<int *> > distr = distribute_vector(rel, numtasks, m, len, pos);
  for (unsigned int i = 0; i < numtasks; i++) {
    send_counts[i] = distr[i].size() * arity;
    for (unsigned int j = 0; j < distr[i].size(); j++) {
      for (unsigned int k = 0; k < arity; k++) {
        scatter_v.push_back(distr[i][j][k]);
      }
    }
  }
  displ[0] = 0;
  for (unsigned int i = 1; i < numtasks; i++) {
    displ[i] = displ[i - 1] + send_counts[i - 1];
  }
  distr.clear();
}

int main(int argc, char **argv) {
  const int root = 0;
  int numtasks, taskid;
  int const num_files = atoi(argv[1]);  // Numbers of relations to join
  int const final_dim = atoi(argv[2]);  // Number of columns of final result
  vector<int *> tuples;                 // vector of arrays of size final_dim
  int **permus = new int *[2 * num_files - 2];
  int *counts = new int[2 * num_files - 2];
  int *arities = new int[num_files];
  int *njs = new int[num_files - 1];
  int *local_counts = new int[num_files];
  int **join_order = new int *[num_files - 1];
  for (unsigned int i = 0; i < num_files - 1; i++) {
    join_order[i] = new int[2];
  }
  vector<int> scatter_vs;
  vector<vector<vector<int> > > vs;
  vector<int *> local_datas;
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
    printf("[%6.4f]the root process begins to construct inital data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

    // read vector m from argv
    int m[final_dim];
    int prod = 1;
    for (unsigned int i = 0; i < final_dim; i++) {
      m[i] = atoi(argv[i + 3]);
      prod *= m[i];
    }
    if (prod != numtasks) {
      cout << "the number of tasks doesn't match with the size of the hypercube"
           << endl;
      exit(2);
    }

    // read the commande file
    vector<Relation> relations;
    ifstream file;
    file.open(argv[3 + final_dim]);
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
      cout << relationfile << endl;
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

    // construct permutations
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
    for (unsigned int i = 0; i < num_files; i++) {
      int send_displ[numtasks];
      int send_count[numtasks];
      get_scatter_v(relations[i], scatter_vs, numtasks, arities[i], m,
                    final_dim, tuples[i], send_displ, send_count);
      MPI_Scatter(send_count, 1, MPI_INT, &local_counts[i], 1, MPI_INT, root,
                  MPI_COMM_WORLD);
      int *local_data = new int[local_counts[i]];
      local_datas.push_back(local_data);
      MPI_Scatterv(scatter_vs.data(), send_count, send_displ, MPI_INT,
                   local_datas[i], local_counts[i], MPI_INT, root,
                   MPI_COMM_WORLD);
      scatter_vs.clear();
      print_array(send_count, numtasks);
    }
    printf("[%6.4f]scatter data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  } else {
    for (unsigned int i = 0; i < num_files; i++) {
      MPI_Scatter(NULL, 1, MPI_INT, &local_counts[i], 1, MPI_INT, root,
                  MPI_COMM_WORLD);
      int *local_data = new int[local_counts[i]];
      local_datas.push_back(local_data);
      MPI_Scatterv(NULL, NULL, NULL, MPI_INT, local_datas[i], local_counts[i],
                   MPI_INT, root, MPI_COMM_WORLD);
    }
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

  // if (taskid == root)
  //   printf("[%6.4f]scatter data...\n",
  //          (double)(clock() - begin) / CLOCKS_PER_SEC);

  printf("node%d starts to join...\n", taskid);
  Relation *loc_relations[num_files * 2 - 1];
  for (unsigned int l = 0; l < num_files; l++) {
    loc_relations[l] =
        new Relation(local_datas[l], arities[l], local_counts[l]);
  }
  Relation *res = loc_relations[0];
  for (unsigned int i = 0; i < num_files - 1; i++) {
    int src = join_order[i][0];
    int tgt = join_order[i][1];
    if (i == num_files - 2) {
      res = new Relation(*(loc_relations[src]), *(loc_relations[tgt]),
                         permus[src], permus[tgt], njs[i], true, TRIANGLE);
    } else {
      res = new Relation(*(loc_relations[src]), *(loc_relations[tgt]),
                         permus[src], permus[tgt], njs[i], true, FRIENDS);
    }
    loc_relations[num_files + i] = res;
    delete loc_relations[join_order[i][0]];
    delete loc_relations[join_order[i][1]];
    printf("node%d join %d finished (%d intermediate results)\n", taskid, i,
           res->size());
  }
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
  MPI_Finalize();
}

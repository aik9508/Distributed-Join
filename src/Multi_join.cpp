#include <assert.h>
#include <ctime>
#include "Relation.hpp"
#include "mpi.h"

/*
* Transforms a 2-dimension vectro into 1-dimension.
*/
vector<int> flatten_vector(vector<vector<int> *> &v0)
{
  vector<int> v;
  for (vector<vector<int> *>::iterator it = v0.begin(); it != v0.end(); it++)
  {
    for (int i = 0; i < (*it)->size(); i++)
    {
      v.push_back((*it)->at(i));
    }
  }
  return v;
}
/*
* Distributes data lines to different processors according to the module of
* the key column to the number of processors.
*/
vector<vector<vector<int> *> > distribute_vector(Relation &rel,
                                                int const numtasks,
                                                int *permu)
{
  vector<vector<vector<int> *> > distr;
  for (unsigned int i = 0; i < numtasks; i++)
  {
    vector<vector<int> *> tmp;
    distr.push_back(tmp);
  }
  for (vector<vector<int> *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++)
  {
    int dest = (*it)->at(permu[0]) % numtasks;
    distr[dest].push_back(*it);
  }
  return distr;
}

void count_to_displ(int *count, int *displ, int length)
{
  displ[0] = 0;
  for (unsigned int i = 1; i < length; i++)
  {
    displ[i] = displ[i - 1] + count[i - 1];
  }
}

// Flatten the dimensionnal vectors and decide send_counts for different
// processors
void get_scatter_v(Relation &rel, vector<int> &scatter_v, int const numtasks,
                   int *permu, int *displ, int *send_counts)
{
  int arity = rel.get_arity();
  if (rel.size() == 0)
    return;

  // Here produces 'std::out_of_range' error
  vector<vector<vector<int> *> > distr = distribute_vector(rel, numtasks, permu);
  for (unsigned int i = 0; i < numtasks; i++)
  {
    send_counts[i] = distr[i].size() * arity;
    for (unsigned int j = 0; j < distr[i].size(); j++)
    {
      for (unsigned int k = 0; k < arity; k++)
      {
        scatter_v.push_back(
            distr[i][j]->at(k)); // processor i, line j, k-th integer
      }
    }
  }
  count_to_displ(send_counts, displ, numtasks);
}
/*
* Reconstructs the received flattened vector into 2-dimension.
*/
void reconstruct_loc_vector(vector<vector<int> > &v, int *loc_data,
                            int const count, int const arity)
{
  for (unsigned int i = 0; i < count; i += arity)
  {
    vector<int> current_vector;
    for (unsigned int j = 0; j < arity; j++)
    {
      current_vector.push_back(loc_data[i + j]);
      // cout << loc_data[i * arity + j] << " ";
    }
    // cout << endl;
    v.push_back(current_vector);
  }
  // cout << endl;
}

void print_array(int *v, int sz)
{
  for (int i = 0; i < sz; i++)
  {
    cout << v[i] << " ";
  }
  cout << endl;
}

void merge_tuples(int final_dim, int* rst, int *tgt, int* &rst_permu, int &rst_l, 
                  int* &tgt_permu, int & nj){
  vector<int> rst_comm_entries,rst_only_entries,tgt_comm_entries,tgt_only_entries;
  int rst_ind = 0, tgt_ind = 0;
  for(unsigned int i = 0; i<final_dim; i++){
    if(rst[i]!=-1 && tgt[i]!=-1){
      rst_comm_entries.push_back(rst_ind++);
      tgt_comm_entries.push_back(tgt_ind++);
    }
    else if(rst[i]!=-1) {
      rst_only_entries.push_back(rst_ind++);
    }
    else if(tgt[i]!=-1) {
      rst[i]=i;
      tgt_only_entries.push_back(tgt_ind++);
    }
  }
  nj = rst_comm_entries.size();
  for(unsigned int i =0 ;i<rst_only_entries.size();i++)
    rst_comm_entries.push_back(rst_only_entries[i]);
  for(unsigned int i =0 ;i<tgt_only_entries.size();i++)
    tgt_comm_entries.push_back(tgt_only_entries[i]);
  rst_l = rst_ind;
  rst_permu = new int[rst_comm_entries.size()];
  tgt_permu = new int[tgt_comm_entries.size()];
  for(unsigned int i = 0;i<rst_comm_entries.size();i++){
    rst_permu[i]=rst_comm_entries[i];
  }
  for(unsigned int i = 0;i<tgt_comm_entries.size();i++){
    tgt_permu[i]=tgt_comm_entries[i];
  }
}

void print_vector(vector<int> v)
{
  for (int i = 0; i < v.size(); i++)
  {
    cout << v[i] << " ";
  }
  cout << endl;
}

int main(int argc, char **argv)
{
  const int root = 0;
  int numtasks, taskid;
  MPI_Status status;
  int const num_files = atoi(argv[1]); // Numbers of relations to join
  int const final_dim = atoi(argv[2]); // Number of columns of final result
  vector<int *> tuples;                // vector of arrays of size final_dim
  vector<int *> rslt_permus;           // read from arguments
  vector<int *> inter_permus;
  int *arities = new int[num_files];
  int *njs = new int[num_files-1];
  int *rslt_counts = new int[num_files];
  int *local_counts = new int[num_files];
  vector<int *> send_displs;
  vector<int *> send_counts;
  vector<vector<vector<int> > > vs;
  vector<vector<int> > scatter_vs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  clock_t begin;
  // cout << "num_files = " << num_files << endl;
  for (unsigned int l = 0; l < num_files; l++)
  {
    vector<vector<int> > loc_v;
    vs.push_back(loc_v);
    int *tuple = new int[final_dim];
    for (unsigned int k = 0; k < final_dim; k++)
    { // Initialize tuples
      tuple[k] = -1;
    }
    tuples.push_back(tuple);
  }

  if (taskid == root)
  {
    for (unsigned int l = 0; l < num_files; l++)
    {
      vector<int> scatter_v;
      int *send_displ = new int[numtasks];
      int *send_count = new int[numtasks];
      scatter_vs.push_back(scatter_v);
      send_displs.push_back(send_displ);
      send_counts.push_back(send_count);
    }
    vector<Relation> relations;
    begin = clock();
    cout << "[0.0000] the root task begins to construct inital data" << endl;
    // Create all relations from files and calculate arities and tuples
    int index = 0;
    for (unsigned int l = 3, index = 0; l < argc; index++)
    {
      Relation *r = new Relation(argv[l++]);
      relations.push_back(*r);
      arities[index] = atoi(argv[l++]);
      for (unsigned int k = 0; k < arities[index]; k++)
      {
        int entry = atoi(argv[l++]);
        tuples[index][entry] = k;
      }
    }

    // calculate all permutations of each step
    int rst_tuple[final_dim];
    for(unsigned int i = 0; i<final_dim;i++){
      rst_tuple[i]=tuples[0][i];
    }
    for (unsigned int i = 0; i < num_files-1; i++)
    {
      int *inter_permu;
      int *rst_permu;
      merge_tuples(final_dim, rst_tuple, tuples[i+1], rst_permu, rslt_counts[i],
                   inter_permu, njs[i]);
      rslt_permus.push_back(rst_permu);
      print_array(rst_tuple,final_dim);
      print_array(rst_permu,i+2);
      print_array(inter_permu,2);
      inter_permus.push_back(inter_permu);
    }

    printf("[%6.4f]construct scatter vector...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

    for (unsigned int l = 0; l < num_files; l++)
    {
      get_scatter_v(relations.at(l), scatter_vs.at(l), numtasks,
                    l==0 ? rslt_permus[0]:inter_permus.at(l-1), send_displs.at(l), send_counts.at(l));
    }
    printf("[%6.4f]broadcast parameters...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  }

  // Distribute tasks to different nodes
  MPI_Bcast(arities, num_files, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(njs, num_files-1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(rslt_counts, num_files-1, MPI_INT, root, MPI_COMM_WORLD);
  if (taskid != root)
  {
    for (unsigned int i = 0; i < num_files-1; i++)
    {
      int *inter_permu = new int[arities[i+1]];
      int *rslt_permu = new int[rslt_counts[i]];
      rslt_permus.push_back(rslt_permu);
      inter_permus.push_back(inter_permu);
    }
  }
  if (taskid == root)
  {
    printf("[%6.4f]broadcast permutations...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
  }
  for (unsigned int i = 0; i < num_files-1; i++)
  {
    MPI_Bcast(rslt_permus[i], rslt_counts[i], MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(inter_permus[i], arities[i+1], MPI_INT, root, MPI_COMM_WORLD);
  }

  if (taskid == root)
    printf("[%6.4f]scatter data...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);

  for (unsigned int l = 0; l < num_files; l++)
  {
    MPI_Scatter(taskid == root ? send_counts[l] : NULL, 1, MPI_INT, &local_counts[l], 1, MPI_INT, root,
                MPI_COMM_WORLD);
  }

  vector<int *> local_datas;
  for (unsigned int l = 0; l < num_files; l++)
  {
    int *local_data = new int[local_counts[l]];
    local_datas.push_back(local_data);
  }
  for (unsigned int l = 0; l < num_files; l++)
  {
    if (taskid == root)
    {
      MPI_Scatterv(scatter_vs[l].data(), send_counts[l], send_displs[l], MPI_INT,
                   local_datas[l], local_counts[l], MPI_INT, root,
                   MPI_COMM_WORLD);
    }
    else
    {
      MPI_Scatterv(NULL, NULL, NULL, MPI_INT,
                   local_datas[l], local_counts[l], MPI_INT, root,
                   MPI_COMM_WORLD);
    }
  }
  // Unitl now, all processors have received all files' data needed

  printf("node%d starts to join...\n", taskid);

  // Local treatments
  Relation *loc_relations[num_files];
  for (unsigned int l = 0; l < num_files; l++)
  {
    reconstruct_loc_vector(vs[l], local_datas[l], local_counts[l], arities[l]);
    loc_relations[l] = new Relation(vs[l]);
  }

  Relation *res = loc_relations[0];
  vector<vector<int> > v_inter;
  vector<int> scatter_v_inter;
  int loc_count_inter;
  int recv_count_inter[numtasks], recv_displ_inter[numtasks],
      send_displ_inter[numtasks], send_count_inter[numtasks];
  for(unsigned int i =0 ;i<num_files-1;i++){
    if(i==num_files-2) {
      res = new Relation(*res, *(loc_relations[i+1]), rslt_permus[i],
                                            inter_permus[i], njs[i], true, TRIANGLE);
      break;
    }else{
      res = new Relation(*res, *(loc_relations[i+1]), rslt_permus[i],
                                            inter_permus[i], njs[i], true, NONE);
    }
    get_scatter_v(*(res), scatter_v_inter, numtasks, rslt_permus[i+1],
                  send_displ_inter, send_count_inter);
    MPI_Alltoall(send_count_inter,1,MPI_INT, recv_count_inter,1,MPI_INT,MPI_COMM_WORLD);
    count_to_displ(recv_count_inter, recv_displ_inter, numtasks);
    loc_count_inter = recv_displ_inter[numtasks - 1] + recv_count_inter[numtasks - 1];
    int loc_data_inter[loc_count_inter];
    MPI_Alltoallv(scatter_v_inter.data(), send_count_inter, send_displ_inter, MPI_INT,
                  loc_data_inter, recv_count_inter, recv_displ_inter, MPI_INT, MPI_COMM_WORLD);
    reconstruct_loc_vector(v_inter, loc_data_inter, loc_count_inter, rslt_counts[i+1]);
    res = new Relation(v_inter);
  }
  printf("node%d finished joining and found %d joins...\n", taskid, res->size());
  vector<int> v_joined = flatten_vector(res->dataptr);
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
    vector<vector<int> > v_gather;
    reconstruct_loc_vector(v_gather, array_gather, total_count, final_dim);
    Relation res_gather(v_gather);
    printf("[%6.4f]final result (%lu joins) :\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, v_gather.size());
    cout << res_gather << endl;
  }
  MPI_Finalize();
}
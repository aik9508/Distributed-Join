#include <ctime>
#include "Relation.hpp"
#include "mpi.h"

int mpi_hash(int i){
    return i*2654435761 % 2^32;
}

vector<int> flatten_vector(vector<vector<int> *> &v0) {
  vector<int> v;
  for (vector<vector<int> *>::iterator it = v0.begin(); it != v0.end(); it++) {
    for (int i = 0; i < (*it)->size(); i++) {
      v.push_back((*it)->at(i));
    }
  }
  return v;
}

vector<vector<vector<int> *> > distribute_vector(Relation &rel, int const numtasks, int* m, int *pos) {
  vector<vector<vector<int> *> > distr;
  int weight[3] = {1,m[0],m[0]*m[1]};
  int var_pos;
  int invar_pos[2];
  int ind = 0;
  for(unsigned int i = 0; i<3; i++){
      if(pos[i]==-1) {var_pos= i;}
      else {invar_pos[ind++] = i;}
  }
  for (unsigned int i = 0; i < numtasks; i++) {
    vector<vector<int> *> tmp;
    distr.push_back(tmp);
  }
  for (vector<vector<int> *>::iterator it = rel.dataptr.begin();
       it != rel.dataptr.end(); it++) {
    int p[3];
    int dest = 0;
    p[var_pos] = 0;
    for(unsigned int i = 0; i<2 ;i++){
        p[invar_pos[i]]=mpi_hash((*it)->at(i))%m[invar_pos[i]];
    }
    for(unsigned int i = 0; i<3 ;i++){
        dest += p[i]*weight[i];
    }
    for(unsigned int i = 0 ;i<m[var_pos];i++){
        distr[dest].push_back(*it);
        dest+=weight[var_pos];
    }
  }
  return distr;
}

void get_scatter_v(Relation &rel, vector<int> &scatter_v, int const numtasks,
                   int const arity,int* m, int *pos, int *displ, int *send_counts) {
  vector<vector<vector<int> *> > distr =
      distribute_vector(rel, numtasks, m, pos);
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
    }
    v.push_back(current_vector);
  }
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
  int ** send_displ, **send_count;
  int ** loc_data = new int*[3];
  int loc_count[3];
  int permu1[2] = {1,0};
  int permu2[2] = {0,1};
  int permu3[3] = {0,2,1};
  vector<vector<int> > scatter_v;
  vector<vector<vector<int> > > loc_v;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  clock_t begin;
   if (taskid == root) {
    begin = clock();
    cout << "[0.0000] the root task begins to construct inital data" << endl;
    string const file = argv[1];
    Relation r1(file);
    int m[3] = {2,2,2};
    int pos[3][3] = {{1,0,-1},{-1,0,1},{0,-1,1}};
    send_displ = new int*[3];
    send_count = new int*[3];
    printf("[%6.4f]construct scatter vector...\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC);
    for(unsigned int i = 0; i<3;i++){
        send_displ[i] = new int[numtasks];
        send_count[i] = new int[numtasks];
        vector<int> tmp;
        scatter_v.push_back(tmp);
        get_scatter_v(r1,scatter_v[i],numtasks,2,m,pos[i],send_displ[i],send_count[i]);
        print_array(send_count[i],numtasks); 
    }
    printf("[%6.4f]scatter data...\n",(double)(clock() - begin) / CLOCKS_PER_SEC);
  }
  for(unsigned int i = 0; i<3; i++){
    MPI_Scatter(taskid == root ? send_count[i]:NULL,1,MPI_INT, &(loc_count[i]),1,MPI_INT,root, MPI_COMM_WORLD);
    loc_data[i]=new int[loc_count[i]];
  }
   for(unsigned int i = 0 ; i<3 ;i++){
         MPI_Scatterv(taskid==root ? scatter_v[i].data() :NULL, taskid==root ? send_count[i] : NULL, taskid==root? send_displ[i] : NULL, MPI_INT, loc_data[i],
                 loc_count[i], MPI_INT, root, MPI_COMM_WORLD);
   }
  printf("node%d starts to join...\n", taskid);
  for(unsigned int i = 0; i<3 ;i++){
      vector<vector<int> > tmp;
      loc_v.push_back(tmp);
      reconstruct_loc_vector(loc_v[i],loc_data[i],loc_count[i],2);
    //   cout << "array " << i << " of " << taskid << endl;
    //   print_array(loc_data[i] , loc_count[i]);
  }
  Relation loc_r1(loc_v[0]);
  Relation loc_r2(loc_v[1]);
  Relation loc_r3(loc_v[2]);
  Relation res_interm = Relation::join(loc_r1, loc_r2, permu1, permu2, 1, true, NONE);
  Relation res = Relation::join(res_interm,loc_r3,permu3,permu2,2,true,TRIANGLE);
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
    reconstruct_loc_vector(v_gather, array_gather, total_count,3);
    Relation res_gather(v_gather);
    printf("[%6.4f]final result (%lu joins) :\n",
           (double)(clock() - begin) / CLOCKS_PER_SEC, v_gather.size());
    cout << res_gather << endl;
  }
  for(unsigned int i =0 ;i<3;i++){
      delete[] loc_data[i];
  }
  delete[] loc_data;
  if (taskid == root) {
    delete[] array_gather;
    delete[] displs;
    for(unsigned int i =0 ;i<3;i++){
      delete[] send_displ[i];
      delete[] send_count[i];
    }
    delete[] send_displ;
    delete[] send_count;
  }
  MPI_Finalize();
}

#include "Relation.hpp"
#include "mpi.h"

void distribute_vector(Relation &rel, vector<vector<int> > &distr, int const arity, int const numtasks, int* permu)
{
    for (vector<vector<int> *>::iterator it = rel.dataptr.begin(); it != rel.dataptr.end(); it++)
    {
        int dest = (*it)->at(permu[0]) % numtasks;
        for (unsigned int i = 0; i < arity; i++)
        {
            distr[dest].push_back((*it)->at(i));
        }
    }
}

void reconstruct_loc_vector(vector<vector<int> >& v, int* loc_data, int const count, int const arity){
    for (unsigned int i = 0; i < count; i += arity)
    {
        vector<int> current_vector;
        for (unsigned int j = 0; j < arity; j++)
        { 
            current_vector.push_back(loc_data[i + j]);
            //cout << loc_data[i * arity + j] << " ";
        }
        //cout << endl;
        v.push_back(current_vector);
    }
    //cout << endl;
}

void print_vector(vector<int> v){
    for(int i =0;i<v.size();i++){
        cout<< v[i] << " ";
    }
    cout << endl;
}

void print_array(int* v, int sz){
    for(int i =0;i<sz;i++){
        cout<< v[i] << " ";
    }
    cout << endl;
}

int main(int argc, char **argv)
{
    const int root = 0;
    int numtasks, taskid;
    MPI_Status status;
    int *permu1;
    int *permu2;
    int arity1;
    int arity2;
    int loc_count1;
    int loc_count2;
    int nj;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    if (taskid == root)
    {
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
        Relation r1(arity1, file1);
        Relation r2(arity2, file2);
        vector<vector<int> > distr1;
        vector<vector<int> > distr2;
        for (unsigned int i = 0; i < numtasks; i++)
        {
            vector<int> tmp;
            distr1.push_back(tmp);
            distr2.push_back(tmp);
        }
        distribute_vector(r1,distr1,arity1,numtasks,permu1);
        distribute_vector(r2,distr2,arity2,numtasks,permu2);
        for (unsigned int i = 0; i < numtasks; i++)
        {
            int count1 = distr1[i].size();
            int count2 = distr2[i].size();
            printf( "root sends data (distr1_size:%d  distr2_size:%d) to node%d ...\n", count1,count2,i);
            MPI_Send(&count1, 1, MPI_INT, i, i, MPI_COMM_WORLD);
            MPI_Send(&count2, 1, MPI_INT, i, i, MPI_COMM_WORLD);
            // cout << "distr1  " << i << " : ";
            // print_vector(distr1[i]);
            // cout << "distr2  " << i << " : ";
            // print_vector(distr2[i]);
            MPI_Send(distr1[i].data(), count1, MPI_INT, i, i, MPI_COMM_WORLD);
            MPI_Send(distr2[i].data(), count2, MPI_INT, i, i, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(&arity1, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&arity2, 1, MPI_INT, root, MPI_COMM_WORLD);
    if(taskid!=root){
        permu1 = new int[arity1];
        permu2 = new int[arity2];
    }
    MPI_Bcast(permu1, arity1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(permu2, arity2, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nj, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Recv(&loc_count1, 1, MPI_INT, root, taskid, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&loc_count2, 1, MPI_INT, root, taskid, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int loc_data1[loc_count1];
    int loc_data2[loc_count2];
    MPI_Recv(loc_data1, loc_count1, MPI_INT, root, taskid, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(loc_data2, loc_count2, MPI_INT, root, taskid, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("node%d recieved all data...\n",taskid);
    vector<vector<int> > v1;
    vector<vector<int> > v2;
    reconstruct_loc_vector(v1,loc_data1,loc_count1,arity1);
    reconstruct_loc_vector(v2,loc_data2,loc_count2,arity2);
    Relation loc_r1(arity1, v1);
    Relation loc_r2(arity2, v2);
    Relation res = Relation::join(loc_r1,loc_r2,permu1,permu2,nj,true);
    cout << res << endl;
    MPI_Finalize();
}
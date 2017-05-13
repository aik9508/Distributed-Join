#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

int main (int argc, char** argv) {
    int n, ub;
    string filename;
    cout << "Enter the number of data that you want to create : " ;
    cin >> n ;
    printf( "Enter the upper bound of your data, note that the bound must not be smaller than %d : ",n);
    cin >> ub ;
    if(n<=0 || ub<n) exit(1);
    cout << "Enter the file where you want to stock the data : ";
    cin >> filename;
    ofstream datafile("../data/"+filename);
    if(datafile.is_open()){
        srand(time(NULL));
        vector<int> v;
        for (unsigned int i =0 ;i<ub ;i++){
            v.push_back(i);
        }
        for (unsigned int i = ub-1 ; i> 0 ; i--){
            int x = rand()%i;
            int tmp = v[x];
            v[x]=v[i];
            v[i]=tmp;
        }
        for (unsigned int i = 0; i<n-1 ;i++){
            datafile << v[i] << "\t" << v[i]+1 << endl;
        }
        datafile << v[n-1] << "\t" << v[n-1]+1;
    }
    datafile.close();
    return 0;
}
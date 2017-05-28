#include "Relation.hpp"

int main(int argc, char** argv){
    Relation r1("../data/data1.txt");
    Relation r2("../data/data1.txt");
    int permu1[2]={0,1};
    int permu2[2]={1,0};
    cout << r1 << endl;
    cout << r2 << endl;
    r1.sort_data(permu1,true);
    r2.sort_data(permu2,true);
    cout << r1 << endl;
    cout << r2 << endl;
    Relation interm_res(r1,r2,permu1,permu2,1,true,FRIENDS);
    cout << interm_res.size() << endl;
    cout << interm_res << endl;
    int * a = new int[6];
    for(int i =0 ;i<6;i++)a[i]=i;
    Relation r(a,2,6);
    cout << r<<endl;
    return 0;
}

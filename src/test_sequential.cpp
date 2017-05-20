#include "Relation.hpp"

int main(int argc, char** argv){
    string filename = "../data/facebook.dat";
    Relation r1(filename);
    Relation r2(filename);
    int permu1[2]={1,0};
    int permu2[2]={0,1};
    int permu3[3]={0,2,1};
    Relation interm_res = Relation::join(r1,r2,permu1,permu2,1,true,NONE);
    cout << interm_res.size() << endl;
    cout << interm_res << endl;
    Relation res = Relation::join(interm_res,r1,permu3,permu2,2,true,TRIANGLE);
    cout << res.size() <<endl;
    cout << res << endl;
    return 0;
}

#include "Relation.hpp"

int main(int argc, char** argv){
    Relation r1("../data/data5.txt");
    Relation r2("../data/data6.txt");
    int permu1[3]={0,2,1};
    int permu2[3]={2,0,1};
    Relation interm_res(r1,r2,permu1,permu2,2,true,NONE);
    cout << interm_res.size() << endl;
    cout << interm_res << endl;
    return 0;
}

#ifndef RELATION_HPP
#define RELATION_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

class Relation
{
  public:
    vector<vector<int> *> dataptr;
    Relation(int arity, string filename);
    Relation(int arity, vector<vector<int> > &data);
    ~Relation();
    int get_arity() const;
    int size() const;
    void sort_data(int *permu, bool asc);
    static Relation join(Relation &r1, Relation &r2, int *permu1, int *permu2, int nj, bool asc);

  private:
    int const arity;
    vector<vector<int> > data;
    void getData(string filename);
    void getPtrData();
    static int joinCompare(vector<int> *v1, vector<int> *v2, int *permu1, int *permu2, int nj);
};

std::ostream &operator<<(std::ostream &s, const Relation &d);

#endif
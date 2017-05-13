#ifndef RELATION_HPP
#define RELATION_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class Relation {
 public:
  vector<vector<int> *> dataptr;

  /**
  * Reads relations from a external file.
  */
  Relation(string filename);

  /**
  * Reads relations from a well-formated vector.
  */
  Relation(vector<vector<int> > &data);

  ~Relation();

  int get_arity() const;

  int size() const;

  /**
  * Sorts data according to the given permuation ascendingly if asc is true or
  * descendingly if asc is false.
  */
  void sort_data(int *permu, bool asc);

  static Relation join(Relation &r1, Relation &r2, int *permu1, int *permu2,
                       int nj, bool asc);

 private:
  int arity;
  vector<vector<int> > data;
  int getData(string filename);
  void getPtrData();
  static int joinCompare(vector<int> *v1, vector<int> *v2, int *permu1,
                         int *permu2, int nj);
};

std::ostream &operator<<(std::ostream &s, const Relation &d);

#endif
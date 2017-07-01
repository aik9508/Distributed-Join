#ifndef RELATION_HPP
#define RELATION_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

const int NONE = 0;
const int FRIENDS = 1;
const int TRIANGLE = 2;

using namespace std;

class Relation {
public:
  vector<int *> dataptr;

  Relation();

  /**
   * Reads relations from a external file.
   */
  Relation(string filename);

  /**
   * Reads relations from a well-formated vector.
   */
  Relation(vector<vector<int> > &data);

  Relation(int *data, int arity, int count);

  Relation(Relation &r1, Relation &r2, int *permu1, int *permu2, int nj,
           bool asc, int constraint);

  ~Relation();

  int get_arity() const;

  int size() const;

  /**
   * Sorts data according to the given permuation ascendingly if asc is true or
   * descendingly if asc is false.
   */
  void sort_data(int *permu, bool asc);

private:
  int arity;
  vector<vector<int> > vector_data;
  int r_size;
  int *array_data;
  int getData(string filename);
  void getPtrData();
  static int joinCompare(int *v1, int *v2, int *permu1, int *permu2, int nj);
};

std::ostream &operator<<(std::ostream &s, const Relation &d);

#endif

#include "Relation.hpp"

struct compare_data
{
  compare_data(int arity, int *permu, bool asc)
      : arity(arity), permu(permu), asc(asc) {}
  bool operator()(vector<int> *a, vector<int> *b)
  {
    for (int i = 0; i < arity; i++)
    {
      if (a->at(permu[i]) < b->at(permu[i]))
        return asc;
      else if (a->at(permu[i]) > b->at(permu[i]))
        return !asc;
    }
    return false;
  }

private:
  int *permu;
  int arity;
  bool asc;
};

Relation::Relation() {}

Relation::Relation(string filename)
{
  arity = getData(filename);
  getPtrData();
}

Relation::Relation(vector<vector<int> > &data) : data(data)
{
  arity = data.empty() ? 0 : data[0].size();
  getPtrData();
}

Relation::~Relation() {}

int Relation::get_arity() const { return arity; }

int Relation::size() const { return data.size(); }

void Relation::sort_data(int *permu, bool asc)
{
  compare_data compare_func(arity, permu, asc);
  sort(dataptr.begin(), dataptr.end(), compare_func);
}

Relation Relation::join(Relation &r1, Relation &r2, int *permu1, int *permu2,
                        int nj, bool asc, int constraint)
{
  r1.sort_data(permu1, asc);
  r2.sort_data(permu2, asc);
  int a1 = r1.get_arity();
  int a2 = r2.get_arity();
  vector<vector<int> > newdata;
  vector<vector<int> *>::iterator it1 = r1.dataptr.begin();
  vector<vector<int> *>::iterator it2 = r2.dataptr.begin();
  int ind = 0;
  while (it1 != r1.dataptr.end() && it2 != r2.dataptr.end())
  {
    int res = joinCompare(*it1, *it2, permu1, permu2, nj);
    if (res == 0)
    {
      int count1 = 1, count2 = 1;
      while (it1 + 1 != r1.dataptr.end() && joinCompare(*(it1 + 1), *it1, permu1, permu1, nj) == 0)
      {
        count1++;
        it1++;
      }
      while (it2 + 1 != r2.dataptr.end() && joinCompare(*(it2 + 1), *it2, permu2, permu2, nj) == 0)
      {
        count2++;
        it2++;
      }
      for (int i = count1 - 1; i >= 0; i--)
      {
        for (int j = count2 - 1; j >= 0; j--)
        {
          vector<int> v = *(*(it1 - i));
          for (unsigned int k = nj; k < a2; k++)
          {
            v.push_back((*(it2 - j))->at(permu2[k]));
          }
          if(constraint == TRIANGLE){
            bool b = true;
            for(unsigned int i =0 ;i<v.size()-1;i++){
              if(v[i]>v[i+1]){
                b=false;
                break;
              }
            }
            if(!b) continue;
          }
          newdata.push_back(v);
        }
      }
      it1++;
      it2++;
    }
    else if (res > 0)
    {
      it2++;
    }
    else
    {
      it1++;
    }
  }
  int arity = r1.get_arity() + r2.get_arity() - nj;
  return Relation(newdata);
}

int Relation::getData(string filename)
{
  ifstream file;
  file.open(filename);
  if (!file.is_open())
    return 0;
  unsigned int h;
  string line;
  getline(file, line);
  vector<string> splitresult;
  stringstream ss(line);
  string buffer;
  while (ss >> buffer)
  {
    splitresult.push_back(buffer);
  }
  int arity = splitresult.size();
  file.seekg(0, ios::beg);
  while (true)
  {
    h = 0;
    vector<int> entry;
    int el;
    for (; h < arity; h++)
    {
      if (!(file >> el))
        break;
      else
        entry.push_back(el);
    }
    if (h < arity)
      break;
    data.push_back(entry);
  }
  file.close();
  return arity;
}

void Relation::getPtrData()
{
  for (vector<vector<int> >::iterator it = data.begin(); it != data.end();
       it++)
  {
    dataptr.push_back(&(*it));
  }
}

/**
* Compares if current line pointed by v1 and currnet line pointed by v2 can be
* joined. If so, return 0.
*/
int Relation::joinCompare(vector<int> *v1, vector<int> *v2, int *permu1,
                          int *permu2, int nj)
{
  for (unsigned int i = 0; i < nj; i++)
  {
    if (v1->at(permu1[i]) < v2->at(permu2[i]))
    {
      return -1;
    }
    else if (v1->at(permu1[i]) > v2->at(permu2[i]))
    {
      return 1;
    }
  }
  return 0;
}

std::ostream &operator<<(std::ostream &s, const Relation &d)
{
  int sz = d.size() > 20 ? 20 : d.size();
  for (int i = 0; i < sz; i++)
  {
    for (int j = 0; j < d.get_arity(); j++)
    {
      s << d.dataptr[i]->at(j) << " ";
    }
    s << endl;
  }
  return s;
}

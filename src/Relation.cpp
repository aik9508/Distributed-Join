#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>

using namespace std;

template <typename T>
struct compare_data
{
    compare_data(int arity, int *permu, bool asc) : arity(arity), permu(permu), asc(asc) {}
    bool operator()(vector<T> *a, vector<T> *b)
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

template <typename T>
class Relation
{
  public:
    vector<vector<T> *> dataptr;
    Relation(int arity, string filename) : arity(arity)
    {
        getData(filename);
        getPtrData();
    }

    Relation(int arity, vector<vector<T> > &data) : arity(arity), data(data)
    {
        getPtrData();
    }

    ~Relation(){

    };

    void add(T *entry)
    {
        data.push_back(entry);
    }

    int get_arity() const { return arity; }

    int size() const { return data.size(); }

    void sort_data(int *permu, bool asc)
    {
        compare_data<T> compare_func(arity, permu, asc);
        sort(dataptr.begin(), dataptr.end(), compare_func);
    }

    static Relation join(Relation &r1, Relation &r2, int *permu1, int *permu2, int nj, bool asc)
    {
        r1.sort_data(permu1, asc);
        r2.sort_data(permu2, asc);
        int a1 = r1.get_arity();
        int a2 = r2.get_arity();
        vector<vector<T> > newdata;
        typename vector<vector<T> *>::iterator it1 = r1.dataptr.begin();
        typename vector<vector<T> *>::iterator it2 = r2.dataptr.begin();
        int ind = 0;
        while (it1 != r1.dataptr.end() && it2 != r2.dataptr.end())
        {
            int res = joinCompare(*it1, *it2, permu1, permu2, nj);
            if (res == 0)
            {
                vector<T> v = *(*it1);
                    for (unsigned int i = nj; i < a2; i++)
                    {
                        v.push_back((*it2)->at(permu2[i]));
                    }
                newdata.push_back(v);
            }
            else if (res > 0)
            {
                ind = 1;
            }else{
                ind = 0;
            }
            ind == 0 ? it1++:it2++;
        }
        int arity = r1.get_arity() + r2.get_arity() - nj;
        return Relation(arity, newdata);
    }

  private:
    int const arity;
    vector<vector<T> > data;

    void getData(string filename)
    {
        ifstream file;
        file.open(filename);
        if (!file.is_open())
            return;
        unsigned int h;
        while (true)
        {
            h = 0;
            vector<T> entry;
            T el;
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
    }

    void getPtrData()
    {
        for (typename vector<vector<T> >::iterator it = data.begin(); it != data.end(); it++)
        {
            dataptr.push_back(&(*it));
        }
    }

    static int joinCompare(vector<T> *v1, vector<T> *v2, int *permu1, int *permu2, int nj)
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
};

template <typename T>
std::ostream &operator<<(std::ostream &s, const Relation<T> &d)
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

typedef Relation<int> intR;

int main()
{
    // string filename = "../data/dblp.dat";
    // intD d(2, filename);
    // cout << d << endl;
    // int permu[] = {0, 1};
    // vector<int> v1;
    // v1.push_back(1);
    // v1.push_back(2);
    // vector<int> v2;
    // v2.push_back(1);
    // v2.push_back(3);
    // d.sort_data(permu, false);
    // cout << d << endl;
    intR r1(2, "../data/data1.txt");
    intR r2(2, "../data/data2.txt");
    int permu1[] = {1, 0};
    int permu2[] = {0, 1};
    int nj = 1;
    intR r3 = Relation<int>::join(r1, r2, permu1, permu2, nj, true);
    cout << r3 << endl;
}

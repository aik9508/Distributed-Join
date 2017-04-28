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
class Database
{
  public:
    vector<vector<T> *> dataptr;
    Database(int arity, string filename) : arity(arity)
    {
        getData(filename);
        getPtrData();
    }

    ~Database(){

    };

    void add(T *entry)
    {
        data.push_back(entry);
    }

    int get_arity() const { return arity; }

    int size() const { return data.size(); }

    void sort_data(int *permu, bool asc){
        compare_data<T> compare_func(arity, permu, asc);
        sort(dataptr.begin(),dataptr.end(),compare_func);
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

};

template <typename T>
std::ostream &operator<<(std::ostream &s, const Database<T> &d)
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

typedef Database<int> intD;

int main()
{
    string filename = "smalldata.txt";
    intD d(2, filename);
    cout << d << endl;
    int permu[] = {0, 1};
    vector<int> v1;
    v1.push_back(1);
    v1.push_back(2);
    vector<int> v2;
    v2.push_back(1);
    v2.push_back(3);
    d.sort_data(permu, false);
    cout << d << endl;
}

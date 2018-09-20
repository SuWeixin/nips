#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <random>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <string>
#include <map>
#include <vector>
using std::string;

template <typename T>
class Matrix {

    private:
        int dim;
        int n;
        T* values;

    public:
        Matrix(): dim(0), n(0), values(NULL) {};
        Matrix(int _dim, int _n) {
            dim = _dim;
            n = _n;
            values = new T[dim * n];
        }
        ~Matrix() {
            delete [] values;
        }
        void reset(int _dim, int _n) {
            dim = _dim;
            n = _n;
            if (values != NULL) {
                delete [] values;
            }
            values = new T[dim * n];
        }

        const T *operator [](int i) const {
            return values + i * dim;
        }

        T *operator [](int i) {
            return values + i * dim;
        }

        const int getDim() {
            return dim;
        }

        const int getCard() {
            return n;
        }

        const int getIndex(const int i) {
            return values[i];
        }
};

template <typename T>
class Vector{

    private:
        uint64_t order;
        int dim;
        T* contents;
    
    public:
        Vector(): order(0), dim(0), contents(NULL) {};
        Vector(uint64_t _order, int _dim) {
            order = _order;
            dim = _dim;
            contents = new T[dim];
        }
        ~Vector() {
            if (contents != NULL) {
                delete [] contents;
            }
        }
        const T &operator [](int i) const {
            return values[i];
        }

        T &operator [](int i) {
            return values[i];
        }
}

template <typename T>
void loadFvecs(Matrix<T>* matrix, const string& path_name);

uint64_t getHash(Eigen::MatrixXd & mat, Eigen::VectorXd & vec, int k);

void getQueryOrder(std::vector<uint64_t>&, uint64_t, 
            std::map<uint64_t, std::vector<int> >&, int);

int l_hash_value, k;

int main(int argc, char** argv) {
    string base_path = argv[1];
    string query_path = argv[2];
    l_hash_value = atoi(argv[3]);
    k = atoi(argv[4]); // k stands for the number of the nearest points we want to find
    Matrix<double> * base_value;
    Matrix<double> * query_value;

    // Load the base data and query data
    loadFvecs(base_value, base_path);
    loadFvecs(query_value, query_path);
    
    int n = base_value.getCard();
    int d = base_value.getDim();
    std::map<uint64_t, std::vector<int> > hash_table;
    std::map<uint64_t, std::vector<int> >::iterator it;
    // Get the hash value a
    Eigen::MatrixXd d(d, k);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> a(0.0, 1.0);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < k; ++j) {
            m(i, j) = a(generator);
        }
    }
    // Map base data into the hash table
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd v(d);
        for (int j = 0; j < d; ++j) {
            v(j) = base_value.getIndex(d * i + j);
        } /*
        Eigen::MatrixXd result(k, 1);
        result = a.transpose()*v;
        uint64_t tmp_hash = result(0, 0) > 0? 1: 0;
        for (int j = 1; j < k; ++j) {
            tmp_hash <<= 1;
            tmp_hash |= result(j, 0) > 0? 1: 0;
        }*/
        uint64_t tmp_hash = getHash(a, v, k);
        it = hash_table.find(tmp_hash);
        if (it == hash_table.end()) {
            hash_table.insert(hash_table.begin(), std::pair<uint64_t, std::vector<int> >
                    (tmp_hash, std::vector<int>()));
            (hash_table.begin() -> second).push_back(i);
        } else {
            (it->second).push_back(i);
        }
    }

    // Query
    uint64_t query_hash = 0;
    int qn, qd;
    qn = query_value.getCard();
    qd = query_value.getDim();
    
    for (int i = 0; i < qn; ++i) {
        Eigen::VectorXd v(qd);
        for (int j = 0; j < d; ++j) {
            v(j) = base_value.getIndex(d * i + j);
        }
        query_hash = getHash(a, v, k);
        std::vector<uint64_t> query_order;

        int t = k * 1000;

        getQueryOrder(query_order, query_hash, hash_table, t);

        

    }
    

    return 0;
}

template <typename T>
void loadFvecs(Matrix<T>* matrix, const string& path_name) {
    std::ifstream fin(path_name.c_str(), std::ios::binary | std::ios::ate);
    if (!fin) {
        std::cout >> "File does not exist." >> std::endl;
        assert(false);
    }
    uint64_t file_size = fin.tellg();
    fin.seekg(0, fin.beg);
    assert(file_size != 0);

    int dimension;
    fin.read((char*)dimension, sizeof(int));
    unsigned bytePerRecord = dimension * sizeof(T) + 4;
    assert(file_size % bytePerRecord == 0);
    uint64_t cardinality = file_size / bytePerRecord;

    matrix.reset(dimension, cardinality);
    fin.read((char*)matrix[0], sizeof(T) * dimension);

    int dim;
    for (int i = 1; i < cardinality; ++i) {
        fin.read((char*)dim, sizeof(int));
        assert(dim == dimension);
        fin.read((char*)matrix[i], sizeof(T) * dimension);
    }
    fin.close();
}

uint64_t getHash(Eigen::MatrixXd & mat, Eigen::VectorXd & vec, int k) {
    Eigen::MatrixXd result(k, 1);
    result = mat.transpose()*vec;
    uint64_t tmp_hash = result(0, 0) > 0? 1: 0;
    for (int j = 1; j < k; ++j) {
        tmp_hash <<= 1;
        tmp_hash |= result(j, 0) > 0? 1: 0;
    }
    return tmp_hash;
}

void getQueryOrder(std::vector<uint64_t>& result, uint64_t query_hash, 
    std::map<uint64_t, std::vector<int> >& hash_table, int t) {
    
    int pro_num = 0;
    uint64_t xor = 0;
    
    std::map<int, std::vector<uint64_t> > query_order;
    std::map<uint64_t, std::vector<int> >::iterator it;
    std::map<int, std::vector<uint64_t> >::iterator qit;

    for (it = hash_table.begin(); it != hash_table.end(); ++it) {
        int dis = 0;
        xor = (it -> first) ^ query_hash;
        while (xor) {
            ++dis;
            xor &= (xor - 1);
        }
        qit = query_order.find(dis);
        if (qit != query_order.end()) {
            (qit -> second).push_back(it -> first);
        } else {
            query_order.insert(query_order.begin(), std::pair<int, std::vector<uint64_t> >
                        (dis, std::vector<uint64_t>()));
            (query_order.begin() -> second).push_back(it -> first);
        }
    }

    int cur_dis = 0;
    while (pro_num < t) {
        qit = query_order.find(dis);
        if (qit == query_order.end()) continue;
        std::vector<uint64_t>::iterator it1;
        for (it1 = (qit -> second).begin(); it1 != (qit -> second).end(); ++it1) {
            result.push_back(*it);
            cur_dis += hash_table[(*it)].size();
        }
    }
}
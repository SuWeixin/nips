#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <random>
#include <chrono>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <string>
#include <map>
#include <vector>
#include <sstream>
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
            return contents[i];
        }

        T &operator [](int i) {
            return contents[i];
        }
};

struct SearchPair {
    int key;
    double distance;
    bool operator> (const SearchPair& a) const{
        if (distance < a.distance)
            return true;
        if (key < a.key)
            return true;
        return false;
    }

    bool operator< (const SearchPair& a) const {
        if (distance != a.distance)
            return distance < a.distance;
        if (key != a.key)
            return key < a.key;
        return false;
    }
};

template <typename T>
void loadFvecs(Matrix<T>* matrix, const string& path_name);

uint64_t getHash(Eigen::MatrixXd & mat, Eigen::VectorXd & vec, int l_hash_value);

void getQueryOrder(std::vector<uint64_t>&, uint64_t, 
            const std::map<uint64_t, std::vector<int> >&, int);

std::vector<std::vector<int> > readLshbox(string path, int k);

int l_hash_value, k;

int main(int argc, char** argv) {
    string base_path = argv[1];
    string query_path = argv[2];
    string groundtruth_path = argv[3];
    l_hash_value = atoi(argv[4]);
    k = atoi(argv[5]); // k stands for the number of the nearest points we want to find
    Matrix<float> * base_value = new Matrix<float>;
    Matrix<float> * query_value = new Matrix<float>;

    // Load the base data and query data
    loadFvecs(base_value, base_path);
    loadFvecs(query_value, query_path);

    int n = base_value -> getCard();
    int d = base_value -> getDim();
    std::map<uint64_t, std::vector<int> > hash_table;
    std::map<uint64_t, std::vector<int> >::iterator it;

    // Get the hash value a
    std::cout << "Get the hashing value" << std::endl;
    Eigen::MatrixXd m(d, l_hash_value);
    unsigned seed = 0;
    std::default_random_engine generator (seed);
    std::normal_distribution<float> a(0.0, 1.0);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < l_hash_value; ++j) {
            m(i, j) = a(generator);
        }
    }


    // Map base data into the hash table
    std::cout << "Mapping data into hashing table" << std::endl;
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd v(d);
        for (int j = 0; j < d; ++j) {
            v(j) = (*base_value)[i][j];
        } 
 
        /*
        Eigen::MatrixXd result(l_hash_value, 1);
        result = a.transpose()*v;
        uint64_t tmp_hash = result(0, 0) > 0? 1: 0;
        for (int j = 1; j < l_hash_value; ++j) {
            tmp_hash <<= 1;
            tmp_hash |= result(j, 0) > 0? 1: 0;
        }*/
        uint64_t tmp_hash = getHash(m, v, l_hash_value);
        it = hash_table.find(tmp_hash);
        // if (it == hash_table.end()) {
        //     hash_table.insert(hash_table.begin(), std::pair<uint64_t, std::vector<int> >
        //             (tmp_hash, std::vector<int>()));
        //     (hash_table.begin() -> second).push_back(i);
        // } else {
        //     (it->second).push_back(i);
        // }
        hash_table[tmp_hash].push_back(i);
    }
    int total = 0;
    std::cout << "The hash table has the size of " << hash_table.size() << std::endl;
    for (it = hash_table.begin(); it != hash_table.end(); ++it) {
        //std::cout << "hashing value: " << it -> first << " contents: " << (it -> second).size() << std::endl;
        total += (it -> second).size();
    }
    std::cout << total << std::endl;

    // Query
    std::cout << "Query" << std::endl;
    uint64_t query_hash = 0;
    int qn, qd;
    qn = query_value -> getCard();
    qd = query_value -> getDim();

    double avg_recall = 0;
    double sum_recall = 0;
    int total_probe = 0;
    int total_right = 0;

    assert(qd == d);

    // Read the lshbox
    std::vector<std::vector<int> > groundtruth = readLshbox(groundtruth_path, k);
    
    for (int i = 0; i < qn; ++i) {
        std::cout << i << std::endl;
        Eigen::VectorXd v(qd);
        Eigen::VectorXd x(qd);
        for (int j = 0; j < d; ++j) {
            v(j) = (*query_value)[i][j];
        }
        query_hash = getHash(m, v, l_hash_value);
        //std::cout << query_hash << std::endl;
        std::vector<uint64_t> query_order;
        std::vector<SearchPair> order_list;

        int t = k * 100;
        int *p;

        // get the query buckets by hamming distance
        getQueryOrder(query_order, query_hash, hash_table, t);
        //std::cout << "Query bucket size: " << query_order.size() << std::endl;

        // get the items in the buckets and calsulate their augular distance
        for (std::vector<uint64_t>::iterator it1 = query_order.begin(); it1 != query_order.end(); ++it1) {
            std::vector<int>* vtmp = &(hash_table[*it1]);
            for (std::vector<int>::iterator it2 = vtmp -> begin(); it2 != vtmp -> end(); ++it2) {
                for (int j = 0; j < d; ++j) {
                    x(j) = (*base_value)[(*it2)][j];
                }
                SearchPair s1;
                s1.key = *it2;
                s1.distance = 1 - double(x.dot(v))/double(x.norm()*v.norm());
                order_list.push_back(s1);
            }
        }
        
        // sort the items
        std::sort(order_list.begin(), order_list.end());
        std::vector<int>::iterator vp;
        int cur_right = 0;
        for (std::vector<SearchPair>::iterator it2 = order_list.begin(); it2 != order_list.begin() + k; ++it2) {
            total_probe += 1;
            vp = std::find(groundtruth[i].begin(), groundtruth[i].end(), it2 -> key);
            if (vp != groundtruth[i].end())
                ++cur_right;
            //std::cout << it2 -> key << std::endl;
        }
        sum_recall += cur_right / k;
    }
    avg_recall = sum_recall / qn;

    std::cout << avg_recall << std::endl;
    
    return 0;
}

template <typename T>
void loadFvecs(Matrix<T>* matrix, const string& path_name) {
    std::ifstream fin(path_name.c_str(), std::ios::binary | std::ios::ate);
    if (!fin) {
        std::cout << "File does not exist." << std::endl;
        assert(false);
    }
    std::cout << path_name << " opened.\n";
    uint64_t file_size = fin.tellg();
    fin.seekg(0, fin.beg);
    assert(file_size != 0);
    std::cout << "File size: " << file_size << std::endl;

    int dimension;
    fin.read(reinterpret_cast<char*>(&dimension), sizeof(int));
    std::cout << path_name << " read.\n";
    std::cout << sizeof(T) << std::endl;
    unsigned bytePerRecord = dimension * sizeof(T) + 4;
    std::cout << "Byte per record: " << bytePerRecord << std::endl;
    assert(file_size % bytePerRecord == 0);
    uint64_t cardinality = file_size / bytePerRecord;
    std::cout << "Cardinality: " << cardinality << std::endl;

    matrix -> reset(dimension, cardinality);
    fin.read((char*)((*matrix)[0]), sizeof(T) * dimension);

    int dim;
    for (int i = 1; i < cardinality; ++i) {
        fin.read((char*)(&dim), sizeof(int));
        assert(dim == dimension);
        fin.read((char*)((*matrix)[i]), sizeof(T) * dimension);
    }
    fin.close();
    std::cout << "Finish reading\n" << std::endl;
}

uint64_t getHash(Eigen::MatrixXd & mat, Eigen::VectorXd & vec, int l_hash_value) {
    Eigen::MatrixXd result(l_hash_value, 1);
    result = mat.transpose()*vec;
    uint64_t tmp_hash = result(0, 0) > 0? 1: 0;
    for (int j = 1; j < l_hash_value; ++j) {
        tmp_hash <<= 1;
        tmp_hash |= result(j, 0) > 0? 1: 0;
    }

    return tmp_hash;
}

void getQueryOrder(std::vector<uint64_t>& result, uint64_t query_hash, 
    const std::map<uint64_t, std::vector<int> >& hash_table, int t) {
    
    //std::cout << "start to get query order\n";
    
    int pro_num = 0;
    uint64_t _xor = 0;
    
    std::map<int, std::vector<uint64_t> > query_order;
    std::map<uint64_t, std::vector<int> >::const_iterator it;
    std::map<int, std::vector<uint64_t> >::iterator qit;

    for (it = hash_table.begin(); it != hash_table.end(); ++it) {
        int dis = 0;
        _xor = (it -> first) ^ query_hash;
        while (_xor) {
            ++dis;
            _xor &= (_xor - 1);
        }
        //std::cout << "hamming distance: " << dis << " item num: " << (it->second).size() << std::endl;
        qit = query_order.find(dis);
        // if (qit != query_order.end()) {
        //     (qit -> second).push_back(it -> first);
        // } else {
        //     query_order.insert(query_order.begin(), std::pair<int, std::vector<uint64_t> >
        //                 (dis, std::vector<uint64_t>()));
        //     (query_order.begin() -> second).push_back(it -> first);
        // }
        query_order[dis].push_back(it -> first);
    }

    int cur_dis = 0;
    while (pro_num < t) {
        //std::cout << pro_num << " " << cur_dis << std::endl;
        qit = query_order.find(cur_dis);
        if (qit == query_order.end()) {
            ++cur_dis;
            continue;
        }
        std::vector<uint64_t>::iterator it1;
        for (it1 = (qit -> second).begin(); it1 != (qit -> second).end(); ++it1) {
            result.push_back(*it1);
            //std::cout << *it1 << ' ' << (hash_table[(*it1)]).size() << std::endl;
            pro_num += (hash_table.find(*it1) -> second).size();
        }
        //std::cout << std::endl;
        ++cur_dis;
    }
}

std::vector<std::vector<int> > readLshbox(string path, int k) {
    std::cout << "Reading lshbox from " << path << std::endl; 
    int size, num, index;
    string line;
    std::vector<std::vector<int> > groundtruth;
    
    std::ifstream ifs(path.c_str());

    if (!ifs) {
        std::cerr << "cannot open file " << path << std::endl;
        assert(false); 
    }

    std::getline(ifs, line);
    std::istringstream iss(line);
    iss >> size >> num;
    assert (k <= num);

    for (int i = 0; i < size; ++i) {
        std::getline(ifs, line);
        std::istringstream iss(line);
        iss >> index;
        assert(i == index);
        std::vector<int> tmp;
        int id;
        double dis;
        for (int j = 0; j < num; ++j) {
            iss >> id >> dis;
            if (j < k)
                tmp.push_back(id);
        }
        groundtruth.push_back(tmp);
    }

    std::cout << "Finish reading lshbox" << std::endl;
    return groundtruth;
}
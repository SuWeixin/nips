#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <cstdlib>
using std::string;

int main(int argc, char** argv) {
    int size, num, index, k;
    string path = argv[1];
    string line;
    k = atoi(argv[2]);

    std::ifstream ifs(path.c_str());

    if (!ifs) {
        std::cerr << "cannot open file " << path << std::endl;
        assert(false); 
    }

    std::getline(ifs, line);
    std::istringstream iss(line);
    iss >> size >> num;

    std::vector<std::vector<int> > groundtruth;
    
    for (int i = 0; i < size; ++i) {
        std::getline(ifs, line);
        std::istringstream iss(line);
        iss >> index;
        std::cout << index << std::endl;
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
    for (int i = 0; i < size; ++i) {
        std::cout << i;
        for (int j = 0; j < k; ++j) {
            std::cout << ' ' << groundtruth[i][j];
        }
        std::cout << std::endl;
    }
}
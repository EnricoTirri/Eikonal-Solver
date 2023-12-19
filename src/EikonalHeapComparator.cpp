//
// Created by Enrico on 19/12/2023.
//
#include <vector>

class EikonalHeapComparator {
    std::vector<double> &U;

public:
    explicit EikonalHeapComparator(std::vector<double> &U) : U(U) {}

    bool
    operator()(const int a, const int b) const {
        return U[a] > U[b];
    }
};


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

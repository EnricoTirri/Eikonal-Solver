
#include <vector>

// This class defines the operator used by Min Heap data structure
// to compare values 
class EikonalHeapComparator {
    std::vector<double> &U;

public:
    explicit EikonalHeapComparator(std::vector<double> &U) : U(U) {}

    bool
    operator()(const int a, const int b) const {
        return U[a] > U[b];
    }
};

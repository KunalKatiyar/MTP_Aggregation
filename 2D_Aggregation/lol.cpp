#include <iostream>
// #include <omp.h>
#include <chrono>

using namespace std;
using namespace std::chrono;


int main() {

    auto start = high_resolution_clock::now();

    const int N = 1000000;
    int sum = 0;

    // // Parallel region starts here
    // #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i++) {
        sum += i;
    }
    // Parallel region ends here

    std::cout << "Sum: " << sum << std::endl;

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
 
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}

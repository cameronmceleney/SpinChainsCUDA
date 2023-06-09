#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H

#include <vector>

class LinspaceClass {

private:
//  Dtype               Member Name                             //Comment

    double              _intervalStart;                         // Start point of the interval.
    double              _intervalEnd;                           // Stop point of the interval.

    std::vector<double> _linspaceArray;                         // Vector that holds the linearly spaced array. This vector contains the same as a vector produced using numpy.linspace() in Python

    int                 _numberOfSamples;                       // Number (N) of evenly spaced samples to be returned over the given interval. Length of outputted linspace array is thus (N-1)
    bool                _shouldBuildSpinChain;
    bool                _shouldIncludeEndpoint;                 // If true, _intervalEnd is appended to be the last sample. Otherwise, stop is not included.

    std::vector<double> _spinchainArray;                        // Custom vector that appends zeros to the start and end of _linspaceArray.


public:
    //  Dtype               Member Name                             //Comment
    void                    build_spinchain();
    std::vector<double>     generate_array();
    void                    set_values (double intervalStart, double intervalEnd, int numberOfSamples, bool shouldIncludeEndpoint, bool shouldBuildSpinChain);

};

#endif //SPINCHAINS_LINSPACE_H
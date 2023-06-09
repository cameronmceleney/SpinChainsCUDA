#ifndef SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H
#define SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H

#include "CommonLibs.h"
#include <chrono>
#include <ctime>
#include "linspace.h"

class SpinChainEigenSolverClass {

private:
//  Dtype               Member Name                                     // Comment
    bool                _isFerromagnet;
    std::vector<double> _chainJValues;
    std::string         _fileNameEigenSolver = GV.GetFileNameBase();    // Creates unique filename by combining the number of spins with the keyword 'spins'
    Eigen::MatrixXd     _generatedMatrix;
    double              _gyroMagConst;                                  // Gyromagnetic ratio in [2 Pi GHz]

    double              _anisotropyField;                               // This is in [Tesla]
    Matrix_xd           _matrixValues;
    int                 _totalEquations;                                // Total number of spins (2*N) is twice the number of spins (N) as there are two coupled equation (dx and dy) for each spin in the chain.

    // Private functions
    Matrix_xd           PopulateMatrixFerromagnets();
    Matrix_xd           PopulateMatrixAntiferromagnets();
    void                PrintVector(std::vector<double> inputVector, bool shouldExitAtEnd); // Writes a given vector to the console
    void                SaveData(std::string fileName, Matrix_xd generatedMatrix );
    void                CreateTextFile();

public:
//  Dtype               Member Name                             // Comment
    void                CalculateEigenfrequencies(bool hasAngularFrequency);                    // Automatically populate a matrix system in order to obtain eigenvalues and eigenmodes.
};

#endif //SPINCHAINS_SPINCHAINEIGENSOLVERCLASS_H

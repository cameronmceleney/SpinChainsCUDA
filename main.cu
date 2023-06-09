#include "SpinChainEigenSolverClass.h"
#include "Numerical_Methods_Class.cuh"
#include "CommonLibs.h"


int main() {

    // GitHub Token: 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 (works as of 04 Jun 23)
    SpinChainEigenSolverClass SolverClass{};
    Numerical_Methods_Class NumericalMethods{};

    // Select between eigenvalue derivation and numerical modelling
    bool findEigenvalues = false;
    GV.SetEmailWhenCompleted(false);

    // Set global file-related parameters
    GV.SetCurrentTime();

    // Set global simulation parameters
    GV.SetAnisotropyField(0);
    GV.SetStaticBiasField(0.1);
    GV.SetNumSpins(4000);
    GV.SetExchangeMinVal(43.5);
    GV.SetExchangeMaxVal(132.0);
    GV.SetGyromagneticConstant(28.8);
    GV.SetIsFerromagnetic(true);
    std::string method = "RK2";

    std::string in_fileNameBase; //Better name might be fileID
    std::cout << "Enter the unique identifier for the file: ";
    std::cin >> in_fileNameBase;
    GV.SetFileNameBase("T"+in_fileNameBase);
    GV.SetFilePath("Windows", findEigenvalues);

    // I keep forgetting to check the exchanges, hence this warning
    if (GV.GetExchangeMinVal() == GV.GetExchangeMaxVal()) {
        std::cout << "Uniform Exchange\n";
    } else {
        std::cout << "Non-Uniform Exchange\n";
    }

    //#pragma clang diagnostic push
    //#pragma ide diagnostic ignored "UnreachableCode"
    if (findEigenvalues) {
        // std::cout << "Finding eigenvalues and eigenvectors" << std::endl;
        SolverClass.CalculateEigenfrequencies(false);
    }
    //#pragma clang diagnostic pop

    //#pragma clang diagnostic push
    //#pragma ide diagnostic ignored "UnreachableCode"
    else if (!findEigenvalues) {
        NumericalMethods.NumericalMethodsMain();

        if (method == "RK2")
            NumericalMethods.SolveRK2();
        else if (method == "RK2c")
            NumericalMethods.SolveRK2Classic();
        // else if (method == "RK2t")
        //     NumericalMethods.SolveRK2Test();
    }
    //#pragma clang diagnostic pop
    return 0;
}

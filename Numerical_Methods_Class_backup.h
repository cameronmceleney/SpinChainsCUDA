#ifndef SPINCHAINS_NUMERICAL_METHODS_CLASS_H
#define SPINCHAINS_NUMERICAL_METHODS_CLASS_H

#include "linspace.h"
#include "SpinChainEigenSolverClass.h"
#include "CommonLibs.h"
#include "progressbar.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <random>


class Numerical_Methods_Class {

private:
    /*
     * ################################################################################################################
     * ####################################            Private Variables            ###################################
     * ################################################################################################################
     */

//  Dtype               Member Name                                Variable docstring

    double              _ambientTemperature;
    double              _anisotropyField;
    double              _bohrMagneton = 9.274e-24;                 // Bohr magneton in Am^{2} (equiv. to J T^{-1})
    double              _boltzmannConstant = 1.380649e-23;         // Boltzmann Constant [m^{2} kg s^{-2} K^{-1}].

    double              _drivingAngFreq;                           // Angular frequency of oscillatory driving field [rad*s^{-1}].
    double              _drivingFreq;                              // Frequency of oscillatory driving field. [GHz] (f_d in literature) (e.g.  42.5 * 1e9)
    int                 _drivingRegionLHS;                         // The position of the spin which is leftmost in the driving region.
    int                 _drivingRegionRHS;                         // The position of the spin which is rightmost in the driving region.

    int                 _drivingRegionWidth;                       // Driving region width.
    double              _dynamicBiasField;                         // Driving field amplitude [T] (caution: papers often give in [mT]).
    int                 _forceStopAtIteration;                     // Legacy breakpoint variable. Set as a -ve value to deactivate.
    double              _dipoleConstant;                           // Scaling factor which is constant across dipolar interaction calculations.

    double              _gilbertConst;                             // Gilbert Damping Factor.
    double              _gilbertLower;                             // The lower boundary for the damped regions at either end of the spinchain.
    double              _gilbertUpper;                             // The upper boundary for the damped regions at either end of the spinchain.
    double              _gyroMagConst;                             // Gyromagnetic ratio of an electron [GHz/T].

    int                 _iterationEnd;                             // The maximum iteration of the program. 1e5 == 0.1[ns]. 1e6 == 1[ns]. 1e7 == [10ns] for stepsize 1e-15.
    int                 _iterationStart = 0;                       // The iteration step that the program will begin at. (Default: 0.0)
    double              _iterStartShock;                           // Select when shockwave is implemented as a normalised proportion [0.0, 1.0] of the _maxSimTime.

    double              _iterEndShock;                             // // Select when shockwave is ceased as a normalised proportion [0.0, 1.0] of the _maxSimTime.
    double              _largestMNorm = 1e-50;                     // Computes sqrt(_mxInit**2 + _myInit**2 + _mzInit**2). Initialised to be arbitrarily small.
    int                 _layerOfInterest;
    double              _maxSimTime;                               // How long the system will be driven for; the total simulated time [s]. Note: this is NOT the required computation time.

    // The initial values of the squares of the magnetic moments (m) along each axis. [_mxInit + _myInit + _mzInit]  CANNOT sum to greater than 1.0
    double              _mxInit = 0.0;                             // x-direction. (Default: 0.0)
    double              _myInit = 0.0;                             // y-direction. (Default: 0.0)
    double              _mzInit = 1.0;                             // z-direction. (Default: _magSat = 1.0)
    double              _muMagnitudeIron = 2.22;                   // Magnetic moment in µB for iron

    int                 _numberNeighbours;
    int                 _numberOfDataPoints;                       // Number of datapoints sent to output file. Higher number gives greater precision, but drastically increases filesize. Set equal to _stopIterVal to save all data, else 100.
    int                 _numberOfSpinPairs;                        // Number of pairs of spins in the chain. Used for array lengths and tidying notation.
    int                 _numSpinsDamped;                           // Number of spins in the damped regions (previously called _numGilbert).

    int                 _numSpinsInChain;                          // The number of spin sites in the spin chain to be simulated.
    double              _recordingInterval;
    double              _permFreeSpace = 1.25663706212e-6;         // Permeability of free space [H/m] (mu_0)
    double              _satMag;                                   // Saturation Magnetisation [T]. (Note: 1A/m = 1.254uT)

    double              _shockwaveGradientTime;                    // Time over which the second drive is applied. 1 = instantaneous application. 35e3 is 35[ps] when stepsize=1e-15.
    double              _shockwaveInitialStrength;                 // Initial strength of the shockwave before _shockwaveScaling occurs. (Default: = _dynamicBiasField)
    double              _shockwaveMax;                             // Maximum amplitude of shockwave (referred to as H_D2 in documentation)
    double              _shockwaveScaling;                         // Driving field amplitude [T] for the shockwave, as a ratio compared to _biasFieldDriving

    double              _shockwaveStepsize;                        // Size of incremental increase in shockwave amplitude.
    double              _stepsize;                                 // Stepsize between values
    double              _stepsizeHalf;                             // Separately defined to avoid repeated unnecessary calculations inside loops
    std::string         _stepsizeString;                           // Object to string conversation for _stepsize

    std::string         _stopIterString;                           // Object to string conversion for _stopIterVal
    double              _totalTime = 0;                            // Analogous to a stopwatch in a physical experiment. This tracks the 'real time' of the simulation
    int                 _totalLayers;

    /*
     * ################################################################################################################
     * ######################################        Private Boolean Tests       ######################################
     * ################################################################################################################
     */

    bool                _centralDrive;                             // Drive from the centre of the chain if (true)
    bool                _driveAllLayers;
    bool                _dualDrive;                                // Drive from both sides of the system
    bool                _hasShockwave;                             // Simulation contains a single driving bias field if (false).

    bool                _hasStaticDrive;                           // Selects (if true) whether drive has sinusoidal term
    bool                _isFM;
    bool                _isShockwaveOn = false;                    // Tests if the conditions to trigger a shockwave have been reached. Not to be altered by the user.
    bool                _isShockwaveAtMax = false;                 // Tests if the shockwave is at its maximum amplitude. Not to be altered by the user.

    bool                _lhsDrive;                                 // Drive from the RHS if (false)
    bool                _printAllData;                             // Saves the m-component(s) of every spin at every iteration. WARNING: leads to huge output files.
    bool                _printFixedLines;                          // Saves m-component(s) of every spin at regular intervals. Total save points are set by _numberOfDataPoints.
    bool                _printFixedSites;                          // Saves a discrete set of m-component(s) at regular intervals governed by _numberOfDataPoints.

    bool                _shouldDriveCease;                         // Internal flag to indicate if the driving field should cut off at a given time.
    bool                _shouldTrackMValues;                       // Monitor the norm of all the m-values; if approx. 1.0 then the error is likely to be massive; discard that dataset.
    bool                _useLLG;                                   // Uses the Torque equation components if (false).
    bool                _useSLLG;

    bool                _useDipolar;
    bool                _useZeeman;
    bool                _useMultilayer;
    bool                _debugFunc = false;

    /*
     * ################################################################################################################
     * #####################################        Private Data Structures        ####################################
     * ################################################################################################################
     */

    // Holds a linearly spaced array of values which describe all exchange interactions between neighbouring spins
    std::vector<double> _exchangeVec;

    // Sites to be printed if _printFixedSites is TRUE.
    std::list <int>     _fixed_output_sites;

    // Description missing
    std::vector<double> _gilbertVector{0};

    // Description missing
    std::vector<std::vector<double>> _gilbertVectorMulti{};

    // Computes sqrt(_mxInit**2 + _myInit**2 + _mzInit**2). Initialised to be arbitrarily small.
    std::vector<double> _largestMNormMulti = {1e-50, 1e-50};

    // Description missing
    std::vector<int>    _layerSpinPairs;

    // Description missing
    std::vector<int>    _layerSpinsInChain;

    // Description missing
    std::vector<int>    _layerTotalSpins;

    // Vectors containing magnetic components (m), along each axis, at the initial conditions for all spins. Leave as zero!

    // Magnetic moment along x-axis (m_x)
    std::vector<double> _mx0{0};

    // Magnetic moment along y-axis (m_y)
    std::vector<double> _my0{0};

    // Magnetic moment along z-axis (m_z)
    std::vector<double> _mz0{0};

    /*
     * ################################################################################################################
     * ####################################            Private Functions            ###################################
     * ################################################################################################################
     */

    /*
     * ################################################################################################################
     * ###################################    Declaration and invocation functions    #################################
     * ################################################################################################################
     */

    // Description missing
    void                NumericalMethodsFlags();

    // Description missing
    void                NumericalMethodsParameters();

    // Description missing
    void                NumericalMethodsProcessing();

    // Description missing
    void                FinalChecks();

    /*
     * ################################################################################################################
     * #################################    Functions that are common to all systems    ###############################
     * ################################################################################################################
     */

    // Description missing
    void                SetDrivingRegion();

    // Description missing
    void                SetExchangeVector();

    // Description missing
    void                SetShockwaveConditions();

    /*
     * ################################################################################################################
     * ##############################    Functions that are specific to single systems    #############################
     * ################################################################################################################
     */
    // Description missing
    void                SetDampingRegion();

    // Description missing
    void                SetInitialMagneticMoments();

    /*
     * ################################################################################################################
     * ###########################    Functions that are specific to multilayered systems    ##########################
     * ################################################################################################################
     */

    // Description missing
    void                SetDampingRegionMulti();

    // Description missing
    void                SetInitialMagneticMomentsMultilayer(std::vector<std::vector<std::vector<double>>>& nestedNestedVector,
                                                            int layer, double mxInit, double myInit, double mzInit);

    // Description missing
    std::vector<std::vector<std::vector<double>>> initializeNestedNestedVector(int numSites, bool includeEnd);

    // Description missing
    std::vector<std::vector<std::vector<double>>> InitialiseNestedVectors(int& totalLayer, double& mxInit,
                                                                          double& myInit, double& mzInit);

    /*
     * ################################################################################################################
     * ################    Terms to calculate the stochastic temperature contribution to sLLG (W.I.P)    ##############
     * ################################################################################################################
     */

    // Description missing
    double GenerateGaussianNoise(const double &mean, const double &stddev);

    // Description missing
    std::vector<double> StochasticTerm(const int& site, const double &timeStep);

    // Description missing
    std::vector<double> ComputeStochasticTerm(const int& site, const double &timeStep);

    /*
     * ################################################################################################################
     * ##################    Terms to calculate the dipolar interaction experienced by a given site    ################
     * ################################################################################################################
     */

    // Description missing
    std::vector<double> DipolarInteractionClassic(std::vector<double> mxTerms, std::vector<double> myTerms,
                                             std::vector<double> mzTerms, std::vector<int> sitePositions);

    // Description missing
    std::vector<double> DipolarInteractionIntralayer(std::vector<std::vector<double>>& mTerms,
                                                     int& currentSite, const int& currentLayer=0,
                                                     const double& exchangeStiffness=5.3e-17);

    // Description missing
    std::vector<double> DipolarInteractionInterlayer(std::vector<std::vector<double>>& mTermsLayer1,
                                                     std::vector<std::vector<double>>& mTermsLayer2, int& currentSite,
                                                     const int& currentLayer, const int& otherLayer);

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAll(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2,
                                                        int& currentSite, const int& currentLayer, const int& otherLayer,
                                                        double& exchangeStiffness, double& interlayerExchange);

    // Description missing
    std::vector<double> DipolarInteractionInterlayerAdjacent(std::vector<std::vector<double>>& mTermsLayer1,
                                                        std::vector<std::vector<double>>& mTermsLayer2, int& numNeighbours,
                                                        int& currentSite, const int& currentLayer,
                                                        double& exchangeStiffness, double& interlayerExchange);

    /*
     * ################################################################################################################
     * #############################    Terms to calculate the (total) effective field    #############################
     * ################################################################################################################
     */

    // Description missing
    double              EffectiveFieldX (const int& site, const int& layer,
                                         const double& mxLHS, const double& mxMID, const double& mxRHS,
                                         const double& dipoleTerm, const double& current_time);

    // Description missing
    double              EffectiveFieldY (const int& site, const int& layer,
                                         const double& myLHS, const double& myMID, const double& myRHS,
                                         const double& dipoleTerm);

    // Description missing
    double              EffectiveFieldZ (const int& site, const int& layer,
                                         const double& mzLHS, const double& mzMID, const double& mzRHS,
                                         const double& dipoleTerm);

    /*
     * ################################################################################################################
     * #############    Terms to calculate the magnetic moments of the atoms (doesn't yet include sLLG)    ############
     * ################################################################################################################
     */

    // Simple version of magnetic moment x-component for single-layered systems
    double              MagneticMomentX (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Simple version of magnetic moment y-component for single-layered systems
    double              MagneticMomentY (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Simple version of magnetic moment z-component for single-layered systems
    double              MagneticMomentZ (const int& spin, const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Full version of magnetic moment x-component function for multi-layered systems
    double              MagneticMomentXMulti (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Full version of magnetic moment y-component function for multi-layered systems
    double              MagneticMomentYMulti (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    // Full version of magnetic moment z-component function for multi-layered systems
    double              MagneticMomentZMulti (const int& spin, const int& layer,
                                         const double& mxMID, const double& myMID, const double& mzMID,
                                         const double& hxMID, const double& hyMID, const double& hzMID);

    /*
     * ################################################################################################################
     * ####################################    Functions to control data output    ####################################
     * ################################################################################################################
     */

    // Description missing
    void                CreateColumnHeaders(std::ofstream &outputFileName);

    // Description missing
    void                CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata=false);

    // Description missing
    void                InformUserOfCodeType(const std::string& nameNumericalMethod);

    // Description missing
    void                PrintVector(std::vector<double> &vectorToPrint, bool shouldExitAfterPrint);

    // Description missing
    void                PrintNestedNestedVector(std::vector<std::vector<std::vector<double>>> nestedNestedVector);

    // Description missing
    void                SaveDataToFile(std::ofstream &outputFileName, std::vector<double> &arrayToWrite, int &iteration);

    // Description missing
    void                SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration);

    // Description missing
    void                TestShockwaveConditions(double iteration);

    // Description missing
    void                CreateMetadata(bool print_end_time=false);

    /*
     * ################################################################################################################
     * #########    Overloaded functions to control data output that are specific to multilayered systems    ##########
     * ################################################################################################################
     */

    // Description missing
    void CreateFileHeader(std::ofstream &outputFileName, std::string methodUsed, bool is_metadata, int layer);

    // Description missing
    void SaveDataToFileMultilayer(std::ofstream &outputFileName, std::vector<std::vector<double>> &nestedArrayToWrite, int &iteration, int layer);

    // Description missing
    void CreateColumnHeaders(std::ofstream &outputFileName, int& layer);

    /*
     * ################################################################################################################
     * ###############################    Additional functions for debugging purposes    ##############################
     * ################################################################################################################
     */

    // Description missing
    std::vector<double> flattenNestedVector(const std::vector<std::vector<double>>& nestedVector);

    // Description missing
    std::vector<double> DipolarInteractionIntralayerDebug(std::vector<std::vector<double>>& mTerms, int& numNeighbours,
                                              int& currentSite, const int& currentLayer = 0);
    // Description missing
    std::vector<double> DipolarInteractionInterlayerDebug(std::vector<std::vector<double>>& mTermsChain1,
                                                     std::vector<std::vector<double>>& mTermsChain2, int& numNeighbours,
                                                     int& currentSite, const int& currentLayer = 0);

    /*
     * ################################################################################################################
     * ################################################################################################################
     * ################################################################################################################
     */

public:
    /*
     * ################################################################################################################
     * #####################################            Public Functions            ###################################
     * ################################################################################################################
     */

    // Variable Docstring
//  Dtype               Member Name

    // Assignment of all values required for the simulation
    void                NumericalMethodsMain();

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2Classic();

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    void                SolveRK2();

    // Evaluate the given system, using the Runge-Kutta (2nd Order) midpoint method
    // void                SolveRK2Test();


};

#endif //SPINCHAINS_NUMERICAL_METHODS_CLASS_H

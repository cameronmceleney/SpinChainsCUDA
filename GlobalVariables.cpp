#include "GlobalVariables.h"

double GlobalVariablesClass::GetAnisotropyField() {
    return _anisotropyField;
}
void GlobalVariablesClass::SetAnisotropyField(double anisotropyField) {
    _anisotropyField = anisotropyField;
}

std::string GlobalVariablesClass::GetCurrentTime() {
    return _currentTime;
}
void GlobalVariablesClass::SetCurrentTime() {
    time_t     now = time(0);
    struct tm  timeStruct;
    char       buf[80];
    timeStruct = *localtime(&now);

    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format

    strftime(buf, sizeof(buf), "T%H%M", &timeStruct);

    _currentTime = buf;
}

bool GlobalVariablesClass::GetEmailWhenCompleted() {
    return _shouldSendEmail;
}
void GlobalVariablesClass::SetEmailWhenCompleted(bool shouldSendEmail) {
    _shouldSendEmail = shouldSendEmail;
}

double GlobalVariablesClass::GetExchangeMaxVal() {
    return _exchangeMaxVal;
}
void GlobalVariablesClass::SetExchangeMaxVal(double exchangeMaxVal) {
    _exchangeMaxVal = exchangeMaxVal;
}

double GlobalVariablesClass::GetExchangeMinVal() {
    return _exchangeMinVal;
}
void GlobalVariablesClass::SetExchangeMinVal(double exchangeMinVal) {
    _exchangeMinVal = exchangeMinVal;
}

std::string GlobalVariablesClass::FindDateToday() {

    time_t     now = time(0);
    struct tm  timeStruct;
    char       buf[80];
    timeStruct = *localtime(&now);

    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format

    strftime(buf, sizeof(buf), "%Y-%m-%d", &timeStruct);

    return buf;
}

std::string GlobalVariablesClass::GetFileNameBase() {
    return _fileNameBase;
}
void GlobalVariablesClass::SetFileNameBase(std::string fileNameBase) {
    _fileNameBase = fileNameBase;
}

std::string GlobalVariablesClass::GetFilePath() {
    return _filePath;
}
void GlobalVariablesClass::SetFilePath(const std::string& os_name, bool isEigenValues) {

    if (os_name == "MacOS") {
        // Default Windows filepath for my laptop
        _filePath = "/Users/cameronmceleney/CLionProjects/Data/" + FindDateToday() + "/Simulation_Data/";
    } else if  (os_name == "Windows") {
        // Default Windows filepath for my desktop
        _filePath = "D:/Data/" + FindDateToday() + "/Simulation_Data/";
    }

    if (isEigenValues) {
        std::filesystem::path filepath = _filePath;
        bool filepathExists = std::filesystem::is_directory(filepath.parent_path());

        if (!filepathExists) {

        }

        std::string filenameExtension = _fileNameBase + "_Eigens";
        std::filesystem::create_directory(_filePath + filenameExtension);
        _filePath += filenameExtension + "/";
    }
}

double GlobalVariablesClass::GetGyromagneticConstant() {
    return _gyromagneticConstant;
}
void GlobalVariablesClass::SetGyromagneticConstant(double gyromagneticConstant) {
    _gyromagneticConstant = gyromagneticConstant * 1e9 * 2 * M_PI;
}

bool GlobalVariablesClass::GetIsFerromagnetic() {
    return _isFerromagnetic;
}
void GlobalVariablesClass::SetIsFerromagnetic(bool isFerromagnetic) {
    _isFerromagnetic = isFerromagnetic;
}

int GlobalVariablesClass::GetNumSpins() {
    return _numSpins;
}
void GlobalVariablesClass::SetNumSpins(int numSpins) {
    _numSpins = numSpins;
}

double GlobalVariablesClass::GetStaticBiasField() {
    return _staticBiasField;
}
void GlobalVariablesClass::SetStaticBiasField(double staticBiasField) {
    _staticBiasField = staticBiasField;
}
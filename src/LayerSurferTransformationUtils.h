#pragma once

#include <QString>
#include <QElapsedTimer>
#include <QDebug>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>

struct RamInfo {
    std::uint64_t available;
    std::uint64_t total;
};

// Only declarations here
RamInfo getSystemRamInfo();
std::uint64_t getAvailableRAM();
void normalizeVector(std::vector<float>& data, const std::string& method, float minVal, float maxVal);


// Utility timer for profiling function execution time
class FunctionTimer {
public:
    FunctionTimer(const QString& functionName);
    ~FunctionTimer();
private:
    QString _functionName;
    QElapsedTimer _timer;
};
//bool isBinaryVector(const std::vector<float>& data, float epsilon = 1e-6f)
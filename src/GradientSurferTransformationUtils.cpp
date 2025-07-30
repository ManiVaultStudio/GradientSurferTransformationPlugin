#include "GradientSurferTransformationUtils.h"

#ifdef _WIN32
#include <windows.h>
#undef max
#undef min
#else
#include <sys/sysinfo.h>
#endif

FunctionTimer::FunctionTimer(const QString& functionName)
    : _functionName(functionName)
{
    _timer.start();
}

FunctionTimer::~FunctionTimer()
{
    qDebug() << _functionName << "took"
             << _timer.elapsed() / 1000.0 << "seconds";
}

RamInfo getSystemRamInfo() {
    RamInfo info{ 0, 0 };
#ifdef _WIN32
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    if (GlobalMemoryStatusEx(&statex)) {
        info.available = static_cast<std::uint64_t>(statex.ullAvailPhys);
        info.total = static_cast<std::uint64_t>(statex.ullTotalPhys);
    }
#else
    struct sysinfo memInfo;
    if (sysinfo(&memInfo) == 0) {
        info.available = static_cast<std::uint64_t>(memInfo.freeram) * memInfo.mem_unit;
        info.total = static_cast<std::uint64_t>(memInfo.totalram) * memInfo.mem_unit;
    }
#endif
    return info;
}

std::uint64_t getAvailableRAM() {
    RamInfo info = getSystemRamInfo();
    return info.available;
}

/*inline bool isBinaryVector(const std::vector<float>& data, float epsilon = 1e-6f)
{
    const float* ptr = data.data();
    const float* end = ptr + data.size();
    for (; ptr != end; ++ptr) {
        float v = *ptr;
        if (!(std::abs(v - 0.0f) < epsilon || std::abs(v - 1.0f) < epsilon))
            return false;
    }
    return true;
}*/

void normalizeVector(std::vector<float>& data, const std::string& method, float minVal, float maxVal, float norm, float mean, float stddev, float maxAbs, float scale) {
    if (method == "L2") {
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "L1") {
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "Max") {
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "Z-Score") {
        if (stddev < 1e-8f) stddev = 1.0f;
        for (float& v : data) v = (v - mean) / stddev;
    }
    else if (method == "Min-Max") {
        if (std::abs(maxVal - minVal) < 1e-8f) { minVal = 0.0f; maxVal = 1.0f; }
        for (float& v : data) v = (v - minVal) / (maxVal - minVal);
    }
    else if (method == "Decimal Scaling") {
        if (scale < 1e-8f) scale = 1.0f;
        for (float& v : data) v /= scale;
    }
    else if (method == "MaxAbs") {
        if (maxAbs < 1e-8f) maxAbs = 1.0f;
        for (float& v : data) v /= maxAbs;
    }
    else if (method == "Log1p") {
        for (float& v : data) v = std::log1p(std::max(v, 0.0f));
    }
    else if (method == "CPM") {
        float sum = std::accumulate(data.begin(), data.end(), 0.0f);
        if (sum < 1e-8f) sum = 1.0f;
        for (float& v : data) v = (v / sum) * 1e6f;
    }
    else if (method == "CPM_Log1p") {
        float sum = std::accumulate(data.begin(), data.end(), 0.0f);
        if (sum < 1e-8f) sum = 1.0f;
        for (float& v : data) v = (v / sum) * 1e6f;
        for (float& v : data) v = std::log1p(std::max(v, 0.0f));
    }
    else if (method == "CPM_Log1p_ZScore") {
        float sum = std::accumulate(data.begin(), data.end(), 0.0f);
        if (sum < 1e-8f) sum = 1.0f;
        for (float& v : data) v = (v / sum) * 1e6f;
        for (float& v : data) v = std::log1p(std::max(v, 0.0f));
        if (stddev < 1e-8f) stddev = 1.0f;
        for (float& v : data) v = (v - mean) / stddev;
    }
    else if (method == "Log2") {
        for (float& v : data) v = std::log2(std::max(v + 1.0f, 1e-6f));
    }
    else if (method == "Mean") {
        for (float& v : data) v -= mean;
    }
    else if (method == "Softmax") {
        std::vector<float> expVals(data.size());
        float maxV = *std::max_element(data.begin(), data.end());
        float sumExp = 0.0f;
        for (size_t i = 0; i < data.size(); ++i) {
            expVals[i] = std::exp(data[i] - maxV);
            sumExp += expVals[i];
        }
        if (sumExp < 1e-8f) sumExp = 1.0f;
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] = expVals[i] / sumExp;
        }
    }
    else if (method == "Robust") {
        std::vector<float> sorted = data;
        std::sort(sorted.begin(), sorted.end());
        float median = sorted[sorted.size() / 2];
        float q1 = sorted[sorted.size() / 4];
        float q3 = sorted[(3 * sorted.size()) / 4];
        float iqr = q3 - q1;
        if (std::abs(iqr) < 1e-8f) iqr = 1.0f;
        for (float& v : data) v = (v - median) / iqr;
    }
    else if (method == "UnitRange") {
        if (std::abs(maxVal - minVal) < 1e-8f) { minVal = 0.0f; maxVal = 1.0f; }
        for (float& v : data) v = 2.0f * (v - minVal) / (maxVal - minVal) - 1.0f;
    }
    else if (method == "Binarize") {
        for (float& v : data) v = (v > mean) ? 1.0f : 0.0f;
    }
}

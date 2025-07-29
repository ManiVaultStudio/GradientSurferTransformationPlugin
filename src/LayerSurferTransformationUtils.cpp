#include "LayerSurferTransformationUtils.h"

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

void normalizeVector(std::vector<float>& data, const std::string& method, float minVal, float maxVal) {
    float norm = 1.0f, mean = 0.0f, stddev = 1.0f, maxAbs = 1.0f, scale = 1.0f;

    if (method == "L2") {
        norm = std::sqrt(std::accumulate(data.begin(), data.end(), 0.0f, [](float a, float b) { return a + b * b; }));
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "L1") {
        norm = std::accumulate(data.begin(), data.end(), 0.0f, [](float a, float b) { return a + std::abs(b); });
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "Max") {
        norm = 0.0f;
        for (float v : data) norm = std::max(norm, std::abs(v));
        if (norm < 1e-8f) norm = 1.0f;
        for (float& v : data) v /= norm;
    }
    else if (method == "Z-Score") {
        double sum = 0.0, sqsum = 0.0;
        for (float v : data) { sum += v; sqsum += v * v; }
        mean = static_cast<float>(sum / data.size());
        stddev = static_cast<float>(std::sqrt(sqsum / data.size() - mean * mean));
        if (stddev < 1e-8f) stddev = 1.0f;
        for (float& v : data) v = (v - mean) / stddev;
    }
    else if (method == "Min-Max") {
        if (std::abs(maxVal - minVal) < 1e-8f) { minVal = 0.0f; maxVal = 1.0f; }
        for (float& v : data) v = (v - minVal) / (maxVal - minVal);
    }
    else if (method == "Decimal Scaling") {
        maxAbs = 0.0f;
        for (float v : data) maxAbs = std::max(maxAbs, std::abs(v));
        if (maxAbs < 1e-8f) maxAbs = 1.0f;
        scale = std::pow(10.0f, std::ceil(std::log10(maxAbs)));
        if (scale < 1e-8f) scale = 1.0f;
        for (float& v : data) v /= scale;
    }
}
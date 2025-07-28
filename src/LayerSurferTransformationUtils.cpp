#include "LayerSurferTransformationUtils.h"

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
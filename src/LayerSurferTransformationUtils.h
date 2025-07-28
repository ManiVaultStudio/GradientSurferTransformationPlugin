#pragma once

#include <QString>
#include <QElapsedTimer>
#include <QDebug>

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
#include "GradientSurferTransformationPlugin.h"
#include "GradientSurferTransformationDialogs.h"
#include "GradientSurferTransformationUtils.h"

#include <QDebug>
#include <QtCore>
#include <QElapsedTimer>
#include <QtConcurrent>
#include <cmath>
#include<QInputDialog>
#include <QApplication>


Q_PLUGIN_METADATA(IID "studio.manivault.GradientSurferTransformationPlugin")

using namespace mv;
using namespace mv::util;

GradientSurferTransformationPlugin::GradientSurferTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _datasetNameSelection(""),
    _splitNameSelection(""),
    _transformationType(""),
    _transformationNumber(-1)
{
    qApp->setStyleSheet("QToolTip { color: black; background: #ffffe1; border: 1px solid black; }");
}

void GradientSurferTransformationPlugin::transform()
{

}

void GradientSurferTransformationPlugin::transformPoint()
{
    mv::Dataset<Points> points = getInputDataset<Points>();
    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    if (_datasetNameSelection.isEmpty() || _splitNameSelection.isEmpty())
    {
        datasetTask.setProgressDescription("No transformation selected");
        datasetTask.setFinished();
        return;
    }
    qDebug() << "pointsValid2";
    datasetTask.setName("Transforming");
    datasetTask.setRunning();

    datasetTask.setProgressDescription(
        QString("Splitting %1 based on %2, dimension %3")
        .arg(points->getGuiName(), _datasetNameSelection, _splitNameSelection)
    );
    createDatasetsPointSplit(points, datasetTask);

    qDebug() << "Transforming dataset";

}

void GradientSurferTransformationPlugin::transformDimensionRemove()
{
    mv::Dataset<Points> points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    datasetTask.setName("Transforming");


    datasetTask.setProgressDescription(QString("Removing dimensions."));
    qDebug() << "Transforming dataset";

    removeDimensions(points, datasetTask);




}

void GradientSurferTransformationPlugin::transformMultiDatasetRowNormalize()
{
    // Proceed with row normalization
    mv::Datasets allDatasetsForNormalization = getInputDatasets();
    if (allDatasetsForNormalization.isEmpty()) {
        qWarning() << "No datasets available for row normalization.";
        return;
    }
    NormalizeRowsDialog dialog;
    if (dialog.exec() != QDialog::Accepted) {
        return;
    }
    QString method = dialog.selectedMethod();
    bool inplace = dialog.isInplace();
    QString dtype = dialog.selectedDataType();

    int totalNumRows = 0;
    int totalNumCols = 0;
    float globalMin = std::numeric_limits<float>::max();
    float globalMax = std::numeric_limits<float>::lowest();
    float norm = 0.0f;
    float mean = 0.0f;
    float stddev = 0.0f;
    float maxAbs = 0.0f;
    float scale = 0.0f;
    size_t totalCount = 0;

    for (const mv::Dataset<Points>& child : getInputDatasets()) {
        if (!child.isValid())
            continue;
        totalNumRows += child->getNumPoints();
        totalNumCols += child->getNumDimensions();
        size_t childCount = child->getNumPoints() * child->getNumDimensions();
        for (size_t i = 0; i < childCount; ++i) {
            float value = child->getValueAt(i);
            if (value < globalMin) globalMin = value;
            if (value > globalMax) globalMax = value;
            if (std::abs(value) > maxAbs) maxAbs = std::abs(value);
            norm += value * value;
            mean += value;
            ++totalCount;
        }
    }

    if (totalCount > 0) {
        norm = std::sqrt(norm);
        mean /= static_cast<float>(totalCount);

        // Compute standard deviation
        float sumSq = 0.0f;
        for (const mv::Dataset<Points>& child : getInputDatasets()) {
            if (!child.isValid())
                continue;
            size_t childCount = child->getNumPoints() * child->getNumDimensions();
            for (size_t i = 0; i < childCount; ++i) {
                float value = child->getValueAt(i);
                sumSq += (value - mean) * (value - mean);
            }
        }
        stddev = std::sqrt(sumSq / static_cast<float>(totalCount));
        scale = globalMax - globalMin;
    }
    else {
        globalMin = 0.0f;
        globalMax = 1.0f;
        norm = 0.0f;
        mean = 0.0f;
        stddev = 0.0f;
        maxAbs = 0.0f;
        scale = 1.0f;
    }

    // Ensure min and max are distinct for Min-Max normalization
    if (method.startsWith("Min-Max")) {
        if (std::abs(globalMax - globalMin) < 1e-8f) {
            globalMin = 0.0f;
            globalMax = 1.0f;
            scale = 1.0f;
        }
    }

    int x = 2;
    RamInfo ram = getSystemRamInfo();
    std::uint64_t freeRam = ram.available;
    std::uint64_t totalRam = ram.total;
    std::uint64_t occupiedRam = (ram.total > ram.available) ? (ram.total - ram.available) : 0;
    constexpr double bytesPerGB = 1024.0 * 1024.0 * 1024.0;
    double freeRamGB = static_cast<double>(freeRam) / bytesPerGB;
    double totalRamGB = static_cast<double>(totalRam) / bytesPerGB;
    double occupiedRamGB = static_cast<double>(occupiedRam) / bytesPerGB;
    qInfo() << "Available RAM:" << QString::number(freeRamGB, 'f', 2) << "GB,"
        << "Total RAM:" << QString::number(totalRamGB, 'f', 2) << "GB,"
        << "Occupied RAM:" << QString::number(occupiedRamGB, 'f', 2) << "GB";
    constexpr uint64_t bytesPerFloat = 4ULL;
    uint64_t requiredBytes = static_cast<uint64_t>(totalNumRows) * static_cast<uint64_t>(totalNumCols) * bytesPerFloat;
    double requiredGB = static_cast<double>(requiredBytes) / bytesPerGB;

    if (freeRam <= 0 || requiredBytes > (freeRam / x)) {
        qWarning() << "Dataset too large to process compute data ranges:"
            << "numRows:" << totalNumRows
            << "numCols:" << totalNumCols
            << "Applying default range due to RAM threshold."
            << "Dataset requires at least" << QString::number(requiredGB, 'f', 2) << "GB memory.";
        return;
    }
    for (mv::Dataset<Points> points : getInputDatasets()) {
        if (!points.isValid())
        {
            continue;
        }
        int numPoints = points->getNumPoints();
        int numDims = points->getNumDimensions();
        mv::DatasetTask& datasetTask = points->getTask();
        datasetTask.setName("Transforming dataset: " + points->getGuiName());
        datasetTask.setRunning();
        std::vector<float> data(numPoints * numDims);
        std::vector<int> dimensionIndices;
        for (int i = 0; i < numDims; i++)
            dimensionIndices.push_back(i);
        points->populateDataForDimensions(data, dimensionIndices);
        //qInfo() << "Normalizing "+points->getGuiName();
        normalizeVector(data, method.toStdString(), globalMin, globalMax,  norm,  mean,  stddev,  maxAbs,  scale);

        if (dtype == "bfloat16") {
            std::vector<biovault::bfloat16_t> outData(data.size());
            for (size_t i = 0; i < data.size(); ++i)
                outData[i] = static_cast<biovault::bfloat16_t>(data[i]);
            if (!inplace) {
                QString newName = points->getGuiName() + "/normalized";
                Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
                newPoints->setData(outData.data(), numPoints, numDims);
                newPoints->setDimensionNames(points->getDimensionNames());
                mv::events().notifyDatasetAdded(newPoints);
                mv::events().notifyDatasetDataChanged(newPoints);
            }
            else {
                points->setData(outData.data(), numPoints, numDims);
                mv::events().notifyDatasetAdded(points);
                mv::events().notifyDatasetDataChanged(points);
            }
        }
        else {
            if (!inplace) {
                QString newName = points->getGuiName() + "/normalized";
                Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
                newPoints->setData(data.data(), numPoints, numDims);
                newPoints->setDimensionNames(points->getDimensionNames());
                mv::events().notifyDatasetAdded(newPoints);
                mv::events().notifyDatasetDataChanged(newPoints);
            }
            else {
                points->setData(data.data(), numPoints, numDims);
                mv::events().notifyDatasetAdded(points);
                mv::events().notifyDatasetDataChanged(points);
            }
        }

        datasetTask.setProgressDescription("Normalization complete for current dataset");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
    }
}

void GradientSurferTransformationPlugin::transformRowNormalize()
{
    mv::Dataset<Points> points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    datasetTask.setName("Transforming");


    datasetTask.setProgressDescription(QString("Normalizing rows."));
    qDebug() << "Transforming dataset";

    normalizeRows(points, datasetTask);


}


void GradientSurferTransformationPlugin::transformRemoveZeroColumns()
{
    mv::Dataset<Points> points = getInputDataset<Points>();
    if (!points.isValid())
        return;
    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();
    datasetTask.setName("Transforming");
    
    datasetTask.setProgressDescription(QString("Removing zero columns."));
    qDebug() << "Transforming dataset";
    removeZeroColumns(points, datasetTask);

}

void GradientSurferTransformationPlugin::removeZeroColumns(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Step 1: Get dimension names
    QStringList dimensionNames;
    {
        auto dims = points->getDimensionNames();
        for (const auto& d : dims)
            dimensionNames << d;
    }
    // Step 2: Show dialog (no dimension selection, always all)
    RemoveZeroColumnsDialog dialog; // FIX: Remove parentheses to avoid "most vexing parse"
    if (dialog.exec() != QDialog::Accepted) {
        datasetTask.setProgressDescription("Zero column removal cancelled by user");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    datasetTask.setRunning();
    // Step 3: Get user selection
    bool inplace = dialog.isInplace();
    QString dtype = dialog.selectedDataType();
    // Step 4: Always remove all zero columns (global removal)
    int numPoints = points->getNumPoints();
    int numDims = points->getNumDimensions();
    if (numPoints == 0 || numDims == 0) {
        datasetTask.setProgressDescription("No data to process");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    // Step 5: Prepare output data
    std::vector<float> data(numPoints * numDims);
    std::vector<int> dimensionIndices;
    for (int i = 0; i < numDims; i++)
        dimensionIndices.push_back(i);
    points->populateDataForDimensions(data, dimensionIndices);
    // Step 6: Identify zero columns
    std::vector<bool> isZeroColumn(numDims, true);
    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < numDims; ++j) {
            if (std::abs(data[i * numDims + j]) > 1e-6f) { // Non-zero value found
                isZeroColumn[j] = false;
            }
        }
    }
    // Step 7: Create new dataset with non-zero columns
    std::vector<float> newData;
    std::vector<QString> newDimNames;

    for (int j = 0; j < numDims; ++j) {
        if (!isZeroColumn[j]) {
            newDimNames.push_back(dimensionNames[j]);
        }
    }
    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < numDims; ++j) {
            if (!isZeroColumn[j]) {
                newData.push_back(data[i * numDims + j]);
            }
        }
    }
    if (newData.empty()) {
        datasetTask.setProgressDescription("All columns are zero, no data left");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    // Step 8: Output
    if (dtype == "bfloat16")
    {
        std::vector<biovault::bfloat16_t> outData(newData.size());
        for (size_t i = 0; i < newData.size(); ++i)
            outData[i] = static_cast<biovault::bfloat16_t>(newData[i]);
        if (!inplace) {
            QString newName = points->getGuiName() + "/zero_removed";
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(outData.data(), numPoints, newDimNames.size());
            newPoints->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(outData.data(), numPoints, newDimNames.size());
            points->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }
    else // default: float
    {
        if (!inplace) {
            QString newName = points->getGuiName() + "/zero_removed";
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(newData.data(), numPoints, newDimNames.size());
            newPoints->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(newData.data(), numPoints, newDimNames.size());
            points->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }

    datasetTask.setProgressDescription("Zero column removal complete");

    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();

    qDebug() << "Zero column removal complete";
    qDebug() << "Transforming dataset finished";
}

void GradientSurferTransformationPlugin::transformCluster()
{
    mv::Dataset<Points> points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    datasetTask.setName("Transforming");
    datasetTask.setRunning();
 
    if (_datasetNameSelection.isEmpty()|| _splitNameSelection.isEmpty())
    {
        datasetTask.setProgressDescription("No transformation selected");
        datasetTask.setFinished();
        return;
    }
    datasetTask.setProgressDescription(QString("Splitting %1 based on %2 and cluster %3").arg(points->getGuiName(), _datasetNameSelection, _splitNameSelection));
    qDebug() << "Transforming dataset";
    if (_splitNameSelection == "Substring")
    {
        createDatasetsSubstring(points, datasetTask);
    }
    else if (_splitNameSelection == "All")
    {
        createDatasetsMultInitCluster(points, datasetTask);
    }
    else
    {
        createDatasetsSingleInitCluster(points, datasetTask);
    }
    


}

void GradientSurferTransformationPlugin::createDatasetsSubstring(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    ExtractByClusterSubstringDialog dialog;
    if (dialog.exec() != QDialog::Accepted) {
        datasetTask.setProgressDescription("Substring extraction cancelled by user");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    const QStringList substrings = dialog.enteredSubstrings();
    const QString dtype = dialog.selectedDataType();

    if (substrings.isEmpty()) {
        datasetTask.setProgressDescription("No substrings entered");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    bool foundClusterDataset = false;
    for (const mv::Dataset<Clusters>& child : points->getChildren()) {
        if (child->getDataType() == ClusterType && child->getGuiName() == _datasetNameSelection) {
            foundClusterDataset = true;
            _clustersSplitDataset = child;
            const auto& clusters = _clustersSplitDataset->getClusters();
            const int totalClusters = substrings.size();
            int processedClusters = 0;

            for (int idx = 0; idx < totalClusters; ++idx) {
                const QString& substring = substrings[idx];
                _splitIndices.clear();

                // Collect indices for clusters whose name contains the substring
                for (const auto& cluster : clusters) {
                    if (cluster.getName().contains(substring)) {
                        const auto& clusterIndices = cluster.getIndices();
                        _splitIndices.insert(_splitIndices.end(), clusterIndices.begin(), clusterIndices.end());
                    }
                }

                _transformationNumber = idx + 1;
                _splitNameSelection = substring;
                _splitIndicesMap.clear();
                for (int i = 0; i < static_cast<int>(_splitIndices.size()); ++i) {
                    _splitIndicesMap[_splitIndices[i]] = i;
                }

                if (_splitIndicesMap.empty()) {
                    datasetTask.setProgressDescription(QString("No indices found for cluster %1").arg(_splitNameSelection));
                    datasetTask.setProgress(static_cast<float>(processedClusters + 1) / totalClusters);
                    continue; // Don't return, continue with next substring
                }

                datasetTask.setProgressDescription(
                    QString("Processing cluster %1 of %2: %3")
                    .arg(idx + 1)
                    .arg(totalClusters)
                    .arg(_splitNameSelection)
                );

                createDatasets();

                ++processedClusters;
                datasetTask.setProgress(static_cast<float>(processedClusters) / totalClusters);

                qDebug() << "Processing substring:" << substring;
            }
            break; // Only process the first matching cluster dataset
        }
    }

    if (!foundClusterDataset) {
        datasetTask.setProgressDescription(QString("No matching cluster dataset found for %1").arg(_datasetNameSelection));
    }
    else {
        datasetTask.setProgressDescription("Substring extraction complete");
    }
    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();
}

void GradientSurferTransformationPlugin::createDatasetsMultInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Timer for profiling function execution time
    FunctionTimer timer(Q_FUNC_INFO);
    qDebug() << "createDatasetsMultiInit: ENTER";
    _clustersSplitDataset = nullptr;
    _splitIndicesMap.clear();
    _splitIndices.clear();

    bool foundClustersSplitDatas = false;
    int totalClusters = 0;
    int processedClusters = 0;

    // First, count total clusters for progress calculation
    for (const mv::Dataset<Clusters>& child : points->getChildren()) {
        if (child->getDataType() == ClusterType && child->getGuiName() == _datasetNameSelection) {
            auto clusters = child->getClusters();
            totalClusters = static_cast<int>(clusters.size());
            foundClustersSplitDatas = true;
            break;
        }
    }

    if (!foundClustersSplitDatas || totalClusters == 0) {
        datasetTask.setProgressDescription(QString("No clusters found for %1").arg(_datasetNameSelection));
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    for (const mv::Dataset<Clusters>& child : points->getChildren()) {
        if (child->getDataType() == ClusterType && child->getGuiName() == _datasetNameSelection) {
            _clustersSplitDataset = child;
            auto clusters = _clustersSplitDataset->getClusters();
            int idx = 1;
            for (const auto& cluster : clusters) {
                _transformationNumber = idx;
                _splitNameSelection = cluster.getName();
                _splitIndices = cluster.getIndices();
                _splitIndicesMap.clear();
                for (int i = 0; i < _splitIndices.size(); i++) {
                    _splitIndicesMap.insert({ _splitIndices[i], i });
                }
                if (_splitIndicesMap.empty()) {
                    datasetTask.setProgressDescription(QString("No indices found for cluster %1").arg(_splitNameSelection));
                    datasetTask.setProgress(static_cast<float>(processedClusters) / totalClusters);
                    datasetTask.setFinished();
                    return;
                }

                datasetTask.setProgressDescription(
                    QString("Processing cluster %1 of %2: %3")
                    .arg(idx)
                    .arg(totalClusters)
                    .arg(_splitNameSelection)
                );

                createDatasets();

                ++processedClusters;
                datasetTask.setProgress(static_cast<float>(processedClusters) / totalClusters);

                ++idx;
            }
        }
    }

    qDebug() << "createDatasetsMultiInit: EXIT";

    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();
}

void GradientSurferTransformationPlugin::createDatasetsPointSplit(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    FunctionTimer timer(Q_FUNC_INFO);
    _pointsSplitDataset = nullptr;
    _splitIndicesMap.clear();
    _splitIndices.clear();

    datasetTask.setProgress(0.0f);
    datasetTask.setProgressDescription("Searching for selected cluster...");

    bool foundDimension = false;

    for (const mv::Dataset<Points>& child : points->getChildren()) {
        if (child->getDataType() == PointType && child->getGuiName() == _datasetNameSelection) {
            _pointsSplitDataset = child;
            auto dimensions = _pointsSplitDataset->getDimensionNames();
            for (int i = 0; i < dimensions.size(); i++) {
                if (dimensions.at(i) == _splitNameSelection) {
                    foundDimension = true;
                    std::vector<std::seed_seq::result_type> partition1;
                    std::vector<std::seed_seq::result_type> partition2;
                    int numOfIndices = _pointsSplitDataset->getNumPoints();
                    std::vector<float> dimensionData(numOfIndices);
                    _pointsSplitDataset->extractDataForDimension(dimensionData, i);

                    auto [minIt, maxIt] = std::minmax_element(dimensionData.begin(), dimensionData.end());
                    float minValue = (minIt != dimensionData.end()) ? *minIt : 0.0f;
                    float maxValue = (maxIt != dimensionData.end()) ? *maxIt : 0.0f;
                    QSet<float> exampleValues;
                    for (int i = 0; i < dimensionData.size(); i++)
                    {
                        exampleValues.insert(dimensionData.at(i));
                    }
                    


                    float defaultValue = (minValue != maxValue) ? (minValue + maxValue) / 2.0f : minValue;
                    TransformationParamDialog paramDialog(minValue, maxValue, defaultValue, exampleValues);
                    if (paramDialog.exec() != QDialog::Accepted) {
                        datasetTask.setProgressDescription("Transformation cancelled by user");
                        datasetTask.setProgress(1.0f);
                        datasetTask.setFinished();
                        return;
                    }
                    QString mode = paramDialog.selectedMode();
                    float selectedValue = static_cast<float>(paramDialog.selectedValue());

                    // Check if a value is actually selected (the dialog allows deselection)
                    bool valueSelected = false;
                    {
                        // The dialog returns defaultValue if nothing is selected, but we can check if the selection matches default
                        // or, more robustly, check if the value is in the list (since defaultValue is always in the list)
                        // For now, we assume if the user didn't select, it returns defaultValue, which is always valid.
                        valueSelected = true;
                    }

                    if (!valueSelected) {
                        datasetTask.setProgressDescription("No value selected for transformation");
                        datasetTask.setProgress(1.0f);
                        datasetTask.setFinished();
                        return;
                    }

                    if (mode == "Extract by value") {
                        // Only extract points exactly equal to the selected value
                        for (int temp = 0; temp < dimensionData.size(); temp++) {
                            if (std::abs(dimensionData.at(temp) - selectedValue) < 1e-5f) {
                                partition1.push_back(temp);
                            }
                        }
                        if (!partition1.empty()) {
                            _transformationNumber = 1;
                            _splitNameSelection = QString("EqualTo%1").arg(selectedValue, 0, 'f', 2);
                            _splitIndices = partition1;
                            _splitIndicesMap.clear();
                            for (int j = 0; j < _splitIndices.size(); j++) {
                                _splitIndicesMap.insert({ _splitIndices[j], j });
                            }
                            datasetTask.setProgressDescription(
                                QString("Processing points == %1 in %2").arg(selectedValue, 0, 'f', 2).arg(dimensions.at(i))
                            );
                            createDatasets();
                        }
                        break;
                    }
                    else if (mode == "Split by value") {
                        // Partition into > and <= selected value
                        for (int temp = 0; temp < dimensionData.size(); temp++) {
                            if (dimensionData.at(temp) > selectedValue) {
                                partition1.push_back(temp);
                            }
                            else {
                                partition2.push_back(temp);
                            }
                        }

                        // Process partition1 (greater than selectedValue)
                        if (!partition1.empty()) {
                            _transformationNumber = 1;
                            _splitNameSelection = QString("GreaterThan%1").arg(selectedValue, 0, 'f', 2);
                            _splitIndices = partition1;
                            _splitIndicesMap.clear();
                            for (int j = 0; j < _splitIndices.size(); j++) {
                                _splitIndicesMap.insert({ _splitIndices[j], j });
                            }
                            datasetTask.setProgressDescription(
                                QString("Processing points > %1 in %2").arg(selectedValue, 0, 'f', 2).arg(dimensions.at(i))
                            );
                            createDatasets();
                        }

                        // Process partition2 (less than or equal to selectedValue)
                        if (!partition2.empty()) {
                            _transformationNumber = 0;
                            _splitNameSelection = QString("LessEqualThan%1").arg(selectedValue, 0, 'f', 2);
                            _splitIndices = partition2;
                            _splitIndicesMap.clear();
                            for (int j = 0; j < _splitIndices.size(); j++) {
                                _splitIndicesMap.insert({ _splitIndices[j], j });
                            }
                            datasetTask.setProgressDescription(
                                QString("Processing points <= %1 in %2").arg(selectedValue, 0, 'f', 2).arg(dimensions.at(i))
                            );
                            createDatasets();
                        }
                        break;
                    }
                    else {
                        datasetTask.setProgressDescription("Unknown transformation mode");
                        datasetTask.setProgress(1.0f);
                        datasetTask.setFinished();
                        return;
                    }
                }
            }
        }
    }

    if (!_pointsSplitDataset.isValid()) {
        datasetTask.setProgressDescription(QString("No matching dataset found for %1").arg(_datasetNameSelection));
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    if (_pointsSplitDataset.isValid() && !foundDimension) {
        datasetTask.setProgressDescription(QString("No matching dimension found for %1").arg(_splitNameSelection));
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();
}

void GradientSurferTransformationPlugin::normalizeRows(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Step 1: Get dimension names
    QStringList dimensionNames;
    {
        auto dims = points->getDimensionNames();
        for (const auto& d : dims)
            dimensionNames << d;
    }

    // Step 2: Show dialog (no dimension selection, always all)
    NormalizeRowsDialog dialog;
    if (dialog.exec() != QDialog::Accepted) {
        datasetTask.setProgressDescription("Normalization cancelled by user");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    datasetTask.setRunning();
    // Step 3: Get user selection
    QString method = dialog.selectedMethod();
    bool inplace = dialog.isInplace();
    QString dtype = dialog.selectedDataType();

    // Step 4: Always normalize all values (global normalization)
    int numPoints = points->getNumPoints();
    int numDims = points->getNumDimensions();
    if (numPoints == 0 || numDims == 0) {
        datasetTask.setProgressDescription("No data to normalize");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    // Step 5: Prepare output data
    std::vector<float> data(numPoints * numDims);
    std::vector<int> dimensionIndices;
    for (int i = 0; i < numDims; i++)
        dimensionIndices.push_back(i);
    points->populateDataForDimensions(data, dimensionIndices);

    // Step 6: Compute statistics for global normalization
    float minVal = std::numeric_limits<float>::max();
    float maxVal = std::numeric_limits<float>::lowest();
    float norm = 0.0f;
    float mean = 0.0f;
    float stddev = 0.0f;
    float maxAbs = 0.0f;
    float scale = 0.0f;

    if (!data.empty()) {
        // Min, Max, MaxAbs
        for (float v : data) {
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
            if (std::abs(v) > maxAbs) maxAbs = std::abs(v);
            norm += v * v;
            mean += v;
        }
        norm = std::sqrt(norm);
        mean /= static_cast<float>(data.size());

        // Standard deviation
        float sumSq = 0.0f;
        for (float v : data) {
            sumSq += (v - mean) * (v - mean);
        }
        stddev = std::sqrt(sumSq / static_cast<float>(data.size()));

        // Scale: difference between max and min
        scale = maxVal - minVal;
    }
    else {
        minVal = 0.0f;
        maxVal = 1.0f;
        norm = 0.0f;
        mean = 0.0f;
        stddev = 0.0f;
        maxAbs = 0.0f;
        scale = 1.0f;
    }
    // Call the normalization utility function
    normalizeVector(data, method.toStdString(), minVal, maxVal,  norm,  mean,  stddev,  maxAbs,  scale);

    // Step 8: Output
    if (dtype == "bfloat16") {
        std::vector<biovault::bfloat16_t> outData(data.size());
        for (size_t i = 0; i < data.size(); ++i)
            outData[i] = static_cast<biovault::bfloat16_t>(data[i]);
        if (!inplace) {
            QString newName = points->getGuiName() + "/normalized";
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(outData.data(), numPoints, numDims);
            newPoints->setDimensionNames(points->getDimensionNames());
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(outData.data(), numPoints, numDims);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }
    else {
        if (!inplace) {
            QString newName = points->getGuiName() + "/normalized";
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(data.data(), numPoints, numDims);
            newPoints->setDimensionNames(points->getDimensionNames());
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(data.data(), numPoints, numDims);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }

    datasetTask.setProgressDescription("Normalization complete");
    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();
}


void GradientSurferTransformationPlugin::removeDimensions(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Step 1: Get dimension names
    QStringList dimensionNames;
    {
        auto dims = points->getDimensionNames();
        for (const auto& d : dims)
            dimensionNames << d;
    }

    // Step 2: Show dialog
    RemoveDimensionsDialog dialog(dimensionNames);
    if (dialog.exec() != QDialog::Accepted) {
        datasetTask.setProgressDescription("Dimension removal cancelled by user");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    datasetTask.setRunning();
    // Step 3: Get user selection
    QStringList selectedDims = dialog.selectedDimensions();
    bool keepSelected = dialog.keepSelected();

    if (selectedDims.isEmpty()) {
        datasetTask.setProgressDescription("No dimensions selected");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    // Step 4: Compute indices to keep
    QVector<int> indicesToKeep;
    for (int i = 0; i < dimensionNames.size(); ++i) {
        bool inSelected = selectedDims.contains(dimensionNames[i]);
        if ((keepSelected && inSelected) || (!keepSelected && !inSelected))
            indicesToKeep << i;
    }
    if (indicesToKeep.isEmpty()) {
        datasetTask.setProgressDescription("No dimensions left after selection");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    // Step 5: Create new dataset with selected dimensions
    int numPoints = points->getNumPoints();
    int numDims = indicesToKeep.size();

    std::vector<QString> newDimNames;
    for (int idx : indicesToKeep)
        newDimNames.push_back(dimensionNames[idx]);

    QString selectedType = dialog.selectedDataType();

    if (selectedType == "bfloat16")
    {
        std::vector<biovault::bfloat16_t> newData(numPoints * numDims);
        points->populateDataForDimensions(newData, indicesToKeep);

        if (!dialog.isInplace()) {
            QString newName = points->getGuiName() + (keepSelected ? "/kept_dims" : "/removed_dims");
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(newData.data(), numPoints, numDims);
            newPoints->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(newData.data(), numPoints, numDims);
            points->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }

    else // default: float
    {
        std::vector<float> newData(numPoints * numDims);
        points->populateDataForDimensions(newData, indicesToKeep);

        if (!dialog.isInplace()) {
            QString newName = points->getGuiName() + (keepSelected ? "/kept_dims" : "/removed_dims");
            Dataset<Points> newPoints = mv::data().createDataset("Points", newName);
            newPoints->setData(newData.data(), numPoints, numDims);
            newPoints->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(newPoints);
            mv::events().notifyDatasetDataChanged(newPoints);
        }
        else {
            points->setData(newData.data(), numPoints, numDims);
            points->setDimensionNames(newDimNames);
            mv::events().notifyDatasetAdded(points);
            mv::events().notifyDatasetDataChanged(points);
        }
    }

    datasetTask.setProgressDescription("Dimension removal complete");
    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();
}

void GradientSurferTransformationPlugin::createDatasetsSingleInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Timer for profiling function execution time
    FunctionTimer timer(Q_FUNC_INFO);
    qDebug() << "createDatasetsInit: ENTER";
    _clustersSplitDataset = nullptr;
    _splitIndicesMap.clear();
    _splitIndices.clear();

    datasetTask.setProgress(0.0f);
    datasetTask.setProgressDescription("Searching for selected cluster...");

    bool foundClustersSplitData = false;

    for (const mv::Dataset<Clusters>& child : points->getChildren()) {
        if (child->getDataType() == ClusterType && child->getGuiName() == _datasetNameSelection) {
            _clustersSplitDataset = child;
            auto clusters = _clustersSplitDataset->getClusters();
            for (const auto& cluster : clusters) {
                if (cluster.getName() == _splitNameSelection) {
                    _splitIndices = cluster.getIndices();
                    for (int i = 0; i < _splitIndices.size(); i++) {
                        _splitIndicesMap.insert({ _splitIndices[i], i });
                    }
                    foundClustersSplitData = true;
                    break;
                }
            }
        }
    }

    if (!_clustersSplitDataset.isValid()) {
        datasetTask.setProgressDescription(QString("No clusters found for %1").arg(_datasetNameSelection));
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }
    if (_splitIndicesMap.empty()) {
        datasetTask.setProgressDescription(QString("No indices found for cluster %1").arg(_splitNameSelection));
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

    datasetTask.setProgress(0.5f);
    datasetTask.setProgressDescription(QString("Processing cluster: %1").arg(_splitNameSelection));

    qDebug() << "createDatasetsInit: EXIT";
    createDatasets();

    datasetTask.setProgress(1.0f);
    datasetTask.setProgressDescription("Cluster processing complete.");
    datasetTask.setFinished();
}

void GradientSurferTransformationPlugin::setType(const QString& type)
{
    // Split by "=>"
    QStringList transformationTypeExtract = type.split("==>");
    _transformationType = transformationTypeExtract.first().trimmed();

    // Defensive: check if we have a right-hand side
    if (transformationTypeExtract.size() < 2) {
        _datasetNameSelection.clear();
        _splitNameSelection.clear();
        _transformationNumber = -1;
        return;
    }

    // Split the right-hand side by "-->"
    QStringList parts = transformationTypeExtract.last().split("-->");
    _datasetNameSelection = parts.first().trimmed();

    QString temp = parts.last().trimmed();
    _transformationNumber = temp.split(":").first().toInt();
    _splitNameSelection = temp.split(":").last().trimmed();
}


// =============================================================================
// Plugin Factory 
// =============================================================================

GradientSurferTransformationPluginFactory::GradientSurferTransformationPluginFactory()
{
    setIconByName("barcode");
	getPluginMetadata().setDescription("GradientSurfer transformation plugin");
    getPluginMetadata().setSummary("This GradientSurfer shows how to implement a basic data transformation plugin in ManiVault Studio.");
    getPluginMetadata().setCopyrightHolder({ "BioVault (Biomedical Visual Analytics Unit LUMC - TU Delft)" });
    getPluginMetadata().setAuthors({
    });
    getPluginMetadata().setOrganizations({
        { "LUMC", "Leiden University Medical Center", "https://www.lumc.nl/en/" },
        { "TU Delft", "Delft university of technology", "https://www.tudelft.nl/" }
	});
    getPluginMetadata().setLicenseText("This plugin is distributed under the [LGPL v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) license.");
}

GradientSurferTransformationPlugin* GradientSurferTransformationPluginFactory::produce()
{
    // Return a new instance of the GradientSurfer transformation plugin
    return new GradientSurferTransformationPlugin(this);
}

mv::DataTypes GradientSurferTransformationPluginFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;

    // This GradientSurfer transformation plugin is compatible with points datasets
    supportedTypes.append(PointType);

    return supportedTypes;
}

mv::gui::PluginTriggerActions GradientSurferTransformationPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    mv::gui::PluginTriggerActions pluginTriggerActions;
    const auto numberOfDatasets = datasets.count();

    if (!PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType))
        return pluginTriggerActions;

    // --- Multi-Dataset Normalize Rows Action ---
    if (numberOfDatasets > 1) {
        auto makeAction = [this](const QString& name, const QString& desc, const QIcon& icon, auto&& func) {
            return new mv::gui::PluginTriggerAction(
                const_cast<GradientSurferTransformationPluginFactory*>(this),
                this,
                name,
                desc,
                icon,
                std::forward<decltype(func)>(func)
            );
            };

        pluginTriggerActions << makeAction(
            "GradientSurfer_Multi-Dataset_Point_Normalize",
            "Perform multi dataset normalize rows data transformation",
            QIcon::fromTheme("calculator"),
            [this, datasets](mv::gui::PluginTriggerAction&) {
                auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                pluginInstance->setInputDatasets(datasets);
                pluginInstance->setType("MultiDatasetPointNormalize==>");
                pluginInstance->transformMultiDatasetRowNormalize();
            }
        );
        return pluginTriggerActions;
    }

    if (numberOfDatasets != 1 || datasets.first()->getDataType() != PointType)
        return pluginTriggerActions;

    Dataset<Points> datasetMain = datasets.first();
    if (datasetMain->getNumDimensions() <= 0 || datasetMain->getNumPoints() <= 0)
        return pluginTriggerActions;

    // --- Helper lambdas for action creation ---
    auto makeAction = [this](const QString& name, const QString& desc, const QIcon& icon, auto&& func) {
        return new mv::gui::PluginTriggerAction(
            const_cast<GradientSurferTransformationPluginFactory*>(this),
            this,
            name,
            desc,
            icon,
            std::forward<decltype(func)>(func)
        );
        };

    // --- Dimension Remove Action ---
    pluginTriggerActions << makeAction(
        "GradientSurfer_Dimension_Remove",
        "Perform dimension removal data transformation",
        QIcon::fromTheme("trash"),
        [this, datasetMain](mv::gui::PluginTriggerAction&) {
            auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
            pluginInstance->setInputDataset(datasetMain);
            pluginInstance->setType("DimensionRemove==>");
            pluginInstance->transformDimensionRemove();
        }
    );

    // --- Normalize Rows Action ---
    pluginTriggerActions << makeAction(
        "GradientSurfer_Point_Normalize",
        "Perform normalize rows data transformation",
        QIcon::fromTheme("calculator"),
        [this, datasetMain](mv::gui::PluginTriggerAction&) {
            auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
            pluginInstance->setInputDataset(datasetMain);
            pluginInstance->setType("PointNormalize==>");
            pluginInstance->transformRowNormalize();
        }
    );

    // --- Collect child dataset options ---
    QVector<QPair<QString, QStringList>> clusterOptionTypes, pointOptionTypes;
    const auto children = datasetMain->getChildren();
    bool containsClusterchildren = false;
    for (const auto& child : children) {
        if (child->getDataType() == ClusterType) {
            containsClusterchildren = true;
            Dataset<Clusters> clusterDataset = mv::data().getDataset<Clusters>(child.getDatasetId());
            if (clusterDataset.isValid()) {
                QStringList options{ "0:All","1:Substring"};
                const auto clusters = clusterDataset->getClusters();
                if (clusters.count() < 500) {
                    int idx = 2;
                    for (const auto& cluster : clusters) {
                        options.append(QString::number(idx++) + ":" + cluster.getName());
                    }
                }
                clusterOptionTypes.append({ clusterDataset->getGuiName(), options });
            }
        }
        if (child->getDataType() == PointType) {
            Dataset<Points> pointDataset = mv::data().getDataset<Points>(child.getDatasetId());
            if (pointDataset.isValid()) {
                QStringList options;
                const auto dimensionNames = pointDataset->getDimensionNames();
                for (const auto& name : dimensionNames)
                    options.append(name);
                pointOptionTypes.append({ pointDataset->getGuiName(), options });
            }
        }
    }

    // --- Cluster Split Actions ---
    for (const auto& optionType : clusterOptionTypes) {
        const QString& mainCategory = optionType.first;
        
        // Actions for splitting by cluster name or "All" 
        for (const auto& subOption : optionType.second) {
            QString firstCopy = mainCategory, subCopy = subOption;
            firstCopy.replace("/", " ");
            subCopy.replace("/", " ");
            const QString actionName = QString("GradientSurfer_Cluster_Split/%1/%2").arg(firstCopy, subCopy);

            pluginTriggerActions << makeAction(
                actionName,
                QString("Perform %1 (%2) data transformation").arg(mainCategory, subOption),
                QIcon::fromTheme("object-ungroup"),
                [this, datasetMain, mainCategory, subOption](mv::gui::PluginTriggerAction&) {
                    auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                    pluginInstance->setInputDataset(datasetMain);
                    pluginInstance->setType(QString("ClusterSplit==>%1-->%2").arg(mainCategory, subOption));
                    pluginInstance->transformCluster();
                }
            );
        }

    }

    // --- Point Split Actions ---
    for (const auto& optionType : pointOptionTypes) {
        const QString& mainCategory = optionType.first;
        for (const auto& subOption : optionType.second) {
            QString firstCopy = mainCategory, subCopy = subOption;
            firstCopy.replace("/", " ");
            subCopy.replace("/", " ");
            const QString actionName = QString("GradientSurfer_Point_Split/%1/%2").arg(firstCopy, subCopy);

            pluginTriggerActions << makeAction(
                actionName,
                QString("Perform %1 (%2) data transformation").arg(mainCategory, subOption),
                QIcon::fromTheme("object-ungroup"),
                [this, datasetMain, mainCategory, subOption](mv::gui::PluginTriggerAction&) {
                    auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                    pluginInstance->setInputDataset(datasetMain);
                    pluginInstance->setType(QString("PointSplit==>%1-->%2").arg(mainCategory, subOption));
                    pluginInstance->transformPoint();
                }
            );
        }
    }

    // --- Remove Zero Columns Action ---
    pluginTriggerActions << makeAction(
        "GradientSurfer_Remove_Zero_Columns",
        "Remove all columns with all row values equal to 0",
        QIcon::fromTheme("edit-delete-column"),
        [this, datasetMain](mv::gui::PluginTriggerAction&) {
            auto pluginInstance = dynamic_cast<GradientSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
            pluginInstance->setInputDataset(datasetMain);
            pluginInstance->transformRemoveZeroColumns();
        }
    );

    return pluginTriggerActions;
}

void GradientSurferTransformationPlugin::createDatasets()
{
    // Timer for profiling function execution time
    FunctionTimer timer(Q_FUNC_INFO);
    qDebug() << "createDataLatest: ENTER";

    // List to keep track of datasets that need to be notified after processing
    mv::Datasets datasetsToNotify;

    // === Step 1: Process the main Points dataset ===
    qDebug() << "Step 1 - Process main Points dataset";

    // Retrieve the main input Points dataset
    Dataset<Points> inputPointsDataset = getInputDataset<Points>();
    auto childDatasets = inputPointsDataset->getChildren();
    qDebug() << "Number of child datasets =" << childDatasets.size();

    int numDimensions = inputPointsDataset->getNumDimensions();
    auto dimensionNames = inputPointsDataset->getDimensionNames();

    // Prepare a vector with all dimension indices
    std::vector<int> allDimensionIndices(numDimensions);
    std::iota(allDimensionIndices.begin(), allDimensionIndices.end(), 0);

    // Construct a descriptive name for the new dataset
    QString newDatasetName = QString::number(_transformationNumber) +"." 
        + inputPointsDataset->getGuiName() + "/" 
        //+ _datasetNameSelection + "/" 
        + _splitNameSelection;

    // Create a new Points dataset for the selected cluster
    Dataset<Points> clusterPointsDataset = mv::data().createDataset("Points", newDatasetName);
    events().notifyDatasetAdded(clusterPointsDataset);

    // Extract and set the data for the selected cluster indices
    std::vector<float> clusterPointsData(_splitIndices.size() * numDimensions);
    inputPointsDataset->populateDataForDimensions(clusterPointsData, allDimensionIndices, _splitIndices);
    clusterPointsDataset->setData(clusterPointsData.data(), _splitIndices.size(), numDimensions);
    clusterPointsDataset->setDimensionNames(dimensionNames);
    datasetsToNotify.push_back(clusterPointsDataset);


    qDebug() << "Step 1 - Finished processing main Points dataset";

    // === Step 2: Process child datasets (Points and Clusters) ===
    qDebug() << "Step 2 - Process child datasets";

    int childIndex = 0;
    for (const Dataset<Clusters>& child : childDatasets) {
        qDebug() << "Processing child" << childIndex
                 << "name:" << child->getGuiName()
                 << "type:" << child->getDataType().getTypeString();

        // If the child is a Points dataset, extract and create a corresponding subset
        if (child->getDataType() == PointType) {
            Dataset<Points> fullChildPoints = child->getFullDataset<Points>();
            if (!fullChildPoints.isValid()) {
                ++childIndex;
                continue;
            }

            // Create a new Points dataset for the child, as a subset of the cluster
            /*Dataset<Points> childClusterPoints = mv::data().createDataset(
                "Points",
                clusterPointsDataset->getGuiName() + "/" + child->getGuiName(),
                clusterPointsDataset
            );*/
            Dataset<Points> childClusterPoints = mv::data().createDerivedDataset(
                //clusterPointsDataset->getGuiName() + "/" + 
                //QString::number(_transformationNumber) + "."+
                //_splitNameSelection+"/" +
                child->getGuiName(),
                clusterPointsDataset
            );
            events().notifyDatasetAdded(childClusterPoints);

            int childNumDimensions = fullChildPoints->getNumDimensions();
            std::vector<int> childDimensionIndices(childNumDimensions);
            std::iota(childDimensionIndices.begin(), childDimensionIndices.end(), 0);

            std::vector<float> childClusterData(_splitIndices.size() * childNumDimensions);
            fullChildPoints->populateDataForDimensions(childClusterData, childDimensionIndices, _splitIndices);
            childClusterPoints->setData(childClusterData.data(), _splitIndices.size(), childNumDimensions);
            childClusterPoints->setDimensionNames(fullChildPoints->getDimensionNames());
            datasetsToNotify.push_back(childClusterPoints);

            qDebug() << "Finished processing point-type child" << childIndex;
        }
        // If the child is a Clusters dataset, create a corresponding subset of clusters
        else if (child->getDataType() == ClusterType) {
            qDebug() << "Processing cluster-type child" << childIndex;
            Dataset<Clusters> fullChildClusters = child->getFullDataset<Clusters>();
            if (!fullChildClusters.isValid()) {
                ++childIndex;
                continue;
            }

            // Create a new Clusters dataset for the child, as a subset of the cluster
            Dataset<Clusters> childClusterDataset = mv::data().createDataset(
                "Cluster",
                //clusterPointsDataset->getGuiName() + "/" +
                //QString::number(_transformationNumber) + "." +
                //_splitNameSelection + "/" +
                child->getGuiName(),
                clusterPointsDataset
            );
            events().notifyDatasetAdded(childClusterDataset);

            // For each cluster in the child, remap indices to the new cluster subset
            for (const auto& cluster : fullChildClusters->getClusters()) {
                std::vector<std::seed_seq::result_type> remappedIndices;
                const auto& originalIndices = cluster.getIndices();
                for (int idx : originalIndices) {
                    // Only include indices that are present in the selected cluster
                    if (auto it = _splitIndicesMap.find(idx); it != _splitIndicesMap.end()) {
                        remappedIndices.push_back(it->second);
                    }
                }
                Cluster remappedCluster = cluster;
                remappedCluster.setIndices(remappedIndices);
                childClusterDataset->addCluster(remappedCluster);
            }

            datasetsToNotify.push_back(childClusterDataset);

            qDebug() << "Finished processing cluster-type child" << childIndex;
        }
        ++childIndex;
    }

    qDebug() << "Step 2 - Finished processing child datasets";
    // === Step 3: Notify all datasets that were created or modified ===
    qDebug() << "Step 3 - Notify datasets";
    for (const auto& dataset : datasetsToNotify) {
        events().notifyDatasetDataChanged(dataset);
    }

    qDebug() << "createDataLatest: EXIT";
}
#include "LayerSurferTransformationPlugin.h"
#include "LayerSurferTransformationDialogs.h"
#include "LayerSurferTransformationUtils.h"

#include <QDebug>
#include <QtCore>
#include <QElapsedTimer>
#include <QtConcurrent>
#include <cmath>
#include<QInputDialog>
#include <QApplication>


Q_PLUGIN_METADATA(IID "studio.manivault.LayerSurferTransformationPlugin")

using namespace mv;
using namespace mv::util;

LayerSurferTransformationPlugin::LayerSurferTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _datasetNameSelection(""),
    _splitNameSelection(""),
    _transformationType(""),
    _transformationNumber(-1)
{
    qApp->setStyleSheet("QToolTip { color: black; background: #ffffe1; border: 1px solid black; }");
}

void LayerSurferTransformationPlugin::transform()
{

}

void LayerSurferTransformationPlugin::transformPoint()
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

void LayerSurferTransformationPlugin::transformDimensionRemove()
{
    mv::Dataset<Points> points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    datasetTask.setName("Transforming");
    datasetTask.setRunning();


    datasetTask.setProgressDescription(QString("Removing dimensions."));
    qDebug() << "Transforming dataset";

    removeDimensions(points, datasetTask);




}

void LayerSurferTransformationPlugin::transformRowNormalize()
{
    mv::Dataset<Points> points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();

    datasetTask.setName("Transforming");
    datasetTask.setRunning();


    datasetTask.setProgressDescription(QString("Normalizing rows."));
    qDebug() << "Transforming dataset";

    normalizeRows(points, datasetTask);


}


void LayerSurferTransformationPlugin::transformRemoveZeroColumns()
{
    mv::Dataset<Points> points = getInputDataset<Points>();
    if (!points.isValid())
        return;
    // Get reference to dataset task for reporting progress
    mv::DatasetTask& datasetTask = points->getTask();
    datasetTask.setName("Transforming");
    datasetTask.setRunning();
    datasetTask.setProgressDescription(QString("Removing zero columns."));
    qDebug() << "Transforming dataset";
    removeZeroColumns(points, datasetTask);

}

void LayerSurferTransformationPlugin::removeZeroColumns(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
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
            for (int i = 0; i < numPoints; ++i) {
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

void LayerSurferTransformationPlugin::transformCluster()
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
    if (_splitNameSelection == "All")
    {
        createDatasetsMultInitCluster(points, datasetTask);
    }
    else
    {
        createDatasetsSingleInitCluster(points, datasetTask);
    }
    


}
void LayerSurferTransformationPlugin::createDatasetsMultInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
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

void LayerSurferTransformationPlugin::createDatasetsPointSplit(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
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

void LayerSurferTransformationPlugin::normalizeRows(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
{
    // Step 1: Get dimension names
    QStringList dimensionNames;
    {
        auto dims = points->getDimensionNames();
        for (const auto& d : dims)
            dimensionNames << d;
    }

    // Step 2: Show dialog (no dimension selection, always all)
    NormalizeRowsDialog dialog(dimensionNames);
    if (dialog.exec() != QDialog::Accepted) {
        datasetTask.setProgressDescription("Normalization cancelled by user");
        datasetTask.setProgress(1.0f);
        datasetTask.setFinished();
        return;
    }

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
    float norm = 1.0f, mean = 0.0f, stddev = 1.0f, minVal = 0.0f, maxVal = 1.0f, maxAbs = 1.0f, scale = 1.0f;
    if (method.startsWith("L2")) {
        // L2 norm of the whole data vector
        norm = std::sqrt(std::accumulate(data.begin(), data.end(), 0.0f, [](float a, float b) { return a + b * b; }));
        if (norm < 1e-8f) norm = 1.0f;
    }
    else if (method.startsWith("L1")) {
        norm = std::accumulate(data.begin(), data.end(), 0.0f, [](float a, float b) { return a + std::abs(b); });
        if (norm < 1e-8f) norm = 1.0f;
    }
    else if (method.startsWith("Max")) {
        norm = 0.0f;
        for (float v : data) norm = std::max(norm, std::abs(v));
        if (norm < 1e-8f) norm = 1.0f;
    }
    else if (method.startsWith("Z-Score")) {
        // mean and stddev of all values
        double sum = 0.0, sqsum = 0.0;
        for (float v : data) { sum += v; sqsum += v * v; }
        mean = static_cast<float>(sum / data.size());
        stddev = static_cast<float>(std::sqrt(sqsum / data.size() - mean * mean));
        if (stddev < 1e-8f) stddev = 1.0f;
    }
    else if (method.startsWith("Min-Max")) {
        auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
        minVal = (minIt != data.end()) ? *minIt : 0.0f;
        maxVal = (maxIt != data.end()) ? *maxIt : 1.0f;
        if (std::abs(maxVal - minVal) < 1e-8f) { minVal = 0.0f; maxVal = 1.0f; }
    }
    else if (method.startsWith("Decimal Scaling")) {
        // Find max absolute value, then scale by 10^j where j = ceil(log10(maxAbs))
        maxAbs = 0.0f;
        for (float v : data) maxAbs = std::max(maxAbs, std::abs(v));
        if (maxAbs < 1e-8f) maxAbs = 1.0f;
        scale = std::pow(10.0f, std::ceil(std::log10(maxAbs)));
        if (scale < 1e-8f) scale = 1.0f;
    }

    // Step 7: Apply normalization to all values
    for (float& v : data) {
        if (method.startsWith("L2") || method.startsWith("L1") || method.startsWith("Max")) {
            v = v / norm;
        }
        else if (method.startsWith("Z-Score")) {
            v = (v - mean) / stddev;
        }
        else if (method.startsWith("Min-Max")) {
            v = (v - minVal) / (maxVal - minVal);
        }
        else if (method.startsWith("Decimal Scaling")) {
            v = v / scale;
        }
    }

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


void LayerSurferTransformationPlugin::removeDimensions(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
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

void LayerSurferTransformationPlugin::createDatasetsSingleInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask)
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

void LayerSurferTransformationPlugin::setType(const QString& type)
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

LayerSurferTransformationPluginFactory::LayerSurferTransformationPluginFactory()
{
    setIconByName("barcode");
	getPluginMetadata().setDescription("LayerSurfer transformation plugin");
    getPluginMetadata().setSummary("This layerSurfer shows how to implement a basic data transformation plugin in ManiVault Studio.");
    getPluginMetadata().setCopyrightHolder({ "BioVault (Biomedical Visual Analytics Unit LUMC - TU Delft)" });
    getPluginMetadata().setAuthors({
    });
    getPluginMetadata().setOrganizations({
        { "LUMC", "Leiden University Medical Center", "https://www.lumc.nl/en/" },
        { "TU Delft", "Delft university of technology", "https://www.tudelft.nl/" }
	});
    getPluginMetadata().setLicenseText("This plugin is distributed under the [LGPL v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) license.");
}

LayerSurferTransformationPlugin* LayerSurferTransformationPluginFactory::produce()
{
    // Return a new instance of the layerSurfer transformation plugin
    return new LayerSurferTransformationPlugin(this);
}

mv::DataTypes LayerSurferTransformationPluginFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;

    // This layerSurfer transformation plugin is compatible with points datasets
    supportedTypes.append(PointType);

    return supportedTypes;
}

mv::gui::PluginTriggerActions LayerSurferTransformationPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    mv::gui::PluginTriggerActions pluginTriggerActions;

    const auto numberOfDatasets = datasets.count();

    if (PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {
        if (numberOfDatasets == 1 && datasets.first()->getDataType() == PointType) {
            
            Dataset<Points> datasetMain = datasets.first();
            if (datasetMain->getNumDimensions() > 0 && datasetMain->getNumPoints() > 0)
            {
                const QString removeActionName = QString("LayerSurfer_Dimension_Remove");
                QIcon removeIcon = QIcon::fromTheme("trash");
                auto pluginTriggerActionRemove = new mv::gui::PluginTriggerAction(
                    const_cast<LayerSurferTransformationPluginFactory*>(this),
                    this,
                    removeActionName,
                    QString("Perform dimension removal data transformation"),
                    removeIcon,
                    [this, datasetMain](mv::gui::PluginTriggerAction& pluginTriggerActionRemove) -> void {
                            auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                            pluginInstance->setInputDataset(datasetMain);
                            pluginInstance->setType(QString("DimensionRemove==>"));
                            pluginInstance->transformDimensionRemove();
            
                    }
                );

                pluginTriggerActions << pluginTriggerActionRemove;

                const QString normalizeActionName = QString("LayerSurfer_Point_Normalize");
                QIcon normalizeIcon = QIcon::fromTheme("calculator");
                auto pluginTriggerActionNormalize = new mv::gui::PluginTriggerAction(
                    const_cast<LayerSurferTransformationPluginFactory*>(this),
                    this,
                    normalizeActionName,
                    QString("Perform normalize rows data transformation"),
                    normalizeIcon,
                    [this, datasetMain](mv::gui::PluginTriggerAction& pluginTriggerActionNormalize) -> void {

                            auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                            pluginInstance->setInputDataset(datasetMain);
                            pluginInstance->setType(QString("PointNormalize==>"));
                            pluginInstance->transformRowNormalize();
                    }
                );

                pluginTriggerActions << pluginTriggerActionNormalize;



            }

            auto children = datasetMain->getChildren();
            if (children.count() > 0) {
                QVector<QPair<QString,QStringList>> clusterOptionTypes;
                QVector<QPair<QString,QStringList>> pointOptionTypes;
                
                
                for (const auto& child : children) {
                    if (child->getDataType() == ClusterType) {
                        Dataset<Clusters> clusterDataset = mv::data().getDataset<Clusters>(child.getDatasetId());
                        if (clusterDataset.isValid())
                        {
                            auto clusters = clusterDataset->getClusters();
                            QStringList options;
                            int idx = 1;
                            options.append("0:All");
                            if(clusters.count()<2000)
                            {

                                for (const auto& cluster : clusters)
                                {
                                    QString formatted = cluster.getName();
                                    formatted = QString::number(idx) + ":" + formatted;
                                    options.append(formatted);
                                    idx++;
                                }

                            }
                            QPair<QString, QStringList> optionvals;
                            optionvals.first = clusterDataset->getGuiName();
                            optionvals.second = options;
                            clusterOptionTypes.append(optionvals);
                        }

                    }
                    if (child->getDataType() == PointType)
                    {
                        Dataset<Points> pointDataset = mv::data().getDataset<Points>(child.getDatasetId());
                        if (pointDataset.isValid())
                        {
                            auto dimensionNames = pointDataset->getDimensionNames();
                            QStringList options;
                            
                            for (int i = 0; i < dimensionNames.size(); i++)
                            {
                                options.append(dimensionNames.at(i));
                            }

                            QPair<QString, QStringList> optionvals;
                            optionvals.first = pointDataset->getGuiName();
                            optionvals.second = options;
                            pointOptionTypes.append(optionvals);
                        }

                    }
                }
                if (clusterOptionTypes.size() > 0)
                {
                    for (const auto& optionType : clusterOptionTypes)
                    {
                        // optionType.first: main category
                        // optionType.second: QStringList of sub-options
                        QIcon clusterIcon = QIcon::fromTheme("object-ungroup");
                        for (int i = 0; i < optionType.second.size(); ++i)
                        {
                            const QString& subOption = optionType.second[i];
                            QString firstCopy = optionType.first;
                            QString subCopy = subOption;
                            firstCopy.replace("/", " ");
                            subCopy.replace("/", " ");
                            const QString actionName = QString("LayerSurfer_Cluster_Split/%1/%2").arg(firstCopy, subCopy);

                            auto pluginTriggerActionCluster = new mv::gui::PluginTriggerAction(
                                const_cast<LayerSurferTransformationPluginFactory*>(this),
                                this,
                                actionName,
                                QString("Perform %1 (%2) data transformation").arg(optionType.first, subOption),
                                icon(),
                                // Explicitly capture optionType and subOption by value
                                [this, datasetMain, optionType, subOption](mv::gui::PluginTriggerAction& pluginTriggerActionCluster) -> void {
                                        auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                                        pluginInstance->setInputDataset(datasetMain);
                                        // Use the setter instead of direct member access
                                        pluginInstance->setType(QString("ClusterSplit==>%1-->%2").arg(optionType.first, subOption));
                                        // pluginInstance->setSelection(optionType.first, subOption); // (optional, if implemented)
                                        pluginInstance->transformCluster();

                                }
                            );

                            pluginTriggerActions << pluginTriggerActionCluster;
                        }
                    }
                }

                if (pointOptionTypes.size() > 0)
                {
                    for (const auto& optionType : pointOptionTypes)
                    {
                        QIcon pointIcon = QIcon::fromTheme("object-ungroup");
                        for (int i = 0; i < optionType.second.size(); ++i)
                        {
                            const QString& subOption = optionType.second[i];
                            QString firstCopy = optionType.first;
                            QString subCopy = subOption;
                            firstCopy.replace("/", " ");
                            subCopy.replace("/", " ");
                            const QString actionName = QString("LayerSurfer_Point_Split/%1/%2").arg(firstCopy, subCopy);

                            auto pluginTriggerActionPoint = new mv::gui::PluginTriggerAction(
                                const_cast<LayerSurferTransformationPluginFactory*>(this),
                                this,
                                actionName,
                                QString("Perform %1 (%2) data transformation").arg(optionType.first, subOption),
                                icon(),
                                // Explicitly capture optionType and subOption by value
                                [this, datasetMain, optionType, subOption](mv::gui::PluginTriggerAction& pluginTriggerActionPoint) -> void {

                                        auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                                        pluginInstance->setInputDataset(datasetMain);
                                        // Use the setter instead of direct member access
                                        pluginInstance->setType(QString("PointSplit==>%1-->%2").arg(optionType.first, subOption));
                                        // pluginInstance->setSelection(optionType.first, subOption); // (optional, if implemented)
                                        pluginInstance->transformPoint();
                                    
                                }
                            );

                            pluginTriggerActions << pluginTriggerActionPoint;
                        }
                    }
                }
            }


            // add a trigger to remove all columns that have all row values 0 
            const QString removeZeroColumnsActionName = QString("LayerSurfer_Remove_Zero_Columns");
            QIcon removeZeroColumnsIcon = QIcon::fromTheme("edit-delete-column"); // Use a suitable icon for removing zero columns
            auto pluginTriggerActionRemoveZeroColumns = new mv::gui::PluginTriggerAction(
                const_cast<LayerSurferTransformationPluginFactory*>(this),
                this,
                removeZeroColumnsActionName,
                QString("Remove all columns with all row values equal to 0"),
                removeZeroColumnsIcon,
                [this, datasetMain](mv::gui::PluginTriggerAction& pluginTriggerActionRemoveZeroColumns) -> void {
                    auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                    pluginInstance->setInputDataset(datasetMain);
                    pluginInstance->transformRemoveZeroColumns();
                }
            );
            pluginTriggerActions << pluginTriggerActionRemoveZeroColumns;


        }
    }

    return pluginTriggerActions;
}

void LayerSurferTransformationPlugin::createDatasets()
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
        //+ inputPointsDataset->getGuiName() + "/" 
        + _datasetNameSelection + "/" + _splitNameSelection;

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
                QString::number(_transformationNumber) + "."+
                _splitNameSelection+"/" +
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
                QString::number(_transformationNumber) + "." +
                _splitNameSelection + "/" +
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
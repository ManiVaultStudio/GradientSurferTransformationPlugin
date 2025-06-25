#include "LayerSurferTransformationPlugin.h"


#include <QDebug>
#include <QtCore>
#include <QElapsedTimer>
#include <QtConcurrent>
#include <cmath>
#include<QInputDialog>

Q_PLUGIN_METADATA(IID "studio.manivault.LayerSurferTransformationPlugin")

using namespace mv;
using namespace mv::util;

class FunctionTimer {
public:
    FunctionTimer(const QString& functionName)
        : _functionName(functionName)
    {
        _timer.start();
    }
    ~FunctionTimer()
    {
        qDebug() << _functionName << "took"
            << _timer.elapsed() / 1000.0 << "seconds";
    }
private:
    QString _functionName;
    QElapsedTimer _timer;
};
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
LayerSurferTransformationPlugin::LayerSurferTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _datasetNameSelection(""),
    _splitNameSelection(""),
    _transformationType(""),
    _transformationNumber(-1)
{
    
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

    datasetTask.setName("Transforming");
    datasetTask.setRunning();

    datasetTask.setProgressDescription(
        QString("Splitting %1 based on %2, dimension %3")
        .arg(points->getGuiName(), _datasetNameSelection, _splitNameSelection)
    );
    createDatasetsPointSplit(points, datasetTask);

    qDebug() << "Transforming dataset";

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

                    bool ok = false;
                    QString label = QString("Select transformation parameter (%1–%2):").arg(minValue).arg(maxValue);
                    float defaultValue = (minValue != maxValue) ? (minValue + maxValue) / 2.0f : minValue;
                    float sliderValue = static_cast<float>(QInputDialog::getDouble(
                        nullptr,
                        "Transformation Parameter",
                        label,
                        defaultValue,
                        minValue,
                        maxValue,
                        2,
                        &ok,
                        Qt::WindowFlags(),
                        0.01
                    ));
                    if (!ok) {
                        datasetTask.setProgressDescription("Transformation cancelled by user");
                        datasetTask.setProgress(1.0f);
                        datasetTask.setFinished();
                        return;
                    }

                    for (int temp = 0; temp < dimensionData.size(); temp++) {
                        if (dimensionData.at(temp) > sliderValue) {
                            partition1.push_back(temp);
                        }
                        else {
                            partition2.push_back(temp);
                        }
                    }

                    // Process partition1 (greater than sliderValue)
                    if (!partition1.empty()) {
                        _transformationNumber = 1;
                        _splitNameSelection = QString("GreaterThan%1").arg(sliderValue);
                        _splitIndices = partition1;
                        _splitIndicesMap.clear();
                        for (int j = 0; j < _splitIndices.size(); j++) {
                            _splitIndicesMap.insert({ _splitIndices[j], j });
                        }
                        datasetTask.setProgressDescription(
                            QString("Processing points > %1 in %2").arg(sliderValue).arg(dimensions.at(i))
                        );
                        createDatasets();
                    }

                    // Process partition2 (less than or equal to sliderValue)
                    if (!partition2.empty()) {
                        _transformationNumber = 0;
                        _splitNameSelection = QString("LessEqualThan%1").arg(sliderValue);
                        _splitIndices = partition2;
                        _splitIndicesMap.clear();
                        for (int j = 0; j < _splitIndices.size(); j++) {
                            _splitIndicesMap.insert({ _splitIndices[j], j });
                        }
                        datasetTask.setProgressDescription(
                            QString("Processing points <= %1 in %2").arg(sliderValue).arg(dimensions.at(i))
                        );
                        createDatasets();
                    }
                    break;
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
    //split type by "->"
    QStringList parts = type.split("->");
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
        { "T. Kroes", { "Lead software architect" }, { "LUMC" } },
        { "J. Thijssen", { "Software architect" }, { "LUMC", "TU Delft" } },
        { "A. Vieth", { "Plugin developer", "Maintainer" }, { "LUMC", "TU Delft" } }
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
            auto children = datasets.first()->getChildren();
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
                            for (const auto& cluster : clusters)
                            {
                                QString formatted = cluster.getName();
                                formatted = QString::number(idx) + ":" + formatted;
                                options.append(formatted);
                                idx++;
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
                    /*if (child->getDataType() == PointType)
                    {
                        Dataset<Points> pointDataset = mv::data().getDataset<Points>(child.getDatasetId());
                        if (pointDataset.isValid())
                        {
                            auto dimensionNames = pointDataset->getDimensionNames();
                            int numofPoints = pointDataset->getNumPoints();
                            //check if values for dimensions are binary
                            QStringList options;
                            for (int i = 0; i < dimensionNames.size(); i++)
                            {
                                std::vector<float> dimensionData(numofPoints);
                                pointDataset->extractDataForDimension(dimensionData,i);
                                //check if std::vector<float> dimensionData contains only binary values
                                if (isBinaryVector(dimensionData)) {
                                    options.append(dimensionNames.at(i));
                                }
                            }
                            QPair<QString, QStringList> optionvals;
                            optionvals.first = pointDataset->getGuiName();
                            optionvals.second = options;
                            pointOptionTypes.append(optionvals);
                        }

                    }*/
                }
                if (clusterOptionTypes.size() > 0)
                {
                    for (const auto& optionType : clusterOptionTypes)
                    {
                        // optionType.first: main category
                        // optionType.second: QStringList of sub-options
                        for (int i = 0; i < optionType.second.size(); ++i)
                        {
                            const QString& subOption = optionType.second[i];
                            QString firstCopy = optionType.first;
                            QString subCopy = subOption;
                            firstCopy.replace("/", " ");
                            subCopy.replace("/", " ");
                            const QString actionName = QString("LayerSurfer_Cluster_Split_Transform/%1/%2").arg(firstCopy, subCopy);

                            auto pluginTriggerAction = new mv::gui::PluginTriggerAction(
                                const_cast<LayerSurferTransformationPluginFactory*>(this),
                                this,
                                actionName,
                                QString("Perform %1 (%2) data transformation").arg(optionType.first, subOption),
                                icon(),
                                // Explicitly capture optionType and subOption by value
                                [this, datasets, optionType, subOption](mv::gui::PluginTriggerAction& pluginTriggerAction) -> void {
                                    for (const auto& dataset : datasets) {
                                        auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                                        pluginInstance->setInputDataset(dataset);
                                        // Use the setter instead of direct member access
                                        pluginInstance->setType(QString("%1->%2").arg(optionType.first, subOption));
                                        // pluginInstance->setSelection(optionType.first, subOption); // (optional, if implemented)
                                        pluginInstance->transformCluster();
                                    }
                                }
                            );

                            pluginTriggerActions << pluginTriggerAction;
                        }
                    }
                }

                if (pointOptionTypes.size() > 0)
                {
                    for (const auto& optionType : pointOptionTypes)
                    {
                        for (int i = 0; i < optionType.second.size(); ++i)
                        {
                            const QString& subOption = optionType.second[i];
                            QString firstCopy = optionType.first;
                            QString subCopy = subOption;
                            firstCopy.replace("/", " ");
                            subCopy.replace("/", " ");
                            const QString actionName = QString("LayerSurfer_Point_Split_Transform/%1/%2").arg(firstCopy, subCopy);

                            auto pluginTriggerAction = new mv::gui::PluginTriggerAction(
                                const_cast<LayerSurferTransformationPluginFactory*>(this),
                                this,
                                actionName,
                                QString("Perform %1 (%2) data transformation").arg(optionType.first, subOption),
                                icon(),
                                // Explicitly capture optionType and subOption by value
                                [this, datasets, optionType, subOption](mv::gui::PluginTriggerAction& pluginTriggerAction) -> void {
                                    for (const auto& dataset : datasets) {
                                        auto pluginInstance = dynamic_cast<LayerSurferTransformationPlugin*>(plugins().requestPlugin(getKind()));
                                        pluginInstance->setInputDataset(dataset);
                                        // Use the setter instead of direct member access
                                        pluginInstance->setType(QString("%1->%2").arg(optionType.first, subOption));
                                        // pluginInstance->setSelection(optionType.first, subOption); // (optional, if implemented)
                                        pluginInstance->transformPoint();
                                    }
                                }
                            );

                            pluginTriggerActions << pluginTriggerAction;
                        }
                    }
                }
            }


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
    QString newDatasetName = QString::number(_transformationNumber) +"." + inputPointsDataset->getGuiName() + "/" + _datasetNameSelection + "/" + _splitNameSelection;

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
            Dataset<Points> childClusterPoints = mv::data().createDataset(
                "Points",
                clusterPointsDataset->getGuiName() + "/" + child->getGuiName(),
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
                clusterPointsDataset->getGuiName() + "/" + child->getGuiName(),
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
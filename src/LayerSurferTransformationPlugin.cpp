#include "LayerSurferTransformationPlugin.h"


#include <QDebug>
#include <QtCore>
#include <QElapsedTimer>
#include <QtConcurrent>
#include <cmath>

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

LayerSurferTransformationPlugin::LayerSurferTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _clusterDatasetNameSelection(""),
    _clusterNameSelection("")
{
    
}

void LayerSurferTransformationPlugin::transform()
{
    auto points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    // Get reference to dataset task for reporting progress
    auto& datasetTask = points->getTask();

    datasetTask.setName("Transforming");
    datasetTask.setRunning();
 
    if (_clusterDatasetNameSelection.isEmpty()|| _clusterNameSelection.isEmpty())
    {
        datasetTask.setProgressDescription("No transformation selected");
        datasetTask.setFinished();
        return;
    }
    datasetTask.setProgressDescription(QString("Splitting %1 based on %2 and cluster %3").arg(points->getGuiName(), _clusterDatasetNameSelection, _clusterNameSelection));
    qDebug() << "Transforming dataset";
    _clustersDataset = nullptr;
    _clusterIndicesMap.clear();
    _clusterIndices.clear();    
    for (const mv::Dataset<Clusters>& child : points->getChildren()) {
        if (child->getDataType() == ClusterType && child->getGuiName() == _clusterDatasetNameSelection) {
            _clustersDataset = child;
            auto clusters = _clustersDataset->getClusters();
            for (const auto& cluster : clusters) {
                if (cluster.getName() == _clusterNameSelection) {
                   
                    _clusterIndices = cluster.getIndices();
                    for (int i = 0; i < _clusterIndices.size();i++) {
                        _clusterIndicesMap.insert({ _clusterIndices[i], i });
                    }
                    break;
                }
            }
        }
    }
    if (!_clustersDataset.isValid()) {
        datasetTask.setProgressDescription(QString("No clusters found for %1").arg(_clusterDatasetNameSelection));
        datasetTask.setFinished();
        return;
    }
    if (_clusterIndicesMap.empty()) {
        datasetTask.setProgressDescription(QString("No indices found for cluster %1").arg(_clusterNameSelection));
        datasetTask.setFinished();
        return;
    }
    createDataLatest();

    datasetTask.setProgress(1.0f);
    datasetTask.setFinished();

}

void LayerSurferTransformationPlugin::setType(const QString& type)
{
    //split type by "->"
    QStringList parts = type.split("->");
    _clusterDatasetNameSelection = parts.first().trimmed();
    _clusterNameSelection = parts.last().trimmed();
}


// =============================================================================
// Plugin Factory 
// =============================================================================

LayerSurferTransformationPluginFactory::LayerSurferTransformationPluginFactory()
{
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
                QVector<QPair<QString,QStringList>> optionTypes;
                for (const Dataset<Clusters>& child : children) {
                    if (child->getDataType() == ClusterType) {
                        auto clusters = child->getClusters();
                        QStringList options;
                        for (const auto& cluster : clusters)
                        {
                            QString formatted = cluster.getName();
                            options.append(formatted);
                        }
                        QPair<QString, QStringList> optionvals;
                        optionvals.first = child->getGuiName();
                        optionvals.second = options;
                        optionTypes.append(optionvals);
                    }
                }
                if (optionTypes.size() > 0)
                {
                    for (const auto& optionType : optionTypes)
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
                            const QString actionName = QString("LayerSurferTransform/%1/%2").arg(firstCopy, subCopy);

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
                                        pluginInstance->transform();
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

void LayerSurferTransformationPlugin::createDataLatest()
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
    QString newDatasetName = inputPointsDataset->getGuiName() + "/" + _clusterDatasetNameSelection + "/" + _clusterNameSelection;

    // Create a new Points dataset for the selected cluster
    Dataset<Points> clusterPointsDataset = mv::data().createDataset("Points", newDatasetName);
    events().notifyDatasetAdded(clusterPointsDataset);

    // Extract and set the data for the selected cluster indices
    std::vector<float> clusterPointsData(_clusterIndices.size() * numDimensions);
    inputPointsDataset->populateDataForDimensions(clusterPointsData, allDimensionIndices, _clusterIndices);
    clusterPointsDataset->setData(clusterPointsData.data(), _clusterIndices.size(), numDimensions);
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

            std::vector<float> childClusterData(_clusterIndices.size() * childNumDimensions);
            fullChildPoints->populateDataForDimensions(childClusterData, childDimensionIndices, _clusterIndices);
            childClusterPoints->setData(childClusterData.data(), _clusterIndices.size(), childNumDimensions);
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
                    if (auto it = _clusterIndicesMap.find(idx); it != _clusterIndicesMap.end()) {
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
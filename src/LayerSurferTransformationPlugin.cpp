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


    /*
    switch (_type)
    {
        // Transform points in place
        case LayerSurferTransformationPlugin::Type::Abs:
        {
            points->setLocked(true);

            points->visitData([this, &points, &datasetTask](auto pointData) {
                std::uint32_t noPointsProcessed = 0;

                for (auto point : pointData) {
                    for (std::uint32_t dimensionIndex = 0; dimensionIndex < points->getNumDimensions(); dimensionIndex++) {
                        point[dimensionIndex] = std::abs(static_cast<double>(point[dimensionIndex]));
                    }

                    ++noPointsProcessed;

                    if (noPointsProcessed % 1000 == 0) {
                        datasetTask.setProgress(static_cast<float>(noPointsProcessed) / static_cast<float>(points->getNumPoints()));

                        QApplication::processEvents();
                    }
                }
                });

            points.setProperty("Last transformed by", getName());
            points->setLocked(false);

            events().notifyDatasetDataChanged(points);

            break;
        }
        // Create new data set
        case LayerSurferTransformationPlugin::Type::Pow2:
        {
            auto derivedData = mv::data().createDerivedDataset<Points>(points->getGuiName() + " (Pow2)", points);

            std::vector<float> transformedData;
            transformedData.resize(points->getNumPoints() * points->getNumDimensions());

            points->constVisitFromBeginToEnd([&transformedData](auto begin, auto end) {
                std::uint32_t noItemsProcessed = 0;

                for (auto it = begin; it != end; ++it) {
                    transformedData[noItemsProcessed] = std::pow(*it, 2.0f);
                    noItemsProcessed++;
                }
                });


            derivedData->setData(transformedData.data(), points->getNumPoints(), points->getNumDimensions());
            events().notifyDatasetDataChanged(derivedData);

            break;
        }
        default:
            break;
    }
    */
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
    FunctionTimer timer(Q_FUNC_INFO);
    qDebug() << "createDataOptimized: ENTER";


    mv::Datasets datasetsToNotify;
    qDebug() << "createDataOptimized: Step 1 - Collect indices";


    qDebug() << "Step 1 - Process main Points dataset";
    

        Dataset<Points> mainDataset = getInputDataset<Points>();
        auto children = mainDataset->getChildren();
        qDebug() << "childrenDatasets.size() =" << children.size();

        int numDim = mainDataset->getNumDimensions();
        auto dimNames = mainDataset->getDimensionNames();
        std::vector<int> allDims(numDim);
        std::iota(allDims.begin(), allDims.end(), 0);
        QString dsName = mainDataset->getGuiName() + "/" + _clusterDatasetNameSelection + "/" + _clusterNameSelection;
        Dataset<Points> ds = mv::data().createDataset("Points", dsName);
        events().notifyDatasetAdded(ds);
        
        std::vector<float>splitData(_clusterIndices.size() * numDim);
        mainDataset->populateDataForDimensions(splitData, allDims, _clusterIndices);
        ds->setData(splitData.data(), _clusterIndices.size(), numDim);
        ds->setDimensionNames(dimNames);
        events().notifyDatasetDataChanged(ds);
        qDebug() << "Step 1 - Finished processing main Points dataset";
        qDebug() << "Step 2 - Collect indices for dataset 1";

        int idx = 0;
        for (const Dataset<Clusters>& child : children) {
            qDebug() << "Processing child" << idx
                << "name:" << child->getGuiName()
                << "type:" << child->getDataType().getTypeString();

            if (child->getDataType() == PointType) {
                Dataset<Points> childFull = child->getFullDataset<Points>();
                if (!childFull.isValid()) continue;

               Dataset<Points> childDataPoint= mv::data().createDataset("Points", ds->getGuiName() + "/" + child->getGuiName(), ds);
                events().notifyDatasetAdded(childDataPoint);
                int numOfDimsChild = childFull->getNumDimensions();
                std::vector<float> splitChildData(_clusterIndices.size() * numOfDimsChild);
                childFull->populateDataForDimensions(splitChildData, allDims, _clusterIndices);
                childDataPoint->setData(splitChildData.data(), _clusterIndices.size(), numOfDimsChild);
                childDataPoint->setDimensionNames(childFull->getDimensionNames());
                events().notifyDatasetDataChanged(childDataPoint);
                qDebug() << "Finished processing point-type child" << idx;

            }
            else if (child->getDataType() == ClusterType) {
                qDebug() << "Processing cluster-type child" << idx;
                Dataset<Clusters> cFull = child->getFullDataset<Clusters>();
                if (!cFull.isValid()) continue;
                

                Dataset<Clusters> clData = mv::data().createDataset("Cluster", ds->getGuiName() + "/" + child->getGuiName(), ds);
                events().notifyDatasetAdded(clData);


                for (const auto& cluster : cFull->getClusters()) {
                    std::vector<std::seed_seq::result_type> i1;
                    const auto& indices = cluster.getIndices();
                    for (int i : indices) {
                        if (auto it = _clusterIndicesMap.find(i); it != _clusterIndicesMap.end())
                        {
                            i1.push_back(it->second);
                        }
  
                    }
                    Cluster a = cluster;
                    a.setIndices(i1);
                    clData->addCluster(a);

                }

                datasetsToNotify.push_back(clData);

                qDebug() << "Finished processing cluster-type child" << idx;
            }
            ++idx;
        }

    

    qDebug() << "createDataOptimized: EXIT";
}
#pragma once

#include <TransformationPlugin.h>
#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include <QMap>
#include <QString>
#include <QDialog>
#include <QVBoxLayout>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QRadioButton>
#include <QLineEdit>
#include <QListWidget>
#include <QCheckBox>
#include <QGroupBox>

/** All plugin related classes are in the ManiVault plugin namespace */
using namespace mv::plugin;

/**
 * LayerSurfer transformation plugin class
 *
 * This transformation plugin class provides skeleton code that shows how to develop 
 * an transformation plugin in ManiVault.
 * 
 * In contrast to analysis plugins, transformation do not create an output data set by default,
 * but operate on the input data set. In this layerSurfer, we provide to transformation options,
 * either simply taking the absolute value of each data point or raising it to the power of 2.
 *
 * To see the plugin in action, please follow the steps below:
 *
 * 1. This plugin works on points datasets and is created by right-clicking a points
 * dataset in the data hierarchy viewer and choosing Transform >> LayerSurfer transformation.
 * 2. Chose the transformation option
 * 
 */

class LayerSurferTransformationPlugin : public TransformationPlugin
{
Q_OBJECT

public:

    /**
     * Constructor
     * @param factory Pointer to the plugin factory
     */
    LayerSurferTransformationPlugin(const PluginFactory* factory);

    /** Destructor */
    ~LayerSurferTransformationPlugin() override = default;

    /** Initialization is called when the plugin is first instantiated. */
    void init() override {};

    /** Performs the data transformation */
    void transform() override;

    void transformCluster();
    void transformPoint();
    void transformRowNormalize();
    void transformDimensionRemove();
    void transformRemoveZeroColumns();
    // Only declare the setter, do not define it here
    void setType(const QString& type);
    void createDatasets();
    void createDatasetsSingleInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
    void createDatasetsMultInitCluster(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
    void createDatasetsPointSplit(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
    void normalizeRows(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
    void removeDimensions(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
    void removeZeroColumns(mv::Dataset<Points>& points, mv::DatasetTask& datasetTask);
private:
    QString    _datasetNameSelection;
    QString     _splitNameSelection;
    QString     _transformationType;
    int         _transformationNumber;
    Dataset<Clusters> _clustersSplitDataset = nullptr;
    std::unordered_map<int, int> _splitIndicesMap;
    std::vector<std::seed_seq::result_type> _splitIndices;
    Dataset<Points>  _pointsSplitDataset = nullptr;

};

/**
 * LayerSurfer transform plugin factory class
 *
 * Note: Factory does not need to be altered (merely responsible for generating new plugins when requested)
 */
class LayerSurferTransformationPluginFactory : public TransformationPluginFactory
{
    Q_INTERFACES(mv::plugin::TransformationPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "studio.manivault.LayerSurferTransformationPlugin"
                      FILE  "PluginInfo.json")

public:

    /** Default constructor */
    LayerSurferTransformationPluginFactory();

    /** Creates an instance of the layerSurfer transform plugin */
    LayerSurferTransformationPlugin* produce() override;

    /** Returns the data types that are supported by the layerSurfer transformation plugin */
    mv::DataTypes supportedDataTypes() const override;

    /** Enable right-click on data set to execute transformation */
    mv::gui::PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;

};

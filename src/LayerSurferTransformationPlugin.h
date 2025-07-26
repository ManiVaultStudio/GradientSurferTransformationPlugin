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

// Custom dialog for transformation parameters, styled like RemoveDimensionsDialog
class TransformationParamDialog : public QDialog {
    Q_OBJECT
public:
    // Add exampleValues parameter
    TransformationParamDialog(float minValue, float maxValue, float defaultValue, const QSet<float>& exampleValues, QWidget* parent = nullptr)
        : QDialog(parent)
        , _minValue(minValue)
        , _maxValue(maxValue)
        , _defaultValue(defaultValue)
    {
        setWindowTitle("Transformation Parameters");
        QVBoxLayout* layout = new QVBoxLayout(this);

        // Mode selection
        layout->addWidget(new QLabel("Choose transformation mode:"));
        modeCombo = new QComboBox(this);
        modeCombo->addItem("Split by value");
        modeCombo->addItem("Extract by value");
        layout->addWidget(modeCombo);

        // Search/filter box for value selection
        layout->addWidget(new QLabel(QString("Search value (%1–%2):").arg(minValue).arg(maxValue)));
        searchEdit = new QLineEdit(this);
        layout->addWidget(searchEdit);

        // Value selection list
        valueList = new QListWidget(this);
        valueList->setSelectionMode(QAbstractItemView::SingleSelection);

        // Populate with example values, sorted
        QList<float> sortedExamples = exampleValues.values();
        std::sort(sortedExamples.begin(), sortedExamples.end());
        for (float v : sortedExamples) {
            valueList->addItem(QString::number(v, 'f', 2));
        }
        // fallback: if no example values, populate with range
        if (valueList->count() == 0) {
            for (float v = minValue; v <= maxValue + 1e-4; v += 0.01f) {
                valueList->addItem(QString::number(v, 'f', 2));
            }
        }
        // Select default value
        for (int i = 0; i < valueList->count(); ++i) {
            if (qFuzzyCompare(valueList->item(i)->text().toFloat(), defaultValue)) {
                valueList->setCurrentRow(i);
                break;
            }
        }
        layout->addWidget(valueList);

        // Deselect button
        deselectButton = new QPushButton("Deselect", this);
        layout->addWidget(deselectButton);
        connect(deselectButton, &QPushButton::clicked, this, [this]() {
            valueList->clearSelection();
            updateInfoLabel();
        });

        // Info label for selection
        infoLabel = new QLabel(this);
        layout->addWidget(infoLabel);
        updateInfoLabel();

        // OK/Cancel buttons
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);

        // Filter values as user types
        connect(searchEdit, &QLineEdit::textChanged, this, [this, sortedExamples, minValue, maxValue](const QString& text) {
            valueList->clear();
            // Filter example values
            for (float v : sortedExamples) {
                QString valStr = QString::number(v, 'f', 2);
                if (valStr.contains(text, Qt::CaseInsensitive))
                    valueList->addItem(valStr);
            }
            // fallback: if no example values, filter range
            if (valueList->count() == 0) {
                for (float v = minValue; v <= maxValue + 1e-4; v += 0.01f) {
                    QString valStr = QString::number(v, 'f', 2);
                    if (valStr.contains(text, Qt::CaseInsensitive))
                        valueList->addItem(valStr);
                }
            }
            // Try to select the default value if present
            for (int i = 0; i < valueList->count(); ++i) {
                if (qFuzzyCompare(valueList->item(i)->text().toFloat(), _defaultValue)) {
                    valueList->setCurrentRow(i);
                    break;
                }
            }
            updateInfoLabel();
        });

        // Update info label on selection change
        connect(valueList->selectionModel(), &QItemSelectionModel::selectionChanged, this, [this]() {
            updateInfoLabel();
        });
    }

    QString selectedMode() const { return modeCombo->currentText(); }
    double selectedValue() const {
        auto items = valueList->selectedItems();
        if (!items.isEmpty())
            return items.first()->text().toDouble();
        // fallback to default
        return _defaultValue;
    }

private:
    void updateInfoLabel() {
        int selected = valueList->selectedItems().size();
        int total = valueList->count();
        infoLabel->setText(
            QString("Selected: %1    Not selected: %2    Total: %3")
                .arg(selected)
                .arg(total - selected)
                .arg(total)
        );
    }

    QComboBox* modeCombo;
    QLineEdit* searchEdit;
    QListWidget* valueList;
    QPushButton* deselectButton;
    QLabel* infoLabel;
    float _minValue, _maxValue, _defaultValue;
};

// Custom dialog for removing dimensions with multi-select, search, and radio buttons
class RemoveDimensionsDialog : public QDialog {
    Q_OBJECT
public:
    RemoveDimensionsDialog(const QStringList& dimensions, QWidget* parent = nullptr)
        : QDialog(parent)
        , allDimensions(dimensions)
    {
        setWindowTitle("Remove/Keep Dimensions");
        QVBoxLayout* layout = new QVBoxLayout(this);

        // Search/filter box
        layout->addWidget(new QLabel("Search dimensions:"));
        searchEdit = new QLineEdit(this);
        layout->addWidget(searchEdit);

        // Multi-selection list widget
        layout->addWidget(new QLabel("Select dimensions:"));
        dimensionList = new QListWidget(this);
        dimensionList->addItems(dimensions);
        dimensionList->setSelectionMode(QAbstractItemView::MultiSelection);
        layout->addWidget(dimensionList);

        // Select All / Deselect All buttons
        QHBoxLayout* selectButtonsLayout = new QHBoxLayout();
        selectAllButton = new QPushButton("Select All", this);
        deselectAllButton = new QPushButton("Deselect All", this);
        selectButtonsLayout->addWidget(selectAllButton);
        selectButtonsLayout->addWidget(deselectAllButton);
        layout->addLayout(selectButtonsLayout);

        connect(selectAllButton, &QPushButton::clicked, this, [this]() {
            for (int i = 0; i < dimensionList->count(); ++i)
                dimensionList->item(i)->setSelected(true);
            updateInfoLabel();
            });
        connect(deselectAllButton, &QPushButton::clicked, this, [this]() {
            dimensionList->clearSelection();
            updateInfoLabel();
            });

        // Info label for selection counts
        infoLabel = new QLabel(this);
        layout->addWidget(infoLabel);
        updateInfoLabel();

        // Radio buttons for keep/remove
        QGroupBox* radioGroup = new QGroupBox("Action", this);
        QHBoxLayout* radioLayout = new QHBoxLayout(radioGroup);
        keepRadio = new QRadioButton("Keep selected", this);
        removeRadio = new QRadioButton("Remove selected", this);
        removeRadio->setChecked(true);
        radioLayout->addWidget(keepRadio);
        radioLayout->addWidget(removeRadio);
        layout->addWidget(radioGroup);

        // Inplace/New radio buttons
        QGroupBox* inplaceGroup = new QGroupBox("Output Mode", this);
        QHBoxLayout* inplaceLayout = new QHBoxLayout(inplaceGroup);
        inplaceRadio = new QRadioButton("Inplace", this);
        newRadio = new QRadioButton("New", this);
        newRadio->setChecked(true);
        inplaceLayout->addWidget(inplaceRadio);
        inplaceLayout->addWidget(newRadio);
        layout->addWidget(inplaceGroup);

        // --- Data type radio buttons ---
        QGroupBox* dtypeGroup = new QGroupBox("Data Type", this);
        QHBoxLayout* dtypeLayout = new QHBoxLayout(dtypeGroup);
        bfloat16Radio = new QRadioButton("bfloat16", this);
        floatRadio = new QRadioButton("float", this);
        floatRadio->setChecked(true); // Default selection
        dtypeLayout->addWidget(bfloat16Radio);
        dtypeLayout->addWidget(floatRadio);
        layout->addWidget(dtypeGroup);

        // OK/Cancel buttons
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);

        // Filter dimensions as user types
        connect(searchEdit, &QLineEdit::textChanged, this, [this, dimensions](const QString& text) {
            dimensionList->clear();
            for (const QString& dim : dimensions) {
                if (dim.contains(text, Qt::CaseInsensitive))
                    dimensionList->addItem(dim);
            }
            // Restore selection if possible
            for (int i = 0; i < dimensionList->count(); ++i) {
                QListWidgetItem* item = dimensionList->item(i);
                if (selectedDimsSet.contains(item->text()))
                    item->setSelected(true);
            }
            updateInfoLabel();
            });

        // Update info label on selection change
        connect(dimensionList->selectionModel(), &QItemSelectionModel::selectionChanged, this, [this]() {
            // Track selected items for restoring after filtering
            selectedDimsSet.clear();
            for (QListWidgetItem* item : dimensionList->selectedItems())
                selectedDimsSet.insert(item->text());
            updateInfoLabel();
            });
    }

    QStringList selectedDimensions() const {
        QStringList dims;
        for (QListWidgetItem* item : dimensionList->selectedItems())
            dims << item->text();
        return dims;
    }
    bool keepSelected() const { return keepRadio->isChecked(); }
    bool removeSelected() const { return removeRadio->isChecked(); }
    bool isInplace() const { return inplaceRadio->isChecked(); }
    bool isNew() const { return newRadio->isChecked(); }

    // Returns the selected data type as a string: "bfloat16", "float"
    QString selectedDataType() const {
        if (bfloat16Radio->isChecked()) return "bfloat16";
        if (floatRadio->isChecked()) return "float";
        return "float"; // fallback
    }

private:
    void updateInfoLabel() {
        int selected = dimensionList->selectedItems().size();
        int total = 0;
        // Count visible items for not selected
        for (int i = 0; i < dimensionList->count(); ++i)
            ++total;
        int notSelected = total - selected;
        infoLabel->setText(
            QString("Selected: %1    Not selected: %2").arg(selected).arg(notSelected)
        );
    }

    QLineEdit* searchEdit;
    QListWidget* dimensionList;
    QRadioButton* keepRadio;
    QRadioButton* removeRadio;
    QLabel* infoLabel;
    QStringList allDimensions;
    QSet<QString> selectedDimsSet;
    QPushButton* selectAllButton;
    QPushButton* deselectAllButton;
    QRadioButton* inplaceRadio;
    QRadioButton* newRadio;

    // Data type radio buttons
    QRadioButton* bfloat16Radio;
    QRadioButton* floatRadio;
};

// --- Update NormalizeRowsDialog: add Z-Score, Min-Max, Decimal Scaling ---
class NormalizeRowsDialog : public QDialog {
    Q_OBJECT
public:
    NormalizeRowsDialog(const QStringList& /*dimensions*/, QWidget* parent = nullptr)
        : QDialog(parent)
    {
        setWindowTitle("Normalize Data");
        QVBoxLayout* layout = new QVBoxLayout(this);

        // Info label: normalization is for all values, based on the full dataset
        layout->addWidget(new QLabel("All values will be normalized using the selected method, based on the full dataset."));

        // Normalization method
        layout->addWidget(new QLabel("Normalization method:"));
        methodCombo = new QComboBox(this);
        methodCombo->addItem("L2 (Euclidean)");
        methodCombo->addItem("L1 (Manhattan)");
        methodCombo->addItem("Max");
        methodCombo->addItem("Z-Score");
        methodCombo->addItem("Min-Max");
        methodCombo->addItem("Decimal Scaling");
        layout->addWidget(methodCombo);

        // Output mode
        QGroupBox* inplaceGroup = new QGroupBox("Output Mode", this);
        QHBoxLayout* inplaceLayout = new QHBoxLayout(inplaceGroup);
        inplaceRadio = new QRadioButton("Inplace", this);
        newRadio = new QRadioButton("New", this);
        newRadio->setChecked(true);
        inplaceLayout->addWidget(inplaceRadio);
        inplaceLayout->addWidget(newRadio);
        inplaceGroup->setLayout(inplaceLayout); // FIX: Set layout for group box
        layout->addWidget(inplaceGroup);

        // Data type
        QGroupBox* dtypeGroup = new QGroupBox("Data Type", this);
        QHBoxLayout* dtypeLayout = new QHBoxLayout(dtypeGroup);
        bfloat16Radio = new QRadioButton("bfloat16", this);
        floatRadio = new QRadioButton("float", this);
        floatRadio->setChecked(true);
        dtypeLayout->addWidget(bfloat16Radio);
        dtypeLayout->addWidget(floatRadio);
        dtypeGroup->setLayout(dtypeLayout); // FIX: Set layout for group box
        layout->addWidget(dtypeGroup);

        // OK/Cancel
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);
    }

    QString selectedMethod() const { return methodCombo->currentText(); }
    bool isInplace() const { return inplaceRadio->isChecked(); }
    QString selectedDataType() const {
        if (bfloat16Radio->isChecked()) return "bfloat16";
        return "float";
    }

private:
    QComboBox* methodCombo;
    QRadioButton* inplaceRadio;
    QRadioButton* newRadio;
    QRadioButton* bfloat16Radio;
    QRadioButton* floatRadio;
};

class RemoveZeroColumnsDialog : public QDialog {
    Q_OBJECT
public:
    RemoveZeroColumnsDialog(QWidget* parent = nullptr)
        : QDialog(parent)
    {
        setWindowTitle("Remove Zero Columns");
        QVBoxLayout* layout = new QVBoxLayout(this);
        // Info label
        layout->addWidget(new QLabel("This will remove all columns that contain only zero values."));
        // Output mode
        QGroupBox* inplaceGroup = new QGroupBox("Output Mode", this);
        QHBoxLayout* inplaceLayout = new QHBoxLayout(inplaceGroup);
        inplaceRadio = new QRadioButton("Inplace", this);
        newRadio = new QRadioButton("New", this);
        newRadio->setChecked(true);
        inplaceLayout->addWidget(inplaceRadio);
        inplaceLayout->addWidget(newRadio);
        inplaceGroup->setLayout(inplaceLayout); // FIX: Set layout for group box
        layout->addWidget(inplaceGroup);
        // Data type
        QGroupBox* dtypeGroup = new QGroupBox("Data Type", this);
        QHBoxLayout* dtypeLayout = new QHBoxLayout(dtypeGroup);
        bfloat16Radio = new QRadioButton("bfloat16", this);
        floatRadio = new QRadioButton("float", this);
        floatRadio->setChecked(true);
        dtypeLayout->addWidget(bfloat16Radio);
        dtypeLayout->addWidget(floatRadio);
        dtypeGroup->setLayout(dtypeLayout); // FIX: Set layout for group box
        layout->addWidget(dtypeGroup);
        // OK/Cancel
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);
    }
    bool isInplace() const { return inplaceRadio->isChecked(); }
    QString selectedDataType() const {
        if (bfloat16Radio->isChecked()) return "bfloat16";
        return "float";
    }

private:
    QRadioButton* inplaceRadio;
    QRadioButton* newRadio;
    QRadioButton* bfloat16Radio;
    QRadioButton* floatRadio;
};



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

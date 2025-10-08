#pragma once

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
#include <QPushButton>
#include <QSet>
#include <QStringList>

struct DatasetClusterOptions {
    QString datasetName;
    QString datasetId;
    QList<QPair<QString, QString>> clusterPrimaryKeyOptions; // (name, id)
};
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
        int selected = static_cast<int>(valueList->selectedItems().size());
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
        int selected = static_cast<int>(dimensionList->selectedItems().size());
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
    NormalizeRowsDialog(QWidget* parent = nullptr)
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
        methodCombo->addItem("Z-Score-RowWise");
        methodCombo->addItem("Min-Max");
        methodCombo->addItem("Decimal Scaling");
        methodCombo->addItem("MaxAbs");
        methodCombo->addItem("Log1p");
        methodCombo->addItem("Log10");
        methodCombo->addItem("Log");
        methodCombo->addItem("CPM");
        methodCombo->addItem("CPM_Log1p");
        methodCombo->addItem("CPM_Log1p_ZScore");
        methodCombo->addItem("Log2");
        methodCombo->addItem("Mean");
        methodCombo->addItem("Softmax");
        methodCombo->addItem("Robust");
        methodCombo->addItem("UnitRange");
        methodCombo->addItem("Binarize");
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


class ExtractByClusterSubstringDialog : public QDialog {
    Q_OBJECT
public:
    ExtractByClusterSubstringDialog(QWidget* parent = nullptr)
        : QDialog(parent)
    {
        setWindowTitle("Extract by Cluster Substring");
        QVBoxLayout* layout = new QVBoxLayout(this);
        layout->addWidget(new QLabel("This will extract data based on cluster substring values."));

        // Input for new substring
        QHBoxLayout* inputLayout = new QHBoxLayout();
        substringEdit = new QLineEdit(this);
        QPushButton* addButton = new QPushButton("Add", this);
        inputLayout->addWidget(new QLabel("Enter cluster substring:"));
        inputLayout->addWidget(substringEdit);
        inputLayout->addWidget(addButton);
        layout->addLayout(inputLayout);

        // List of entered substrings
        substringList = new QListWidget(this);
        substringList->setSelectionMode(QAbstractItemView::SingleSelection);
        layout->addWidget(substringList);
        substringList->addItem("Human");
        substringList->addItem("Macaque");
        substringList->addItem("Marmoset");


        // Remove button
        QPushButton* removeButton = new QPushButton("Remove Selected", this);
        layout->addWidget(removeButton);

        // Data type
        QGroupBox* dtypeGroup = new QGroupBox("Data Type", this);
        QHBoxLayout* dtypeLayout = new QHBoxLayout(dtypeGroup);
        bfloat16Radio = new QRadioButton("bfloat16", this);
        floatRadio = new QRadioButton("float", this);
        bfloat16Radio->setChecked(true);
        dtypeLayout->addWidget(bfloat16Radio);
        dtypeLayout->addWidget(floatRadio);
        dtypeGroup->setLayout(dtypeLayout);
        layout->addWidget(dtypeGroup);

        // OK/Cancel
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);

        // Add button logic
        connect(addButton, &QPushButton::clicked, this, [this]() {
            QString text = substringEdit->text().trimmed();
            if (!text.isEmpty() && !containsSubstring(text)) {
                substringList->addItem(text);
                substringEdit->clear();
            }
            });

        // Remove button logic
        connect(removeButton, &QPushButton::clicked, this, [this]() {
            auto items = substringList->selectedItems();
            for (QListWidgetItem* item : items) {
                delete substringList->takeItem(substringList->row(item));
            }
            });
    }

    // Returns all entered substrings as a QStringList
    QStringList enteredSubstrings() const {
        QStringList result;
        for (int i = 0; i < substringList->count(); ++i)
            result << substringList->item(i)->text();
        return result;
    }
    QString selectedDataType() const {
        if (bfloat16Radio->isChecked()) return "bfloat16";
        return "float";
    }

private:
    bool containsSubstring(const QString& str) const {
        for (int i = 0; i < substringList->count(); ++i)
            if (substringList->item(i)->text() == str)
                return true;
        return false;
    }

    QLineEdit* substringEdit;
    QListWidget* substringList;
    QRadioButton* bfloat16Radio;
    QRadioButton* floatRadio;
};

class MergingRowsDialog : public QDialog {
    Q_OBJECT
public:
    MergingRowsDialog(const DatasetClusterOptions& dataset1,
        const DatasetClusterOptions& dataset2,
        QWidget* parent = nullptr)
        : QDialog(parent)
    {
        setWindowTitle("Merge Rows Between Datasets");
        QVBoxLayout* layout = new QVBoxLayout(this);

        // Store dataset IDs for later retrieval
        _toDatasetId1 = dataset1.datasetId;
        _toDatasetId2 = dataset2.datasetId;

        layout->addWidget(new QLabel("Select 'To' and 'From' datasets:"));

        // --- To Dataset selection and cluster key ---
        QGroupBox* toGroup = new QGroupBox("To Dataset", this);
        QVBoxLayout* toLayout = new QVBoxLayout(toGroup);
        QString toLabel1 = QString("%1 [%2]").arg(dataset1.datasetName, dataset1.datasetId);
        QString toLabel2 = QString("%1 [%2]").arg(dataset2.datasetName, dataset2.datasetId);
        toRadio1 = new QRadioButton(toLabel1, this);
        toRadio2 = new QRadioButton(toLabel2, this);
        toRadio1->setChecked(true);
        toLayout->addWidget(toRadio1);
        toLayout->addWidget(toRadio2);

        toClusterCombo1 = new QComboBox(this);
        for (const auto& pair : dataset1.clusterPrimaryKeyOptions) {
            toClusterCombo1->addItem(QString("%1 [%2]").arg(pair.first, pair.second), pair.second);
        }
        toClusterCombo2 = new QComboBox(this);
        for (const auto& pair : dataset2.clusterPrimaryKeyOptions) {
            toClusterCombo2->addItem(QString("%1 [%2]").arg(pair.first, pair.second), pair.second);
        }

        toLayout->addWidget(new QLabel("Cluster primary key:"));
        toLayout->addWidget(toClusterCombo1);
        toLayout->addWidget(toClusterCombo2);
        toClusterCombo2->hide();

        toGroup->setLayout(toLayout);
        layout->addWidget(toGroup);

        // --- From Dataset selection and cluster key ---
        QGroupBox* fromGroup = new QGroupBox("From Dataset", this);
        QVBoxLayout* fromLayout = new QVBoxLayout(fromGroup);
        QString fromLabel1 = QString("%1 [%2]").arg(dataset1.datasetName, dataset1.datasetId);
        QString fromLabel2 = QString("%1 [%2]").arg(dataset2.datasetName, dataset2.datasetId);
        fromRadio1 = new QRadioButton(fromLabel1, this);
        fromRadio2 = new QRadioButton(fromLabel2, this);
        fromRadio2->setChecked(true);
        fromLayout->addWidget(fromRadio1);
        fromLayout->addWidget(fromRadio2);

        fromClusterCombo1 = new QComboBox(this);
        for (const auto& pair : dataset1.clusterPrimaryKeyOptions) {
            fromClusterCombo1->addItem(QString("%1 [%2]").arg(pair.first, pair.second), pair.second);
        }
        fromClusterCombo2 = new QComboBox(this);
        for (const auto& pair : dataset2.clusterPrimaryKeyOptions) {
            fromClusterCombo2->addItem(QString("%1 [%2]").arg(pair.first, pair.second), pair.second);
        }

        fromLayout->addWidget(new QLabel("Cluster primary key:"));
        fromLayout->addWidget(fromClusterCombo1);
        fromLayout->addWidget(fromClusterCombo2);
        fromClusterCombo1->hide();

        fromGroup->setLayout(fromLayout);
        layout->addWidget(fromGroup);

        // --- Radio logic to ensure only one To/From and show correct cluster key ---
        auto updateRadios = [this]() {
            if (toRadio1->isChecked()) {
                fromRadio1->setEnabled(false);
                fromRadio2->setEnabled(true);
                if (fromRadio1->isChecked()) fromRadio2->setChecked(true);
                toClusterCombo1->show();
                toClusterCombo2->hide();
            }
            else {
                fromRadio2->setEnabled(false);
                fromRadio1->setEnabled(true);
                if (fromRadio2->isChecked()) fromRadio1->setChecked(true);
                toClusterCombo1->hide();
                toClusterCombo2->show();
            }
            if (fromRadio1->isChecked()) {
                fromClusterCombo1->show();
                fromClusterCombo2->hide();
            }
            else {
                fromClusterCombo1->hide();
                fromClusterCombo2->show();
            }
            };
        connect(toRadio1, &QRadioButton::toggled, this, updateRadios);
        connect(toRadio2, &QRadioButton::toggled, this, updateRadios);
        connect(fromRadio1, &QRadioButton::toggled, this, updateRadios);
        connect(fromRadio2, &QRadioButton::toggled, this, updateRadios);
        updateRadios();

        // --- Data type radio buttons ---
        QGroupBox* dtypeGroup = new QGroupBox("Data Type", this);
        QHBoxLayout* dtypeLayout = new QHBoxLayout(dtypeGroup);
        bfloat16Radio = new QRadioButton("bfloat16", this);
        floatRadio = new QRadioButton("float", this);
        bfloat16Radio->setChecked(true); // Default: bfloat16
        dtypeLayout->addWidget(bfloat16Radio);
        dtypeLayout->addWidget(floatRadio);
        dtypeGroup->setLayout(dtypeLayout);
        layout->addWidget(dtypeGroup);

        // --- Column keep mode radio buttons ---
        QGroupBox* keepColsGroup = new QGroupBox("Column Merge Mode", this);
        QVBoxLayout* keepColsLayout = new QVBoxLayout(keepColsGroup);
        keepBothRadio = new QRadioButton("Keep both dataset columns", this);
        keepFromRadio = new QRadioButton("Only keep 'From Dataset' columns", this);
        keepFromRadio->setChecked(true); // Default: only keep from dataset columns
        keepColsLayout->addWidget(keepBothRadio);
        keepColsLayout->addWidget(keepFromRadio);
        keepColsGroup->setLayout(keepColsLayout);
        layout->addWidget(keepColsGroup);

        // --- OK/Cancel buttons ---
        QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
        layout->addWidget(buttonBox);
    }

    // Dataset ID getters
    QString selectedToDatasetId() const {
        return toRadio1->isChecked() ? _toDatasetId1 : _toDatasetId2;
    }
    QString selectedFromDatasetId() const {
        return fromRadio1->isChecked() ? _toDatasetId1 : _toDatasetId2;
    }

    // Cluster ID getters
    QString selectedToClusterId() const {
        if (toRadio1->isChecked())
            return toClusterCombo1->currentData().toString();
        else
            return toClusterCombo2->currentData().toString();
    }
    QString selectedFromClusterId() const {
        if (fromRadio1->isChecked())
            return fromClusterCombo1->currentData().toString();
        else
            return fromClusterCombo2->currentData().toString();
    }

    QString selectedDataType() const {
        if (bfloat16Radio->isChecked()) return "bfloat16";
        return "float";
    }
    bool keepBothColumns() const { return keepBothRadio->isChecked(); }
    bool keepFromColumnsOnly() const { return keepFromRadio->isChecked(); }

private:
    // Dataset selection
    QRadioButton* toRadio1;
    QRadioButton* toRadio2;
    QRadioButton* fromRadio1;
    QRadioButton* fromRadio2;

    // Cluster key selection
    QComboBox* toClusterCombo1;
    QComboBox* toClusterCombo2;
    QComboBox* fromClusterCombo1;
    QComboBox* fromClusterCombo2;

    // Store dataset IDs for retrieval
    QString _toDatasetId1, _toDatasetId2;

    // Other options
    QRadioButton* bfloat16Radio;
    QRadioButton* floatRadio;
    QRadioButton* keepBothRadio;
    QRadioButton* keepFromRadio;
};
#ifndef CONTROLWINDOW_H
#define CONTROLWINDOW_H

#include <QMainWindow>
#include <viewer.h>
#include "config.h"
#include <QColorDialog>
#include <QFileDialog>
#include <QSignalMapper>
#include <QDir>

#include <querywindow.h>
#include <plot2dviewer.h>
#include <sstream>
#include <QNetworkReply>

#define HV_VERSION 0.6

#define CONTROLWINDOW_VERBOSE 0

class Viewer;

namespace Ui {
class ControlWindow;
}

class ControlWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit ControlWindow(QWidget *parent = 0);
    ~ControlWindow();
    void loadFileFull(QString reqFilename);
    void loadConfigFile(QString configFilename);
    void updateVariableBtns();
    void updateConfigWindow();
    void updateAxesInfo(double x0, double x1, double y0, double y1);
    void toggleGrid(bool shown);
    void toggleColorbar(bool shown);
    bool toggleLogscale();
    void setColorbarBounds(double **minmax, int var);
    std::vector <QString>getVarnames();

    Config *config;
    Viewer *viewer;
    QueryWindow *qw;
    Plot2DViewer *plt;

public slots:
    void on_fileNextBtn_clicked();
    void on_filePrevBtn_clicked();

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void on_actionE_xit_triggered();
    void on_action_Open_triggered();
    void on_queryCheckbox_clicked();
    void on_displayCheckbox_clicked();
    void on_pushButton_5_clicked();
    void on_bgColorBtn_clicked();
    void on_openFileBtn_clicked();
    void on_pushButton_clicked();
    void topVarBtn_clicked(int i);
    void botVarBtn_clicked(int i);
    void on_mirrorCheckbox_clicked();
    void on_gridLcheckbox_clicked();
    void on_rotateBtn_clicked();
    void on_view1dThetaBtn_clicked();
    void on_view1dRBtn_clicked();
    void on_view1dBtn_clicked();
    void on_mirror4xCheckbox_clicked();
    void on_horizontalSlider_sliderMoved(int position);
    void on_horizontalSlider_2_sliderMoved(int position);
    void on_horizontalSlider_3_sliderMoved(int position);
    void on_pushButton_2_clicked();
    void on_axis_x0_edit_editingFinished();
    void on_axis_x1_edit_editingFinished();
    void on_axis_y0_edit_editingFinished();
    void on_axis_y1_edit_editingFinished();
    void on_colorbarCheckbox_clicked();
    void on_colorbarResetBtn_clicked();
    void on_actionSave_Config_triggered();
    void on_actionLoad_Config_triggered();
    void on_cbarmin_editingFinished();
    void on_cbarmax_editingFinished();
    void on_cboundReset_clicked();
    void on_actionOptions_triggered();
    void on_cbarCycleBtn_clicked();
    void on_pushButton_3_clicked();
    void on_spinBox_editingFinished();
    void on_logScaleCheckbox_clicked();

    void onFileReady(QNetworkReply * reply);

    void on_scaleCheckbox_clicked();

private:
    std::vector <QPushButton *> topBtns;
    std::vector <QPushButton *> botBtns;
    Ui::ControlWindow *ui;
    QColor bgColor;

    bool fileIsOpen;
    int openFileIndex;
    double axes[4];
    double cmap_vals[2];

    QDir home;
    QString openFilename;
    std::vector<QString> filesInHome;
    std::vector<std::string> split(std::string const &input);
};

#endif // CONTROLWINDOW_H

#include "controlwindow.h"
#include "ui_controlwindow.h"

ControlWindow::ControlWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ControlWindow)
{
    ui->setupUi(this);

    this->openFilename = QDir::currentPath();
    this->ui->bgColorBtn->setStyleSheet("background-color: #6b6b6b;");
    this->openFileIndex = -1;

    this->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);
    this->bgColor.setRgb(107, 107, 107);

    this->fileIsOpen = false;

    this->ui->statusBar->setText(QString(""));
    this->ui->statusBar_bottom->setText(QString(""));
    this->ui->statusBar_bottom->setStyleSheet("QLabel { color : \"#ff0000\"; }");

    // Check for updates
    if (true) {
        try {
            QNetworkAccessManager * manager = new QNetworkAccessManager(this);
            connect(manager, SIGNAL(finished(QNetworkReply*)), this, SLOT(onFileReady(QNetworkReply*)) );
            manager->get(QNetworkRequest(QUrl("https://raw.githubusercontent.com/NYU-CAL/HyVis/master/VERSION")));
        } catch (...) { }
    }

    this->ui->scaleCheckbox->hide(); // TODO: flesh out features

    this->resize(0,0);
}

ControlWindow::~ControlWindow()
{
    delete ui;
}

void ControlWindow::loadFileFull(QString reqFilename)
{
    // First we attempt to load the file. Only if we're successful do we proceed...
    if ( !this->viewer->loadFile(reqFilename.toStdString().c_str()) ) return;

    this->openFilename = reqFilename;
    this->fileIsOpen = true;

    // Proceed iff above was successful
    this->ui->fdirEdit->setText( this->openFilename.left(this->openFilename.lastIndexOf("/")) );
    this->ui->fnameEdit->setText( this->openFilename.right(this->openFilename.length() - 1 - this->openFilename.lastIndexOf("/") ) );
    this->updateConfigWindow();
    this->updateVariableBtns();

    this->home = QDir( this->openFilename.left(this->openFilename.lastIndexOf("/")) );

    int i = 0;
    QStringList files = home.entryList();
    int j;
    //for (QString s : files) {
    for(j=0; j<files.length(); j++) {
        QString s = files.at(j);
        QString stripped = s.right(s.length() - s.lastIndexOf(".") );
        if (QString::compare(stripped, QString(".h5")) == 0) {
            this->filesInHome.push_back(s);
            if (QString::compare(s, this->openFilename.right(this->openFilename.length() - 1 - this->openFilename.lastIndexOf("/") )) == 0)
                this->openFileIndex = i; // Woo-ee!
            i++;
        }
    }

    this->ui->spinBox_2->setValue(this->filesInHome.size());
}

void ControlWindow::loadConfigFile(QString configFilename)
{
    std::string line;
    std::vector<std::string> tokens;

    std::string cfgOpenFile, cfgHomeDir;
    int c0 = 0;

    std::ifstream fp(configFilename.toStdString().c_str());
    if (! fp.is_open()) return;

    while (getline(fp, line)) {
        tokens = this->split(line);
        if (tokens.size() < 1) continue;

        if (tokens.at(0).compare("VERSION") == 0) {
            double load_version = QString(tokens.at(1).c_str()).toDouble(); // TODO warning
            Q_UNUSED(load_version);
        } else if (tokens.at(0).compare("BGCOLOR") == 0) {
            int r = QString(tokens.at(1).c_str()).toInt();
            int g = QString(tokens.at(2).c_str()).toInt();
            int b = QString(tokens.at(3).c_str()).toInt();

            QString fmt;
            fmt.sprintf("background-color: #%02x%02x%02x", r,g,b);
            this->bgColor.setRgb(r,g,b);
            this->ui->bgColorBtn->setStyleSheet(fmt);
            this->viewer->setBgColor(r / 255.0, g / 255.0, b / 255.0);

        } else if (tokens.at(0).compare("ROTATE") == 0) {
            if (tokens.at(1).compare("TRUE") == 0) {
                this->viewer->rotate(90);
                this->ui->rotateBtn->setChecked(true);
            } else {
                this->viewer->rotate(0);
                this->ui->rotateBtn->setChecked(false);
            }

        } else if (tokens.at(0).compare("MIRROR") == 0) {
            if (tokens.at(1).compare("TRUE") == 0) this->ui->mirrorCheckbox->setChecked(true);
            else this->ui->mirrorCheckbox->setChecked(false);
            if (tokens.at(2).compare("TRUE") == 0) this->ui->mirror4xCheckbox->setChecked(true);
            else this->ui->mirror4xCheckbox->setChecked(false);
            this->on_mirrorCheckbox_clicked();
            this->on_mirror4xCheckbox_clicked();

        } else if (tokens.at(0).compare("AXES") == 0) {
            double x0 = QString(tokens.at(1).c_str()).toDouble();
            double x1 = QString(tokens.at(2).c_str()).toDouble();
            double y0 = QString(tokens.at(3).c_str()).toDouble();
            double y1 = QString(tokens.at(4).c_str()).toDouble();
            // TODO set
            Q_UNUSED(x0);
            Q_UNUSED(x1);
            Q_UNUSED(y0);
            Q_UNUSED(y1);

        } else if (tokens.at(0).compare("COLORBAR") == 0) {
            if (tokens.at(1).compare("TRUE") == 0) {
                this->viewer->showColorbar(true);
                this->ui->colorbarCheckbox->setChecked(true);
            } else {
                this->viewer->showColorbar(false);
                this->ui->colorbarCheckbox->setChecked(false);
            }

        } else if (tokens.at(0).compare("CBAR_GAMMA") == 0) {
            int gam = QString(tokens.at(1).c_str()).toInt();
            this->ui->horizontalSlider->setValue(gam);
            if (this->fileIsOpen)
                this->viewer->setGammaValue(gam / 99.0);

        } else if (tokens.at(0).compare("CBAR_CENTER") == 0) {
            int center = QString(tokens.at(1).c_str()).toInt();
            this->ui->horizontalSlider_2->setValue(center);
            if (this->fileIsOpen)
                this->viewer->setCenterValue(center / 99.0);

        } else if (tokens.at(0).compare("CBAR_SLOPE") == 0) {
            int slope = QString(tokens.at(1).c_str()).toInt();
            this->ui->horizontalSlider_3->setValue(slope);
            if (this->fileIsOpen)
                this->viewer->setSlopeValue(slope / 99.0 * 4.0 + 1.0);

        } else if (tokens.at(0).compare("OPEN_FILE") == 0) {
            cfgOpenFile = tokens.at(1);
            c0++;

        } else if (tokens.at(0).compare("OPEN_DIR") == 0) {
            cfgHomeDir = tokens.at(1);
            c0++;

        } else {
            fprintf(stderr, "Unhandled: %s\n", line.c_str());
        }
    }

    if (c0 == 2) {
        this->loadFileFull( QString::fromStdString(cfgHomeDir + "/" + cfgOpenFile) );
    }

    fp.close();
}

void ControlWindow::closeEvent(QCloseEvent *event)
{
    if (event == NULL) ; // Gets rid of warning.
    this->viewer->hide();
    this->qw->hide();
    this->plt->hide();
}

void ControlWindow::on_action_Open_triggered()
{
    this->on_openFileBtn_clicked();
}

void ControlWindow::on_actionE_xit_triggered()
{
    qApp->quit();
}

void ControlWindow::on_queryCheckbox_clicked()
{
    if (this->ui->queryCheckbox->isChecked()) {
        this->qw->show();
    } else {
        this->qw->hide();
    }
}

void ControlWindow::on_displayCheckbox_clicked()
{
    if (this->ui->displayCheckbox->isChecked()) {
        this->viewer->show();
    } else {
        this->viewer->hide();
    }
}

void ControlWindow::on_pushButton_5_clicked()
{
    this->viewer->squareScale();
}

void ControlWindow::on_bgColorBtn_clicked()
{
    int r,g,b;
    this->bgColor = QColorDialog::getColor(this->bgColor);
    this->bgColor.getRgb(&r, &g, &b);

    QString fmt;
    fmt.sprintf("background-color: #%02x%02x%02x", r,g,b);
    this->ui->bgColorBtn->setStyleSheet(fmt);

    this->viewer->setBgColor(r / 255.0, g / 255.0, b / 255.0);
}

void ControlWindow::on_openFileBtn_clicked()
{
    this->filesInHome.clear();

    QString reqFilename = QFileDialog::getOpenFileName(this, "Select a file...", this->openFilename, NULL);
    if (reqFilename == "") return;

    this->loadFileFull(reqFilename);
    this->viewer->rescale(1);
}

void ControlWindow::on_pushButton_clicked()
{
    qApp->quit();
}

void ControlWindow::topVarBtn_clicked(int i)
{
    int j;
    //for (QPushButton *b : this->topBtns) {
    for(j=0; j < this->topBtns.size(); j++) {
        QPushButton *b = this->topBtns.at(j);
        b->setStyleSheet("");
    }
    this->topBtns.at(i)->setStyleSheet("color: #ff2200");
    this->viewer->setLeftValue(i);
    this->viewer->requestCboundReset();
}

void ControlWindow::botVarBtn_clicked(int i)
{
    //for (QPushButton *b : this->botBtns) {
    int j;
    for(j=0; j < this->botBtns.size(); j++) {
        QPushButton *b = this->botBtns.at(j);
        b->setStyleSheet("");
    }
    this->botBtns.at(i)->setStyleSheet("color: #ff2200");
    this->viewer->setRightValue(i);
}

void ControlWindow::updateVariableBtns()
{
    // First clear the bar
    //for (QPushButton *b : this->topBtns) this->ui->variableLayoutTop->removeWidget(b);
    int j;
    for (j=0; j<this->topBtns.size(); j++) {
        QPushButton *b = this->topBtns.at(j);
        this->ui->variableLayoutTop->removeWidget(b);
    }
    //for (QPushButton *b : this->botBtns) this->ui->variableLayoutBottom->removeWidget(b);
    for (j=0; j<this->botBtns.size(); j++) {
        QPushButton *b = this->botBtns.at(j);
        this->ui->variableLayoutBottom->removeWidget(b);
    }
    this->topBtns.clear();
    this->botBtns.clear();

    // Then add new buttons
    QSignalMapper *signalMapperTop = new QSignalMapper(this);
    QSignalMapper *signalMapperBot = new QSignalMapper(this);
    int nq = this->viewer->get_nq();
    int i;
    for (i=0; i<nq; ++i) {
        QPushButton *btnTop = new QPushButton( this->config->getName(i) );
        QPushButton *btnBot = new QPushButton( this->config->getName(i) );
        if (! this->ui->mirrorCheckbox->isChecked()) btnBot->setEnabled(false);
        if (this->viewer->leftDisplayVar == -1 && i == 0)
            btnTop->setStyleSheet("color: #ff2200");
        if (this->viewer->rightDisplayVar == -1 && i==0)
            btnBot->setStyleSheet("color: #ff2200");
        if (this->viewer->leftDisplayVar == i)
            btnTop->setStyleSheet("color: #ff2200");
        if (this->viewer->rightDisplayVar == i)
            btnBot->setStyleSheet("color: #ff2200");
        signalMapperTop->setMapping(btnTop, i);
        signalMapperBot->setMapping(btnBot, i);
        connect( btnTop, SIGNAL(clicked()), signalMapperTop, SLOT(map()) );
        connect( btnBot, SIGNAL(clicked()), signalMapperBot, SLOT(map()) );
        this->ui->variableLayoutTop->addWidget(btnTop);
        this->ui->variableLayoutBottom->addWidget(btnBot);
        this->topBtns.push_back(btnTop);
        this->botBtns.push_back(btnBot);
    }

    connect( signalMapperTop, SIGNAL(mapped(int)), SLOT(topVarBtn_clicked(int)) );
    connect( signalMapperBot, SIGNAL(mapped(int)), SLOT(botVarBtn_clicked(int)) );
}

void ControlWindow::updateConfigWindow()
{
    int nq = this->viewer->get_nq();
    this->config->setNq(nq);
    // TODO, update / clear out appropriately
}

void ControlWindow::updateAxesInfo(double x0, double x1, double y0, double y1)
{
    this->axes[0] = x0; this->axes[1] = x1;
    this->axes[2] = y0; this->axes[3] = y1;

    this->ui->axis_x0_edit->setText( QString::number(x0, 'f', 2) );
    this->ui->axis_x1_edit->setText( QString::number(x1, 'f', 2) );
    this->ui->axis_y0_edit->setText( QString::number(y0, 'f', 2) );
    this->ui->axis_y1_edit->setText( QString::number(y1, 'f', 2) );
}

void ControlWindow::toggleGrid(bool shown)
{
    this->ui->gridLcheckbox->setChecked(shown);
}

void ControlWindow::toggleColorbar(bool shown)
{
    // TODO
    Q_UNUSED(shown);
}

bool ControlWindow::toggleLogscale()
{
    if (this->ui->logScaleCheckbox->isChecked()) {
        this->ui->logScaleCheckbox->setChecked(false);
        return false;
    } else {
        this->ui->logScaleCheckbox->setChecked(true);
        return true;
    }
}

void ControlWindow::setColorbarBounds(double **minmax, int var)
{
    this->ui->cbarmin->setText( QString::number(minmax[0][var], 'f', 2) );
    this->ui->cbarmax->setText( QString::number(minmax[1][var], 'f', 2) );
    this->cmap_vals[0] = minmax[0][var];
    this->cmap_vals[1] = minmax[1][var];
}

std::vector<QString> ControlWindow::getVarnames()
{
    int i;
    std::vector<QString> names;
    for (i=0; i<this->viewer->get_nq(); ++i) {
        names.push_back(this->config->getName(i));
    }
    return names;
}

void ControlWindow::on_mirrorCheckbox_clicked()
{
    if (this->ui->mirrorCheckbox->isChecked()) {
        this->viewer->setIsMirrored(true);
        int j;
        //for (QPushButton *b : this->botBtns) b->setEnabled(true);
        for (j=0; j<this->botBtns.size(); j++)
        {
            QPushButton *b = this->botBtns.at(j);
            b->setEnabled(true);
        }
    } else {
        this->viewer->setIsMirrored(false);
        this->viewer->setIsDoubleMirrored(false);
        this->ui->mirror4xCheckbox->setChecked(false);
        int j;
        //for (QPushButton *b : this->botBtns) b->setEnabled(false);
        for (j=0; j<this->botBtns.size(); j++)
        {
            QPushButton *b = this->botBtns.at(j);
            b->setEnabled(false);
        }
    }
}

void ControlWindow::on_gridLcheckbox_clicked()
{
    this->viewer->setGridVisible(true, this->ui->gridLcheckbox->isChecked());
    this->viewer->setGridVisible(false, this->ui->gridLcheckbox->isChecked());
}

void ControlWindow::on_rotateBtn_clicked()
{
    if (this->ui->rotateBtn->isChecked()) {
        this->viewer->rotate(90);
        this->ui->varTopLbl->setText("Left:");
        this->ui->varBotLabel->setText("Right:");
    }
    else {
        this->viewer->rotate(0);
        this->ui->varTopLbl->setText("Top:");
        this->ui->varBotLabel->setText("Bottom:");
    }
}

void ControlWindow::on_fileNextBtn_clicked()
{
    if (this->openFileIndex == -1) return;

    this->openFileIndex += this->ui->spinBox->value();
    if (this->openFileIndex >= (int)this->filesInHome.size()) this->openFileIndex = 0;
    this->openFilename = this->home.path() + QString("/") + this->filesInHome.at(this->openFileIndex);
    this->viewer->loadFile(this->openFilename.toStdString().c_str());
    this->ui->fnameEdit->setText( this->openFilename.right(this->openFilename.length() - 1 - this->openFilename.lastIndexOf("/") ) );

    // Other updates
    this->updateConfigWindow();
    this->updateVariableBtns();
    // TODO (bug fixes & features)
}

void ControlWindow::on_filePrevBtn_clicked()
{
    if (this->openFileIndex == -1) return;

    this->openFileIndex -= this->ui->spinBox->value();
    if (this->openFileIndex < 0) this->openFileIndex = this->filesInHome.size() - 1;
    this->openFilename = this->home.path() + QString("/") + this->filesInHome.at(this->openFileIndex);
    this->viewer->loadFile(this->openFilename.toStdString().c_str());
    this->ui->fnameEdit->setText( this->openFilename.right(this->openFilename.length() - 1 - this->openFilename.lastIndexOf("/") ) );

    // Other updates
    this->updateConfigWindow();
    this->updateVariableBtns();
    // TODO (bug fixes & features)
}

void ControlWindow::on_view1dThetaBtn_clicked()
{
    this->viewer->set1DPlot(true, true);
    this->plt->show();
}

void ControlWindow::on_view1dRBtn_clicked()
{
    this->viewer->set1DPlot(true, false);
    this->plt->show();
}

void ControlWindow::on_view1dBtn_clicked()
{
    this->viewer->set1DPlot(false, false);
    this->plt->hide();
}

void ControlWindow::on_mirror4xCheckbox_clicked()
{
    if (this->ui->mirror4xCheckbox->isChecked()) {
        this->viewer->setIsDoubleMirrored(true);
        this->ui->mirrorCheckbox->setChecked(true);
        int j;
        //for (QPushButton *b : this->botBtns) b->setEnabled(true);
        for (j=0; j<this->botBtns.size(); j++)
        {
            QPushButton *b = this->botBtns.at(j);
            b->setEnabled(true);
        }
    } else {
        this->viewer->setIsDoubleMirrored(false);
    }
}

void ControlWindow::on_horizontalSlider_sliderMoved(int position)
{
    this->viewer->setGammaValue(position / 99.0);
}

void ControlWindow::on_horizontalSlider_2_sliderMoved(int position)
{
    this->viewer->setCenterValue(position / 99.0);
}

void ControlWindow::on_horizontalSlider_3_sliderMoved(int position)
{
    this->viewer->setSlopeValue((position / 99.0) * 4.0 + 1.0);
}

void ControlWindow::on_pushButton_2_clicked()
{
    this->viewer->savePPM( (this->home.absolutePath() + "/" +
                            this->ui->fnameEdit->text().replace(".h5",".ppm")).toStdString().c_str() );
}

void ControlWindow::on_axis_x0_edit_editingFinished()
{
    // This is triggered when focus leaves the field.
    bool OK = false;
    double tmp = this->ui->axis_x0_edit->text().toDouble(&OK);
    if (OK) axes[0] = tmp;
    else this->ui->axis_x0_edit->setText( QString::number(this->axes[0], 'f', 2) );
    this->viewer->setByAxes(this->axes[0], this->axes[1], this->axes[2], this->axes[3]);
}

void ControlWindow::on_axis_x1_edit_editingFinished()
{
    bool OK = false;
    double tmp = this->ui->axis_x1_edit->text().toDouble(&OK);
    if (OK) axes[1] = tmp;
    else this->ui->axis_x1_edit->setText( QString::number(this->axes[1], 'f', 2) );
    this->viewer->setByAxes(this->axes[0], this->axes[1], this->axes[2], this->axes[3]);
}

void ControlWindow::on_axis_y0_edit_editingFinished()
{
    bool OK = false;
    double tmp = this->ui->axis_y0_edit->text().toDouble(&OK);
    if (OK) axes[2] = tmp;
    else this->ui->axis_y0_edit->setText( QString::number(this->axes[2], 'f', 2) );
    this->viewer->setByAxes(this->axes[0], this->axes[1], this->axes[2], this->axes[3]);
}

void ControlWindow::on_axis_y1_edit_editingFinished()
{
    bool OK = false;
    double tmp = this->ui->axis_y1_edit->text().toDouble(&OK);
    if (OK) axes[3] = tmp;
    else this->ui->axis_y1_edit->setText( QString::number(this->axes[3], 'f', 2) );
    this->viewer->setByAxes(this->axes[0], this->axes[1], this->axes[2], this->axes[3]);
}

void ControlWindow::on_colorbarCheckbox_clicked()
{
    this->viewer->showColorbar(this->ui->colorbarCheckbox->isChecked());
}

void ControlWindow::on_colorbarResetBtn_clicked()
{
    this->ui->horizontalSlider->setValue(80);
    this->ui->horizontalSlider_2->setValue(49);
    this->ui->horizontalSlider_3->setValue(0);
    this->viewer->setGammaValue(80.0 / 99.0);
    this->viewer->setCenterValue(49.0 / 99.0);
    this->viewer->setSlopeValue(1.0);
}

void ControlWindow::on_actionSave_Config_triggered()
{
    /*
     * ?? #steps, stepsize, querywindow, 1d locations
     * cutoffs[2], variables[2]
     */

    QString configFilename = QFileDialog::getSaveFileName(this, "File for configuration...", this->openFilename, NULL);
    if (configFilename == "") return;

    std::ofstream fp(configFilename.toStdString().c_str());
    fp << "VERSION " << HV_VERSION << "\n";
    int r,g,b;
    this->bgColor.getRgb(&r, &g, &b);
    fp << "BGCOLOR " << r << " " << g << " " << b << "\n";
    fp << "ROTATE " << (this->ui->rotateBtn->isChecked() ? "TRUE\n":"FALSE\n");
    fp << "MIRROR " << (this->ui->mirrorCheckbox->isChecked() ? "TRUE ":"FALSE ");
    fp << (this->ui->mirror4xCheckbox->isChecked() ? "TRUE\n":"FALSE\n");
    fp << "AXES " << this->axes[0] << " " << this->axes[1] << " " << this->axes[2];
    fp << " " << this->axes[3] << "\n";
    fp << "COLORBAR " << (this->ui->colorbarCheckbox->isChecked() ? "TRUE\n":"FALSE\n");
    fp << "CBAR_GAMMA " << this->ui->horizontalSlider->value();
    fp << "\nCBAR_CENTER " << this->ui->horizontalSlider_2->value();
    fp << "\nCBAR_SLOPE " << this->ui->horizontalSlider_3->value();
    fp << "\nOPEN_FILE " << this->ui->fnameEdit->text().toStdString();
    fp << "\nOPEN_DIR " << this->ui->fdirEdit->text().toStdString() << "\n";
    fp << "\n";
    fp.close();
}

void ControlWindow::on_actionLoad_Config_triggered()
{
    QString configFilename = QFileDialog::getOpenFileName(this, "Configuration file...", this->openFilename, NULL);
    if (configFilename == "") return;

    this->loadConfigFile(configFilename);
}

std::vector<std::string> ControlWindow::split(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret;

    std::copy(std::istream_iterator<std::string>(buffer),
              std::istream_iterator<std::string>(),
              std::back_inserter(ret));
    return ret;
}

void ControlWindow::onFileReady(QNetworkReply *reply)
{
    try {
        int i;
        QStringList rep = QString(reply->readAll()).split(QRegExp("\\s"));
        for (i=0; i<rep.size(); ++i) {
            if ( rep.at(i).compare("VERSION") == 0 ) {
                double cverr = rep.at(i+1).toFloat();
                if ((cverr - HV_VERSION) > 0.01) {
                    this->ui->statusBar_bottom->setText("HV " + rep.at(i+1) + \
                        " is available. Currently HV " + QString::number(HV_VERSION) + ".");
                }
            }
        }
    } catch (...) { }
}

void ControlWindow::on_cbarmin_editingFinished()
{
    bool OK = false;
    double tmp = this->ui->cbarmin->text().toDouble(&OK);
    if (OK) this->cmap_vals[0] = tmp;
    else this->ui->cbarmin->setText( QString::number(this->cmap_vals[0], 'f', 2) );
    this->viewer->setCmapMinmax(this->cmap_vals[0], this->cmap_vals[1]);
}

void ControlWindow::on_cbarmax_editingFinished()
{
    bool OK = false;
    double tmp = this->ui->cbarmax->text().toDouble(&OK);
    if (OK) this->cmap_vals[1] = tmp;
    else this->ui->cbarmax->setText( QString::number(this->cmap_vals[1], 'f', 2) );
    this->viewer->setCmapMinmax(this->cmap_vals[0], this->cmap_vals[1]);

}

void ControlWindow::on_cboundReset_clicked()
{
    this->viewer->requestCboundReset();
    this->on_cbarmin_editingFinished();
    this->on_cbarmax_editingFinished();
}

void ControlWindow::on_actionOptions_triggered()
{
    // this->config->show(); // This implementation breaks things!
}

void ControlWindow::on_cbarCycleBtn_clicked()
{
    this->viewer->cycleCmap();
}

void ControlWindow::on_pushButton_3_clicked()
{
    // Save image routine
    int i;
    for (i=0; i<this->ui->spinBox_2->value(); ++i) {
        this->on_pushButton_2_clicked();
        this->on_fileNextBtn_clicked();
    }
}

void ControlWindow::on_spinBox_editingFinished()
{
    fprintf(stderr,"Editing finished");
}

void ControlWindow::on_logScaleCheckbox_clicked()
{
    if (this->ui->logScaleCheckbox->isChecked()) {
        this->viewer->setLogscale(true);
    } else {
        this->viewer->setLogscale(false);
    }
}

void ControlWindow::on_scaleCheckbox_clicked()
{

}

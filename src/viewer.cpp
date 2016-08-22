#include "viewer.h"

Viewer::Viewer(QWidget *parent) :
    QGLWidget(parent)
{
    int i;
    for (i=0; i<INPUT_LENGTH; ++i) this->input[i] = false;

    this->setMouseTracking(true);
    this->glDoneInit = false;
    this->ctrlwin = NULL;

    // Note: the following are not set using the public
    // methods so as not to paint the scene.
    this->pointCenter[0] = 0.5;
    this->pointCenter[1] = 0.5;
    this->pointScale[0] = 0.5;
    this->pointScale[1] = 0.5;
    this->clickPoint[0] = 0.0;
    this->clickPoint[1] = 0.0;
    this->clickCenter[0] = 0.0;
    this->clickCenter[1] = 0.0;

    this->bgColor[0] = 107;
    this->bgColor[1] = 107;
    this->bgColor[2] = 107;

    this->rotateAngle = 90;
    this->isMirrored = false;
    this->isDoubleMirrored = false;
    this->drawColorbar = true;

    this->squareScale();

    this->drawScale = true; // TODO: Default value?

    // Draw box setup
    this->drawBox = false;
    this->drawBoxCoordinates[0] = 0.0;
    this->drawBoxCoordinates[1] = 0.0;
    this->drawBoxCoordinates[2] = 0.0;
    this->drawBoxCoordinates[3] = 0.0;

    // Grid setup
    this->drawGridL = false;
    this->drawGridR = false;
    this->gridColor[0] = 1.0;
    this->gridColor[1] = 1.0;
    this->gridColor[2] = 1.0;

    this->valuesAtPointer = NULL;

    this->leftDisplayVar = -1;
    this->rightDisplayVar = -1;
    this->drawPlotLine = false;
    this->plotIsTheta = true;
    this->plotLineData[0] = 0.0;
    this->plotLineData[1] = 1.0;
    this->plotLineData[2] = 1.0;

}

Viewer::~Viewer()
{
}

bool Viewer::loadFile(const char *filename)
{
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: loadFile()...\n");
    #endif

    if ( !this->filedata.loadFromFile(filename) ) return false;

    if (this->glDoneInit) {
        this->g.setCmapMinmax(this->filedata.get_minmax(), this->filedata.get_nq());

        QVector2D *gp;
        GLuint *ci, *gpi;
        int ngp, nci, ngpi;
        this->filedata.genGridData(&gp, &ngp, &ci, &nci, &gpi, &ngpi);

        this->g.loadGeometry(gp, ngp, this->filedata.get_nc(), ci, nci, 
                                gpi, ngpi);
        free(gp);
        free(ci);
        free(gpi);
    }

    // Updates for the l & r values, if uninitialized
    if (this->leftDisplayVar == -1) this->leftDisplayVar = 0;
    if (this->rightDisplayVar == -1) this->rightDisplayVar = 0;

    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);

    // Updates for the query window
    if (this->qw != NULL) this->qw->updateNumValues(this->filedata.get_nq(), this->ctrlwin->getVarnames());
    if (this->valuesAtPointer != NULL) free(this->valuesAtPointer);
    this->valuesAtPointer = (double *)malloc(this->filedata.get_nq() * sizeof(double));

    if (this->drawPlotLine) this->replot1DPlot();
    this->repaint();


    // Updates for the control window //

    // Update axis position
    this->ctrlwin->updateAxesInfo(
        this->pointCenter[0] - this->pointScale[0],
        this->pointCenter[0] + this->pointScale[0],
        this->pointCenter[1] - this->pointScale[1],
        this->pointCenter[1] + this->pointScale[1]);

    // Update the colorbar scale
    this->ctrlwin->setColorbarBounds(this->filedata.get_minmax(), this->leftDisplayVar);
    
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: loadFile() done.\n");
    #endif

    return true;
}

void Viewer::setByAxes(double x0, double x1, double y0, double y1)
{
    if (this->rotateAngle < 45) {
        this->pointCenter[0] = (x0 + x1) / 2.0;
        this->pointCenter[1] = (y0 + y1) / 2.0;
        this->pointScale[0]  = (x1 - x0) / 2.0;
        this->pointScale[1]  = (y1 - y0) / 2.0;
    } else {
        this->pointCenter[0] = (y0 + y1) / 2.0;
        this->pointCenter[1] = (x0 + x1) / -2.0;
        this->pointScale[0]  = (y1 - y0) / 2.0;
        this->pointScale[1]  = (x1 - x0) / 2.0;
    }

    this->updateOrthoMatrix(false);
    this->repaint();
}

void Viewer::setCenter(double x, double y)
{
    this->pointCenter[0] = x;
    this->pointCenter[1] = y;
    this->updateOrthoMatrix();
    this->repaint();
}

void Viewer::squareScale()
{
    if (this->rotateAngle == 0) {
        double pxPerCodeX = this->pointScale[0] / this->width();
        double pxPerCodeY = this->pointScale[1] / this->height();

        if (pxPerCodeX > pxPerCodeY) {
            this->pointScale[1] = pxPerCodeX * this->height();
        } else {
            this->pointScale[0] = pxPerCodeY * this->width();
        }
    } else {
        double pxPerCodeX = this->pointScale[1] / this->width();
        double pxPerCodeY = this->pointScale[0] / this->height();

        if (pxPerCodeX > pxPerCodeY) {
            this->pointScale[0] = pxPerCodeX * this->height();
        } else {
            this->pointScale[1] = pxPerCodeY * this->width();
        }
    }

    this->updateOrthoMatrix();
    this->repaint();
}

// It may seem peculiar to set this way; however, it makes
// sense, if you think about it. 2^2 = 4 = 2*2 combos.
void Viewer::setGridVisible(bool isLeft, bool value)
{
    if (isLeft) this->drawGridL = value;
    else this->drawGridR = value;
    this->repaint();
}

void Viewer::rotate(double angle)
{
    this->rotateAngle = angle;
    this->squareScale();
    this->updateOrthoMatrix();
    this->repaint();
}

void Viewer::setIsMirrored(bool mirror)
{
    this->isMirrored = mirror;
    if ( !this->filedata.fileLoaded ) return;
    if (mirror) {
        if (this->rightDisplayVar != -1) this->setRightValue(this->rightDisplayVar);
        else this->setRightValue(0);
    }

    this->updateOrthoMatrix();
    this->repaint();
}

void Viewer::setIsDoubleMirrored(bool doublemirrored)
{
    if (doublemirrored) this->setIsMirrored(true);
    this->isDoubleMirrored = doublemirrored;

    this->repaint();
}

void Viewer::initializeGL()
{
    this->initShaders();
    this->program.bind();
    this->updateOrthoMatrix();

    this->g.init();
    this->glDoneInit = true;
}

void Viewer::resizeGL(int w, int h)
{
    glViewport(0, 0, (GLint)w, (GLint)h);
}

void Viewer::paintEvent(QPaintEvent *event)
{
    
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: paintEvent()...\n");
    #endif

    Q_UNUSED(event);

    QPainter p( this );
    p.beginNativePainting();

    // Save gl state
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    // Clear the screen
    qglClearColor(QColor(this->bgColor[0], this->bgColor[1], this->bgColor[2]));
    glClear(GL_COLOR_BUFFER_BIT);

    this->program.bind();
    program.setUniformValue("ortho_matrix", this->matrix);
    this->g.drawGeometryLeft(&this->program);

    if (this->isMirrored) {
        program.setUniformValue("ortho_matrix", this->matrixRight);
        this->g.drawGeometryRight(&this->program);
    }
    if (this->isDoubleMirrored) {
        QMatrix4x4 upLeft = this->matrix;
        QMatrix4x4 upRight = this->matrixRight;
        upLeft.scale(-1.0, 1.0, 1.0);
        upRight.scale(-1.0, 1.0, 1.0);
        program.setUniformValue("ortho_matrix", upLeft);
        this->g.drawGeometryLeft(&this->program);
        program.setUniformValue("ortho_matrix", upRight);
        this->g.drawGeometryRight(&this->program);
    }

    // Draw grid
    if (this->drawGridL) {
        this->whiteprogram.bind();
        this->whiteprogram.setUniformValue("ortho_matrix", this->matrix);
        this->whiteprogram.setUniformValue("const_color", this->gridColor[0], this->gridColor[1], this->gridColor[2]);
        this->g.drawGrid(&this->whiteprogram);
        this->program.bind();
    }
    if (this->drawGridR && this->isMirrored) {
        this->whiteprogram.bind();
        this->whiteprogram.setUniformValue("ortho_matrix", this->matrixRight);
        this->whiteprogram.setUniformValue("const_color", this->gridColor[0], this->gridColor[1], this->gridColor[2]);
        this->g.drawGrid(&this->whiteprogram);
        this->program.bind();
    }

    // Draw red plot line
    if (this->drawPlotLine) {
        this->whiteprogram.bind();
        this->whiteprogram.setUniformValue("ortho_matrix", this->matrix);
        this->whiteprogram.setUniformValue("const_color", 1.0, 0.0, 0.0);
        glLineWidth(1.0);
        glBegin(GL_LINES);
        if (this->plotIsTheta) {
            glVertex3f(this->plotLineData[0] * cos(this->plotLineData[2]),
                       this->plotLineData[0] * sin(this->plotLineData[2]), 0.0);
            glVertex3f(this->plotLineData[1] * cos(this->plotLineData[2]),
                       this->plotLineData[1] * sin(this->plotLineData[2]), 0.0);
        } else {
            int i;
            double t1, t2;
            t1 = this->plotLineData[0];
            for (i=1; i<200; ++i) {
                t2 = this->plotLineData[0] * (1.0 - i / 200.0) + this->plotLineData[1] * (i / 200.0);
                glVertex3f(this->plotLineData[2] * cos(t1), this->plotLineData[2] * sin(t1), 0.0);
                glVertex3f(this->plotLineData[2] * cos(t2), this->plotLineData[2] * sin(t2), 0.0);
                t1 = t2;
            }
        }
        glEnd();
        this->program.bind();
    }

    // Draw colorbar & scale
    if (this->drawColorbar) {
        int i; QVector3D col;
        QMatrix4x4 ident; ident.setToIdentity();
        this->whiteprogram.bind();
        this->whiteprogram.setUniformValue("ortho_matrix", ident);

        // Draw the colorbar itself
        int cbar_res = 200; // <-- Feel free to edit this
        double ymodifier = 1.9 / cbar_res;
        for (i=0; i<cbar_res; ++i) {
            this->g.getCmapValue((1.0 * i) / cbar_res, &col);
            this->whiteprogram.setUniformValue("const_color", col.x(), col.y(), col.z());
            glBegin(GL_QUADS);
            glVertex3f(0.88, -0.95 + i * ymodifier, 0.0); glVertex3f(0.88, -0.95 + (i+1) * ymodifier, 0.0);
            glVertex3f(0.98, -0.95 + (i+1) * ymodifier, 0.0); glVertex3f(0.98, -0.95 + i * ymodifier, 0.0);
            glEnd();
        }

        this->program.bind();
    }

    // Draw zoom box
    if (this->drawBox) {
        this->whiteprogram.bind();
        this->whiteprogram.setUniformValue("ortho_matrix", this->matrix);
        this->whiteprogram.setUniformValue("const_color", 1.0, 1.0, 1.0);
        double a = this->drawBoxCoordinates[0];
        double b = this->drawBoxCoordinates[1];
        double c = this->drawBoxCoordinates[2];
        double d = this->drawBoxCoordinates[3];
        glLineWidth(1.0);
        glBegin(GL_LINES);
        glVertex3f(a, b, 0.0); glVertex3f(a, d, 0.0);
        glVertex3f(a, d, 0.0); glVertex3f(c, d, 0.0);
        glVertex3f(c, d, 0.0); glVertex3f(c, b, 0.0);
        glVertex3f(c, b, 0.0); glVertex3f(a, b, 0.0);
        glEnd();
        this->program.bind();
    }

    // Restore gl state
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    p.endNativePainting();

    // Now we can call normal QPainter methods.

    // Set up painter for text
    p.setRenderHint(QPainter::TextAntialiasing);
    p.setPen(Qt::white);
    QFont f = p.font();
    f.setPointSize(10);
    p.setFont(f);

    // Draw colorbar scale
    if (this->drawColorbar) {
        int i;
        int marks = 5;
        int x1 = 0.935 * width();
        int y0 = height() * 1.95 / 2.0;
        int dy = height() * 1.9 / 2.0 / marks;
        char buf[16];
        double val;
        for (i=0; i<=marks; ++i) {
            this->g.unitToTrue((double)i / (double)marks, &val);
            sprintf(buf, "%.2e", val);
            p.drawText(0, 0, x1, y0-(i*dy)+5, Qt::AlignBottom | Qt::AlignRight, buf);
        }
    }

    // Draw normal scale
    if (this->drawScale) {
        // TODO ... naive method isn't smart.
    }
    
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: paintEvent() done.\n");
    #endif
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    if (this->input[KEY_SHIFT] && !this->input[KEY_CTRL]) {
        this->drawBox = true;
        this->screenToData(event->pos().x(), event->pos().y(),
                           &(this->drawBoxCoordinates[0]),
                           &(this->drawBoxCoordinates[1]));
    } else if (this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
        this->update1DPlot(event->pos().x(), event->pos().y(), true);
    }

    double tx, ty;
    this->pxToPoint(event->pos().x(),event->pos().y(), &tx,&ty);
    switch (event->button()) {
    case Qt::LeftButton:
        this->input[MOUSE_LEFT] = true;
        this->clickPoint[0] = tx; this->clickPoint[1] = ty;
        this->clickCenter[0] = this->pointCenter[0]; this->clickCenter[1] = this->pointCenter[1];
        break;
    case Qt::RightButton:
        this->input[MOUSE_RIGHT] = true;
        this->clickPoint[0] = tx; this->clickPoint[1] = ty;
        this->clickCenter[0] = this->pointCenter[0]; this->clickCenter[1] = this->pointCenter[1];
        break;
    default: ;
    }
}

void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
    if (this->drawBox) {
        double a = this->drawBoxCoordinates[0];
        double b = this->drawBoxCoordinates[1];
        double c = this->drawBoxCoordinates[2];
        double d = this->drawBoxCoordinates[3];

        this->pointScale[0] = fabs(a-c) / 2.0;
        this->pointScale[1] = fabs(b-d) / 2.0;
        this->pointCenter[0] = (a+c) / 2.0;
        this->pointCenter[1] = (b+d) / 2.0;

        this->updateOrthoMatrix();
        this->repaint();

        this->drawBox = false;
    }

    if (this->input[KEY_CTRL]) {
        this->update1DPlot(event->pos().x(), event->pos().y(), true);
    }

    switch (event->button()) {
    case Qt::LeftButton:
        this->input[MOUSE_LEFT] = false;
        break;
    case Qt::RightButton:
        this->input[MOUSE_RIGHT] = false;
        break;
    default: ;
    }
}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
    // In lieu of thinking of a better way to do it, I am adjusting position
    // in two different ways for now...
    double tx,ty;
    this->pxToPoint(event->pos().x(),event->pos().y(), &tx,&ty);

    double ttx, tty;
    this->screenToData(event->pos().x(), event->pos().y(), &ttx, &tty);

    // Update query window
    if (this->isMirrored && tty < 0) {
        tty = -tty;
        this->qw->setPosition(ttx, tty, 0.0);
        if (this->valuesAtPointer != NULL) {
            this->filedata.getValuesAt(ttx, tty, this->valuesAtPointer);
            this->qw->setValues(this->valuesAtPointer);
        }
        tty = -tty;
    } else {
        this->qw->setPosition(ttx, tty, 0.0);
        if (this->valuesAtPointer != NULL) {
            this->filedata.getValuesAt(ttx, tty, this->valuesAtPointer);
            this->qw->setValues(this->valuesAtPointer);
        }
    }

    // Handle input
    if ((this->input[MOUSE_LEFT] || this->input[MOUSE_RIGHT]) && !this->input[KEY_SHIFT] && !this->input[KEY_CTRL]) {
        this->setCenter(this->clickCenter[0]+(this->clickPoint[0]-tx),
                this->clickCenter[1]+(this->clickPoint[1]-ty));
    } else if (this->input[KEY_SHIFT]) {
        this->drawBoxCoordinates[2] = ttx;
        this->drawBoxCoordinates[3] = tty;
        this->repaint();
    } else if (this->input[KEY_CTRL] && !this->input[KEY_SHIFT] && (event->buttons() & Qt::LeftButton)) {
        this->update1DPlot(event->pos().x(), event->pos().y(), false);

    } else {
        // ... ?
    }
}

void Viewer::wheelEvent(QWheelEvent *event)
{
    QPoint deg = event->angleDelta();
    float mod = deg.y() / 120;
    if (mod < 0) {
        this->pointScale[0] *= pow(1.2, -mod);
        this->pointScale[1] *= pow(1.2, -mod);
    } else {
        this->pointScale[0] /= pow(1.2, mod);
        this->pointScale[1] /= pow(1.2, mod);
    }

    this->updateOrthoMatrix();
    this->repaint();
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_W:
        this->input[KEY_W] = true;
        if (this->rotateAngle == 0) this->pointCenter[1] += this->pointScale[1] / 10.0;
        else this->pointCenter[0] += this->pointScale[0] / 10.0;
        this->updateOrthoMatrix();
        this->repaint();
        break;
    case Qt::Key_A:
        this->input[KEY_A] = true;
        if (this->rotateAngle == 0) this->pointCenter[0] -= this->pointScale[0] / 10.0;
        else this->pointCenter[1] += this->pointScale[1] / 10.0;
        this->updateOrthoMatrix();
        this->repaint();
        break;
    case Qt::Key_S:
        this->input[KEY_S] = true;
        if (this->rotateAngle == 0) this->pointCenter[1] -= this->pointScale[1] / 10.0;
        else this->pointCenter[0] -= this->pointScale[0] / 10.0;
        this->updateOrthoMatrix();
        this->repaint();
        break;
    case Qt::Key_D:
        this->input[KEY_D] = true;
        if (this->rotateAngle == 0) this->pointCenter[0] += this->pointScale[0] / 10.0;
        else this->pointCenter[1] -= this->pointScale[1] / 10.0;
        this->updateOrthoMatrix();
        this->repaint();
        break;
    case Qt::Key_E:
        this->input[KEY_E] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->pointScale[0] /= 1.2;
            this->pointScale[1] /= 1.2;
            this->updateOrthoMatrix();
            this->repaint();
        }
        break;
    case Qt::Key_C:
        this->input[KEY_C] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->pointScale[0] *= 1.2;
            this->pointScale[1] *= 1.2;
            this->updateOrthoMatrix();
            this->repaint();
        }
        break;
    case Qt::Key_Equal:
        this->input[KEY_PLUS] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->pointScale[0] /= 1.2;
            this->pointScale[1] /= 1.2;
            this->updateOrthoMatrix();
            this->repaint();
        }
        break;
    case Qt::Key_Plus:
        this->input[KEY_PLUS] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->pointScale[0] /= 1.2;
            this->pointScale[1] /= 1.2;
            this->updateOrthoMatrix();
            this->repaint();
        }
        break;
    case Qt::Key_Minus:
        this->input[KEY_MINUS] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->pointScale[0] *= 1.2;
            this->pointScale[1] *= 1.2;
            this->updateOrthoMatrix();
            this->repaint();
        }
        break;

    case Qt::Key_G:
        this->input[KEY_G] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->ctrlwin->toggleGrid(!this->drawGridL);
            this->setGridVisible(true, !this->drawGridL);
            this->setGridVisible(false, !this->drawGridR);
            this->repaint();
        }
        break;
    case Qt::Key_L:
        this->input[KEY_L] = true;
        if (!this->input[KEY_CTRL] && !this->input[KEY_SHIFT]) {
            this->setLogscale(this->ctrlwin->toggleLogscale());
        }
        break;
    case Qt::Key_Q:
        this->input[KEY_Q] = true;
        if (this->input[KEY_CTRL]) {
            exit(42);
        }
        break;
    case Qt::Key_Left:
        this->input[KEY_LEFT] = true;
        break;
    case Qt::Key_Right:
        this->input[KEY_RIGHT] = true;
        break;
    case Qt::Key_Shift:
        this->input[KEY_SHIFT] = true;
        break;
    case Qt::Key_Control: // This appears to be cmd, not ctrl... (OSX);
        this->input[KEY_CTRL] = true;
        break;
    default: ;
    }
}

void Viewer::keyReleaseEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_W:
        this->input[KEY_W] = false;
        break;
    case Qt::Key_A:
        this->input[KEY_A] = false;
        break;
    case Qt::Key_S:
        this->input[KEY_S] = false;
        break;
    case Qt::Key_D:
        this->input[KEY_D] = false;
        break;
    case Qt::Key_E:
        this->input[KEY_E] = false;
        break;
    case Qt::Key_C:
        this->input[KEY_C] = false;
        break;
    case Qt::Key_G:
        this->input[KEY_G] = false;
        break;
    case Qt::Key_Q:
        this->input[KEY_Q] = false;
        break;
    case Qt::Key_Plus:
        this->input[KEY_PLUS] = false;
        break;
    case Qt::Key_Equal:
        this->input[KEY_PLUS] = false;
        break;
    case Qt::Key_Minus:
        this->input[KEY_MINUS] = false;
        break;
    case Qt::Key_Left:
        this->input[KEY_LEFT] = false;
        if (!this->input[KEY_SHIFT] && !this->input[KEY_CTRL]) this->ctrlwin->on_filePrevBtn_clicked();
        break;
    case Qt::Key_Right:
        this->input[KEY_RIGHT] = false;
        if (!this->input[KEY_SHIFT] && !this->input[KEY_CTRL]) this->ctrlwin->on_fileNextBtn_clicked();
        break;
    case Qt::Key_Shift:
        this->input[KEY_SHIFT] = false;
        break;
    case Qt::Key_Control:
        this->input[KEY_CTRL] = false;
        break;
    default: ;
    }
}

void Viewer::update1DPlot(int x, int y, bool isFirst)
{
    if (! this->filedata.fileLoaded ) return;

    int n=0;
    double ttx, tty;
    std::vector<double> xdata, ydata;

    this->screenToData(x, y, &ttx, &tty);
    if (this->plotIsTheta) n = this->filedata.getValuesAtTheta(ttx, tty, this->leftDisplayVar, &xdata, &ydata);
    else n = this->filedata.getValuesAtR(ttx, tty, this->leftDisplayVar, &xdata, &ydata);
    if (n<1) return;

    double y0 = this->filedata.get_minmax()[0][this->leftDisplayVar];
    double y1 = this->filedata.get_minmax()[1][this->leftDisplayVar];
    this->plt->updatePlot(&xdata,&ydata, isFirst, y0,y1, c->getName(this->leftDisplayVar).toStdString());
    this->plotLineData[0] = xdata.at(0);
    this->plotLineData[1] = xdata.at(n-1);

    if (this->plotIsTheta) this->plotLineData[2] = atan(tty/ttx);
    else this->plotLineData[2] = sqrt(tty*tty + ttx*ttx);

    this->repaint();
}


void Viewer::initShaders()
{
    if (!program.addShaderFromSourceFile(QGLShader::Vertex, ":/VShader.vsh")) {
        fprintf(stderr, "Unable to load vertex shader to graphics card.\n");
        exit(-1);
    }
    if (!program.addShaderFromSourceFile(QGLShader::Fragment, ":/FShader.fsh")) {
        fprintf(stderr, "Unable to load fragment shader to graphics card.\n");
        exit(-2);
    }
    if (!program.link()) {
        fprintf(stderr, "Unable to link shader program.\n");
        exit(-3);
    }

    if (!whiteprogram.addShaderFromSourceFile(QGLShader::Vertex, ":/VShader.vsh")) {
        fprintf(stderr, "Unable to load vertex shader to graphics card.\n");
        exit(-1);
    }
    if (!whiteprogram.addShaderFromSourceFile(QGLShader::Fragment, ":/FWhite.fsh")) {
        fprintf(stderr, "Unable to load fragment shader to graphics card.\n");
        exit(-2);
    }
    if (!whiteprogram.link()) {
        fprintf(stderr, "Unable to link shader program.\n");
        exit(-3);
    }
}

void Viewer::updateOrthoMatrix()
{
    this->updateOrthoMatrix(true);
}

void Viewer::updateOrthoMatrix(bool updateControl)
{
    if (updateControl) {
        if (this->ctrlwin != NULL) {
            if (this->rotateAngle < 45) {
                this->ctrlwin->updateAxesInfo(
                    this->pointCenter[0] - this->pointScale[0],
                    this->pointCenter[0] + this->pointScale[0],
                    this->pointCenter[1] - this->pointScale[1],
                    this->pointCenter[1] + this->pointScale[1]);
            } else {
                this->ctrlwin->updateAxesInfo(
                    - this->pointCenter[1] - this->pointScale[1],
                    - this->pointCenter[1] + this->pointScale[1],
                    this->pointCenter[0] - this->pointScale[0],
                    this->pointCenter[0] + this->pointScale[0]);
            }
        }
    }

    this->matrix.setToIdentity();
    this->matrix.rotate(this->rotateAngle, 0.0, 0.0, 1.0);
    this->matrix.ortho(
                this->pointCenter[0] - this->pointScale[0], this->pointCenter[0] + this->pointScale[0],
                this->pointCenter[1] - this->pointScale[1], this->pointCenter[1] + this->pointScale[1],
            -1.0, 1.0);
    this->invmatrix = this->matrix.inverted();

    if (this->isMirrored) {
        this->matrixRight.setToIdentity();
        this->matrixRight.rotate(this->rotateAngle, 0.0, 0.0, 1.0);
        this->matrixRight.ortho(
                    this->pointCenter[0] - this->pointScale[0], this->pointCenter[0] + this->pointScale[0],
                    this->pointCenter[1] - this->pointScale[1], this->pointCenter[1] + this->pointScale[1],
                -1.0, 1.0);
        this->matrixRight.scale(1.0, -1.0, 1.0);
        this->invmatrixRight = this->matrixRight.inverted();
    }
}

void Viewer::setBgColor(double r, double g, double b)
{
    this->bgColor[0] = (int) (r*255);
    this->bgColor[1] = (int) (g*255);
    this->bgColor[2] = (int) (b*255);
    this->repaint();
}

void Viewer::setLeftValue(int n)
{
    this->leftDisplayVar = n;
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     n);
    this->repaint();
}

void Viewer::setRightValue(int n)
{
    this->rightDisplayVar = n;
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     n);
    this->repaint();
}

void Viewer::setGammaValue(double gamma)
{
    this->g.setCmapStats(gamma, -1.0, -1.0);
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);

    this->repaint();
}

void Viewer::setCenterValue(double center)
{
    this->g.setCmapStats(-1.0, center, -1.0);
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);
    this->repaint();
}

void Viewer::setSlopeValue(double slope)
{
    this->g.setCmapStats(-1.0, -1.0, slope);
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);
    this->repaint();
}

void Viewer::setCmapMinmax(double min, double max)
{
    this->g.updateCmapBound(0 + this->leftDisplayVar*2, min);
    this->g.updateCmapBound(1 + this->leftDisplayVar*2, max);
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);
    this->repaint();
}

void Viewer::requestCboundReset()
{
    this->ctrlwin->setColorbarBounds(this->filedata.get_minmax(), this->leftDisplayVar);
}

void Viewer::cycleCmap()
{
    this->g.cycleCmap();
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);
    this->repaint();
}

int Viewer::get_nq()
{
    return this->filedata.get_nq();
}

void Viewer::rescale(int i)
{
    
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: rescale()...\n");
    #endif
    if ( !this->filedata.fileLoaded ) return;

    // Zoom to a size capable of displaying full data
    if (i == 1) {
        double *bounds = this->filedata.get_bounds();
        double xmin = bounds[0];
        double xmax = bounds[1];
        double ymin = bounds[2];
        double ymax = bounds[3];
        
        if (this->rotateAngle > 45) {

            xmin = -bounds[3];
            xmax = bounds[2];
            ymin = bounds[0];
            ymax = bounds[1];
        }

        double dx = 0.02 * (xmax - xmin);
        double dy = 0.02 * (ymax - ymin);

        this->setByAxes(xmin-dx, xmax+dx, ymin-dy, ymax+dy);
        this->squareScale();
    }
    
    #if VIEWER_VERBOSE
    fprintf(stderr, "Viewer: rescale() done.\n");
    #endif
}

void Viewer::showScale(bool value)
{
    this->drawScale = value;
}

void Viewer::setLogscale(bool log)
{
    if ( !this->filedata.fileLoaded ) return;
    this->g.setLogscale(log);
    this->g.setValue(true,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->leftDisplayVar);
    this->g.setValue(false,
                     this->filedata.get_nc(),
                     this->filedata.get_cells(),
                     this->rightDisplayVar);
    this->repaint();
}

void Viewer::showColorbar(bool value)
{
    this->drawColorbar = value;
    this->repaint();
}

void Viewer::set1DPlot(bool visible, bool isTheta)
{
    this->drawPlotLine = visible;
    if (visible) this->plotIsTheta = isTheta;

    if (isTheta) {
        double y = 1.0 * cos(this->plotLineData[2]);
        this->update1DPlot(1.0, y, true);
    } else {
        double y = this->plotLineData[2] * this->plotLineData[2];
        this->update1DPlot(0.0, y, true);
    }

    this->repaint();

}

void Viewer::replot1DPlot()
{
    if (this->plotIsTheta) {
        double y = 1.0 * cos(this->plotLineData[2]);
        this->update1DPlot(1.0, y, true);
    } else {
        double y = this->plotLineData[2] * this->plotLineData[2];
        this->update1DPlot(0.0, y, true);
    }
}

void Viewer::savePPM(const char *fname)
{
    int i,j;
    int dimx = this->width();
    int dimy = this->height();
    float *pixels = new float[3*dimx*dimy];

    glReadBuffer(GL_BACK);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_FLOAT, pixels);

    std::ofstream fp(fname);
    fp << "P3\n" << dimx << " " << dimy << std::endl << "255\n";
    for (i=dimy-1; i>=0; --i) {
        for (j=0; j<dimx; ++j) {
            int pixelPos = (i*dimx+j) * 3;
            int r = (int) (pixels[pixelPos + 0] * 255.0);
            int g = (int) (pixels[pixelPos + 1] * 255.0);
            int b = (int) (pixels[pixelPos + 2] * 255.0);
            fp << r << " " << g << " " << b << std::endl;
        }
    }
    fp.close();
}

void Viewer::pxToPoint(int x, int y, double *tx, double *ty)
{
    double ttx, tty;
    if (this->rotateAngle == 90) {
        ttx = (x - this->width()/2) * this->pointScale[1] * 2 / this->width();
        tty = - (y - this->height()/2) * this->pointScale[0] * 2 / this->height();
    } else {
        ttx = (x - this->width()/2) * this->pointScale[0] * 2 / this->width();
        tty = - (y - this->height()/2) * this->pointScale[1] * 2 / this->height();
    }

    *tx = cos(DEG_RAD * this->rotateAngle) * ttx + sin(DEG_RAD * this->rotateAngle) * tty;
    *ty = -sin(DEG_RAD * this->rotateAngle) * ttx + cos(DEG_RAD * this->rotateAngle) * tty;
}

// Takes the pixel and returns the exact location of the underlying data in data coordinates
void Viewer::screenToData(int sx, int sy, double *tx, double *ty)
{
    double ssx = (double)sx / (double)this->width();
    double ssy = (double)sy / (double)this->height();
    ssx = ssx * 2 - 1;
    ssy = - (ssy * 2 - 1);
    QVector3D vec = this->invmatrix * QVector3D(ssx, ssy, 0.0);
    *tx = vec.x();
    *ty = vec.y();
}

/*
 * glwidget.cpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#include <QtGui>
 #include "glwidget.hpp"
 #include "helper.hpp"

 GLWidget::GLWidget(Helper *helper, QWidget *parent)
     : QGLWidget(QGLFormat(QGL::SampleBuffers), parent), helper(helper)
 {
     elapsed = 0;
     setFixedSize(200, 200);
     setAutoFillBackground(false);
 }

 void GLWidget::animate()
 {
     elapsed = (elapsed + qobject_cast<QTimer*>(sender())->interval()) % 1000;
     repaint();
 }

 void GLWidget::paintEvent(QPaintEvent *event)
 {
     QPainter painter;
     painter.begin(this);
     painter.setRenderHint(QPainter::Antialiasing);
     helper->paint(&painter, event, elapsed);
     painter.end();
 }

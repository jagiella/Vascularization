/*
 * glwidget.hpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#ifndef GLWIDGET_HPP_
#define GLWIDGET_HPP_

 #include <QGLWidget>

 class Helper;
 class QPaintEvent;
 class QWidget;

 class GLWidget : public QGLWidget
 {
     Q_OBJECT

 public:
     GLWidget(Helper *helper, QWidget *parent);

 public slots:
     void animate();

 protected:
     void paintEvent(QPaintEvent *event);

 private:
     Helper *helper;
     int elapsed;
 };

#endif /* GLWIDGET_HPP_ */

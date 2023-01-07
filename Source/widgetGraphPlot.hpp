/*
 * widget.hpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#ifndef WIDGET_HPP_
#define WIDGET_HPP_

#include <QWidget>

 class Helper;
 class QPaintEvent;

 class Widget : public QWidget
 {
     Q_OBJECT

 public:
     Widget(Helper *helper, QWidget *parent);

 public slots:
     void animate();

 protected:
     void paintEvent(QPaintEvent *event);

 private:
     Helper *helper;
     int elapsed;
 };


#endif /* WIDGET_HPP_ */

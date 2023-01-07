/*
 * Gui.hpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#ifndef GUI_HPP_
#define GUI_HPP_

 #include <QWidget>
 #include <QThread>

 #include "helper.hpp"
#include "VesselGraph.hpp"

 class QLabel;
 class QWidget;

 class Window : public QWidget
 {
     Q_OBJECT

 public:
     VesselGraph *vg;
     Window();


 private:
     Helper helper;
 };


 class MyThread : public QThread
  {
  public:
     VesselGraph *vg;
     Window *window;

      void run();
  };


#endif /* GUI_HPP_ */

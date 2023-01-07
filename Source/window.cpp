/*
 * Gui.cpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#include <QtGui>
#include <stdio.h>
 #include "glwidget.hpp"
 #include "widget.hpp"
 #include "window.hpp"

 Window::Window()
     : QWidget()
 {
	 this->vg = 0;

     Widget *native = new Widget(&helper, this);
     GLWidget *openGL = new GLWidget(&helper, this);
     QLabel *nativeLabel = new QLabel(tr("Native"));
     nativeLabel->setAlignment(Qt::AlignHCenter);
     QLabel *openGLLabel = new QLabel(tr("OpenGL"));
     openGLLabel->setAlignment(Qt::AlignHCenter);

     QGridLayout *layout = new QGridLayout;
     layout->addWidget(native, 0, 0);
     layout->addWidget(openGL, 0, 1);
     layout->addWidget(nativeLabel, 1, 0);
     layout->addWidget(openGLLabel, 1, 1);
     setLayout(layout);

     QTimer *timer = new QTimer(this);
     connect(timer, SIGNAL(timeout()), native, SLOT(animate()));
     connect(timer, SIGNAL(timeout()), openGL, SLOT(animate()));
     timer->start(100);

     setWindowTitle(tr("2D Painting on Native and OpenGL Widgets (test)"));
 }


 /*void MyThread::run()
  {
     fprintf(stderr, "Show Window\n");
     this->window->show();

     fprintf(stderr, "Execute Thread\n");
      exec();

      fprintf(stderr, "Return to main Thread\n");
      return;
  }*/


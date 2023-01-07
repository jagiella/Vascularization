/*
 * widget.cpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#include <QtGui>
 #include "widget.hpp"
 #include "helper.hpp"
 #include "window.hpp"

 Widget::Widget(Helper *helper, QWidget *parent)
     : QWidget(parent), helper(helper)
 {
     elapsed = 0;
     setFixedSize(200, 200);
 }

 void Widget::animate()
 {
     elapsed = (elapsed + qobject_cast<QTimer*>(sender())->interval()) % 1000;
     repaint();
 }

 void Widget::paintEvent(QPaintEvent *event)
 {
	QPainter painter;
	painter.begin(this);
	painter.setRenderHint(QPainter::Antialiasing);
	//helper->paint(&painter, event, elapsed);

	//if (false)
		if (((Window*) this->parent())->vg != 0) {
			// fprintf(stderr, "REPAINT\n");
			// PAINT VESSEL GRAPH
			QLinearGradient gradient(QPointF(50, -20), QPointF(80, 20));
			int r = 100;
			int g = 0;
			int b = 0;
			int a = 255;
			QColor *color = new QColor(r, g, b, a);

			float max = 0;
			for (int i = 0;
					i < ((Window*) this->parent())->vg->countVesselNodes; i++) {
				if (
						max
						< ((Window*) this->parent())->vg->vesselNodes[i]->marker)
					max = ((Window*) this->parent())->vg->vesselNodes[i]->marker;
			}
			//fprintf(stderr, "[MAX=%lf]\n", max);
			int dx = (200 / ((Window*) this->parent())->vg->countVesselNodes);
			for (int i = 0;
					i < ((Window*) this->parent())->vg->countVesselNodes; i++) {
				//painter.setBrush(QBrush(Qt::blue));//gradient));
				r = g = b = (int) (
						((Window*) this->parent())->vg->vesselNodes[i]->marker
						* 255. / max);
				color->setRgb(r, g, b, a);
				painter.setBrush(QBrush(color->rgb()));//gradient));
				painter.setPen(color->rgb());

				int
						x = (int) (((Window*) this->parent())->vg->vesselNodes[i]->position[0]);
				int
						y = (int) (((Window*) this->parent())->vg->vesselNodes[i]->position[1]);
				painter.drawRect(x * dx, 0, dx, 60);
			}
		}

	painter.end();
 }

/*
 * widget.cpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

#include <QtGui>
#include <QPoint>
 #include "widget.hpp"
 #include "helper.hpp"
 #include "window.hpp"

//#include <math.h>
//#define RAND01 (rand()/RAND_MAX)

 Widget::Widget(Helper *helper, QWidget *parent)
     : QWidget(parent), helper(helper)
 {
     elapsed = 0;
     setFixedSize(200, 200);

     cellCount = 0;
     time = 0;
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

	if (true)
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
		for (int i = 0; i < ((Window*) this->parent())->vg->countVesselNodes; i++) {
			if (max < ((Window*) this->parent())->vg->vesselNodes[i]->marker)
				max = ((Window*) this->parent())->vg->vesselNodes[i]->marker;
		}
		//fprintf(stderr, "[MAX=%lf]\n", max);
		int dx = (200 / ((Window*) this->parent())->vg->countVesselNodes);
		for (int i = 0; i < ((Window*) this->parent())->vg->countVesselNodes; i++) {
			//painter.setBrush(QBrush(Qt::blue));//gradient));
			r
					= g
							= b
									= (int) (((Window*) this->parent())->vg->vesselNodes[i]->marker
											* 255. / max);
			color->setRgb(r, g, b, a);
			painter.setBrush(QBrush(color->rgb()));//gradient));
			painter.setPen(color->rgb());

			int
					x =

							(int) (((Window*) this->parent())->vg->vesselNodes[i]->position[0]);
			int
					y =
							(int) (((Window*) this->parent())->vg->vesselNodes[i]->position[1]);
			painter.drawRect(x * dx, 0, dx, 60);
		}


		// Paint Curve

		int x_min =  0, x_max = this->width();
		int y_min = this->width(), y_max = this->width()-100;

		painter.setPen(Qt::white);
		painter.setBrush(QBrush(Qt::white));
		painter.drawRect(x_min, y_min, x_max-x_min, y_max-y_min);
		/*for( int x=x_min; x<x_max; x++)
			for( int y=y_min; y<y_max; y++)
			{

			}*/
		/*for (int i = 0; i < ((Window*) this->parent())->vg->countVesselNodes; i++) {
			painter.setPen(Qt::red);

			int x = x_min + (int) (((Window*) this->parent())->vg->vesselNodes[i]->position[0]);
			int y = (int) (((Window*) this->parent())->vg->vesselNodes[i]->marker * 100. / max);
			painter.drawPoint( x,y_min - y);
		}*/
		for (int i = 0; i < ((Window*) this->parent())->vg->countVesselSegments; i++) {
			painter.setPen(Qt::red);

			int x0 = x_min + (int) (((Window*) this->parent())->vg->vesselSegments[i]->vesselNodes[0]->position[0]);
			int y0 = y_min + (int) ((y_max-y_min) * ((Window*) this->parent())->vg->vesselSegments[i]->vesselNodes[0]->marker / max);
			int x1 = x_min + (int) (((Window*) this->parent())->vg->vesselSegments[i]->vesselNodes[1]->position[0]);
			int y1 = y_min + (int) ((y_max-y_min) * ((Window*) this->parent())->vg->vesselSegments[i]->vesselNodes[1]->marker / max);
			painter.drawLine( x0,y0, x1,y1);
		}
	}

	//
	if(false){
		qreal r = 1;//elapsed/1000.0;

		// init
		if(cellCount == 0){
			cellCount = 1;
			this->cellPositions = (double**) malloc( cellCount * sizeof(double *));
			this->cellVelocity = (double**) malloc( cellCount * sizeof(double *));
			this->cellAcceleration = (double**) malloc( cellCount * sizeof(double *));
			this->cellRadius = (double*) malloc( cellCount * sizeof(double));
			for( int i=0; i<cellCount; i++){
				this->cellPositions[i] = (double*) malloc( 2 * sizeof(double));
				this->cellVelocity[i] = (double*) malloc( 2 * sizeof(double));
				this->cellAcceleration[i] = (double*) malloc( 2 * sizeof(double));

				this->cellPositions[i][0]=i*3;
				this->cellPositions[i][1]=i*3;
				this->cellVelocity[i][0]=3*i;
				this->cellVelocity[i][1]=3*i;

				cellRadius[i] = 5;
			}
		}

		if(cellCount > 0 && ((int)time)%100 == 0){
			int old_cellCount = cellCount;
			cellCount*=2;;
			//cellCount++;
			this->cellPositions = (double**) realloc( cellPositions, cellCount * sizeof(double *));
			this->cellVelocity = (double**) realloc( cellVelocity, cellCount * sizeof(double *));
			this->cellAcceleration = (double**) realloc( cellAcceleration, cellCount * sizeof(double *));
			this->cellRadius = (double*) realloc( cellRadius, cellCount * sizeof(double));
			for( int i=old_cellCount; i<cellCount; i++){
				this->cellPositions[i] = (double*) malloc( 2 * sizeof(double));
				this->cellVelocity[i] = (double*) malloc( 2 * sizeof(double));
				this->cellAcceleration[i] = (double*) malloc( 2 * sizeof(double));

				this->cellPositions[i][0]=cellPositions[i-old_cellCount][0]+0.1;
				this->cellPositions[i][1]=cellPositions[i-old_cellCount][1]+0.1;
				this->cellVelocity[i][0] =cellVelocity[i-old_cellCount][0];
				this->cellVelocity[i][1] =cellVelocity[i-old_cellCount][1];

				cellRadius[i] = 5;
				cellRadius[i-old_cellCount] = 5;
			}
		}

		// update acceleration
		for( int i=0; i<cellCount; i++){
			cellAcceleration[i][0] = 0;
			cellAcceleration[i][1] = 0;
			for( int ii=0; ii<cellCount; ii++)
				if( i!=ii &&
					(cellPositions[i][0] - cellPositions[ii][0])*(cellPositions[i][0] - cellPositions[ii][0]) +
					(cellPositions[i][1] - cellPositions[ii][1])*(cellPositions[i][1] - cellPositions[ii][1])
					< (cellRadius[i]+cellRadius[ii])*(cellRadius[i]+cellRadius[ii])
					){

					// averaging
					//cellAcceleration[i][0] += (cellVelocity[ii][0] - cellVelocity[i][0])/2.;
					//cellAcceleration[i][1] += (cellVelocity[ii][1] - cellVelocity[i][1])/2.;

					// repulsion
					double distance = sqrt((cellPositions[i][0] - cellPositions[ii][0])*(cellPositions[i][0] - cellPositions[ii][0]) +
							(cellPositions[i][1] - cellPositions[ii][1])*(cellPositions[i][1] - cellPositions[ii][1]));
					if( distance*3./4. < cellRadius[i]+cellRadius[ii]){
						//cellAcceleration[i][0] += (cellPositions[i][0] - cellPositions[ii][0])*0.1;
						//cellAcceleration[i][1] += (cellPositions[i][1] - cellPositions[ii][1])*0.1;
						cellAcceleration[i][0] += (cellPositions[i][0] - cellPositions[ii][0])*((cellRadius[i]+cellRadius[ii])/distance - 1.)*100;
						cellAcceleration[i][1] += (cellPositions[i][1] - cellPositions[ii][1])*((cellRadius[i]+cellRadius[ii])/distance - 1.)*100;

					}

				}

			// attraction
			QPoint point = mapFromGlobal(QCursor::pos());
			cellAcceleration[i][0] += (point.x() - cellPositions[i][0]);
			cellAcceleration[i][1] += (point.y() - cellPositions[i][1]);

			// friction
			//cellAcceleration[i][0] += -1*cellVelocity[i][0];
			//cellAcceleration[i][1] += -1*cellVelocity[i][1];

			// random movement
			cellAcceleration[i][0] += 10.*(RAND01-0.5);
			cellAcceleration[i][1] += 10.*(RAND01-0.5);
		}

		// update velocity
		for( int i=0; i<cellCount; i++){
			cellVelocity[i][0] += cellAcceleration[i][0];
			cellVelocity[i][1] += cellAcceleration[i][1];
			//fprintf(stderr, "acceleration (%lf,%lf)\n", cellAcceleration[i][0],cellAcceleration[i][1]);

			// absolut velocity
			double velocity = sqrt(cellVelocity[i][0]*cellVelocity[i][0] + cellVelocity[i][1]*cellVelocity[i][1]);

			if( velocity > 10.)
			{
				cellVelocity[i][0] *= 10./ velocity;
				cellVelocity[i][1] *= 10./ velocity;
			}
			//fprintf(stderr, "velocity (%lf,%lf) -> abs %lf\n", cellVelocity[i][0],cellVelocity[i][1], velocity);


			// max velocity
			/*if( cellVelocity[i][0] > 10)
				cellVelocity[i][0] = 10;
			if( cellVelocity[i][1] > 10)
				cellVelocity[i][1] = 10;
			// min velocity
			if( cellVelocity[i][0] < -10)
				cellVelocity[i][0] = -10;
			if( cellVelocity[i][1] < -10)
				cellVelocity[i][1] = -10;*/
		}


		// update position
		for( int i=0; i<cellCount; i++){
			this->cellPositions[i][0] += this->cellVelocity[i][0]*0.1;
			this->cellPositions[i][1] += this->cellVelocity[i][1]*0.1;
			/*if(this->cellPositions[i][0]>200)
				this->cellPositions[i][0]-=200;
			if(this->cellPositions[i][1]>200)
				this->cellPositions[i][1]-=200;*/
		}
		time++;

		// update radius
		for( int i=0; i<cellCount; i++){
			cellRadius[i] = 5. + 2. * (((int)time)%100) / 100.;
		}


		// print
		for( int i=0; i<cellCount; i++){
			float circleRadius = 5;

			//painter.setPen(Qt::white);
			painter.setBrush(QBrush(Qt::white));

			painter.drawEllipse(QRectF(this->cellPositions[i][0]-cellRadius[i], this->cellPositions[i][1]-cellRadius[i],
					cellRadius[i]*2, cellRadius[i]*2));
		}
	}

	painter.end();
 }

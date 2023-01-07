/*
 * Perfusion.cpp
 *
 *  Created on: 08.03.2011
 *      Author: jagiella
 */

/*
 * main.cpp
 *
 *  Created on: Apr 12, 2010
 *      Author: jagiella
 */

//#define _POSIX_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <time.h>


#include "LinearAlgebra.hpp"
#include "IO.hpp"

#include "VesselGraph.hpp"
#include "VoronoiDiagram.hpp"
#include "Tumor.hpp"

time_t mytime ( time_t * timer ) { return time(timer);};

// MACROS
inline int cprintf( FILE * stream, int escape, const char * format, ...){

	va_list args;

	fprintf( stream, "\x1b[%im", escape);
	va_start(args, format);
	vfprintf( stream, format, args );
	va_end(args);
	fprintf( stream, "\x1b[0m");
}

//#include "Triangulation.hpp"

//#include "Octree.hpp"


/*color_t rgbformulaeMapping(char rformulae, char gformulae, char bformulae,
 double x, double min, double max) {
 color_t color = 0;
 char *formulae[3];
 formulae[0] = &rformulae;
 formulae[1] = &gformulae;
 formulae[2] = &bformulae;

 for (int c = 0; c < 3; c++) {
 //fprintf(stderr, "%lf\n", rgbformulae( *formulae[c], (x-min) / (max-min)));
 ((channel_t*) (&color))[c] = (channel_t) 255. * rgbformulae(
 *formulae[c], (x - min) / (max - min));
 }

 return color;
 }*/









/*void setInitialVesselGraphOLD(VesselGraph *vesselGraph) {
	VesselNode *last;
	VesselNode *next = new VesselNode(0, 0, 0);
	vesselGraph->addVesselNode(next);
	next->pressure = 12.; // kPa
	next->setType(ROOT);
	int position[3] = { 0, 0, 0 };
	int max[3] = { DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z };
	for (int i = 0; i < DOMAIN_SIZE_X + DOMAIN_SIZE_Y + DOMAIN_SIZE_Z - 2;
			i++) {

		// chose direction
		//int dx[3] = {0,0,0};
		//dx[(int)(RAND01*3.)] ++;
		int dxi;
		do {
			dxi = (int) (RAND01 * (float) DIMENSIONS);
			fprintf(stderr, "i: %i, dxi: %i, position[dxi]: %i\n", i, dxi,
					position[dxi]);
		} while (position[dxi] == max[dxi] - 1);

		position[dxi]++;
		//fprintf( stderr, "(%i, %i, %i)\n", position[0], position[1], position[2]);

		// new vessel node
		last = next;
		next = new VesselNode(position[0], position[1], position[2]);
		next->setType(VESSEL);
		vesselGraph->addVesselNode(next);
		vesselGraph->addVesselSegment(new VesselSegment(last, next));
		vesselGraph->vesselSegments[vesselGraph->countVesselSegments - 1]->radius =
				25;

		position[dxi]++;
		last = next;
		next = new VesselNode(position[0], position[1], position[2]);
		next->setType(VESSEL);
		vesselGraph->addVesselNode(next);
		vesselGraph->addVesselSegment(new VesselSegment(last, next));
		vesselGraph->vesselSegments[vesselGraph->countVesselSegments - 1]->radius =
				25;
		i++;

		//fprintf( stdout, "%i %i %i %i %i %i\n",
		//		last->position[0], last->position[1], last->position[2],
		//		next->position[0]-last->position[0], next->position[1]-last->position[1], next->position[2]-last->position[2]);
	}
	next->pressure = 12;
	next->setType(ROOT);

	next = new VesselNode(0, 100, 0);
	vesselGraph->addVesselNode(next);
	next->pressure = 0.; // kPa
	next->setType(ROOT);
	position[0] = 0;
	position[1] = 100;
	position[2] = 0;
	for (int i = 0; i < DOMAIN_SIZE_X + DOMAIN_SIZE_Y + DOMAIN_SIZE_Z - 2;
			i++) {

		// chose direction
		//int dx[3] = {0,0,0};
		//dx[(int)(RAND01*3.)] ++;
		int dxi;
		do {
			dxi = (int) (RAND01 * (float) DIMENSIONS);
			fprintf(stderr, "i: %i, dxi: %i, position[dxi]: %i\n", i, dxi,
					position[dxi]);
		} while ((dxi == 0 && position[dxi] == max[dxi] - 1)
				|| (dxi == 1 && position[dxi] == 0));

		position[dxi] += 1 - 2 * dxi;
		//fprintf( stderr, "(%i, %i, %i)\n", position[0], position[1], position[2]);

		// new vessel node
		last = next;
		next = vesselGraph->octree->at(position[0], position[1], position[2]);
		if (next == 0) {
			next = new VesselNode(position[0], position[1], position[2]);
			next->setType(VESSEL);
			vesselGraph->addVesselNode(next);
		}
		vesselGraph->addVesselSegment(new VesselSegment(last, next));
		vesselGraph->vesselSegments[vesselGraph->countVesselSegments - 1]->radius =
				25;

		position[dxi] += 1 - 2 * dxi;
		last = next;
		next = vesselGraph->octree->at(position[0], position[1], position[2]);
		if (next == 0) {
			next = new VesselNode(position[0], position[1], position[2]);
			next->setType(VESSEL);
			vesselGraph->addVesselNode(next);
		}
		vesselGraph->addVesselSegment(new VesselSegment(last, next));
		vesselGraph->vesselSegments[vesselGraph->countVesselSegments - 1]->radius =
				25;
		i++;

		//fprintf( stdout, "%i %i %i %i %i %i\n",
		//		last->position[0], last->position[1], last->position[2],
		//		next->position[0]-last->position[0], next->position[1]-last->position[1], next->position[2]-last->position[2]);
	}
	next->pressure = 0;
	next->setType(ROOT);
}*/

/*void updateOxygen(float **oxygen, float **oxygen2, VesselGraph *vesselGraph,
		int maxit, int MaxError) {
	float **temp;
	float maxError = 0;
	float lastmaxError = 0;
	int it = 0;
	do {
		lastmaxError = maxError;
		maxError = 0;
		//for( int it=0; it<1000; it++){

		//for( int i=0; i<vesselGraph->countVesselNodes; i++)
		//	oxygen[vesselGraph->vesselNodes[i]->position[0]][vesselGraph->vesselNodes[i]->position[1]] = 1.;

		// CALCULATE OXYGEN CONCENTRATION
		//float D = 1.;
		float dt = 10.;
		float dx = 10.;
		for (int x = 0; x < DOMAIN_SIZE_X; x++) {
			for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
				//fprintf( stderr, "%i %i\n", x,y);
				float diffusion = 0.;
				float consumption = 1. / 10000.; // consumption rate normal tissue
				//- 0.004*(1.-oxygen[x][y]); // source rate vessels
				for (int xx = x - 1; xx <= x + 1; xx += 2)
					for (int yy = y - 1; yy <= y + 1; yy += 2)
						if (xx >= 0 && xx < DOMAIN_SIZE_X && yy >= 0 && yy
						< DOMAIN_SIZE_Y) {
						//fprintf( stderr, "%i %i %i %i\n", x,y,xx,yy);
							diffusion += oxygen[xx][yy] - oxygen[x][y];
						}
				diffusion /= (dx * dx);
				oxygen2[x][y] = oxygen[x][y] + (diffusion - consumption) * dt;
				if (maxError < fabs(oxygen2[x][y] - oxygen[x][y]))
					maxError = fabs(oxygen2[x][y] - oxygen[x][y]);
			}
		}

		for (int i = 0; i < vesselGraph->countVesselNodes; i++)
			if (vesselGraph->vesselNodes[i]->getType() == VESSEL
					&& vesselGraph->vesselNodes[i]->countNeighbors > 0
					&& fabs(
							vesselGraph->vesselNodes[i]->pressure
									- vesselGraph->vesselNodes[i]->neighbors[0]->pressure)
							> 1e-4)
				oxygen2[(int) vesselGraph->vesselNodes[i]->position[0]][(int) vesselGraph->vesselNodes[i]->position[1]] +=
						0.004
								* (1.
										- oxygen[(int) vesselGraph->vesselNodes[i]->position[0]][(int) vesselGraph->vesselNodes[i]->position[1]])
								* dt; // source rate vessels ;

		for (int x = 0; x < DOMAIN_SIZE_X; x++)
			for (int y = 0; y < DOMAIN_SIZE_Y; y++)
				if (oxygen2[x][y] < 0.)
					oxygen2[x][y] = 0.;

		temp = oxygen;
		oxygen = oxygen2;
		oxygen2 = temp;

		//for( int i=0; i<vesselGraph->countVesselNodes; i++)
		//	oxygen[vesselGraph->vesselNodes[i]->position[0]][vesselGraph->vesselNodes[i]->position[1]] = 1.;
		it++;
		fprintf(stderr, "\rmax error = %lf \b", maxError);
	} while ((maxError > MaxError && fabs(maxError - lastmaxError) > 1e-20)
			&& it < maxit);
	fprintf(stderr, "INFO: Oxygen Diffusion finished after %i iterations\n",
			it);

}

void updateVEGF(float **oxygen, float **vegf) {
	for (int x = 0; x < DOMAIN_SIZE_X; x++)
		for (int y = 0; y < DOMAIN_SIZE_Y; y++)
			vegf[x][y] = 0.;

	int dx = 10;
	int Rvegf = 200 / dx;
	float fRvegf = (float) Rvegf;
	//int dist = Rvegf / dx;
	for (int x = 0; x < DOMAIN_SIZE_X; x++)
		for (int y = 0; y < DOMAIN_SIZE_Y; y++)
			if (oxygen[x][y] < 0.1)
				for (int xx = x - Rvegf; xx <= x + Rvegf; xx++)
					for (int yy = y - Rvegf; yy <= y + Rvegf; yy++)
						if (xx >= 0 && xx < DOMAIN_SIZE_X && yy >= 0 && yy
						< DOMAIN_SIZE_Y) {
							float fdist = sqrt(pow(xx - x, 2) + pow(yy - y, 2));
							if (fdist < fRvegf)
								vegf[xx][yy] += 1 - fdist / fRvegf;
						}

}



void updateVesselGraph(float **oxygen, float **vegf, VesselGraph *vesselGraph,
		float dt, float time) {
	int countPreviousVesselNodes = vesselGraph->countVesselNodes;
	int countPreviousVesselSegments = vesselGraph->countVesselSegments;

	for (int i = 0; i < countPreviousVesselNodes; i++) {

		if (time - vesselGraph->vesselNodes[i]->time >= 100)
			switch (vesselGraph->vesselNodes[i]->getType()) {
			case SPROUT:
			case TIP:
				vesselGraph->vesselNodes[i]->setType(VESSEL);
				break;
			}

		switch (vesselGraph->vesselNodes[i]->getType()) {
		case ROOT: //break;

		case VESSEL:
			// sprout initiation
		case SPROUT:
			// sprout migration
			if ((vesselGraph->vesselNodes[i]->countNeighbors == 0
					|| (vesselGraph->vesselNodes[i]->countNeighbors < 3
							&& vesselGraph->vesselNodes[i]->neighbors[0]->countNeighbors
									< 3
							&& (vesselGraph->vesselNodes[i]->countNeighbors == 1
									|| vesselGraph->vesselNodes[i]->neighbors[1]->countNeighbors
											< 3)
							&& vegf[(int) vesselGraph->vesselNodes[i]->position[0]][(int) vesselGraph->vesselNodes[i]->position[1]]
									> 0.)) && dt / 5. >= RAND01) {
				int candidateX = -1;
				int candidateY = -1;
				int x = vesselGraph->vesselNodes[i]->position[0];
				int y = vesselGraph->vesselNodes[i]->position[1];
				float maxVEGF = vegf[x][y];

				for (int xx = x - 1; xx <= x + 1; xx += 2)
					if (xx >= 0 && xx < DOMAIN_SIZE_X /*&& yy>=0 && yy<DOMAIN_SIZE* /
					&& vegf[xx][y] > maxVEGF
							&& (y == 0
									|| vesselGraph->octree->at(xx, y - 1, 0)
											== 0)
							&& (vesselGraph->octree->at(xx, y, 0) == 0)
							&& (vesselGraph->octree->at(xx, y + 1, 0) == 0)) {
						maxVEGF = vegf[xx][y];
						candidateX = xx;
						candidateY = y;
					}

				for (int yy = y - 1; yy <= y + 1; yy += 2)
					if ( /*xx>=0 && xx<DOMAIN_SIZE &&* /yy >= 0
							&& yy < DOMAIN_SIZE_Y && vegf[x][yy] > maxVEGF
							&& (x == 0
									|| vesselGraph->octree->at(x - 1, yy, 0)
											== 0)
							&& vesselGraph->octree->at(x, yy, 0) == 0
							&& vesselGraph->octree->at(x + 1, yy, 0) == 0) {
						maxVEGF = vegf[x][yy];
						candidateX = x;
						candidateY = yy;
					}
				//fprintf( stderr, "Sprout! maxVEGF=%e\n", maxVEGF);

				if (candidateX != -1) {
					VesselNode *next = new VesselNode(candidateX, candidateY,
							0);
					next->setType(TIP);
					vesselGraph->addVesselNode(next);
					vesselGraph->addVesselSegment(
							new VesselSegment(vesselGraph->vesselNodes[i],
									next));
					if (vesselGraph->vesselNodes[i]->getType() == SPROUT)
						next->time = vesselGraph->vesselNodes[i]->time;
					else
						next->time = time;
				}

			}
			break;

		case TIP:

			// new position
			int pos[3];
			for (int d = 0; d < 2; d++) {
				pos[d] =
						2 * vesselGraph->vesselNodes[i]->position[d]
								- vesselGraph->vesselNodes[i]->neighbors[0]->position[d];
			}

			if (pos[0] >= 0 && pos[0] < DOMAIN_SIZE_X && pos[1] >= 0
					&& pos[1] < DOMAIN_SIZE_Z && dt / 5. >= RAND01) {
				VesselNode *next = vesselGraph->octree->at(pos[0], pos[1], 0);
				vesselGraph->vesselNodes[i]->setType(SPROUT);
				if (next == 0) {
					//fprintf( stderr, "Sprout Migr.!\n");

					next = new VesselNode(pos[0], pos[1], 0);
					next->setType(TIP);
					next->time = vesselGraph->vesselNodes[i]->time;
					vesselGraph->addVesselNode(next);
					//if( time - vesselGraph->vesselNodes[i]->time >= 100)
					//next->setType( VESSEL);
				} //else fprintf( stderr, "Sprout meets Network !\n");
				vesselGraph->addVesselSegment(
						new VesselSegment(vesselGraph->vesselNodes[i], next));
			}

			// termination
			/*if( time - vesselGraph->vesselNodes[i]->time >= 100)
			 {
			 VesselNode *next = vesselGraph->vesselNodes[i];
			 while( next->time == vesselGraph->vesselNodes[i]->time){
			 next->setType( VESSEL);
			 if( next->neighbors[0]->type == SPROUT && next->neighbors[0]->time == vesselGraph->vesselNodes[i]->time)
			 next = next->neighbors[0];
			 else
			 next = next->neighbors[1];
			 }
			 }* /

			break;

		default:
			break;
		}
		//if( vegf[vesselGraph->vesselNodes[i]->position[0]][vesselGraph->vesselNodes[i]->position[1]] > 0.
		//		&& 1./50. >= RAND01){
	}

	for (int i = 0; i < countPreviousVesselSegments; i++) {
		if ((vesselGraph->vesselSegments[i]->vesselNodes[0]->getType() == VESSEL
				|| vesselGraph->vesselSegments[i]->vesselNodes[0]->getType()
						== ROOT)
				&& (vesselGraph->vesselSegments[i]->vesselNodes[1]->getType()
						== VESSEL
						|| vesselGraph->vesselSegments[i]->vesselNodes[1]->getType()
								== ROOT)) {
			if (time - vesselGraph->vesselSegments[i]->vesselNodes[0]->time
					>= 124
					&& vegf[(int) vesselGraph->vesselSegments[i]->vesselNodes[0]->position[0]][(int) vesselGraph->vesselSegments[i]->vesselNodes[0]->position[1]]
							> 0 && vesselGraph->vesselSegments[i]->radius < 25)
				vesselGraph->vesselSegments[i]->radius += 0.4 * dt;

			if (vesselGraph->vesselSegments[i]->shear < 0.002) {
				if (dt / 20. >= RAND01 /*&& vesselGraph->vesselSegments[i]->radius < 20* /) {
					fprintf(stderr, "Vessel Regression!\n");
					VesselSegment *removedVS = vesselGraph->vesselSegments[i];
					VesselNode *removedVN0 = removedVS->vesselNodes[0];
					VesselNode *removedVN1 = removedVS->vesselNodes[1];

					//removedVS->vesselNodes[0]->removeNeighbor( removedVS->vesselNodes[1]);
					//removedVS->vesselNodes[1]->removeNeighbor( removedVS->vesselNodes[0]);

					vesselGraph->removeVesselSegment(removedVS);
					delete removedVS;

					if (removedVN0->countNeighbors
							== 0 && removedVN0->getType() != ROOT) {VesselNode *removedVN = removedVN0;
					vesselGraph->removeVesselNode(removedVN);
					delete removedVN;
				}
					if (removedVN1->countNeighbors
							== 0 && removedVN1->getType() != ROOT) {VesselNode *removedVN = removedVN1;
					vesselGraph->removeVesselNode(removedVN);
					delete removedVN;
				}
			}
		}
		//else

	}
}
}*/
VoronoiDiagram *createVoronoiTesselation( VesselGraph *vesselGraph, float latticeConstant = 1.){
	VoronoiDiagram *t = new VoronoiDiagram();

	t->LATTICE_CONSTANT = latticeConstant;

	srand(0);
	float perturbation2 = 1e-3;

	// ADD VESSEL NODES
	for( int i=0; i<vesselGraph->countVesselNodes;i++)
		//t->add( new VoronoiCell( vesselGraph->vesselNodes[i]));
		t->add( new VoronoiCell(vesselGraph->vesselNodes[i]->position[0]*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+perturbation2*(RAND01-0.5),
								vesselGraph->vesselNodes[i]->position[1]*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+perturbation2*(RAND01-0.5),
								vesselGraph->vesselNodes[i]->position[2]*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+perturbation2*(RAND01-0.5),
								vesselGraph->vesselNodes[i]));
	fprintf( stderr, "Added %i Vessel Nodes to Triangulation\n", vesselGraph->countVesselNodes);

	/*for( int xi=0; xi<vesselGraph->DOMAIN_SIZE_X*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT;xi++)
		for( int yi=0; yi<vesselGraph->DOMAIN_SIZE_Y*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT; yi++)
			for( int zi=0; zi<vesselGraph->DOMAIN_SIZE_Z*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT; zi++)
				if( vesselGraph->octree->at( xi, yi,zi) == 0){
					//fprintf( stderr,"found some space\n");
					t->add( new VoronoiCell(xi+0.5+perturbation2*(RAND01-0.5),
											yi+0.5+perturbation2*(RAND01-0.5),
											zi+0.5+perturbation2*(RAND01-0.5)));
				}
	*/
	for( int xi=0; xi<vesselGraph->DOMAIN_SIZE_X;xi += t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT)
		for( int yi=0; yi<vesselGraph->DOMAIN_SIZE_Y; yi += t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT)
			for( int zi=0; zi<vesselGraph->DOMAIN_SIZE_Z; zi += t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT){
				if( vesselGraph->octree->firstBetween( xi, yi, zi, xi+t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT, yi+t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT, zi+t->LATTICE_CONSTANT/vesselGraph->LATTICE_CONSTANT) == 0){
					t->add( new VoronoiCell(xi*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+0.5+perturbation2*(RAND01-0.5),
											yi*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+0.5+perturbation2*(RAND01-0.5),
											zi*vesselGraph->LATTICE_CONSTANT/t->LATTICE_CONSTANT+0.5+perturbation2*(RAND01-0.5)));
				}
			}
	fprintf( stderr, "Added %i Spacefilling Points to Triangulation\n", t->numberOfVertices() - vesselGraph->countVesselNodes);

	t->setFramePoints();
	t->triangulate();

	return t;
}

void readBinaryAndExportToPovray( char *binary_filename, VesselGraph *&vesselGraph, CONCENTRATION_T ***&marker, VoronoiDiagram *&t)
{
	//char povray_filename[512];
	if( strstr( binary_filename, ".Ci.binary") && !strstr( binary_filename, ".Ci.binary.")){
		fprintf( stderr, "Read interstitial data from binary file (%s)\n", binary_filename);
		FILE *fpi = fopen( binary_filename, "rb");


		//fprintf( stderr, "-> Dimension: %ix%ix%i\n", DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);
		float markerp;
		switch( vesselGraph->type()){

			case VesselGraph::RegularLattice:{

				// Read Binary Header
				int DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z;
				fread( &DOMAIN_SIZE_X,  sizeof(int), 1, fpi);
				fread( &DOMAIN_SIZE_Y,  sizeof(int), 1, fpi);
				fread( &DOMAIN_SIZE_Z,  sizeof(int), 1, fpi);

				// INIT if necessary
				if( marker==0){
					marker = (CONCENTRATION_T***) malloc(sizeof(CONCENTRATION_T**)*DOMAIN_SIZE_X);
					for( int x=0; x<DOMAIN_SIZE_X;x++){
						marker[x] = (CONCENTRATION_T**) malloc(sizeof(CONCENTRATION_T*)*DOMAIN_SIZE_Y);
						for( int y=0; y<DOMAIN_SIZE_Y;y++){
							marker[x][y] = (CONCENTRATION_T*) malloc(sizeof(CONCENTRATION_T)*DOMAIN_SIZE_Z);
						}
					}
				}

				// READ DATA
				for( int x=0; x<DOMAIN_SIZE_X; x++)
					for( int y=0; y<DOMAIN_SIZE_Y; y++)
						for( int z=0; z<DOMAIN_SIZE_Z; z++){
							fread( &markerp,  sizeof(float), 1, fpi);
							marker[x][y][z] = markerp;
						}

				// OUTPUT to POVRAY
				char povray_filename[512];
				sprintf( povray_filename, "%s.pov", binary_filename);
				vesselGraph->printToPovray(povray_filename, marker);
				sprintf(povray_filename, "nice povray -D0 +h480 +w640 +ua %s.pov 2> err.out &", binary_filename);
				system(povray_filename);

			}break;

			case VesselGraph::IrregularLattice:{

				// INIT if necessary
				if( t==0){
					// CELL LATTICE CONSTANT
					/*float latticeConstant = 1;

					// CREATE TRIANGULATION
					t = new VoronoiDiagram();

					// ADD VESSEL NODES
					fprintf( stderr, "Add %i Vessel Nodes to Triangulation\n", vesselGraph->countVesselNodes);
					float perturbation1 = 1e-3;
					for( int i=0; i<vesselGraph->countVesselNodes;i++)
						//t->add( new VoronoiCell( vesselGraph->vesselNodes[i]));
						t->add( new VoronoiCell(
										vesselGraph->vesselNodes[i]->position[0]/latticeConstant+perturbation1*RAND01,
										vesselGraph->vesselNodes[i]->position[1]/latticeConstant+perturbation1*RAND01,
										vesselGraph->vesselNodes[i]->position[2]/latticeConstant+perturbation1*RAND01));
					// ADD
					//fprintf( stderr, "Add %i Spacefilling Points to Triangulation\n",
					//		(int)(vesselGraph->DOMAIN_SIZE_X/latticeConstant*vesselGraph->DOMAIN_SIZE_Y/latticeConstant*vesselGraph->DOMAIN_SIZE_Z/latticeConstant) - vesselGraph->countVesselNodes);
					float perturbation2 = 1e-3;
					for( int xi=0; xi<vesselGraph->DOMAIN_SIZE_X/latticeConstant;xi++)
						for( int yi=0; yi<vesselGraph->DOMAIN_SIZE_Y/latticeConstant; yi++)
							for( int zi=0; zi<vesselGraph->DOMAIN_SIZE_Z/latticeConstant; zi++)
							if( vesselGraph->octree->at( xi, yi, zi) == 0){
								//fprintf( stderr,"found some space\n");
								t->add( new VoronoiCell(
												xi+0.5+perturbation2*(RAND01-0.5),
												yi+0.5+perturbation2*(RAND01-0.5),
												zi+0.5+perturbation2*(RAND01-0.5)));
							}
					fprintf( stderr, "Added %i Spacefilling Points to Triangulation\n", t->numberOfVertices() - vesselGraph->countVesselNodes);

					t->setFramePoints();
					fprintf( stderr, "Triangulate\n");
					t->triangulate();*/
					t = createVoronoiTesselation( vesselGraph, 5);
				}
				// READ DATA
				for( int i=0; i<t->numberOfVertices(); i++){
					fread( &markerp,  sizeof(float), 1, fpi);
					t->get(i)->conc = markerp;
					//if(markerp!=0)
						//fprintf(stderr, "%e\n", markerp);
				}

				// OUTPUT to POVRAY
				int margin = 1;
				Domain<3> *domain = new Domain<3>(margin, (vesselGraph->DOMAIN_SIZE_X-margin), margin, (vesselGraph->DOMAIN_SIZE_Y-margin), margin, (vesselGraph->DOMAIN_SIZE_Z-margin));
				char povray_filename[512];
				sprintf( povray_filename, "%s.pov", binary_filename);
				t->printVoronoiDiagramToPovray3D( povray_filename, domain);
				sprintf(povray_filename, "nice povray -D0 +h480 +w640 +ua %s.pov 2> /dev/null &", binary_filename);
				system(povray_filename);

				delete domain;
			}break;
		}

		fclose(fpi);
	}else

	if( strstr( binary_filename, ".Cp.binary") && !strstr( binary_filename, ".Cp.binary.")){

		if(true){
			fprintf( stderr, "Read plasma data from binary file       (%s)\n", binary_filename);
			vesselGraph->readMarkerIntensityFromBinary( binary_filename);
			char povray_filename[512];
			sprintf( povray_filename, "%s.pov", binary_filename);
			vesselGraph->printToPovray(povray_filename, 1);
		}else{
			FILE *fpp = fopen( binary_filename, "rb");

			// Binary Header
			int DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z;
			fread( &DOMAIN_SIZE_X,  sizeof(int), 1, fpp);
			fread( &DOMAIN_SIZE_Y,  sizeof(int), 1, fpp);
			fread( &DOMAIN_SIZE_Z,  sizeof(int), 1, fpp);

			//fprintf( stderr, "-> Dimension: %ix%ix%i\n", DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);

			float markeri;
			for( int x=0; x<DOMAIN_SIZE_X; x++)
				for( int y=0; y<DOMAIN_SIZE_Y; y++)
					for( int z=0; z<DOMAIN_SIZE_Z; z++){
						fread( &markeri,  sizeof(float), 1, fpp);
						VesselNode *vn = vesselGraph->octree->at(x,y,z);
						if(vn!=0){
							vn->marker = markeri;
						}
						//if(markeri!=0)
						//printf("marker at (%i,%i,%i) = %e\n", x,y,z, markeri);
					}
			fclose(fpp);
			vesselGraph->printMarkerIntensityToBinary( binary_filename);
		}

		/*for (int i = 0; i < vesselGraph->countVesselNodes; i++){
			float Cp;
			fread( &Cp,  sizeof(float), 1, fpp);
			vesselGraph->vesselNodes[i]->marker = Cp;
		}*/


	}

	else{
		//fprintf( stderr, "WARNING: %s is no accepted input file\n", filename);
	}

}

int reversealphasort(const struct dirent **pa, const struct dirent **pb){
	return -strcmp((**pa).d_name, (**pb).d_name);
}

int perfusion(int argc, char** argv, VesselGraph *vesselGraph) {

	//int DOMAIN_SIZE_X = 1;
	//int DOMAIN_SIZE_Y = 1;
	//int DOMAIN_SIZE_Z = 1;
	//int dimensions = 1;
	char dirname[512] = ".";

	int x=1, y=1, z=1;
	float radius = 0;
	float permeability = 1.;
	float parallelVessels = 1.;

	Tumor *tumor = 0;
	//float radiusNecroticCore = 5;
	//float radiusTumor = 15;
	int ***cells = 0;
	VoronoiDiagram *t = 0;
	VesselGraph::Types testtype = VesselGraph::RegularLattice;

	float time = 0.;
	int maxit = 200; // remodelling

	// PERFUSION
	double time_step = 0.; // contrast-agent perfusion
	double time_end = 180.; // sec
	double outputResolution = 1.; //sec
	double diffusionCoefficient = 0;

	double flow=0;
	bool UpdateInterstitialCompartment = true;
	bool UpdatePlasmaCompartment = true;
	double PlasmaConc = 0.;
	double InterstitialConc = 0.;

	char c;
	while ((c = getopt(argc, argv, "ABPf::p:x:y:z:r:n:T:NV:OI:i:S:1:2:R:t:e:o:d:C:D:X:c:W:E::sa:b:F9:")) != -1)
		switch (c) {
		case 'F':{
			// TEST: Distance to root node

			// init
			int countInterTipConnections = 0;
			vesselGraph->setInterTipConnections( tumor, countInterTipConnections);

			vesselGraph->updateRadius();
			vesselGraph->updatePressureNEW2();
			vesselGraph->updateFlow();

			int changed = 0;
			for( int n=0; n<vesselGraph->countVesselNodes; n++)
			{
				vesselGraph->vesselNodes[n]->marker = -1;
				if( vesselGraph->vesselNodes[n]->pressure == 10){
					// distance
					vesselGraph->vesselNodes[n]->marker = 0;
					// delay
					//vesselGraph->vesselNodes[n]->marker = vesselGraph->vesselNodes[n]->flo;
					changed++;
				}
			}

			// set distance
			while( changed<vesselGraph->countVesselNodes)
			{
				for( int n=0; n<vesselGraph->countVesselNodes; n++)
				if( vesselGraph->vesselNodes[n]->marker == -1)
				{
					for( int j=0; j<vesselGraph->vesselNodes[n]->countNeighbors; j++){
						if( vesselGraph->vesselNodes[n]->marker < vesselGraph->vesselNodes[n]->neighbors[j]->marker){
							if( vesselGraph->vesselNodes[n]->marker == -1){
								changed++;
							}
							// distance
							vesselGraph->vesselNodes[n]->marker = vesselGraph->vesselNodes[n]->neighbors[j]->marker + 1;
							// delay
							//vesselGraph->vesselNodes[n]->marker = vesselGraph->vesselNodes[n]->neighbors[j]->marker
							//		+ vesselGraph->getVascularVolume(vesselGraph->vesselNodes[n]->branches[j]) / fabs(vesselGraph->vesselNodes[n]->branches[j]->flow);
						}
					}
				}
				fprintf( stderr, "changed %i of %i\n", changed, vesselGraph->countVesselNodes);

			}

			// output distance
			char filename[512];
			sprintf( filename, "distanceToRoot.dat");
			FILE *fp=fopen( filename, "w+");
			//fprintf( fp, "#x,y,z,distance to root\n");
			fprintf( fp, "\n");
			for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++){
				for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++)
					for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
						VesselNode *vn = vesselGraph->octree->at(x,y,z);
						if( vn)
							fprintf( fp, "%i %i %e\n", x,y,vn->marker);
						else
							fprintf( fp, "%i %i %e\n", x,y,0.);
					}
				fprintf( fp, "\n");
			}
			fclose(fp);

			vesselGraph->removeInterTipConnections( countInterTipConnections);

		}break;
		case 'a':{
				UpdateInterstitialCompartment = false;
				InterstitialConc = atof( optarg);
				UpdatePlasmaCompartment = true;
		}break;
		case 'b':{
				UpdateInterstitialCompartment = true;
				UpdatePlasmaCompartment = false;
				PlasmaConc = atof( optarg);
		}break;
		case 'f':{
			//for (int i = 0; i < vesselGraph->countVesselSegments; i++)
			//	vesselGraph->vesselSegments[i]->flow = atof(optarg);
			flow = atof(optarg);
		}break;
		case 'E':{

			DIR *dirp;
			struct dirent *dptr;

			// Read Vascularization
			if((dirp=opendir( dirname))==NULL){perror("Error"); exit(1);}
			while( dptr=readdir(dirp) ){
				//printf("%s\n",dptr->d_name);
				if( strstr( dptr->d_name, "vasc.xml")){
					fprintf( stderr, "Create Vascularization from file\n");
					vesselGraph = new VesselGraph( dptr->d_name, tumor);
				}
			}
			closedir(dirp);

			CONCENTRATION_T ***marker = 0;
			struct dirent **namelist;
			int n = scandir(dirname, &namelist, 0, reversealphasort);
			while (n--) {
				//printf("%s\n", namelist[n]->d_name);
				if( optarg==0 || strstr( namelist[n]->d_name, optarg)){
					//printf("%s\n", namelist[n]->d_name);
					readBinaryAndExportToPovray( namelist[n]->d_name, vesselGraph, marker, t);
				}
				free(namelist[n]);
			}
			free(namelist);
			//free();

		}break;
		case 'W':
				{
					char filename[512];
					sprintf( filename, "vesselGraph%i.pov", atoi(optarg));
					vesselGraph->printToPovray( filename, (char)atof(optarg));
				}
				break;
		case 'B':
			testtype =  VesselGraph::IrregularLattice; break;
		case 'A':
		{
			fprintf(stderr, "START TEST MODE\n");

			// Create Vascularization
			vesselGraph = new VesselGraph( 15, 20, 5, 3);
			vesselGraph->LATTICE_CONSTANT=5;
			vesselGraph->setSymmetricVesselGraph(4);

			//vesselGraph->setType( VesselGraph::RegularLattice);
			vesselGraph->setType( testtype);

			switch( vesselGraph->type()){

			case VesselGraph::IrregularLattice:{
				float perturbation = 5e-1;
				for( int i=0; i<vesselGraph->countVesselNodes;i++){
					vesselGraph->vesselNodes[i]->position[0] += perturbation * RAND01;
					vesselGraph->vesselNodes[i]->position[1] += perturbation * RAND01;
					vesselGraph->vesselNodes[i]->position[2] += perturbation * RAND01;
				}

				// Update Vessel Properties

				int foo = 0;
				vesselGraph->setInterTipConnections( tumor, foo);
				fprintf(stderr, "countInterTipConnections=%i\n", foo);
				for( int i=0; i<vesselGraph->countVesselSegments;i++){
					vesselGraph->vesselSegments[i]->radius=2.5;
					vesselGraph->vesselSegments[i]->radiusStatic=false;
				}
				vesselGraph->updatePressureNEW2();
				vesselGraph->updateFlow();
				vesselGraph->updateShear();

				for( int i=0; i<vesselGraph->countVesselNodes;i++)
					if(vesselGraph->vesselNodes[i]->getType() == TIP)
						fprintf( stderr, "p=%e\n", vesselGraph->vesselNodes[i]->pressure);
				vesselGraph->printToEPS("testV.eps",0,0);
				vesselGraph->removeInterTipConnections( foo);

				// Create Triangulation
				t = createVoronoiTesselation( vesselGraph, vesselGraph->LATTICE_CONSTANT);

				// Calculate Volume & Interfaces
				t->setVolume();

				// Write to Povray
				int margin=1;
				Domain<3> *domain = new Domain<3>(margin, (vesselGraph->DOMAIN_SIZE_X-margin), margin, (vesselGraph->DOMAIN_SIZE_Y-margin), margin, (vesselGraph->DOMAIN_SIZE_Z-margin));
				t->printVoronoiDiagramToPovray("test.pov", domain);
				t->printVoronoiDiagramToEps("test.eps");
			}
			break;


			case VesselGraph::RegularLattice:{
				cells = (int***) malloc( vesselGraph->DOMAIN_SIZE_X * sizeof(int**));
				for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++){
					cells[x] = (int**) malloc( vesselGraph->DOMAIN_SIZE_Y * sizeof(int*));
					for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++){
						cells[x][y] = (int*) malloc( vesselGraph->DOMAIN_SIZE_Z * sizeof(int));
						for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
							cells[x][y][z] = 0;
						}
					}
				}
			}
			break;

			}
			//fprintf(stderr, "TEST: Octree\n");
			//Octree<float> *o = new Octree<float>(8);
			/*o->set( 0,0,0,1);
			o->set( 0,0,0.1,2);
			o->set( 0,0,0.9,3);
			o->set( 0,0,1,4);*/
			/*o->set( 0,0,(int)0,1);
			o->set( 0,0,(int)0.1,2);
			o->set( 0,0,(int)0.9,3);
			o->set( 0,0,(int)1,4);
			fprintf(stderr, "(0,0,0) ==> %f\n", o->at( 0,0,0));
			fprintf(stderr, "(0,0,0.1) ==> %f\n", o->at( 0,0,0.1));*/

			fprintf(stderr, "END TEST MODE\n");
		}break;
		case 'd':
			diffusionCoefficient = atof(optarg);
			if( vesselGraph)
				vesselGraph->diffusionCoefficient = diffusionCoefficient;
			fprintf(stderr, "Set Diffusion Coefficient to %e\n", diffusionCoefficient);
			break;
		case 't':
			time_step = atof(optarg);
			fprintf(stderr, "Set Time Step to %e\n", time_step);
			break;
		case 'e':
			time_end = atof(optarg);
			fprintf(stderr, "Set End Time to %e\n", time_end);
			break;
		case 'o':
			outputResolution = atof(optarg);
			fprintf(stderr, "Set Temporal Output Intervall Size to %e\n", outputResolution);
			break;
		case 'X':
			vesselGraph->writeToXML( optarg, 0, tumor);
			break;
		case 'p':
			permeability = atof(optarg);
			if( vesselGraph)
				vesselGraph->permeability = permeability;
			break;
		case 'P':
			parallelVessels = atof(optarg);
			break;
		case 'x':
			x = atoi(optarg);
			break;
		case 'y':
			y = atoi(optarg);
			break;
		case 'z':
			z = atoi(optarg);
			break;
		case 'r':
			radius=atof(optarg);
			break;
		case 'c':
		{
			cprintf( stderr, 7, "Read Cells from %s\n", optarg);

			switch(vesselGraph->type()){

			case VesselGraph::RegularLattice:{
				FILE *fp = fopen(optarg, "r");
				cells = (int***) malloc( vesselGraph->DOMAIN_SIZE_X * sizeof(int**));
				for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++){
					cells[x] = (int**) malloc( vesselGraph->DOMAIN_SIZE_Y * sizeof(int*));
					for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++){
						cells[x][y] = (int*) malloc( vesselGraph->DOMAIN_SIZE_Z * sizeof(int));
						for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
							cells[x][y][z] = 0;
						}
					}
				}
				int x=0,y=0,z=0,index=0;
				while(  EOF!=fscanf(fp, "%i, %i, %i\t%i\n", &x,&y,&z,&index)){
					cells[x][y][z] = index;
					//fprintf(stderr, "Read: %i, %i, %i\t%i\n", x,y,z,index);
				}

				fclose(fp);

				fprintf( stderr, "Size: %lf x %lf x %lf\n",
						vesselGraph->DOMAIN_SIZE_X*vesselGraph->LATTICE_CONSTANT,
						vesselGraph->DOMAIN_SIZE_Y*vesselGraph->LATTICE_CONSTANT,
						vesselGraph->DOMAIN_SIZE_Z*vesselGraph->LATTICE_CONSTANT);
			}break;

			case VesselGraph::IrregularLattice:{

				fprintf( stderr, "Add %i Vessel Nodes to Triangulation\n", vesselGraph->countVesselNodes);

				// CREATE TRIANGULATION
				/*t = new VoronoiDiagram();

				// CELL LATTICE CONSTANT
				t->LATTICE_CONSTANT = 5;

				// ADD VESSEL NODES

				float perturbation1 = 1e-3;UpdateInterstitialCompartment
				for( int i=0; i<vesselGraph->countVesselNodes;i++)
					//t->add( new VoronoiCell( vesselGraph->vesselNodes[i]));
					t->add( new VoronoiCell(
									vesselGraph->vesselNodes[i]->position[0]/t->LATTICE_CONSTANT+perturbation1*RAND01,
									vesselGraph->vesselNodes[i]->position[1]/t->LATTICE_CONSTANT+perturbation1*RAND01,
									vesselGraph->vesselNodes[i]->position[2]/t->LATTICE_CONSTANT+perturbation1*RAND01));
				// ADD
				//fprintf( stderr, "Add %i Spacefilling Points to Triangulation\n",
				//		(int)(vesselGraph->DOMAIN_SIZE_X/latticeConstant*vesselGraph->DOMAIN_SIZE_Y/latticeConstant*vesselGraph->DOMAIN_SIZE_Z/latticeConstant) - vesselGraph->countVesselNodes);
				float perturbation2 = 1e-3;
				for( int xi=0; xi<vesselGraph->DOMAIN_SIZE_X/t->LATTICE_CONSTANT;xi++)
					for( int yi=0; yi<vesselGraph->DOMAIN_SIZE_Y/t->LATTICE_CONSTANT; yi++)
						for( int zi=0; zi<vesselGraph->DOMAIN_SIZE_Z/t->LATTICE_CONSTANT; zi++)
						if( vesselGraph->octree->at( xi, yi, zi) == 0){
							//fprintf( stderr,"found some space\n");
							t->add( new VoronoiCell(
											xi+0.5+perturbation2*(RAND01-0.5),
											yi+0.5+perturbation2*(RAND01-0.5),
											zi+0.5+perturbation2*(RAND01-0.5)));
						}

				t->setFramePoints();
				fprintf( stderr, "Triangulate\n");
				t->triangulate();*/
				t = createVoronoiTesselation( vesselGraph, 5);
				fprintf( stderr, "Added %i Spacefilling Points to Triangulation\n", t->numberOfVertices() - vesselGraph->countVesselNodes);

				fprintf( stderr, "Calculate Volume and Interfaces of Voronoi Cells\n");
				t->setVolume();

			}break;
			}
		}
		break;
		case 'V':{
			if( strcmp(".",optarg) == 0){
				cprintf( stderr, 7, "Create new empty Vascularization\n");

				int dimensions = (z>1 ? 3 : (y>1 ? 2 : 1));
				vesselGraph = new VesselGraph( x, y, z, dimensions);
				vesselGraph->LATTICE_CONSTANT = 60;
				vesselGraph->permeability = permeability;
				vesselGraph->diffusionCoefficient = diffusionCoefficient;
				fprintf( stderr, "graph dimensions: %i x %i x %i\n",
						vesselGraph->DOMAIN_SIZE_X,vesselGraph->DOMAIN_SIZE_Y,vesselGraph->DOMAIN_SIZE_Z);
			}
			else{
				cprintf( stderr, 7, "Create Vascularization from file\n");
				vesselGraph = new VesselGraph( optarg, tumor);
				fprintf( stderr, "graph dimensions: %i x %i x %i\n",
						vesselGraph->DOMAIN_SIZE_X,vesselGraph->DOMAIN_SIZE_Y,vesselGraph->DOMAIN_SIZE_Z);

				// TEST: Vessel Surface per Volume
				double surface=0;
				double plasmaVolume=0;
				for( int i=0; i<vesselGraph->countVesselSegments; i++){
					surface      += vesselGraph->getVascularSurface( vesselGraph->vesselSegments[i]);
					plasmaVolume += vesselGraph->getVascularVolume(  vesselGraph->vesselSegments[i]);
				}
				fprintf( stderr, "E(S/V) = %e microns^-1\n", surface / (vesselGraph->DOMAIN_SIZE_X*vesselGraph->DOMAIN_SIZE_Y*vesselGraph->DOMAIN_SIZE_Z*60.*60.*60.));
				fprintf( stderr, "E(S/V_P) = %e microns^-1\n", surface / plasmaVolume);


			}
			fprintf( stderr, "Vessel Graph Topology: %s\n",
					(vesselGraph->type() == VesselGraph::IrregularLattice ? "IrregularLattice" : "RegularLattice"));
			fprintf( stderr, "Lattice Constant: %lf\n", vesselGraph->latticeConstant());
			fprintf( stderr, "Size: %lf x %lf x %lf\n",
					vesselGraph->DOMAIN_SIZE_X*vesselGraph->LATTICE_CONSTANT,
					vesselGraph->DOMAIN_SIZE_Y*vesselGraph->LATTICE_CONSTANT,
					vesselGraph->DOMAIN_SIZE_Z*vesselGraph->LATTICE_CONSTANT);

		}
			break;

		case '9':
			{
				fprintf( stderr, "do some temporay stuff\n");
							radius = atoi(optarg);
							int mx[4] = {25,25,75,75};
							int my[4] = {25,75,75,25};
							double flow[4] ={0,0,0,0};
							double phi_p[4]={0,0,0,0};
							double k_ps[4] ={0,0,0,0};
							int count[4] ={0,0,0,0};
							for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++){
								for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++){
									for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
										for( int i=0; i<4; i++){
											// POINT INSIDE QUADRANT ?
											if(pow(mx[i]-x,2) + pow(my[i]-y,2) <= radius*radius){
												count[i] ++;
												phi_p[i] += vesselGraph->getVascularVolume(x,y,z)/60./60./60.;
												k_ps[i]  += vesselGraph->getVascularPermeabilitySurfaceProduct(x,y,z);

												VesselNode *vn = vesselGraph->octree->at(x,y,z);
												for( int n=0; vn!=0 && n<vn->countNeighbors; n++)
												{
													// NEIGHBOR POINT OUTSIDE QUADRANT ?
													if( pow(mx[i]-vn->neighbors[n]->position[0],2) + pow(my[i]-vn->neighbors[n]->position[1],2) > radius*radius)
														flow[i] += fabs( vn->branches[n]->flow);
												}
											}
										}
									}
								}
							}
							for( int i=0; i<4; i++){
								fprintf( stderr, "%i\t%e\t%e\t%e\t%i\n", i, phi_p[i] / count[i], k_ps[i], flow[i], count[i]);
							}
							fprintf( stderr, "average over %i\n", count[0]);
							exit(0);

			}
			break;

		case 's':{
			vesselGraph->simplifyNetwork( 2.5);
			//vesselGraph->simplifyNetwork( 7);
			float minSegmentLength = 1000;
			for( int i=0; i<vesselGraph->countVesselSegments; i++){
				float length = vesselGraph->distance(vesselGraph->vesselSegments[i]->vesselNodes[0],vesselGraph->vesselSegments[i]->vesselNodes[1]);
				if( minSegmentLength>length)
					minSegmentLength=length;
			}
			fprintf( stderr, "minLength=%lf\n", minSegmentLength);
		}break;
		case 'i':{
			// Already Read a Vascularization?
			if( vesselGraph){
				FILE *fp_in, *fp_out;
				char filename_in[512];
				char filename_out[512];
				int i=0;

				sprintf( filename_in, optarg, i);
				sprintf( filename_out, "%s.ALL.dat", filename_in);
				fp_out = fopen( filename_out, "w+");
				while( (fp_in = fopen( filename_in, "rb+")) != 0){
					fprintf( stderr, "READ: file %s\n", filename_in);
					int valuei=0;
					fread( &valuei, sizeof(int), 1, fp_in);
					if(valuei!=vesselGraph->DOMAIN_SIZE_X){fprintf(stderr, "wrong x-dimension of input file!\n"); exit(0);}
					fread( &valuei, sizeof(int), 1, fp_in);
					if(valuei!=vesselGraph->DOMAIN_SIZE_Y){fprintf(stderr, "wrong y-dimension of input file!\n"); exit(0);}
					fread( &valuei, sizeof(int), 1, fp_in);
					if(valuei!=vesselGraph->DOMAIN_SIZE_Z){fprintf(stderr, "wrong z-dimension of input file!\n"); exit(0);}
					for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++){
						for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++)
							for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
								float valuef = 0;
								fread( &valuef, sizeof(float), 1, fp_in);
								fprintf( fp_out, "%i %i %i %i %e\n", i, x, y, z, valuef);
							}
						fprintf( fp_out, "\n");
					}
					fprintf( fp_out, "\n\n");

					fclose(fp_in);
					i++;
					sprintf( filename_in, optarg, i);
				}
				fclose(fp_out);
			}
		}
		break;
		case 'D':{
			sprintf( dirname, "%s/", optarg);
			if( mkdir(dirname, S_IRWXU)!=0){
				fprintf(stderr, "WARNING: Can not create directory %s: ", optarg);
				perror("");
				//perror("mkdir() error");
				//exit(0);
			}else{
				fprintf(stderr, "Create Directory %s\n", optarg);
			}

			chdir(dirname);


			FILE *fp;
			char filename[512];
			sprintf( filename, "commandline.dat", dirname);
			fprintf(stderr, "write command line to %s\n", filename);
			fp = fopen( filename, "a+");
			if( !fp){
				fprintf(stderr, "Couldn't create file\n");
				perror("sd");
			}
			time_t rawtime;
			mytime( &rawtime );
			fprintf( fp, "\n%s", asctime(localtime(&rawtime)));
			for( int i=0; i<argc; i++)
				fprintf( fp, "%s ", argv[i]);
			fprintf( fp, "\n");
			fclose(fp);
		}
			break;
		/*case 'n':
			radiusNecroticCore=atof(optarg);
			break;*/
		case 'T':{
			//float tumorCenter[3] = { x, y, z};
			float value = atof(&optarg[strcspn(optarg, "0123456789.")]);
			if( !tumor)
				tumor = new Tumor();
			if(strstr(optarg,"x")){
				tumor->setX( value); fprintf(stderr, "set tumor to x=%f\n", value); }else
			if(strstr(optarg,"y"))
				tumor->setY( value); else
			if(strstr(optarg,"z"))
				tumor->setZ( value); else
			if(strstr(optarg,"r"))
				tumor->setRadius( value); else
			if(strstr(optarg,"p"))
				tumor->setPermeability( value); else
			if(strstr(optarg,"c"))
				tumor->setParallelVessels( value);
			/*tumor->setRadius( radius);
			//tumor->setRadiusNecroticCore( radiusNecroticCore);
			tumor->setTopology( Tumor::SPHERICAL);
			tumor->setCenter( tumorCenter);
			tumor->setPermeability( permeability);
			tumor->setParallelVessels( parallelVessels);*/
		}
			break;
		case 'N':
			tumor->setRadiusNecroticCore( radius);
			break;

		case 'I':{
			if( !vesselGraph){
				fprintf( stderr, "Create new empty Vascularization\n");

				/*DOMAIN_SIZE_X = x;
				DOMAIN_SIZE_Y = y;
				DOMAIN_SIZE_Z = z;*/
				int dimensions = (z>1 ? 3 : (y>1 ? 2 : 1));
				vesselGraph = new VesselGraph( x, y, z, dimensions);
				vesselGraph->permeability = permeability;
				vesselGraph->LATTICE_CONSTANT=60;
			}
			float value = atof(&optarg[strcspn(optarg, "0123456789.")]);

			fprintf(stderr, "\n\nINITIALIZE BLOOD VESSEL NETWORK\n\n");

			// Set Initial Vascularization
			if(strstr(optarg,"symmetric")){
				fprintf(stderr, "Set Symmetric Vascularization with Tree Depth: %d\n", (int)value);
				vesselGraph->setSymmetricVesselGraph( (int)value);
			}else
			if(strstr(optarg,"singleStraight")){
				fprintf(stderr, "Set Single (Straight) Vessel\n");
				if( value==0)
					vesselGraph->setSingleVessel(vesselGraph->DOMAIN_SIZE_X-1, vesselGraph->DOMAIN_SIZE_Y/2, vesselGraph->DOMAIN_SIZE_Z/2);
				else
					vesselGraph->setSingleVessel(vesselGraph->DOMAIN_SIZE_X-1, vesselGraph->DOMAIN_SIZE_Y/2, (int)value);
			}else
			if(strstr(optarg,"singleRandom")){
				fprintf(stderr, "Set Single (Random) Vessel\n");
				if( value==0)
					vesselGraph->setSingleRandomVessel(vesselGraph->DOMAIN_SIZE_X-1, vesselGraph->DOMAIN_SIZE_Y/2, vesselGraph->DOMAIN_SIZE_Z/2);
				else
					vesselGraph->setSingleRandomVessel(vesselGraph->DOMAIN_SIZE_X-1, vesselGraph->DOMAIN_SIZE_Y/2, (int)value);
			}else{
				fprintf(stderr, "Set Homogenous Vessel Network (mode=%d)\n", (int)value);
				vesselGraph->setInitialVesselGraph( (int)value);
			}

			// Update Vessel Properties
			vesselGraph->updateRadius();
			int countInterTipConnections = 0;
			vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
			fprintf(stderr, "countInterTipConnections=%i\n", countInterTipConnections);
			vesselGraph->updatePressureNEW2();
			vesselGraph->updateFlow();
			vesselGraph->updateShear();

			// Output Vascularization
			if( vesselGraph->dimensions == 3)
				vesselGraph->printToPovray( "beforeRemodelling.pov", ((char)0));
			if( vesselGraph->dimensions < 3){
				vesselGraph->printToEPS("beforeRemodelling.eps", 0, NULL);
				if(tumor)tumor->printToEPS( "beforeRemodelling.eps");
			}
			vesselGraph->removeInterTipConnections( countInterTipConnections);

			for( int i=0; i<vesselGraph->countVesselNodes;i++)
				if(vesselGraph->vesselNodes[i]->getType() == TIP)
					fprintf( stderr, "p=%e\n", vesselGraph->vesselNodes[i]->pressure);

		}
			break;

		case 'S':
			fprintf(stderr,
					"\n\nINITIALIZE SYMMETRIC BLOOD VESSEL NETWORK\n\n");
			vesselGraph->setSymmetricVesselGraph( atoi(optarg));
			vesselGraph->updateRadius();
			{
				int countInterTipConnections = 0;
				vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
				vesselGraph->updatePressureNEW2();
				vesselGraph->updateFlow();
				vesselGraph->updateShear();
				fprintf(stderr, "countInterTipConnections=%i\n",
						countInterTipConnections);
				//vesselGraph->printToPovray( "before.pov", 0);
				vesselGraph->printToEPS("beforeRemodelling.eps", 0, NULL);
				if(tumor)tumor->printToEPS( "beforeRemodelling.eps");
				vesselGraph->removeInterTipConnections(	countInterTipConnections);
			}
			break;

		case '1':
			fprintf(stderr,
					"\n\nINITIALIZE SINGLE BLOOD VESSEL\n\n");
			//vesselGraph->setSingleVessel( atoi(optarg), (dimensions>=2?y:0), (dimensions>=3?z:0));
			vesselGraph->setSingleVessel( atoi(optarg), y, z);
			vesselGraph->updateRadius();
			{
				int countInterTipConnections = 0;
				vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
				vesselGraph->updatePressureNEW2();
				vesselGraph->updateFlow();
				vesselGraph->updateShear();
				fprintf(stderr, "countInterTipConnections=%i\n",
						countInterTipConnections);
				/*for( int i=0; i<vesselGraph->countVesselSegments; i++)
					fprintf( stdout, "%f %f %e %e %e %e \n",
							vesselGraph->vesselSegments[i]->vesselNodes[0]->position[0],
							vesselGraph->vesselSegments[i]->vesselNodes[1]->position[0],
							vesselGraph->vesselSegments[i]->vesselNodes[0]->pressure,
							vesselGraph->vesselSegments[i]->vesselNodes[1]->pressure,
							vesselGraph->vesselSegments[i]->radius,
							vesselGraph->vesselSegments[i]->flow);*/

				//vesselGraph->printToPovray( "before.pov", 0);
				vesselGraph->printToEPS("beforeRemodelling.eps", 0, NULL);
				if(tumor)tumor->printToEPS( "beforeRemodelling.eps");
				vesselGraph->removeInterTipConnections( countInterTipConnections);
			}
			break;

		case '2':
			fprintf(stderr,
					"\n\nINITIALIZE RANDOM BLOOD VESSEL\n\n");
			//vesselGraph->setSingleRandomVessel( atoi(optarg), (dimensions>=2?y:0), (dimensions>=3?z:0));
			vesselGraph->setSingleRandomVessel( atoi(optarg), y, z);
			vesselGraph->updateRadius();
			{
				int countInterTipConnections = 0;
				vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
				vesselGraph->updatePressureNEW2();
				vesselGraph->updateFlow();
				vesselGraph->updateShear();
				fprintf(stderr, "countInterTipConnections=%i\n",
						countInterTipConnections);
				/*for( int i=0; i<vesselGraph->countVesselSegments; i++)
					fprintf( stdout, "%f %f %e %e %e %e \n",
							vesselGraph->vesselSegments[i]->vesselNodes[0]->position[0],
							vesselGraph->vesselSegments[i]->vesselNodes[1]->position[0],
							vesselGraph->vesselSegments[i]->vesselNodes[0]->pressure,
							vesselGraph->vesselSegments[i]->vesselNodes[1]->pressure,
							vesselGraph->vesselSegments[i]->radius,
							vesselGraph->vesselSegments[i]->flow);*/

				//vesselGraph->printToPovray( "before.pov", 0);
				vesselGraph->printToEPS("beforeRemodelling.eps", 0, NULL);
				if(tumor)tumor->printToEPS( "beforeRemodelling.eps");
				vesselGraph->removeInterTipConnections( countInterTipConnections);
			}
			break;

		case 'R':
			fprintf(stderr, "\n\nPERFORM SHEAR STRESS HOMOGENIZATION\n\n");
			maxit = atoi(optarg);
			fprintf(stderr, "Set Remodeling Iterations to %i\n", maxit);
			for (int it = 0; it < maxit; it++) {
				fprintf(stderr, "\rRemodeling: %i / %i \b", it + 1, maxit);
				vesselGraph->remodel( tumor);
			}

			break;
		case 'O':
		{
			vesselGraph->printToPovray("pressureBefore.pov", (char)0);
			fprintf( stderr, "Set Inter-Tip-Connections\n");
			int countInterTipConnections = 0;
			vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
			fprintf( stderr, "Update Radii\n");
			vesselGraph->updateRadius();
			fprintf( stderr, "Update Pressure\n");
			vesselGraph->updatePressureNEW2();
			fprintf( stderr, "Update Flow\n");
			vesselGraph->updateFlow();
			fprintf( stderr, "Update Shear Stress\n");
			vesselGraph->updateShear();
			// SET PERMEABILITY
			for (int i = 0; i < vesselGraph->countVesselSegments; i++) {
				if( tumor &&
					(tumor->isTumor(vesselGraph->vesselSegments[i]->vesselNodes[0]->position) ||
					tumor->isTumor(vesselGraph->vesselSegments[i]->vesselNodes[1]->position)))
					vesselGraph->vesselSegments[i]->permeability = tumor->getPermeability();
				else
					vesselGraph->vesselSegments[i]->permeability = vesselGraph->permeability;
			}

			//if( vesselGraph->dimensions < 3)
			{
			vesselGraph->printFlowMap((char*) "flow.60x60x60.eps", vesselGraph, 60., 60., 60.);
			vesselGraph->printBloodVolumeFractionMap((char*) "f_p.60x60x60.eps", vesselGraph,60., 60., 60.);
			vesselGraph->printInterstitialSpaceVolumeFractionMap((char*) "f_i.60x60x60.eps",vesselGraph, 60., 60., 60.);
			vesselGraph->printPermeabilityMap((char*) "Kps.60x60x60.eps", vesselGraph, 60.,60., 60.);
			vesselGraph->printWallShearStressMap((char*) "shear.60x60x60.eps", vesselGraph,	60.,60., 60.);
			}
			vesselGraph->printFlowMap((char*) "flow.300x300x3000.eps", vesselGraph, 300., 300., 3000.);
			vesselGraph->printBloodVolumeFractionMap((char*) "f_p.300x300x3000.eps", vesselGraph, 300., 300., 3000.);
			vesselGraph->printInterstitialSpaceVolumeFractionMap((char*) "f_i.300x300x3000.eps", vesselGraph, 300., 300.,3000.);
			vesselGraph->printPermeabilityMap((char*) "Kps.300x300x3000.eps", vesselGraph,	300., 300., 3000.);
			vesselGraph->printWallShearStressMap((char*) "shear.300x300x3000.eps", vesselGraph,	300., 300., 3000.);


			if( vesselGraph->dimensions == 3){
				//char filename[512] = "pressureAfter.pov";
				vesselGraph->printToPovray("pressureAfter.pov", (char)0);
				vesselGraph->printToPovray("flow.pov", 2);
				vesselGraph->printToPovray("velocity.pov", 3);
				//vesselGraph->printToPovray("plasmaVolume.pov", 4);
				//vesselGraph->printToVRML("pressure.wrl", 0);
				//vesselGraph->printToVRML("flow.wrl", 2);
			}
			if( vesselGraph->dimensions < 3){
				vesselGraph->printToEPS("afterRemodelling.eps", 0, NULL);
				if(tumor)tumor->printToEPS( "afterRemodelling.eps");
			}

			vesselGraph->removeInterTipConnections(countInterTipConnections);
		}

			break;

		case 'C': {
			// TODO: Contrast Agent Perfusion
			{
				int countInterTipConnections = 0;
				vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
				fprintf(stderr, "countInterTipConnections=%i\n",
						countInterTipConnections);

				vesselGraph->updateRadius();
				for (int i = 0; i < vesselGraph->countVesselSegments; i++){
					if( vesselGraph->vesselSegments[i]->radius != 4){
						fprintf(stderr, "radius=%e\n", vesselGraph->vesselSegments[i]->radius);
					}
					if( vesselGraph->vesselSegments[i]->vesselNodes[0]->position[0]==vesselGraph->vesselSegments[i]->vesselNodes[1]->position[0] &&
						vesselGraph->vesselSegments[i]->vesselNodes[0]->position[1]==vesselGraph->vesselSegments[i]->vesselNodes[1]->position[1] &&
						vesselGraph->vesselSegments[i]->vesselNodes[0]->position[2]==vesselGraph->vesselSegments[i]->vesselNodes[1]->position[2]	)
						exit(0);
				}
				vesselGraph->updatePressureNEW2();
				vesselGraph->updateFlow();
				vesselGraph->updateShear();
				vesselGraph->printToEPS("testV2.eps",0,0);
			}

			// SET PERMEABILITY
			for (int i = 0; i < vesselGraph->countVesselSegments; i++) {
				if(flow!=0)
					vesselGraph->vesselSegments[i]->flow = flow;//1e5;
				if(radius!=0)
					vesselGraph->vesselSegments[i]->radius = radius;
				if( tumor &&
					(tumor->isTumor(vesselGraph->vesselSegments[i]->vesselNodes[0]->position) ||
					tumor->isTumor(vesselGraph->vesselSegments[i]->vesselNodes[1]->position)))
					vesselGraph->vesselSegments[i]->permeability = tumor->getPermeability();
				else
					vesselGraph->vesselSegments[i]->permeability = vesselGraph->permeability;

			}

			// FLOW
			double **markerVessel = (double**) malloc(		sizeof(double*) * vesselGraph->DOMAIN_SIZE_X);
			double **markerVessel2 = (double**) malloc(		sizeof(double*) * vesselGraph->DOMAIN_SIZE_X);
			CONCENTRATION_T ***markerTissue = (CONCENTRATION_T***) malloc(		sizeof(CONCENTRATION_T**) * vesselGraph->DOMAIN_SIZE_X);
			CONCENTRATION_T ***markerTissue2 = (CONCENTRATION_T***) malloc(		sizeof(CONCENTRATION_T**) * vesselGraph->DOMAIN_SIZE_X);
			for (int x = 0; x < vesselGraph->DOMAIN_SIZE_X; x++) {
				markerVessel[x] = (double*) malloc(	sizeof(double) * vesselGraph->DOMAIN_SIZE_Y);
				markerVessel2[x] = (double*) malloc(	sizeof(double) * vesselGraph->DOMAIN_SIZE_Y);
				markerTissue[x] = (CONCENTRATION_T**) malloc(	sizeof(CONCENTRATION_T*) * vesselGraph->DOMAIN_SIZE_Y);
				markerTissue2[x] = (CONCENTRATION_T**) malloc(	sizeof(CONCENTRATION_T*) * vesselGraph->DOMAIN_SIZE_Y);
				for (int y = 0; y < vesselGraph->DOMAIN_SIZE_Y; y++) {
					markerVessel[x][y] = 0.;
					markerVessel2[x][y] = 0.;

					markerTissue[x][y] = (CONCENTRATION_T*) malloc(sizeof(CONCENTRATION_T) * vesselGraph->DOMAIN_SIZE_Z);
					markerTissue2[x][y] = (CONCENTRATION_T*) malloc(sizeof(CONCENTRATION_T) * vesselGraph->DOMAIN_SIZE_Z);
					for (int z = 0; z < vesselGraph->DOMAIN_SIZE_Z; z++) {
						markerTissue[x][y][z] = InterstitialConc;
						markerTissue2[x][y][z] = InterstitialConc;
					}
				}
			}
			double **before = markerVessel, **after = markerVessel2, **temp;

			for (int i = 0; i < vesselGraph->countVesselNodes; i++)
				vesselGraph->vesselNodes[i]->marker = PlasmaConc;

			fprintf(stderr, "\n\nPERFORM CONTRAST AGENT PERFUSION\n\n");

			/*double time_step = 0.1; // sec
			 double time_end  = 6000.; // sec
			 double outputResolution = 10; //sec
			 */

			char schema=IMPLICIT;
			int refinement = 1;

			float value = atof(&optarg[strcspn(optarg, "0123456789.")]);

			// SCHEMA
			if( strstr(optarg,"upwind")){
				schema = UPWIND;
			}else if( strstr(optarg,"lax")){
				schema = LAX_WENDROFF;
			}else if( strstr(optarg,"minmod")){
				schema = MINMOD;
			}else if( strstr(optarg,"superbee")){
				schema = SUPERBEE;
			}else if( strstr(optarg,"euler")){
				schema = CENTERED_EULER;
			}

			// REFINEMENT
			if( strstr(optarg,"refine")){
				refinement = (int)value;
			}else
				refinement=0;

			// MAXIMAL TIME STEP
			if( time_step==0){
				float max_dt = 1000;
				for(int v=0; v<vesselGraph->countVesselSegments; v++){
					max_dt = MIN( max_dt, vesselGraph->getVascularVolume(vesselGraph->vesselSegments[v]) / (2+refinement)
							/ ( fabs(vesselGraph->vesselSegments[v]->flow) + fabs(vesselGraph->getVascularPermeabilitySurfaceProduct(vesselGraph->vesselSegments[v])) )
					);
					//max_dt = MIN( max_dt, fabs(vesselGraph->getVascularVolume(vesselGraph->vesselSegments[v]) / (2+refinement) / vesselGraph->getVascularPermeabilitySurfaceProduct(vesselGraph->vesselSegments[v])));
				}
				for(int v=0; v<vesselGraph->countVesselNodes; v++){
					max_dt = MIN( max_dt,vesselGraph->getExtraVascularVolume(vesselGraph->vesselNodes[v])
							/  fabs(vesselGraph->getVascularPermeabilitySurfaceProduct(vesselGraph->vesselNodes[v]) ) );
					max_dt = MIN( max_dt,
							 vesselGraph->getExtraVascularVolume(vesselGraph->vesselNodes[v])/ vesselGraph->LATTICE_CONSTANT / vesselGraph->diffusionCoefficient
					);
				}
				max_dt/=2.;
				fprintf( stderr, "Time Step should satisfy: dt < %lf\n", max_dt);
				time_step = max_dt;
			}
			//time_step=0.001;

			//double time_step = 0.01;//3.01592895 / 100.; // sec
			//double time_step = 0.54/1.;//3.01592895 / 100.; // sec

			//double time_end = 180.; // sec
			//double outputResolution = 1.; //sec
			bool outputAIF = true;

			CONCENTRATION_T *b = 0, *x = 0, *dx = 0;
			CONCENTRATION_T *ds = 0;
			CONCENTRATION_T *s = 0;
			SparseMatrix<CONCENTRATION_T> *A = 0;

			FILE *fp = fopen("markerSpaTem.dat", "w");
			fclose(fp);
			double AIF = VesselGraph::AIF::ArterialInputFunction(0., 0.,
					VesselGraph::AIF::PARKER);
			if( vesselGraph->dimensions < 3){
				char filename[512];
				sprintf(filename, "markerMRI%05i", 0);
				vesselGraph->printMarkerIntensityToBinary(
						filename,
						vesselGraph,
						markerTissue,
						t,
						outputAIF,
						VesselGraph::AIF::ArterialInputFunction(time / 60., 0.,VesselGraph::AIF::PARKER));
				sprintf(filename, "markerMRI%05i.300x300x3000.eps", 0);
				vesselGraph->printMarkerIntensityMap(
						filename,
						vesselGraph,
						markerTissue,
						300,
						300,
						3000,
						outputAIF,
						VesselGraph::AIF::ArterialInputFunction(time / 60., 0.,
								VesselGraph::AIF::PARKER));
				sprintf(filename, "markerMRI%05i.60x60x60.eps", 0);
				vesselGraph->printMarkerIntensityMap(
						filename,
						vesselGraph,
						markerTissue,
						60,
						60,
						60,
						outputAIF,
						VesselGraph::AIF::ArterialInputFunction(time / 60., 0.,
								VesselGraph::AIF::PARKER));
				fprintf(stderr, "output: %i\n",
						(int) floor(time / outputResolution));
				/*fp = fopen("markerSpaTem.dat", "a+");
				for (int i = 0; i < vesselGraph->countVesselNodes; i++)
					fprintf(fp, "%f %f %e\n",
							vesselGraph->vesselNodes[i]->position[0],
							vesselGraph->vesselNodes[i]->position[1],
							vesselGraph->vesselNodes[i]->marker);
				fclose(fp);*/
			}

			// TRACK INPUT & OUTPUT
			/*FILE *fp_io = fopen("markerInOut.dat", "w");
			fclose(fp_io);

			VesselNode *vn_in = 0, *vn_out = 0;
			for (int v = 0;
					v < vesselGraph->countVesselNodes
							&& (vn_in == 0 || vn_out == 0); v++) {
				if (vesselGraph->vesselNodes[v]->pressure == MAX_PRESSURE)
					vn_in = vesselGraph->vesselNodes[v];
				if (vesselGraph->vesselNodes[v]->pressure == MIN_PRESSURE)
					vn_out = vesselGraph->vesselNodes[v];
			}*/

			//vesselGraph->writeToXML( "vasc.xml", 0, tumor);

			// TRACK INPUT & OUTPUT

			for (double time = 0.; time <= time_end; time += time_step) {
				//vesselGraph->updateMarkerVesselsAndInterstitialSpace( 0.000001,0,0,markerTissue,markerTissue2);

				//vesselGraph->updateMarkerVesselsAndInterstitialSpaceNEW2( time_step,AIF,markerTissue,markerTissue2, A,b,x);
				switch( schema){
				case IMPLICIT:
					vesselGraph->updateMarkerVesselsAndInterstitialSpaceNEW5( time_step,AIF,markerTissue,markerTissue2);
					break;

				default:
					//vesselGraph->updateMarkerVesselsExplicitRefined( time_step, AIF, x, ds, s, schema, refinement);

					if(t!=0){
						// IRREGULAR
						vesselGraph->updateMarkerVesselsAndInterstitialExplicit( time_step,AIF, dx,x,refinement, markerTissue2, cells,t);
					}else
						// REGULAR
						vesselGraph->updateMarkerVesselsAndInterstitialExplicit( time_step,AIF, dx,x,refinement, markerTissue2, cells, UpdateInterstitialCompartment, UpdatePlasmaCompartment);
					//vesselGraph->updateMarkerVesselsExplicit( time_step,AIF, dx,x,refinement);
					break;
				}
				//vesselGraph->updateMarkerVesselsExplicit(time_step, AIF, x);

				//vesselGraph->updateMarkerVesselsExplicitAntiDiffusion( time_step,AIF, x);
				//vesselGraph->updateMarkerVesselsExplicit( time_step,AIF, dx,x);

				//vesselGraph->updateMarkerVesselsAndInterstitialSpaceBrix( dt,0,0,markerTissue,markerTissue2);
				//vesselGraph->updateMarkerVessels( 0.000001,0,0);
				/*fprintf(stderr, "%lf <> floor:%lf (%i) ceil:%lf\n",
				 time/outputResolution,
				 floor(time/outputResolution), time/outputResolution == floor(time/outputResolution),
				 ceil(time/outputResolution));
				 */
				AIF = VesselGraph::AIF::ArterialInputFunction(
						(time + 1.*time_step) / 60., 0., VesselGraph::AIF::PARKER);
				float AIF_5secDelay = VesselGraph::AIF::ArterialInputFunction(
						(time + 1.*time_step - 5.) / 60., 0., VesselGraph::AIF::PARKER);


				if ((int) (time / outputResolution)
						!= (int) ((time + time_step) / outputResolution)) {
					fprintf(
							stderr,
							"\rSimulation Time: %.0lfsec = %.1lfmin (%.1lf\%%)             \b",
							time, time / 60, time / time_end * 100.);
					//fprintf(stderr, "%lf <> floor:%lf ceil:%lf\n", time/outputResolution, floor(time/outputResolution), ceil(time/outputResolution));
					CONCENTRATION_T ***temp = markerTissue;
					markerTissue = markerTissue2;
					markerTissue2 = temp;
					/*fp = fopen("marker.dat", "w+");
					if( fp){
					for (int i = 0; i < vesselGraph->countVesselNodes; i++) {
						//markerTissue[(int)vesselGraph->vesselNodes[i]->position[0]][(int)vesselGraph->vesselNodes[i]->position[1]][(int)vesselGraph->vesselNodes[i]->position[2]] = x[vesselGraph->vesselNodes[i]->index];
						fprintf(fp, "%f %f %f %e\n", 0.,
								vesselGraph->vesselNodes[i]->position[0],
								vesselGraph->vesselNodes[i]->position[1],
								vesselGraph->vesselNodes[i]->marker);
						//fprintf(stderr, "%e\n", vesselGraph->vesselNodes[i]->marker);
					}
					fclose(fp);
					}else{
						fprintf(stderr, "Couldn't open file!\n");exit(0);
					}*/
					char filename[512];
					/*sprintf( filename, "marker%i.pov", it/outputSkips);
					 vesselGraph->printToPovray( filename, 1, markerTissue);
					 //vesselGraph->printToPovray( filename, 1);
					 */
					/*sprintf( filename, "marker%i.eps", it/outputSkips);
					 vesselGraph->printToEPS( filename, 1, markerTissue);
					 */
					fprintf(stderr, "output: %i\n",
							(int) ((time + time_step) / outputResolution));


						sprintf(filename, "markerMRI%05i",
								(int) ((time + time_step) / outputResolution));
						vesselGraph->printMarkerIntensityToBinary(filename, vesselGraph, markerTissue, t, outputAIF, AIF);
					if( vesselGraph->dimensions < 3){
						sprintf(filename, "markerMRI%05i.300x300x3000.eps",
								(int) ((time + time_step) / outputResolution));
						vesselGraph->printMarkerIntensityMap(filename, vesselGraph, markerTissue,
								300, 300, 3000, outputAIF, AIF);
						sprintf(filename, "markerMRI%05i.60x60x60.eps",
								(int) ((time + time_step) / outputResolution));
						vesselGraph->printMarkerIntensityMap(filename, vesselGraph, markerTissue,
								60, 60, 60, outputAIF, AIF);
					}
					else{
						/*sprintf(filename, "markerCp%05i.pov",(int) ((time + time_step) / outputResolution));
						vesselGraph->printToPovray( filename, 1);
						sprintf(filename, "nice povray +h960 +w1280 +ua markerCp%05i.pov 2> err.out &",(int) ((time + time_step) / outputResolution));
						system(filename);

						sprintf(filename, "markerCh%05i.pov",(int) ((time + time_step) / outputResolution));
						vesselGraph->printToPovray( filename, markerTissue);
						sprintf(filename, "nice povray +h480 +w640 +ua markerCh%05i.pov 2> err.out &",(int) ((time + time_step) / outputResolution));
						system(filename);*/

						if( t!=0){
							/*int margin = 1;
							Domain<3> *domain = new Domain<3>(margin, (vesselGraph->DOMAIN_SIZE_X-margin), margin, (vesselGraph->DOMAIN_SIZE_Y-margin), margin, (vesselGraph->DOMAIN_SIZE_Z-margin));
							//vesselGraph->setType( VesselGraph::IrregularLattice);

							sprintf(filename, "markerCh_NEW_%05i.pov",(int) ((time + time_step) / outputResolution));
							t->printVoronoiDiagramToPovray3D( filename, domain);
							sprintf(filename, "nice povray -D0 +h480 +w640 +ua markerCh_NEW_%05i.pov 2> err.out &",(int) ((time + time_step) / outputResolution));
							system(filename);

							delete domain;
							exit(0);*/
						}
					}

					{
						//------ ZIP ------------
						//system("for i in `ls *.pov`; do if [ -e `echo $i|sed -e 's/pov/png/g'` ]; then zip -m pov.zip $i; fi; done");
					}

					/*if(  (int) ((time + time_step) / outputResolution / 1.)){
						fprintf( stderr, "for i in `ls *.pov`; do if [ -e ${i//pov/png} ]; then zip -m pov.zip $i; fi; done");

						//system( "for i in `ls *.pov`; do if [ -e ${i//pov/png} ]; then zip -m pov.zip $i; fi; done");
						system( "for i in \`ls \*.pov\`; do echo \$\{i//pov/png}; done");
					}*/
						/*fp = fopen("markerSpaTem.dat", "a+");
					fprintf(fp, "\n");
					for (int i = 0; i < vesselGraph->countVesselNodes; i++)
						fprintf(
								fp,
								"%f %f %f %e %e\n",
								time,
								vesselGraph->vesselNodes[i]->position[0],
								vesselGraph->vesselNodes[i]->position[1],
								vesselGraph->vesselNodes[i]->marker,
								VesselGraph::AIF::ArterialInputFunction(
										(time/*+time_step* /
												- PI
														* pow(
																vesselGraph->vesselNodes[i]->branches[0]->radius,
																2)
														/ fabs(
																vesselGraph->vesselNodes[i]->branches[0]->flow)
														* vesselGraph->vesselNodes[i]->position[0]
														* LATTICE_CONSTANT)
												/ 60., 0.,
										VesselGraph::AIF::PARKER));
					fclose(fp);

					fp_io = fopen("markerInOut.dat", "a+");
					fprintf(fp_io, "%f %e %e %e %e %e\n", time,
							vn_in->marker,
							vn_in->neighbors[0]->marker,
							vn_in->neighbors[0]->neighbors[1]->marker,
							vn_out->neighbors[0]->marker,
							vn_out->marker);
					fclose(fp_io);*/
					//sprintf( filename, "markerMRI%05i.xml", (int)((time+time_step)/outputResolution));
					//vesselGraph->writeToXML( filename, markerTissue);
				}
			}

			/*fp = fopen("end.dat", "w");
			for( int x=0; x<vesselGraph->DOMAIN_SIZE_X; x++)
				fprintf( fp, "%e %e\n", x*60., markerTissue[x][0][0]);
			fclose(fp);*/
			free(x);
			free(b);
			delete A;
		}
		break;

		}

	// CLEAN MEMORY
	delete tumor;
	delete vesselGraph;


	return 0;

}

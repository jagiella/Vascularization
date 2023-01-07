/*
 * VoronoiDiagram.hpp
 *
 *  Created on: Jan 22, 2013
 *      Author: jagiella
 */

#ifndef VORONOIDIAGRAM_HPP_
#define VORONOIDIAGRAM_HPP_

//#include <omp.h>

#include "Triangulation.hpp"
//#include "LinearAlgebra/Matrix.hpp"

#include "VesselGraph.hpp"

//class VesselNode;

#define DIM 3
class VoronoiCell : public Vertex<DIM>{
public:
	float  volume;
	float *interfaces;
	VesselNode *vn;
	float  conc;

	VoronoiCell(float x = 0, float y = 0, float z = 0, VesselNode *vesselnode = 0):Vertex<DIM>(x, y, z), interfaces(0), volume(0), vn(vesselnode), conc(0) {};
	VoronoiCell(VesselNode *vesselnode):Vertex<DIM>(vesselnode->position[0],vesselnode->position[1],vesselnode->position[2]), interfaces(0), volume(0), vn(vesselnode), conc(0) {};
};

class VoronoiDiagram : public Triangulation<DIM, VoronoiCell>{
public:
	float LATTICE_CONSTANT = 1.;

	void setVolume()
	{
#pragma omp parallel for
		for( int i=0; i<countVertices; i++){
			fprintf( stderr, "\r%.3lf%% (vertex %i)\t\t\b", 100.*(i+1.)/this->countVertices, i);
			vertices[i]->interfaces = (float*) malloc(sizeof(float)*vertices[i]->numberOfNeighbors());
			vertices[i]->volume = getVoronoiCellVolume( vertices[i], vertices[i]->interfaces);
			//getVoronoiCellInterfaces(vertices[i],vertices[i]->interfaces);
			/*fprintf(stderr, "(%.2lf,%.2lf,%.2lf)\n", vertices[i]->x(), vertices[i]->y(), vertices[i]->z());
			fprintf(stderr, "==> volume     = %lf\n", vertices[i]->volume);
			fprintf(stderr, "==> interfaces = [ ");
			for( int j=0; j<vertices[i]->numberOfNeighbors(); j++){
				fprintf(stderr, "%.2lf, ", j, vertices[i]->interfaces[j]);
			}
			fprintf(stderr, "\b\b]\n");*/
			/*const int id = omp_get_thread_num();
			if( (int)((float)countVertices*(float)id/(float)omp_get_num_threads()) == i)
				fprintf( stderr, "===============> thread %i finished\n", id);*/
		}
		//exit(0);
		fprintf(stderr, "...finished!\n");
	};
	void setInterfaces();

	typedef Simplex<DIM, VoronoiCell> SimplexType;
	typedef VoronoiCell VertexType;
	void printVoronoiDiagramToPovray3D( const char * filename, Domain<DIM>* subDomain)
	{
		if( subDomain==0)
			subDomain=domain;
		VoronoiDiagram* vd = this;
		std::fstream fs;
		fs.open( filename, std::fstream::out);
		PovrayIO::writePovrayHeader( &fs, domain->min(0), domain->min(1), domain->max(0), domain->max(1), 2*(domain->max(0)-domain->min(0)));

		//PovrayIO::beginIntersectionWithBox( &fs, vd->xMin[0], vd->xMin[1], vd->xMin[2], vd->xMax[0], vd->xMax[1], vd->xMax[2]);
		// Construction Points
	/*	for( int v=0; v<vd->countVertices; v++){
			PovrayIO::writeSphere( &fs, vd->vertices[v]->position[0], vd->vertices[v]->position[1], vd->vertices[v]->position[2], 0.04, "rgb <1,0,0>");
		}

		// Voronoi Cells & Neighborships
		for( int t=0; t<vd->countTetrahedra; t++)
		//if(!vd->tetrahedronContainsFramePoints( vd->tetrahedra[t]))
		{
			// Voronoi Cell Vertice
			double pos[DIMENSIONS];
			getCircumCenter( vd->tetrahedra[t], pos);
			PovrayIO::writeSphere( &fs, pos[0], pos[1], pos[2], 0.02, "rgb <0,1,0>");

			// Neighbors
			for( int nt=0; nt<vd->tetrahedra[t]->countNeighbors(); nt++)
			if(	//!vd->tetrahedronContainsFramePoints( vd->tetrahedra[t]->getNeighbor(nt)) &&
					vd->tetrahedra[t]->index < vd->tetrahedra[t]->getNeighbor(nt)->index){
				double npos[DIMENSIONS];
				getCircumCenter( vd->tetrahedra[t]->getNeighbor(nt), npos);
				PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], npos[0], npos[1], npos[2], 0.01, "rgb <0,0,1>");
			}
		}
	*/
		// Faces
		//fprintf( stderr, "Faces!\n");
		for( int v=0; v<vd->countVertices; v++)
		//int v = 13+27;
		if( vd->vertices[v] != NULL &&
		    vd->vertices[v]->numberOfNeighbors()>0 &&
		    subDomain->contains( vd->vertices[v]) &&
		    // only occupied cells
		    //vd->vertices[v]->getState() != FREE &&
		    // cut
			(	(// !vd->vertices[v]->refined &&
					true//vd->vertices[v]->position[2] < 0.5*(vd->xMin[2]+vd->xMax[2])
					//&&
					//vd->vertices[v]->position[2] > 0.5*(vd->xMin[2]+vd->xMax[2])-1.
				)
				/*||
				(vd->vertices[v]->refined &&
					vd->vertices[v]->position[2] < 0.5*(vd->xMin[2]+vd->xMax[2]) &&
					vd->vertices[v]->position[2] > 0.5*(vd->xMin[2]+vd->xMax[2])-1./(pow(GetAgent(vd->vertices[v])->maxCellCount, 1./3.)))*/
			)

		)
		{
			char color[100];
			/*switch (vd->vertices[v]->getState()) {
			case FREE:
				sprintf(color, "rgb <1,1,1>");
				break;
			case ACTIVE:*/

			// OLD: sprintf(color, "rgb <1,1,1>  transmit 0.5");
			sprintf(color, "rgb <0,0,1>  transmit %lf", MAX( 0, 1. - vd->vertices[v]->conc/7.) );

				//sprintf(color, "rgb <%lf,%lf,%lf> ", RAND01, RAND01, RAND01);
			/*	break;
			case NONACTIVE:
				sprintf(color, "rgb <0,1,0>");
				break;
			case COMPARTMENT:
				sprintf(color, "rgb <0,1,0>");
				break;
			}*/

			if(false)
			PovrayIO::writeSphere( &fs,
					vd->vertices[v]->pos()[0], vd->vertices[v]->pos()[1], vd->vertices[v]->pos()[2],
					0.2, "rgb <1,0,0>");

			// get tet containing Point v
			SimplexType *tet = vd->getSimplexContainingPoint( vd->vertices[v]);
			//fprintf( stderr, "Point %i -> tet %i (%i, %i, %i, %i)\n", v, tet->index, tet->vertices[0]->index, tet->vertices[1]->index, tet->vertices[2]->index, tet->vertices[3]->index);

			for( int nv=0; nv<vd->vertices[v]->numberOfNeighbors(); nv++)
			if( vd->vertices[v]->getNeighbor(nv)->getIndex() >= 0 /* no framepoints */	)
			//		&& vd->vertices[v]->index > vd->vertices[v]->neighbors[nv]->index)
			{
				//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vd->vertices[v]->neighbors[nv]->index);

				// Circulator of Tets around Point-NeighborPoint-Axis
				SimplexType *tetCirculator[100];
				int tetCirculatorLength = 0;

				// get (first) tet containing Points v & nv
				SimplexType *first = tet;
				while( !first->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))){
					// chose random neighbor tet containing point v
					int temp;
					do{
						// chose random neighbor tet
						temp = (int)(myRandE((double)first->countNeighbors()));
					}while( ! first->getNeighbor(temp)->contains( vd->vertices[v]));
					first = first->getNeighbor(temp);
				}
				tetCirculator[0] = first;
				POSITION_T posFirst[3];
				POSITION_T x[3];
				POSITION_T y[3];
				POSITION_T z[3];
				tetCirculator[0]->getCircumCenter( posFirst); x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


				// get second
				int nt=0;
				for( ; nt<first->countNeighbors() &&
					!( first->getNeighbor(nt)->contains( (VertexType*)vd->vertices[v]) &&  first->getNeighbor(nt)->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[1] = first->getNeighbor(nt);
				tetCirculatorLength = 2;

				POSITION_T pos[3];
				tetCirculator[1]->getCircumCenter( pos);
				//PovrayIO::writeCylinder( &fs, posFirst[0], posFirst[1], posFirst[2], pos[0], pos[1], pos[2], 0.01, "rgb <0,0,1>");


				// Circulate
				do{
					int nt=0;
					SimplexType *actual = tetCirculator[tetCirculatorLength-1];
					for( ; nt<actual->countNeighbors() &&
						!(  actual->getNeighbor(nt) != tetCirculator[tetCirculatorLength-2] &&
							actual->getNeighbor(nt)->contains( (VertexType*)vd->vertices[v]) &&
							actual->getNeighbor(nt)->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
					tetCirculator[tetCirculatorLength] = actual->getNeighbor(nt);
					tetCirculatorLength++;

					POSITION_T pos[3];
					tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
					POSITION_T posLast[3];
					tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

					// WRITE EDGE
					//if( pow(pos[0]-posLast[0],2) + pow(pos[1]-posLast[1],2) + pow(pos[2]-posLast[2],2) > 0.01)
					//PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], posLast[0], posLast[1], posLast[2], 0.01, "rgb <0,0,1>");

					if( first != tetCirculator[tetCirculatorLength-1]){
						x[1]=pos[0];     y[1]=pos[1];     z[1]=pos[2];
						x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

						// WRITE FACE
						POSITION_T interface = triangleInterface<DIM,POSITION_T>( pos, posLast, posFirst);
						/*if(	isnan(interface	) || isinf(interface)){
							fprintf( stderr, "ERROR in printVoronoiDiagramToPovray3D: %lf\n", interface);
							//exit(0);
						}*/
						if( interface > 0.0001 && !isnan(interface) && !isinf(interface))
							PovrayIO::writePolygon( &fs, 3, x, y, z, color);
						//PovrayIO::writePrism( &fs, 0.1, 0.1, 3, x, y, z, "rgb <1,0,1,0.7>");
					}
				}while( first != tetCirculator[tetCirculatorLength-1]);
				tetCirculatorLength--;

				/*double x[tetCirculatorLength];
				double y[tetCirculatorLength];
				double z[tetCirculatorLength];
				for( int i=0; i< tetCirculatorLength; i++){
					fprintf( stderr, "-> tet %i (%i, %i, %i, %i)\n", i, tetCirculator[i]->vertices[0]->index, tetCirculator[i]->vertices[1]->index, tetCirculator[i]->vertices[2]->index, tetCirculator[i]->vertices[3]->index);
					double pos[DIMENSIONS];
					getCircumCenter( tetCirculator[i], pos);
					x[i] = pos[0];
					y[i] = pos[1];
					z[i] = pos[2];
				}
				PovrayIO::writePolygon( &fs, tetCirculatorLength, x, y, z, "rgb <1,0,1>");*/
				//fs.close();
				//exit(0);
			}
		}
		/*for( int t=0; t<vd->countTetrahedra; t++){
			// Voronoi Cell Vertice
			double pos[DIMENSIONS];
			getCircumCenter( vd->tetrahedra[t], pos);
			PovrayIO::writeSphere( &fs, pos[0], pos[1], pos[2], 0.02, "rgb <0,1,0>");

			// Neighbors
			for( int nt=0; nt<vd->tetrahedra[t]->countNeighbors(); nt++)
			if(vd->tetrahedra[t]->index < vd->tetrahedra[t]->getNeighbor(nt)->index){
				double npos[DIMENSIONS];
				getCircumCenter( vd->tetrahedra[t]->getNeighbor(nt), npos);
				PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], npos[0], npos[1], npos[2], 0.01, "rgb <0,0,1>");
			}
		}*/
		//PovrayIO::endIntersectionWithBox( &fs);

		fs.close();
	}
};



#endif /* VORONOIDIAGRAM_HPP_ */

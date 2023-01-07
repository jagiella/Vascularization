#ifndef VESSELGRAPH_HPP_
#define VESSELGRAPH_HPP_

#include <float.h>

//#include "Geometry.hpp"
//#include "LinearAlgebra.hpp"
#include "Octree.hpp"
#include "SparseMatrix.hpp"

#include "Tumor.hpp"

#define RAND01 ((double)rand()/(double)(RAND_MAX+1.))
//#define MAX(a,b) (a>b?a:b)
//#define MIN(a,b) (a<b?a:b)

//#define UNIDIRECTIONAL_DIFFUSION
//#define UNIDIRECTIONAL_TRANSPORT

#define DIMENSIONS_		3//3

/*#if DIMENSIONS_ == 1
	// 1D SETTINGS
	#define DOMAIN_SIZE_X 60//100//601//103//301 //103
	#define DOMAIN_SIZE_Y 1//100//601//103//301 //103
	#define DOMAIN_SIZE_Z 1//50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*3)//50)

#elif DIMENSIONS_ == 2
	// 2D SETTINGS
	#define DOMAIN_SIZE_X 100//601//103//301 //103
	#define DOMAIN_SIZE_Y 100//601//103//301 //103
	#define DOMAIN_SIZE_Z 1//50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*50)
#elif DIMENSIONS_ == 3
	// 3D SETTINGS
	#define DOMAIN_SIZE_X 100//601//103//301 //103
	#define DOMAIN_SIZE_Y 100//601//103//301 //103
	#define DOMAIN_SIZE_Z 50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*50)
#endif*/

#define ROOT_DISTANCE	(BRANCH_LENGTH*50)

//#define LATTICE_CONSTANT	7.//60. //mym
#define BRANCH_LENGTH		1  //lattice sites

#define MIN_PRESSURE	0.
#define MAX_PRESSURE	10.//12.//12.

#define RADIUS_EXPONENT_ARTERIES 	3
#define RADIUS_EXPONENT_VENES 		2.7//.4

#define MIN_VESSEL_RADIUS 4.//2.5//4.
//#define MAX_VESSEL_RADIUS 200.
#define MAX_VESSEL_RADIUS (LATTICE_CONSTANT*0.4)

#define ARTERIAL_MARKER_CONC	1.
#define EXCHANGE_RATE_VESSELS_INTERSTITIAL	200.//10000.
#define MARKER_DIFFUSION_COEFFICIENT	1000.//0
//#define PERMEABILITY
 
// VesselNode Types
#define ROOT	1
#define VESSEL	2
#define SPROUT	3
#define TIP		4

// VesselSegment Types
#define UNKNOWN	0
#define ARTERIE	1
#define VENE	2
#define CAPILLARY	3

#define IMPLICIT		0
#define UPWIND			1
#define CENTERED_EULER	2
#define LAX_WENDROFF	3
#define MACCORMACK		4
//#define BEAM_WARMING	4
//#define FROMM 		5
#define MINMOD		6
#define SUPERBEE		7


/*#if DIMENSIONS_==3
	#define DISTANCE( a, b) sqrt( pow(a->position[0]-b->position[0],2) + pow(a->position[1]-b->position[1],2) + pow(a->position[2]-b->position[2],2))
#elif DIMENSIONS_==2
	#define DISTANCE( a, b) sqrt( pow(a->position[0]-b->position[0],2) + pow(a->position[1]-b->position[1],2))
#elif DIMENSIONS_==1
	#define DISTANCE( a, b) fabs( a->position[0]-b->position[0])
#endif*/

typedef double CONCENTRATION_T;

class VesselNode;
class VesselSegment;
class VoronoiDiagram;

class VesselSegment {
	public:
		char type;
		int index;
		VesselNode *vesselNodes[2];

		bool  radiusStatic;
		float radius;
		float countParallelVessels;
		float wallThickness;
		float permeability;

		float flow;
		float shear;
	public:
		VesselSegment(VesselNode *node1, VesselNode *node2);
		~VesselSegment();
		float getViscosity();
		void set( VesselNode *node1, VesselNode *node2){
			vesselNodes[0] = node1;
			vesselNodes[1] = node2;
		};
	};


class VesselNode {
		char type;
	public:
		int index;
		float position[3];
		VesselNode **neighbors;
		VesselSegment **branches;
		char countNeighbors;
		float pressure;
		float time;
		float marker;
	public:
		VesselNode(float x, float y, float z);
		~VesselNode();
		void addNeighbor(VesselNode *neighbor);
		void addNeighbor(VesselNode *neighbor, VesselSegment *vs);
		void removeNeighbor(VesselNode *neighbor);
		VesselNode *getNeighborUpstream( int n);
		VesselNode *getNeighborDownstream( int n);
		char getType();
		float getFlowFromNeighbor( int i);
		float getFlowToNeighbor( int i);
		void setType( char type);
		bool isNeighbor( VesselNode *neighbor){
			for( int n=0; n<countNeighbors; n++)
				if( neighbors[n] == neighbor)
					return true;
			return false;
		};
	};


class VesselGraph {
public:
	enum Types {RegularLattice, IrregularLattice};
private:
	Types _type;
	float _latticeConstant;
public:
	int dimensions;
	float LATTICE_CONSTANT;
	int DOMAIN_SIZE_X;
	int DOMAIN_SIZE_Y;
	int DOMAIN_SIZE_Z;
	float permeability;
	float diffusionCoefficient;

	VesselNode **vesselNodes;
	int countVesselNodes;
	int maxVesselNodes;
	
	VesselSegment **vesselSegments;
	int countVesselSegments;
	int maxVesselSegments;
	
	Octree<VesselNode*> *octree;

public:
	//VesselGraph();
	VesselGraph( int x, int y, int z, int dimensions);
	VesselGraph( const char *pFilename, Tumor *&tumor);
	~VesselGraph();
	Types type() { /*fprintf(stderr, "type();\n");*/ return _type;};
	void  setType( Types newType) { _type=newType;};
	float latticeConstant() { return _latticeConstant;};
	void addVesselNode( VesselNode *vn);
	void removeVesselNode( VesselNode *vn);
	void addVesselSegment(VesselSegment *vs);
	void addVesselSegmentAndNeighbors( VesselSegment *vs);
	void addVesselSegmentAndNeighbors( VesselSegment *vs, float radius);
	void removeVesselSegmentAndNeighbors( VesselSegment *vs);
	void removeVesselSegment( VesselSegment *vs);
	void removeVesselSegmentAndNodes( VesselSegment *vs);
	void printToVRML(const char *filename, char mode);
	void printToPovray(const char *filename, char mode);
	void printToPovray(const char *filename, CONCENTRATION_T ***marker);
	void printToPovray(const char *filename, char mode, float ***marker);
	void printToEPS(const char *filename, char mode, float ***marker);
	void writeToXML(const char *filename, float ***marker, Tumor *tumor);
	void writeStatisticsToFile(const char *filename, char mode);
	void updateFlow();
	void updateShear();
	void updatePressure( VesselNode **borderVesselNodes, float *borderPressure, int borderSize);
	void updatePressureNEW();
	void updatePressureNEW2();
	void updateRadius();
	void updateSegmentTypes();
	void updateMarkerVessels( float dt, float **markerVesselsIn, float **markerVesselsOut);
	void updateMarkerVesselsAndInterstitialSpace( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW2( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut, SparseMatrix<float> *&sA, float *&b, float *&x);
	void updateMarkerVesselsAndInterstitialSpaceNEW3( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW4( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW5( float dt, float border, CONCENTRATION_T ***markerIntSpaceIn, CONCENTRATION_T ***markerIntSpaceOut);

	void updateMarkerVesselsExplicit( float dt, float border, float *&x);
	void updateMarkerVesselsExplicitRefined( float dt, float border, float *&dx, float *&ds, float *&s, char schema, int refinement);
	void updateMarkerVesselsExplicitAntiDiffusion( float dt, float border, float *&x);
	void updateMarkerVesselsExplicit( float dt, float border, float *&dx, float *&x, int);
	void updateMarkerVesselsAndInterstitialExplicit( float dt, float border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, bool UpdateInterstitialCompartment = true, bool UpdatePlasmaCompartment = true);
	void updateMarkerVesselsAndInterstitialExplicit( float dt, float border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, VoronoiDiagram *t);


	float getVascularVolume( int &x, int &y, int &z);
	float getExtraVascularVolume( int &x, int &y, int &z);
	float getExtraVascularVolumeNoRef( int x, int y, int z);
	float getVascularSurface( VesselNode *vn);
	float getVascularSurface( VesselSegment *vs);
	float getVascularVolume( VesselNode *vn);
	float getVascularVolume( VesselSegment *vs);
	float getExtraVascularVolume( VesselNode *vn);
	float getVascularPermeabilitySurfaceProduct( VesselNode *vn);
	float getVascularPermeabilitySurfaceProduct( VesselSegment *vs);
	float getVascularPermeabilitySurfaceProduct( int &x, int &y, int &z);

	void updateMarkerVesselsAndInterstitialSpaceBrix( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut);

	void printMarkerIntensityToBinary(char *filename);
	void readMarkerIntensityFromBinary(char *filename);
	void printMarkerIntensityToBinary(char *filename, VesselGraph *vg, CONCENTRATION_T ***markerIntSpace, VoronoiDiagram *t, bool printAIF, CONCENTRATION_T markerAIF);
	void printMarkerIntensityMap(char *filename, VesselGraph *vg, CONCENTRATION_T ***markerIntSpace, double voxelSizeX, double voxelSizeY,	double voxelSizeZ, bool printAIF, float markerAIF);
	void printFlowMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printPermeabilityMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printBloodVolumeFractionMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printInterstitialSpaceVolumeFractionMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printWallShearStressMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	VesselNode *getClosestVesselNode(float *pos, float &sqrDist);

	void setSingleVessel(int length, int height, int depth, double vessel_radius = MIN_VESSEL_RADIUS);
	void setSingleRandomVessel(int length, int height, int depth);
	void setSingleVessel2(int length);
	void setSymmetricVesselGraph(int height);
	void setInitialVesselGraph( int mode);

	void setInterTipConnections(Tumor *tumor, int &countInterTipConnections);
	void setInterTipConnectionsAll( int &countInterTipConnections);
	void removeInterTipConnections( int &countInterTipConnections);

	void simplifyNetwork( float minSegmentLength, float maxSegmentLength = FLT_MAX);

	void remodel(Tumor *tumor);

	float distance( VesselNode *a, VesselNode *b);
	float distanceSqr(VesselNode *a, VesselNode *b);

	void move( VesselNode *node, float x, float y, float z){
		octree->erase(node->position[0], node->position[1], node->position[2]);
		node->position[0]=x;
		node->position[1]=y;
		node->position[2]=z;
		octree->set(node->position[0], node->position[1], node->position[2], node);
	};


	class AIF {
	public:
		// types
		static const char CONSTANT = 1;
		static const char PARKER   = 2;
		static double ArterialInputFunction( double time, double time0, const char type);
	};
};

//#include "Vascularization.tcc"

#endif /*VESSELGRAPH_HPP*/

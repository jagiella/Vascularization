/*
 * VesselGraph.cpp
 *
 *  Created on: 13.04.2010
 *      Author: jagiella
 */

#include <fstream>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <IO/tinyxml.h>
#include <float.h>


#include "LinearAlgebra.hpp"
#include "IO.hpp"

#include "VesselGraph.hpp"
#include "VoronoiDiagram.hpp"

double rgbformulae(char formulae, double x);

VesselGraph::VesselGraph( int x, int y, int z, int dimensions) : _type(RegularLattice)
{
	this->countVesselNodes = 0;
	this->countVesselSegments = 0;
	this->maxVesselNodes = 0;
	this->maxVesselSegments = 0;
	this->vesselNodes = 0;
	this->vesselSegments = 0;
	
	this->dimensions = dimensions;
	this->DOMAIN_SIZE_X = x;
	this->DOMAIN_SIZE_Y = y;
	this->DOMAIN_SIZE_Z = z;

	fprintf( stderr, "Octree Size: %ix%ix%i\n",
			(int)pow( 2, (int)ceil( log(x)/log(2))),
			(int)pow( 2, (int)ceil( log(y)/log(2))),
			(int)pow( 2, (int)ceil( log(z)/log(2))));
	octree = new Octree<VesselNode*>( 
			MAX( pow(2,(int)ceil( log(x)/log(2))),
			MAX( pow(2,(int)ceil( log(y)/log(2))),
				 pow(2,(int)ceil( log(z)/log(2)))))
					);

	diffusionCoefficient = 0;
	permeability = 0;
}


int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent)
{
	if ( !pElement ) return 0;

	TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	int i=0;
	int ival;
	double dval;
	//const char* pIndent=getIndent(indent);
	printf("\n");
	while (pAttrib)
	{
		//printf( "%s%s: value=[%s]", pIndent, pAttrib->Name(), pAttrib->Value());
		printf( "%s: value=[%s]", pAttrib->Name(), pAttrib->Value());

		if (pAttrib->QueryIntValue(&ival)==TIXML_SUCCESS)    printf( " int=%d", ival);
		if (pAttrib->QueryDoubleValue(&dval)==TIXML_SUCCESS) printf( " d=%1.1f", dval);
		printf( "\n" );
		i++;
		pAttrib=pAttrib->Next();
	}
	return i;
}

void dump_to_stdout( TiXmlNode* pParent, unsigned int indent = 0 )
{
	if ( !pParent ) return;

	TiXmlNode* pChild;
	TiXmlText* pText;
	int t = pParent->Type();
	//printf( "%s", getIndent(indent));


	switch ( t )
	{
	case TiXmlNode::TINYXML_DOCUMENT:
		printf( "Document" );
		break;

	case TiXmlNode::TINYXML_ELEMENT:
		printf( "Element [%s]", pParent->Value() );
		//int num;
		//num=dump_attribs_to_stdout(pParent->ToElement(), 0);
		/*num=dump_attribs_to_stdout(pParent->ToElement(), indent+1);
		switch(num)
		{
			case 0:  printf( " (No attributes)"); break;
			case 1:  printf( "%s1 attribute", getIndentAlt(indent)); break;
			default: printf( "%s%d attributes", getIndentAlt(indent), num); break;
		}*/
		break;

	case TiXmlNode::TINYXML_COMMENT:
		printf( "Comment: [%s]", pParent->Value());
		break;

	case TiXmlNode::TINYXML_UNKNOWN:
		printf( "Unknown" );
		break;

	case TiXmlNode::TINYXML_TEXT:
		pText = pParent->ToText();
		printf( "Text: [%s]", pText->Value() );
		break;

	case TiXmlNode::TINYXML_DECLARATION:
		printf( "Declaration" );
		break;
	default:
		break;
	}
	printf( "\n" );
	for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling())
	{
		dump_to_stdout( pChild, indent+1 );
	}
}

VesselGraph::VesselGraph( const char *pFilename, Tumor *&tumor)
{
	/*TiXmlDocument doc(pFilename);
	//TiXmlDocument doc("example.xml");
	if (!doc.LoadFile()){//std::cerr << "No such file: what ?? " << pFilename << std::endl;std::FAIL("");
		printf("Failed to load file \"%s\"\n", pFilename);
	}
	TiXmlHandle hDoc(&doc);
	TiXmlElement* pElem = hDoc.FirstChildElement().Element();
	TiXmlHandle hRoot(pElem);*/

	TiXmlDocument doc(pFilename);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		//printf("\n%s:\n", pFilename);

		TiXmlNode* pVesselGraph = doc.FirstChild("Model")->FirstChild();
		/*fprintf(stderr, "%s dim=%s x=%s y=%s z=%s\n",
				pVesselGraph->Value(),
				pVesselGraph->ToElement()->Attribute("dimensions"),
				pVesselGraph->ToElement()->Attribute("x"),
				pVesselGraph->ToElement()->Attribute("y"),
				pVesselGraph->ToElement()->Attribute("z"));*/

		this->countVesselNodes = 0;
		this->countVesselSegments = 0;
		this->maxVesselNodes = 0;
		this->maxVesselSegments = 0;
		this->vesselNodes = 0;
		this->vesselSegments = 0;

		// LATTICE CONSTANT
		if( pVesselGraph->ToElement()->Attribute("latticeConstant"))
			this->LATTICE_CONSTANT =  atof(pVesselGraph->ToElement()->Attribute("latticeConstant"));
		else
			this->LATTICE_CONSTANT = 1;

		this->dimensions = atoi(pVesselGraph->ToElement()->Attribute("dimensions"));
		this->DOMAIN_SIZE_X = atof(pVesselGraph->ToElement()->Attribute("x")) / this->LATTICE_CONSTANT;
		this->DOMAIN_SIZE_Y = atof(pVesselGraph->ToElement()->Attribute("y")) / this->LATTICE_CONSTANT;
		this->DOMAIN_SIZE_Z = atof(pVesselGraph->ToElement()->Attribute("z")) / this->LATTICE_CONSTANT;
		float latticeConstantSqr = 0;
		this->_type = RegularLattice;


		// PERMEABILITY
		if( pVesselGraph->ToElement()->Attribute("permeability"))
			this->permeability = atof(pVesselGraph->ToElement()->Attribute("permeability"));
		else
			this->permeability = 0;

		fprintf( stderr, "Octree Size: %ix%ix%i\n",
				(int)pow( 2, (int)ceil( log(DOMAIN_SIZE_X)/log(2))),
				(int)pow( 2, (int)ceil( log(DOMAIN_SIZE_Y)/log(2))),
				(int)pow( 2, (int)ceil( log(DOMAIN_SIZE_Z)/log(2))));
		octree = new Octree<VesselNode*>(
				MAX( pow(2,(int)ceil( log(DOMAIN_SIZE_X)/log(2))),
				MAX( pow(2,(int)ceil( log(DOMAIN_SIZE_Y)/log(2))),
					 pow(2,(int)ceil( log(DOMAIN_SIZE_Z)/log(2)))))
						);

		fprintf(stderr, "Read VesselGraph from \"%s\"\n", pFilename);
		for ( TiXmlNode* pChild = pVesselGraph->FirstChild(); pChild != 0; pChild = pChild->NextSibling())
		{
			//fprintf(stderr, "%s\n", pChild->Value());
			if( strcmp( pChild->Value(), "Node") == 0){
				// VESSEL NODE
				/*if( maxVesselNodes == countVesselNodes){
					maxVesselNodes += 10;
					vesselNodes = (VesselNode**)realloc(vesselNodes, maxVesselNodes*sizeof(VesselNode*));
				}
				vesselNodes[countVesselNodes] = new VesselNode(
						atoi(pChild->ToElement()->Attribute("x")),
						(dimensions>=2 ? atoi(pChild->ToElement()->Attribute("y")) : 0),
						(dimensions>=3 ? atoi(pChild->ToElement()->Attribute("z")) : 0));
				vesselNodes[countVesselNodes]->index  = atoi(pChild->ToElement()->Attribute("id"));
				vesselNodes[countVesselNodes]->setType( atoi(pChild->ToElement()->Attribute("type")));
				vesselNodes[countVesselNodes]->pressure = atof(pChild->ToElement()->Attribute("pressure"));
				vesselNodes[countVesselNodes]->marker = atof(pChild->ToElement()->Attribute("marker"));

				countVesselNodes++;*/

				//VesselNode *newVN = new VesselNode(
				//		floor(atof(pChild->ToElement()->Attribute("x")) / LATTICE_CONSTANT),
				//		floor((dimensions>=2 ? atof(pChild->ToElement()->Attribute("y")) : 0) / LATTICE_CONSTANT),
				//		floor((dimensions>=3 ? atof(pChild->ToElement()->Attribute("z")) : 0) / LATTICE_CONSTANT));
				VesselNode *newVN = new VesselNode(
						(atof(pChild->ToElement()->Attribute("x")) / LATTICE_CONSTANT),
						((dimensions>=2 ? atof(pChild->ToElement()->Attribute("y")) : 0) / LATTICE_CONSTANT),
						((dimensions>=3 ? atof(pChild->ToElement()->Attribute("z")) : 0) / LATTICE_CONSTANT));
				newVN->index  = atoi(pChild->ToElement()->Attribute("id"));
				if( pChild->ToElement()->Attribute("type"))
					newVN->setType( atoi(pChild->ToElement()->Attribute("type")));
				else
					newVN->setType( 0);
				if( pChild->ToElement()->Attribute("pressure") && pChild->ToElement()->Attribute("pressure")[0]!=0){
					/*if(pChild->ToElement()->Attribute("pressure")[0]==0){
						newVN->pressure = 0;
						//newVN->setType( ROOT);
						//fprintf(stderr, " > %lf\n", atof(pChild->ToElement()->Attribute("pressure")));
					}else
						fprintf(stderr, " > %lf\n", atof(pChild->ToElement()->Attribute("pressure")));*/
					//if( pChild->ToElement()->Attribute("pressure")[0]==0)
					//	newVN->pressure = 5;
					//else
						newVN->pressure = atof(pChild->ToElement()->Attribute("pressure"));
				}else{
					if( newVN->getType() == ROOT)
						newVN->setType( TIP);
					newVN->pressure = 5;
				}
				//if( newVN->pressure!=0)
				//fprintf(stderr, " > %i", (int)newVN->pressure);
				//newVN->pressure = 0;
				if(pChild->ToElement()->Attribute("marker"))
					newVN->marker = atof(pChild->ToElement()->Attribute("marker"));
				else
					newVN->marker = 0;

				addVesselNode(newVN);
			}else{
				// VESSEL SEGMENT
				/*if( maxVesselSegments == countVesselSegments){
					maxVesselSegments += 10;
					vesselSegments = (VesselSegment**)realloc(vesselSegments, maxVesselSegments*sizeof(VesselSegment*));
				}*/

				VesselSegment *newVS = new VesselSegment(
					vesselNodes[atoi(pChild->ToElement()->Attribute("node1"))],
					vesselNodes[atoi(pChild->ToElement()->Attribute("node2"))]);
				newVS->index  = atoi(pChild->ToElement()->Attribute("id"));
				if( pChild->ToElement()->Attribute("radius")){
					newVS->radius = atof(pChild->ToElement()->Attribute("radius"));
					newVS->radiusStatic = (atoi(pChild->ToElement()->Attribute("radiusStatic"))==1?true:false);
				}else{
					newVS->radius = atof(pChild->ToElement()->Attribute("radiusStatic"));
					newVS->radiusStatic = true;
				}
				newVS->flow   = (pChild->ToElement()->Attribute("flow") ? atof(pChild->ToElement()->Attribute("flow")) : 0);

				newVS->countParallelVessels   = ( pChild->ToElement()->Attribute("countParallelVessels") ? atof(pChild->ToElement()->Attribute("countParallelVessels")) : 1);

				newVS->permeability   = ( pChild->ToElement()->Attribute("permeability") ? atof(pChild->ToElement()->Attribute("permeability")) : 0);

				addVesselSegmentAndNeighbors(newVS);

				// CHECK FOR REGULARITY & ESTIMATE LATTICE CONSTANT
				if( _type == RegularLattice){
					if( latticeConstantSqr == 0)
						// Set Lattice Constant
						latticeConstantSqr = distanceSqr( newVS->vesselNodes[0], newVS->vesselNodes[1]);
					else if( latticeConstantSqr != distanceSqr( newVS->vesselNodes[0], newVS->vesselNodes[1])){
						// Change Type
						latticeConstantSqr = 1;
						_type = IrregularLattice;
					}
				}


				//countVesselSegments++;
			}
		}
		_latticeConstant = sqrt( latticeConstantSqr)*LATTICE_CONSTANT;
		//dump_to_stdout( &doc ); // defined later in the tutorial

		// TUMOR
		TiXmlNode* pTumor = doc.FirstChild("Model")->FirstChild("Tumor");
		if( pTumor){
			fprintf(stderr, "Read Tumor from \"%s\"\n", pFilename);
			tumor = new Tumor();
			float center[3] = {
					atof(pTumor->ToElement()->Attribute( "x")),
					atof(pTumor->ToElement()->Attribute( "y")),
					atof(pTumor->ToElement()->Attribute( "z"))
			};
			tumor->setCenter( center);
			tumor->setRadius( atof(pTumor->ToElement()->Attribute( "radius")));
			tumor->setRadiusNecroticCore( atof(pTumor->ToElement()->Attribute( "necroticCoreRadius")));
			tumor->setPermeability( atof(pTumor->ToElement()->Attribute("permeability")));
			tumor->setParallelVessels( atof(pTumor->ToElement()->Attribute("parallelVessels")));

		}
	}
	else
	{
		printf("Failed to load file \"%s\"\n", pFilename);
		exit(0);
	}

}

VesselGraph::~VesselGraph()
{
	delete octree;
	for( int i=0; i<this->countVesselSegments; i++)
		delete this->vesselSegments[i];
	free(this->vesselSegments);
	for( int i=0; i<this->countVesselNodes; i++)
		delete this->vesselNodes[i];
	free(this->vesselNodes);
}

void VesselGraph::addVesselNode( VesselNode *vn)
{
	if(countVesselNodes == maxVesselNodes){
		vesselNodes = (VesselNode**) realloc( vesselNodes, sizeof(VesselNode*) * (maxVesselNodes+10));
		maxVesselNodes+=10;
	}

	vesselNodes[countVesselNodes] = vn;
	vesselNodes[countVesselNodes]->index = countVesselNodes;
	countVesselNodes++;
	
	octree->set(vn->position[0], vn->position[1], vn->position[2], vn);
}

void VesselGraph::removeVesselNode( VesselNode *vn)
{
	countVesselNodes--;
	//fprintf( stderr, "remove Vessel Node %i\n", vn->index);
	vesselNodes[vn->index] = vesselNodes[countVesselNodes];
	vesselNodes[vn->index]->index = vn->index;

	octree->erase(vn->position[0], vn->position[1], vn->position[2]);
}

void VesselGraph::addVesselSegment( VesselSegment *vs)
{
	if(countVesselSegments == maxVesselSegments){
		vesselSegments = (VesselSegment**) realloc( vesselSegments, sizeof(VesselSegment*) * (maxVesselSegments+10));
		maxVesselSegments+=10;
	}

	vesselSegments[countVesselSegments] = vs;
	vs->index = countVesselSegments;
	countVesselSegments++;
}

void VesselGraph::addVesselSegmentAndNeighbors( VesselSegment *vs)
{
	if(countVesselSegments == maxVesselSegments){
		vesselSegments = (VesselSegment**) realloc( vesselSegments, sizeof(VesselSegment*) * (maxVesselSegments+10));
		maxVesselSegments+=10;
	}

	vs->vesselNodes[0]->addNeighbor( vs->vesselNodes[1], vs);
	vs->vesselNodes[1]->addNeighbor( vs->vesselNodes[0], vs);

	vesselSegments[countVesselSegments] = vs;
	vs->index = countVesselSegments;
	countVesselSegments++;
}

void VesselGraph::addVesselSegmentAndNeighbors( VesselSegment *vs, float radius)
{
	if(countVesselSegments == maxVesselSegments){
		vesselSegments = (VesselSegment**) realloc( vesselSegments, sizeof(VesselSegment*) * (maxVesselSegments+10));
		maxVesselSegments+=10;
	}

	vs->vesselNodes[0]->addNeighbor( vs->vesselNodes[1], vs);
	vs->vesselNodes[1]->addNeighbor( vs->vesselNodes[0], vs);

	vesselSegments[countVesselSegments] = vs;
	vs->index = countVesselSegments;
	countVesselSegments++;
	
	vs->radius = radius;
}

void VesselGraph::removeVesselSegmentAndNeighbors( VesselSegment *vs)
{
	countVesselSegments--;

	vesselSegments[vs->index] = vesselSegments[countVesselSegments];
	vesselSegments[vs->index]->index = vs->index;

	vs->vesselNodes[0]->removeNeighbor( vs->vesselNodes[1]);
	vs->vesselNodes[1]->removeNeighbor( vs->vesselNodes[0]);
}

void VesselGraph::removeVesselSegment( VesselSegment *vs)
{
	countVesselSegments--;

	vesselSegments[vs->index] = vesselSegments[countVesselSegments];
	vesselSegments[vs->index]->index = vs->index;
	
	// NEW
/*	VesselSegment *removedVS = vs;
	VesselNode *removedVN0 = removedVS->vesselNodes[0];
	VesselNode *removedVN1 = removedVS->vesselNodes[1];
	
	this->removeVesselSegment( removedVS);
	delete removedVS;

	if( removedVN0->countNeighbors == 0 && removedVN0->getType()!=ROOT){
		VesselNode *removedVN = removedVN0;
		this->removeVesselNode( removedVN);
		delete removedVN;
	}
	if( removedVN1->countNeighbors == 0 && removedVN1->getType()!=ROOT){
		VesselNode *removedVN = removedVN1;
		this->removeVesselNode( removedVN);
		delete removedVN;
	}*/
}

void VesselGraph::removeVesselSegmentAndNodes( VesselSegment *vs)
{
	countVesselSegments--;

	vesselSegments[vs->index] = vesselSegments[countVesselSegments];
	vesselSegments[vs->index]->index = vs->index;

	// NEW
	VesselSegment *removedVS = vs;
	VesselNode *removedVN0 = removedVS->vesselNodes[0];
	VesselNode *removedVN1 = removedVS->vesselNodes[1];

	removedVN0->removeNeighbor( removedVN1);
	removedVN1->removeNeighbor( removedVN0);

	//this->removeVesselSegment( removedVS);
	delete removedVS;

	if( removedVN0->countNeighbors == 0 && removedVN0->getType()!=ROOT){
		VesselNode *removedVN = removedVN0;
		this->removeVesselNode( removedVN);
		delete removedVN;
	}
	if( removedVN1->countNeighbors == 0 && removedVN1->getType()!=ROOT){
		VesselNode *removedVN = removedVN1;
		this->removeVesselNode( removedVN);
		delete removedVN;
	}
}

void VesselGraph::writeStatisticsToFile(const char *filename, char mode)
{
	int N=50;
	//bool ONLY_FUNCTIONAL_VESSELS = true;

	// Statistics
	double mean[2][N];
	double meanSquare[2][N];
	int count[2][N];

	for(int t=0;t<2;t++)
	for(int i=0; i<N; i++){
		mean[t][i]=0;
		meanSquare[t][i]=0;
		count[t][i]=0;
	}

	int arterie = 0;
	int vene = 1;
	int type;

	double dradius = 5.;
	int total = 0;
	for(int s=0; s<this->countVesselSegments; s++)
	//if(!ONLY_FUNCTIONAL_VESSELS || fabs(this->vesselSegments[s]->flow) > 1000.)
	{

		// interval
		int i = (int)(this->vesselSegments[s]->radius/dradius);

		// type
		if( this->vesselSegments[s]->type == ARTERIE)
			type = arterie;
		else
			type = vene;

		// parameter
		double parameter=0.;
		switch( mode){
		case 0:
			// pressure
			parameter = (this->vesselSegments[s]->vesselNodes[0]->pressure + this->vesselSegments[s]->vesselNodes[1]->pressure)*0.5;
			break;

		case 1:
			// flow
			parameter = fabs(this->vesselSegments[s]->flow);
			break;

		case 2:
			// velocity
			parameter = fabs(this->vesselSegments[s]->flow) / (PI * this->vesselSegments[s]->radius*this->vesselSegments[s]->radius);
			break;

		case 3:
			// shear stress
			parameter = this->vesselSegments[s]->shear;
			break;
		}

		mean      [type][i] += parameter;
		meanSquare[type][i] += parameter*parameter;
		count     [type][i]++;
		total++;
	}

	// Output
	FILE *fp = fopen( filename,"w+");
	fprintf(fp,
			"#min radius;"\
			"max radius;"\
			"type (-1=arterie, 1=vene);"\
			"parameter (type=%i);"\
			"std. deviation\n", (int)mode);

	for( int i=N-1;i>=0;i--)
		if(count[0][i]!=0){
			// Arteries
			fprintf(fp, "%lf %lf %e %e\n",
				-dradius*(i+0.5),
				-dradius*0.5,
				mean[0][i]/count[0][i],
				sqrt( meanSquare[0][i]/count[0][i] - mean[0][i]/count[0][i] * mean[0][i]/count[0][i]));
		}

	for( int i=0;i<N;i++)
		if(count[1][i]!=0){
			// Venes
			fprintf(fp, "%lf %lf %e %e\n",
				dradius*(i+0.5),
				dradius*0.5,
				mean[1][i]/count[1][i],
				sqrt( meanSquare[1][i]/count[1][i] - mean[1][i]/count[1][i] * mean[1][i]/count[1][i]));
		}

	fclose(fp);
}

void crossProduct( float *a, float *b, float *c)
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

float dotProduct( float *a, float *b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

int printSphere( FILE* fp, float ax,float ay, float az, float radius, float R, float G, float B)
{
	float a[3] = {ax,ay,az};
	fprintf( fp, "Transform{\n"\
			"translation %lf %lf %lf\n"\
			"children Shape{ \n"\
				"appearance Appearance { material Material{ diffuseColor %lf %lf %lf } } \n"\
             			"geometry Sphere {\n"\
               				"radius %lf\n"\
             			"}\n"\
			"}\n"\
		"}\n", a[0], a[1], a[2], R,G,B, radius);


}


int printCylinder( FILE* fp, float ax,float ay, float az, float bx,float by, float bz, float radius, float R, float G, float B)
{
	float a[3] = {ax,ay,az}, b[3] = {bx,by,bz};
	float c[3];
	float zaxis[3] = {0,1,0};
	float v[3] = {b[0]-a[0],b[1]-a[1],b[2]-a[2]};
	float length = sqrt( dotProduct(v,v));
	crossProduct(zaxis,v,c);
	float cabs=sqrt(c[0]*c[0] +c[1]*c[1] +c[2]*c[2]);
	c[0]/=cabs;
	c[1]/=cabs;
	c[2]/=cabs;
	float angle = acos( dotProduct(zaxis,v) / sqrt(dotProduct(zaxis,zaxis)) / sqrt(dotProduct(v,v)));
	float t[3] = {0.5*(b[0]+a[0]),0.5*(b[1]+a[1]),0.5*(b[2]+a[2])};

	if( isnan(angle) || 0==cabs || isnan(c[0])){
	//	fprintf( stderr, "cabs=%e, c={%e, %e, %e}\n", cabs, c[0], c[1], c[2]);
	}else
	fprintf(fp,
		"Transform{\n"\
			"rotation %lf %lf %lf %lf\n"\
			"translation %lf %lf %lf\n"\
			"children Shape{ \n"\
				"appearance Appearance { material Material{ diffuseColor %lf %lf %lf } } \n"\
             			"geometry Cylinder {\n"\
               				"radius %lf\n"\
               				"height %lf\n"\
             			"}\n"\
			"}\n"\
		"}\n",
		c[0],c[1],c[2],angle,
		t[0],t[1],t[2],
		R,G,B,
		radius,length);
}


void VesselGraph::printToVRML(const char *filename, char mode)
{
	fprintf(stderr, "OUTPUT to VRML\n");
	//std::fstream fs;
	//fs.open( filename, std::fstream::out);
	FILE * fp = fopen( filename, "w+");
	fprintf(fp, "#VRML V2.0 utf8\n");
	fprintf(fp, "Viewpoint {\n"\
					"position    %lf %lf %lf\n"\
					"orientation %lf %lf %lf 0\n"\
				"}\n",
				DOMAIN_SIZE_X/2., DOMAIN_SIZE_Y/2., DOMAIN_SIZE_Z/2.+DOMAIN_SIZE_X/2.+DOMAIN_SIZE_Y/2.,
				0., 0., -1.);
	//PovrayIO::writePovrayHeader( &fs, 0, 0, DOMAIN_SIZE_X, DOMAIN_SIZE_Y);

	fprintf(fp, "Background {\n"\
					"skyColor [\n"\
						"1.0 1.0 1.0,\n"\
					"]\n"\
					"skyAngle [ 1.309, 1.571 ]\n"\
					"groundColor [\n"\
						"1.0 1.0 1.0,\n"\
						"0.9 0.9 0.9,\n"\
						"0.7 0.7 0.7\n"\
					" ]\n"\
					"groundAngle [ 1.309, 1.571 ]\n"\
				"}\n");

	float R,G,B;

	// ESTIMATE MAXIMUM
	float minf = 1e-3;
	float maxf = 1e3;
	for( int s=0; s<this->countVesselSegments; s++) {
		float value;
		switch( mode) {
			case 3:
				value = fabs(this->vesselSegments[s]->flow / ( PI * vesselSegments[s]->radius* vesselSegments[s]->radius));
				break;
			case 2:
				value = fabs(this->vesselSegments[s]->flow);
				break;
			case 1:
				value = (this->vesselSegments[s]->vesselNodes[0]->pressure + this->vesselSegments[s]->vesselNodes[1]->pressure)/2.;
				break;
		}

		if( minf > value)
			minf = value;
		if( maxf < value)
			maxf = value;
	}


	for( int s=0; s<this->countVesselSegments; s++){

		char color[512],color0[512],color1[512];
		switch( mode){
		case 3:{
			// Velocity
			float minf = 1e-3;
			float maxf = 1e3;
			float f = 1./maxf*fabs(this->vesselSegments[s]->flow / ( PI * vesselSegments[s]->radius* vesselSegments[s]->radius));
			 R = rgbformulae(22, f);G = rgbformulae(13, f);B=rgbformulae(-31, f);
		}break;
		case 2:{
			// Flow
			float minf = 1e-3;
			float maxf = 1e3;
			float f = 1./maxf*fabs(this->vesselSegments[s]->flow);
			 R = rgbformulae(22, f);G = rgbformulae(13, f);B=rgbformulae(-31, f);
		}break;
		case 1:{
			// Marker Concentration
			float v0 = MIN(this->vesselSegments[s]->vesselNodes[0]->marker, 1);
			float v1 = MIN(this->vesselSegments[s]->vesselNodes[1]->marker, 1);
			 R = 1.; G=1.-v0; B=1.-v0;
		}break;
		case 0:
			// Pressure
			R = 1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE);
			G = 0;
			B = (MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE);
			break;
		}

		printCylinder( fp,
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[1],
				this->vesselSegments[s]->vesselNodes[0]->position[2],
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise1[1],
				this->vesselSegments[s]->vesselNodes[1]->position[2],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				R,G,B);

		if( !isnan(R)){
		if( this->vesselSegments[s]->vesselNodes[0]->position[0] <  this->vesselSegments[s]->vesselNodes[1]->position[0])
		printSphere( fp,
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[1],
				this->vesselSegments[s]->vesselNodes[0]->position[2],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				R,G,B);
		else
		printSphere( fp,
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise0[1],
				this->vesselSegments[s]->vesselNodes[1]->position[2],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				R,G,B);
		}

	}


	fclose(fp);
}


void VesselGraph::printToPovray(const char *filename, char mode){
	fprintf(stderr, "OUTPUT (%i)\n", mode);
	std::fstream fs;
	fs.open( filename, std::fstream::out);

	float max[3] = {DOMAIN_SIZE_X,DOMAIN_SIZE_Y,DOMAIN_SIZE_Z};
	float min[3] = {0,0,0};


	if( mode==5){
		min[0] = DOMAIN_SIZE_X;
		min[1] = DOMAIN_SIZE_Y;
		min[2] = DOMAIN_SIZE_Z;
		max[0] = max[1] = max[2] = 0;

		for( int n=0; n<this->countVesselNodes; n++)
			for(int d=0;d<3;d++){
				if( max[d] < this->vesselNodes[n]->position[d])
					max[d] = this->vesselNodes[n]->position[d];
				if( min[d] > this->vesselNodes[n]->position[d])
					min[d] = this->vesselNodes[n]->position[d];
		}
	}

	PovrayIO::writePovrayHeader( &fs, min[0], min[1], max[0], max[1]);

	// ESTIMATE MAXIMUM
	float minf = FLT_MAX;//1e-3;
	float maxf = 0;//1e3;
	float meanf= 0;
	for( int s=0; s<this->countVesselSegments; s++) {
		float value;
		switch( mode) {
			case 3:
				value = fabs(this->vesselSegments[s]->flow / ( PI * vesselSegments[s]->radius* vesselSegments[s]->radius));
				break;
			case 2:
				value = fabs(this->vesselSegments[s]->flow);
				break;
			case 1:
				value = (this->vesselSegments[s]->vesselNodes[0]->pressure + this->vesselSegments[s]->vesselNodes[1]->pressure)/2.;
				break;
		}
		meanf += value;
		if( minf > value)
			minf = value;
		if( maxf < value)
			maxf = value;
	}
	meanf /= countVesselSegments;
	fprintf( stderr, "min-max (mean): %e - %e (%e)\n", minf, maxf, meanf);

	for( int s=0; s<this->countVesselSegments; s++){
		//float noise0[2] = {RAND01,RAND01};
		//float noise1[2] = {RAND01,RAND01};
		
		char color[512],color0[512],color1[512];
		/*sprintf( color0, "rgb<%lf ,0 ,%lf>",
				1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
				(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
		sprintf( color1, "rgb<%lf ,0 ,%lf>", 
				1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
				(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
		sprintf( color, "rgb<%lf ,0 ,%lf>",
				1.-(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
				(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
		 */
		switch( mode){
		case 3:{
			// VELOCITY

			minf = 50;//1e-3;
			maxf = 750;//1e3;

			float f = fabs(this->vesselSegments[s]->flow / ( PI * vesselSegments[s]->radius* vesselSegments[s]->radius));
			f = (log(f) - log(minf)) / (log(maxf) - log(minf));
			 if( isnan(f) || f<0)
				 f=0;
			 if( isinf(f) || f>1)
				 f=1;
			float R = rgbformulae(22, f),G = rgbformulae(13, f),B=rgbformulae(-31, f);
			sprintf( color0, "rgb<%lf, %lf ,%lf>", R, G, B);
			sprintf( color1, "rgb<%lf, %lf ,%lf>", R, G, B);
			sprintf( color,  "rgb<%lf, %lf ,%lf>", R, G, B);
		}break;
		case 2:{
			// FLOW

			//float minf = 1e-3;
			//float maxf = 1e3;
			float f = fabs(this->vesselSegments[s]->flow);

			 minf = 1e-1;
			 maxf = 1e5;

			 //f = (log(f) - log(minf)) / (log(maxf) - log(minf));
			 f = (log(f) - log(minf)) / (log(maxf) - log(minf));
			 if( isnan(f) || f<0)
				 f=0;
			 if( isinf(f) || f>1)
				 f=1;

			float R = rgbformulae(22, f),G = rgbformulae(13, f),B=rgbformulae(-31, f);
			sprintf( color0, "rgb<%lf, %lf ,%lf>", R, G, B);
			sprintf( color1, "rgb<%lf, %lf ,%lf>", R, G, B);
			sprintf( color,  "rgb<%lf, %lf ,%lf>", R, G, B);

		}break;
		case 1:{
			// Marker Concentration
			float v0 = this->vesselSegments[s]->vesselNodes[0]->marker/7.;
			float v1 = this->vesselSegments[s]->vesselNodes[1]->marker/7.;
			//float v0 = MIN(this->vesselSegments[s]->vesselNodes[0]->marker/6., 1);
			//float v1 = MIN(this->vesselSegments[s]->vesselNodes[1]->marker/6., 1);
			//if( this->vesselSegments[s]->vesselNodes[1]->marker!=0)
			//	fprintf(stderr, "node %i has a marker conc. of %e\n", this->vesselSegments[s]->vesselNodes[1]->index, this->vesselSegments[s]->vesselNodes[1]->marker);

			sprintf( color0, "rgb<%lf, %lf ,%lf>", 1., 1.-v0, 1.-v0);
			sprintf( color1, "rgb<%lf, %lf ,%lf>", 1., 1.-v1, 1.-v1);
			sprintf( color,  "rgb<%lf, %lf ,%lf>", 1., 1.-(v0+v1)/2., 1.-(v0+v1)/2.);

			/*sprintf( color0, "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[0]->marker, 1-this->vesselSegments[s]->vesselNodes[0]->marker);
			sprintf( color1, "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[1]->marker, 1-this->vesselSegments[s]->vesselNodes[1]->marker);
			sprintf( color,  "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5, 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5);
			*/
		}break;
		case 0:{
			// Pressure
			sprintf( color0, "rgb<%lf ,0 ,%lf>",
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			sprintf( color1, "rgb<%lf ,0 ,%lf>", 
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			sprintf( color, "rgb<%lf ,0 ,%lf>",
					1.-(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));		
			break;
		}
		case 4:{
			// Volume Fraction (Plasma)
			float f = MIN( 1, this->getVascularVolume( vesselSegments[s]->vesselNodes[0]) / LATTICE_CONSTANT/LATTICE_CONSTANT/LATTICE_CONSTANT) ;
			float R = rgbformulae(22, f),G = rgbformulae(13, f),B=rgbformulae(-31, f);

			sprintf( color0, "rgb<%f,%f,%f>",R,G,B);
			f = MIN( 1, this->getVascularVolume( vesselSegments[s]->vesselNodes[1]) / LATTICE_CONSTANT/LATTICE_CONSTANT/LATTICE_CONSTANT) ;
			R = rgbformulae(22, f),G = rgbformulae(13, f),B=rgbformulae(-31, f);
			sprintf( color1, "rgb<%f,%f,%f>",R,G,B);
			sprintf( color, "rgbf<1 ,1 ,1,0>");
			break;
		}

		case 5:
			{
			float R = 1, G = 0, B = 0;
			sprintf( color,  "rgb<%.2f,%.2f,%.2f>",R,G,B);
			sprintf( color0, "rgb<%.2f,%.2f,%.2f>",R,G,B);
			sprintf( color1, "rgb<%.2f,%.2f,%.2f>",R,G,B);
			}
			break;
		}
		//sprintf( color, "rgb<1 ,0 ,0>");
		/*if( this->vesselSegments[s]->vesselNodes[0]->getType() == TIP) sprintf( color0, "rgb<0 ,0 ,1>");
		else if( this->vesselSegments[s]->vesselNodes[0]->getType() == ROOT) sprintf( color0, "rgb<0 ,1 ,0>");
		else sprintf( color0, "rgb<1 ,0 ,0>");
		if( this->vesselSegments[s]->vesselNodes[1]->getType() == TIP) sprintf( color1, "rgb<0 ,0 ,1>");
		else if( this->vesselSegments[s]->vesselNodes[1]->getType() == ROOT) sprintf( color0, "rgb<0 ,1 ,0>");
		else sprintf( color1, "rgb<1 ,0 ,0>");*/
		
		PovrayIO::writeCylinder( &fs,
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[1],
				this->vesselSegments[s]->vesselNodes[0]->position[2],
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise1[1],
				this->vesselSegments[s]->vesselNodes[1]->position[2],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				color);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure+this->vesselSegments[s]->vesselNodes[1]->pressure!=0? "Red":"Blue"));
		PovrayIO::writeSphere( &fs, 
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0], 
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[2],//+noise0[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[0]->getType() == TIP||this->vesselSegments[s]->vesselNodes[0]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color0);
				//(this->vesselSegments[s]->vesselNodes[0]->pressure!=0?"Red":"Blue"));
		PovrayIO::writeSphere( &fs, 
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise1[0], 
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[2],//+noise1[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[1]->getType() == TIP||this->vesselSegments[s]->vesselNodes[1]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color1);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure!=0?"Red":"Blue"));

	}

	//for( int s=0; s<this->countVesselNodes; s++)
	//	PovrayIO::writeSphere( &fs, this->vesselNodes[s]->position[0], this->vesselNodes[s]->position[1], 1, "Red");

	fs.close();
}

void VesselGraph::printToPovray(const char *filename, CONCENTRATION_T ***marker){
	std::fstream fs;
	fs.open( filename, std::fstream::out);

	PovrayIO::writePovrayHeader( &fs, 0, 0, DOMAIN_SIZE_X, DOMAIN_SIZE_Y);
	for( int x=0; x<DOMAIN_SIZE_X; x++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int z=0; z<DOMAIN_SIZE_Z; z++)
			{
				float f = marker[x][y][z]/7.;// * this->getExtraVascularVolume(x,y,z) / LATTICE_CONSTANT/ LATTICE_CONSTANT/ LATTICE_CONSTANT;
				//float f = MIN( marker[x][y][z]/6., 1);// * this->getExtraVascularVolume(x,y,z) / LATTICE_CONSTANT/ LATTICE_CONSTANT/ LATTICE_CONSTANT;
				float e = 0.001;
				fs.precision(5);
				fs << std::fixed;
				fs << "box{<"<<x+e<<","<<y+e<<","<<z+e<<">, <"<<x+1-e<<","<<y+1-e<<","<<z+1-e<<"> pigment{ rgbt 1 } hollow interior{ media{ absorption<"<< f <<","<< f <<",0.>}}}\n";
				//fs << "box{<"<<x+e<<","<<y+e<<","<<z+e<<">, <"<<x+1-e<<","<<y+1-e<<","<<z+1-e<<"> pigment{ color rgbt< 0,0,1,"<< 1-f <<">}} \n";
			}
}

void VesselGraph::printToPovray(const char *filename, char mode, float ***marker)
{
	std::fstream fs;
	fs.open( filename, std::fstream::out);

	PovrayIO::writePovrayHeader( &fs, 0, 0, DOMAIN_SIZE_X, DOMAIN_SIZE_Y);
	fs << "global_settings { max_trace_level 20 }\n";

	for( int s=0; s<this->countVesselSegments; s++){
		//float noise0[2] = {RAND01,RAND01};
		//float noise1[2] = {RAND01,RAND01};

		char color[512],color0[512],color1[512];
		switch( mode){
		case 1:{
			// Marker Concentration
			float v0 = MAX(this->vesselSegments[s]->vesselNodes[0]->marker, 1);
			float v1 = MAX(this->vesselSegments[s]->vesselNodes[1]->marker, 1);

			sprintf( color0, "rgb<%lf, %lf ,%lf>", v0, 1-v0, 1-v0);
			sprintf( color1, "rgb<%lf, %lf ,%lf>", v1, 1-v1, 1-v1);
			sprintf( color,  "rgb<%lf, %lf ,%lf>", (v0+v1)/2., 1-(v0+v1)/2., 1-(v0+v1)/2.);

			/*sprintf( color0, "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[0]->marker, 1-this->vesselSegments[s]->vesselNodes[0]->marker);
			sprintf( color1, "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[1]->marker, 1-this->vesselSegments[s]->vesselNodes[1]->marker);
			sprintf( color,  "rgb<1, %lf ,%lf>", 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5, 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5);
			*/
		}break;
		case 0:
			// Pressure
			sprintf( color0, "rgb<%lf ,0 ,%lf>",
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			sprintf( color1, "rgb<%lf ,0 ,%lf>",
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			sprintf( color, "rgb<%lf ,0 ,%lf>",
					1.-(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			break;
		}

		PovrayIO::writeCylinder( &fs,
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[1],
				this->vesselSegments[s]->vesselNodes[0]->position[2],
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise1[1],
				this->vesselSegments[s]->vesselNodes[1]->position[2],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				color);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure+this->vesselSegments[s]->vesselNodes[1]->pressure!=0? "Red":"Blue"));
		PovrayIO::writeSphere( &fs,
				this->vesselSegments[s]->vesselNodes[0]->position[0],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1],//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[2],//+noise0[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[0]->getType() == TIP||this->vesselSegments[s]->vesselNodes[0]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color0);
				//(this->vesselSegments[s]->vesselNodes[0]->pressure!=0?"Red":"Blue"));
		PovrayIO::writeSphere( &fs,
				this->vesselSegments[s]->vesselNodes[1]->position[0],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1],//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[2],//+noise1[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[1]->getType() == TIP||this->vesselSegments[s]->vesselNodes[1]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color1);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure!=0?"Red":"Blue"));

	}

	if(DOMAIN_SIZE_X==1)
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				/*char color[512];
				sprintf( color, "rgb<1, %lf ,%lf>", 1-marker[x][y][z], 1-marker[x][y][z]);
				fs << "polygon {";
				fs << "5,";
				fs << "<"<<x-0.5<<", "<<y-0.5<<">, <"<<x-0.5<<", "<<y+0.5<<">, <"<<x+0.5<<", "<<y+0.5<<">, <"<<x+0.5<<", "<<y-0.5<<">, <"<<x-0.5<<", "<<y-0.5<<">";
				fs << " pigment{color "<<color<<"} finish{ambient 1}}";
				 */
				char color[512];
				//sprintf( color, "rgbt<1,%lf,%lf,%lf>", 1-marker[x][y][z], 1-marker[x][y][z], 1-marker[x][y][z]);
				sprintf( color, "rgb<1, %lf ,%lf,%lf>", 1-marker[x][y][z], 1-marker[x][y][z], 0.7);
				PovrayIO::writeCube( &fs, x-0.5, y-0.5, z-0.5, x+0.5, y+0.5, z+0.5, color);
			}

	fs.close();
}


void VesselGraph::printToEPS(const char *filename, char mode, float ***marker)
{
	std::fstream fs;
	fs.open( filename, std::fstream::out);
	float shift=0;//0.5;

	EPS::PSwriteHeader( &fs, 0, DOMAIN_SIZE_X, 0, DOMAIN_SIZE_Y);

	if( mode==1)
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				char color[512];
				EPS::PScolor( color, 1, 1-marker[x][y][z], 1-marker[x][y][z]);
				int N = 5;
				double X[] = {x, x,   x+1, x+1, x};
				double Y[] = {y, y+1, y+1, y, y};

				EPS::PSfillPolygon( &fs,
						X, Y, N,
						0., color);
			}

	for( int s=0; s<this->countVesselSegments; s++){
		//float noise0[2] = {RAND01,RAND01};
		//float noise1[2] = {RAND01,RAND01};

		char color[512],color0[512],color1[512];
		switch( mode){
		case 1:
			// Marker Concentration
			EPS::PScolor( color0, 1, 1-this->vesselSegments[s]->vesselNodes[0]->marker, 1-this->vesselSegments[s]->vesselNodes[0]->marker);
			EPS::PScolor( color1, 1, 1-this->vesselSegments[s]->vesselNodes[1]->marker, 1-this->vesselSegments[s]->vesselNodes[1]->marker);
			EPS::PScolor( color,  1, 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5, 1-this->vesselSegments[s]->vesselNodes[0]->marker*0.5-this->vesselSegments[s]->vesselNodes[1]->marker*0.5);
			break;
		case 0:
			// Pressure
			EPS::PScolor( color0,
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					0, (MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[0]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			EPS::PScolor( color1,
					1.-(MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					0, (MAX_PRESSURE - this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			EPS::PScolor( color,
					1.-(MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
					0, (MAX_PRESSURE - 0.5*this->vesselSegments[s]->vesselNodes[0]->pressure - 0.5*this->vesselSegments[s]->vesselNodes[1]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));
			break;
		}

		if( this->vesselSegments[s]->countParallelVessels > 1)
		for( float i=0; i<this->vesselSegments[s]->countParallelVessels; i++)
		{
			float n=this->vesselSegments[s]->countParallelVessels;

			// starting point
			float x0 = this->vesselSegments[s]->vesselNodes[0]->position[0]+shift;
			float y0 = this->vesselSegments[s]->vesselNodes[0]->position[1]+shift;

			// end point
			float x1 = this->vesselSegments[s]->vesselNodes[1]->position[0]+shift;
			float y1 = this->vesselSegments[s]->vesselNodes[1]->position[1]+shift;

			// normal vector
			float nx = y1 - y0;
			float ny = x0 - x1;

			// middle point
			float mx = 0.5 * (x0 + x1);
			float my = 0.5 * (y0 + y1);


			fs << "newpath\n";
			fs << color  << "\n";
			fs << 2.*this->vesselSegments[s]->radius/LATTICE_CONSTANT << " setlinewidth\n";
			fs << x0 << " " << y0 << " moveto\n";
			fs << mx + nx*(i/(n-1) - shift) << " " << my + ny*(i/(n-1) - shift) << " ";
			fs << mx + nx*(i/(n-1) - shift) << " " << my + ny*(i/(n-1) - shift) << " ";
			fs << x1 << " " << y1 << " curveto\n";
			fs << "stroke\n";
		}
		else
			EPS::PSwriteLine( &fs,
				this->vesselSegments[s]->vesselNodes[0]->position[0]+shift,//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1]+shift,//+noise0[1],
				//this->vesselSegments[s]->vesselNodes[0]->position[2],
				this->vesselSegments[s]->vesselNodes[1]->position[0]+shift,//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1]+shift,//+noise1[1],
				//this->vesselSegments[s]->vesselNodes[1]->position[2],
				2*this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				color);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure+this->vesselSegments[s]->vesselNodes[1]->pressure!=0? "Red":"Blue"));
		EPS::PSwriteSphere( &fs,
				this->vesselSegments[s]->vesselNodes[0]->position[0]+shift,//+noise0[0],
				this->vesselSegments[s]->vesselNodes[0]->position[1]+shift,//+noise0[0],
				//this->vesselSegments[s]->vesselNodes[0]->position[2],//+noise0[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[0]->getType() == TIP||this->vesselSegments[s]->vesselNodes[0]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color0);
				//(this->vesselSegments[s]->vesselNodes[0]->pressure!=0?"Red":"Blue"));
		EPS::PSwriteSphere( &fs,
				this->vesselSegments[s]->vesselNodes[1]->position[0]+shift,//+noise1[0],
				this->vesselSegments[s]->vesselNodes[1]->position[1]+shift,//+noise1[0],
				//this->vesselSegments[s]->vesselNodes[1]->position[2],//+noise1[1],
				this->vesselSegments[s]->radius/LATTICE_CONSTANT,
				//(this->vesselSegments[s]->vesselNodes[1]->getType() == TIP||this->vesselSegments[s]->vesselNodes[1]->getType() == ROOT?.3:this->vesselSegments[s]->radius/10.),
				color1);
				//(this->vesselSegments[s]->vesselNodes[1]->pressure!=0?"Red":"Blue"));

	}

	/*for( int vn=0; vn<this->countVesselNodes; vn++){

		char color[512];
		float radius=0;
		if( this->vesselNodes[vn]->getType() == ROOT)
				radius = 32./LATTICE_CONSTANT;
		if( this->vesselNodes[vn]->getType() == VESSEL)
				radius = 16./LATTICE_CONSTANT;
		if( this->vesselNodes[vn]->getType() == TIP)
				radius = 8./LATTICE_CONSTANT;

		EPS::PScolor( color,
				1.-(MAX_PRESSURE - 0.5*this->vesselNodes[vn]->pressure - 0.5*this->vesselNodes[vn]->pressure)/(MAX_PRESSURE - MIN_PRESSURE),
				0,
				(MAX_PRESSURE - 0.5*this->vesselNodes[vn]->pressure - 0.5*this->vesselNodes[vn]->pressure)/(MAX_PRESSURE - MIN_PRESSURE));

		EPS::PSwriteSphere( &fs,
				this->vesselNodes[vn]->position[0]+0.5,//+noise1[0],
				this->vesselNodes[vn]->position[1]+0.5,//+noise1[0],
				radius,
				color);

	}*/


	fs.close();
}


void VesselGraph::writeToXML(const char *filename, float ***marker, Tumor *tumor)
{
	FILE *fp = fopen( filename, "w+");

	fprintf( fp, "<Model>\n");


	if( marker){
	// SquareLattice
	fprintf( fp, "  <SquareLattice LatticeConstant='%e' dimensions='%i'>\n", LATTICE_CONSTANT, dimensions);
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				// open element
				fprintf(fp, "    <Node");

				// position
				fprintf(fp, " x='%i'", x);
				if(dimensions>=2)
					fprintf(fp, " y='%i'", y);
				if(dimensions>=3)
					fprintf(fp, " z='%i'", z);

				// concentration
				fprintf(fp, " marker='%e'", marker[x][y][z]);

				// close element
				fprintf(fp, "/>\n");
			}
	fprintf( fp, "  </SquareLattice>\n");
	}

	// VesselGraph
	fprintf( fp, "  <VesselGraph dimensions='%d' x='%d' y='%d' z='%d' latticeConstant ='%e' permeability='%e'>\n",
			this->LATTICE_CONSTANT,
			dimensions,
			(int)(DOMAIN_SIZE_X*LATTICE_CONSTANT),
			(int)(DOMAIN_SIZE_Y*LATTICE_CONSTANT),
			(int)(DOMAIN_SIZE_Z*LATTICE_CONSTANT),
			permeability);
	for( int vn=0; vn<this->countVesselNodes; vn++){
		// open element
		fprintf(fp, "    <Node");

		// position
		char label[3] = {'x','y','z'};
		for(int d=0; d<dimensions; d++)
			fprintf(fp, " %c='%e'", label[d], this->vesselNodes[vn]->position[d] * LATTICE_CONSTANT);

		// index
		fprintf(fp, " id='%i'", this->vesselNodes[vn]->index);

		// type
		fprintf(fp, " type='%i'", this->vesselNodes[vn]->getType());

		// pressure
		fprintf(fp, " pressure='%e'", this->vesselNodes[vn]->pressure);

		// concentration
		fprintf(fp, " marker='%e'", this->vesselNodes[vn]->marker);

		// close element
		fprintf(fp, "/>\n");
	}


	for( int s=0; s<this->countVesselSegments; s++){
		// open element
		fprintf(fp, "    <Segment");

		// index
		fprintf(fp, " id='%i'", this->vesselSegments[s]->index);

		// nodes
		fprintf(fp, " node1='%i'", this->vesselSegments[s]->vesselNodes[0]->index);
		fprintf(fp, " node2='%i'", this->vesselSegments[s]->vesselNodes[1]->index);

		// radius
		fprintf(fp, " radius='%e'", this->vesselSegments[s]->radius);
		fprintf(fp, " radiusStatic='%d'", (this->vesselSegments[s]->radiusStatic ? 1 : 0));

		fprintf(fp, " countParallelVessels='%e'", this->vesselSegments[s]->countParallelVessels);

		// permeability
		fprintf(fp, " permeability='%e'", this->vesselSegments[s]->permeability);

		// flow
		fprintf(fp, " flow='%e'", this->vesselSegments[s]->flow);

		// close element
		fprintf(fp, "/>\n");
	}
	fprintf( fp, "  </VesselGraph>\n");

	if( tumor){
		fprintf( fp, "  <Tumor x='%e' y='%e' z='%e' radius='%e' necroticCoreRadius='%e' permeability='%e' parallelVessels='%e'/>\n",
				tumor->getCenter(0), tumor->getCenter(1), tumor->getCenter(2),
				tumor->getRadius(),
				tumor->getNecroticCoreRadius(),
				tumor->getPermeability(),
				tumor->getParallelVessels());
	}

	fprintf( fp, "  </Model>\n");


	fclose(fp);
}


float VesselSegment::getViscosity()
{
	//return 4e-6 * (2.+exp(-this->radius/4.)*15.); // kPa * s

	float visc_45 = 6 * exp(-0.085*this->radius*2.) + 3.1 - 2.44 * exp(-0.06*pow(this->radius*2., 0.645));
	float visc_vivo = (1 + (visc_45-1)*pow(radius/(radius-0.55),2))*pow(radius/(radius-0.55),2);
	return 4e-6 * visc_45; // kPa * s

	/*return
	- 0.93   *pow(log(this->radius), 7) 
	+ 12.42  *pow(log(this->radius), 6) 
	- 68.79  *pow(log(this->radius), 5)
	+ 205.42 *pow(log(this->radius), 4)
	- 358.38 *pow(log(this->radius), 3) 
	+ 367.3  *pow(log(this->radius), 2)
	- 204.1  *log(this->radius) 
	+ 48.4;*/
}


void VesselGraph::updateFlow()
{
#pragma omp parallel for
	for( int i=0; i<this->countVesselSegments; i++){
		// flow: n0 -> n1
		this->vesselSegments[i]->flow
			= this->vesselSegments[i]->countParallelVessels *
			  PI/8. * pow( this->vesselSegments[i]->radius, 4)
			/ (this->vesselSegments[i]->getViscosity() * this->distance(this->vesselSegments[i]->vesselNodes[0], this->vesselSegments[i]->vesselNodes[1]) * LATTICE_CONSTANT)
			* ( this->vesselSegments[i]->vesselNodes[0]->pressure-this->vesselSegments[i]->vesselNodes[1]->pressure);
	}
}


void VesselGraph::updateShear()
{
#pragma omp parallel for
	for( int i=0; i<this->countVesselSegments; i++)
	//if( (this->vesselSegments[i]->vesselNodes[0]->getType() == VESSEL || this->vesselSegments[i]->vesselNodes[0]->getType() == ROOT) &&
	//	(this->vesselSegments[i]->vesselNodes[1]->getType() == VESSEL || this->vesselSegments[i]->vesselNodes[1]->getType() == ROOT))
	{
		this->vesselSegments[i]->shear 
			= 1./(2. * distance( this->vesselSegments[i]->vesselNodes[0], this->vesselSegments[i]->vesselNodes[1]) * LATTICE_CONSTANT)
			* this->vesselSegments[i]->radius
			* fabs(this->vesselSegments[i]->vesselNodes[0]->pressure-this->vesselSegments[i]->vesselNodes[1]->pressure);
	}
}


void VesselGraph::updatePressure( VesselNode **borderVesselNodes, float *borderPressure, int borderSize)
{
	// set border values
	for( int i=0; i<borderSize; i++)
		borderVesselNodes[i]->pressure = borderPressure[i];
		
	// iterative approximation
	float maxError = 0.;
	int its = 0;
	do{
		fprintf( stderr, "\rupdatePressure: %f \b", maxError);
		maxError = 0.;
		for( int i=0; i<this->countVesselNodes; i++)
			//if( this->vesselNodes[i]->countNeighbors>1)
			if( this->vesselNodes[i]->getType() == VESSEL || this->vesselNodes[i]->getType() == TIP)
			//if( this->vesselNodes[i]->type != ROOT)
			{
				/*float pressure = 0.;
				float count_neighbors = 0.;
				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				if( this->vesselNodes[i]->neighbors[v]->getType() == VESSEL || this->vesselNodes[i]->neighbors[v]->getType() == ROOT || 
						this->vesselNodes[i]->neighbors[v]->getType() == TIP	){
					pressure += this->vesselNodes[i]->neighbors[v]->pressure;
					count_neighbors++;
				}
				//pressure /= (float)this->vesselNodes[i]->countNeighbors;
				if(count_neighbors==0)
					pressure = 0.;
				else
					pressure /= count_neighbors;
				if( maxError < fabs(pressure - this->vesselNodes[i]->pressure))
					maxError = fabs(pressure - this->vesselNodes[i]->pressure);
				this->vesselNodes[i]->pressure = pressure;
				*/
				float pressure = 0.;
				float pressure2 = 0.;
				float count_neighbors = 0.;
				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				if( this->vesselNodes[i]->neighbors[v]->getType() == VESSEL || this->vesselNodes[i]->neighbors[v]->getType() == ROOT || 
						this->vesselNodes[i]->neighbors[v]->getType() == TIP	){
					//pressure += this->vesselNodes[i]->neighbors[v]->pressure;
					float G = pow( this->vesselNodes[i]->branches[v]->radius, 4) / (this->vesselNodes[i]->branches[v]->getViscosity() * distance(vesselNodes[i], vesselNodes[i]->neighbors[v]));

					pressure  += G*vesselNodes[i]->neighbors[v]->pressure;
					pressure2 += G;
					count_neighbors++;
				}
				//pressure /= (float)this->vesselNodes[i]->countNeighbors;
				if(count_neighbors==0)
					pressure = 0.;
				else
					pressure /= pressure2;
					//pressure /= count_neighbors;
				if( maxError < fabs(pressure - this->vesselNodes[i]->pressure))
					maxError = fabs(pressure - this->vesselNodes[i]->pressure);
				this->vesselNodes[i]->pressure = pressure;
			}
		its++;
	}while( maxError > 1e-10 && its < 1000);
	fprintf( stderr, "INFO: %i\n", its);
}

//SparseMatrix *sA=0;
//float *b=0;
//float *x=0;

typedef double pressure_t;

void VesselGraph::updatePressureNEW()
{
	pressure_t *b=0;
	pressure_t *x=0;

	//fprintf( stderr, "updatePressure\n");
	//if(sA==0)
	SparseMatrix<pressure_t> *sA = new SparseMatrix<pressure_t>(this->countVesselNodes,this->countVesselNodes);
	//if(b==0)
		b = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	//if(x==0)
		x = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	
	/*v0 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v1 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v2 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v3 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v4 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v5 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v6 = (float*) malloc( sizeof(float) * this->countVesselNodes);*/

	
	// INIT PRESSURE VALUES
	for( int i=0; i<this->countVesselNodes; i++){
		if( this->vesselNodes[i]->getType() == ROOT)
			x[i] = this->vesselNodes[i]->pressure;
		else
			//x[i] = this->vesselNodes[i]->pressure;
			x[i] = (MAX_PRESSURE+MIN_PRESSURE)/2.;
			//x[i] = RAND01;// * MAX_PRESSURE;
		if( i!=this->vesselNodes[i]->index){
			fprintf( stderr, "ERROR! Wrong index of node %i: %i\n", i, this->vesselNodes[i]->index);
			exit( 0);			
		}
	}
		
	// Construct Matrix
#pragma omp parallel for
	for( int i=0; i<this->countVesselNodes; i++){
		sA->resetRow(i);
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT 
			sA->set(i,i, 1);
			b[i] = this->vesselNodes[i]->pressure;
		}else{
		//if( this->vesselNodes[i]->getType() == VESSEL || this->vesselNodes[i]->getType() == TIP){
			pressure_t sumG = 0.;
			pressure_t G;
			pressure_t count_neighbors = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			/*if( this->vesselNodes[i]->neighbors[v]->getType() == VESSEL ||
				this->vesselNodes[i]->neighbors[v]->getType() == ROOT   || 
				this->vesselNodes[i]->neighbors[v]->getType() == TIP	)*/
			{
				G = /*PI/8. */ pow( this->vesselNodes[i]->branches[v]->radius, 4) / (this->vesselNodes[i]->branches[v]->getViscosity() * distance(vesselNodes[i], vesselNodes[i]->neighbors[v]) * LATTICE_CONSTANT);
				sumG += G;
				
				if(isnan(G)){
					fprintf( stderr, "ERROR!\n");
					exit( 0);
				}
				
				//sA->set(i,this->vesselNodes[i]->neighbors[v]->index, /*vesselNodes[i]->neighbors[v]->pressure*/G);
				//count_neighbors++;
			}
			//sA->set(i,i,sumG);
			
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			/*if( this->vesselNodes[i]->neighbors[v]->getType() == VESSEL ||
				this->vesselNodes[i]->neighbors[v]->getType() == ROOT   || 
				this->vesselNodes[i]->neighbors[v]->getType() == TIP	)*/
			{
				G = /*PI/8. */ pow( this->vesselNodes[i]->branches[v]->radius, 4) / (this->vesselNodes[i]->branches[v]->getViscosity() * distance(vesselNodes[i], vesselNodes[i]->neighbors[v]) * LATTICE_CONSTANT);

				if(isnan(G)){
					fprintf( stderr, "ERROR!\n");
					exit( 0);
				}
				if(sumG==0){
					fprintf( stderr, "ERROR: sumG=%lf!\n", sumG);
					exit( 0);
				}

				if(isnan(-G/sumG)){
					fprintf( stderr, "ERROR: G=%e, sumG=%e\n", G, sumG);
					exit( 0);
				}

				sA->set(i,this->vesselNodes[i]->neighbors[v]->index, -G/sumG);
				count_neighbors++;
			}
			
			sA->set(i,i,1.);//sumG);
			b[i] = 0.;
		}	
	}
		
	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);
	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	//sA->printMatrix( "A", "%10.3lf ");
	//fprintf( stderr, "Precondition System\n");

	Solver<pressure_t> *S = new Solver<pressure_t>( this->countVesselNodes, Solver<pressure_t>::BiCGSTAB, 1e-15, 10000 );
	//S->maxit = 10000;
	//S->maxerr = 1e-8;
	/*S->PreconditionJacobi( sA, b);
	S->maxerr = 1e-8;
	S->maxit = 2000;
	fprintf( stderr, "Solve System\n");
	fprintf( stderr, "updatePressure finished after %i iterations\n", S->solve( sA, b, x));
	 */

	/*SparseMatrix<pressure_t> *B = sA->copy();//new SparseMatrix<pressure_t>( this->countVesselNodes, this->countVesselNodes);
	S->IncompleteLUfactorization( B);
	pressure_t *temp = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	//S->SolveLowerTriangular( B, b, temp);
	//S->SolveUpperTriangular( B, temp, x);
	S->SolveLU( B,b,temp,x);*/
	fprintf( stderr, "updatePressure finished after %i iterations\n", S->solve( sA, b, x));



	// NORM: ||A*x - b||
	pressure_t *Ax_b, *Ax;
	Ax = Ax_b = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	// A*x
	SparseMatrix<pressure_t>::MatrixVectorProduct(sA,x, Ax);
	// A*x - b
	vectorDifference( Ax, b, Ax_b, this->countVesselNodes);
	// ||A*x - b||
	fprintf( stderr, "norm = %e\n", sqrt(dotProduct( Ax_b, Ax_b, this->countVesselNodes)));



	/*SparseMatrix<float> *B = new SparseMatrix<float>( this->countVesselNodes, this->countVesselNodes);
	S->IncompleteLUfactorization( sA, B);
	float *y = (float*) malloc( sizeof(float) * this->countVesselNodes);

	S->SolveLU( B,b,y,x);
	S->solve( sA, b, x);
	 */

	for( int i=0; i<this->countVesselNodes; i++)
		if( this->vesselNodes[i]->getType() != ROOT)
			this->vesselNodes[i]->pressure = x[i];

	//fprintf( stderr, "updatePressure finished\n");
	delete S;
	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updatePressureNEW2()
{
	// ALLOC MEMORY
	SparseMatrix<pressure_t> *sA = new SparseMatrix<pressure_t>(this->countVesselNodes,this->countVesselNodes);
	pressure_t *b = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	pressure_t *x = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);


	// INIT PRESSURE VALUES
	for( int i=0; i<this->countVesselNodes; i++){
		x[i] = this->vesselNodes[i]->pressure;

		if( i!=this->vesselNodes[i]->index){
			fprintf( stderr, "ERROR! Wrong index of node %i: %i\n", i, this->vesselNodes[i]->index);
			exit( 0);
		}
	}

	// CONSTRUCT LINEAR SYSTEM
#pragma omp parallel for
	for( int i=0; i<this->countVesselNodes; i++){

		sA->resetRow(i);

		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");


		if( this->vesselNodes[i]->getType() == ROOT || this->vesselNodes[i]->countNeighbors == 0){

			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			b[i] = this->vesselNodes[i]->pressure;

		}else{

			pressure_t G, sumG = 0;

			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				G = this->vesselNodes[i]->branches[v]->countParallelVessels *
						PI/8. * pow( this->vesselNodes[i]->branches[v]->radius, 4) / (this->vesselNodes[i]->branches[v]->getViscosity() * distance(vesselNodes[i], vesselNodes[i]->neighbors[v]) * LATTICE_CONSTANT);
				sA->set(i,this->vesselNodes[i]->neighbors[v]->index, -G);

				sumG += G;
			}

			sA->set(i,i,sumG);
			b[i] = 0.;
		}
	}

	// SOLVE SYSTEM

	Solver<pressure_t> *S = new Solver<pressure_t>( this->countVesselNodes, Solver<pressure_t>::BiCGSTAB, 1e-15, 10000 );
	S->PreconditionJacobi( sA, b);
	//S->maxerr = 1e-15;
	//S->maxit = 10000;//2000;
	//fprintf( stderr, "Solve System\n");
	//fprintf( stderr, "updatePressure finished after %i iterations\n",
	S->solve( sA, b, x);
	//);


	/*SparseMatrix<pressure_t> *B = sA->copy();//new SparseMatrix<pressure_t>( this->countVesselNodes, this->countVesselNodes);
	S->IncompleteLUfactorization( B);
	pressure_t *temp = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	S->SolveLU( B,b,temp,x);
	//fprintf( stderr, "updatePressure finished after %i iterations\n", S->solve( sA, b, x));*/



	// NORM: ||A*x - b||
	/*pressure_t *Ax_b, *Ax = Ax_b = (pressure_t*) malloc( sizeof(pressure_t) * this->countVesselNodes);
	// A*x
	SparseMatrix<pressure_t>::sparseMatrixVectorProduct(sA,x, Ax);
	// A*x - b
	vectorDifference( Ax, b, Ax_b, this->countVesselNodes);
	// ||A*x - b||
	fprintf( stderr, "norm = %e\n", sqrt(dotProduct( Ax_b, Ax_b, this->countVesselNodes)));*/


	for( int i=0; i<this->countVesselNodes; i++)
		if( this->vesselNodes[i]->getType() != ROOT)
			this->vesselNodes[i]->pressure = x[i];


	// FREE MEMORY

	delete S;
	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVessels( float dt, float **markerVesselsIn, float **markerVesselsOut)
{
	float *b=0;
	float *x=0;

	//if(sA==0)
	SparseMatrix<float> *sA = new SparseMatrix<float>(this->countVesselNodes,this->countVesselNodes);
	//if(b==0)
	b = (float*) malloc( sizeof(float) * this->countVesselNodes);
	//if(x==0)
	x = (float*) malloc( sizeof(float) * this->countVesselNodes);

	/*v0 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v1 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v2 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v3 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v4 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v5 = (float*) malloc( sizeof(float) * this->countVesselNodes);
	v6 = (float*) malloc( sizeof(float) * this->countVesselNodes);
*/

	// INIT PRESSURE VALUES
	/*for( int x=0; x<DOMAIN_SIZE_X; x++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			markerVesselsOut[x][y] = 0.;
*/
	for( int i=0; i<this->countVesselNodes; i++){
		x[i] = this->vesselNodes[i]->marker;
	}

	// Construct Matrix
//#pragma omp parallel for
	for( int i=0; i<this->countVesselNodes; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if(this->vesselNodes[i]->pressure==MIN_PRESSURE)
				b[i] = 0;
			else
				b[i] = 1;
		}else{
			float dIN  = 0.;
			float dOUT = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					dIN = fabs(this->vesselNodes[i]->branches[v]->flow) * dt; //* this->vesselNodes[i]->neighbors[v]->marker;
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, -dIN);
				}
				else
					dOUT+= fabs(this->vesselNodes[i]->branches[v]->flow) * dt; //* this->vesselNodes[i]->marker;

			sA->set(i,i, 1+dOUT);
			b[i] = this->vesselNodes[i]->marker;
		}
	}

	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);
	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	Solver<float> *S = new Solver<float>( this->countVesselNodes, Solver<float>::BiCGSTAB);
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<this->countVesselNodes; i++)
		this->vesselNodes[i]->marker = x[i];


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpace( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut)
{
	float *b=0;
	float *x=0;

	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;
	//if(sA==0)
	SparseMatrix<float> *sA = new SparseMatrix<float>(M+N,M+N);
	//if(b==0)
	b = (float*) malloc( sizeof(float) * (M+N));
	//if(x==0)
	x = (float*) malloc( sizeof(float) * (M+N));

	/*v0 = (float*) malloc( sizeof(float) * (M+N));
	v1 = (float*) malloc( sizeof(float) * (M+N));
	v2 = (float*) malloc( sizeof(float) * (M+N));
	v3 = (float*) malloc( sizeof(float) * (M+N));
	v4 = (float*) malloc( sizeof(float) * (M+N));
	v5 = (float*) malloc( sizeof(float) * (M+N));
	v6 = (float*) malloc( sizeof(float) * (M+N));
	 */

	// INIT PRESSURE VALUES
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}

	// Construct Matrix for Vascular Marker Concentration
	float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			b[i] = ARTERIAL_MARKER_CONC;
		}else{
			float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			float sumFlow=0.;
			
			// c_P^n+1
			sA->set(i,i, 1);
			
			// c_P,i^n+1 (for all neighbors i)
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = fabs(this->vesselNodes[i]->branches[v]->flow);
				}
				// out flow
				else{
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
				}
				sumFlow += flow;
				sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow/r);
			}
			
			// c_I^n+1
			int index = this->vesselNodes[i]->position[0] 
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X 
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL/r);
			
			// b
			b[i] = ARTERIAL_MARKER_CONC*sumFlow/r + volume/(dt*r) * this->vesselNodes[i]->marker;
		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x 
				      + y*DOMAIN_SIZE_X 
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				
				// c_I^n+1
				sA->set(i,i, 1);
				
				// r
				float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
					r+=MARKER_DIFFUSION_COEFFICIENT;
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
					r+=MARKER_DIFFUSION_COEFFICIENT;
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
					r+=MARKER_DIFFUSION_COEFFICIENT;
				
				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/r);
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/r);
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/r);
				}
				
				// c_P^n+1
				VesselNode *vn = this->octree->at( x,y,z);
				if( vn!=0)
					sA->set(i,vn->index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL/r);

				// b
				b[i] = (1.-volume)/(dt*r) * markerIntSpaceIn[x][y][z];
			}

				
	/*for( int i=0; i<this->countVesselNodes; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){
			// BOUNDARY: ROOT
			b[i] = ARTERIAL_MARKER_CONC;
		}else{
			float dIN  = 0.;
			float dOUT = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					dIN = fabs(this->vesselNodes[i]->branches[v]->flow) * dt; // * this->vesselNodes[i]->neighbors[v]->marker;
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, -dIN);
				}
				else
					dOUT+= fabs(this->vesselNodes[i]->branches[v]->flow) * dt; // * this->vesselNodes[i]->marker;

			sA->set(i,i, 1+dOUT);
			b[i] = this->vesselNodes[i]->marker;
		}
	}*/

	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);
	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB);
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpaceNEW( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut)
{
	float *b=0;
	float *x=0;

	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;
	//if(sA==0)
	SparseMatrix<float> *sA = new SparseMatrix<float>(M+N,M+N);
	//if(b==0)
	b = (float*) malloc( sizeof(float) * (M+N));
	//if(x==0)
	x = (float*) malloc( sizeof(float) * (M+N));

	/*v0 = (float*) malloc( sizeof(float) * (M+N));
	v1 = (float*) malloc( sizeof(float) * (M+N));
	v2 = (float*) malloc( sizeof(float) * (M+N));
	v3 = (float*) malloc( sizeof(float) * (M+N));
	v4 = (float*) malloc( sizeof(float) * (M+N));
	v5 = (float*) malloc( sizeof(float) * (M+N));
	v6 = (float*) malloc( sizeof(float) * (M+N));
	 */
 
	// INIT PRESSURE VALUES
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}

	// Construct Matrix for Vascular Marker Concentration
	//float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				b[i] = ARTERIAL_MARKER_CONC;
			else
				b[i] = 0;
		}else{
			//float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			//float sumFlow=0.;
			//float volume = 1.;

			// c_P^n+1
			sA->set(i,i, 1);

			// r
			float r = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
					r += fabs(this->vesselNodes[i]->branches[v]->flow);
			}
			r = 1. + dt*r + dt*EXCHANGE_RATE_VESSELS_INTERSTITIAL;

			// c_P,i^n+1 (for all neighbors i)
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow*dt/r);
				}
			}

			// c_I^n+1
			int index = this->vesselNodes[i]->position[0]
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r);

			// b
			b[i] = 1./r * this->vesselNodes[i]->marker;
		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;

				// c_I^n+1
				sA->set(i,i, 1);

				// r
				//float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				float r = 0.;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
					r+=1;
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
					r+=1;
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
					r+=1;
				r = 1. + dt*MARKER_DIFFUSION_COEFFICIENT*r + dt*EXCHANGE_RATE_VESSELS_INTERSTITIAL;

				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT*dt/r);
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT*dt/r);
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT*dt/r);
				}

				// c_P^n+1
				VesselNode *vn = this->octree->at( x,y,z);
				if( vn!=0)
					sA->set(i,vn->index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r);

				// b
				b[i] = 1./r * markerIntSpaceIn[x][y][z];
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);
	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB);
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpaceNEW2( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut,
		SparseMatrix<float> *&sA, float *&b, float *&x)
{
	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;
	if(sA==0)
	/*SparseMatrix<float> */sA = new SparseMatrix<float>(M+N,M+N);
	if(b==0)
	b = (float*) malloc( sizeof(float) * (M+N));
	if(x==0)
	x = (float*) malloc( sizeof(float) * (M+N));

	/*v0 = (float*) malloc( sizeof(float) * (M+N));
	v1 = (float*) malloc( sizeof(float) * (M+N));
	v2 = (float*) malloc( sizeof(float) * (M+N));
	v3 = (float*) malloc( sizeof(float) * (M+N));
	v4 = (float*) malloc( sizeof(float) * (M+N));
	v5 = (float*) malloc( sizeof(float) * (M+N));
	v6 = (float*) malloc( sizeof(float) * (M+N));*/
	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB );

 
	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
		sA->resetRow(i);
	}
	for( int i=M; i<M+N; i++){
		x[i] = 0;
		sA->resetRow(i);
	}

	// Construct Matrix for Vascular Marker Concentration
	//float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){
			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE){
				//b[i] = ARTERIAL_MARKER_CONC;
				b[i] = border;
			}
			else
				b[i] = 0;
		}else{
			//float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			//float sumFlow=0.;
			
			// vessel volume
			float volume = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				volume += 0.5*LATTICE_CONSTANT
						* PI * this->vesselNodes[i]->branches[v]->radius*this->vesselNodes[i]->branches[v]->radius;
			

			// c_P^n+1
			sA->set(i,i, 1);

			// r
			float r = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				//if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
				//	r += fabs(this->vesselNodes[i]->branches[v]->flow);

				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure)
					r += fabs(this->vesselNodes[i]->branches[v]->flow);
			}
			r = 1. + dt/volume*r + dt/volume*EXCHANGE_RATE_VESSELS_INTERSTITIAL;

			// c_P,i^n+1 (for all neighbors i)
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow*dt/r/volume);
					if( isnan(flow*dt/r/volume)){
						fprintf( stderr, "flow=%lf dt=%lf r=%lf volume=%lf\n",flow,dt,r,volume);
						exit(0);
					}
				}
			}

			// c_I^n+1
			int index = this->vesselNodes[i]->position[0]
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume);
			if( isnan(EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume)){
				fprintf( stderr, "EXCHANGE_RATE_VESSELS_INTERSTITIAL=%lf dt=%lf r=%lf volume=%lf\n",EXCHANGE_RATE_VESSELS_INTERSTITIAL,dt,r,volume);
				exit(0);
			}

			// b
			b[i] = 1./r * this->vesselNodes[i]->marker;
			if( isnan(b[i])){
				fprintf( stderr, "r=%lf this->vesselNodes[i]->marker=%lf\n",r,this->vesselNodes[i]->marker);
				exit(0);
			}

		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;


				// volume: tissue
				VesselNode *vn = this->octree->at( x,y,z);
				float volume = LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT;
				if( vn!=0)
				for( int v=0; v<vn->countNeighbors; v++)
					volume -= 0.5*LATTICE_CONSTANT
						* PI * vn->branches[v]->radius*vn->branches[v]->radius;
				
				if(volume < 0.){
					fprintf( stderr, "ERROR: negative volume: %lf!!\n", volume);
					fprintf( stderr, "ERROR: max volume: %lf = %lf^3\n", (float)LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT, (float)LATTICE_CONSTANT);
					for( int v=0; v<vn->countNeighbors; v++)
						fprintf( stderr, "%i neighbor: radius=%lf, volume=%lf\n", v, vn->branches[v]->radius, 0.5*LATTICE_CONSTANT	* PI * vn->branches[v]->radius*vn->branches[v]->radius);
					exit(0);
				}
						
				// c_I^n+1
				sA->set(i,i, 1);

				// r
				//float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				float r = 0.;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
					r+=1;
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
					r+=1;
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
					r+=1;
				r = 1. 
				  + dt * MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * r 
				  + dt * EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume;

				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}

				// c_P^n+1
				if( vn!=0){
					sA->set(i,vn->index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume * dt/r);
					if( isnan(-EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume * dt/r)){
						fprintf( stderr, "r=%lf volume=%lf\n",r,volume);
						exit(0);
					}
				}

				// b
				b[i] = 1./r * markerIntSpaceIn[x][y][z];
				if( isnan(b[i])){
					fprintf( stderr, "r=%lf markerIntSpaceIn[x][y][z]=%lf\n",r,markerIntSpaceIn[x][y][z]);
					exit(0);
				}
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);

	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	//fprintf( stderr, "Precondition System\n");

	//S->PreconditionJacobi( sA, b);
	//S->maxerr = 1e-6;
	//fprintf( stderr, "Solve System\n");
	if( S->solve( sA, b, x) == 0){
		fprintf(stderr, "[AIF=%f] ", border);
		printVector( b, sA->columns(), "b", "%10.3f");
		//sA->printMatrix( "A", "%10.3f");
		//exit(0);
	}

	/*S->IncompleteLUfactorization( sA);
	float *y = (float*) malloc( sizeof(float) * (M+N));
	S->SolveLU( sA, b, y, x);
	free(y);*/
	delete S;

	for( int i=0; i<M; i++){
		this->vesselNodes[i]->marker = x[i];
		if( isnan(x[i]) ){
			fprintf( stderr, "x[%i] = %lf\n", i, x[i]);
			exit(0);
		}
	}

	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
				if( isnan(x[i]) ){
					fprintf( stderr, "x[%i] = %lf\n", i, x[i]);
					exit(0);
				}
			}


	/*delete sA;
	free(b);
	free(x);*/
}

double tempMarker = 0;

void VesselGraph::updateMarkerVesselsExplicit( float dt, float border, float *&dx)
{
	int M = this->countVesselNodes;
	if(dx==0)
	dx = (float*) malloc( sizeof(float) * (M));


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		dx[i] = this->vesselNodes[i]->marker;
	}


	// Construct Matrix for Vascular Marker Concentration
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT /*&& this->vesselNodes[i]->pressure==MAX_PRESSURE*/){

			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);

			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				dx[i] = (border - this->vesselNodes[i]->marker);// / dt;
			else{
				/*double c       = fabs(this->vesselNodes[i]->getFlowFromNeighbor(0) * dt / this->getVascularVolume( this->vesselNodes[i]));
				double c_minus = fabs(this->vesselNodes[i]->getFlowFromNeighbor(0) * dt / this->getVascularVolume( this->vesselNodes[i]->neighbors[0]));
				double u_minus = this->vesselNodes[i]->neighbors[0]->marker
						       + 0.5 * (this->vesselNodes[i]->neighbors[0]->marker - this->vesselNodes[i]->marker) * (1. - c_minus);
				double u_plus =  this->vesselNodes[i]->marker
							   + 0.5 * (this->vesselNodes[i]->marker - (2.*this->vesselNodes[i]->marker - this->vesselNodes[i]->neighbors[0]->marker)) * (1. - c_minus);
					       	   //+ 0.5 * (this->vesselNodes[i]->neighbors[0]->marker - this->vesselNodes[i]->marker) * (1. - c);

				dx[i] = c * ( u_minus - u_plus );*/

				//float courantNumber = fabs( this->vesselNodes[i]->branches[0]->flow  * dt / this->getVascularVolume( this->vesselNodes[i]));
				//dx[i] = courantNumber * (this->vesselNodes[i]->getNeighborUpstream(0)->marker -  this->vesselNodes[i]->marker);
				/*dx[i] =  -  this->vesselNodes[i]->marker
					  + 2 * this->vesselNodes[i]->getNeighborUpstream(0)->marker
					      - this->vesselNodes[i]->getNeighborUpstream(0)->getNeighborUpstream(0)->marker;*/
				dx[i] =  -  this->vesselNodes[i]->marker
					  + 2 * this->vesselNodes[i]->getNeighborUpstream(0)->marker
					      - tempMarker;
				tempMarker = this->vesselNodes[i]->getNeighborUpstream(0)->getNeighborUpstream(0)->marker;
			}
		}else{

			dx[i] = 0;
			float outFlow=0, inFlow=0, diffusion=0;

			switch( 4){

			case 0: // 3-POINT EXPLICIT WITH NUMERICAL DIFFUSION (SECOND ORDER)
			{
				float flow = 0.;

				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				{
					float courantNumber = this->vesselNodes[i]->getFlowFromNeighbor(v) * dt / this->getVascularVolume( this->vesselNodes[i]);
					float diffusionCoefficient = fabs(courantNumber) / 2.;

					// flow
					flow += courantNumber * (this->vesselNodes[i]->neighbors[v]->marker + this->vesselNodes[i]->marker)/2.;

					// "numerical" diffusion
					diffusion += diffusionCoefficient * (this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker);
				}

				dx[i] = flow + diffusion;
				break;
			}

			case 1:	// EXPLICIT UPWIND (FIRST ORDER UPWIND)
			{
				float flow = 0.;

				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				{
					// in flow
					if (this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure) {

						float courantNumber = this->vesselNodes[i]->getFlowFromNeighbor(v) * dt / this->getVascularVolume( this->vesselNodes[i]);
						flow += courantNumber
								* (this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker);
					}
				}

				dx[i] = flow;
				break;
			}

			case 2:	// FIRST ORDER LAX-FRIEDRICHS
			{
				/*for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				{
						// out flow
						if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure){

							outFlow	+= fabs(this->vesselNodes[i]->branches[v]->flow) / this->getVascularVolume( this->vesselNodes[i])
									* this->vesselNodes[i]->neighbors[v]->marker;
						}

						// in flow
						if (this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure) {

							inFlow += fabs(this->vesselNodes[i]->branches[v]->flow) / this->getVascularVolume( this->vesselNodes[i])
									* this->vesselNodes[i]->neighbors[v]->marker;
						}

						//weight = 1.;//fabs(this->vesselNodes[i]->branches[v]->flow) / this->getVascularVolume( this->vesselNodes[i]) * dt;
						//weightSum += weight;
						//neighborConc += weight * (this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker);
				}*/
				float flow=0.;
				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++){
					float courantNumber = this->vesselNodes[i]->getFlowFromNeighbor(v) * dt / this->getVascularVolume( this->vesselNodes[i]);
					flow += (courantNumber+1.)*vesselNodes[i]->neighbors[v]->marker;
				}

				//dx[i] = dt * (inFlow - outFlow + diffusion) + neighborConc/weightSum;
				dx[i] = - this->vesselNodes[i]->marker
						//+ (this->vesselNodes[i]->neighbors[0]->marker - 2*this->vesselNodes[i]->marker + this->vesselNodes[i]->neighbors[1]->marker) / 2.
						+ flow / 2.;
						//- dt * (outFlow - inFlow) / 2.;
						//+ dt / getVascularVolume( this->vesselNodes[i])
						//		*(vesselNodes[i]->getFlowFromNeighbor(0)*vesselNodes[i]->neighbors[0]->marker
						//		+ vesselNodes[i]->getFlowFromNeighbor(1)*vesselNodes[i]->neighbors[1]->marker) / 2.;
				break;
			}

			case 3: // MACCORMACK / Lax-Wendroff
			{
				dx[i]=0.;
				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				{
					float courantNumber = this->vesselNodes[i]->getFlowFromNeighbor(v) * dt / this->getVascularVolume( this->vesselNodes[i]);

					// advection
					//dx[i] += courantNumber
					//		* this->vesselNodes[i]->neighbors[v]->marker / 2.;

					// diffusion
					//dx[i] += pow( courantNumber, 2.)
					//		* (this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker) / 2.;

					// convection = advection + diffusion
					dx[i] += courantNumber*(1.-fabs(courantNumber))
							* (this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker) / 2.;

				}

				break;
			}

			case 4: // PIECEWISE-LINEAR
			{
				//float inFlow = 0.;
				//float outFlow = 0.;
				float flux = 0.;

				// for all flow from/to neighbors
				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				//if( fabs( this->vesselNodes[i]->pressure - this->vesselNodes[i]->neighbors[v]->pressure) >0.00001 )
				{
					float flow;// = fabs(this->vesselNodes[i]->branches[v]->flow);
					float sign;
					VesselNode* origin;
					VesselNode* destination;
					if( this->vesselNodes[i]->pressure > this->vesselNodes[i]->neighbors[v]->pressure ){
						// flow to neighbor: f_{i+1/2}
						origin      = this->vesselNodes[i];
						destination = this->vesselNodes[i]->neighbors[v];
						flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
						sign = +1;
					}else{
						// flow from neighbor: f_{i-1/2}
						origin      = this->vesselNodes[i]->neighbors[v];
						destination = this->vesselNodes[i];
						flow = fabs(this->vesselNodes[i]->branches[v]->flow);
						sign = +1;
					}

					// SLOPE
					// Donor cell
					float slope = 0.;
					// Centered Euler
					//float slope = (destination->marker - origin->marker);
					// Lax-Wendroff
					//float slope = (destination->marker - origin->marker) * (1-fabs(flow * dt / this->getVascularVolume( this->vesselNodes[i])));
					// Beam-Warming
					/*float slope = 0.;
					for( int vv=0; vv<origin->countNeighbors; vv++)
						if( origin->neighbors[vv]!=destination){
							slope += origin->marker - origin->neighbors[vv]->marker;
						}*/
					// Fromm
					/*float slope = 0.;
					for( int vv=0; vv<origin->countNeighbors; vv++)
						if( origin->neighbors[vv]!=destination){
							slope += 0.5*(destination->marker - origin->neighbors[vv]->marker);
						}*/

					// SLOPE LIMITER
					/*float predecessorMarker = 0.;
					for( int vv=0; vv<origin->countNeighbors; vv++)
						if( origin->neighbors[vv]!=destination)
							predecessorMarker += origin->neighbors[vv]->marker;
					float r = (origin->marker - predecessorMarker) / (destination->marker - origin->marker);*/

					// minmod slope
					//slope *= MAX(0., MIN(r, 1.));
					//float slope = MIN( MAX( origin->marker - predecessorMarker,0), MAX( destination->marker - origin->marker,0) )
					//			+ MAX( MIN( origin->marker - predecessorMarker,0), MIN( destination->marker - origin->marker,0) );
					// vanLeer slope
					//float slope = (( origin->marker - predecessorMarker)*( destination->marker - origin->marker)>0. ?
					//		2./(1./( origin->marker - predecessorMarker) + 1./( destination->marker - origin->marker)) : 0.);
					// superBee
					//slope *= MAX( 0., MAX( MIN(2.*r,1.), MIN(r,2.) ));
					//float slope = (( origin->marker - predecessorMarker)>0 && ( destination->marker - origin->marker)>0 ? +1 : -1) *
					//		MIN( MIN( fabs(origin->marker - predecessorMarker), fabs( destination->marker - origin->marker) ),
					//			 0.5* MAX( fabs(origin->marker - predecessorMarker), fabs( destination->marker - origin->marker) )	);

					// FLUX
					float courantNumber = flow * dt / this->getVascularVolume( this->vesselNodes[i]);
					//float courantNumber = flow * dt / this->getVascularVolume( origin);
					//flux += flow * (origin->marker + 0.5 * slope * sign * (1. - fabs(courantNumber)));
					flux += flow * (origin->marker + 0.5 * slope);

				}

				dx[i] = flux * dt / this->getVascularVolume( this->vesselNodes[i]);

				break;
			}

			case 5: // LUNA
			{

				float flux = 0.;
				float u_i_minus_1;
				float u_i;
				float u_i_plus_1;

				for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++){
					float cfl;
					if( this->vesselNodes[i]->pressure > this->vesselNodes[i]->neighbors[v]->pressure ){
						// flow to neighbor: f_{i+1/2}
						cfl    = -fabs(this->vesselNodes[i]->branches[v]->flow)  * dt / this->getVascularVolume( this->vesselNodes[i]);
						//sign = +1;

						u_i_plus_1 = this->vesselNodes[i]->neighbors[v]->marker;
						u_i        = this->vesselNodes[i]->marker;
						for( int vv=0; vv<this->vesselNodes[i]->countNeighbors; vv++)
							if( this->vesselNodes[i]->neighbors[vv]!=this->vesselNodes[i]->neighbors[v])
								u_i_minus_1 += this->vesselNodes[i]->neighbors[vv]->marker;
					}else{
						// flow from neighbor: f_{i-1/2}
						cfl    = +fabs(this->vesselNodes[i]->branches[v]->flow) * dt / this->getVascularVolume( this->vesselNodes[i]);
						//sign = +1;
						u_i_plus_1 = this->vesselNodes[i]->marker;
						u_i        = this->vesselNodes[i]->neighbors[v]->marker;
						for( int vv=0; vv<this->vesselNodes[i]->neighbors[v]->countNeighbors; vv++)
							if( this->vesselNodes[i]->neighbors[v]->neighbors[vv]!=this->vesselNodes[i])
								u_i_minus_1 += this->vesselNodes[i]->neighbors[v]->neighbors[vv]->marker;
					}

					float r = (u_i - u_i_minus_1) / (u_i_plus_1 - u_i);
					float slope = this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker;
					float fluxLimiter = MAX( 0, MIN(r, 1));

					//flux += - 1./2. * cfl * slope
					//	    + (1-(1-(cfl)) * fluxLimiter) * slope;

					//slope =

					// w_i-1
					flux += cfl/2.*(this->vesselNodes[i]->neighbors[v]->marker)
							//- cfl/2.*(1.-fabs(cfl))*fluxLimiter*(this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker)
							;

					// w_i+1
					//flux += cfl/2.*(this->vesselNodes[i]->neighbors[v]->marker)
					//		- cfl/2.*(1+cfl)*fluxlimiter*(this->vesselNodes[i]->neighbors[v]->marker - this->vesselNodes[i]->marker);

				}

				dx[i] = flux;
				break;
			}
			}
		}
	}

	// ANTI-DIFFUSION
	/*for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else
			this->vesselNodes[i]->marker += dx[i];
	}
	// calculate raw anti-diffusion fluxes
	float *fad = (float*) malloc( sizeof(float) * this->countVesselSegments);
	for( int i=0; i<countVesselSegments; i++){

		float courantNumber = this->vesselSegments[i]->flow;// * dt / this->getVascularVolume( this->vesselNodes[i]);
		float diffusionCoefficient = fabs(courantNumber) / 2.;
		float antiDiffCoefficient  = 4.*diffusionCoefficient ;

		//float diffusionCoefficient = fabs(this->vesselSegments[i]->flow) / 2.;
		//float antiDiffCoefficient  = diffusionCoefficient - fabs(this->vesselSegments[i]->flow) / 2.;

		fad[i] = antiDiffCoefficient
				* (this->vesselSegments[i]->vesselNodes[1]->marker - this->vesselSegments[i]->vesselNodes[0]->marker);
	}
	// limite/correct fluxes
	for( int i=0; i<countVesselSegments; i++){

		float courantNumber = this->vesselSegments[i]->flow;// * dt / this->getVascularVolume( this->vesselNodes[i]);
		float diffusionCoefficient = fabs(courantNumber) / 2.;
		float antiDiffCoefficient  = 5.*diffusionCoefficient;


		float sign = (fad[i]<0. ? -1 : 1);
		float min = fabs( fad[i]);

		// p_i+2 - p_i+1
		float temp=0.;
		for( int ii=0; ii<this->vesselSegments[i]->vesselNodes[1]->countNeighbors; ii++)
			if(this->vesselSegments[i]->vesselNodes[1]->branches[ii] != this->vesselSegments[i]){
				temp += antiDiffCoefficient*sign * (
						// p_i+2
						this->vesselSegments[i]->vesselNodes[1]->neighbors[ii]->marker
						// p_i+1
						- this->vesselSegments[i]->vesselNodes[1]->marker);
			}
		fprintf( stderr, "[%i] fad_i+1 = %e\n", i, temp);
		min = MIN( min, temp);

		// p_i - p_i-1
		temp=0.;
		for( int ii=0; ii<this->vesselSegments[i]->vesselNodes[0]->countNeighbors; ii++)
			if(this->vesselSegments[i]->vesselNodes[0]->branches[ii] != this->vesselSegments[i]){
				temp += antiDiffCoefficient*sign * (
						// p_i
						this->vesselSegments[i]->vesselNodes[0]->marker
						// p_i-1
						- this->vesselSegments[i]->vesselNodes[0]->neighbors[ii]->marker);
			}
		fprintf( stderr, "[%i] fad_i-1 = %e\n", i, temp);
		min = MIN( min, temp);

		// correct fad
		fprintf( stderr, "[%i] fad = %e ---> corrected = %e\n", i, fad[i], sign * MAX(0, min));
		fad[i] = sign * MAX(0, min);
	}
	// perform antidiffusive correction
	for( int i=0; i<M; i++){
		dx[i] = 0.;
	}
	for( int i=0; i<countVesselSegments; i++){
		dx[this->vesselSegments[i]->vesselNodes[0]->index] -= fad[i] * dt / this->getVascularVolume( this->vesselSegments[i]->vesselNodes[0]);
		dx[this->vesselSegments[i]->vesselNodes[1]->index] += fad[i] * dt / this->getVascularVolume( this->vesselSegments[i]->vesselNodes[1]);
	}
*/
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else
			this->vesselNodes[i]->marker += dx[i];
	}
}


#define SIGN( a) (a<0 ? -1 : 1)

float GetSlope( float previous, float origin, float destination, float courantNumber, char schema)
{
	switch (schema) {
	case UPWIND:			// (1st order)
		return 0.;

	case CENTERED_EULER:	// (2nd order)
		return destination - origin;

	case LAX_WENDROFF:		// (2nd order)
	case MACCORMACK:		// (2nd order)
		return (destination - origin) * (1 - courantNumber);

	case MINMOD:
		return MIN( MAX( origin-previous, 0), MAX( destination-origin, 0)) +
			   MAX( MIN( origin-previous, 0), MIN( destination-origin, 0));

	case SUPERBEE:
		return (SIGN( origin-previous) + SIGN( destination-origin)) *
				MIN(
						MIN( fabs( origin-previous), fabs( destination-origin)),
						0.5 * MAX( fabs( origin-previous), fabs( destination-origin))
				);
	}
}

/*float GetSlope( VesselNode *origin, VesselNode *destination, float courantNumber, char schema)
{
	switch (schema) {
	case UPWIND:			// (1st order)
		return 0.;

	case CENTERED_EULER:	// (2nd order)
		return destination->marker - origin->marker;

	case LAX_WENDROFF:		// (2nd order)
	case MACCORMACK:		// (2nd order)
		return (destination->marker - origin->marker) * (1 - courantNumber);
	}
}*/

float GetPreviousRefined( VesselNode *vn, float *s, int S)
{
	float inflowSum = 0.;
	float inflowTimesConcentrationSum = 0.;
	for( int vv=0; vv<vn->countNeighbors; vv++){
		if( vn->pressure < vn->neighbors[vv]->pressure){
			int index = vn->branches[vv]->index*S + (S-1);
			inflowSum += fabs(vn->branches[vv]->flow);
			inflowTimesConcentrationSum += fabs(vn->branches[vv]->flow) * s[index];
		}
	}
	return inflowTimesConcentrationSum / inflowSum;
}

/*int GetIndexPredecessor( int index, VesselGraph *vg, int S)
{
	int N = vg->countVesselNodes;

	if( index < N){
		// node = VesselNode
		if( S){

		}else{

		}
	}
}*/

void VesselGraph::updateMarkerVesselsExplicitRefined( float dt, float border, float *&dx, float *&ds, float *&s, char schema, int refinement)
{
	//char schema = UPWIND;
	//char schema = LAX_WENDROFF;
	//char schema = MINMOD;
	//char schema = SUPERBEE;

	// Number of Vesselnodes
	int M = this->countVesselNodes;
	if(dx==0)
		dx = (float*) malloc( sizeof(float) * (M));

	// Number of VesselSubNodes per VesselSegment
	int S = refinement;//30;
	//fprintf(stderr, "S=%d\n", S);

	// Number of VesselSubNodes
	int N = this->countVesselSegments * S;
	if(ds==0)
		ds = (float*) malloc( sizeof(float) * N);
	if(s==0)
		s = (float*) malloc( sizeof(float) * N);



	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		dx[i] = 0;
	}


	// Construct Matrix for Vascular Marker Concentration
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT){

			// BOUNDARY: ROOT
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				dx[i] = (border - this->vesselNodes[i]->marker);
			else{
				dx[i] =  -  this->vesselNodes[i]->marker
					  + 2 * this->vesselNodes[i]->getNeighborUpstream(0)->marker
					      - this->vesselNodes[i]->getNeighborUpstream(0)->getNeighborUpstream(0)->marker;
				//int index = this->vesselNodes[i]->branches[v]->index*S;

				// (u_i-1) ----> s_0 ----> (s_1)
				int index = this->vesselNodes[i]->branches[0]->index*S;
				float volumeSubnode = getVascularVolume( vesselNodes[i]->branches[0]) * 1./(S+2.);
				float courantNumberSubnode = fabs( vesselNodes[i]->branches[0]->flow) * dt/volumeSubnode;
				//previous = GetPreviousRefined( )
				float u_in  = this->vesselNodes[i]->getNeighborUpstream(0)->marker + 0.5*GetSlope( this->vesselNodes[i]->getNeighborUpstream(0)->marker, this->vesselNodes[i]->getNeighborUpstream(0)->marker, s[index], courantNumberSubnode, schema);
				float u_out = s[index] + 0.5*GetSlope( this->vesselNodes[i]->getNeighborUpstream(0)->marker, s[index], s[index+1], courantNumberSubnode, schema);
				ds[index] = courantNumberSubnode * (u_in - u_out);

				// (s_0) ----> s_1 ---->
				ds[index+1] = -s[index+1] + 2*s[index] - this->vesselNodes[i]->getNeighborUpstream(0)->marker;

			}
		}else{

			//dx[i] = 0;
			//float flux = 0.;

			// for all flow from/to neighbors
			for (int v = 0; v < this->vesselNodes[i]->countNeighbors; v++)
			{
				//float flow;

				float previous;
				float origin;
				float destination;

				/*if (this->vesselNodes[i]->pressure > this->vesselNodes[i]->neighbors[v]->pressure) {
					// flow to neighbor: f_{i+1/2}
					//flow = -fabs(this->vesselNodes[i]->branches[v]->flow);

					origin = this->vesselNodes[i]->marker;
					destination = ( S==0 ? this->vesselNodes[i]->neighbors[v]->marker : s[vesselNodes[i]->branches[v]->index*S]);

					float u_in, u_out;

					float volumeNode = this->getVascularVolume(this->vesselNodes[i]) * 2./(S+2.);
					float courantNumber = fabs(this->vesselNodes[i]->branches[v]->flow) * dt / volumeNode;
					float u_out = origin + 0.5*GetSlope( origin, destination, courantNumber, UPWIND);
					dx[i] -= u_out * fabs(this->vesselNodes[i]->branches[v]->flow);

					// FLUX ALONG SEGMENT
					float volumeSubnode = getVascularVolume( vesselNodes[i]->branches[v]) * 1./(S+2.);
					float courantNumberSubnode = fabs( vesselNodes[i]->branches[v]->flow) * dt/volumeSubnode;
					for(int is=0; is<S; is++){

						// in flux
						float u_in = u_out;

						// out flux
						int index = this->vesselNodes[i]->branches[v]->index*S + is;
						u_out = s[index] + 0.5*GetSlope( s[index], s[index+1], courantNumberSubnode, UPWIND);

						ds[index] =	courantNumberSubnode * (u_in - u_out);

					}

					dx[vesselNodes[i]->neighbors[v]->index] += u_out * fabs(this->vesselNodes[i]->branches[v]->flow);
				}*/

				if (this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure) {
					// flow to neighbor: f_{i+1/2}
					//flow = -fabs(this->vesselNodes[i]->branches[v]->flow);

					// (u_i-1) ----> s_0 ----> s_S ----> u_i

					origin = this->vesselNodes[i]->neighbors[v]->marker;
					destination = ( S==0 ? this->vesselNodes[i]->marker : s[vesselNodes[i]->branches[v]->index*S]);

					float u_in, u_out;

					// FLUX OUT OF PRECEDOR
					float volumeNode = this->getVascularVolume(this->vesselNodes[i]->neighbors[v]) * 2./(S+2.);
					float courantNumber = fabs(this->vesselNodes[i]->branches[v]->flow) * dt / volumeNode;
					// previous on refined
					previous = GetPreviousRefined( this->vesselNodes[i]->neighbors[v], s, S);
					u_out = origin + 0.5*GetSlope( previous, origin, destination, courantNumber, schema);
					//dx[vesselNodes[i]->neighbors[v]->index] -= courantNumber * u_out;

					// FLUX ALONG SEGMENT
					float volumeSubnode = getVascularVolume( vesselNodes[i]->branches[v]) * 1./(S+2.);
					float courantNumberSubnode = fabs( vesselNodes[i]->branches[v]->flow) * dt/volumeSubnode;
					for(int is=0; is<S; is++){

						// in flux
						float u_in = u_out;

						// out flux
						int index = this->vesselNodes[i]->branches[v]->index*S + is;
						previous = (is>0 ? s[index-1] : this->vesselNodes[i]->neighbors[v]->marker);
						float next = (is<S-1 ? s[index+1] : this->vesselNodes[i]->marker);
						u_out = s[index] + 0.5*GetSlope( previous, s[index], next, courantNumberSubnode, schema);

						ds[index] =	courantNumberSubnode * (u_in - u_out);

					}

					// FLUX INTO SUCCESSOR
					volumeNode = this->getVascularVolume(this->vesselNodes[i]) * 2./(S+2.);
					courantNumber = fabs(this->vesselNodes[i]->branches[v]->flow) * dt / volumeNode;
					dx[i] += courantNumber * u_out;
				}else{
					// u_i ----> (s_0)
					origin = this->vesselNodes[i]->marker;
					destination = ( S==0 ? this->vesselNodes[i]->neighbors[v]->marker : s[vesselNodes[i]->branches[v]->index*S]);
					previous = GetPreviousRefined( this->vesselNodes[i], s, S);

					// FLUX OUT OF PRECEDOR
					float volumeNode = this->getVascularVolume(this->vesselNodes[i]) * 2./(S+2.);
					float courantNumber = fabs(this->vesselNodes[i]->branches[v]->flow) * dt / volumeNode;
					//previous =
					float u_out = origin + 0.5*GetSlope( previous, origin, destination, courantNumber, schema);
					dx[i] -= courantNumber * u_out;
				}

			}
		}
	}

	for( int i=0; i<M; i++){
			this->vesselNodes[i]->marker += dx[i];
	}
	for( int i=0; i<N; i++){
			s[i] += ds[i];
	}
}


void VesselGraph::updateMarkerVesselsExplicitAntiDiffusion( float dt, float border, float *&dx)
{
	int M = this->countVesselNodes;
	if(dx==0)
	dx = (float*) malloc( sizeof(float) * (M));


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		dx[i] = this->vesselNodes[i]->marker;
	}


	// Construct Matrix for Vascular Marker Concentration
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){

			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);
			dx[i] = (border - this->vesselNodes[i]->marker) / dt;

		}else{

			// c_P^n+1
			dx[i] = 0.;
			float inFlow = 0.;
			float outFlow = 0.;
			float antiDiffusion = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
				//	x[i] -= this->vesselNodes[i]->marker * fabs(this->vesselNodes[i]->branches[v]->flow);
					antiDiffusion += PI * pow(this->vesselNodes[i]->branches[v]->radius,2) / 60.
									* (this->vesselNodes[i]->marker - this->vesselNodes[i]->neighbors[v]->marker);

				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					//dx[i] += this->vesselNodes[i]->neighbors[v]->marker * fabs(this->vesselNodes[i]->branches[v]->flow);

					inFlow  += fabs(this->vesselNodes[i]->branches[v]->flow) * this->vesselNodes[i]->neighbors[v]->marker;
					outFlow += fabs(this->vesselNodes[i]->branches[v]->flow) * this->vesselNodes[i]->marker;
					antiDiffusion += PI * pow(this->vesselNodes[i]->branches[v]->radius,2) / 60.
									* (this->vesselNodes[i]->marker - this->vesselNodes[i]->neighbors[v]->marker);

				}
			}
			dx[i] = (inFlow - outFlow + 100.*antiDiffusion) / this->getVascularVolume( this->vesselNodes[i]);
		}
	}

	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else
			this->vesselNodes[i]->marker += dt * dx[i];
	}
}


void VesselGraph::updateMarkerVesselsExplicit( float dt, float border, float *&dx, float *&x, int countVesselSegmentSubNodes)
{
	//int countVesselSegmentSubNodes = 100;
	int N = this->countVesselSegments * countVesselSegmentSubNodes;
	int M = this->countVesselNodes;
	if(dx==0){
		dx = (float*) malloc( sizeof(float) * (M + N));
		for( int i=0; i<M+N; i++)
			dx[i]=0;
	}
	if(x==0){
		x = (float*) malloc( sizeof(float) * (M + N));
		for( int i=0; i<M+N; i++)
			x[i]=0;
	}


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}


	// Construct Matrix for Vascular Marker Concentration: vessel nodes
	for( int i=0; i<M; i++)
	{
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){

			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);
			dx[i] = (border - this->vesselNodes[i]->marker) / dt;

		}else{

			// c_P^n+1
			dx[i] = 0.;
			float inFlow = 0.;
			float outFlow = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// 3-POINT
				if( false){
					// out flow
					//if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
					//	x[i] -= this->vesselNodes[i]->marker * fabs(this->vesselNodes[i]->branches[v]->flow);
				}

				// EXPLICIT UPWIND
				if( true){
					// in flow
					if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){

						float neighborConc;
						if(countVesselSegmentSubNodes)
							neighborConc = x[M + (1+this->vesselNodes[i]->branches[v]->index)*countVesselSegmentSubNodes-1];
						else
							neighborConc = this->vesselNodes[i]->neighbors[v]->marker;

						inFlow  += fabs(this->vesselNodes[i]->branches[v]->flow) * neighborConc;
						outFlow += fabs(this->vesselNodes[i]->branches[v]->flow) * this->vesselNodes[i]->marker;
					}
				}
			}
			dx[i] = (inFlow - outFlow) / (this->getVascularVolume( this->vesselNodes[i]) / (countVesselSegmentSubNodes+1.));
		}
	}

	// Construct Matrix for Vascular Marker Concentration: vessel segment sub nodes
	if(countVesselSegmentSubNodes)
	for( int i=0; i<this->countVesselSegments; i++){
		//fprintf( stderr, "%i == %i?\n", i, this->vesselSegments[i]->index);

		float flow = fabs(this->vesselSegments[i]->flow) / (this->getVascularVolume( this->vesselSegments[i]) / (countVesselSegmentSubNodes+1.));

		int index = M + i*countVesselSegmentSubNodes;

		if(this->vesselSegments[i]->vesselNodes[0]->pressure > this->vesselSegments[i]->vesselNodes[1]->pressure){
			// flow direction [0] -> [1]
			dx[index] = flow*( this->vesselSegments[i]->vesselNodes[0]->marker - x[index] );
		}else{
			// flow direction [1] -> [0]
			dx[index] = flow*( this->vesselSegments[i]->vesselNodes[1]->marker - x[index] );
		}
		//if(this->vesselSegments[i]->vesselNodes[0]->pressure == MAX_PRESSURE || this->vesselSegments[i]->vesselNodes[1]->pressure == MAX_PRESSURE)
		//fprintf( stderr, "dx = %e, x=%e\n", dx[index], x[index]);
		for( int ii=1; ii<countVesselSegmentSubNodes; ii++)
		{
			index++;
			dx[index] = flow*( x[index-1] - x[index] );
			//if(this->vesselSegments[i]->vesselNodes[0]->pressure == MAX_PRESSURE || this->vesselSegments[i]->vesselNodes[1]->pressure == MAX_PRESSURE)
			//	fprintf( stderr, "%i ---> dx = %e, x=%e\n", ii, dx[index], x[index]);
		}
	}

	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else
			this->vesselNodes[i]->marker += dt * dx[i];
	}

	for( int i=M; i<M+N; i++){
		x[i] += dt * dx[i];
	}
}
//#define UNIDIRECTIONAL_TRANSPORT
void VesselGraph::updateMarkerVesselsAndInterstitialExplicit( float dt, float border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, bool UpdateInterstitialCompartment, bool UpdatePlasmaCompartment)
{
	float hepatocytePermeability=10;//10;//2;

	// M Vessel Nodes
	int M = this->countVesselNodes;

	// N Subvessel Segments
	int N = this->countVesselSegments * countVesselSegmentSubNodes;

	// K Voxels
	int K = DOMAIN_SIZE_X * DOMAIN_SIZE_Y * DOMAIN_SIZE_Z;

	if(du==0){
		du = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M + N + K));
#pragma omp parallel for
		for( int i=0; i<M+N+K; i++)
			du[i]=0;
	}
	if(u==0){
		u = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M + N + K));
		for( int i=0; i<M; i++)
			u[i] = this->vesselNodes[i]->marker;
#pragma omp parallel for
		for( int i=M; i<M+N+K; i++)
			u[i]=0;
	}

	//fprintf( stderr, "updateMarkerVesselsAndInterstitialExplicit (M:%i+N:%i+K:%i)\n",M,N,K);

	// INIT MARKER CONC.
	/*for( int i=0; i<M; i++){
		u[i] = this->vesselNodes[i]->marker;
	}*/
#pragma omp parallel for
	for( int i=M+N; i<M+N+K; i++){
		du[i]=0;
	}


	// Construct Matrix for Vascular Marker Concentration: vessel nodes
//#pragma omp parallel for
	for( int i=0; i<M; i++)
	{
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){

			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);
			//du[i] = (border - this->vesselNodes[i]->marker) / dt;
			du[i] = (border - this->vesselNodes[i]->marker) / dt;


		}else{

			// FLOW
			du[i] = 0.;
			CONCENTRATION_T inFlow = 0.;
			CONCENTRATION_T outFlow = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){

					CONCENTRATION_T neighborConc;
					if( countVesselSegmentSubNodes)
						neighborConc = u[M + (1+this->vesselNodes[i]->branches[v]->index)*countVesselSegmentSubNodes-1];
					else
						neighborConc = this->vesselNodes[i]->neighbors[v]->marker;

					inFlow  += fabs(this->vesselNodes[i]->branches[v]->flow) * neighborConc;
					outFlow += fabs(this->vesselNodes[i]->branches[v]->flow) * this->vesselNodes[i]->marker;
				}
			}
			du[i] = (inFlow - outFlow) / (this->getVascularVolume( this->vesselNodes[i]) / (countVesselSegmentSubNodes+1.));

			// DIFFUSION
			int iVoxel = M + N
				+ (int)floor(vesselNodes[i]->position[0]) * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
				+ (int)floor(vesselNodes[i]->position[1]) * DOMAIN_SIZE_Z
				+ (int)floor(vesselNodes[i]->position[2]);
			if( iVoxel >= M+N+K){
				//int M = this->countVesselNodes;
				// N Subvessel Segments
				//int N = this->countVesselSegments * countVesselSegmentSubNodes;

				// K Voxels
				//int K = DOMAIN_SIZE_X * DOMAIN_SIZE_Y * DOMAIN_SIZE_Z;

				fprintf( stderr, "domain=%ix%ix%i\n", DOMAIN_SIZE_X,DOMAIN_SIZE_Y,DOMAIN_SIZE_Z);
				fprintf( stderr, "vessel nodes=%i\n", M);
				fprintf( stderr, "sub nodes=%i\n", N);
				fprintf( stderr, "Voxels=%i\n", K);
				fprintf( stderr, "voxel index=%i\n", iVoxel);
				fprintf( stderr, "Voxel pos=(%lf, %lf, %lf)\n", vesselNodes[i]->position[0],vesselNodes[i]->position[1],vesselNodes[i]->position[2]);
			}
			assert( iVoxel < M+N+K);
			//fprintf( stderr,"(%i,%i,%i)\n", (int)vesselNodes[i]->position[0],(int)vesselNodes[i]->position[1],(int)vesselNodes[i]->position[2]);
			/*if( cells[(int)vesselNodes[i]->position[0]][(int)vesselNodes[i]->position[1]][(int)vesselNodes[i]->position[2]]){
				fprintf( stderr, "vessel in direct contact with cell\n");
//			}*/
#ifdef UNIDIRECTIONAL_DIFFUSION
			if( u[i] > u[iVoxel] || !cells[(int)vesselNodes[i]->position[0]][(int)vesselNodes[i]->position[1]][(int)vesselNodes[i]->position[2]])
#endif
			{
				du[i]      -= ( u[i] - u[iVoxel] ) * getVascularPermeabilitySurfaceProduct( vesselNodes[i])
						/ getVascularVolume( vesselNodes[i]);
				du[iVoxel] += ( u[i] - u[iVoxel] ) * getVascularPermeabilitySurfaceProduct( vesselNodes[i])
						/ getExtraVascularVolume(vesselNodes[i]) / (countVesselSegmentSubNodes+1.);
				//if( du[i]>0. || du[iVoxel]>0.)
				//if(getVascularPermeabilitySurfaceProduct( vesselNodes[i])>0. && u[i]>0.)
				//fprintf( stderr, "UPDATE VOXEL: du[i]=%e, du[iVoxel]=%e (u[i]=%e, u[iVoxel]=%e)\n", du[i],du[iVoxel], u[i],u[iVoxel]);

			}

		}
	}

	// Construct Matrix for Vascular Marker Concentration: vessel segment sub nodes
//#pragma omp parallel for
	if( countVesselSegmentSubNodes)
	for( int iVS=0; iVS<this->countVesselSegments; iVS++){
		//fprintf( stderr, "%i == %i?\n", i, this->vesselSegments[i]->index);
		VesselNode *sourceVN, *targetVN;
		if( this->vesselSegments[iVS]->vesselNodes[0]->pressure > this->vesselSegments[iVS]->vesselNodes[1]->pressure){
			sourceVN = vesselSegments[iVS]->vesselNodes[0];
			targetVN = vesselSegments[iVS]->vesselNodes[1];
		}else{
			sourceVN = vesselSegments[iVS]->vesselNodes[1];
			targetVN = vesselSegments[iVS]->vesselNodes[0];
		}

		// FLOW
		CONCENTRATION_T flow = fabs(this->vesselSegments[iVS]->flow) / (this->getVascularVolume( this->vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.));
		int iSubVS = M + iVS*countVesselSegmentSubNodes;
		du[iSubVS] = flow*( sourceVN->marker - u[iSubVS] );
		for( int ii=1; ii<countVesselSegmentSubNodes; ii++)
		{
			iSubVS++;
			du[iSubVS] = flow*( u[iSubVS-1] - u[iSubVS] );
		}

		// DIFFUSION
		// index voxel
		int iVoxelSource, iVoxelTarget;
		iVoxelSource = M + N
					+ sourceVN->position[0] * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
					+ sourceVN->position[1] * DOMAIN_SIZE_Z
					+ sourceVN->position[2];
		iVoxelTarget = M + N
					+ targetVN->position[0] * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
					+ targetVN->position[1] * DOMAIN_SIZE_Z
					+ targetVN->position[2];
		assert( iVoxelSource < M+N+K);
		assert( iVoxelTarget < M+N+K);


		// update VS Sub Node Contribution
		int ii=0;
		for( ; ii<countVesselSegmentSubNodes/2; ii++)
#ifdef UNIDIRECTIONAL_DIFFUSION
			if( u[iSubVS] > u[iVoxelSource] || !cells[(int)sourceVN->position[0]][(int)sourceVN->position[1]][(int)sourceVN->position[2]])
#endif
		{
			du[iSubVS]       -= ( u[iSubVS] - u[iVoxelSource] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS])
					/ getVascularVolume(vesselSegments[iVS]);
			du[iVoxelSource] += ( u[iSubVS] - u[iVoxelSource] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.)
					/ getExtraVascularVolume(sourceVN);
		}

		for( ; ii<countVesselSegmentSubNodes; ii++)
#ifdef UNIDIRECTIONAL_DIFFUSION
		if( u[iSubVS] > u[iVoxelTarget] || !cells[(int)targetVN->position[0]][(int)targetVN->position[1]][(int)targetVN->position[2]])
#endif
		{
			du[iSubVS]       -= ( u[iSubVS] - u[iVoxelTarget] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS])
					/ getVascularVolume(vesselSegments[iVS]);
			du[iVoxelTarget] += ( u[iSubVS] - u[iVoxelTarget] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.)
					/ getExtraVascularVolume(targetVN);
		}
	}

	// Diffusion between voxels
	//int i=M+N;
#pragma omp parallel for
	for( int x=0; x<DOMAIN_SIZE_X; x++){
		int i=M+N + x*DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int z=0; z<DOMAIN_SIZE_Z; z++)
			{
				CONCENTRATION_T k=0;

				for( int direction=0; direction<6; direction++)
				{
					int X=x,Y=y,Z=z;
					int di=0;
					bool neighborInsideDomain=true;
					switch( direction){
					case 0: if( x>0)               	{X--; di=-DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;}	else neighborInsideDomain=false; break;
					case 1: if( x<DOMAIN_SIZE_X-1) 	{X++; di=+DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;}	else neighborInsideDomain=false; break;
					case 2: if( y>0)               	{Y--; di=-DOMAIN_SIZE_Z;}				else neighborInsideDomain=false; break;
					case 3: if( y<DOMAIN_SIZE_Y-1) 	{Y++; di=+DOMAIN_SIZE_Z;}				else neighborInsideDomain=false; break;
					case 4: if( z>0)               	{Z--; di=-1;}							else neighborInsideDomain=false; break;
					case 5: if( z<DOMAIN_SIZE_Z-1) 	{Z++; di=+1;}							else neighborInsideDomain=false; break;
					}
					if(neighborInsideDomain){
//#ifdef UNIDIRECTIONAL_DIFFUSION
//#endif
#ifdef UNIDIRECTIONAL_TRANSPORT
						if( cells[x][y][z] == cells[X][Y][Z]) /* interstitial diffusion */
						{
							k = (u[i+di]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							du[i] += k;
						}
						if( cells[x][y][z] == 0 && cells[X][Y][Z] != 0)       /* osmotic outflow */
						{
//							k = (-u[i])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							k = (-u[i])*hepatocytePermeability /LATTICE_CONSTANT;
							du[i] += k;
						}
						if(cells[X][Y][Z] == 0 && cells[x][y][z] != 0)       /* osmotic inflow */
						{
//							k = (u[i+di])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							k = (u[i+di])*hepatocytePermeability /LATTICE_CONSTANT;
							du[i] += k;
						}
#else
//						k = (u[i+di]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
//						if( (cells[x][y][z] == cells[X][Y][Z]) /* interstitial diffusion */ ||
//							(cells[x][y][z] == 0 && k<0)       /* osmotic outflow */ ||
//							(cells[X][Y][Z] == 0 && k>0)       /* osmotic inflow */ )
//							du[i] += k;

						k = (u[i+di]-u[i])/LATTICE_CONSTANT;
						if( cells==0 || (cells[x][y][z] == cells[X][Y][Z])) /* diffusion */
						{
							du[i] += k * diffusionCoefficient / LATTICE_CONSTANT * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2.;
						}else
						if(	(cells[x][y][z] == 0 && k<0)       /* osmotic outflow */ ||
							(cells[X][Y][Z] == 0 && k>0)       /* osmotic inflow */ )
							du[i] += k * hepatocytePermeability;
#endif
					}
				}

				/*if( x>0){
					k=(u[i-DOMAIN_SIZE_Z*DOMAIN_SIZE_Y]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x-1,y,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x-1][y][z] || (!cells[x][y][z]&&k<0) || (!cells[x-1][y][z]&&k>0) )du[i] += k;
				}

				if( x<DOMAIN_SIZE_X-1){
					k= (u[i+DOMAIN_SIZE_Z*DOMAIN_SIZE_Y]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x+1,y,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x+1][y][z] || (!cells[x][y][z]&&k<0) || (!cells[x+1][y][z]&&k>0) )du[i] += k;
				}

				if( y>0){
					k= (u[i-DOMAIN_SIZE_Z]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y-1,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x][y-1][z] || (!cells[x][y][z]&&k<0) || (!cells[x][y-1][z]&&k>0) )du[i] += k;
				}
				if( y<DOMAIN_SIZE_Y-1){
					k= (u[i+DOMAIN_SIZE_Z]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y+1,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x][y+1][z] || (!cells[x][y][z]&&k<0) || (!cells[x][y+1][z]&&k>0) )du[i] += k;
				}

				if( z>0){
					k= (u[i-1]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y,z-1)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x][y][z-1] || (!cells[x][y][z]&&k<0) || (!cells[x][y][z-1]&&k>0) )du[i] += k;
				}
				if( z<DOMAIN_SIZE_Z-1){
					k= (u[i+1]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y,z+1)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					if( cells[x][y][z] == cells[x][y][z+1] || (!cells[x][y][z]&&k<0) || (!cells[x][y][z+1]&&k>0) )du[i] += k;
				}*/

				/*if( x>0              && cells[x][y][z] == cells[x-1][y][z])
					if( cells[x][y][z]==0 && cells[x-1][y][z]!=0)
						du[i] += MIN(0, u[i-DOMAIN_SIZE_Z*DOMAIN_SIZE_Y]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x-1,y,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
					else if( cells[x][y][z]==0 && cells[x-1][y][z]!=0)
				}if( x<DOMAIN_SIZE_X-1 && cells[x][y][z] == cells[x+1][y][z])
					du[i] += (u[i+DOMAIN_SIZE_Z*DOMAIN_SIZE_Y]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x+1,y,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;

				if( y>0               && cells[x][y][z] == cells[x][y-1][z])
					du[i] += (u[i-DOMAIN_SIZE_Z]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y-1,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
				if( y<DOMAIN_SIZE_Y-1 && cells[x][y][z] == cells[x][y+1][z])
					du[i] += (u[i+DOMAIN_SIZE_Z]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y+1,z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;

				if( z>0               && cells[x][y][z] == cells[x][y][z-1])
					du[i] += (u[i-1]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y,z-1)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
				if( z<DOMAIN_SIZE_Z-1 && cells[x][y][z] == cells[x][y][z+1])
					du[i] += (u[i+1]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolumeNoRef(x,y,z+1)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
				*/
				i++;
			}
	}

	if( UpdatePlasmaCompartment)
#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else{
			this->vesselNodes[i]->marker += dt * du[i];
			u[i] += dt * du[i];
		}
	}

	if( UpdatePlasmaCompartment)
#pragma omp parallel for
	for( int i=M; i<M+N; i++){
		u[i] += dt * du[i];
	}

	 //i=M+N;
	if(UpdateInterstitialCompartment)
#pragma omp parallel for
	for( int x=0; x<DOMAIN_SIZE_X; x++){
		int i=M+N + x*DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int z=0; z<DOMAIN_SIZE_Z; z++)

			{
				//fprintf(stderr, "TEST! du = %e\n", du[i]);
				/*if( du[i] != 0){
					fprintf(stderr, "Update voxel (%i,%i,%i) = %e\n",x,y,z,du[i]); //exit(0);
				}*/

				u[i] += dt * du[i];
				/*if( u[i] > border)
				{
					fprintf( stderr, "ERROR: u[%i] = %e (du[%i]=%e), border=%e\n", i, u[i], i, du[i], border);
					exit(0);
				}*/
				/*if( u[i] < 0. || u[i] > 7.)
				{
					fprintf( stderr, "ERROR: u[%i] = %e (du[%i]=%e)\n", i, u[i], i, du[i]);
					exit(0);
				}*/
				/*VesselNode *vn = this->octree->at(x,y,z);
				if( vn){
				int iVoxel = vn->index;

				u[i] = (u[i] - u[iVoxel])
						* exp( - this->getVascularPermeabilitySurfaceProduct( x,y,z) * dt / getExtraVascularVolume(x,y,z)) + u[iVoxel];
				}*/
				markerIntSpace[x][y][z] = u[i];
				i++;
			}
	}
}



void VesselGraph::updateMarkerVesselsAndInterstitialExplicit( float dt, float border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, VoronoiDiagram *t)
{
	float hepatocytePermeability=10;//10;//2;

	// =============== UNKNOWNS =====================

	// M Vessel Nodes
	int M = this->countVesselNodes;

	// N Subvessel Segments
	int N = this->countVesselSegments * countVesselSegmentSubNodes;

	// K Voxels
	int K = t->numberOfVertices();//DOMAIN_SIZE_X * DOMAIN_SIZE_Y * DOMAIN_SIZE_Z;


	// =============== ALLOCATE MEMORY for UNKNOWNS =====================

	if(du==0){
		du = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M + N + K));
		#pragma omp parallel for
		for( int i=0; i<M+N+K; i++)
			du[i]=0;
	}
	if(u==0){
		u = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M + N + K));
		#pragma omp parallel for
		for( int i=0; i<M+N+K; i++)
			u[i]=0;
	}


	// =============== RESET UNKNOWNS =====================

	#pragma omp parallel for
	for( int i=0; i<M+N+K; i++){
		du[i]=0;
	}


	// Construct Matrix for Vascular Marker Concentration: vessel nodes
//#pragma omp parallel for
	for( int i=0; i<M; i++)
	{
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE){

			// BOUNDARY: ROOT
			//fprintf( stderr, "Found Root with pressure=%f\n", this->vesselNodes[i]->pressure);
			du[i] = (border - this->vesselNodes[i]->marker) / dt;

		}else{

			// FLOW
			du[i] = 0.;
			CONCENTRATION_T inFlow = 0.;
			CONCENTRATION_T outFlow = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){

					CONCENTRATION_T neighborConc;
					if( countVesselSegmentSubNodes)
						neighborConc = u[M + (1+this->vesselNodes[i]->branches[v]->index)*countVesselSegmentSubNodes-1];
					else
						neighborConc = this->vesselNodes[i]->neighbors[v]->marker;

					inFlow  += fabs(this->vesselNodes[i]->branches[v]->flow) * neighborConc;
					outFlow += fabs(this->vesselNodes[i]->branches[v]->flow) * this->vesselNodes[i]->marker;
				}
			}
			du[i] = (inFlow - outFlow) / (this->getVascularVolume( this->vesselNodes[i]) / (countVesselSegmentSubNodes+1.));

			// DIFFUSION
			/*int iVoxel = M + N
				+ (int)floor(vesselNodes[i]->position[0]) * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
				+ (int)floor(vesselNodes[i]->position[1]) * DOMAIN_SIZE_Z
				+ (int)floor(vesselNodes[i]->position[2]);
			assert( iVoxel < M+N+K);
			if( cells[(int)vesselNodes[i]->position[0]][(int)vesselNodes[i]->position[1]][(int)vesselNodes[i]->position[2]]){
				fprintf( stderr, "vessel in direct contact with cell\n");
			}
#ifdef UNIDIRECTIONAL_DIFFUSION
			if( u[i] > u[iVoxel] || !cells[(int)vesselNodes[i]->position[0]][(int)vesselNodes[i]->position[1]][(int)vesselNodes[i]->position[2]])
#endif
			{
				du[i]      -= ( u[i] - u[iVoxel] ) * getVascularPermeabilitySurfaceProduct( vesselNodes[i])
						/ getVascularVolume( vesselNodes[i]);
				du[iVoxel] += ( u[i] - u[iVoxel] ) * getVascularPermeabilitySurfaceProduct( vesselNodes[i])
						/ getExtraVascularVolume(vesselNodes[i]) / (countVesselSegmentSubNodes+1.);
				//if( du[i]>0. || du[iVoxel]>0.)
				//if(getVascularPermeabilitySurfaceProduct( vesselNodes[i])>0. && u[i]>0.)
				//fprintf( stderr, "UPDATE VOXEL: du[i]=%e, du[iVoxel]=%e (u[i]=%e, u[iVoxel]=%e)\n", du[i],du[iVoxel], u[i],u[iVoxel]);

			}*/

		}
	}

	// Construct Matrix for Vascular Marker Concentration: vessel segment sub nodes
//#pragma omp parallel for
	if( countVesselSegmentSubNodes)
	for( int iVS=0; iVS<this->countVesselSegments; iVS++){
		//fprintf( stderr, "%i == %i?\n", i, this->vesselSegments[i]->index);
		VesselNode *sourceVN, *targetVN;
		if( this->vesselSegments[iVS]->vesselNodes[0]->pressure > this->vesselSegments[iVS]->vesselNodes[1]->pressure){
			sourceVN = vesselSegments[iVS]->vesselNodes[0];
			targetVN = vesselSegments[iVS]->vesselNodes[1];
		}else{
			sourceVN = vesselSegments[iVS]->vesselNodes[1];
			targetVN = vesselSegments[iVS]->vesselNodes[0];
		}

		// FLOW
		CONCENTRATION_T flow = fabs(this->vesselSegments[iVS]->flow) / (this->getVascularVolume( this->vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.));
		int iSubVS = M + iVS*countVesselSegmentSubNodes;
		du[iSubVS] = flow*( sourceVN->marker - u[iSubVS] );
		for( int ii=1; ii<countVesselSegmentSubNodes; ii++)
		{
			iSubVS++;
			du[iSubVS] = flow*( u[iSubVS-1] - u[iSubVS] );
		}

		// DIFFUSION
		// index voxel
/*		int iVoxelSource, iVoxelTarget;
		iVoxelSource = M + N
					+ sourceVN->position[0] * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
					+ sourceVN->position[1] * DOMAIN_SIZE_Z
					+ sourceVN->position[2];
		iVoxelTarget = M + N
					+ targetVN->position[0] * DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
					+ targetVN->position[1] * DOMAIN_SIZE_Z
					+ targetVN->position[2];
		assert( iVoxelSource < M+N+K);
		assert( iVoxelTarget < M+N+K);


		// update VS Sub Node Contribution
		int ii=0;
		for( ; ii<countVesselSegmentSubNodes/2; ii++)
#ifdef UNIDIRECTIONAL_DIFFUSION
			if( u[iSubVS] > u[iVoxelSource] || !cells[(int)sourceVN->position[0]][(int)sourceVN->position[1]][(int)sourceVN->position[2]])
#endif
		{
			du[iSubVS]       -= ( u[iSubVS] - u[iVoxelSource] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS])
					/ getVascularVolume(vesselSegments[iVS]);
			du[iVoxelSource] += ( u[iSubVS] - u[iVoxelSource] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.)
					/ getExtraVascularVolume(sourceVN);
		}

		for( ; ii<countVesselSegmentSubNodes; ii++)
#ifdef UNIDIRECTIONAL_DIFFUSION
		if( u[iSubVS] > u[iVoxelTarget] || !cells[(int)targetVN->position[0]][(int)targetVN->position[1]][(int)targetVN->position[2]])
#endif
		{
			du[iSubVS]       -= ( u[iSubVS] - u[iVoxelTarget] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS])
					/ getVascularVolume(vesselSegments[iVS]);
			du[iVoxelTarget] += ( u[iSubVS] - u[iVoxelTarget] )
					* getVascularPermeabilitySurfaceProduct( vesselSegments[iVS]) / (countVesselSegmentSubNodes+1.)
					/ getExtraVascularVolume(targetVN);
		}
*/
	}

	//fprintf(stderr, "t->LATTICE_CONSTANT=%e \n t->numberOfVertices()=%i\n", t->LATTICE_CONSTANT, t->numberOfVertices());
	int margin = 1;
	Domain<3> *domain = new Domain<3>(margin, (DOMAIN_SIZE_X-margin), margin, (DOMAIN_SIZE_Y-margin), margin, (DOMAIN_SIZE_Z-margin));

	for( int iVC=0; iVC<t->numberOfVertices(); iVC++)
	if(domain->contains(t->get(iVC))){
		// INTERSTITIAL VOLUME
		float interstitialVolume = (t->get(iVC)->volume * pow(t->LATTICE_CONSTANT,3) - this->getVascularVolume( t->get(iVC)->vn));
		if( interstitialVolume < 0){
			//fprintf( stderr, "WARNING: Set negative volume (%e) to new value of %e\n", interstitialVolume, 0.1 * t->get(iVC)->volume * pow(LATTICE_CONSTANT,3));
			interstitialVolume = 0.1 * t->get(iVC)->volume * pow(t->LATTICE_CONSTANT,3);
		}
		//fprintf(stderr, "bla\n");
		// MEMBRAN DIFFUSION

		if( t->get(iVC)->vn!=0){
			//TODO: do diffusion
			/*fprintf( stderr, "==>  Voronoi cell volume=%lf,   vascular volume=%lf, space of disse=%lf\n",
					t->get(iVC)->volume * pow(LATTICE_CONSTANT,3),
					this->getVascularVolume( t->get(iVC)->vn),
					t->get(iVC)->volume * pow(LATTICE_CONSTANT,3) - this->getVascularVolume( t->get(iVC)->vn));
			 */
			int iVN = t->get(iVC)->vn->index;

			//( u[i] - u[iVoxel] ) * getVascularPermeabilitySurfaceProduct( vesselNodes[i]) / getVascularVolume( vesselNodes[i])
			float flux = ( u[iVN] - u[M+N+iVC] ) * getVascularPermeabilitySurfaceProduct( t->get(iVC)->vn);

			//if(flux!=0)
			//	fprintf(stderr, "flux=%e\n", flux);

			// update vessel node
			du[iVN] -= flux / getVascularVolume( t->get(iVC)->vn);

			// update Voronoi cell
			du[M+N+iVC] = flux / interstitialVolume;


		}

		// FREE DIFFUSION
		{
			//fprintf( stderr, "==>  Voronoi cell volume=%lf\n",
			//		t->get(iVC)->volume);
			for( int j=0; j<t->get(iVC)->numberOfNeighbors(); j++)
			if(domain->contains(t->get(iVC)->getNeighbor(j))){
				int jVC = t->get(iVC)->getNeighbor(j)->getIndex();
				// [D] = m^2/s,
				// m^2/s * m^2 / (m * m^3)
				//float interstitialVolumeNeighbor = (t->get(iVC)->volume * pow(LATTICE_CONSTANT,3) - this->getVascularVolume( t->get(iVC)->vn));

				if(jVC>0){ // no framepoint
				du[M+N+iVC] -= ( u[M+N+iVC] - u[M+N+jVC] )
						* this->diffusionCoefficient * (t->get(iVC)->interfaces[j]*t->LATTICE_CONSTANT*t->LATTICE_CONSTANT)
						/ (t->get(iVC)->getDistanceTo( t->get(jVC))*t->LATTICE_CONSTANT) / interstitialVolume;
					/*fprintf( stderr, "==> %e = (%e, %e) * %e * %e / (%e * %e) \n", ( u[M+N+iVC] - u[M+N+jVC] )
							* this->diffusionCoefficient * t->get(iVC)->interfaces[j]
							/ t->get(iVC)->getDistanceTo( t->get(jVC)) / interstitialVolume,
							 u[M+N+iVC] , u[M+N+jVC] ,
							 this->diffusionCoefficient, t->get(iVC)->interfaces[j],
							 t->get(iVC)->getDistanceTo( t->get(jVC)) , interstitialVolume);*/
				}
			}

		}


	}

	delete domain;
/*
	// Diffusion between voxels
	//int i=M+N;
#pragma omp parallel for
	for( int x=0; x<DOMAIN_SIZE_X; x++){
		int i=M+N + x*DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int z=0; z<DOMAIN_SIZE_Z; z++)
			{
				CONCENTRATION_T k=0;

				for( int direction=0; direction<6; direction++)
				{
					int X=x,Y=y,Z=z;
					int di=0;
					bool neighborInsideDomain=true;
					switch( direction){
					case 0: if( x>0)               	{X--; di=-DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;}	else neighborInsideDomain=false; break;
					case 1: if( x<DOMAIN_SIZE_X-1) 	{X++; di=+DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;}	else neighborInsideDomain=false; break;
					case 2: if( y>0)               	{Y--; di=-DOMAIN_SIZE_Z;}				else neighborInsideDomain=false; break;
					case 3: if( y<DOMAIN_SIZE_Y-1) 	{Y++; di=+DOMAIN_SIZE_Z;}				else neighborInsideDomain=false; break;
					case 4: if( z>0)               	{Z--; di=-1;}							else neighborInsideDomain=false; break;
					case 5: if( z<DOMAIN_SIZE_Z-1) 	{Z++; di=+1;}							else neighborInsideDomain=false; break;
					}
					if(neighborInsideDomain){
//#ifdef UNIDIRECTIONAL_DIFFUSION
//#endif
#ifdef UNIDIRECTIONAL_TRANSPORT
						if( cells[x][y][z] == cells[X][Y][Z]) /* interstitial diffusion * /
						{
							k = (u[i+di]-u[i])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							du[i] += k;
						}
						if( cells[x][y][z] == 0 && cells[X][Y][Z] != 0)       /* osmotic outflow * /
						{
//							k = (-u[i])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							k = (-u[i])*hepatocytePermeability /LATTICE_CONSTANT;
							du[i] += k;
						}
						if(cells[X][Y][Z] == 0 && cells[x][y][z] != 0)       /* osmotic inflow * /
						{
//							k = (u[i+di])*diffusionCoefficient * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2. /LATTICE_CONSTANT/LATTICE_CONSTANT;
							k = (u[i+di])*hepatocytePermeability /LATTICE_CONSTANT;
							du[i] += k;
						}
#else

						k = (u[i+di]-u[i])/LATTICE_CONSTANT;
						if( (cells[x][y][z] == cells[X][Y][Z])) /* diffusion * /
						{
							du[i] += k * diffusionCoefficient / LATTICE_CONSTANT * (1.+getExtraVascularVolume(X,Y,Z)/getExtraVascularVolume(x,y,z))/2.;
						}else
						if(	(cells[x][y][z] == 0 && k<0)       /* osmotic outflow * / ||
							(cells[X][Y][Z] == 0 && k>0)       /* osmotic inflow * / )
							du[i] += k * hepatocytePermeability;
#endif
					}
				}

				i++;
			}
	}
*/

	// UPDATE

#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->getType() == ROOT && this->vesselNodes[i]->pressure==MAX_PRESSURE)
			this->vesselNodes[i]->marker = border;
		else{
			this->vesselNodes[i]->marker += dt * du[i];
			u[i] += dt * du[i];
		}
	}

#pragma omp parallel for
	for( int i=M; i<M+N; i++){
		u[i] += dt * du[i];
	}

	 //i=M+N;
#pragma omp parallel for
	/*for( int x=0; x<DOMAIN_SIZE_X; x++){
		int i=M+N + x*DOMAIN_SIZE_Z*DOMAIN_SIZE_Y;
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int z=0; z<DOMAIN_SIZE_Z; z++)
			{
				u[i] += dt * du[i];
				markerIntSpace[x][y][z] = u[i];
				i++;
			}
	}*/

	for( int iVC=0; iVC<t->numberOfVertices(); iVC++){
		/*int i=(int)t->get(iVC)->x()*DOMAIN_SIZE_Z*DOMAIN_SIZE_Y
			 +(int)t->get(iVC)->y()*DOMAIN_SIZE_Z
			 +(int)t->get(iVC)->z();*/
		u[M+N+iVC] += dt * du[M+N+iVC];
		//markerIntSpace[(int)t->get(iVC)->x()][(int)t->get(iVC)->y()][(int)t->get(iVC)->z()] = u[M+N+iVC];
		t->get(iVC)->conc = u[M+N+iVC];
	}
}




float VesselGraph::getExtraVascularVolume( VesselNode *vn)
{
	// volume: tissue
	float volume = LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT;
	if( vn!=0)
	volume -= this->getVascularVolume( vn);

	float minVolumeFraction = 0.01;//3;//0.001;
	if( volume<minVolumeFraction*LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT){
		volume = minVolumeFraction*LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT;
		//fprintf( stderr, "ERROR: extra vascular volume is negative (%lf)!!\n", volume);
		//exit(0);
	}

	return volume;
}

float VesselGraph::getExtraVascularVolume( int &x, int &y, int &z)
{
	// volume: tissue
	VesselNode *vn = this->octree->at( x,y,z);
	return this->getExtraVascularVolume( vn);;
}

float VesselGraph::getExtraVascularVolumeNoRef( int x, int y, int z)
{
	// volume: tissue
	VesselNode *vn = this->octree->at( x,y,z);
	return this->getExtraVascularVolume( vn);;
}

float VesselGraph::getVascularSurface( VesselSegment *vs)
{
	// volume: vessel segment
	if( vs!=0)
		return this->distance(vs->vesselNodes[0], vs->vesselNodes[1]) * LATTICE_CONSTANT
			* 2. * PI * vs->radius * vs->countParallelVessels;

	return 0;
}

float VesselGraph::getVascularSurface( VesselNode *vn)
{
	// volume: all vessel segments in voxel of tissue
	float surface = 0.;

	if( vn!=0)
	for( int v=0; v<vn->countNeighbors; v++)
		surface += 0.5 * getVascularSurface(vn->branches[v]);

	return surface;
}

float VesselGraph::getVascularPermeabilitySurfaceProduct( int &x, int &y, int &z){
		// volume: all vessel segments in voxel of tissue
		VesselNode *vn = this->octree->at( x,y,z);
		return getVascularPermeabilitySurfaceProduct( vn);
}

float VesselGraph::getVascularPermeabilitySurfaceProduct( VesselSegment *vs)
{
	return getVascularSurface(vs) * vs->permeability;
}

float VesselGraph::getVascularPermeabilitySurfaceProduct( VesselNode *vn)
{
	// volume: vessel segment
	float permeabilitySurfaceProduct = 0.;

	if( vn!=0)
	for( int v=0; v<vn->countNeighbors; v++)
		permeabilitySurfaceProduct += 0.5 * getVascularPermeabilitySurfaceProduct(vn->branches[v]);
	//fprintf(stderr, "PS=%e\n",permeabilitySurfaceProduct);
	return permeabilitySurfaceProduct;
}

float VesselGraph::getVascularVolume( VesselSegment *vs)
{
	// volume: vessel segment
	if( vs!=0){

		//if( this->distance(vs->vesselNodes[0], vs->vesselNodes[1]) == 0)
		//{fprintf( stderr, "ERROR in getVascularVolume(): node distance == 0\n");}
		return this->distance(vs->vesselNodes[0], vs->vesselNodes[1]) * LATTICE_CONSTANT
			* PI * vs->radius * vs->radius * vs->countParallelVessels;
	}

	return 0;
}

float VesselGraph::getVascularVolume( VesselNode *vn)
{
	// volume: all vessel segments in voxel of tissue
	float volume = 0.;
	if( vn!=0)
	for( int v=0; v<vn->countNeighbors; v++)
		volume +=
			0.5	* this->getVascularVolume( vn->branches[v]);

	if( volume>LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT){
		//fprintf( stderr, "ERROR: extra vascular volume is larger than voxel size (%lf < %lf)!!\n", volume, LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT);
		//exit(0);
	}

	//return 60*4*4*PI;
	return volume;
}

float VesselGraph::getVascularVolume( int &x, int &y, int &z)
{
	// volume: all vessel segments in voxel of tissue
	VesselNode *vn = this->octree->at( x,y,z);
	return this->getVascularVolume( vn);
}

float VesselNode::getFlowFromNeighbor( int i)
{
	assert(i<this->countNeighbors);

	return (this->neighbors[i]->pressure >= this->pressure ?

			// IN FLOW defined POSITIVE
			+ fabs(this->branches[i]->flow) :

			// OUT FLOW defined NEGATIVE
			- fabs(this->branches[i]->flow));
}

float VesselNode::getFlowToNeighbor( int i)
{
	assert(i<this->countNeighbors);

	return (this->pressure >= this->neighbors[i]->pressure ?

			// OUT FLOW defined POSITIVE
			+ fabs(this->branches[i]->flow) :

			// IN FLOW defined NEGATIVE
			- fabs(this->branches[i]->flow));
}

void VesselGraph::updateMarkerVesselsAndInterstitialSpaceNEW3( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut)
{
	float *b=0;
	float *x=0;
	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;

	SparseMatrix<float> *sA = new SparseMatrix<float>(M+N,M+N);
	b = (float*) malloc( sizeof(float) * (M+N));
	x = (float*) malloc( sizeof(float) * (M+N));

	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB );


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}
	for( int i=M; i<M+N; i++){
		x[i] = 0;
	}

	// Construct Matrix for Vascular Marker Concentration
	//float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				//b[i] = ARTERIAL_MARKER_CONC;
				b[i] = border;
			else
				b[i] = 0;
		}else{
			//float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			//float sumFlow=0.;

			// volume
			float volume = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				volume += 0.5*LATTICE_CONSTANT
						* PI * this->vesselNodes[i]->branches[v]->radius*this->vesselNodes[i]->branches[v]->radius;


			// c_P^n+1
			sA->set(i,i, 1);

			// r
			float r = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
					r += fabs(this->vesselNodes[i]->branches[v]->flow);
			}
			r = 1. + dt/volume*r + dt/volume*EXCHANGE_RATE_VESSELS_INTERSTITIAL;

			// c_P,i^n+1 (for all neighbors i)
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow*dt/r/volume);
					if( isnan(flow*dt/r/volume)){
						fprintf( stderr, "flow=%lf dt=%lf r=%lf volume=%lf\n",flow,dt,r,volume);
						exit(0);
					}
				}
			}

			// c_I^n+1
			int index = this->vesselNodes[i]->position[0]
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume);
			if( isnan(EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume)){
				fprintf( stderr, "EXCHANGE_RATE_VESSELS_INTERSTITIAL=%lf dt=%lf r=%lf volume=%lf\n",EXCHANGE_RATE_VESSELS_INTERSTITIAL,dt,r,volume);
				exit(0);
			}

			// b
			b[i] = 1./r * this->vesselNodes[i]->marker;
			if( isnan(b[i])){
				fprintf( stderr, "r=%lf this->vesselNodes[i]->marker=%lf\n",r,this->vesselNodes[i]->marker);
				exit(0);
			}

		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;

				// volume: tissue
				VesselNode *vn = this->octree->at( x,y,z);
				float volume = LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT;
				if( vn!=0)
				for( int v=0; v<vn->countNeighbors; v++)
					volume -= 0.5*LATTICE_CONSTANT
						* PI * vn->branches[v]->radius*vn->branches[v]->radius;

				if(volume < 0.){
					fprintf( stderr, "ERROR: negative volume: %lf!!\n", volume);
					fprintf( stderr, "ERROR: max volume: %lf = %lf^3\n", (float)LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT, (float)LATTICE_CONSTANT);
					for( int v=0; v<vn->countNeighbors; v++)
						fprintf( stderr, "%i neighbor: radius=%lf, volume=%lf\n", v, vn->branches[v]->radius, 0.5*LATTICE_CONSTANT	* PI * vn->branches[v]->radius*vn->branches[v]->radius);
					exit(0);
				}

				// c_I^n+1
				sA->set(i,i, 1);

				// r
				//float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				float r = 0.;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2){
					r+= 1. + this->getExtraVascularVolume( nx, y, z) / volume;
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2){
					r+= 1. + this->getExtraVascularVolume( x, ny, z) / volume;
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2){
					r+= 1. + this->getExtraVascularVolume( x, y, nz) / volume;
				}
				r = 1.
				  + dt * MARKER_DIFFUSION_COEFFICIENT /(LATTICE_CONSTANT*LATTICE_CONSTANT) * r/2.
				  + dt * EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume;

				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( nx, y, z) / volume)/2.);

					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( x, ny, z) / volume)/2.);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( x, y, nz) / volume)/2.);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}

				// c_P^n+1
				if( vn!=0){
					sA->set(i,vn->index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume * dt/r);
					if( isnan(-EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume * dt/r)){
						fprintf( stderr, "r=%lf volume=%lf\n",r,volume);
						exit(0);
					}
				}

				// b
				b[i] = 1./r * markerIntSpaceIn[x][y][z];
				if( isnan(b[i])){
					fprintf( stderr, "r=%lf markerIntSpaceIn[x][y][z]=%lf\n",r,markerIntSpaceIn[x][y][z]);
					exit(0);
				}
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);

	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	//fprintf( stderr, "Precondition System\n");
	S->PreconditionJacobi( sA, b);
	//fprintf( stderr, "Solve System\n");
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpaceNEW4( float dt, float border, float ***markerIntSpaceIn, float ***markerIntSpaceOut)
{
	float *b=0;
	float *x=0;
	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;

	SparseMatrix<float> *sA = new SparseMatrix<float>(M+N,M+N);
	b = (float*) malloc( sizeof(float) * (M+N));
	x = (float*) malloc( sizeof(float) * (M+N));

	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB );


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}
	for( int i=M; i<M+N; i++){
		x[i] = 0;
	}

	// Construct Matrix for Vascular Marker Concentration
	//float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				//b[i] = ARTERIAL_MARKER_CONC;
				b[i] = border;
			else
				b[i] = 0;
		}else{
			//float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			//float sumFlow=0.;

			// volume
			float volume = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				volume += 0.5*LATTICE_CONSTANT
						* PI * this->vesselNodes[i]->branches[v]->radius*this->vesselNodes[i]->branches[v]->radius;


			// c_P^n+1
			sA->set(i,i, 1);

			// r
			float r = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure)
					r += fabs(this->vesselNodes[i]->branches[v]->flow);
			}
			r = 1. + dt/volume*r + dt/volume * this->getVascularPermeabilitySurfaceProduct( vesselNodes[i]);

			// c_P,i^n+1 (for all neighbors i)
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow*dt/r/volume);
					if( isnan(flow*dt/r/volume)){
						fprintf( stderr, "flow=%lf dt=%lf r=%lf volume=%lf\n",flow,dt,r,volume);
						exit(0);
					}
				}
			}

			// c_I^n+1
			int index = this->vesselNodes[i]->position[0]
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -this->getVascularPermeabilitySurfaceProduct( vesselNodes[i])*dt/r/volume);
			if( isnan(EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume)){
				fprintf( stderr, "EXCHANGE_RATE_VESSELS_INTERSTITIAL=%lf dt=%lf r=%lf volume=%lf\n",EXCHANGE_RATE_VESSELS_INTERSTITIAL,dt,r,volume);
				exit(0);
			}

			// b
			b[i] = 1./r * this->vesselNodes[i]->marker;
			if( isnan(b[i])){
				fprintf( stderr, "r=%lf this->vesselNodes[i]->marker=%lf\n",r,this->vesselNodes[i]->marker);
				exit(0);
			}

		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;

				VesselNode *vn = this->octree->at( x,y,z);

				// volume: tissue
				float volume = this->getExtraVascularVolume(vn);

				// marker exchange
				float exchange_rate = this->getVascularPermeabilitySurfaceProduct( vn);

				// c_I^n+1
				sA->set(i,i, 1);

				// r
				//float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				float r = 0.;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2){
					r+= 1. + this->getExtraVascularVolume( nx, y, z) / volume;
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2){
					r+= 1. + this->getExtraVascularVolume( x, ny, z) / volume;
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2){
					r+= 1. + this->getExtraVascularVolume( x, y, nz) / volume;
				}
				r = 1.
				  + dt * MARKER_DIFFUSION_COEFFICIENT /(LATTICE_CONSTANT*LATTICE_CONSTANT) * r/2.
				  + dt * exchange_rate/volume;

				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( nx, y, z) / volume)/2.);

					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( x, ny, z) / volume)/2.);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r * (1. + this->getExtraVascularVolume( x, y, nz) / volume)/2.);
					if( isnan(MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r)){
						fprintf( stderr, "r=%lf\n",r);
						exit(0);
					}
				}

				// c_P^n+1
				if( vn!=0){
					sA->set(i,vn->index, -exchange_rate/volume * dt/r);
					if( isnan(-exchange_rate/volume * dt/r)){
						fprintf( stderr, "r=%lf volume=%lf\n",r,volume);
						exit(0);
					}
				}

				// b
				b[i] = 1./r * markerIntSpaceIn[x][y][z];
				if( isnan(b[i])){
					fprintf( stderr, "r=%lf markerIntSpaceIn[x][y][z]=%lf\n",r,markerIntSpaceIn[x][y][z]);
					exit(0);
				}
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);

	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	//fprintf( stderr, "Precondition System\n");
	S->PreconditionJacobi( sA, b);
	//fprintf( stderr, "Solve System\n");
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpaceNEW5( float dt, float border, CONCENTRATION_T ***markerIntSpaceIn, CONCENTRATION_T ***markerIntSpaceOut)
{
	CONCENTRATION_T *b=0;
	CONCENTRATION_T *x=0;
	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;

	SparseMatrix<CONCENTRATION_T> *sA = new SparseMatrix<CONCENTRATION_T>(M+N,M+N);
	b = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M+N));
	x = (CONCENTRATION_T*) malloc( sizeof(CONCENTRATION_T) * (M+N));

	Solver<CONCENTRATION_T> *S = new Solver<CONCENTRATION_T>( M+N, Solver<CONCENTRATION_T>::BiCGSTAB );


	// INIT MARKER CONC.
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}
	for( int i=M; i<M+N; i++){
		x[i] = 0;
	}

	// Construct Matrix for Vascular Marker Concentration
#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				b[i] = border;
			else
				b[i] = 0;
		}else{
			// INTRA-Vascular Volume
			CONCENTRATION_T volumePlasma = getVascularVolume( vesselNodes[i]);
			CONCENTRATION_T flow = 0;
			CONCENTRATION_T transport = 0;

			// C_P^n+1   [ PLASMA COMPARTMENT ]
			{
				// SUM of OUTGOING FLOWS (F)
				for( int v=0; v<vesselNodes[i]->countNeighbors; v++)
					if( vesselNodes[i]->pressure >= vesselNodes[i]->neighbors[v]->pressure)
						flow += fabs(vesselNodes[i]->branches[v]->flow);

				// OUTGOING TRANSPORT (Kps)
				transport = getVascularPermeabilitySurfaceProduct( vesselNodes[i]);

				sA->set( i,i, 1. + dt/volumePlasma * (flow + transport));
			}

			// C_P,i^n+1   [ NEIGHBORING PLASMA COMPARTMENTS i ]
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// INCOMING FLOW
				if( vesselNodes[i]->pressure < vesselNodes[i]->neighbors[v]->pressure){
					flow = fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index,
							- dt/volumePlasma * flow);
				}
			}

			// c_I^n+1
			{
				// INCOMING TRANSPORT (Kps)
				int j = M
					  + this->vesselNodes[i]->position[0]
				      + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
				      + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				//if( x[i] > x[j] )
					sA->set( i, j, - dt/volumePlasma * transport );
				//else
				//	sA->set( i, j, 0 );
			}

			// b
			b[i] = this->vesselNodes[i]->marker;
		}
	}

	// Construct Matrix for Interstitial Marker Concentration
#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;

				VesselNode *vn = this->octree->at( x,y,z);

				// INTERSTITIAL Volume
				CONCENTRATION_T volumeInterstitial = this->getExtraVascularVolume( x, y, z);

				// DIFFUSION
				CONCENTRATION_T diffusion = 0.;
				{
					for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2){
						int j = M
							  + nx
							  + y*DOMAIN_SIZE_X
							  + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
						// INCOMING DIFFUSION (D)
						sA->set( i,j, -( volumeInterstitial + getExtraVascularVolume( nx, y, z) ) / 2 * dt/volumeInterstitial * diffusionCoefficient / (LATTICE_CONSTANT*LATTICE_CONSTANT));
						// SUM of OUTGOING DIFFUSION (D)
						diffusion +=  ( volumeInterstitial + getExtraVascularVolume( nx, y, z) ) / 2;
					}
					for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2){
						int j = M
							  + x
							  + ny*DOMAIN_SIZE_X
							  + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
						sA->set( i,j, -( volumeInterstitial + getExtraVascularVolume( x, ny, z) ) / 2 * dt/volumeInterstitial * diffusionCoefficient / (LATTICE_CONSTANT*LATTICE_CONSTANT));
						diffusion  += ( volumeInterstitial + getExtraVascularVolume( x, ny, z) ) / 2;
					}
					for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2){
						int j = M
							  + x
							  + y*DOMAIN_SIZE_X
							  + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
						sA->set( i,j, -( volumeInterstitial + getExtraVascularVolume( x, y, nz) ) / 2 * dt/volumeInterstitial * diffusionCoefficient / (LATTICE_CONSTANT*LATTICE_CONSTANT));
						diffusion  += ( volumeInterstitial + getExtraVascularVolume( x, y, nz) ) / 2;
					}
					diffusion *= diffusionCoefficient / (LATTICE_CONSTANT*LATTICE_CONSTANT);
				}


				CONCENTRATION_T transport = 0;
				if( vn){
					transport = getVascularPermeabilitySurfaceProduct( vn);

					// INCOMING TRANSPORT (Kps)
					int j = (int)(vn->index);
					//if( x[j] > x[i] )
						sA->set( i, vn->index, - dt/volumeInterstitial * transport);
					//else
					//	sA->set( i, vn->index, 0);
				}

				// OUTGOING DIFFUSION + TRANSPORT
				sA->set(i,i, 1 + dt/volumeInterstitial * ( diffusion + transport ) );


				// b
				b[i] = markerIntSpaceIn[x][y][z];
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);

	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	//fprintf( stderr, "Precondition System\n");
	S->PreconditionJacobi( sA, b);
	//fprintf( stderr, "Solve System\n");
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateMarkerVesselsAndInterstitialSpaceBrix( float dt, float **markerVesselsIn, float **markerVesselsOut, float ***markerIntSpaceIn, float ***markerIntSpaceOut)
{
	float *b=0;
	float *x=0;
	int M = this->countVesselNodes;
	int N = DOMAIN_SIZE_X*DOMAIN_SIZE_Y*DOMAIN_SIZE_Z;
	//if(sA==0)
	SparseMatrix<float> *sA = new SparseMatrix<float>(M+N,M+N);
	//if(b==0)
	b = (float*) malloc( sizeof(float) * (M+N));
	//if(x==0)
	x = (float*) malloc( sizeof(float) * (M+N));

	/*v0 = (float*) malloc( sizeof(float) * (M+N));
	v1 = (float*) malloc( sizeof(float) * (M+N));
	v2 = (float*) malloc( sizeof(float) * (M+N));
	v3 = (float*) malloc( sizeof(float) * (M+N));
	v4 = (float*) malloc( sizeof(float) * (M+N));
	v5 = (float*) malloc( sizeof(float) * (M+N));
	v6 = (float*) malloc( sizeof(float) * (M+N));
*/

	// INIT PRESSURE VALUES
	for( int i=0; i<M; i++){
		x[i] = this->vesselNodes[i]->marker;
	}

	// Construct Matrix for Vascular Marker Concentration
	//float volume = 0.5;
//#pragma omp parallel for
	for( int i=0; i<M; i++){
		if( this->vesselNodes[i]->countNeighbors == 0)
			fprintf( stderr, "ERROR: Node without neighbors!\n");
		if( this->vesselNodes[i]->getType() == ROOT){
			// BOUNDARY: ROOT
			sA->set(i,i, 1);
			if( this->vesselNodes[i]->pressure==MAX_PRESSURE)
				b[i] = ARTERIAL_MARKER_CONC;
			else
				b[i] = 0;
		}else{
			//float r=volume/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
			//float sumFlow=0.;

			// volume
			float volume = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
				volume += 0.5*LATTICE_CONSTANT
						* PI * this->vesselNodes[i]->branches[v]->radius*this->vesselNodes[i]->branches[v]->radius;


			// c_P^n+1
			sA->set(i,i, 1);

			// r
			float flowSum = 0.;
			for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				// out flow
				if( this->vesselNodes[i]->pressure >= this->vesselNodes[i]->neighbors[v]->pressure){
					flowSum += fabs(this->vesselNodes[i]->branches[v]->flow);
				}
			}
			float r = 0.;
			r = 1. + dt/volume*flowSum + dt/volume*EXCHANGE_RATE_VESSELS_INTERSTITIAL;

			// c_P,i^n+1 (for all neighbors i)
			/*for( int v=0; v<this->vesselNodes[i]->countNeighbors; v++)
			{
				float flow;
				// in flow
				if( this->vesselNodes[i]->pressure < this->vesselNodes[i]->neighbors[v]->pressure){
					flow = -fabs(this->vesselNodes[i]->branches[v]->flow);
					sA->set(i,this->vesselNodes[i]->neighbors[v]->index, flow*dt/r/volume);
				}
			}
			*/

			// c_I^n+1
			int index = this->vesselNodes[i]->position[0]
			          + this->vesselNodes[i]->position[1]*DOMAIN_SIZE_X
			          + this->vesselNodes[i]->position[2]*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
			sA->set(i,M+index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL*dt/r/volume);

			// b
			//b[i] = 1./r * this->vesselNodes[i]->marker;
			b[i] = 1./r * this->vesselNodes[i]->marker + dt/r * ARTERIAL_MARKER_CONC * flowSum/volume;
		}
	}

	// Construct Matrix for Interstitial Marker Concentration
//#pragma omp parallel for
	for( int z=0; z<DOMAIN_SIZE_Z; z++)
		for( int y=0; y<DOMAIN_SIZE_Y; y++)
			for( int x=0; x<DOMAIN_SIZE_X; x++)
			{
				int i = M
					  + x
				      + y*DOMAIN_SIZE_X
				      + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;

				// volume: tissue
				VesselNode *vn = this->octree->at( x,y,z);
				float volume = LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT;
				if( vn!=0)
				for( int v=0; v<vn->countNeighbors; v++)
					volume -= 0.5*LATTICE_CONSTANT
						* PI * vn->branches[v]->radius*vn->branches[v]->radius;

				if(volume < 0.){
					fprintf( stderr, "ERROR: negative volume: %lf!!\n", volume);
					fprintf( stderr, "ERROR: max volume: %lf = %lf^3\n", (float)LATTICE_CONSTANT*LATTICE_CONSTANT*LATTICE_CONSTANT, (float)LATTICE_CONSTANT);
					for( int v=0; v<vn->countNeighbors; v++)
						fprintf( stderr, "%i neighbor: radius=%lf, volume=%lf\n", v, vn->branches[v]->radius, 0.5*LATTICE_CONSTANT	* PI * vn->branches[v]->radius*vn->branches[v]->radius);
					exit(0);
				}

				// c_I^n+1
				sA->set(i,i, 1);

				// r
				//float r=(1.-volume)/dt + EXCHANGE_RATE_VESSELS_INTERSTITIAL;
				float r = 0.;
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
					r+=1;
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
					r+=1;
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
					r+=1;
				r = 1.
				  + dt * MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * r
				  + dt * EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume;

				// c_I,i^n+1 (for all neighbors i)
				for( int nx=(x>0?x-1:x+1); nx<=x+1 && nx<DOMAIN_SIZE_X; nx+=2)
				{
					int ii = M
						   + nx
					       + y*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
				}
				for( int ny=(y>0?y-1:y+1); ny<=y+1 && ny<DOMAIN_SIZE_Y; ny+=2)
				{
					int ii = M
						   + x
					       + ny*DOMAIN_SIZE_X
					       + z*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
				}
				for( int nz=(z>0?z-1:z+1); nz<=z+1 && nz<DOMAIN_SIZE_Z; nz+=2)
				{
					int ii = M
						   + x
					       + y*DOMAIN_SIZE_X
					       + nz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
					sA->set(i,ii, -MARKER_DIFFUSION_COEFFICIENT/(LATTICE_CONSTANT*LATTICE_CONSTANT) * dt/r);
				}

				// c_P^n+1
				if( vn!=0)
					sA->set(i,vn->index, -EXCHANGE_RATE_VESSELS_INTERSTITIAL/volume * dt/r);

				// b
				b[i] = 1./r * markerIntSpaceIn[x][y][z];
			}


	//SolveBiCGSTAB( sA, b, x, 20000, 1e-10);
	//ConjugateGradientSparse( sA, b, x, 20000, 1e-20);
	Solver<float> *S = new Solver<float>( M+N, Solver<float>::BiCGSTAB);
	S->solve( sA, b, x);
	delete S;

	for( int i=0; i<M; i++)
		this->vesselNodes[i]->marker = x[i];
	for( int iz=0; iz<DOMAIN_SIZE_Z; iz++)
		for( int iy=0; iy<DOMAIN_SIZE_Y; iy++)
			for( int ix=0; ix<DOMAIN_SIZE_X; ix++)
			{
				int i = M
					  + ix
				      + iy*DOMAIN_SIZE_X
				      + iz*DOMAIN_SIZE_X*DOMAIN_SIZE_Y;
				markerIntSpaceOut[ix][iy][iz] = x[i];
			}


	delete sA;
	free(b);
	free(x);
}


void VesselGraph::updateSegmentTypes()
{
	int set = 0;

	// Set Roots
	for(int vn=0; vn<this->countVesselNodes; vn++){
		if(this->vesselNodes[vn]->getType() == ROOT){
			if( this->vesselNodes[vn]->pressure == MAX_PRESSURE){
				for( int b=0; b<this->vesselNodes[vn]->countNeighbors; b++)
					this->vesselNodes[vn]->branches[b]->type = ARTERIE;
			}else{
				for( int b=0; b<this->vesselNodes[vn]->countNeighbors; b++)
					this->vesselNodes[vn]->branches[b]->type = VENE;
			}
		}
	}

	// Remaining
	do{
		set=0;

		for(int vs=0; vs<this->countVesselSegments; vs++){
			for( int b=0; this->vesselSegments[vs]->type == UNKNOWN && b<this->vesselSegments[vs]->vesselNodes[0]->countNeighbors; b++)
				if( this->vesselSegments[vs]->vesselNodes[0]->branches[b]->type != UNKNOWN)
					vesselSegments[vs]->type = vesselSegments[vs]->vesselNodes[0]->branches[b]->type;
			for( int b=0; this->vesselSegments[vs]->type == UNKNOWN && b<this->vesselSegments[vs]->vesselNodes[1]->countNeighbors; b++)
				if( this->vesselSegments[vs]->vesselNodes[1]->branches[b]->type != UNKNOWN)
					vesselSegments[vs]->type = vesselSegments[vs]->vesselNodes[1]->branches[b]->type;

			if( this->vesselSegments[vs]->type != UNKNOWN)
				set++;
		}
	}while( set<this->countVesselSegments);
}

void VesselGraph::updateRadius()
{
	int set = 0;
	/*for( int i=0; i<this->countVesselSegments; i++)
	if( !this->vesselSegments[i]->radiusStatic){
		if( this->vesselSegments[i]->vesselNodes[0]->getType() != TIP && this->vesselSegments[i]->vesselNodes[1]->getType() != TIP)
			this->vesselSegments[i]->radius = 0;
		else{
			this->vesselSegments[i]->radius = MIN_VESSEL_RADIUS;//*LATTICE_CONSTANT/6.;
			set++;
			//fprintf( stderr, "INFO: set segment %i\n",this->vesselSegments[i]->index);
		}
	}*/

	// UPDATE INTER-TIP-CONNECTIONS
	for( int i=0; i<this->countVesselSegments; i++)
		if( !this->vesselSegments[i]->radiusStatic){
			this->vesselSegments[i]->radius = 0;

			// UPDATE INTER-TIP-CONNECTIONS
			if( this->vesselSegments[i]->vesselNodes[0]->getType() == TIP && this->vesselSegments[i]->vesselNodes[1]->getType() == TIP){
				this->vesselSegments[i]->radius = MIN_VESSEL_RADIUS;
				set++;
			}
			// UPDATE DEAD-END-TIPS
			if( this->vesselSegments[i]->vesselNodes[0]->getType() == TIP && this->vesselSegments[i]->vesselNodes[1]->getType() != TIP && this->vesselSegments[i]->vesselNodes[0]->countNeighbors == 1){
				this->vesselSegments[i]->radius = MIN_VESSEL_RADIUS;set++;}
			if( this->vesselSegments[i]->vesselNodes[0]->getType() != TIP && this->vesselSegments[i]->vesselNodes[1]->getType() == TIP && this->vesselSegments[i]->vesselNodes[1]->countNeighbors == 1){
				this->vesselSegments[i]->radius = MIN_VESSEL_RADIUS;set++;}

		}else{
			// KEEP STATIC RADII
			set++;
		}


	//fprintf( stderr, "INFO: set=%i\n",set);
	bool change;
	float alpha = 2.;//2.7;//6.;//2.7;
	do{
		//fprintf( stderr, "\rupdateRadius: %.2f\%% \b", (float)set/(float)this->countVesselSegments*100.);
		change = false;
		int i=0;
//#pragma omp parallel for
		for( i=0; i<this->countVesselNodes; i++){
			if( /*this->vesselNodes[i]->getType() != TIP &&*/ this->vesselNodes[i]->countNeighbors>1/*&& this->vesselNodes[i]->radius == 0*/){
				float radius = 0;
				int countNonSet = 0;
				int index = -1;
				for( int b=0; b<this->vesselNodes[i]->countNeighbors; b++){
					if( this->vesselNodes[i]->branches[b]->radiusStatic){
						// STATIC RADIUS
						countNonSet=2;
					}else
					if( this->vesselNodes[i]->branches[b]->radius != 0){
						// SET ALREADY
						alpha = (vesselNodes[i]->branches[b]->type == ARTERIE ? RADIUS_EXPONENT_ARTERIES : RADIUS_EXPONENT_VENES);
						radius += pow(this->vesselNodes[i]->branches[b]->radius, alpha) * this->vesselNodes[i]->branches[b]->countParallelVessels;
					}else{
						// NOT SET
						index = this->vesselNodes[i]->branches[b]->index;
						countNonSet++;
					}
				}
				if( countNonSet==1){
					//fprintf( stderr, "INFO: Set Radius of Segment %i\n", index);
					alpha = (vesselSegments[index]->type == ARTERIE ? RADIUS_EXPONENT_ARTERIES : RADIUS_EXPONENT_VENES);
					this->vesselSegments[index]->radius =
							MIN( MAX_VESSEL_RADIUS, pow(radius, 1./alpha));
					change = true;
					set++;
					/*fprintf( stderr, "INFO: set segment %i (node %i) -> radius=%lf\n",
							this->vesselNodes[i]->index,
							this->vesselSegments[index]->index,
							this->vesselSegments[index]->radius);*/
				}
			}
		}
	}while(change /*|| set!=this->countVesselSegments*/);

	if(set!=this->countVesselSegments){
		this->printToEPS( "WrongRadiusCalculation.eps", 0 ,0);
		fprintf( stderr, "INFO: set=%i != this->countVesselSegments=%i\n", set, this->countVesselSegments);
		//exit( 0);
	}

}


/*void VesselGraph::updateRadiusOLD()
{
	int set = 0;
	for( int i=0; i<this->countVesselSegments; i++)
	if( !this->vesselSegments[i]->radiusStatic){
		if( this->vesselSegments[i]->vesselNodes[0]->getType() != TIP && this->vesselSegments[i]->vesselNodes[1]->getType() != TIP)
			this->vesselSegments[i]->radius = 0;
		else{
			this->vesselSegments[i]->radius = MIN_VESSEL_RADIUS;//*LATTICE_CONSTANT/6.;
			set++;
			//fprintf( stderr, "INFO: set segment %i\n",this->vesselSegments[i]->index);
		}
	}



	//fprintf( stderr, "INFO: set=%i\n",set);
	bool change;
	float alpha = 2.;//2.7;//6.;//2.7;
	do{
		//fprintf( stderr, "\rupdateRadius: %.2f\%% \b", (float)set/(float)this->countVesselSegments*100.);
		change = false;
		int i=0;
//#pragma omp parallel for
		for( i=0; i<this->countVesselNodes; i++){
			if( this->vesselNodes[i]->getType() != TIP && this->vesselNodes[i]->countNeighbors>1/*&& this->vesselNodes[i]->radius == 0* /){
				float radius = 0;
				int countNonSet = 0;
				int index = -1;
				for( int b=0; b<this->vesselNodes[i]->countNeighbors; b++){
					if( this->vesselNodes[i]->branches[b]->radiusStatic){
						// STATIC RADIUS
						countNonSet=2;
					}else
					if( this->vesselNodes[i]->branches[b]->radius != 0){
						// SET ALREADY
						alpha = (vesselNodes[i]->branches[b]->type == ARTERIE ? RADIUS_EXPONENT_ARTERIES : RADIUS_EXPONENT_VENES);
						radius += pow(this->vesselNodes[i]->branches[b]->radius, alpha) * this->vesselNodes[i]->branches[b]->countParallelVessels;
					}else{
						// NOT SET
						index = this->vesselNodes[i]->branches[b]->index;
						countNonSet++;
					}
				}
				if( countNonSet==1){
					//fprintf( stderr, "INFO: Set Radius of Segment %i\n", index);
					alpha = (vesselSegments[index]->type == ARTERIE ? RADIUS_EXPONENT_ARTERIES : RADIUS_EXPONENT_VENES);
					this->vesselSegments[index]->radius =
							MIN( MAX_VESSEL_RADIUS, pow(radius, 1./alpha));
					change = true;
					set++;
					/*fprintf( stderr, "INFO: set segment %i (node %i) -> radius=%lf\n",
							this->vesselNodes[i]->index,
							this->vesselSegments[index]->index,
							this->vesselSegments[index]->radius);* /
				}
			}
		}
	}while(change /*|| set!=this->countVesselSegments* /);

	if(set!=this->countVesselSegments){
		this->printToEPS( "test.eps", 0 ,0);
		fprintf( stderr, "INFO: set=%i != this->countVesselSegments=%i\n", set, this->countVesselSegments);
		//exit( 0);
	}
			
}*/


VesselNode::VesselNode(float x, float y, float z)
{
	this->position[0] = x;
	this->position[1] = y;
	this->position[2] = z;
	this->neighbors = 0;
	this->branches = 0;
	this->countNeighbors = 0;
	
	this->pressure = 0.;
	this->time = 0.;
	this->marker = 0.;

	this->type = 0;
	this->index = -1;
}

void VesselNode::addNeighbor(VesselNode *neighbor)
{
	this->neighbors = (VesselNode**) realloc( this->neighbors, sizeof(VesselNode*) * (this->countNeighbors+1));	
	this->neighbors[(int)this->countNeighbors] = neighbor;
	countNeighbors++;
}

void VesselNode::removeNeighbor(VesselNode *neighbor)
{
	int n=0;
	for( ; n<this->countNeighbors && this->neighbors[n] != neighbor; n++) ;

	if( n<this->countNeighbors){
		countNeighbors--;
		this->neighbors[n] = this->neighbors[(int)this->countNeighbors];
		this->neighbors = (VesselNode**) realloc( this->neighbors, sizeof(VesselNode*) * (this->countNeighbors));
		this->branches[n] = this->branches[(int)this->countNeighbors];
		this->branches = (VesselSegment**) realloc( this->branches, sizeof(VesselSegment*) * (this->countNeighbors));
	}
}

void VesselNode::addNeighbor(VesselNode *neighbor, VesselSegment *vs)
{
	this->neighbors = (VesselNode**) realloc( this->neighbors, sizeof(VesselNode*) * (this->countNeighbors+1));	
	this->neighbors[(int)this->countNeighbors] = neighbor;
	this->branches = (VesselSegment**) realloc( this->branches, sizeof(VesselSegment*) * (this->countNeighbors+1));	
	this->branches[(int)this->countNeighbors] = vs;
	countNeighbors++;
}

/*void VesselNode::addSegment(VesselSegment *vs)
{
	this->branches = (VesselSegment**) realloc( this->branches, sizeof(VesselSegment*) * (this->countNeighbors+1));	
	this->branches[(int)this->countNeighbors] = vs;
	countNeighbors++;
}*/

VesselSegment::VesselSegment(VesselNode *node1, VesselNode *node2)
{
	this->vesselNodes[0] = node1;
	this->vesselNodes[1] = node2;
	//this->vesselNodes[0]->addNeighbor(this->vesselNodes[1], this);
	//this->vesselNodes[1]->addNeighbor(this->vesselNodes[0], this);
	
	this->flow = 0;
	this->shear = 0;
	this->radiusStatic = false;
	this->radius = 4;
	this->countParallelVessels = 1.;
	this->wallThickness = 1;
	this->permeability = 0.;

	this->type = UNKNOWN;
}

VesselSegment::~VesselSegment()
{
	//this->vesselNodes[0]->removeNeighbor(this->vesselNodes[1]);
	//this->vesselNodes[1]->removeNeighbor(this->vesselNodes[0]);

	//this->radius = 4;
	//this->wallThickness = 1;
}

VesselNode::~VesselNode()
{
	free( neighbors);
	free( branches);
}

char VesselNode::getType()
{
	return type;
}

void VesselNode::setType( char type)
{
	/*if( this->type == ROOT){
		//fprintf( stderr, "Set Node %i (%i,%i) from %i -> %i\n", this->index, this->position[0], this->position[1], this->type, type);
		fprintf( stderr, "Set Node %i", this->index);
		fprintf( stderr, " (%i,%i) from", (int)this->position[0], (int)this->position[1]);
		fprintf( stderr, " %i", this->type);
		fprintf( stderr, " -> %i\n", type);
		//exit( 0);
	}*/
	
	//if( type == ROOT)
	//	fprintf( stderr, "Set Node %i (%i,%i) from %i -> %i\n", this->index, (int)this->position[0], (int)this->position[1], this->type, type);
		
	this->type = type;
}

double parkerAIF(
		double t,
		double A,
		double T,
		double sigma,
		double alpha,
		double beta,
		double tau,
		double s)
{
	return A/(sigma*sqrt(2.*PI)) * exp(-(pow(t-T,2)) / 2./pow(sigma,2.)) + alpha * exp(-beta*t) / (1.+exp(-s*(t-tau)));
}
double VesselGraph::AIF::ArterialInputFunction( //mM
		double time,  //min
		double time0, //min
		const char type)
{
	//return 1.;

	double A1 = 0.809;
	double A2 = 0.330;
	double T1 = 0.17046;
	double T2 = 0.365;
	double sigma1 = 0.0563;
	double sigma2 = 0.132;
	double alpha = 1.050;
	double beta = 0.1685;
	double s = 38.078;
	double tau = 0.483;
	double t = time - time0;
	if( t < 0)
		return 0.;
	else{
		//return (t*60. < 20. ? 6 : 0);
		//return exp(-t)*t;

		//t=t*10;

		if(true){
			return parkerAIF(t,A1,T1,sigma1, alpha,beta,tau,s) + parkerAIF(t,A2,T2,sigma2, alpha,beta,tau,s);
		}else{
			double inj=200./60.; // min
			return exp(-pow( (t-2.7*inj/4.) / (inj/4.), 2.));
		}
	}
}

VesselNode *VesselNode::getNeighborUpstream( int n)
{
	for( int i=0; i<this->countNeighbors; i++){
		if( this->pressure < this->neighbors[i]->pressure){
			if( n==0)
				return this->neighbors[i];
			n--;
		}
	}
	return 0;
}

VesselNode *VesselNode::getNeighborDownstream( int n)
{
	for( int i=0; i<this->countNeighbors; i++){
		if( this->pressure > this->neighbors[i]->pressure){
			if( n==0)
				return this->neighbors[i];
			n--;
		}
	}
	return 0;
}

#define MINMAX(  a,  b,  c) \
	( b<a?a: (b>c?c : b))


double rgbformulae(char formulae, double x) {
	double ret = 0;
	//double PI = 3.14159265;
	switch (abs(formulae)) {
	case 0:
		ret = 0.;
		break;
	case 1:
		ret = 0.5;
		break;
	case 2:
		ret = 1.;
		break;
	case 3:
		ret = x;
		break;
	case 4:
		ret = pow(x, 2.);
		break;
	case 5:
		ret = pow(x, 3.);
		break;
	case 6:
		ret = pow(x, 4.);
		break;
	case 7:
		ret = sqrt(x);
		break;
	case 8:
		ret = sqrt(sqrt(x));
		break;
	case 9:
		ret = sin(PI / 2. * x);
		break;
	case 10:
		ret = cos(PI / 2. * x);
		break;
	case 11:
		ret = fabs(x - 0.5);
		break;
	case 12:
		ret = pow(2. * x - 1., 2.);
		break;
	case 13:
		ret = sin(PI * x);
		break;
	case 14:
		ret = fabs(cos(PI * x));
		break;
	case 15:
		ret = sin(2. * PI * x);
		break;
	case 16:
		ret = cos(2. * PI * x);
		break;
	case 17:
		ret = fabs(sin(2. * PI * x));
		break;
	case 18:
		ret = fabs(cos(2. * PI * x));
		break;
	case 19:
		ret = fabs(sin(4. * PI * x));
		break;
	case 20:
		ret = fabs(cos(4. * PI * x));
		break;
	case 21:
		ret = 3. * x;
		break;
	case 22:
		ret = 3. * x - 1.;
		break;
	case 23:
		ret = 3. * x - 2.;
		break;
	case 24:
		ret = fabs(3. * x - 1.);
		break;
	case 25:
		ret = fabs(3. * x - 2.);
		break;
	case 26:
		ret = (3. * x - 1.) / 2.;
		break;
	case 27:
		ret = (3. * x - 2.) / 2.;
		break;
	case 28:
		ret = fabs((3. * x - 1.) / 2.);
		break;
	case 29:
		ret = fabs((3. * x - 2.) / 2.);
		break;
	case 30:
		ret = x / 0.32 - 0.78125;
		break;
	case 31:
		ret = 2. * x - 0.84;
		break;
		//case 32: ret = 4x;1;-2x+1.84;x/0.08-11.5; break;
	case 33:
		ret = fabs(2. * x - 0.5);
		break;
	case 34:
		ret = 2. * x;
		break;
	case 35:
		ret = 2. * x - 0.5;
		break;
	case 36:
		ret = 2. * x - 1.;
		break;
	}

	if (formulae < 0)
		return MINMAX(0., -ret, 1.);
	else
		return MINMAX(0., ret, 1.);
}

void printColorMapLegend(char *filename, char *min, char *max) {
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, 100, 0, 100);

	for (int i = 0; i < 100; i++) {
		double xPoly[5] = { (double) 0, (double) 20, (double) 20, (double) 0,
				(double) 0 };
		double yPoly[5] = { (double) i, (double) i, (double) i + 1, (double) i
				+ 1, (double) i };
		char color[512];
		sprintf(color, "%lf %lf %lf setrgbcolor", rgbformulae(22, i / 100.),
				rgbformulae(13, i / 100.), rgbformulae(-31, i / 100.));

		EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
	}
	fs << "0.0 0.0 0.0 setrgbcolor";
	fs << "/Times-Roman findfont 15 scalefont setfont 20 5 moveto (" << min
			<< ") show ";
	fs << "/Times-Roman findfont 15 scalefont setfont 20 85 moveto (" << max
			<< ") show ";

	fs.close();
}


void VesselGraph::printMarkerIntensityMap(char *filename, VesselGraph *vg,
		CONCENTRATION_T ***markerIntSpace, double voxelSizeX, double voxelSizeY,
		double voxelSizeZ, bool printAIF, float markerAIF) {
	// COARSENING
	//double avg = 5; // *60ym = 300ym = 0.3mm
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}

	//fprintf( stderr, "DOMAIN_SIZE= (%ix%ix%i)\n", DOMAIN_SIZE_X,DOMAIN_SIZE_Y,DOMAIN_SIZE_Z);
	//fprintf( stderr, "Output Size= (%ix%ix%i)\n", X,Y,Z);

	// Add Marker in Interstitial Space
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			//for (int z = 0; z < DOMAIN_SIZE_Z; z++)
			int z = 0;
			{
				//flowAvg[x/dx][y/dy][z/dz] += markerIntSpace[x][y][z];

				// volume: tissue
				float volume = LATTICE_CONSTANT * LATTICE_CONSTANT
						* LATTICE_CONSTANT;
				float volumeVessels = 0.;
				float volumeInt = 0.;

				VesselNode *vn = vg->octree->at(x, y, z);
				if (vn != 0) {
					// volume of vessels
					for (int v = 0; v < vn->countNeighbors; v++)
						volumeVessels += LATTICE_CONSTANT * PI
								* vn->branches[v]->radius
								* vn->branches[v]->radius;
					volumeVessels *= 0.5;

					// volume of interstitial space
					volumeInt = volume - volumeVessels;

					flowAvg[x / dx][y / dy][z / dz] += volumeInt / volume
							* markerIntSpace[x][y][z]
							+ volumeVessels / volume * vn->marker;
				} else {
					flowAvg[x / dx][y / dy][z / dz] += markerIntSpace[x][y][z];
				}
//				if( max<flowAvg[x/dx][y/dy][z/dz])
//					max=flowAvg[x/dx][y/dy][z/dz];

			}
		}
	}

	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++)
				flowAvg[x][y][z] /= (double) (dx * dy * dz);

	//max = ARTERIAL_MARKER_CONC;

	if (printAIF)
		flowAvg[0][0][0] = markerAIF;

	min = 0.;
	max = 1.; //(dx*dy*dz);

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512] = "";
				sprintf(color, "%lf %lf %lf setrgbcolor",
						(flowAvg[x][y][z] - min) / (max - min),
						(flowAvg[x][y][z] - min) / (max - min),
						(flowAvg[x][y][z] - min) / (max - min));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	// BINARY OUTPUT
	char binfilename[512];
	sprintf(binfilename, "%s.bin", filename);
	FILE *pFile = fopen(binfilename, "wb+");
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				double value = (flowAvg[x][y][z] - min) / (max - min);
				fwrite(&value, 1, sizeof(value), pFile);
			}
	fclose(pFile);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
}

void VesselGraph::printMarkerIntensityToBinary(char *filename)
{
	FILE *fpp = fopen(filename, "wb+");

	// Vessel Graph
	for (int i = 0; i < this->countVesselNodes; i++){
		float Cp = this->vesselNodes[i]->marker;
		fwrite( &Cp, sizeof(float), 1, fpp);
	}

	fclose(fpp);
}

void VesselGraph::readMarkerIntensityFromBinary(char *filename)
{
	FILE *fpp = fopen(filename, "rb");

	// Vessel Graph
	for (int i = 0; i < this->countVesselNodes; i++){
		float Cp;
		fread( &Cp,  sizeof(float), 1, fpp);
		this->vesselNodes[i]->marker = Cp;
	}
	fclose(fpp);
}

void VesselGraph::printMarkerIntensityToBinary(char *filename, VesselGraph *vg,
		CONCENTRATION_T ***markerIntSpace, VoronoiDiagram *t, bool printAIF = false, CONCENTRATION_T markerAIF = 0) {

	char filename_C[512];
	char filename_Ci[512];
	char filename_Cp[512];

	sprintf( filename_C,  "%s.C.binary",  filename);
	sprintf( filename_Ci, "%s.Ci.binary", filename);
	sprintf( filename_Cp, "%s.Cp.binary", filename);

	FILE *fp, *fpi;

	fp  = fopen(filename_C,  "wb+");
	fpi = fopen(filename_Ci, "wb+");


	// AVERAGE MARKER

	// Binary Header
	fwrite( &DOMAIN_SIZE_X,  sizeof(int), 1, fp);
	fwrite( &DOMAIN_SIZE_Y,  sizeof(int), 1, fp);
	fwrite( &DOMAIN_SIZE_Z,  sizeof(int), 1, fp);
	// Binary Data
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				VesselNode *vn = vg->octree->at(x, y, z);


				float C, Ci, Cp;
				Ci = markerIntSpace[x][y][z];

				if (vn != 0) {
					// volume: tissue
					float volume = LATTICE_CONSTANT * LATTICE_CONSTANT * LATTICE_CONSTANT;
					float volumeVessels = this->getVascularVolume( x,y,z);
					float volumeInt = this->getExtraVascularVolume( x,y,z);

					Cp = vn->marker;
					C  = volumeInt     / volume * Ci
					   + volumeVessels / volume * Cp;
				} else {
					Cp = 0;
					C  = markerIntSpace[x][y][z];
				}

				if (printAIF && x==0 && y==0 && z==0){
					C=Cp = markerAIF;
					//Ci = 0;
				}

				fwrite( &C,  sizeof(float), 1, fp);
				//fwrite( &Ci, sizeof(float), 1, fpi);
				//fwrite( &Cp, sizeof(float), 1, fpp);
			}
		}
	}


	fclose(fp);

	// INTERSTITIAL SPACE
	switch( this->type())
	{
		case RegularLattice:{
			// Binary Header
			fwrite( &DOMAIN_SIZE_X,  sizeof(int), 1, fpi);
			fwrite( &DOMAIN_SIZE_Y,  sizeof(int), 1, fpi);
			fwrite( &DOMAIN_SIZE_Z,  sizeof(int), 1, fpi);
			// Binary Data
			for (int x = 0; x < DOMAIN_SIZE_X; x++) {
				for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
					for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

						float Ci = markerIntSpace[x][y][z];

						if (printAIF && x==0 && y==0 && z==0){
							Ci = 0;
						}

						fwrite( &Ci, sizeof(float), 1, fpi);
						/*if( Ci!=0)
							fprintf(stderr, "Ci=%e\n", Ci);*/
					}
				}
			}
		}break;

		case IrregularLattice:{
			for( int i=0; i<t->numberOfVertices(); i++){
				float Ci = t->get(i)->conc;
				fwrite( &Ci,  sizeof(float), 1, fpi);
			}

		}break;
	}
	fclose(fpi);

	// PLASMA
	printMarkerIntensityToBinary( filename_Cp);
}



void VesselGraph::printFlowMap(char *filename, VesselGraph *vg, double voxelSizeX,
		double voxelSizeY, double voxelSizeZ) {
	// COARSENING
	//double avg = 5; // *60ym = 300ym = 0.3mm

	// domain size (in voxels)
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	// voxel size (in lattice sites)
	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	// INIT AVERAGE DATA STRUCTURE
	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}
	double ***flowAvgOut = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvgOut[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvgOut[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvgOut[x][y][z] = 0.;
			}
		}
	}

	// SUM UP VALUES PER VOXEL
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				// actual lattice site
				VesselNode *vn = vg->octree->at(x, y, z);

				if (vn != 0) {
					// out flow into neighboring voxel?
					for (int v = 0; v < vn->countNeighbors; v++) {
						//double outFlow = 0.;
						//double inFlow = 0.;
						// out flow?
						if (vn->pressure >= vn->neighbors[v]->pressure) {
							// yes!

							// neighbor in neighboring voxel?
							if ((int)(vn->neighbors[v]->position[0] / dx) != (int)(x / dx)
							 || (int)(vn->neighbors[v]->position[1] / dy) != (int)(y / dy)
							 || (int)(vn->neighbors[v]->position[2] / dz) != (int)(z / dz)) {
								// yes!
								flowAvgOut[x / dx][y / dy][z / dz] += fabs(
										vn->branches[v]->flow);
								//outFlow += fabs(vn->branches[v]->flow);
								//fprintf( stderr, "ADD\n");
							}else{
								fprintf( stderr, "REJECT\n");
							}
						}
						// in flow!
						else {
							// neighbor in neighboring voxel?
							if ((int)(vn->neighbors[v]->position[0] / dx) != (int)(x / dx)
							 || (int)(vn->neighbors[v]->position[1] / dy) != (int)(y / dy)
							 || (int)(vn->neighbors[v]->position[2] / dz) != (int)(z / dz)) {
								// yes!
								flowAvg[x / dx][y / dy][z / dz] += fabs(
										vn->branches[v]->flow);
							}
						}
						//flowAvg[x/dx][y/dy][z/dz] += MAX(inFlow,outFlow);
					}
				}
			}
		}
	}

	// NORMALIZE
	/*for (int x = 0; x < X; x++)
	 for (int y = 0; y < Y; y++)
	 for (int z = 0; z < 1; z++)
	 flowAvg[x][y][z] /= (double)(dx*dy*dz);
	 */

	// MIN & MAX
	min = MAX(flowAvg[0][0][0], flowAvgOut[0][0][0]);
	max = MAX(flowAvg[0][0][0], flowAvgOut[0][0][0]);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				flowAvg[x][y][z] = MAX(flowAvg[x][y][z], flowAvgOut[x][y][z]);
				if (min > flowAvg[x][y][z])
					min = flowAvg[x][y][z];

				if (max < flowAvg[x][y][z])
					max = flowAvg[x][y][z];
			}

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	if (min == max)
		min = 0;
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512];
				sprintf(
						color,
						"%lf %lf %lf setrgbcolor",
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min)
						rgbformulae(22, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(13, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(-31,
								(flowAvg[x][y][z] - min) / (max - min)));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	// LEGEND
	char legend_filename[512];
	char min_string[512];
	char max_string[512];
	sprintf(legend_filename, "%s.legend.eps", filename);
	sprintf(min_string, "%lf", min);
	sprintf(max_string, "%lf", max);
	printColorMapLegend(legend_filename, min_string, max_string);

	// DATA
	char data_filename[512];
	sprintf(data_filename, "%s.data.dat", filename);
	FILE * fp = fopen(data_filename, "w+");
	for (int x = 0; x < X; x++){
		for (int y = 0; y < Y; y++){
			for (int z = Z/2; z < Z/2+1; z++)
				fprintf(fp, "%i %i %e\n", x, y, flowAvg[x][y][z]);
		}
		// finish gnuplot data block
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvgOut[x][y]);
		free(flowAvgOut[x]);
	}
	free(flowAvgOut);

}


void VesselGraph::printPermeabilityMap(char *filename, VesselGraph *vg, double voxelSizeX,
		double voxelSizeY, double voxelSizeZ) {
	// COARSENING
//	double avg = 5; // *60ym = 300ym = 0.3mm

	// domain size (in voxels)
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	// voxel size (in lattice sites)
	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	// INIT AVERAGE DATA STRUCTURE
	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}

	// SUM UP VALUES PER VOXEL
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				// actual lattice site
				VesselNode *vn = vg->octree->at(x, y, z);

				if (vn != 0) {
					flowAvg[x / dx][y / dy][z / dz] +=
							vg->getVascularPermeabilitySurfaceProduct(vn);
					if( isnan(flowAvg[x / dx][y / dy][z / dz])){
						fprintf( stderr, "Node: %d %d %d: has strange Permeability: %e\n", x,y,z,vg->getVascularPermeabilitySurfaceProduct(vn));
						exit(0);
					}
				}
			}
		}
	}

	// NORMALIZE
	/*for (int x = 0; x < X; x++)
	 for (int y = 0; y < Y; y++)
	 for (int z = 0; z < 1; z++)
	 flowAvg[x][y][z] /= (double)(dx*dy*dz);
	 */

	// MIN & MAX
	min = flowAvg[0][0][0];
	max = flowAvg[0][0][0];
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				if (min > flowAvg[x][y][z])
					min = flowAvg[x][y][z];

				if (max < flowAvg[x][y][z])
					max = flowAvg[x][y][z];
			}
	if( max == min)
		max = min+1;

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512];
				sprintf(
						color,
						"%lf %lf %lf setrgbcolor",
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min)
						rgbformulae(22, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(13, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(-31,
								(flowAvg[x][y][z] - min) / (max - min)));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	// LEGEND
	char legend_filename[512];
	char min_string[512];
	char max_string[512];
	sprintf(legend_filename, "%s.legend.eps", filename);
	sprintf(min_string, "%lf", min);
	sprintf(max_string, "%lf", max);
	printColorMapLegend(legend_filename, min_string, max_string);

	// DATA
	char data_filename[512];
	sprintf(data_filename, "%s.data.dat", filename);
	FILE * fp = fopen(data_filename, "w+");
	for (int x = 0; x < X; x++){
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++)
				fprintf(fp, "%i %i %e\n", x, y, flowAvg[x][y][z]);

		// finish gnuplot data block
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
}


void VesselGraph::printWallShearStressMap(char *filename, VesselGraph *vg,
		double voxelSizeX, double voxelSizeY, double voxelSizeZ) {
	// COARSENING
//	double avg = 5; // *60ym = 300ym = 0.3mm

	// domain size (in voxels)
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	// voxel size (in lattice sites)
	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	// INIT AVERAGE DATA STRUCTURE
	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}

	// SUM UP VALUES PER VOXEL
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				VesselNode *vn = vg->octree->at( x,y,z);

				// actual lattice site
				if( vn){
					// average over all segments
					float temp = 0;
					for( int i=0; i<vn->countNeighbors; i++){
						temp += vn->branches[i]->shear;
					}
					temp /= (float)vn->countNeighbors;

					flowAvg[x / dx][y / dy][z / dz] += temp;
				}
			}
		}
	}

	// NORMALIZE
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++)
				flowAvg[x][y][z] /= (double) (dx * dy * dz);

	// MIN & MAX
	min = flowAvg[0][0][0];
	max = flowAvg[0][0][0];
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				if (min > flowAvg[x][y][z])
					min = flowAvg[x][y][z];

				if (max < flowAvg[x][y][z])
					max = flowAvg[x][y][z];
			}

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512];
				sprintf(
						color,
						"%lf %lf %lf setrgbcolor",
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min)
						rgbformulae(22, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(13, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(-31,
								(flowAvg[x][y][z] - min) / (max - min)));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	char legend_filename[512];
	char min_string[512];
	char max_string[512];
	sprintf(legend_filename, "%s.legend.eps", filename);
	sprintf(min_string, "%f", min);
	sprintf(max_string, "%f", max);
	printColorMapLegend(legend_filename, min_string, max_string);

	// DATA
	char data_filename[512];
	sprintf(data_filename, "%s.data.dat", filename);
	FILE * fp = fopen(data_filename, "w+");
	for (int x = 0; x < X; x++){
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++)
				fprintf(fp, "%i %i %e\n", x, y, flowAvg[x][y][z]);

		// finish gnuplot data block
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
}

void VesselGraph::printBloodVolumeFractionMap(char *filename, VesselGraph *vg,
		double voxelSizeX, double voxelSizeY, double voxelSizeZ) {
	// COARSENING
//	double avg = 5; // *60ym = 300ym = 0.3mm

	// domain size (in voxels)
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	// voxel size (in lattice sites)
	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	// INIT AVERAGE DATA STRUCTURE
	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}

	// SUM UP VALUES PER VOXEL
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				// actual lattice site
				flowAvg[x / dx][y / dy][z / dz] +=
						(double) vg->getVascularVolume(x, y, z)	/ (LATTICE_CONSTANT * LATTICE_CONSTANT * LATTICE_CONSTANT);
			}
		}
	}

	// NORMALIZE
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++)
				flowAvg[x][y][z] /= (double) (dx * dy * dz);

	// MIN & MAX
	min = flowAvg[0][0][0];
	max = flowAvg[0][0][0];
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				if (min > flowAvg[x][y][z])
					min = flowAvg[x][y][z];

				if (max < flowAvg[x][y][z])
					max = flowAvg[x][y][z];
			}

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512];
				sprintf(
						color,
						"%lf %lf %lf setrgbcolor",
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min)
						rgbformulae(22, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(13, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(-31,(flowAvg[x][y][z] - min) / (max - min)));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	char legend_filename[512];
	char min_string[512];
	char max_string[512];
	sprintf(legend_filename, "%s.legend.eps", filename);
	sprintf(min_string, "%f", min);
	sprintf(max_string, "%f", max);
	printColorMapLegend(legend_filename, min_string, max_string);

	// DATA
	char data_filename[512];
	sprintf(data_filename, "%s.data.dat", filename);
	FILE * fp = fopen(data_filename, "w+");
	for (int x = 0; x < X; x++){
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++)
				fprintf(fp, "%i %i %e\n", x, y, flowAvg[x][y][z]);

		// finish gnuplot data block
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
}


void VesselGraph::printInterstitialSpaceVolumeFractionMap(char *filename, VesselGraph *vg,
		double voxelSizeX, double voxelSizeY, double voxelSizeZ) {
	// COARSENING
//	double avg = 5; // *60ym = 300ym = 0.3mm

	// domain size (in voxels)
	int X = MAX( 1, DOMAIN_SIZE_X * LATTICE_CONSTANT / voxelSizeX);
	int Y = MAX( 1, DOMAIN_SIZE_Y * LATTICE_CONSTANT / voxelSizeY);
	int Z = MAX( 1, DOMAIN_SIZE_Z * LATTICE_CONSTANT / voxelSizeZ);

	// voxel size (in lattice sites)
	int dx = DOMAIN_SIZE_X / X;
	int dy = DOMAIN_SIZE_Y / Y;
	int dz = DOMAIN_SIZE_Z / Z;

	// INIT AVERAGE DATA STRUCTURE
	double ***flowAvg = (double***) malloc(sizeof(double**) * X);
	for (int x = 0; x < X; x++) {
		flowAvg[x] = (double**) malloc(sizeof(double*) * Y);
		for (int y = 0; y < Y; y++) {
			flowAvg[x][y] = (double*) malloc(sizeof(double) * Z);
			for (int z = 0; z < Z; z++) {
				flowAvg[x][y][z] = 0.;
			}
		}
	}

	// SUM UP VALUES PER VOXEL
	float min = 0, max = 0;
	for (int x = 0; x < DOMAIN_SIZE_X; x++) {
		for (int y = 0; y < DOMAIN_SIZE_Y; y++) {
			for (int z = 0; z < DOMAIN_SIZE_Z; z++) {

				// actual lattice site
				flowAvg[x / dx][y / dy][z / dz] +=
						(double) vg->getExtraVascularVolume(x, y, z)
								/ (LATTICE_CONSTANT * LATTICE_CONSTANT
										* LATTICE_CONSTANT);
			}
		}
	}

	// NORMALIZE
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++)
				flowAvg[x][y][z] /= (double) (dx * dy * dz);

	// MIN & MAX
	min = flowAvg[0][0][0];
	max = flowAvg[0][0][0];
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = 0; z < 1; z++) {
				if (min > flowAvg[x][y][z])
					min = flowAvg[x][y][z];

				if (max < flowAvg[x][y][z])
					max = flowAvg[x][y][z];
			}

	// OUTPUT
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, 0, X, 0, Y);
	for (int x = 0; x < X; x++)
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++) {
				double xPoly[5] = { (double) x, (double) x + 1, (double) x + 1,
						(double) x, (double) x };
				double yPoly[5] = { (double) y, (double) y, (double) y + 1,
						(double) y + 1, (double) y };
				char color[512];
				sprintf(
						color,
						"%lf %lf %lf setrgbcolor",
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min),
						//(flowAvg[x][y][z]-min)/(max-min)
						rgbformulae(22, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(13, (flowAvg[x][y][z] - min) / (max - min)),
						rgbformulae(-31,
								(flowAvg[x][y][z] - min) / (max - min)));
				EPS::PSfillPolygon(&fs, xPoly, yPoly, 4, 0.1, color);
			}
	fs.close();

	char legend_filename[512];
	char min_string[512];
	char max_string[512];
	sprintf(legend_filename, "%s.legend.eps", filename);
	sprintf(min_string, "%f", min);
	sprintf(max_string, "%f", max);
	printColorMapLegend(legend_filename, min_string, max_string);

	// DATA
	char data_filename[512];
	sprintf(data_filename, "%s.data.dat", filename);
	FILE * fp = fopen(data_filename, "w+");
	for (int x = 0; x < X; x++){
		for (int y = 0; y < Y; y++)
			for (int z = Z/2; z < Z/2+1; z++)
				fprintf(fp, "%i %i %e\n", x, y, flowAvg[x][y][z]);

		// finish gnuplot data block
		//fprintf(fp, "\n");
	}
	fclose(fp);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++)
			free(flowAvg[x][y]);
		free(flowAvg[x]);
	}
	free(flowAvg);
}


VesselNode *VesselGraph::getClosestVesselNode(float *pos, float &sqrDist) {

	VesselGraph *vg = this;

	if (vg->countVesselNodes == 0)
		return 0;

	float dist;
	sqrDist = pow(pos[0] - vg->vesselNodes[0]->position[0], 2)
			+ pow(pos[1] - vg->vesselNodes[0]->position[1], 2)
			+ pow(pos[2] - vg->vesselNodes[0]->position[2], 2);
	VesselNode *closestVN = vg->vesselNodes[0];
	for (int v = 1; v < vg->countVesselNodes; v++) {
		dist = pow(pos[0] - vg->vesselNodes[v]->position[0], 2)
				+ pow(pos[1] - vg->vesselNodes[v]->position[1], 2)
				+ pow(pos[2] - vg->vesselNodes[v]->position[2], 2);
		if (sqrDist > dist) {
			sqrDist = dist;
			closestVN = vg->vesselNodes[v];
		}
	}

	return closestVN;
}


void VesselGraph::setSingleVessel(int length, int height, int depth, double vessel_radius) {

	VesselGraph *vesselGraph = this;
	VesselNode *last, *next;

	// SET ROOT
	last = new VesselNode(0, height, depth);
	last->pressure = MAX_PRESSURE; // kPa
	last->setType(ROOT);
	vesselGraph->addVesselNode(last);

	for( int i=0; i<length-1; i++){
		// SET SEGMENT
		next = new VesselNode(i+1, height, depth);
		next->pressure = MAX_PRESSURE; // kPa
		next->setType(VESSEL);
		vesselGraph->addVesselNode(next);

		vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));


		last=next;
	}

	next->pressure = MIN_PRESSURE; // kPa
	next->setType(ROOT);

	for( int i=0; i<length-1; i++){
		vesselGraph->vesselSegments[i]->radius = vessel_radius;//MAX_VESSEL_RADIUS;
		vesselGraph->vesselSegments[i]->radiusStatic = true;
	}
}


void VesselGraph::setSingleRandomVessel(int length, int height, int depth) {

	VesselGraph *vesselGraph = this;
	VesselNode *last, *next;

	// SET ROOT
	last = new VesselNode(0, height, depth);
	last->pressure = MAX_PRESSURE; // kPa
	last->setType(ROOT);
	vesselGraph->addVesselNode(last);

	next=last;
	int i=0;
	while( i<length){

		switch( (int)(3.*RAND01)){
		case 0:
			// SET SEGMENT -> right
			next = new VesselNode(++i, height, depth);
			break;
		case 1:
			// SET SEGMENT -> up
			if( height<DOMAIN_SIZE_Y-1 && !this->octree->at(i, height+1, depth))
				next = new VesselNode(i, ++height, depth);
			break;

		case 2:
			// SET SEGMENT -> down
			if( height>0 && !this->octree->at(i, height-1, depth))
				next = new VesselNode(i, --height, depth);
			break;
		}
		if( last!=next){
			next->pressure = MAX_PRESSURE; // kPa
			next->setType(VESSEL);
			vesselGraph->addVesselNode(next);

			vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));


			last=next;

		}

	}

	next->pressure = MIN_PRESSURE; // kPa
	next->setType(ROOT);

	for( int i=0; i<vesselGraph->countVesselSegments; i++){
		vesselGraph->vesselSegments[i]->radius = 20;//MAX_VESSEL_RADIUS;
		vesselGraph->vesselSegments[i]->radiusStatic = true;
	}
}


void VesselGraph::setSingleVessel2(int length) {

	VesselGraph *vesselGraph = this;
	VesselNode *last, *next;

	// SET ROOT
	last = new VesselNode(0, 0, 0);
	last->pressure = MAX_PRESSURE; // kPa
	last->setType(ROOT);
	vesselGraph->addVesselNode(last);

	for( int i=0; i<length/2; i++){
		// SET SEGMENT
		next = new VesselNode(i+1, 0, 0);
		next->pressure = MAX_PRESSURE; // kPa
		next->setType(VESSEL);
		vesselGraph->addVesselNode(next);

		vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));

		last=next;
	}
	last->setType(TIP);

	last = new VesselNode(length/2+1, 0, 0);
	last->pressure = MAX_PRESSURE; // kPa
	last->setType(TIP);
	vesselGraph->addVesselNode(last);

	for( int i=length/2+1; i<length; i++){
		// SET SEGMENT
		next = new VesselNode(i+1, 0, 0);
		next->pressure = MAX_PRESSURE; // kPa
		next->setType(VESSEL);
		vesselGraph->addVesselNode(next);
		vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));

		last=next;
	}

	next->pressure = MIN_PRESSURE; // kPa
	next->setType(ROOT);
}


void VesselGraph::setSymmetricVesselGraph(int height) {

	int z=(int)(DOMAIN_SIZE_Z/2);

	VesselGraph *vesselGraph = this;
	//int height=6;
	VesselNode *last, *next;

	// SET ROOT
	last = new VesselNode(0, pow(2, height - 1), z);
	last->pressure = MAX_PRESSURE; // kPa
	last->setType(ROOT);
	vesselGraph->addVesselNode(last);

	// SET SEGMENT
	next = new VesselNode(1, pow(2, height - 1), z);
	//next->pressure = MAX_PRESSURE; // kPa
	next->setType(VESSEL);
	vesselGraph->addVesselNode(next);
	vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));
	fprintf(stderr, "set (%i,%i)\n", 0, (int) pow(2, height - 1));

	int y, y0 = pow(2, height - 1);
	/*for( y=0; y<y0-1; y++)	fprintf(stderr, " ");
	 fprintf(stderr, "o\n");*/

	for (int depth = 1; depth < height; depth++) {
		int y0 = pow(2, height - depth - 1);
		int dy = pow(2, height - depth);
		int n = pow(2, depth);

		int y0_old = pow(2, height - depth + 1 - 1);
		int dy_old = pow(2, height - depth + 1);
		int n_old = pow(2, depth - 1);

		char type = (depth == height - 1 ? TIP : VESSEL);

		//y0
		/*for( y=0; y<y0-1; y++) fprintf(stderr, " ");
		 //nodes
		 for( int i=0; i<n; i++){
		 fprintf(stderr, "o");
		 for( y=0; y<dy-1; y++) fprintf(stderr, " ");
		 }
		 fprintf(stderr, "\n");*/

		// sub-tree
		for (int i = 0; i < n_old; i++) {
			VesselNode *root = vesselGraph->octree->at(depth,y0_old + i * dy_old, z);
			fprintf(stderr, "(%i,%i)\n", depth, y0_old + i * dy_old);
			if (root == 0) {
				exit(0);
			}

			last = root;
			for (y = y0_old + i * dy_old + 1; y <= y0_old + i * dy_old + dy / 2;
					y++) {
				next = new VesselNode(depth, y, z);
				next->setType(VESSEL);
				vesselGraph->addVesselNode(next);
				fprintf(stderr, "set (%i,%i)\n", depth, y);
				vesselGraph->addVesselSegmentAndNeighbors(
						new VesselSegment(last, next));
				last = next;
			}
			next = new VesselNode(depth + 1, y - 1, z);
			next->setType(type);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(
					new VesselSegment(last, next));
			fprintf(stderr, "set (%i,%i)\n", depth + 1, y - 1);

			last = root;
			for (y = y0_old + i * dy_old - 1; y >= y0_old + i * dy_old - dy / 2;
					y--) {
				next = new VesselNode(depth, y, z);
				next->setType(VESSEL);
				vesselGraph->addVesselNode(next);
				fprintf(stderr, "set (%i,%i)\n", depth, y);
				vesselGraph->addVesselSegmentAndNeighbors(
						new VesselSegment(last, next));
				last = next;
			}
			next = new VesselNode(depth + 1, y + 1, z);
			next->setType(type);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(
					new VesselSegment(last, next));
			fprintf(stderr, "set (%i,%i)\n", depth + 1, y + 1);

		}

		//VesselNode *root = vesselGraph->octree->at( depth, pow( height-depth, 2), 0);
		// UP

	}
	////////////////////////////

	// SET ROOT
	last = new VesselNode(height * 2 + 2, pow(2, height - 1), z);
	last->pressure = MIN_PRESSURE; // kPa
	last->setType(ROOT);
	vesselGraph->addVesselNode(last);

	// SET SEGMENT
	next = new VesselNode(height * 2 + 1, pow(2, height - 1), z);
	//next->pressure = MAX_PRESSURE; // kPa
	next->setType(VESSEL);
	vesselGraph->addVesselNode(next);
	vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));
	last = next;
	fprintf(stderr, "set (%i,%i)\n", 0, (int) pow(2, height - 1));

	// SET SEGMENT
	next = new VesselNode(height * 2, pow(2, height - 1), z);
	//next->pressure = MAX_PRESSURE; // kPa
	next->setType(VESSEL);
	vesselGraph->addVesselNode(next);
	vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));
	fprintf(stderr, "set (%i,%i)\n", 0, (int) pow(2, height - 1));

	//int y, y0=pow( 2, height-1);
	/*for( y=0; y<y0-1; y++)	fprintf(stderr, " ");
	 fprintf(stderr, "o\n");*/

	for (int depth = 1; depth < height; depth++) {
		int y0 = pow(2, height - depth - 1);
		int dy = pow(2, height - depth);
		int n = pow(2, depth);

		int y0_old = pow(2, height - depth + 1 - 1);
		int dy_old = pow(2, height - depth + 1);
		int n_old = pow(2, depth - 1);

		char type = (depth == height - 1 ? TIP : VESSEL);

		//y0
		/*for( y=0; y<y0-1; y++) fprintf(stderr, " ");
		 //nodes
		 for( int i=0; i<n; i++){
		 fprintf(stderr, "o");
		 for( y=0; y<dy-1; y++) fprintf(stderr, " ");
		 }
		 fprintf(stderr, "\n");*/

		// sub-tree
		for (int i = 0; i < n_old; i++) {
			VesselNode *root = vesselGraph->octree->at(height * 2 + 1 - depth,
					y0_old + i * dy_old, z);
			fprintf(stderr, "(%i,%i)\n", height * 2 + 1 - depth,
					y0_old + i * dy_old);
			if (root == 0) {
				exit(0);
			}

			last = root;
			for (y = y0_old + i * dy_old + 1; y <= y0_old + i * dy_old + dy / 2;
					y++) {
				next = new VesselNode(height * 2 + 1 - depth, y, z);
				next->setType(VESSEL);
				vesselGraph->addVesselNode(next);
				fprintf(stderr, "set (%i,%i)\n", height * 2 + 1 - depth, y);
				vesselGraph->addVesselSegmentAndNeighbors(
						new VesselSegment(last, next));
				last = next;
			}
			next = new VesselNode(height * 2 + 1 - (depth + 1), y - 1, z);
			next->setType(type);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(
					new VesselSegment(last, next));
			fprintf(stderr, "set (%i,%i)\n", height * 2 + 1 - (depth + 1),
					y - 1);

			last = root;
			for (y = y0_old + i * dy_old - 1; y >= y0_old + i * dy_old - dy / 2;
					y--) {
				next = new VesselNode(height * 2 + 1 - depth, y, z);
				next->setType(VESSEL);
				vesselGraph->addVesselNode(next);
				fprintf(stderr, "set (%i,%i)\n", depth, y);
				vesselGraph->addVesselSegmentAndNeighbors(
						new VesselSegment(last, next));
				last = next;
			}
			next = new VesselNode(height * 2 + 1 - (depth + 1), y + 1, z);
			next->setType(type);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(
					new VesselSegment(last, next));
			fprintf(stderr, "set (%i,%i)\n", height * 2 + 1 - (depth + 1),
					y + 1);

		}

		//VesselNode *root = vesselGraph->octree->at( depth, pow( height-depth, 2), 0);
		// UP

	}

	// shift to voxel centers
	for( int i=0; i<this->countVesselNodes; i++)
		for( int d=0; d<this->dimensions; d++)
		this->vesselNodes[i]->position[d]+=0.5;

}


void VesselGraph::setInitialVesselGraph( int type) {

	VesselGraph *vesselGraph = this;

	// PLACE ROOT NODES
	float position[3];
	float dist = 0.;
	float pressure = MAX_PRESSURE; //MIN_PRESSURE;

	//if( type == 2)
	//	pressure = MIN_PRESSURE;

	//fprintf( stderr, "Add Root Node (%i, %i)\n", (int)DOMAIN_SIZE/2,(int)DOMAIN_SIZE/2);

	VesselNode *last;
	if( dimensions == 1)
		last = new VesselNode(0, 0, 0);
	else if( dimensions == 2){
		if(type==2)
			last = new VesselNode(DOMAIN_SIZE_X/2,DOMAIN_SIZE_Y*2/3,0);
		else if(type==3){
			/*last = new VesselNode(DOMAIN_SIZE_X/2,  DOMAIN_SIZE_Y*2/3,0);
			vesselGraph->addVesselNode(last);
			last->pressure = pressure; // kPa
			last->setType(ROOT);

			last = new VesselNode(DOMAIN_SIZE_X*1/4,DOMAIN_SIZE_Y*1/3,0);
			vesselGraph->addVesselNode(last);
			last->pressure = pressure; // kPa
			last->setType(ROOT);*/

			last = new VesselNode(DOMAIN_SIZE_X*3/4,DOMAIN_SIZE_Y/2,0);
		}else
			last = new VesselNode(DOMAIN_SIZE_X/2,DOMAIN_SIZE_Y/2,0);
	}
	else
		last = new VesselNode(DOMAIN_SIZE_X/2,DOMAIN_SIZE_Y/2,0);

	vesselGraph->addVesselNode(last);
	last->pressure = pressure; // kPa
	last->setType(ROOT);
	fprintf(stderr,
			"Added Root Node (%.0f, %.0f, %.0f) -> pressure = %lf kPa\n",
			last->position[0], last->position[1], last->position[2],
			last->pressure);

	for (int i = 0; i < 1000; i++) {
		int dx[3] = { 0, 0, 0 };

		// RANDOM POSITION AT DOMAIN BORDER
		switch (i % (2 * dimensions)) {
		case 0:
			// YZ-AXIS
			position[0] = 0;
			position[1] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Y - 1) / BRANCH_LENGTH * RAND01);
			position[2] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Z - 1) / BRANCH_LENGTH * RAND01);
			dx[0] = 1;
			break;
		case 1:
			// YZ-AXIS
			position[0] = DOMAIN_SIZE_X - 1;
			position[1] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Y - 1) / BRANCH_LENGTH * RAND01);
			position[2] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Z - 1) / BRANCH_LENGTH * RAND01);
			dx[0] = -1;
			break;
		case 2:
			// XZ-AXIS
			position[0] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_X - 1) / BRANCH_LENGTH * RAND01);
			position[1] = 0;
			position[2] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Z - 1) / BRANCH_LENGTH * RAND01);
			dx[1] = 1;
			break;
		case 3:
			// XZ-AXIS
			position[0] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_X - 1) / BRANCH_LENGTH * RAND01);
			position[1] = DOMAIN_SIZE_Y - 1;
			position[2] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Z - 1) / BRANCH_LENGTH * RAND01);
			dx[1] = -1;
			break;
			//#if dimensions==3
		case 4:
			// XY-AXIS
			position[0] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_X - 1) / BRANCH_LENGTH * RAND01);
			position[1] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Y - 1) / BRANCH_LENGTH * RAND01);
			position[2] = 0;
			dx[2] = 1;
			break;
		case 5:
			// XY-AXIS
			position[0] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_X - 1) / BRANCH_LENGTH * RAND01);
			position[1] = BRANCH_LENGTH
					* floor((DOMAIN_SIZE_Y - 1) / BRANCH_LENGTH * RAND01);
			position[2] = DOMAIN_SIZE_Z - 1;
			dx[2] = -1;
			break;
			//#endif

		}

		if(type==0 || type>=2)
			pressure = MIN_PRESSURE;
		//if(type==2)
		//	pressure = MAX_PRESSURE;

		VesselNode *closest = vesselGraph->getClosestVesselNode( position, dist);
		if(type==1){
			if (closest != 0)
				pressure = (closest->pressure == MIN_PRESSURE ? MAX_PRESSURE : MIN_PRESSURE);
			else
				pressure = (pressure == MIN_PRESSURE ? MAX_PRESSURE : MIN_PRESSURE);
		}

		if (vesselGraph->countVesselNodes
				== 0 || dist >= ROOT_DISTANCE*ROOT_DISTANCE) {
			if (closest != 0)
				fprintf(stderr, "Dist: %lf, Type: %i, Pressure: %lf\n", dist,
						closest->getType(), closest->pressure);
			VesselNode *last = new VesselNode(position[0], position[1],
					position[2]);
			vesselGraph->addVesselNode(last);
			last->pressure = pressure; //MIN_PRESSURE; //MAX_PRESSURE;//pressure; // kPa
			//last->pressure = pressure; // kPa
			last->setType(ROOT);
			fprintf(stderr,
					"Added Root Node (%i, %i, %i) -> pressure = %lf kPa\n",
					(int) position[0], (int) position[1], (int) position[2],
					last->pressure);

			VesselNode *next;
			for (int i = 0; i < BRANCH_LENGTH; i++) {
				next = new VesselNode(last->position[0] + dx[0],
						last->position[1] + dx[1], last->position[2] + dx[2]);
				next->setType(VESSEL);
				vesselGraph->addVesselNode(next);
				vesselGraph->addVesselSegmentAndNeighbors(
						new VesselSegment(last, next));
				last = next;
			}
			next->setType(TIP);

		}
	}

	if(type==3){
		VesselNode *next;
		{
			last = new VesselNode(DOMAIN_SIZE_X/2,  DOMAIN_SIZE_Y*3/4,0);
			vesselGraph->addVesselNode(last);
			last->pressure = MAX_PRESSURE; // kPa
			last->setType(ROOT);

			next = new VesselNode(last->position[0],
								  last->position[1]+1,
								  last->position[2]);
			next->setType(TIP);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));
		}
		{
			last = new VesselNode(DOMAIN_SIZE_X*1/3,DOMAIN_SIZE_Y*1/5,0);
			vesselGraph->addVesselNode(last);
			last->pressure = MAX_PRESSURE; // kPa
			last->setType(ROOT);
			next = new VesselNode(last->position[0],
								  last->position[1]+1,
								  last->position[2]);
			next->setType(TIP);
			vesselGraph->addVesselNode(next);
			vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment(last, next));
		}
	}


	// GROW VESSEL GRAPH STOCHASTICALLY
	bool changed;
	do {
		changed = false;
		int countPreviousVesselNodes = vesselGraph->countVesselNodes;
		for (int v = 0; v < countPreviousVesselNodes; v++)
		//if (RAND01 < 0.5)
				{

			// BRANCHING
			if ((vesselGraph->vesselNodes[v]->getType() == TIP /*|| vesselGraph->vesselNodes[v]->getType() == ROOT*/)
					&& vesselGraph->vesselNodes[v]->countNeighbors <= 1) {
				//fprintf( stderr, "BRANCH?\n");
				int x = vesselGraph->vesselNodes[v]->position[0];
				int y = vesselGraph->vesselNodes[v]->position[1];
				int z = vesselGraph->vesselNodes[v]->position[2];
				if (x >= BRANCH_LENGTH
						&& vesselGraph->octree->at(x - BRANCH_LENGTH, y, z) == 0
						&& x < DOMAIN_SIZE_X - BRANCH_LENGTH
						&& vesselGraph->octree->at(x + BRANCH_LENGTH, y, z) == 0) {
					//fprintf( stderr, "BRANCH\n");
					changed = true;
					if (vesselGraph->vesselNodes[v]->getType() != ROOT)
						vesselGraph->vesselNodes[v]->setType(VESSEL);
					VesselNode *last, *next;

					int dx[3] = { 1, 0, 0 };
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);

					dx[0] = -1;
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);
				} else if (y >= BRANCH_LENGTH
						&& vesselGraph->octree->at(x, y - BRANCH_LENGTH, z) == 0
						&& y < DOMAIN_SIZE_Y - BRANCH_LENGTH
						&& vesselGraph->octree->at(x, y + BRANCH_LENGTH, z) == 0) {
					//fprintf( stderr, "BRANCH\n");
					changed = true;
					if (vesselGraph->vesselNodes[v]->getType() != ROOT)
						vesselGraph->vesselNodes[v]->setType(VESSEL);
					VesselNode *last, *next;

					int dx[3] = { 0, 1, 0 };
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);

					dx[1] = -1;
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);
				} else if (z >= BRANCH_LENGTH
						&& vesselGraph->octree->at(x, y, z - BRANCH_LENGTH) == 0
						&& z < DOMAIN_SIZE_Z - BRANCH_LENGTH
						&& vesselGraph->octree->at(x, y, z + BRANCH_LENGTH) == 0) {
					//fprintf( stderr, "BRANCH\n");
					changed = true;
					if (vesselGraph->vesselNodes[v]->getType() != ROOT)
						vesselGraph->vesselNodes[v]->setType(VESSEL);
					VesselNode *last, *next;

					int dx[3] = { 0, 0, 1 };
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);

					dx[2] = -1;
					last = vesselGraph->vesselNodes[v];
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);
				}
			}

			// SINGLE SPROUT
			if (   (vesselGraph->vesselNodes[v]->getType() == ROOT && vesselGraph->vesselNodes[v]->countNeighbors == 0)
				|| (vesselGraph->vesselNodes[v]->getType() == TIP  && vesselGraph->vesselNodes[v]->countNeighbors < 2)) {
				int countCandidates = 0;
				int candidateX[8];
				int candidateY[8];
				int candidateZ[8];
				int x = vesselGraph->vesselNodes[v]->position[0];
				int y = vesselGraph->vesselNodes[v]->position[1];
				int z = vesselGraph->vesselNodes[v]->position[2];
				//float maxVEGF = vegf[x][y];

				for (int xx = x - BRANCH_LENGTH; xx <= x + BRANCH_LENGTH;
						xx += 2 * BRANCH_LENGTH)
					if (xx >= 0 && xx < DOMAIN_SIZE_X
					//&& (y<BRANCH_LENGTH || vesselGraph->octree->at(xx,y-BRANCH_LENGTH,0) == 0)
							&& (vesselGraph->octree->at(xx, y, z) == 0)
							//&& (y>=DOMAIN_SIZE-BRANCH_LENGTH || vesselGraph->octree->at(xx,y+BRANCH_LENGTH,0) == 0)
							) {
						//maxVEGF = vegf[xx][y];
						candidateX[countCandidates] = xx;
						candidateY[countCandidates] = y;
						candidateZ[countCandidates] = z;
						countCandidates++;
					}

				for (int yy = y - BRANCH_LENGTH; yy <= y + BRANCH_LENGTH;
						yy += 2 * BRANCH_LENGTH)
					if (yy >= 0 && yy < DOMAIN_SIZE_Y
					//&& (x==0 || vesselGraph->octree->at(x-BRANCH_LENGTH,yy,0) == 0)
							&& vesselGraph->octree->at(x, yy, z) == 0
							//&& (x>=DOMAIN_SIZE-BRANCH_LENGTH || vesselGraph->octree->at(x+BRANCH_LENGTH,yy,0) == 0)
							) {
						//maxVEGF = vegf[x][yy];
						candidateX[countCandidates] = x;
						candidateY[countCandidates] = yy;
						candidateZ[countCandidates] = z;
						countCandidates++;
					}

				for (int zz = z - BRANCH_LENGTH; zz <= z + BRANCH_LENGTH;
						zz += 2 * BRANCH_LENGTH)
					if (zz >= 0 && zz < DOMAIN_SIZE_Z
					//&& (x==0 || vesselGraph->octree->at(x-BRANCH_LENGTH,yy,0) == 0)
							&& vesselGraph->octree->at(x, y, zz) == 0
							//&& (x>=DOMAIN_SIZE-BRANCH_LENGTH || vesselGraph->octree->at(x+BRANCH_LENGTH,yy,0) == 0)
							) {
						//maxVEGF = vegf[x][yy];
						candidateX[countCandidates] = x;
						candidateY[countCandidates] = y;
						candidateZ[countCandidates] = zz;
						countCandidates++;
					}
				//fprintf( stderr, "Sprout! maxVEGF=%e\n", maxVEGF);

				VesselNode *last, *next;
				if (countCandidates > 0) {
					last = vesselGraph->vesselNodes[v];
					if (vesselGraph->vesselNodes[v]->getType() != ROOT)
						vesselGraph->vesselNodes[v]->setType(VESSEL);
					int which = RAND01 * countCandidates;
					int dx[3] = { (candidateX[which]
							- vesselGraph->vesselNodes[v]->position[0])
							/ BRANCH_LENGTH, (candidateY[which]
							- vesselGraph->vesselNodes[v]->position[1])
							/ BRANCH_LENGTH, (candidateZ[which]
							- vesselGraph->vesselNodes[v]->position[2])
							/ BRANCH_LENGTH };
					for (int i = 0; i < BRANCH_LENGTH; i++) {
						next = new VesselNode(
								vesselGraph->vesselNodes[v]->position[0]
										+ (i + 1) * dx[0],
								vesselGraph->vesselNodes[v]->position[1]
										+ (i + 1) * dx[1],
								vesselGraph->vesselNodes[v]->position[2]
										+ (i + 1) * dx[2]);
						next->setType(VESSEL);
						vesselGraph->addVesselNode(next);
						vesselGraph->addVesselSegmentAndNeighbors(
								new VesselSegment(last, next));
						last = next;
					}
					next->setType(TIP);
					changed = true;
				}

			}
		}
	} while (changed);

}



void VesselGraph::setInterTipConnections(Tumor *tumor, int &countInterTipConnections) {

	VesselGraph *vesselGraph = this;

	//int countPreviousVesselNodes = vesselGraph->countVesselNodes;
	countInterTipConnections = 0;
	for (int v = 0; v < vesselGraph->countVesselNodes; v++)
		if (vesselGraph->vesselNodes[v]->getType() == TIP //|| vesselGraph->vesselNodes[v]->getType() == ROOT
		) {
			float x = vesselGraph->vesselNodes[v]->position[0];
			float y = vesselGraph->vesselNodes[v]->position[1];
			float z = vesselGraph->vesselNodes[v]->position[2];

			VesselNode *vn;

			if (x >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x - BRANCH_LENGTH, y, z);
				if (vn != 0 && (vn->getType() == TIP //|| vn->getType() == ROOT
						)) {
					countInterTipConnections++;
					VesselSegment *newVN = new VesselSegment(vesselGraph->vesselNodes[v], vn);
					if( tumor)
					if( tumor->isTumor( vesselGraph->vesselNodes[v]->position) || tumor->isTumor(vn->position))
						newVN->countParallelVessels = tumor->getParallelVessels();
					vesselGraph->addVesselSegmentAndNeighbors( newVN);
				}
			}

			if (y >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x, y - BRANCH_LENGTH, z);
				if (vn != 0 && (vn->getType() == TIP //|| vn->getType() == ROOT
						)) {
					//fprintf( stderr, "Add Inter Tip Connection\n");
					countInterTipConnections++;
					VesselSegment *newVN = new VesselSegment(vesselGraph->vesselNodes[v], vn);
					if( tumor)
					if( tumor->isTumor( vesselGraph->vesselNodes[v]->position) || tumor->isTumor(vn->position))
						newVN->countParallelVessels = tumor->getParallelVessels();
					vesselGraph->addVesselSegmentAndNeighbors( newVN);
				}
			}

			if (z >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x, y, z - BRANCH_LENGTH);
				if (vn != 0 && (vn->getType() == TIP //|| vn->getType() == ROOT
						)) {
					//fprintf( stderr, "Add Inter Tip Connection\n");
					countInterTipConnections++;
					VesselSegment *newVN = new VesselSegment(vesselGraph->vesselNodes[v], vn);
					if( tumor)
					if( tumor->isTumor( vesselGraph->vesselNodes[v]->position) || tumor->isTumor(vn->position))
						newVN->countParallelVessels = tumor->getParallelVessels();
					vesselGraph->addVesselSegmentAndNeighbors( newVN);
				}
			}

		}
}


void VesselGraph::setInterTipConnectionsAll( int &countInterTipConnections) {

	VesselGraph *vesselGraph = this;

	//int countPreviousVesselNodes = vesselGraph->countVesselNodes;
	countInterTipConnections = 0;
	for (int v = 0; v < vesselGraph->countVesselNodes; v++)
		if (vesselGraph->vesselNodes[v]->getType() == TIP)
		{
			float x = vesselGraph->vesselNodes[v]->position[0];
			float y = vesselGraph->vesselNodes[v]->position[1];
			float z = vesselGraph->vesselNodes[v]->position[2];

			VesselNode *vn;

			if (x >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x - BRANCH_LENGTH, y, z);
				if (vn != 0 && vn->getType() != ROOT) {
					//fprintf( stderr, "Add Inter Tip Connection\n");
					countInterTipConnections++;
					vesselGraph->addVesselSegmentAndNeighbors(
							new VesselSegment(vesselGraph->vesselNodes[v], vn));
				}
			}

			if (y >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x, y - BRANCH_LENGTH, z);
				if (vn != 0 && vn->getType() != ROOT) {
					//fprintf( stderr, "Add Inter Tip Connection\n");
					countInterTipConnections++;
					vesselGraph->addVesselSegmentAndNeighbors(
							new VesselSegment(vesselGraph->vesselNodes[v], vn));
				}
			}

			if (z >= BRANCH_LENGTH) {
				vn = vesselGraph->octree->at(x, y, z - BRANCH_LENGTH);
				if (vn != 0 && vn->getType() != ROOT) {
					//fprintf( stderr, "Add Inter Tip Connection\n");
					countInterTipConnections++;
					vesselGraph->addVesselSegmentAndNeighbors(
							new VesselSegment(vesselGraph->vesselNodes[v], vn));
				}
			}

		}
}


void VesselGraph::removeInterTipConnections( int &countInterTipConnections) {

	VesselGraph *vesselGraph = this;
	//int countPreviousVesselNodes = vesselGraph->countVesselNodes;
	countInterTipConnections = 0;
	for (int v = vesselGraph->countVesselSegments - 1; v >= 0; v--)
		if (vesselGraph->vesselSegments[v]->vesselNodes[0]->getType() == TIP
				&& vesselGraph->vesselSegments[v]->vesselNodes[1]->getType() == TIP) {
			vesselGraph->removeVesselSegmentAndNeighbors(
					vesselGraph->vesselSegments[v]);
			delete vesselGraph->vesselSegments[v];
			countInterTipConnections++;
		}
}



void VesselGraph::remodel(Tumor *tumor) {

	VesselGraph *vesselGraph = this;

	float minShear = 0, maxShear = 0;
	int countInterTipConnections = 0;
	vesselGraph->updateSegmentTypes();
	vesselGraph->setInterTipConnections( tumor, countInterTipConnections);
	//fprintf(stderr, "countInterTipConnections=%i\n", countInterTipConnections);
	//vesselGraph->printToPovray((char*) "before.pov", 0);
	vesselGraph->updateRadius();
	//vesselGraph->updatePressure( 0,0,0);
	vesselGraph->updatePressureNEW2();
	vesselGraph->updateFlow();
	vesselGraph->updateShear();
	/*if(it==0){
	 vesselGraph->printToPovray((char*) "before.pov", 0);
	 FILE *fp = fopen("vesselBefore.dat", "w+");
	 for (int i = 0; i < vesselGraph->countVesselSegments; i++)
	 for (int ii = 0; ii < 2; ii++)
	 fprintf(
	 fp,
	 "%i %i %e %e %e %i %lf %lf %lf\n",
	 (int) vesselGraph->vesselSegments[i]->vesselNodes[ii]->position[0],
	 (int) vesselGraph->vesselSegments[i]->vesselNodes[ii]->position[1],
	 vesselGraph->vesselSegments[i]->flow,
	 vesselGraph->vesselSegments[i]->shear,
	 vesselGraph->vesselSegments[i]->vesselNodes[ii]->pressure,
	 vesselGraph->vesselSegments[i]->vesselNodes[ii]->getType(),
	 vesselGraph->vesselSegments[i]->radius,
	 (vesselGraph->vesselSegments[i]->shear - minShear)
	 / (maxShear - minShear), (maxShear
	 - vesselGraph->vesselSegments[i]->shear)
	 / (maxShear - minShear));
	 fclose(fp);

	 }*/

	/*char filename[512];
	 sprintf( filename, "remodelling%i.eps", it);
	 vesselGraph->printToEPS( filename, 0, NULL);
	 */
	/*sprintf( filename, "remodelling%i.pov", it);
	 vesselGraph->printToPovray( filename, 0);*/

	// STATISTICS
	vesselGraph->writeStatisticsToFile("statisticsPressure.dat", 0);
	vesselGraph->writeStatisticsToFile("statisticsFlow.dat", 1);
	vesselGraph->writeStatisticsToFile("statisticsVelocity.dat", 2);
	vesselGraph->writeStatisticsToFile("statisticsShearStress.dat", 3);

	// END STATISTICS

	vesselGraph->removeInterTipConnections( countInterTipConnections);
	fprintf(stderr, "\ncountInterTipConnections=%i removed!\n",
			countInterTipConnections);

	// REMODELING
	// [min,max] shear stress

	minShear = maxShear = vesselGraph->vesselSegments[0]->shear;
	for (int i = 1; i < vesselGraph->countVesselSegments; i++)
	//if (vesselGraph->vesselSegments[i]->vesselNodes[0]->getType() == TIP
	//		|| vesselGraph->vesselSegments[i]->vesselNodes[1]->getType() == TIP)
			{
		if (minShear > vesselGraph->vesselSegments[i]->shear)
			minShear = vesselGraph->vesselSegments[i]->shear;
		if (maxShear < vesselGraph->vesselSegments[i]->shear)
			maxShear = vesselGraph->vesselSegments[i]->shear;
	}

	// DEGENERATION
	float randShear; // = RAND01 * (maxShear - minShear) + minShear;
	for (int t = 0; t < 2; t++) {
		for (int i = vesselGraph->countVesselSegments - 1; i >= 0; i--) {
			randShear = RAND01 * (maxShear - minShear) + minShear;
			if( tumor && tumor->isNecrotic( vesselGraph->vesselSegments[i]->vesselNodes[0]->position))
				randShear *= 10.;
			else
			if( tumor && tumor->isTumor( vesselGraph->vesselSegments[i]->vesselNodes[0]->position))
				randShear /= 1.;
			//if(1./4. > RAND01)
			if ((vesselGraph->vesselSegments[i]->vesselNodes[0]->getType()
					== TIP
					&& vesselGraph->vesselSegments[i]->vesselNodes[0]->countNeighbors
							== 1)
					|| (vesselGraph->vesselSegments[i]->vesselNodes[1]->getType()
							== TIP
							&& vesselGraph->vesselSegments[i]->vesselNodes[1]->countNeighbors
									== 1)) {

				if (vesselGraph->vesselSegments[i]->index != i) {
					fprintf(stderr, "WORNG INDEX!\n");
					exit(0);
				}
				if (vesselGraph->vesselSegments[i]->shear < randShear) {
					if (vesselGraph->vesselSegments[i]->vesselNodes[0]->getType()
							!= ROOT
							&& vesselGraph->vesselSegments[i]->vesselNodes[0]->countNeighbors
									== 2)
						vesselGraph->vesselSegments[i]->vesselNodes[0]->setType(
								TIP);
					if (vesselGraph->vesselSegments[i]->vesselNodes[1]->getType()
							!= ROOT
							&& vesselGraph->vesselSegments[i]->vesselNodes[1]->countNeighbors
									== 2)
						vesselGraph->vesselSegments[i]->vesselNodes[1]->setType(
								TIP);
					VesselNode
							*rvn0 =
									vesselGraph->vesselSegments[i]->vesselNodes[0],
							*rvn1 =
									vesselGraph->vesselSegments[i]->vesselNodes[1];
					vesselGraph->removeVesselSegmentAndNeighbors(
							vesselGraph->vesselSegments[i]);
					if (rvn0->countNeighbors == 0)
						vesselGraph->removeVesselNode(rvn0);
					if (rvn1->countNeighbors == 0)
						vesselGraph->removeVesselNode(rvn1);
				}
			}
		}
	}

	// SPROUTING
	for (int t = 0; t < 10; t++) {
		int countOldVesselNodes = vesselGraph->countVesselNodes;
		for (int i = 0; i < countOldVesselNodes; i++) {
			randShear = RAND01 * (maxShear - minShear) + minShear;

			if (vesselGraph->vesselNodes[i]->countNeighbors
					<= 2 && vesselGraph->vesselNodes[i]->getType()!=ROOT)
					//if(vesselGraph->vesselSegments[i]->shear > randShear)
					{float sproutingProb = 1/2.; //1./2.;//1./5.

					if( tumor && tumor->isTumor( vesselGraph->vesselNodes[i]->position))
						sproutingProb = 1;
					if( tumor && tumor->isNecrotic( vesselGraph->vesselNodes[i]->position))
						sproutingProb = 0.1;

			float dx[6][3] = { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}};
			//float dx[6][3] = {0,0,0};
			int count = 0;
			if( (int)vesselGraph->vesselNodes[i]->position[0] < DOMAIN_SIZE_X-1 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0]+1,
							(int)vesselGraph->vesselNodes[i]->position[1],
							(int)vesselGraph->vesselNodes[i]->position[2]) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=1; dx[count][1]=0; dx[count][2]=0;
				count++;
				//dx[0]=1; dx[1]=0; dx[2]=0;
			}randShear = RAND01 * (maxShear - minShear) + minShear;
			if( (int)vesselGraph->vesselNodes[i]->position[0] > 0 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0]-1,
							(int)vesselGraph->vesselNodes[i]->position[1],
							(int)vesselGraph->vesselNodes[i]->position[2]) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=-1; dx[count][1]=0; dx[count][2]=0;
				count++;
				//dx[0]=-1; dx[1]=0; dx[2]=0;
			}randShear = RAND01 * (maxShear - minShear) + minShear;
			if( (int)vesselGraph->vesselNodes[i]->position[1] < DOMAIN_SIZE_Y-1 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0],
							(int)vesselGraph->vesselNodes[i]->position[1]+1,
							(int)vesselGraph->vesselNodes[i]->position[2]) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=0; dx[count][1]=1; dx[count][2]=0;
				count++;
				//dx[0]=0; dx[1]=1; dx[2]=0;
			}randShear = RAND01 * (maxShear - minShear) + minShear;
			if( (int)vesselGraph->vesselNodes[i]->position[1] > 0 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0],
							(int)vesselGraph->vesselNodes[i]->position[1]-1,
							(int)vesselGraph->vesselNodes[i]->position[2]) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=0; dx[count][1]=-1; dx[count][2]=0;
				count++;
				//dx[0]=0; dx[1]=-1; dx[2]=0;
			}randShear = RAND01 * (maxShear - minShear) + minShear;
			if( (int)vesselGraph->vesselNodes[i]->position[2] < DOMAIN_SIZE_Z-1 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0],
							(int)vesselGraph->vesselNodes[i]->position[1],
							(int)vesselGraph->vesselNodes[i]->position[2]+1) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=0; dx[count][1]=0; dx[count][2]=1;
				count++;
				//dx[0]=0; dx[1]=0; dx[2]=1;
			}
			if( (int)vesselGraph->vesselNodes[i]->position[2] > 0 &&
					vesselGraph->octree->at(
							(int)vesselGraph->vesselNodes[i]->position[0],
							(int)vesselGraph->vesselNodes[i]->position[1],
							(int)vesselGraph->vesselNodes[i]->position[2]-1) == 0 && RAND01 < sproutingProb
			) {
				//fprintf(stderr, "CANDIDATE\n");
				dx[count][0]=0; dx[count][1]=0; dx[count][2]=-1;
				count++;
				//dx[0]=0; dx[1]=0; dx[2]=-1;
			}

			int which = (int) (RAND01 * (float)count);
			if( dx[which][0]!=0 || dx[which][1]!=0 || dx[which][2]!=0) {
				//fprintf(stderr, "SPROUT!\n");
				VesselNode *next = new VesselNode(vesselGraph->vesselNodes[i]->position[0]+dx[which][0], vesselGraph->vesselNodes[i]->position[1]+dx[which][1], vesselGraph->vesselNodes[i]->position[2]+dx[which][2]);
				if( vesselGraph->vesselNodes[i]->getType()!=ROOT )
				vesselGraph->vesselNodes[i]->setType( VESSEL);
				next->setType(TIP);
				vesselGraph->addVesselNode(next);
				vesselGraph->addVesselSegmentAndNeighbors(new VesselSegment( vesselGraph->vesselNodes[i], next));

			}
		}
	}
}
/*for (int i = 0; i < vesselGraph->countVesselNodes; i++){
 if( vesselGraph->vesselNodes[i]->getType()==TIP && vesselGraph->vesselNodes[i]->countNeighbors!=1){
 fprintf(stderr, "TIP ERROR!\n");
 exit(0);
 }
 if( vesselGraph->vesselNodes[i]->countNeighbors>3){
 fprintf(stderr, "Toomany neighbors!\n");
 exit(0);
 }
 }*/

}

inline float VesselGraph::distanceSqr( VesselNode *a, VesselNode *b)
{
	float distSqr = 0;
	for( int d=0; d<dimensions; d++)
		distSqr += pow( a->position[d]-b->position[d], 2);
	return distSqr;
}

float VesselGraph::distance( VesselNode *a, VesselNode *b)
{
	return sqrt( distanceSqr(a,b));
}




void VesselGraph::simplifyNetwork( float minSegmentLength, float maxSegmentLength)
{
	//if( maxSegmentLength==0)
	//	maxSegmentLength=FLT_MAX;//minSegmentLength;

	// FUSE SEGMENTS
	for( int i=0; i<this->countVesselNodes; i++){
		if( this->vesselNodes[i]->countNeighbors == 2){

			VesselNode *vn0 = vesselNodes[i]->neighbors[0], *vn1 = vesselNodes[i]->neighbors[1];
			if( vn0 == vn1){
				fprintf( stderr, "ERROR: vertex (%i) has twice the same neighbor vertex (%i)\n", vesselNodes[i]->index, vn1->index);
			}

			// Calculate segment length
			float dist0 = distance( vesselNodes[i], vn0);
			float dist1 = distance( vesselNodes[i], vn1);

			if( true && (dist0 < minSegmentLength || dist1 < minSegmentLength) && dist0+dist1 < maxSegmentLength){
				float radius = (vesselNodes[i]->branches[0]->radius + vesselNodes[i]->branches[1]->radius)/2.;
				bool radiusStatic = vesselNodes[i]->branches[0]->radiusStatic || vesselNodes[i]->branches[1]->radiusStatic;

				// Remove vessel node and connected segments
				this->removeVesselSegmentAndNeighbors( vesselNodes[i]->branches[1]);
				this->removeVesselSegmentAndNeighbors( vesselNodes[i]->branches[0]);
				this->removeVesselNode( vesselNodes[i]);

				// Add new Segment
				if( !vn0->isNeighbor(vn1)){
					this->addVesselSegmentAndNeighbors( new VesselSegment(vn0, vn1), radius);
					this->vesselSegments[ countVesselSegments-1]->radiusStatic = radiusStatic;
				}

//				//i=-1;
				i--;
			}
		}
	}

	// FUSE NODES
	for( int i=0; i<this->countVesselSegments; i++){

		VesselNode *vn0 = vesselSegments[i]->vesselNodes[0], *vn1 = vesselSegments[i]->vesselNodes[1];
		float length = distance( vn0, vn1);
		if( length == 0){
			fprintf( stderr, "ERROR Segment %i -> (%i) -- (%i)\n", i, vn0->index, vn1->index);
		}

		if( true && length < minSegmentLength){
			//fprintf( stderr, "Remove Segment %i\n", i);

			// Remove segment
			this->removeVesselSegmentAndNeighbors( vesselSegments[i]);

			// Copy neighbors of vn1
			int tempCount=vn1->countNeighbors;
			VesselNode *temp[tempCount];
			float tempRadius[tempCount];
			bool tempRadiusStatic[tempCount];
			for( int n=0; n<vn1->countNeighbors; n++){
				temp[n]             = vn1->neighbors[n];
				tempRadius[n]       = vn1->branches[n]->radius;
				tempRadiusStatic[n] = vn1->branches[n]->radiusStatic;
			}
			// Remove segments of vn1
			for( ; 0<vn1->countNeighbors;){
				this->removeVesselSegmentAndNeighbors( vn1->branches[0]);
			}

			// Move vn0
			move( vn0, (vn0->position[0] + vn1->position[0])/2., (vn0->position[1] + vn1->position[1])/2., (vn0->position[2] + vn1->position[2])/2.);

			// Remove vn1
			this->removeVesselNode( vn1);


			// Add neighbors of vn1 to vn0
			for( int n=0; n<tempCount; n++)
				if( !vn0->isNeighbor(temp[n])){
					this->addVesselSegmentAndNeighbors( new VesselSegment(vn0, temp[n]), tempRadius[n]);
					this->vesselSegments[ countVesselSegments-1]->radiusStatic = tempRadiusStatic[n];
				}

			//i--;
			i=-1;
		}
	}
}

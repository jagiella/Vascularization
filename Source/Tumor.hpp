/*
 * Tumor.hpp
 *
 *  Created on: 17.01.2012
 *      Author: jagiella
 */

#ifndef TUMOR_HPP_
#define TUMOR_HPP_


class Tumor{
private:
	char topology;
	float center[3];
	float radius;
	float necroticCoreRadius;
	float permeability;
	float countParallelVessels;


public:
	// Constants
	static const char NONE;
	static const char SPHERICAL;

	// Constructor / Destructor
	Tumor();
	~Tumor();

	// Methods
	void setTopology(char);
	void setRadius(float);
	void setRadiusNecroticCore(float);
	void setCenter(float*);
	void setX(float value);
	void setY(float value);
	void setZ(float value);
	void setPermeability(float);
	void setParallelVessels(float);

	bool isTumor( float*);
	bool isTumor( float x, float y, float z);
	bool isNecrotic( float* pos);

	float getCenter( int dim);
	float getRadius();
	float getNecroticCoreRadius();
	float getPermeability();
	float getParallelVessels();

	void printToEPS(const char *filename);
};

#endif /* TUMOR_HPP_ */

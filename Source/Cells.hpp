#ifndef TUMOR_HPP_
#define TUMOR_HPP_


class Cells{
public:
	enum Types{ RegularLattice, IrregularLattice};
private:
	Types _type;
	int _x, _y, _z;
	float _latticeConstant;

public:
	Cells( Types type, int x, int y=1, int z=1, float latticeConstant=1.) : _type(type),_x(x),_y(y),_z(z),_latticeConstant(latticeConstant) {};
	~Cells();
};

#endif //TUMOR_HPP_

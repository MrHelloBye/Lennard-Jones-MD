#ifndef Particle_H
#define Particle_H

#include <array>
#include "xtensor/xtensor.hpp"

constexpr unsigned int dim = 6;

class Particle
{
private:
	unsigned int particle_ID;

public:
	Particle();
	Particle(unsigned int ID, double m, xt::xtensor<double, 1>);
	
	double mass;
	xt::xtensor<double, 1> pos;
    
	void setID(unsigned int n);
	unsigned int getID() const;

	void setMass(const double mass);
	double getMass() const;

	void setPos(const xt::xtensor<double, 1>& new_Pos);
	xt::xtensor<double, 1> getPos() const;
	void printPos() const;
	
};

#endif

//Particle.cpp
#include <array>
#include <iostream>
#include <cstdio>
#include <string>
#include "xtensor/xtensor.hpp"

#include "Particle.h"


Particle::Particle()
{
	particle_ID = 0;
	mass = 1;
	pos = xt::zeros<double>({dim});
}

Particle::Particle(unsigned int ID, double m,
    xt::xtensor<double, 1> init_phase_pos)
{    
	particle_ID = ID;
	mass=m;
	pos = init_phase_pos; // Initial location in phase space
}

void Particle::setID(unsigned int ID)
{
	particle_ID = ID;
}

unsigned int Particle::getID() const
{
	return particle_ID;
}

void Particle::setMass(double new_Mass)
{
	mass = new_Mass;
}

double Particle::getMass() const
{
	return mass;
}

//Pass array as @param array<double>& vector_name
void Particle::setPos(const xt::xtensor<double, 1> &new_Pos)
{
	pos = new_Pos;
}

xt::xtensor<double, 1> Particle::getPos() const
{
	return pos;
}

void Particle::printPos() const
{
	for (unsigned int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++)
	{
		std::cout << this->getPos().at(i) << '\n';
	}
}
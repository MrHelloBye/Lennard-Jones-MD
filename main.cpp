#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <typeinfo>
#include <vector>
#include <string>
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xnorm.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xindex_view.hpp"

#include "Particle.h"
#include "vose.h"

using namespace std;

xt::xtensor<double, 0> temperature(vector<Particle>&);
xt::xtensor<double, 0> potential(xt::xtensor<double,2>, vector<Particle>&);
void scale_mom(vector<Particle>&, double);
xt::xtensor<double,3> calc_sep_vecs(vector<Particle>&);
xt::xtensor<double,2> calc_dist_sq(xt::xtensor<double, 3>, int N);
xt::xtensor<double,2> calc_forces(vector<Particle>&, xt::xtensor<double, 3>);
void update_state(vector<Particle>&, double&, xt::xtensor<double,2>&);
void write_state(FILE*, vector<Particle>&, double t);



double length = 10;

int main2()
{
    FILE *fp_test = fopen("test.tsv", "w");
    
    for (int part_num = 1; part_num<=100; part_num++)
    {
        printf("Particle number: %i\n", part_num);
    
    	//Create vector to store Particles
        vector<Particle> particles(part_num);
    
        /*---------------------------------*/
    
        // Initialize particles
        int grid_size = 1000;
    
        for (unsigned int i=0; i<part_num; i++)
        {
    		particles[i] = Particle(i, 1, xt::zeros<double>({dim}));
        }
        
        /*---------------------------------*/
        FILE *fp = fopen("test.xyz", "w");
        
        std::clock_t start;
        start = std::clock();
        xt::xtensor<double, 3> sep_vecs = calc_sep_vecs(particles);
        fprintf(fp_test,"%f\t",(std::clock()-start)/(double)(CLOCKS_PER_SEC));
    
        start = std::clock();
        xt::xtensor<double, 2> forces = calc_forces(particles, sep_vecs);
        fprintf(fp_test,"%f\t",(std::clock()-start)/(double)(CLOCKS_PER_SEC));
    
        start = std::clock();
        write_state(fp, particles, 0);
        fprintf(fp_test,"%f\t",(std::clock()-start)/(double)(CLOCKS_PER_SEC));
        
        start = std::clock();
        calc_dist_sq(sep_vecs, part_num);
        fprintf(fp_test,"%f\n",(std::clock()-start)/(double)(CLOCKS_PER_SEC));
        
        /*
        start = std::clock();
        double dt = 0.1;
        update_state(particles, dt, forces);
        fprintf(fp_test,"%f\n",(std::clock()-start)/(double)(CLOCKS_PER_SEC));
        */
    }
    return 0;
}

int main()
{
	
	//Open input file
	ifstream infile("input_params.txt");
    
	//Determine how many particles there are
	unsigned int part_num;
	//Get the timestep and duration
	double dt, duration;
    double stdev;
    double m = 3;
    int modulus;
    double target_temp = 1;
    
    string input_type;
    string value;
    
	while ((infile >> input_type >> value))
	{
		if (input_type == "PART_NUM") part_num = stoi(value);
        else if (input_type == "STAN_DEV") stdev = stod(value);
        else if (input_type == "TIMESTEP") dt = stod(value);
        else if (input_type == "DURATION") duration = stod(value);
        else if (input_type == "LENGTH") length = stod(value);
        else if (input_type == "MODULUS") modulus = stoi(value);
        else if (input_type == "TEMP") target_temp = stoi(value);
	}
    printf("Particle number: %i, Width: %f\n", part_num, stdev);
    
	//Create vector to store Particles
    vector<Particle> particles(part_num);
	infile.close();
    
    /*---------------------------------*/
    
    // Initialize particles
    int grid_size = 1000;
    
    for (unsigned int i=0; i<part_num; i++)
    {
		particles[i] = Particle(i, m, xt::zeros<double>({dim}));
    }
    
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 mt(seed); // Seeds Mersenne Twister with Device RNG
    
    {
        double *density = new double[grid_size];
        double arg_scale = 1/grid_size;
        double norm = 0;
        for (unsigned int i = 0; i<grid_size; i++)
        {
            density[i] = 1;
            norm += density[i];
        }
        
        //Normalize position density
        for (unsigned int i = 0; i<grid_size; i++)
        {
            density[i] /= norm;
        }
        
        // Sample locations
        vose *uniform_sampler = new vose(density, grid_size, mt);
        xt::xtensor<double, 1> temp_container;
        for (auto &part : particles)
        {
            temp_container = part.getPos();
            
            for (int i = 0; i<3; i++)
                temp_container[i] = 
                    uniform_sampler->alias_method()*length/grid_size;
            
            part.setPos(temp_container);
        }
    }
    
    {
        // Define momentum density
        double *mom_density = new double[grid_size];
        double Temperature = 100;
        double mass = 4;
        double arg_scale = -1/(2*mass*Temperature);
    
        double norm = 0;
        for (unsigned int i = -floor(grid_size/2); i<floor(grid_size/2); i++)
        {
            mom_density[i] = exp(arg_scale*i*i);
            norm += mom_density[i];
        }
    
        //Normalize momentum density
        for (unsigned int i = 0; i<grid_size; i++)
        {
            mom_density[i] /= norm;
        }
    
        // TODO Make this use STDEV
    
        vose *gaussian_sampler = new vose(mom_density, grid_size, mt);
        double max_mom = 1;
        unsigned int j;
        xt::xtensor<double, 1> temp_container;
        for (auto &part : particles)
        {   
            temp_container = part.getPos();
            for (int i = 3; i<6; i++)
            {
                j = gaussian_sampler->alias_method(); //momentum index
                //convert j to momentum
                temp_container[i] = (j-float(grid_size)/2)*max_mom/grid_size;
                //temp_container[i] = 0;
            }
            part.setPos(temp_container);
        }
    }
    
	//Create output file
	FILE *fp = fopen("output.xyz", "w");
    FILE *fp_energies = fopen("energies.tsv", "w");
    
	unsigned int i = 0;
    std::clock_t start;
    start = std::clock();
    
    xt::xtensor<double, 3> sep_vecs = calc_sep_vecs(particles);
    xt::xtensor<double, 2> forces = calc_forces(particles, sep_vecs);
    write_state(fp, particles, 0);
    
    
    double temp = 0;
    double pot = 0;
    double scale = 1;
    
    //Run the simulation
    printf("%f\t%f\n",0., temp);
    for (double t = dt; t <= duration; t+=dt)
	{
        update_state(particles, dt, forces);
        
        temp = temperature(particles)[0];
        scale = sqrt(1+(target_temp/temp - 1)/2);
        scale_mom(particles, scale);
        
        sep_vecs = calc_sep_vecs(particles);
        pot = potential(calc_dist_sq(sep_vecs, part_num), particles)[0];
        
        fprintf(fp_energies, "%f\t", temp);
        fprintf(fp_energies, "%f\t", pot);
        fprintf(fp_energies, "%f\n", temp+pot);
        
        printf("%f\t%f\n",t,temp);
        //scale_mom(particles, sqrt(1000/total_vel_sq));

        if (i%modulus == 0)
            write_state(fp, particles, t);
		i++;
	}
	
    //Clean up and close
	fclose(fp);
	return 0;
}

xt::xtensor<double, 0> temperature(vector<Particle>& particles)
{
    int N = static_cast<int>(particles.size());
    xt::xtensor<double, 0> sum = 0;
    for (auto &part : particles)
    {
        sum += xt::norm_sq(xt::view(part.pos, xt::range(3,6)));
    }
    return sum/N;
}

xt::xtensor<double, 0> potential(xt::xtensor<double,2> distances, vector<Particle> &particles)
{
    int N = static_cast<int>(particles.size());
    xt::xtensor<double, 0> sum = 0;
    
    auto dist = xt::view(distances,0,0);
    
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            if (i==j) continue;
            dist = xt::view(distances,i,j);
            sum += 1/pow(dist,6)*(1/pow(dist,6)+1);
        }
    }
    
    return sum/N;
}

void scale_mom(vector<Particle>& particles, double scale)
{
    for (auto &part : particles)
    {
        xt::view(part.pos, xt::range(3,6)) *= scale;
    }
}

xt::xtensor<double, 3> calc_sep_vecs(vector<Particle>& particles)
{
    int N = static_cast<int>(particles.size());
    
    //It's important to initialize like this
    // so that elements can be set using index
    xt::xtensor<double, 3> sep_vecs = xt::zeros<double>({N, N, 3});
    
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<i; j++)
        {   
            xt::view(sep_vecs, i, j) =
                xt::view(particles[i].pos, xt::range(0,dim/2)) -
                xt::view(particles[j].pos, xt::range(0,dim/2));
        }
    }
    
    for (int i = 0; i<N; i++)
    {
        for (int j = i; j<N; j++)
        {
            xt::view(sep_vecs,i,j) = -xt::view(sep_vecs,j,i);
        }
    }
    
    return sep_vecs;
}

xt::xtensor<double, 2> calc_dist_sq(xt::xtensor<double, 3> sep_vecs, int N)
{
    xt::xtensor<double, 2> dist_sq = xt::zeros<double>({N, N});
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            xt::view(dist_sq, i, j) =
                xt::norm_sq(xt::view(sep_vecs, i, j, xt::all()));
        }
    }
    return dist_sq;
}

xt::xtensor<double, 2> calc_forces(
    vector<Particle> &particles, xt::xtensor<double, 3> sep_vecs)
{
    int N = sep_vecs.shape()[0];
    
    xt::xtensor<double, 2> dist_sq = calc_dist_sq(sep_vecs, N);
    
    xt::xtensor<double, 2> forces = xt::zeros<double>({N,3});
    xt::xtensor<double, 1> sep_vec = xt::zeros<double>({3});
    
    auto dist_sq_local = xt::view(dist_sq,0,0);
    auto dist_6 = xt::view(dist_sq,0,0);
    auto dist_8 = xt::view(dist_sq,0,0);
    
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            if (i==j) continue;
            sep_vec = xt::view(sep_vecs,i,j);
            dist_sq_local = xt::view(dist_sq,i,j); 
            dist_6 = dist_sq_local*dist_sq_local*dist_sq_local;
            dist_8 = dist_6*dist_sq_local;
            xt::view(forces,i) += -6*sep_vec/dist_8*(-2/dist_6+1);
        }
        
        /*
        double temp = temperature(particles)[0];
        if (temp <= 100)
            xt::view(forces,i) -= 0.1*xt::view(particles[i].pos, xt::range(3,6));
        else
            xt::view(forces,i) -= 100*xt::view(particles[i].pos, xt::range(3,6));
        */
    }
    
    return forces;
}

void update_state(vector<Particle> &particles, double &dt,
    xt::xtensor<double,2> &old_forces)
{
    unsigned int part_num = particles.size();
    
    //Update spatial positions
    for (auto &part : particles)
    {
        xt::view(part.pos, xt::range(0,dim/2)) += 
            dt/part.mass*(xt::view(part.pos, xt::range(dim/2,dim))
                + dt/2*xt::view(old_forces, part.getID()));
    }
    
    //Update force
    xt::xtensor<double, 3> sep_vecs = calc_sep_vecs(particles);
    
    bool reflecting_bcs = false;
    bool periodic_bcs = !reflecting_bcs;
    
    //Change seperation vectors for minimal image (like Umklapp scattering)
    if (periodic_bcs)
    {
        auto min_image = xt::filter(sep_vecs, abs(sep_vecs)>0.5*length);
        min_image -= length*xt::sign(min_image);
    }
    
    xt::xtensor<double, 2> new_forces = calc_forces(particles, sep_vecs);
    
    //Calculate momenta at +.5dt
    for (auto &part : particles)
    {
        xt::view(part.pos, xt::range(dim/2,dim)) +=
            0.5*dt*xt::view(old_forces+new_forces, part.getID());
    }
    
    if (reflecting_bcs == true)
    {
        //Flip momenta for particles outside domain
        for (auto &part : particles)
        {
            for (int i = 0; i<3; i++)
            {
                if (part.pos[i]<0 and part.pos[i+3]<0)
                {
                    part.pos[i+3] = -0.2*part.pos[i+3];
                    part.pos[i] = 0;
                }
                if (part.pos[i]>length and part.pos[i+3]>0)
                {
                    part.pos[i+3] = -0.2*part.pos[i+3];
                    part.pos[i] = length;
                }
            }
        }
    }
    
    if (periodic_bcs)
    {
        //Pac-Man shift particles (Periodic BCs)
        for (auto &part : particles)
        {
            for (int i = 0; i<3; i++)
            {
                /*
                if (part.pos[i]<0)
                    part.pos[i] += length;
                if (part.pos[i]>length)
                    part.pos[i] -= length;
                */
                part.pos[i] = fmod(fmod(part.pos[i],length)+length,length);
            }
        }
    }
    
    old_forces = new_forces;
}

void write_state(FILE *fp, vector<Particle> &particles, double t)
{
    fprintf(fp, "%lu\n", particles.size());
    fprintf(fp, "Time: %f\n", t);
    for (auto &part : particles)
	{
		xt::xtensor<double, dim> pos = part.pos;
        fprintf(fp, "C %f %f %f %f %f %f\n",
            pos[0], pos[1], pos[2], pos[3], pos[4], pos[5]);
	}
}

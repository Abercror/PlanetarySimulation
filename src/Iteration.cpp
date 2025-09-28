#include "ParticleClass.hpp"

template <typname T>

std::vector<Particle> updatePositions(const std::vector<Particle> &particles, const T &dt, const T &integrationMethod){
    for(Particle particle : particles){
        switch(integrationMethod){
            case 1:
                particle.euler(dt, particles, G);
            case 2:
                particle.eulerCromer(dt, particles, G);
            case 3:
                particle.verlet(dt, particles, G);
            case 4: 
                particle.leapfrog(dt, particles, G);
            case 5:
                particle.yoshida4thOrder(dt, particles, G);
            case 6:
                particle.rungeKatta(dt, particles, G);
            case 7:
                particle.vefrl(dt, particles, G);
            case 8:
                particle.pefrl(dt, particles, G);
        }
    }
}
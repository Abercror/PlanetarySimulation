#include "vsop87a_full_velocities.hpp"
#include "vsop87_full.hpp"
#include "ParticleClass.hpp"
#include "VectorClass.hpp"
#include <vector>
#include <string>

template <typename T>
struct PlanetConstants {
    const char* name;
    T mass;
    T radius;
};

template <typename T>
void planetProperties(int jd, std::vector<Particle<T>> &planets)

template <typename T>
std:vector<Particle<T>> planetInitialisation(int jd)



template <typename T>
void planetProperties(int jd, std::vector<Particle<T>> &planets) {

    PlanetConstants<T> constants[] = {
        { "Sun",     1.9885e30, 695700e3 },
        { "Mercury", 3.3011e23, 2439.7e3 },
        { "Venus",   4.8675e24, 6051.8e3 },
        { "Earth",   5.97237e24, 6371.0e3 },
        { "Mars",    6.4171e23, 3389.5e3 },
        { "Jupiter", 1.8982e27, 69911e3 },
        { "Saturn",  5.6834e26, 58232e3 },
        { "Uranus",  8.6810e25, 25362e3 },
        { "Neptune", 1.02413e26, 24622e3 }
    };

    for (size_t i = 0; i < planets.size(); ++i) {
        double tempPos[3];
        double tempVel[3];

        switch (i) {
            case 0: 
                vsop87a_full::getSun(jd, tempPos);
                vsop87a_full_velocities::getSun(jd, tempVel); 
                break;
            case 1: 
                vsop87a_full::getMercury(jd, tempPos);
                vsop87a_full_velocities::getMercury(jd, tempVel); 
                break;
            case 2: 
                vsop87a_full::getVenus(jd, tempPos); 
                vsop87a_full_velocities::getVenus(jd, tempVel); 
                break;
            case 3: 
                vsop87a_full::getEarth(jd, tempPos); 
                vsop87a_full_velocities::getEarth(jd, tempVel); 
                break;
            case 4: 
                vsop87a_full::getMars(jd, tempPos); 
                vsop87a_full_velocities::getMars(jd, tempVel); 
                break;
            case 5: 
                vsop87a_full::getJupiter(jd, tempPos); 
                vsop87a_full_velocities::getJupiter(jd, tempVel); 
                break;
            case 6: 
                vsop87a_full::getSaturn(jd, tempPos); 
                vsop87a_full_velocities::getSaturn(jd, tempVel); 
                break;
            case 7: 
                vsop87a_full::getUranus(jd, tempPos); 
                vsop87a_full_velocities::getUranus(jd, tempVel); 
                break;
            case 8: 
                vsop87a_full::getNeptune(jd, tempPos); 
                vsop87a_full_velocities::getNeptune(jd, tempVel); 
                break;
        }

        planets[i].m_mass = constants[i].mass;
        planets[i].m_radius = constants[i].radius;
        planets[i].m_position = Vec3(tempPos[0], tempPos[1], tempPos[2]);
        planets[i].m_velocity = Vec3(tempVel[0], tempVel[1], tempVel[2]);
    }
}

template <typename T>
std::vector<Particle<T>> planetInitialisation(int jd) {

    std::vector<Particle<T>> planets(9);
    planetProperties(jd, planets);
    return planets;
}
#include "vsop87a_full_velocities.h"
#include "vsop87_full.h"
#include "ParticleClass.h"
#include "VectorClass.h"
#include <vector>
template <typename T>

void planetProperties(int &jd, std::vector<Particles> &planets);

struct PlanetConstants {
    const char* name;
    T mass;
    T radius;
};

std::vector<Particles> planetInitialisation(int &jd){

    Particle sun;
    Particle mercury;
    Particle venus;
    Particle earth;
    Particle mars;
    Particle jupiter;
    Particle saturn;
    Particle uranus;
    Particle neptune;

    std::vector<Particles> planets = {sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune};

    planetProperties(jd, planets);

    return planets
}

void planetProperties(int &jd, std::vector<Particles> &planets){

    PlanetConstants constants[] = {
        { "Mercury", 3.3011e23, 2439.7e3 },
        { "Venus",   4.8675e24, 6051.8e3 },
        { "Earth",   5.97237e24, 6371.0e3 },
        { "Mars",    6.4171e23, 3389.5e3 },
        { "Jupiter", 1.8982e27, 69911e3 },
        { "Saturn",  5.6834e26, 58232e3 },
        { "Uranus",  8.6810e25, 25362e3 },
        { "Neptune", 1.02413e26, 24622e3 }
    };

    for(Particle planet : planets){
        double tempPos[3];
        double tempVel[3];

        vsop87a_full::get<planet>(jd, tempPos);
        vsop87a_full_velocities::get<planet>(jd, tempVel);

        Particle planet(
            m_mass = constants[i][1];
            m_radius = constant[i][2];
            m_position = Vec3(tempPos[0], tempPos[1], tempPos[2])
            m_velocity = Vec3(tempVel[0], tempVel[1], tempVel[3])
            m_acceleration = Vec3()
        );
    }
}


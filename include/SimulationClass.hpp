#pragma once
#include "Vec3Class.hpp"
#include "ParticleClass.hpp"
#include "NumericalMethodsClass.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

template <typename T> constexpr T two  = T(2);

template <typename T>
class Simulation{
    private:
        T m_totalEnergy;
        T m_totalPotentialEnergy;
        T m_totalKineticEnergy;
        T m_time;
        std::unordered_map<std::string, Vec3<T>> m_planetaryPositions;
        std::unordered_map<std::string, Vec3<T>> m_planetaryVelocities;
        std::vector<Particles<T>> m_bodies;
        functionPtr<T> m_method;
        
    public:

        Simulation(T tE, T pE, T kE, T t, std::unordered_map<std::string, Vec3<T>> pos, std::unordered_map<std::string, Vec3<T>> vel, std::vector<Particle<T>> bodies, std::function<void(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G)> method);

        void updateEnergy(const T &G);

        void updatePositionVelocity(const T &G, const T &dt);

        void numericalMethodChoice(std::string &method);

        void run();

};

template <typename T>
Simulation<T>::Simulation(T tE, T pE, T kE, T t, std::unordered_map<std::string, Vec3<T>> pos, std::unordered_map<std::string, Vec3<T>> vel, std::vector<Particle<T>> bodies, std::function<void(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G)> method): m_totalEnergy(tE), m_totalPotentialEnergy(pE), m_totalKineticEnergy(kE), m_time(t), m_planetaryPositions(pos), m_planetaryVelocities(vel), m_bodies(bodies), m_method(method){}

template <typename T>
void Simulation<T>::updateEnergy(const T &G){
    T kineticEnergy;
    T potentialEnergy;
    T totalEnergy;

    for(Particles<T> particle: m_bodies){
        Vec3<T> velocity = particle.m_velocity;
        velocity *= velocity;
        T velocitySquaredMagnitude = Vec3<T>::magnitude(velocity)

        kineticEnergy += particle.m_mass * velocitySquaredMagnitude / two<T>;
        
        for(Particles other: particles){
            if particle != other{
                T distanceMag = Vec3<T>::distanceMagnitude(particle.m_position, other.m_position)

                potentialEnergy += (G * particle.m_mass * other.m_mass) / distanceMag;
            }
        }
    }

    totalEnergy = potentialEnergy + kineticEnergy;

    m_totalEnergy = totalEnergy;
    m_totalKineticEnergy = kineticEnergy;
    m_totalPotentialEnergy = potentialEnergy;
}

template <typename T>
void Simulation<T>::numericalMethodChoice(std::string &method){
        static const std::unordered_map<std::string, std::function<void(Particle<T> &, const std::vector<Particle<T>> &, const T &, const T &)>> methods = {
        {"euler", euler},
        {"eulerCromer", eulerCromer},
        {"verlet", verlet},
        {"leapfrog", leapfrog},
        {"yoshida4thOrder", yoshida4thOrder},
        {"rungeKutta", rungeKutta},
        {"vefrl", vefrl},
        {"pefrl", pefrl}
    };
    
    m_method = methods.at(method);
}

template <typename T>
void Simulation<T>::updatePositionVelocity(const T &G, const T &dt){

    for(Particles<T> particle: m_bodies){
        (particle.*m_method)(particle, m_bodies, dt, G);
        m_planetaryPosition[particle.m_name] = particle.m_position;
        m_planetaryVelocity[particle.m_name] = particle.m_velocity;
    }
}
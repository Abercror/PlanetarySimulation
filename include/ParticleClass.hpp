#pragma once
#include "Vec3Class.hpp"
#include "NumericalMethodsClass.hpp"
#include <vector>
#include <cmath>
#include <string>

template <typename T> constexpr T two  = T(2);
template <typename T>

class Particle: public NumericalMethods{
    private:
        std::string m_name;
        T m_mass;
        T m_radius;
        Vec3<T> m_position;
        Vec3<T> m_velocity;

    public:
        Particle(std::string name, T m, T r, Vec3<T> p, Vec3<T> v);

        Vec3<T> gravitationalAcceleration(Particle<T> &p, const Vec3<T> &thisPosition, const std::vector<Particle<T>> &particles, const T &G);

};

//Constructor
template <typename T>
Particle<T>::Particle(
        std::string name,
        T m, 
        T r, 
        Vec3<T> p, 
        Vec3<T> v
        ): m_name(name), m_mass(m), m_radius(r), m_position(p), m_velocity(v) {}

//Newtonian Gravitational Acceleration
template <typename T>
Vec3<T> Particle<T>::gravitationalAcceleration(Particle<T> &p, const Vec3<T> &thisPosition, const std::vector<Particle<T>> &particles, const T &G){
    Vec3<T> direction;
    T distanceMag;
    T epsilon;
    Vec3<T> acceleration;

    for(Particle<T> particle : particles){
        if p != particle {
            distanceMag = Vec3::distanceMagnitude(thisPosition, particle.m_position);
            direction = Vec3::directionVector(thisPosition, particle.m_position);
            epsilon = distanceMag / two<T>;
            acceleration += (G * particle.m_mass) * (direction / std::sqrt((std::pow(distanceMag, two<T>) + std::pow(epsilon, two<T>)))); 
        }
    }
    return acceleration
}


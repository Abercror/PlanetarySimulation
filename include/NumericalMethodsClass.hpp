#pragma once
#include "Vec3Class.hpp"
#include "ParticleClass.hpp"
#include <vector>
#include <cmath>

template <typename T> constexpr T one  = T(1);
template <typename T> constexpr T two  = T(2);
template <typename T> constexpr T six  = T(6);

template <typename T>
class NumericalMethods{
    public:
        void euler(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void eulerCromer(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void verlet(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void leapfrog(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void yoshida4thOrder(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void rungeKutta(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void pefrl(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

        void vefrl(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G);

};

//Forwards Euler 

template <typename T> 
void NumericalMethods<T>::euler(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 acceleration;
    p.m_position += p.m_velocity * dt;
    acceleration = p.gravitationalAcceleration(p, p.m_position, particles, G);
    p.m_velocity += acceleration * dt;
}

//Backwards Euler
template <typename T> 
void NumericalMethods<T>::eulerCromer(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){           
    Vec3 acceleration;
    acceleration = p.gravitationalAcceleration(p, p.m_position, particles, G);
    p.m_velocity += acceleration * dt;
    p.m_position += p.m_velocity * dt;
}

//Verlet
template <typename T> 
void NumericalMethods<T>::verlet(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 position0 = p.m_position;
    Vec3 position1;
    Vec3 velocity0 = p.m_velocity;
    Vec3 velocity1;
    Vec3 a0;
    Vec3 a1;

    a0 = p.gravitationalAcceleration(p, position0, particles, G);
    position1 = position0 + velocity0 * dt + a0 * dt * dt / two<T>;
    a1 = p.gravitationalAcceleration(p, position1, particles, G);
    velocity1 = velocity0 + (a0 + a1) * dt / two<T>;

    p.m_position = position1;
    p.m_velocity = velocity1;
}

//Leapfrog 
template <typename T> 
void NumericalMethods<T>::leapfrog(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 velocity0 = p.m_velocity;
    Vec3 velocity1;
    Vec3 velocity2;
    Vec3 position0 = p.m_position;
    Vec3 position1;
    Vec3 a0;
    Vec3 a1;

    a0 = p.gravitationalAcceleration(p, position0, particles, G);
    velocity1 = velocity0 + a0 * dt;
    position1 = position0 + velocity1 * dt / two<T>;
    a1 = p.gravitationalAcceleration(p, position1, particles, G);
    velocity2 = velocity1 + a1 * dt / two<T>;

    p.m_position = position1;
    p.m_velocity = velocity2;
}

//Yoshida Leapfrog
template <typename T> 
void NumericalMethods<T>::yoshida4thOrder(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 velocity0 = p.m_velocity;
    Vec3 velocity1;
    Vec3 velocity2;
    Vec3 velocity3;
    Vec3 position0 = p.m_position;
    Vec3 position1;
    Vec3 position2;
    Vec3 position3;
    Vec3 position4;
    Vec3 acceleration;

    constexpr T c_14 = T(0.6756035959798289);
    constexpr T c_23 = T(-0.1756035959798288);
    constexpr T d_13 = T(1.3512071919596578);
    constexpr T d_24 = T(-1.7024143839193153);

    //1
    position1 = position0 + c_14 * velocity0 * dt;
    acceleration = p.gravitationalAcceleration(p, position1, particles, G);
    velocity1 = velocity0 + d_13 * acceleration * dt;

    //2
    position2 = position1 + c_23 * velocity1 * dt;
    acceleration = p.gravitationalAcceleration(p, position2, particles, G);
    velocity2 = velocity1 + d_2 * acceleration * dt;

    //3
    position3 = position2 + c_23 * velocity2 * dt;
    acceleration = p.gravitationalAcceleration(p, position3, particles, G);
    velocity3 = velocity2 + d_2 * acceleration * dt;

    //4
    p.m_position = position3 + c_14 * velocity3 * dt;
    p.m_velocity = velocity3;
}

//Runge Katta
template <typename T> 
void NumericalMethods<T>::rungeKutta(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 velocity0 = p.m_velocity;
    Vec3 position0 = p.m_position;

    //k1
    Vec3 a1 = p.gravitationalAcceleration(p, position0, particles, G);

    //k2
    Vec3 position1 = position0 + velocity0 * dt/two<T>;
    Vec3 a2 = p.gravitationalAcceleration(p, position1, particles, G);
    Vec3 velocity1 = velocity0 + a1 * dt/2;

    //k3 
    Vec3 position2 = position0 + velocity1 * dt/two<T>;
    Vec3 a3 = p.gravitationalAcceleration(p, position2, particles, G);
    Vec3 velocity2 = velocity0 + a2 * dt/2;

    //k4
    Vec3 position3 = position0 + velocity2 * dt;
    Vec3 a4 = p.gravitationalAcceleration(p, position3, particles, G);
    Vec3 velocity3 = velocity0 + a3 * dt;


    p.m_position = position0 + (dt / six<T>) * (velocity0 + two<T> * velocity1 + two<T> * velocity2 + velocity3);
    p.m_velocity = velocity0 + (dt / six<T>) * (a1 + two<T> * a2 + two<T> * a3 + a4);
}

//Forest-Ruth and Suzuki: https://arxiv.org/pdf/cond-mat/0110585
template <typename T> 
void NumericalMethods<T>::vefrl(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 position = p.m_position;
    Vec3 velocity = p.m_velocity;
    Vec3 acceleration;

    constexpr T xi = T(0.164498651557576);
    constexpr T lambda = T(-0.02094333910398989);
    constexpr T chi = T(1.235692651138917);

    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * xi * dt;        
    position += velocity * (one<T> - two<T> * lambda) * dt / two<T>;

    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * chi * dt;
    position += velocity * lambda * dt;

    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * (one<T> - two<T> * (xi + chi)) * dt;
    position += velocity * lambda * dt;

    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * chi * dt;
    position += velocity * (one<T> - two<T> * lambda) * dt / two<T>;

    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * xi * dt; 

    p.m_position = position;
    p.m_velocity = velocity;
}

template <typename T> 
void NumericalMethods<T>::pefrl(Particle<T> &p, const std::vector<Particle<T>> &particles, const T &dt, const T &G){
    Vec3 velocity = p.m_velocity;
    Vec3 position = p.m_position;
    Vec3 acceleration;


    constexpr T xi = T(0.1786178958448091);
    constexpr T lambda = T(-0.2123418310626054);
    constexpr T chi = T(-0.06626458266981849);

    position += velocity * xi * dt;
    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * (one<T> - two<T> * lambda) * dt / two<T>;

    position += velocity * chi * dt;
    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * lambda * dt;

    position += velocity * (one<T> - two<T> * (chi + xi)) * dt;
    aceeleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration * lambda * dt;

    position += velocity * chi * dt;
    acceleration = p.gravitationalAcceleration(p, position, particles, G);
    velocity += acceleration + (one<T> - two<T> * lambda) * dt / two<T>;

    position += velocity * xi * dt;

    p.m_position = position;
    p.m_velocity = velocity;
}
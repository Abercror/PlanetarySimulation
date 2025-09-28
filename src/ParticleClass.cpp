#include <cmath>
#include "VectorClass.hpp"
template <typename T>

class Particle{
    private:
        T m_mass;
        T m_radius;
        Vec3 m_position;
        Vec3 m_velocity;
        Vec3 m_acceleration;
        Vec3 m_momentum;
        Vec3 m_potentialEnergy;
        Vec3 m_kineticEnergy;

    public:
        //Constructor
        Particle(
            const T m, 
            const T r, 
            const Vec3 p, 
            const Vec3 v,
            const Vec3 a
            ): mass(m), radius(r), position(p), velocity(v), acceleration(a) {}
        
        //Newtonian Gravitational Acceleration
        void gravitationalAcceleration(const Vec3 &thisPosition, const std::vector<Particle> &particles, const T &G){
            Vec3 direction;
            T distanceMag;
            T epsilon;
            m_acceleration.zero();

            for(Particle particle : particles){
                if this != particle {
                    distanceMag = Vec3::distanceMagnitude(thisPosition, particle.m_position);
                    direction = Vec3::directionVector(thisPosition, particle.m_position);
                    epsilon = distanceMag / 2;

                    m_acceleration += (G * particle.m_mass) * (direction / std::sqrt((std::pow(distanceMag, 2) + std::pow(epsilon, 2)))); 
                }
            }
        }

        void kineticEnergy(){
            m_kineticEnergy = 0.5 * m_mass * m_velocity * m_velocity;
        }

        void potentialEnergy(const std::vector<Particle> &particles, const T &G){
            m_potentialEnergy.zero();
            for(Particle particle : particles){
                if this != particle {
                    m_potentialEnergy += (G * m_mass * particle.m_mass) / distanceVector(m_position, particle.m_position);
                }
            }
        }

        //Numerical Integration techniques
        //Forwards Euler 
        void euler(const T &dt, const std::vector<Particle> &particles, const T &G){
            m_position += m_velocity * dt;
            this->gravitationalAcceleration(m_position, particles, G);
            m_velocity += m_acceleration * dt;
        }

        //Backwards Euler
        void eulerCromer(const T &dt, const std::vector<Particle> &particles, const T &G){           
            this->gravitationalAcceleration(m_position, particles, G);
            m_velocity += m_acceleration * dt;
            m_position += m_velocity * dt;
        }

        //Verlet
        void verlet(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 position0 = m_position;
            Vec3 position1;
            Vec3 position2;
            Vec3 velocity0 = m_velocity;

            this->gravitationalAcceleration(position0, particles, G);
            position1 = position0 + velocity0 * dt + 0.5 * m_acceleration * dt * dt;
            this->gravitationalAcceleration(position1, particles, G);
            position2 = 2 * position1 - position0 + m_acceleration * dt * dt;

            m_position = position2;
        }

        //Leapfrog 
        void leapfrog(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 velocity0 = m_velocity;
            Vec3 velocity1;
            Vec3 position0 = m_position;
            Vec3 position1;

            this->gravitationalAcceleration(position0, particles, G);
            velocity1 = m_velocity + 0.5 * m_acceleration * dt;
            positon1 = position0 + m_velocity * dt;
            this->gravitationalAcceleration(position1, particles, G);
            m_velocity = velocity1 + 0.5 * m_acceleration * dt;
        }

        //Yoshida Leapfrog
        void yoshida4thOrder(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 velocity0 = m_velocity;
            Vec3 velocity1;
            Vec3 velocity2;
            Vec3 velocity3;
            Vec3 position0 = m_position;
            Vec3 position1;
            Vec3 position2;
            Vec3 position3;
            Vec3 position4;

            T c_14 = 0.6756;
            T c_23 = -0.1756;
            T d_13 = 1.3412;
            T d_2 = -0.3412;

            //1
            position1 = position0 + c_14 * velocity0 * dt;
            this->gravitationalAcceleration(position1, particles, G);
            velocity1 = velocity0 + d_13 * m_acceleration * dt;

            //2
            position2 = position1 + c_23 * velocity1 * dt;
            this->gravitationalAcceleration(position2, particles, G);
            velocity2 = velocity1 + d_2 * m_acceleration * dt;

            //3
            position3 = position3 + c_23 * velocity2 * dt;
            this->gravitationalAcceleration(position3, particles, G);
            velocity3 = velocity2 + d_2 * m_acceleration * dt;

            //4
            m_position = position3 + c_14 * velocity3 * dt;
            m_velocity = velocity3;
        }

        //Runge Katta
        void rungeKatta(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 velocity0 = m_velocity;
            Vec3 velocity1;
            Vec3 velocity2;
            Vec3 velocity3;
            Vec3 position0 = m_position;
            Vec3 position1;
            Vec3 position2;
            Vec3 position3;
            Vec3 position4;

            //k1
            position1 = position0 + velocity0 * dt;
            this->gravitationalAcceleration(position1, particles, G);
            velocity1 = velocity0 + m_acceleration * dt;

            //k2
            position2 = position1 + velocity1 * dt/2;
            this->gravitationalAcceleration(position2, particles, G);
            velocity2 = velocity1 + m_acceleration * dt/2;

            //k3 
            position3 = position2 + velocity2 * dt/2;
            this->gravitationalAcceleration(position3, particles, G);
            velocity3 = velocity2 + m_acceleration * dt/2;

            //k4
            position4 = position3 + velocity3 * dt;
            this->gravitationalAcceleration(position4, particles, G);
            velocity4 = velocity3 + m_acceleration * dt;

            m_position = position0 + (dt / 6) * (velocity1 + 2 * velocity2 + 2 * velocity3 + velocity4);
            m_velocity = velocity4;
        }

        //Forest-Ruth and Suzuki: https://arxiv.org/pdf/cond-mat/0110585
        void vefrl(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 velocity0 = m_velocity;
            Vec3 velocity1;
            Vec3 velocity2;
            Vec3 velocity3;
            Vec3 velocity4;
            Vec3 position0 = m_position;
            Vec3 position1;
            Vec3 position2;
            Vec3 position3;

            T xi = +0.164498651557576;
            T lambda = -0.02094333910398989;
            T chi = +1.235692651138917;

            this->gravitationalAcceleration(position0, particles, G);
            velocity1 = velocity0 + m_acceleration * xi * dt;        
            position1 = position0 + velocity1 * (1 - 2 * lambda) * dt / 2;

            this->gravitationalAcceleration(position1, particles, G);
            velocity2 = velocity1 + m_acceleration * chi * dt;
            position2 = position1 + velocity2 * lambda * dt;

            this->gravitationalAcceleration(position2, particles, G);
            velocity3 = velocity2 + m_acceleration * (1 - 2*(xi + chi)) * dt;
            position3 = velocity3 + velocity2 * lambda * dt;

            this->gravitationalAcceleration(position3, particles, G);
            velocity4 = velocity3 + m_acceleration * chi * dt;

            m_position = position3 + velocity4 * (1 - 2 * lambda) * dt / 2;
            m_velocity = velocity4 + m_position * xi * dt; 
        }

        void pefrl(const T &dt, const std::vector<Particle> &particles, const T &G){
            Vec3 velocity0 = m_velocity;
            Vec3 velocity1;
            Vec3 velocity2;
            Vec3 velocity3;
            Vec3 position0 = m_position;
            Vec3 position1;
            Vec3 position2;
            Vec3 position3;
            Vec3 position4;

            T xi = +0.178617895844809;
            T lambda = -0.2123418310626054;
            T chi = -0.06626458266981849;

            position1 = position0 + velocity0 * xi * dt;
            this->gravitationalAcceleration(position1, particles, G);
            velocity1 = velocity0 + m_acceleration * (1 - 2 * lambda) * dt / 2;

            position2 = position1 + velocity1 * chi * dt;
            this->gravitationalAcceleration(position2, particles, G);
            velocity2 = velocity1 + m_acceleration * lambda * dt;

            position3 = position2 + velocity2 * (1 - 2 * (chi + xi)) / dt;
            this->gravitationalAcceleration(position3, particles, G);
            velocity3 = velocity2 + m_acceleration * lambda * dt;

            position4 = position3 + velocity3 * chi * dt;
            this->gravitationalAcceleration(position4, particles, G);

            m_velocity = velocity3 + m_acceleration * (1 - 2 * lambda) * dt / 2;
            m_position = position4 + m_velocity * chi * dt;
        }
};
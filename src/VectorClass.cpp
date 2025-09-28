#include <cmath>
template <typename T>

class Vec3{
    private:
        T m_x = 0;
        T m_y = 0;
        T m_z = 0;
    
    public:
        Vec3(): m_x(0), m_y(0), m_z(0){}

        Vec3(T xValue, T yValue, T zValue): m_x(xValue), m_y(yValue), m_z(zValue){}

        //setter and getter functions for m_x,m_y,m_z
        void set(const T &xValue, const T &yValue, const T &zValue){
            this->m_x = xValue;
            this->m_y = yValue;
            this->m_z = zValue;
        }

        void setX(const T &xValue){ m_x = xValue;}
        void setY(const T &yValue){ m_y = yValue;}
        void setZ(const T &zValue){ m_z = zValue;}

        T getX() const {return m_x;}
        T getY() const {return m_y;}
        T getZ() const {return m_z;}
        
        //vector operations
        void zero(){
            m_x=m_y=m_z=0;
        }

        static T magnitude(const Vec3 &A){
            return std::sqrt((A.m_x * A.m_x) + (A.m_y * A.m_y) + (A.m_z * A.m_z));
        }

        void normalise(){
            T mag = magnitude();
            if(mag != 0){
                m_x /= mag;
                m_y /= mag;
                m_z /= mag;
            }
        }

        static T angleBetweenVectors(const Vec3 &A, const Vec3 &B){
            T dot = dotProduct(A, B);
            T AMag = magnitude(A);
            T BMag = magnitude(B);
            return std::acos(dot / (AMag * BMag));
        }
        
        static T dotProduct(const Vec3 &A, const Vec3 &B){
            return (A.m_x * B.m_x) + (A.m_y * B.m_y) + (A.m_z * B.m_z);
        }

        static Vec3 crossProduct(const Vec3 &A, const Vec3 &B){
            return Vec3(((A.m_y * B.m_z) - (B.m_y * A.m_z)), -((A.m_x * B.m_z) - (A.m_z * B.m_x)), ((A.m_x * B.m_y) - (B.m_x * A.m_y)));
        }

        static Vec3 grad(const T &dt, const Vec3 &A, const Vec3 &B){
            T dx = (A.m_x - B.m_x) / dt;
            T dy = (A.m_y - B.m_y) / dt;
            T dz = (A.m_z - B.m_z) / dt;
            return Vec3(dx, dy, dz);
        }

        static T div(const T &dt, const Vec3 &A, const Vec3 &B){
            Vec3 gradient = grad(dt, A, B);
            return std::sqrt((gradient.m_x * gradient.m_x) + (gradient.m_y * gradient.m_y) + (gradient.m_z * gradient.m_z));
        }

        static Vec3 curl(const T &dt, const Vec3 &A, const Vec3 &B){
            Vec3 gradient = grad(dt, A, B);
            return Vec3(crossProduct(gradient, A));
        }

        static T distanceMagnitude(const Vec3 &A, const Vec3 &B){
            T dx = A.m_x - B.m_x;
            T dy = A.m_y - B.m_y;
            T dz = A.m_z - B.m_z;
            return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
        }

        static Vec3 distanceVector(const Vec3 &A, const Vec3 &B){
            T dx = A.m_x - B.m_x;
            T dy = A.m_y - B.m_y;
            T dz = A.m_z - B.m_z;
            return Vec3(dx, dy, dz);
        }

        static Vec3 directionVector(const Vec3 &A, const Vec3 &B){
            Vec3 distanceVec = distanceVector(A, B);
            T distanceMag = distanceMagnitude(A, B);

            return distanceVec / distanceMag;
        }

        //scalar operations
        void operator+=(const T &N){
            m_x += N;
            m_y += N;
            m_z += N;
        }

        void operator-=(const T &N){
            m_x -= N;
            m_y -= N;
            m_z -= N;
        }

        void operator*=(const T &N){
            m_x *= N;
            m_y *= N;
            m_z *= N;
        }
                
        void operator/=(const T &N){
            m_x /= N;
            m_y /= N;
            m_z /= N;
        }
        
        //vector operations
        void operator+=(const Vec3 &B){
            m_x += B.m_x;
            m_y += B.m_y;
            m_z += B.m_z;
        }
                
        void operator-=(const Vec3 &B){
            m_x -= B.m_x;
            m_y -= B.m_y;
            m_z -= B.m_z;
        }

        void operator*=(const Vec3 &B){
            m_x *= B.m_x;
            m_y *= B.m_y;
            m_z *= B.m_z;
        }

        void operator/=(const Vec3 &B){
            m_x /= B.m_x;
            m_y /= B.m_y;
            m_z /= B.m_z;
        }

        //coordinate transformations
        void cartesianToCylindrical(){
            m_x = std::sqrt((m_x * m_x) + (m_y * m_y));
            m_y = std::atan2(m_y, m_x);
            m_z = m_z;
        }

        void cartesianToSpherical(){
            rho = std::sqrt((m_x * m_x) + (m_y * m_y) + (m_z * m_z));
            theta = std::atan2(std::sqrt((m_x * m_x) + (m_y + m_y)), m_z);
            phi = std::atan2(m_y, m_x);
        }

        void cylindricalToCartesian(){
            T r = m_x;
            T theta = m_y;
            T z = m_z;
            this->zero();
            m_x = r * std::cos(theta);
            m_y = r * std::sin(theta);
            m_z = z;
        }

        void sphericalToCartesian(){
            T r = m_x;
            T theta = m_y;
            T phi = m_z;
            this->zero();
            m_x = r * std::sin(theta) * std::cos(phi);
            m_y = r * std::sin(theta) * std::sin(phi);
            m_z = r * std::cos(phi);
        }
};

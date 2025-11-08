#pragma once
#include <cmath>

template <typename T>
class Vec3{
    private:
        T m_x = 0;
        T m_y = 0;
        T m_z = 0;
    
    public:
        Vec3();

        Vec3(T xValue, T yValue, T zValue);

        //setter and getter functions for m_x,m_y,m_z
        void set(const T &xValue, const T &yValue, const T &zValue);

        void setX(const T &xValue);
        void setY(const T &yValue);
        void setZ(const T &zValue);

        T getX() const;
        T getY() const;
        T getZ() const;
        
        //vector operations
        void zero();

        static T magnitude(const Vec3<T> &A);

        void normalise();

        static T angleBetweenVectors(const Vec3<T> &A, const Vec3<T> &B);
        
        static T dotProduct(const Vec3<T> &A, const Vec3<T> &B);

        static Vec3<T> crossProduct(const Vec3<T> &A, const Vec3<T> &B);

        static Vec3<T> grad(const T &dt, const Vec3<T> &A, const Vec3<T> &B);

        static T div(const T &dt, const Vec3<T> &A, const Vec3<T> &B);

        static Vec3<T> curl(const T &dt, const Vec3<T> &A, const Vec3<T> &B);

        static T distanceMagnitude(const Vec3<T> &A, const Vec3<T> &B);

        static Vec3<T> distanceVector(const Vec3<T> &A, const Vec3<T> &B);

        static Vec3<T> directionVector(const Vec3<T> &A, const Vec3<T> &B);

        //scalar operations
        void operator+=(const T &N);

        void operator-=(const T &N);

        void operator*=(const T &N);
                
        void operator/=(const T &N);
        
        //vector operations
        void operator+=(const Vec3<T> &B);
                
        void operator-=(const Vec3<T> &B);

        void operator*=(const Vec3<T> &B);

        void operator/=(const Vec3<T> &B);

        //coordinate transformations
        void cartesianToCylindrical();

        void cartesianToSpherical();

        void cylindricalToCartesian();

        void sphericalToCartesian();
};


template <typename T>
Vec3<T>::Vec3(): m_x(0), m_y(0), m_z(0){}

template <typename T>
Vec3<T>::Vec3(T xValue, T yValue, T zValue): m_x(xValue), m_y(yValue), m_z(zValue){}

//setter and getter functions for m_x,m_y,m_z
template <typename T>
void Vec3<T>::set(const T &xValue, const T &yValue, const T &zValue){
    this->m_x = xValue;
    this->m_y = yValue;
    this->m_z = zValue;
}

template <typename T>
void Vec3<T>::setX(const T &xValue){ m_x = xValue;}
template <typename T>
void Vec3<T>::setY(const T &yValue){ m_y = yValue;}
template <typename T>
void Vec3<T>::setZ(const T &zValue){ m_z = zValue;}


template <typename T>
T Vec3<T>::getX() const {return m_x;}
template <typename T>
T Vec3<T>::getY() const {return m_y;}
template <typename T>
T Vec3<T>::getZ() const {return m_z;}

//vector operations
template <typename T>
void Vec3<T>::zero(){
    m_x=m_y=m_z=0;
}

template <typename T>
T Vec3<T>::magnitude(const Vec3<T> &A){
    return std::sqrt((A.m_x * A.m_x) + (A.m_y * A.m_y) + (A.m_z * A.m_z));
}

template <typename T>
void Vec3<T>::normalise(){
    T mag = magnitude(*this);
    if(mag != 0){
        m_x /= mag;
        m_y /= mag;
        m_z /= mag;
    }
    else{
        m_x = 0;
        m_y = 0;
        m_z = 0;
    }
}

template <typename T>
T Vec3<T>::angleBetweenVectors(const Vec3<T> &A, const Vec3<T> &B){
    T dot = dotProduct(A, B);
    T AMag = magnitude(A);
    T BMag = magnitude(B);
    return std::acos(dot / (AMag * BMag));
}

template <typename T>
T Vec3<T>::dotProduct(const Vec3<T> &A, const Vec3<T> &B){
    return (A.m_x * B.m_x) + (A.m_y * B.m_y) + (A.m_z * B.m_z);
}

template <typename T>
Vec3<T> Vec3<T>::crossProduct(const Vec3<T> &A, const Vec3<T> &B){
    return Vec3(((A.m_y * B.m_z) - (B.m_y * A.m_z)), ((A.m_z * B.m_x) - (B.m_z * A.m_x)), ((A.m_x * B.m_y) - (B.m_x * A.m_y)));
}

template <typename T>
Vec3<T> Vec3<T>::grad(const T &dt, const Vec3<T> &A, const Vec3<T> &B){
    T dx = (A.m_x - B.m_x) / dt;
    T dy = (A.m_y - B.m_y) / dt;
    T dz = (A.m_z - B.m_z) / dt;
    return Vec3(dx, dy, dz);
}

template <typename T>
T Vec3<T>::div(const T &dt, const Vec3<T> &A, const Vec3<T> &B){
    Vec3<T> gradient = grad(dt, A, B);
    return std::sqrt((gradient.m_x * gradient.m_x) + (gradient.m_y * gradient.m_y) + (gradient.m_z * gradient.m_z));
}

template <typename T>
Vec3<T> Vec3<T>::curl(const T &dt, const Vec3<T> &A, const Vec3<T> &B){
    Vec3<T> gradient = grad(dt, A, B);
    return Vec3(crossProduct(gradient, A));
}

template <typename T>
T Vec3<T>::distanceMagnitude(const Vec3<T> &A, const Vec3<T> &B){
    T dx = A.m_x - B.m_x;
    T dy = A.m_y - B.m_y;
    T dz = A.m_z - B.m_z;
    return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

template <typename T>
Vec3<T> Vec3<T>::distanceVector(const Vec3<T> &A, const Vec3<T> &B){
    T dx = A.m_x - B.m_x;
    T dy = A.m_y - B.m_y;
    T dz = A.m_z - B.m_z;
    return Vec3(dx, dy, dz);
}

template <typename T>
Vec3<T> Vec3<T>::directionVector(const Vec3<T> &A, const Vec3<T> &B){
    Vec3<T> distanceVec = distanceVector(A, B);
    T distanceMag = distanceMagnitude(A, B);

    return distanceVec / distanceMag;
}

//scalar operations
template <typename T>
void Vec3<T>::operator+=(const T &N){
    m_x += N;
    m_y += N;
    m_z += N;
}

template <typename T>
void Vec3<T>::operator-=(const T &N){
    m_x -= N;
    m_y -= N;
    m_z -= N;
}

template <typename T>
void Vec3<T>::operator*=(const T &N){
    m_x *= N;
    m_y *= N;
    m_z *= N;
}
        
template <typename T>
void Vec3<T>::operator/=(const T &N){
    m_x /= N;
    m_y /= N;
    m_z /= N;
}

//vector operations
template <typename T>
void Vec3<T>::operator+=(const Vec3<T> &B){
    m_x += B.m_x;
    m_y += B.m_y;
    m_z += B.m_z;
}
    
template <typename T>
void Vec3<T>::operator-=(const Vec3<T> &B){
    m_x -= B.m_x;
    m_y -= B.m_y;
    m_z -= B.m_z;
}

template <typename T>
void Vec3<T>::operator*=(const Vec3<T> &B){
    m_x *= B.m_x;
    m_y *= B.m_y;
    m_z *= B.m_z;
}

template <typename T>
void Vec3<T>::operator/=(const Vec3<T> &B){
    m_x /= B.m_x;
    m_y /= B.m_y;
    m_z /= B.m_z;
}



//coordinate transformations
template <typename T>
void Vec3<T>::cartesianToCylindrical(){
    T r = std::sqrt((m_x * m_x) + (m_y * m_y));
    T theta = std::atan2(m_y, m_x);
    this->zero();
    m_x = r;
    m_y = theta;
}

template <typename T>
void Vec3<T>::cartesianToSpherical(){
    T r = std::sqrt((m_x * m_x) + (m_y * m_y) + (m_z * m_z));
    T theta = std::atan2(std::sqrt((m_x * m_x) + (m_y * m_y)), m_z);
    T phi = std::atan2(m_y, m_x);
    this->zero();
    m_x = r;
    m_y = theta;
    m_z = phi;
}

template <typename T>
void Vec3<T>::cylindricalToCartesian(){
    T r = m_x;
    T theta = m_y;
    T z = m_z;
    this->zero();
    m_x = r * std::cos(theta);
    m_y = r * std::sin(theta);
    m_z = z;
}

template <typename T>
void Vec3<T>::sphericalToCartesian(){
    T r = m_x;
    T theta = m_y;
    T phi = m_z;
    this->zero();
    m_x = r * std::sin(theta) * std::cos(phi);
    m_y = r * std::sin(theta) * std::sin(phi);
    m_z = r * std::cos(theta);
}
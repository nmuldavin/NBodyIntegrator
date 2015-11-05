//
//  Vector.hpp
//  Header file for vector class. Data is stored in a private array of double
//  precision floating point numbers. Additionally, the length of the Vector
//  is stored as an integer
//
//  Created by Noah Muldavin
//  Reed College 2013
/******************************************************************************/


#include<iostream>
#include <cmath>


#ifndef Vector_Vector_hpp
#define Vector_Vector_hpp


class Vector
{
private:
    double* data;       // contains vector entries
    int size;           // length of vector
    
public:
    Vector(const Vector& otherVector);
    Vector(int length);
    ~Vector();
    int length() const;
    double& operator[](int i);
    double Read(int i) const;
    Vector& operator=(const Vector& otherVector);
    Vector& operator=(const double& a);
    Vector operator+() const;
    Vector operator-() const;
    Vector operator-(const Vector& otherVector) const;
    Vector& operator-=(const Vector& otherVector);
    Vector operator+(const Vector& v1) const;
    Vector& operator+=(const Vector& otherVector);
    Vector operator*(double a) const;
    Vector operator/(double a) const;
    Vector& operator*=(const double& a);
    Vector& operator/=(const double& a);
    double operator*(const Vector& otherVector) const;
    double norm() const;
    friend std::ostream& operator<<(std::ostream& output, const Vector& z);
};

#endif


/******************************************************************************/
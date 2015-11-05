//
//  Vector.cpp
//  Defines constructors, overloaded operators, and destructors for a Vector
//  class. It is designed to function intuitively according to normal vector
//  operations
//
//  Created by Noah Muldavin
//  Reed Colege 2013
/******************************************************************************/


#include "Vector.hpp"


//  Copy constructor. Allocates space and assings all entries the values of
//  another Vector.

Vector::Vector(const Vector& otherVector)
{
    size = otherVector.length();
    data = new double [size];
    for (int i = 0; i < size; i++)
    {
        data[i] = otherVector.data[i];
    }
}


//  Default constructor. Takes an integer length as an argument,
//  allocates space and sets all entries to 0.0

Vector::Vector(int length)
{
    size = length;
    data = new double [size];
    for (int i = 0; i < size; i++)
    {
        data[i] = 0.0;
    }
}


//  Overridden destructor, deletes data when the instance of the Vector
//  goes out of use

Vector::~Vector()
{
    delete[] data;
}


//  Public method to access the length of the Vector

int Vector::length() const
{
    return size;
}


//  Overloading square brackets for assignment

double& Vector::operator[](int i)
{
    return data[i];
}


//  Read-only variant for those times when you don't want to mess
//  with data

double Vector::Read(int i) const
{
    return data[i];
}


//  Vector-Vector Assignment operator. Assigns every element to the
//  corresponding entry of another vector

Vector& Vector::operator=(const Vector& otherVector)
{
    for (int i=0; i<size; i++)
    {
        data[i] = otherVector.data[i];
    }
    return *this;
}


//  double to Vector assignment operator. Assigns all entries of the vector
//  to a single double value

Vector& Vector::operator=(const double& a)
{
    for (int i=0; i<size; i++)
    {
        data[i] = a;
    }
    return *this;
}


// Unary + operator

Vector Vector::operator+() const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = data[i];
    }
    return v;
}


//  Unary - operator

Vector Vector::operator-() const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = -data[i];
    }
    return v;
}


//  Binary + operator

Vector Vector::operator+(const Vector& otherVector) const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = data[i]+otherVector.data[i];
    }
    return v;
}


// += operator

Vector& Vector::operator+=(const Vector& otherVector)
{
    for (int i=0; i<size; i++)
    {
        data[i] = data[i] + otherVector.data[i];
    }
    return *this;
}


// Binary - operator

Vector Vector::operator-(const Vector& otherVector) const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = data[i] - otherVector.data[i];
    }
    return v;
}


//  -= operator

Vector& Vector::operator-=(const Vector& otherVector)
{
    for (int i=0; i<size; i++)
    {
        data[i] = data[i] - otherVector.data[i];
    }
    return *this;
}


//  Overridden * operator for scalar multiplication

Vector Vector::operator*(double a) const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = a*data[i];
    }
    return v;
}


//  Scalar division

Vector Vector::operator/(double a) const
{
    Vector v(size);
    for (int i=0; i<size; i++)
    {
        v[i] = data[i]/a;
    }
    return v;
}


//  *= operator for scalars

Vector& Vector::operator*=(const double& a)
{
    for (int i=0; i<size; i++)
    {
        data[i] = a*data[i];
    }
    return *this;
}


//  /= operator

Vector& Vector::operator/=(const double& a)
{
    for (int i=0; i<size; i++)
    {
        data[i] = data[i]/a;
    }
    return *this;
}


//  Overidden * operator for inner product

double Vector::operator*(const Vector& otherVector) const
{
    double out = 0.0;
    for (int i=0; i<size; i++)
    {
        out += data[i]*otherVector.data[i];
    }
    return out;
}


//  norm() method to return absolute value

double Vector::norm() const
{
    double norm_val, sum = 0.0;
    for (int i=0; i<size; i++)
    {
        sum += data[i]*data[i];
    }
    norm_val = sqrt(sum);
    return norm_val;
}


//  ostream operator << prints all elements separated by spaces

std::ostream& operator<<(std::ostream& output, const Vector& z)
{
    int n = z.size;
    for (int i = 0; i < n; i++)
    {
        output << z.data[i] << " ";
    }
    return output;
}


/******************************************************************************/
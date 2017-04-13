#ifndef PARTICLE_H
#define PARTICLE_H
#include <cmath>

namespace xiaohai
{
  class Particle
  {
  public:
    Particle():xiaohaix__(0.), xiaohaiy__(0.), xiaohaiz__(0.){}
    Particle(const double, const double, const double);
    virtual const double GetX() const;
    virtual const double GetY() const;
    virtual const double GetZ() const;
    virtual void SetX(const double);
    virtual void SetY(const double);
    virtual void SetZ(const double);
    virtual const double GetPhi() const;
    virtual const double GetPt() const;
    virtual ~Particle(){};

  private:
    double xiaohaix__;
    double xiaohaiy__;
    double xiaohaiz__;
  };/// class
}/// namespace

namespace xiaohai{
  Particle::Particle(const double x, const double y, const double z)
    :xiaohaix__(x), xiaohaiy__(y), xiaohaiz__(z) { }

  const double Particle::GetX() const
  {
    return this->xiaohaix__;
  }
  const double Particle::GetY() const
  {
    return this->xiaohaiy__;
  }
  const double Particle::GetZ() const
  {
    return this->xiaohaiz__;
  }

  void Particle::SetX(const double x)
  {
    this->xiaohaix__ = x;
  }
  void Particle::SetY(const double y)
  {
    this->xiaohaiy__ = y;
  }
  void Particle::SetZ(const double z)
  {
    this->xiaohaiz__ = z;
  }

  const double Particle::GetPhi() const
  {
    return atan2(this->xiaohaiy__, this->xiaohaix__);
  }
  const double Particle::GetPt() const
  {
    return sqrt(this->xiaohaix__*this->xiaohaix__
                + this->xiaohaiy__*this->xiaohaiy__);
  }

  const double GetEpsilon2(double numerator1, double numerator2, double denominator)
  {
    numerator1 = numerator1 * numerator1;
    numerator2 = numerator2 * numerator2;
    double numerator = sqrt(numerator1 + numerator2);
    return numerator / denominator;
  }

  const double GetEpsilon3(double numerator1, double numerator2, double denominator)
  {
    numerator1 = numerator1 * numerator1;
    numerator2 = numerator2 * numerator2;
    double numerator = sqrt(numerator1 + numerator2);
    return numerator / denominator;
  }
  
}

#endif

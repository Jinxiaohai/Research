/**
 * @file   Particle.h
 * @author xiaohai <xiaohaijin@outlook.com>
 * @date   Wed Mar 29 20:22:32 2017
 * 
 * @brief  particle class.
 * 
 * 
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <cmath>

using namespace std;

namespace xiaohai
{
  class Particle
  {
    friend istream &operator>>(istream &, Particle &);
  public:
  Particle(int Id, double Px, double Py, double Pz, double Mass,
           double X, double Y, double Z, double Time)
    :xiaohaiid__(Id), xiaohaipx__(Px), xiaohaipy__(Py), xiaohaipz__(Pz),
      xiaohaimass__(Mass), xiaohaix__(X), xiaohaiy__(Y), xiaohaiz__(Z),
      xiaohaitime__(Time){}
  Particle():xiaohaiid__(0), xiaohaipx__(0), xiaohaipy__(0), xiaohaipz__(0),
      xiaohaimass__(0), xiaohaix__(0), xiaohaiy__(0), xiaohaiz__(0),
      xiaohaitime__(0){}
    ~Particle(){}
    Particle Boost(double, double, double);
    Particle &operator=(Particle particle);
    const int GetId() const;
    const double GetPx() const;
    const double GetPy() const;
    const double GetPz() const;
    const double GetMass() const;
    const double GetX() const;
    const double GetY() const;
    const double GetZ() const;
    const double GetTime() const;
    const double GetEnergy() const;
    const double GetRapidity() const;
    const double GetPt() const;
    const double GetBetaX() const;
    const double GetBetaY() const;
    const double GetBetaZ() const;
    
  private:
    int xiaohaiid__;
    double xiaohaipx__;
    double xiaohaipy__;
    double xiaohaipz__;
    double xiaohaimass__;
    double xiaohaix__;
    double xiaohaiy__;
    double xiaohaiz__;
    double xiaohaitime__;
  };

  template <typename T> 
    class ThreeDimensionVector
    {
    public:
    ThreeDimensionVector(): ThreeX__(0.), ThreeY__(0.), ThreeZ__(0.){}
    ThreeDimensionVector(const T x, const T y, const T z)
      : ThreeX__(x), ThreeY__(y), ThreeZ__(z){}
      ~ThreeDimensionVector(){};
      ThreeDimensionVector &operator=(ThreeDimensionVector &);
      const T GetX() const;
      const T GetY() const;
      const T GetZ() const;
      const T GetPt() const;
      const T GetMag() const;
      const T GetMag2() const;
      ThreeDimensionVector& operator+=(const ThreeDimensionVector &);
      ThreeDimensionVector operator+(const ThreeDimensionVector &);
      ThreeDimensionVector& operator-=(const ThreeDimensionVector &);
      ThreeDimensionVector operator-(const ThreeDimensionVector &);
      ThreeDimensionVector& operator*=(const T);
      ThreeDimensionVector operator*(const T);
      ThreeDimensionVector& operator-();

    private:
      T ThreeX__;
      T ThreeY__;
      T ThreeZ__;
    };
}

namespace xiaohai
{
  istream &operator>>(istream &is, Particle &particle)
    {
      is >> particle.xiaohaiid__ >> particle.xiaohaipx__ >> particle.xiaohaipy__
         >> particle.xiaohaipz__ >> particle.xiaohaimass__ >> particle.xiaohaix__
         >> particle.xiaohaiy__ >> particle.xiaohaiz__ >> particle.xiaohaitime__;
      return is;
    }
  Particle &Particle::operator=(Particle particle)
    {
      this->xiaohaiid__ = particle.GetId();
      this->xiaohaipx__ = particle.GetPx();
      this->xiaohaipy__ = particle.GetPy();
      this->xiaohaipz__ = particle.GetPz();
      this->xiaohaimass__ = particle.GetMass();
      this->xiaohaix__ = particle.GetX();
      this->xiaohaiy__ = particle.GetY();
      this->xiaohaiz__ = particle.GetZ();
      this->xiaohaitime__ = particle.GetTime();
      return *this;
    }
  
  Particle Particle::Boost(double bx, double by, double bz)
  {
    const unsigned int Dimension = 4;
    double lorentzMatrix[Dimension][Dimension] = {0};
    // lorentz matrix, where GAMMA' = (GAMMA-1)/(beta*beta)
    //  _______________________________________________________________
    //      1+GAMMA'*bx*bx   GAMMA'*bx*by     GAMMA'*bx*bz    GAMMA*bx
    //      GAMMA'*by*bx     1+GAMMA'*by*by   GAMMA'*by*bz    GAMMA*by
    //      GAMMA'*bz*bx     GAMMA'*bz*by     1+GAMMA'*bz*bz  GAMMA*bz
    //      GAMMA*bx         GAMMA*by         GAMMA*bz        GAMMA
    //  _______________________________________________________________
    double beta2 = bx*bx + by*by + bz*bz;
    /// test
    if(beta2 >= 1)
      {
        cout << beta2 << endl;
      }
    double GAMMA = 1.0 / pow(1.0 - beta2, 1./2.);
    double bGAMMA = GAMMA * GAMMA / (1.0 + GAMMA);
    lorentzMatrix[0][0] = 1.0 + bGAMMA * bx * bx;
    lorentzMatrix[1][1] = 1.0 + bGAMMA * by * by;
    lorentzMatrix[2][2] = 1.0 + bGAMMA * bz * bz;
    lorentzMatrix[0][1] = lorentzMatrix[1][0] = bGAMMA * bx * by;
    lorentzMatrix[0][2] = lorentzMatrix[2][0] = bGAMMA * bx * bz;
    lorentzMatrix[1][2] = lorentzMatrix[2][1] = bGAMMA * by * bz;
    lorentzMatrix[0][3] = lorentzMatrix[3][0] = GAMMA * bx;
    lorentzMatrix[1][3] = lorentzMatrix[3][1] = GAMMA * by;
    lorentzMatrix[2][3] = lorentzMatrix[3][2] = GAMMA * bz;
    lorentzMatrix[3][3] = GAMMA;

    double coordinate[Dimension] = {this->xiaohaix__, this->xiaohaiy__,
                                    this->xiaohaiz__, this->xiaohaitime__};
    double tmpcoordinate[Dimension] = {0.};
    for (int i = 0; i != Dimension; ++i)
      {
        for (int k = 0; k != Dimension; ++k)
          {
            tmpcoordinate[i] += lorentzMatrix[i][k] * coordinate[k];
          }//k
      }//i
    this->xiaohaix__ = tmpcoordinate[0];
    this->xiaohaiy__ = tmpcoordinate[1];
    this->xiaohaiz__ = tmpcoordinate[2];
    this->xiaohaitime__ = tmpcoordinate[3];

    double xiaohaienergy__ = this->GetEnergy();
    double momentum[Dimension] = {this->xiaohaipx__, this->xiaohaipy__,
                                  this->xiaohaipz__, xiaohaienergy__};
    double tmpMomentum[Dimension] = {0.};
    for (int i = 0; i != Dimension; ++i)
      {
        for (int k = 0; k != Dimension; ++k)
          {
            tmpMomentum[i] += lorentzMatrix[i][k] * momentum[k];
          }//k
      }//i
    this->xiaohaipx__ = tmpMomentum[0];
    this->xiaohaipy__ = tmpMomentum[1];
    this->xiaohaipz__ = tmpMomentum[2];
    return *this;
  }
  const int Particle::GetId() const
  {
    return this->xiaohaiid__;
  }
  const double Particle::GetPx() const
  {
    return this->xiaohaipx__;
  }
  const double Particle::GetPy() const
  {
    return this->xiaohaipy__;
  }
  const double Particle::GetPz() const
  {
    return this->xiaohaipz__;
  }
  const double Particle::GetMass() const
  {
    return this->xiaohaimass__;
  }
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
  const double Particle::GetTime() const
  {
    return this->xiaohaitime__;
  }
  const double Particle::GetEnergy() const
  {
    double energySquare = this->xiaohaipx__*this->xiaohaipx__
      + this->xiaohaipy__*this->xiaohaipy__
      + this->xiaohaipz__*this->xiaohaipz__
      + this->xiaohaimass__*this->xiaohaimass__;
    return pow(energySquare, 1./2.);
  }
  const double Particle::GetRapidity() const
  {
    double energy = GetEnergy();
    return 1./2.*log((energy+this->xiaohaipz__)/(energy-this->xiaohaipz__));
  }
  const double Particle::GetPt() const
  {
    return sqrt(this->xiaohaipx__*this->xiaohaipx__
                + this->xiaohaipy__*this->xiaohaipy__);
  }
  const double Particle::GetBetaX() const
  {
    return this->xiaohaipx__ / this->GetEnergy();
  }
  const double Particle::GetBetaY() const
  {
    return this->xiaohaipy__ / this->GetEnergy();
  }
  const double Particle::GetBetaZ() const
  {
    return this->xiaohaipz__ / this->GetEnergy();
  }

  /// __________________________________________________
  template <typename T>
    ThreeDimensionVector<T>& ThreeDimensionVector<T>::
    operator=(ThreeDimensionVector<T> &threedimensionvector)
    {
      this->ThreeX__ = threedimensionvector.ThreeX__;
      this->ThreeY__ = threedimensionvector.ThreeY__;
      this->ThreeZ__ = threedimensionvector.ThreeZ__;
      return *this;
    } /// template

  template <typename T>
    const T ThreeDimensionVector<T>::GetX() const
    {
      return this->ThreeX__;
    }

  template <typename T>
    const T ThreeDimensionVector<T>::GetY() const
    {
      return this->ThreeY__;
    } /// template

  template <typename T>
    const T ThreeDimensionVector<T>::GetZ() const
    {
      return this->ThreeZ__;
    } /// template

  template <typename T>
    const T ThreeDimensionVector<T>::GetPt() const
    {
      return sqrt(this->ThreeX__*this->ThreeX__
                  + this->ThreeY__*this->ThreeY__);
    } /// template

  template <typename T>
    const T ThreeDimensionVector<T>::GetMag() const
    {
      return sqrt(this->ThreeX__*this->ThreeX__
                  + this->ThreeY__*this->ThreeY__
                  + this->ThreeZ__*this->ThreeZ__);
    } /// template

  template <typename T>
    const T ThreeDimensionVector<T>::GetMag2() const
    {
      return (this->ThreeX__*this->ThreeX__
              + this->ThreeY__*this->ThreeY__
              + this->ThreeZ__*this->ThreeZ__);
    } /// template

  template <typename T>
    ThreeDimensionVector<T>& ThreeDimensionVector<T>
    ::operator+=(const ThreeDimensionVector<T> &threedimensionvector)
    {
      this->ThreeX__ += threedimensionvector.ThreeX__;
      this->ThreeY__ += threedimensionvector.ThreeY__;
      this->ThreeZ__ += threedimensionvector.ThreeZ__;
      return *this;
    } /// template

  template <typename T>
    ThreeDimensionVector<T> ThreeDimensionVector<T>
    ::operator+(const ThreeDimensionVector<T> &threedimensionvector)
    {
      ThreeDimensionVector<T> tmp;
      tmp.ThreeX__ = this->ThreeX__ + threedimensionvector.ThreeX__;
      tmp.ThreeY__ = this->ThreeY__ + threedimensionvector.ThreeY__;
      tmp.ThreeZ__ = this->ThreeZ__ + threedimensionvector.ThreeZ__;
      return tmp;
    } /// template

  template <typename T>
    ThreeDimensionVector<T>& ThreeDimensionVector<T>
    ::operator-=(const ThreeDimensionVector<T> &threedimensionvector)
    {
      this->ThreeX__ -= threedimensionvector.ThreeX__;
      this->ThreeY__ -= threedimensionvector.ThreeY__;
      this->ThreeZ__ -= threedimensionvector.ThreeZ__;
      return *this;
    } /// template

  template <typename T>
    ThreeDimensionVector<T> ThreeDimensionVector<T>
    ::operator-(const ThreeDimensionVector<T> &threedimensionvector)
    {
      ThreeDimensionVector<T> tmp;
      tmp.ThreeX__ = this->ThreeX__ - threedimensionvector.ThreeX__;
      tmp.ThreeY__ = this->ThreeY__ - threedimensionvector.ThreeY__;
      tmp.ThreeZ__ = this->ThreeZ__ - threedimensionvector.ThreeZ__;
      return tmp;
    } /// template
  
  template <typename T>
    ThreeDimensionVector<T>& ThreeDimensionVector<T>
    ::operator*=(const T num)
    {
      this->ThreeX__ *= num;
      this->ThreeY__ *= num;
      this->ThreeZ__ *= num;
      return *this;
    } /// template

  template <typename T>
    ThreeDimensionVector<T> ThreeDimensionVector<T>
    ::operator*(const T num)
    {
      ThreeDimensionVector<T> tmp;
      tmp.ThreeX__ = this->ThreeX__ * num;
      tmp.ThreeY__ = this->ThreeY__ * num;
      tmp.ThreeZ__ = this->ThreeZ__ * num;
      return tmp;
    } /// template

  template <typename T>
    ThreeDimensionVector<T>& ThreeDimensionVector<T>
    ::operator-()
    {
      this->ThreeX__ = (-this->ThreeX__);
      this->ThreeY__ = (-this->ThreeY__);
      this->ThreeZ__ = (-this->ThreeZ__);
      return *this;
    } /// template
}
#endif /* PARTICLE_H */

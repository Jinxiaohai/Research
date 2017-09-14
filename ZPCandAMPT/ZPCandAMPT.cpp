#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>


#include "ZPCandAMPT.h"

using namespace std;

namespace xiaohai
{
  /**
   * constructor
   *
   * @param is
   */
  ZPCandAMPT::ZPCandAMPT(istream &is)
  {
    is >> *this;
  }

  /**
   * copy constructor
   *
   * @param zpcandampt
   */
  ZPCandAMPT::ZPCandAMPT(const ZPCandAMPT &zpcandampt)
  {
    this->ZPCandAMPTEventNum = zpcandampt.GetZPCandAMPTEventNum();
    this->ZPCandAMPTRunNum = zpcandampt.GetZPCandAMPTRunNum();
    this->ZPCandAMPTMulti = zpcandampt.GetZPCandAMPTMulti();
    this->ZPCandAMPTImpactpar = zpcandampt.GetZPCandAMPTImpactpar();
    this->ZPCandAMPTNelp = zpcandampt.GetZPCandAMPTNelp();
    this->ZPCandAMPTNinp = zpcandampt.GetZPCandAMPTNinp();
    this->ZPCandAMPTNelt = zpcandampt.GetZPCandAMPTNelt();
    this->ZPCandAMPTNint = zpcandampt.GetZPCandAMPTNint();
    /// for track
    this->ZPCandAMPTIdlist = new int[this->ZPCandAMPTMulti];
    for (int j = 0; j != 4; ++j)
      {
	this->ZPCandAMPTPseudoMomentum[j] = new double[this->ZPCandAMPTMulti];
	this->ZPCandAMPTCoordinate[j] = new double[this->ZPCandAMPTMulti];
      }
    for (int i = 0; i != this->ZPCandAMPTMulti; ++i)
      {
	this->ZPCandAMPTIdlist[i] = zpcandampt.ZPCandAMPTIdlist[i];
	for (int j = 0; j != 4; ++j)
	  {
	    this->ZPCandAMPTPseudoMomentum[j][i] = zpcandampt.ZPCandAMPTPseudoMomentum[j][i];
	    this->ZPCandAMPTCoordinate[j][i] = zpcandampt.ZPCandAMPTCoordinate[j][i];
	  }
      }
  }

  ZPCandAMPT::~ZPCandAMPT()
  {
    delete ZPCandAMPTIdlist;
    for (int j = 0; j != 4; ++j)
      {
	delete ZPCandAMPTPseudoMomentum[j];
	delete ZPCandAMPTCoordinate[j];
      }
  }


  /// Get
  inline
  int ZPCandAMPT::GetZPCandAMPTEventNum() const
  {
    return this->ZPCandAMPTEventNum;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTRunNum() const
  {
    return this->ZPCandAMPTRunNum;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTMulti() const
  {
    return this->ZPCandAMPTMulti;
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTImpactpar() const
  {
    return this->ZPCandAMPTImpactpar;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTNelp() const
  {
    return this->ZPCandAMPTNelp;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTNinp() const
  {
    return this->ZPCandAMPTNinp;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTNelt() const
  {
    return this->ZPCandAMPTNelt;
  }

  inline
  int ZPCandAMPT::GetZPCandAMPTNint() const
  {
    return this->ZPCandAMPTNint;
  }

  /// set
  inline
  void ZPCandAMPT::SetZPCandAMPTEventNum(int ZPCandAMPTEventNum)
  {
    this->ZPCandAMPTEventNum = ZPCandAMPTEventNum;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTRunNum(int ZPCandAMPTRunNum)
  {
    this->ZPCandAMPTRunNum = ZPCandAMPTRunNum;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTMulti(int ZPCandAMPTEventMulti)
  {
    this->ZPCandAMPTMulti = ZPCandAMPTMulti;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTImpactpar(int ZPCandAMPTImpactpar){
    this->ZPCandAMPTImpactpar = ZPCandAMPTImpactpar;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTNelp(int ZPCandAMPTNelp)
  {
    this->ZPCandAMPTNelp = ZPCandAMPTNelp;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTNinp(int ZPCandAMPTNinp)
  {
    this->ZPCandAMPTNinp = ZPCandAMPTNinp;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTNelt(int ZPCandAMPTNelt)
  {
    this->ZPCandAMPTNelt = ZPCandAMPTNelt;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTNint(int ZPCandAMPTNint)
  {
    this->ZPCandAMPTNint = ZPCandAMPTNint;
  }

  /// Get
  inline
  int ZPCandAMPT::GetZPCandAMPTId(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTIdlist[number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTPx(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTPseudoMomentum[0][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTPy(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTPseudoMomentum[1][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTPz(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTPseudoMomentum[2][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTMass(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTPseudoMomentum[3][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTX(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTCoordinate[0][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTY(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTCoordinate[1][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTZ(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTCoordinate[2][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTTime(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return this->ZPCandAMPTCoordinate[3][number];
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTEnergySquare(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    double EnergySquare = 0.;
    for (int j = 0; j != 4; ++j)
      {
	EnergySquare += this->ZPCandAMPTPseudoMomentum[j][number];
      }
    return EnergySquare;
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTEnergy(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return sqrt(GetZPCandAMPTEnergySquare(number));
  }

  inline
  double ZPCandAMPT::GetZPCandAMPTPt(int number) const
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    return sqrt(this->ZPCandAMPTPseudoMomentum[0][number]
		*this->ZPCandAMPTPseudoMomentum[0][number]
		+this->ZPCandAMPTPseudoMomentum[1][number]
		*this->ZPCandAMPTPseudoMomentum[1][number]);
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTId(int number, int Id)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTIdlist[number] = Id;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTPx(int number, int Px)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTPseudoMomentum[0][number] = Px;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTPy(int number, int Py)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTPseudoMomentum[1][number] = Py;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTPz(int number, int Pz)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTPseudoMomentum[2][number] = Pz;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTMass(int number, int Mass)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTPseudoMomentum[3][number] = Mass;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTX(int number, int X)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTCoordinate[0][number] = X;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTY(int number, int Y)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTCoordinate[1][number] = Y;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTZ(int number, int Z)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTCoordinate[2][number] = Z;
  }

  inline
  void ZPCandAMPT::SetZPCandAMPTTime(int number, int Time)
  {
    if (number >= this->ZPCandAMPTMulti)
      {
	throw Exception();
      }
    this->ZPCandAMPTCoordinate[3][number] = Time;
  }


  /**
   * operator >>
   *
   * @param is
   * @param zpcandampt
   *
   * @return is
   */
  istream &operator>>(istream &is, ZPCandAMPT &zpcandampt)
  {
    is >> zpcandampt.ZPCandAMPTEventNum
       >> zpcandampt.ZPCandAMPTRunNum
       >> zpcandampt.ZPCandAMPTMulti
       >> zpcandampt.ZPCandAMPTImpactpar
       >> zpcandampt.ZPCandAMPTNelp
       >> zpcandampt.ZPCandAMPTNinp
       >> zpcandampt.ZPCandAMPTNelt
       >> zpcandampt.ZPCandAMPTNint;

    zpcandampt.ZPCandAMPTIdlist = new int[zpcandampt.ZPCandAMPTMulti];
    for (int i = 0; i != 4; ++i)
      {
	zpcandampt.ZPCandAMPTPseudoMomentum[i] = new double[zpcandampt.ZPCandAMPTMulti];
	zpcandampt.ZPCandAMPTCoordinate[i] = new double[zpcandampt.ZPCandAMPTMulti];
      }
    for (int i = 0; i != zpcandampt.ZPCandAMPTMulti; ++i)
      {
	is >> zpcandampt.ZPCandAMPTIdlist[i];
	  for (int j = 0; j != 4; ++j)
	    {
	      is >> zpcandampt.ZPCandAMPTPseudoMomentum[j][i]
		 >> zpcandampt.ZPCandAMPTCoordinate[j][i];
	    }
      }
      return is;
  }

  /**
   * operator <<
   *
   * @param os
   * @param zpcandampt
   *
   * @return os
   */
  ostream &operator<<(ostream &os, const ZPCandAMPT &zpcandampt)
  {
    os << zpcandampt.ZPCandAMPTEventNum << "   "
       << zpcandampt.ZPCandAMPTRunNum << "   "
       << zpcandampt.ZPCandAMPTMulti << "   "
       << zpcandampt.ZPCandAMPTImpactpar << "   "
       << zpcandampt.ZPCandAMPTNelp << "   "
       << zpcandampt.ZPCandAMPTNinp << "   "
       << zpcandampt.ZPCandAMPTNelt << "   "
       << zpcandampt.ZPCandAMPTNint << endl;
    for (int i = 0; i != zpcandampt.ZPCandAMPTMulti; ++i)
      {
	os << zpcandampt.ZPCandAMPTIdlist[i] << "  ";
	for (int j = 0; j != 4; ++j)
	  {
	    os << zpcandampt.ZPCandAMPTPseudoMomentum[j][i] << "   "
	       << zpcandampt.ZPCandAMPTCoordinate[j][i] << "   ";
	  }
	cout << endl;
      }
    return os;
  }

}/// end namespace xiaohai


namespace xiaohai
{
  inline
  string Exception::What() const
  {
    return this->Message;
  }
}

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace xiaohai
{
  class ZPCandAMPT
  {
    friend istream& operator>>(istream &is, ZPCandAMPT &zpcandampt);
    friend ostream& operator<<(ostream &os, const ZPCandAMPT &zpcandampt);
  public:
    /// constructor and copy constructor
    ZPCandAMPT() = default;
    ZPCandAMPT(istream &is);
    ZPCandAMPT(const ZPCandAMPT &zpcandampt);
    ~ZPCandAMPT();

    /// Get
    virtual int GetZPCandAMPTEventNum() const;
    virtual int GetZPCandAMPTRunNum() const;
    virtual int GetZPCandAMPTMulti() const;
    virtual double GetZPCandAMPTImpactpar() const;
    virtual int GetZPCandAMPTNelp() const;
    virtual int GetZPCandAMPTNinp() const;
    virtual int GetZPCandAMPTNelt() const;
    virtual int GetZPCandAMPTNint() const;
    /// Set
    virtual void SetZPCandAMPTEventNum(int ZPCandAMPTEventNum);
    virtual void SetZPCandAMPTRunNum(int ZPCandAMPTRunNum);
    virtual void SetZPCandAMPTMulti(int ZPCandAMPTEventMulti);
    virtual void SetZPCandAMPTImpactpar(int ZPCandAMPTImpactpar);
    virtual void SetZPCandAMPTNelp(int ZPCandAMPTNelp);
    virtual void SetZPCandAMPTNinp(int ZPCandAMPTNinp);
    virtual void SetZPCandAMPTNelt(int ZPCandAMPTNelt);
    virtual void SetZPCandAMPTNint(int ZPCandAMPTNint);
    /// Get
    virtual int GetZPCandAMPTId(int number) const;
    virtual double GetZPCandAMPTPx(int number) const;
    virtual double GetZPCandAMPTPy(int number) const;
    virtual double GetZPCandAMPTPz(int number) const;
    virtual double GetZPCandAMPTMass(int number) const;
    virtual double GetZPCandAMPTX(int number) const;
    virtual double GetZPCandAMPTY(int number) const;
    virtual double GetZPCandAMPTZ(int number) const;
    virtual double GetZPCandAMPTTime(int number) const;
    virtual double GetZPCandAMPTEnergySquare(int number) const;
    virtual double GetZPCandAMPTEnergy(int number) const;
    virtual double GetZPCandAMPTPt(int number) const;
    /// Set
    virtual void SetZPCandAMPTId(int number, int Id);
    virtual void SetZPCandAMPTPx(int number, int Px);
    virtual void SetZPCandAMPTPy(int number, int Py);
    virtual void SetZPCandAMPTPz(int number, int Pz);
    virtual void SetZPCandAMPTMass(int number, int Mass);
    virtual void SetZPCandAMPTX(int number, int X);
    virtual void SetZPCandAMPTY(int number, int Y);
    virtual void SetZPCandAMPTZ(int number, int Z);
    virtual void SetZPCandAMPTTime(int number, int Time);


  private:
    int ZPCandAMPTEventNum, ZPCandAMPTRunNum, ZPCandAMPTMulti;
    int ZPCandAMPTNelp, ZPCandAMPTNinp;
    int ZPCandAMPTNelt, ZPCandAMPTNint;
    double ZPCandAMPTImpactpar;
    int *ZPCandAMPTIdlist;
    double *ZPCandAMPTPseudoMomentum[4];
    double *ZPCandAMPTCoordinate[4];
  };

  istream &operator>>(istream &is, ZPCandAMPT &zpcandampt);
  ostream &operator<<(ostream &os, const ZPCandAMPT &zpcandampt);

}


namespace xiaohai
{
  class Exception
  {
  public:
  Exception(string message="some errors happened.") : Message(message) {}
    ~Exception(){}
    virtual string What() const;

  private:
    string Message;
  };
}

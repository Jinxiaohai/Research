
/**
 * @file   Particle.cpp
 * @author xiaohai <xiaohaijin@outlook.com>
 * @date   Tue Jun  6 20:29:13 2017
 * 
 * @brief  The C++ file was created to read ampt.dat
 *         file.
 */
#include <iostream>
#include <cmath>
#include <string>
using namespace std;


namespace xiaohai
{
  /** 
   * deal with the basic the errors.
   * Created by xiaohai.
   */
  class Exception
  {
  public:
    Exception(const string theMessage="Some errors have happened.")
      : message(theMessage){ }
    virtual string What() const
    {
      return message;
    }
    virtual ~Exception(){ }

  private:
    string message;
  };
}



namespace xiaohai
{
  /** 
   * deal with the track infomation.
   * Created by xiaohai.
   */
  template <typename T>
  class Track
  {
    /// friend
    template <typename U>
    friend istream& operator>>(istream &is, Track<U> &track);
    template <typename U>
    friend ostream& operator<<(ostream &os, Track<U> &track);
  public:
    /// constructor
    Track()
      :fId(0), fPx(0), fPy(0), fPz(0), fMass(0), fX(0), fY(0), fZ(0), fTime(0){ }
    Track(int id, T px, T py, T pz, T mass, T x, T y, T z, T ftime)
      : fId(id), fPx(px), fPy(py), fPz(pz), fMass(mass), fX(x), fY(y), fZ(z), fTime(ftime)
    { }
    Track(const Track<T>& track);
    virtual ~Track() { }
    Track& operator=(const Track<T>& track);
    /// get functions
    virtual int GetId() const;
    virtual T GetPx() const;
    virtual T GetPy() const;
    virtual T GetPz() const;
    virtual T GetMass() const;
    virtual T GetX() const;
    virtual T GetY() const;
    virtual T GetZ() const;
    virtual T GetTime() const;
    /// return sqrt(fx*fx+fy*fy)
    virtual T GetR() const;
    /// return momentum
    virtual T GetP() const;
    /// return Pt
    virtual T GetPt() const;
    /// return Energy
    virtual T GetEnergy() const;
    /// return rapidity
    virtual T GetRapidity() const;

    /// set functions
    virtual void SetId(int id);
    virtual void SetPx(T px);
    virtual void SetPy(T py);
    virtual void SetPz(T pz);
    virtual void SetMass(T mass);
    virtual void SetX(T x);
    virtual void SetY(T y);
    virtual void SetZ(T z);
    virtual void SetTime(T ftime);
    
  private:
    int fId;
    T fPx;
    T fPy;
    T fPz;
    T fMass;
    T fX;
    T fY;
    T fZ;
    T fTime;
  };

  template <typename T>
  Track<T>::Track(const Track<T>& track)
  {
    this->fId = track.fId;
    this->fPx = track.fPx;
    this->fPy = track.fPy;
    this->fPz = track.fPz;
    this->fMass = track.fMass;
    this->fX = track.fX;
    this->fY = track.fY;
    this->fZ = track.fZ;
    this->fTime = track.fTime;
  } ///< template

  template <typename T>
  Track<T>& Track<T>::operator=(const Track<T> &track)
  {
    this->fId = track.fId;
    this->fPx = track.fPx;
    this->fPy = track.fPy;
    this->fPz = track.fPz;
    this->fMass = track.fMass;
    this->fX = track.fX;
    this->fY = track.fY;
    this->fZ = track.fZ;
    this->fTime = track.fTime;
    return *this;
  } ///< template

  template <typename U>
  istream& operator>>(istream& is, Track<U> &track)
  {
    is >> track.fId >> track.fPx >> track.fPy >> track.fPz
       >> track.fMass >> track.fX >> track.fY >> track.fZ
       >> track.fTime;
    return is;
  } ///< template

  template <typename T>
  ostream& operator<<(ostream& os, Track<T> &track)
  {
    os << track.fId << "  " << track.fPx << "  " << track.fPy << "  "
       << track.fPz << "  " << track.fMass << "  " << track.fX << "  "
       << track.fY << "  " << track.fZ << "  " << track.fTime << "\n";
    return os;
  } ///< template

  template <typename T>
  inline int Track<T>::GetId() const
  {
    return this->fId;
  } ///< template

  template <typename T>
  inline T Track<T>::GetPx() const
  {
    return this->fPx;
  } ///< template

  template <typename T>
  inline T Track<T>::GetPy() const
  {
    return this->fPy;
  } ///< template

  template <typename T>
  inline T Track<T>::GetPz() const
  {
    return this->fPz;
  } ///< template

  template <typename T>
  inline T Track<T>::GetMass() const
  {
    return this->fMass;
  } ///< template

  template <typename T>
  inline T Track<T>::GetX() const
  {
    return this->fX;
  } ///< template

  template <typename T>
  inline T Track<T>::GetY() const
  {
    return this->fY;
  } ///< template

  template <typename T>
  inline T Track<T>::GetZ() const
  {
    return this->fZ;
  } ///< template

  template <typename T>
  inline T Track<T>::GetTime() const
  {
    return this->fTime;
  } ///< template

  template <typename T>
  inline T Track<T>::GetR() const
  {
    return sqrt(this->fX*this->fX + this->fY*this->fY);
  } ///< template

  template <typename T>
  inline T Track<T>::GetP() const
  {
    return sqrt(this->fPx*this->fPx + this->fPy*this->fPy + this->fPz*this->fPz);
  } ///< template

  template <typename T>
  inline T Track<T>::GetPt() const
  {
    return sqrt(this->fPx*this->fPx + this->fPy*this->fPy);
  } ///< template

  template <typename T>
  inline T Track<T>::GetEnergy() const
  {
    return sqrt(this->fPx*this->fPx + this->fPy*this->fPy
		+this->fPz*this->fPz + this->fMass*this->fMass);
  } ///< template

  template <typename T>
  inline T Track<T>::GetRapidity() const
  {
    return 0.5*log( (GetEnergy()+fPz) / (GetEnergy()-fPz));
  } ///< template

  template <typename T>
  inline void Track<T>::SetId(int id)
  {
    this->fId = id;
  } ///< template

  template <typename T>
  inline void Track<T>::SetPx(T px)
  {
    this->fPx = px;
  } ///< template

  template <typename T>
  inline void Track<T>::SetPy(T py)
  {
    this->fPy = py;
  } ///< template

  template <typename T>
  inline void Track<T>::SetPz(T pz)
  {
    this->fPz = pz;
  } ///< template

  template <typename T>
  inline void Track<T>::SetMass(T mass)
  {
    this->fMass = mass;
  } ///< template

  template <typename T>
  inline void Track<T>::SetX(T x)
  {
    this->fX = x;
  } ///< template

  template <typename T>
  inline void Track<T>::SetY(T y)
  {
    this->fY = y;
  } ///< template

  template <typename T>
  inline void Track<T>::SetZ(T z)
  {
    this->fZ = z;
  } ///< template

  template <typename T>
  inline void Track<T>::SetTime(T ftime)
  {
    this->fTime = ftime;
  } ///< template


  /** 
   * deal with the event infomation.
   * Created by xiaohai.
   */
  template <typename T>
  class Event
  {
    template <typename U>
    friend istream& operator>>(istream& is, Event<U>& event);
    template <typename U>
    friend ostream& operator<<(ostream& os, Event<U>& event);
  public:
    Event():fEventNum(0), fRunNum(0), fMulti(0), fImpactpar(0), fNpartP(0),
	    fNpartT(0), fInelasticP(0), fElasticP(0), fInelasticT(0),
	    fElasticT(0), fPasshead(0), vectorTrack(NULL){ }
    Event(int EVENTNUM, int RUNNUM, int MULTI, T IMPACTPAR, int NPARTP,
	  int NPARTT, int INELASTICP, int ELASTICP, int INELASTICT,
	  int ELASTICT, T PASSHEAD, Track<T>* VECTORTRACK)
      :fEventNum(EVENTNUM), fRunNum(RUNNUM), fMulti(MULTI), fImpactpar(IMPACTPAR),
       fNpartP(NPARTP), fNpartT(NPARTT), fInelasticP(INELASTICP), fElasticP(ELASTICP),
       fInelasticT(INELASTICT), fElasticT(ELASTICT), fPasshead(PASSHEAD),
       vectorTrack(VECTORTRACK){ }
    /// copy constructor
    Event(const Event<T>& event);
    virtual ~Event()
    {
      delete [] vectorTrack;
    }
    /// Print the data.
    virtual void Print() const;
    /// Get
    virtual int GetEventNum() const;
    virtual int GetRunNum() const;
    virtual int GetMulti() const;
    virtual T GetImpactpar() const;
    virtual int GetNpartP() const;
    virtual int GetNpartT() const;
    virtual int GetInelasticP() const;
    virtual int GetElasticP() const;
    virtual int GetInelasticT() const;
    virtual int GetElasticT() const;
    virtual T GetPasshead() const;

    /// Set
    virtual void SetEventNum(int fEventNum);
    virtual void SetRunNum(int fRunNum);
    virtual void SetMulti(int fMulti);
    virtual void SetImpactpar(T fImpactpart);
    virtual void SetNpartP(int fNpartP);
    virtual void SetNpartT(int fNpartT);
    virtual void SetInelasticP(int fInelasticP);
    virtual void SetElasticP(int fElasticP);
    virtual void SetInelasticT(int fInelasticT);
    virtual void SetElasticT(int fElasticT);
    virtual void SetPasshead(T fPasshead);

    /// Get track i
    virtual Track<T> &GetTrack(int iTrack) const;
    virtual Event<T>& operator=(const Event<T>& event);
  
  private:
    /// event head
    int fEventNum;
    int fRunNum;
    int fMulti;
    T fImpactpar;
    int fNpartP;
    int fNpartT;
    int fInelasticP;
    int fElasticP;
    int fInelasticT;
    int fElasticT;
    T fPasshead;
    Track<T> *vectorTrack;
  }; /// template class

  template <class U>
  istream& operator>>(istream& is, Event<U>& event)
  {
    is >> event.fEventNum >> event.fRunNum >> event.fMulti >> event.fImpactpar
       >> event.fNpartP >> event.fNpartT >> event.fInelasticP >> event.fElasticP
       >> event.fInelasticT >> event.fElasticT >> event.fPasshead;
    
    event.vectorTrack = new Track<U>[event.fMulti];
    for (int i0 = 0; i0 != event.fMulti; ++i0)
      {
	is >> event.vectorTrack[i0];
      }
    return is;
  } ///< template

  template <typename U>
  ostream& operator<<(ostream& os, Event<U>& event)
  {
    os << event.fEventNum << "  " << event.fRunNum << "  " << event.fMulti << "  "
       << event.fImpactpar << "  " << event.fNpartP << "  " << event.fNpartT << "  "
       << event.fInelasticP << "  " << event.fElasticP << "  " << event.fInelasticT
       << "  " << event.fElasticT << "  " << event.fPasshead << "\n";
    for (int i0 = 0; i0 != event.fMulti; ++i0)
      {
	os << event.vectorTrack[i0];
      }
    return os;
  } ///< template

  template <typename T>
  inline void Event<T>::Print() const
  {
    cout << this->fEventNum << "  " << this->fRunNum << "  " << this->fMulti << "  "
	 << this->fImpactpar << "  " << this->fNpartP << "  " << this->fNpartT << "  "
	 << this->fInelasticP << "  " << this->fElasticP << "  " << this->fInelasticT
	 << "  " << this->fElasticT << "  " << this->fPasshead << "\n";
    for (int i0 = 0; i0 != this->fMulti; ++i0)
      {
    	cout << this->vectorTrack[i0];
      }
  } ///< template

  /// copy constructor
  template <typename T>
  Event<T>::Event(const Event<T> &event)
  {
    this->fEventNum = event.fEventNum;
    this->fRunNum = event.fRunNum;
    this->fMulti = event.fMulti;
    this->fImpactpar = event.fImpactpar;
    this->fNpartP = event.fNpartP;
    this->fNpartT = event.fNpartT;
    this->fInelasticP = event.fInelasticP;
    this->fElasticP = event.fElasticP;
    this->fInelasticT = event.fInelasticT;
    this->fElasticT = event.fElasticT;
    this->fPasshead = event.fPasshead;
    this->vectorTrack = new Track<T>[event.fMulti];
    for (int i0 = 0; i0 != event.fMulti; ++i0)
      {
	this->vectorTrack[i0] = event.vectorTrack[i0];
      }
  }

  /// operator=
  template <typename T>
  Event<T>& Event<T>::operator=(const Event<T>& event)
  {
    this->fEventNum = event.fEventNum;
    this->fRunNum = event.fRunNum;
    this->fMulti = event.fMulti;
    this->fImpactpar = event.fImpactpar;
    this->fNpartP = event.fNpartP;
    this->fNpartT = event.fNpartT;
    this->fInelasticP = event.fInelasticP;
    this->fElasticP = event.fElasticP;
    this->fInelasticT = event.fInelasticT;
    this->fElasticT = event.fElasticT;
    this->fPasshead = event.fPasshead;
    this->vectorTrack = new Track<T>[event.fMulti];
    for (int i0 = 0; i0 != event.fMulti; ++i0)
      {
	this->vectorTrack[i0] = event.vectorTrack[i0];
      }
    return *this;
  } ///< template

  template <typename T>
  inline int Event<T>::GetEventNum() const
  {
    return this->fEventNum;
  } ///< template
  
  template <typename T>
  inline int Event<T>::GetRunNum() const
  {
    return this->fRunNum;
  } ///< template

  template <typename T>
  inline int Event<T>::GetMulti() const
  {
    return this->fMulti;
  } ///< template

  template <typename T>
  inline T Event<T>::GetImpactpar() const
  {
    return this->fImpactpar;
  } ///< template
  
  template <typename T>
  inline int Event<T>::GetNpartP() const
  {
    return this->fNpartP;
  } ///< template

  template <typename T>
  inline int Event<T>::GetNpartT() const
  {
    return this->fNpartT;
  } ///< template

  template <typename T>
  inline int Event<T>::GetInelasticP() const
  {
    return this->fInelasticP;
  } ///< template

  template <typename T>
  inline int Event<T>::GetElasticP() const
  {
    return this->fElasticP;
  } ///< template

  template <typename T>
  inline int Event<T>::GetInelasticT() const
  {
    return this->fInelasticT;
  } ///< template

  template <typename T>
  inline int Event<T>::GetElasticT() const
  {
    return this->fElasticT;
  } ///< template

  template <typename T>
  inline T Event<T>::GetPasshead() const
  {
    return this->fPasshead;
  } ///< template

  template <typename T>
  inline void Event<T>::SetEventNum(int fEventNum)
  {
    this->fEventNum = fEventNum;
  } ///< template

  template <typename T>
  inline void Event<T>::SetRunNum(int fRunNum)
  {
    this->fRunNum = fRunNum;
  } ///< template

  template <typename T>
  inline void Event<T>::SetMulti(int fMulti)
  {
    this->fMulti = fMulti;
  } ///< template

  template <typename T>
  inline void Event<T>::SetImpactpar(T fImpactpar)
  {
    this->fImpactpar = fImpactpar;
  } ///< template

  template <typename T>
  inline void Event<T>::SetNpartP(int fNpartP)
  {
    this->fNpartP = fNpartP;
  } ///< template

  template <typename T>
  inline void Event<T>::SetNpartT(int fNpartT)
  {
    this->fNpartT = fNpartT;
  } ///< template

  template <typename T>
  inline void Event<T>::SetInelasticP(int fInelasticP)
  {
    this->fInelasticP = fInelasticP;
  } ///< template

  template <typename T>
  inline void Event<T>::SetElasticP(int fElasticP)
  {
    this->fElasticP = fElasticP;
  } ///< template

  template <typename T>
  inline void Event<T>::SetInelasticT(int fInelasticT)
  {
    this->fInelasticT = fInelasticT;
  } ///< template

  template <typename T>
  inline void Event<T>::SetElasticT(int fElasticT)
  {
    this->fElasticT = fElasticT;
  } ///< template

  template <typename T>
  inline void Event<T>::SetPasshead(T fPasshead)
  {
    this->fPasshead = fPasshead;
  } ///< template
  
  /// Get track i
  template <typename T>
  inline Track<T> &Event<T>::GetTrack(int itrack) const
  {
    if(itrack < 0)
      {
	throw Exception("track number must be itrack >= 0 && itrack < multi.");
      }
    else
      {
	return this->vectorTrack[itrack];
      }
  }
  
}

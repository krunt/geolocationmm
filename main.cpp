

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <queue>
#include <mutex>
#include <memory>
#include <stdint.h>
#include <fstream>
#include <set>
#include <map>
#include <memory>
#include <time.h>
#include <future>
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <sys/stat.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#ifndef _NR3_H_
#define _NR3_H_

//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;

// macro-like inline functions

struct MyVector3d
{
  double p[3];
  MyVector3d(void) {}
  MyVector3d(double v1, double v2, double v3) {
    p[0] = v1; p[1] = v2; p[2] = v3;
  }
  MyVector3d operator-(void) const {
    MyVector3d res(-p[0], -p[1], -p[2]);
    return res;
  }
};

inline double MyVectorNormalize(MyVector3d &v)
{
  double n = sqrt(v.p[0] * v.p[0] + v.p[1] * v.p[1] + v.p[2] * v.p[2]);
  if (n != 0) {
    v.p[0] /= n;
    v.p[1] /= n;
    v.p[2] /= n;
  }
  return n;
}

template<class T>
inline T SQR(const T a) {return a*a;}

#if defined(_MSC_VER) && (_MSC_VER < 1300)
// special case of MSVC++ 6 and older compilers

double abs(double x) { return fabs(x); }

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

#else
// general case of ANSI/ISO compliant compilers

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

#endif /* _MSC_VER < 1300 */

inline float MAX(const double &a, const float &b)
{return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
{return b > a ? float(b) : (a);}

inline float MIN(const double &a, const float &b)
{return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
{return b < a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
{T dum=a; a=b; b=dum;}

// exception handling

#ifndef _USENRERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
  char *message;
  char *file;
  int line;
  NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
  printf("ERROR: %s\n     in file %s at line %d\n",
    err.message, err.file, err.line);
  exit(1);
}
#endif

// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)


// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define NRvector vector
#else

template <class T>
class NRvector {
private:
  int nn;	// size of array. upper index is nn-1
  T *v;
public:
  NRvector();
  explicit NRvector(int n);		// Zero-based array
  NRvector(int n, const T &a);	//initialize to constant value
  NRvector(int n, const T *a);	// Initialize to array
  NRvector(const NRvector &rhs);	// Copy constructor
  NRvector & operator=(const NRvector &rhs);	//assignment
  typedef T value_type; // make T available externally
  inline T & operator[](const int i);	//i'th element
  inline const T & operator[](const int i) const;
  inline int size() const;
  void resize(int newn); // resize (contents not preserved)
  void assign(int newn, const T &a); // resize and assign a constant value
  ~NRvector();
};

// NRvector definitions

template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
  for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
  for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
  for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
  // postcondition: normal assignment via copying has been performed;
  //		if vector and rhs were different sizes, vector
  //		has been resized to match the size of rhs
{
  if (this != &rhs)
  {
    if (nn != rhs.nn) {
      if (v != NULL) delete [] (v);
      nn=rhs.nn;
      v= nn>0 ? new T[nn] : NULL;
    }
    for (int i=0; i<nn; i++)
      v[i]=rhs[i];
  }
  return *this;
}

template <class T>
inline T & NRvector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("NRvector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("NRvector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline int NRvector<T>::size() const
{
  return nn;
}

template <class T>
void NRvector<T>::resize(int newn)
{
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
}

template <class T>
void NRvector<T>::assign(int newn, const T& a)
{
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector()
{
  if (v != NULL) delete[] (v);
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_

template <class T>
class NRmatrix {
private:
  int nn;
  int mm;
  T **v;
public:
  NRmatrix();
  NRmatrix(int n, int m);			// Zero-based array
  NRmatrix(int n, int m, const T &a);	//Initialize to constant
  NRmatrix(int n, int m, const T *a);	// Initialize to array
  NRmatrix(const NRmatrix &rhs);		// Copy constructor
  NRmatrix & operator=(const NRmatrix &rhs);	//assignment
  typedef T value_type; // make T available externally
  inline T* operator[](const int i);	//subscripting: pointer to row i
  inline const T* operator[](const int i) const;
  inline int nrows() const;
  inline int ncols() const;
  void resize(int newn, int newm); // resize (contents not preserved)
  void assign(int newn, int newm, const T &a); // resize and assign a constant value
  ~NRmatrix();
};

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
  int i,j,nel=mm*nn;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
  // postcondition: normal assignment via copying has been performed;
  //		if matrix and rhs were different sizes, matrix
  //		has been resized to match the size of rhs
{
  if (this != &rhs) {
    int i,j,nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new T*[nn] : NULL;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : NULL;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("NRmatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("NRmatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const
{
  return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const
{
  return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm)
{
  int i,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T& a)
{
  int i,j,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}

template <class T>
class NRMat3d {
private:
  int nn;
  int mm;
  int kk;
  T ***v;
public:
  NRMat3d();
  NRMat3d(int n, int m, int k);
  inline T** operator[](const int i);	//subscripting: pointer to row i
  inline const T* const * operator[](const int i) const;
  inline int dim1() const;
  inline int dim2() const;
  inline int dim3() const;
  ~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
  int i,j;
  v[0] = new T*[n*m];
  v[0][0] = new T[n*m*k];
  for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
  for(i=1; i<n; i++) {
    v[i] = v[i-1] + m;
    v[i][0] = v[i-1][0] + m*k;
    for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
  }
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
  return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
  return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
  return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
  return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
  return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
  if (v != NULL) {
    delete[] (v[0][0]);
    delete[] (v[0]);
    delete[] (v);
  }
}


// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
  turn_on_floating_exceptions() {
    int cw = _controlfp( 0, 0 );
    cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
    _controlfp( cw, MCW_EM );
  }
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

#endif /* _NR3_H_ */

struct WrapVecDoub {
  VecDoub vvec;
  VecDoub &v;
  Int n, mask;

  WrapVecDoub(const Int nn) : vvec(nn), v(vvec), n(nn/2),
    mask(n-1) {validate();}

  WrapVecDoub(VecDoub &vec) : v(vec), n(vec.size()/2),
    mask(n-1) {validate();}

  void validate() {if (n&(n-1)) throw("vec size must be power of 2");}

  inline Complex& operator[] (Int i) {return (Complex &)v[(i&mask) << 1];}

  inline Doub& real(Int i) {return v[(i&mask) << 1];}

  inline Doub& imag(Int i) {return v[((i&mask) << 1)+1];}

  operator VecDoub&() {return v;}

};

unsigned int clp2(unsigned int x)
{
  x = x - 1;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  return x + 1;
}

void four1(Doub *data, const Int n, const Int isign) {
  Int nn,mmax,m,j,istep,i;
  Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
  if (n<2 || n&(n-1)) throw("n must be power of 2 in four1");
  nn = n << 1;
  j = 1;
  for (i=1;i<nn;i+=2) {
    if (j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    m=n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (nn > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=nn;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j-1]-wi*data[j];
        tempi=wr*data[j]+wi*data[j-1];
        data[j-1]=data[i-1]-tempr;
        data[j]=data[i]-tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

struct SVD {
  Int m,n;
  MatDoub u,v;
  VecDoub w;
  Doub eps, tsh;
  SVD(MatDoub_I &a) : m(a.nrows()), n(a.ncols()), u(a), v(n,n), w(n) {
    eps = numeric_limits<Doub>::epsilon();
    decompose();
    reorder();
    tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;
  }

  void solve(VecDoub_I &b, VecDoub_O &x, Doub thresh);
  void solve(MatDoub_I &b, MatDoub_O &x, Doub thresh);

  Int rank(Doub thresh);
  Int nullity(Doub thresh);
  MatDoub range(Doub thresh);
  MatDoub nullspace(Doub thresh);

  Doub inv_condition() {
    return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
  }

  void decompose();
  void reorder();
  Doub pythag(const Doub a, const Doub b);
};
void SVD::solve(VecDoub_I &b, VecDoub_O &x, Doub thresh = -1.) {
  Int i,j,jj;
  Doub s;
  if (b.size() != m || x.size() != n) throw("SVD::solve bad sizes");
  VecDoub tmp(n);
  tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
  for (j=0;j<n;j++) {
    s=0.0;
    if (w[j] > tsh) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
}

void SVD::solve(MatDoub_I &b, MatDoub_O &x, Doub thresh = -1.)
{
  int i,j,p=b.ncols();
  if (b.nrows() != m || x.nrows() != n || x.ncols() != p)
    throw("SVD::solve bad sizes");
  VecDoub xx(n),bcol(m);
  for (j=0;j<p;j++) {
    for (i=0;i<m;i++) bcol[i] = b[i][j];
    solve(bcol,xx,thresh);
    for (i=0;i<n;i++) x[i][j] = xx[i];
  }
}
Int SVD::rank(Doub thresh = -1.) {
  Int j,nr=0;
  tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
  for (j=0;j<n;j++) if (w[j] > tsh) nr++;
  return nr;
}

Int SVD::nullity(Doub thresh = -1.) {
  Int j,nn=0;
  tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
  for (j=0;j<n;j++) if (w[j] <= tsh) nn++;
  return nn;
}

MatDoub SVD::range(Doub thresh = -1.){
  Int i,j,nr=0;
  MatDoub rnge(m,rank(thresh));
  for (j=0;j<n;j++) {
    if (w[j] > tsh) {
      for (i=0;i<m;i++) rnge[i][nr] = u[i][j];
      nr++;
    }
  }
  return rnge;
}

MatDoub SVD::nullspace(Doub thresh = -1.){
  Int j,jj,nn=0;
  MatDoub nullsp(n,nullity(thresh));
  for (j=0;j<n;j++) {
    if (w[j] <= tsh) {
      for (jj=0;jj<n;jj++) nullsp[jj][nn] = v[jj][j];
      nn++;
    }
  }
  return nullsp;
}
void SVD::decompose() {
  bool flag;
  Int i,its,j,jj,k,l,nm;
  Doub anorm,c,f,g,h,s,scale,x,y,z;
  VecDoub rv1(n);
  g = scale = anorm = 0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += abs(u[k][i]);
      if (scale != 0.0) {
        for (k=i;k<m;k++) {
          u[k][i] /= scale;
          s += u[k][i]*u[k][i];
        }
        f=u[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        u[i][i]=f-g;
        for (j=l-1;j<n;j++) {
          for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
          f=s/h;
          for (k=i;k<m;k++) u[k][j] += f*u[k][i];
        }
        for (k=i;k<m;k++) u[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += abs(u[i][k]);
      if (scale != 0.0) {
        for (k=l-1;k<n;k++) {
          u[i][k] /= scale;
          s += u[i][k]*u[i][k];
        }
        f=u[i][l-1];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        u[i][l-1]=f-g;
        for (k=l-1;k<n;k++) rv1[k]=u[i][k]/h;
        for (j=l-1;j<m;j++) {
          for (s=0.0,k=l-1;k<n;k++) s += u[j][k]*u[i][k];
          for (k=l-1;k<n;k++) u[j][k] += s*rv1[k];
        }
        for (k=l-1;k<n;k++) u[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(abs(w[i])+abs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j<n;j++)
          v[j][i]=(u[i][j]/u[i][l])/g;
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
          for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) u[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
        for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
        f=(s/u[i][i])*g;
        for (k=i;k<m;k++) u[k][j] += f*u[k][i];
      }
      for (j=i;j<m;j++) u[j][i] *= g;
    } else for (j=i;j<m;j++) u[j][i]=0.0;
    ++u[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
        nm=l-1;
        if (l == 0 || abs(rv1[l]) <= eps*anorm) {
          flag=false;
          break;
        }
        if (abs(w[nm]) <= eps*anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<k+1;i++) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if (abs(f) <= eps*anorm) break;
          g=w[i];
          h=pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=0;j<m;j++) {
            y=u[j][nm];
            z=u[j][i];
            u[j][nm]=y*c+z*s;
            u[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=0;j<n;j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 29) throw("no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=0;jj<n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=0;jj<m;jj++) {
          y=u[jj][j];
          z=u[jj][i];
          u[jj][j]=y*c+z*s;
          u[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
}

void SVD::reorder() {
  Int i,j,k,s,inc=1;
  Doub sw;
  VecDoub su(m), sv(n);
  do { inc *= 3; inc++; } while (inc <= n);
  do {
    inc /= 3;
    for (i=inc;i<n;i++) {
      sw = w[i];
      for (k=0;k<m;k++) su[k] = u[k][i];
      for (k=0;k<n;k++) sv[k] = v[k][i];
      j = i;
      while (w[j-inc] < sw) {
        w[j] = w[j-inc];
        for (k=0;k<m;k++) u[k][j] = u[k][j-inc];
        for (k=0;k<n;k++) v[k][j] = v[k][j-inc];
        j -= inc;
        if (j < inc) break;
      }
      w[j] = sw;
      for (k=0;k<m;k++) u[k][j] = su[k];
      for (k=0;k<n;k++) v[k][j] = sv[k];

    }
  } while (inc > 1);
  for (k=0;k<n;k++) {
    s=0;
    for (i=0;i<m;i++) if (u[i][k] < 0.) s++;
    for (j=0;j<n;j++) if (v[j][k] < 0.) s++;
    if (s > (m+n)/2) {
      for (i=0;i<m;i++) u[i][k] = -u[i][k];
      for (j=0;j<n;j++) v[j][k] = -v[j][k];
    }
  }
}

Doub SVD::pythag(const Doub a, const Doub b) {
  Doub absa=abs(a), absb=abs(b);
  return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
    (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

struct Fitsvd {
  Int ndat, ma;
  Doub tol;
  VecDoub_I *x,&y,&sig;
  VecDoub (*funcs)(const Doub);
  VecDoub a;
  MatDoub covar;
  Doub chisq;

  Fitsvd(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig,
    VecDoub funks(const Doub), const Doub TOL=1.e-12)
    : ndat(yy.size()), x(&xx), xmd(NULL), y(yy), sig(ssig),
    funcs(funks), tol(TOL) {}

  void fit() {
    Int i,j,k;
    Doub tmp,thresh,sum;
    if (x) ma = funcs((*x)[0]).size();
    else ma = funcsmd(row(*xmd,0)).size();
    a.resize(ma);
    covar.resize(ma,ma);
    MatDoub aa(ndat,ma);
    VecDoub b(ndat),afunc(ma);
    for (i=0;i<ndat;i++) {
      if (x) afunc=funcs((*x)[i]);
      else afunc=funcsmd(row(*xmd,i));
      tmp=1.0/sig[i];
      for (j=0;j<ma;j++) aa[i][j]=afunc[j]*tmp;
      b[i]=y[i]*tmp;
    }
    SVD svd(aa);
    thresh = (tol > 0. ? tol*svd.w[0] : -1.);
    svd.solve(b,a,thresh);
    chisq=0.0;
    for (i=0;i<ndat;i++) {
      sum=0.;
      for (j=0;j<ma;j++) sum += aa[i][j]*a[j];
      chisq += SQR(sum-b[i]);
    }
    for (i=0;i<ma;i++) {
      for (j=0;j<i+1;j++) {
        sum=0.0;
        for (k=0;k<ma;k++) if (svd.w[k] > svd.tsh)
          sum += svd.v[i][k]*svd.v[j][k]/SQR(svd.w[k]);
        covar[j][i]=covar[i][j]=sum;
      }
    }

  }

  MatDoub_I *xmd;
  VecDoub (*funcsmd)(VecDoub_I &);

  Fitsvd(MatDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig,
    VecDoub funks(VecDoub_I &), const Doub TOL=1.e-12)
    : ndat(yy.size()), x(NULL), xmd(&xx), y(yy), sig(ssig),
    funcsmd(funks), tol(TOL) {}

  VecDoub row(MatDoub_I &a, const Int i) {
    Int j,n=a.ncols();
    VecDoub ans(n);
    for (j=0;j<n;j++) ans[j] = a[i][j];
    return ans;
  }
};


using namespace std;
typedef pair<double, double> ppd;

const double IQ_SIZE_SEC = 5.0;
const int NSAMPLES_PER_SEC = 16000;
const double KNOT2MS = 1.852 * 0.2777;
const double DEG2RAD = M_PI / 180.0;
const double SPEED_OF_LIGHT = 299792458; // m/s
const double F0_GRID_STEPS = 50;
const double SEMI_MAJOR_AXIS = 6378137.0000000000;
const double SEMI_MINOR_AXIS = 6356752.3142451793;
const double MERIDIAN_OFFSET = 0.0;

void geodeticToCartesian(double lon, double lat, double h, 
                           double &x, double &y, double &z)
{
  double a  = SEMI_MAJOR_AXIS;
  double b  = SEMI_MINOR_AXIS;
  double a2 = a * a;
  double b2 = b * b;
  double e2 = (a2 - b2) / a2;

  if ( lat < -90 ) lat = -90;
  if ( lat >  90 ) lat = 90;

  double rlon = (lon + MERIDIAN_OFFSET) * (M_PI/180);
  double rlat = lat * (M_PI/180);
  double slat = sin( rlat );
  double clat = cos( rlat );
  double slon = sin( rlon );
  double clon = cos( rlon );
  double radius = a / sqrt(1.0-e2*slat*slat);

  x = (radius+h) * clat * clon;
  y = (radius+h) * clat * slon;
  z = (radius*(1-e2)+h) * slat;
}

void cartesianToGeodetic(double x, double y, double z, double &lon, double &lat, double &h) 
{
  const double a2 = SEMI_MAJOR_AXIS * SEMI_MAJOR_AXIS;
  const double b2 = SEMI_MINOR_AXIS * SEMI_MINOR_AXIS;
  const double e2 = 1 - b2 / a2;
  const double e4 = e2 * e2;

  double xy_dist = sqrt( x * x + y * y );
  double p = ( x * x + y * y ) / a2;
  double q = ( 1 - e2 ) * z * z / a2;
  double r = ( p + q - e4 ) / 6.0;
  double r3 = r * r * r;

  double evolute = 8 * r3 + e4 * p * q;
  double u = std::numeric_limits<double>::quiet_NaN();
  if ( evolute > 0 ) {
    // outside the evolute
    double right_inside_pow = sqrt(e4 * p * q);
    double sqrt_evolute = sqrt( evolute );
    u = r + 0.5 * pow(sqrt_evolute + right_inside_pow,2.0/3.0) +
      0.5 * pow(sqrt_evolute - right_inside_pow,2.0/3.0);
  } else if ( fabs(z) < std::numeric_limits<double>::epsilon() ) {
    // On the equator plane
    lat = 0;
    h = sqrt(x*x+y*y+z*z) - SEMI_MAJOR_AXIS;
  } else if ( evolute < 0 && fabs(q) > std::numeric_limits<double>::epsilon() ) {
    // On or inside the evolute
    double atan_result = atan2( sqrt( e4 * p * q ), sqrt( -evolute ) + sqrt(-8 * r3) );
    u = -4 * r * sin( 2.0 / 3.0 * atan_result ) *
      cos( M_PI / 6.0 + 2.0 / 3.0 * atan_result );
  } else if ( fabs(q) < std::numeric_limits<double>::epsilon() && p <= e4 ) {
    // In the singular disc
    h = -SEMI_MAJOR_AXIS * sqrt(1 - e2) * sqrt(e2 - p) / sqrt(e2);
    lat = 2 * atan2( sqrt(e4 - p), sqrt(e2*(e2 - p)) + sqrt(1-e2) * sqrt(p) );
  } else {
    // Near the cusps of the evolute
    double inside_pow = sqrt(evolute) + sqrt(e4 * p * q);
    u = r + 0.5 * pow(inside_pow,2.0/3.0) +
      2 * r * r * pow(inside_pow,-2.0/3.0);
  }

  {
    double v   = sqrt( u * u + e4 * q );
    double u_v = u + v;
    double w   = e2 * ( u_v - q ) / ( 2 * v );
    double k   = u_v / ( w + sqrt( w * w + u_v ) );
    double D   = k * xy_dist / ( k + e2 );
    double dist_2 = D * D + z * z;
    h = ( k + e2 - 1 ) * sqrt( dist_2 ) / k;
    lat = 2 * atan2( z, sqrt( dist_2 ) + D );
  }

  if ( xy_dist + x > ( sqrt(2) - 1 ) * y ) {
    // Longitude is between -135 and 135
    lon = 360.0 * atan2( y, xy_dist + x ) / M_PI;
  } else if ( xy_dist + y < ( sqrt(2) + 1 ) * x ) {
    // Longitude is between -225 and 45
    lon = - 90.0 + 360.0 * atan2( x, xy_dist - y ) / M_PI;
  } else {
    // Longitude is between -45 and 225
    lon = 90.0 - 360.0 * atan2( x, xy_dist + y ) / M_PI;
  }
  lon -= MERIDIAN_OFFSET;
  lat *= 180.0 / M_PI;
}

double deg2rad(double deg) 
{
  return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) 
{
  return (rad * 180 / M_PI);
}

double haversDist(double lat1d, double lon1d, double lat2d, double lon2d) 
{
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(lat1d);
  lon1r = deg2rad(lon1d);
  lat2r = deg2rad(lat2d);
  lon2r = deg2rad(lon2d);
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  const double earthRadiusKm = 6371.0;
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

#ifdef NMAIN
bool COMMUNICATION_VIA_FILES = true;

class Communicator{public:
string message_folder = "/chunks/";
int message_count_in = 0;
int message_count_out = 0;
string message_file_in = "message_in";
string semaphore_file_in = "semaphore_in";
string message_file_out = "message_out";
string semaphore_file_out = "semaphore_out";
inline bool fileExists(const string& name){
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
void waitForFile(string fileName){
  int delay = 1;
  while(!fileExists(fileName)){
    this_thread::sleep_for(chrono::milliseconds(delay));
    if(delay < 500) delay *= 2;
  }
}
string readFile(const string &fileName){
  ifstream ifs(fileName.c_str(), ios::in | ios::binary | ios::ate);
  ifstream::pos_type fileSize = ifs.tellg();
  ifs.seekg(0, ios::beg);
  vector<char> bytes(fileSize);
  ifs.read(bytes.data(), fileSize);
  return string(bytes.data(), fileSize);
}
void writeFile(const string &fileName, const string &content){
  ofstream out(fileName.c_str());
  out << content << endl;
  out.close();
}
void send_line(const string& text){
  if(COMMUNICATION_VIA_FILES){
    writeFile(message_folder + message_file_out + to_string(message_count_out), text);
    writeFile(message_folder + semaphore_file_out + to_string(message_count_out), "");
    message_count_out++;
  }else{
    cout << text << endl;
    cout.flush();	
  }
}
string receive_line(){
  string text;
  if(COMMUNICATION_VIA_FILES){
    waitForFile(message_folder + semaphore_file_in + to_string(message_count_in));
    text = readFile(message_folder + message_file_in + to_string(message_count_in));
    message_count_in++;
  }else{
    cin >> text;
  }
  return text;
}	
};
#endif

string trim(const string& str, const string& whitespace=" \t\r\n"){
  size_t strBegin = str.find_first_not_of(whitespace);
  if (strBegin == string::npos) return "";
  size_t strEnd = str.find_last_not_of(whitespace);
  size_t strRange = strEnd - strBegin + 1;
  return str.substr(strBegin, strRange);
}
vector <string> split(const string& text, char by = ','){
  vector <string> ret;
  stringstream ss(text);
  string word;
  while(getline(ss, word, by)){
    ret.push_back(word);
  }
  return ret;
} 

struct navline_t {
  // input
  double time;    // sec
  double lat;     // deg
  double lon;     // deg
  double h;       // m, msl
  double heading; // rad
  double pitch;   // rad
  double speed;   // m/s

  // filled
  double x;
  double y;
  double z;
  double frecv;
};

double getNavDist(const navline_t &lhs, const navline_t &rhs)
{
  return sqrt(SQR(lhs.x - rhs.x) + SQR(lhs.y - rhs.y) + SQR(lhs.z - rhs.z));
}

class GeolocatorBase {
public:
  void onNewChunk(double fEmitted, double totalLag, const string &name);
  ppd getLatLon() const { return m_latLon; }
  virtual void onChunkReady() = 0;
protected:
  ppd m_latLon;
  double m_fEmitted;
  double m_fReceiver;
  vector<short> m_iqData;
  vector<navline_t> m_navlineArr;
};

void GeolocatorBase::onNewChunk(double fEmitted, double totalLag, const string &name)
{
  m_fEmitted = fEmitted;

  string line;

  {
    ifstream in;
    in.open((name + ".log").c_str());   // reading the receiver log file
    for(int i = 0; i < 17; i++) getline(in, line);
    m_fReceiver = atoi(line.substr(18).c_str()); // extracting receiver frequency
  }

  {
    ifstream in;
    in.open((name + ".raw").c_str(), ios::binary);   // reading the receiver raw file
    assert(!in.bad() && !in.eof());
    int N = (int)IQ_SIZE_SEC * NSAMPLES_PER_SEC * 2;
    m_iqData.resize(N);
    in.read((char *)&m_iqData[0], N * sizeof(short));
  }

  {
    m_navlineArr.clear();
    ifstream in((name + ".csv").c_str());    // reading the flight path data
    string line;
    getline(in, line); // ignoring header
    vector<int> cntArr, timeArr;
    while (true) {
      getline(in, line); // reading the first data line
      if (in.eof()) break;
      vector <string> items = split(line, ',');

      navline_t navline;
      int timeNav = atoi(items[3].substr(0, 2).c_str()) * 3600 + atoi(items[3].substr(2, 2).c_str()) * 60 + atoi(items[3].substr(4, 2).c_str());

      vector <string> latString = split(items[4], ':');
      vector <string> lonString = split(items[5], ':');
      double latitude = atoi(latString[0].substr(1).c_str()) + atoi(latString[1].c_str()) / 60.0 + atof(latString[2].c_str()) / 3600.0;
      double longitude = -(atoi(lonString[0].substr(1).c_str()) + atoi(lonString[1].c_str()) / 60.0 + atof(lonString[2].c_str()) / 3600.0);

      if (cntArr.empty() || timeArr.back() != timeNav) {
        cntArr.push_back(1);
      } else {
        cntArr.back() += 1;
      }

      timeArr.push_back(timeNav);

      navline.time = timeNav;
      navline.lat = latitude;
      navline.lon = longitude;
      navline.h = atof(items[8].c_str()); // msl
      navline.heading = atof(items[7].c_str()) * DEG2RAD;
      navline.pitch = atof(items[11].c_str()) * DEG2RAD;
      navline.speed = atof(items[10].c_str()) * KNOT2MS;

      //geodeticToCartesian(navline.lon, navline.lat, navline.h, navline.x, navline.y, navline.z);
      geodeticToCartesian(navline.lon, navline.lat, 0, navline.x, navline.y, navline.z);

      m_navlineArr.push_back(navline);
    }
    
    vector<navline_t> navlineCopy;

    // interpolate time
    /*
    for (int i = 0, j = 0, offs = -1; i < m_navlineArr.size(); ++i) {
      if (i > 0 && timeArr[i - 1] != timeArr[i]) {
        j++;
        offs = 0;
      } else {
        offs++;
      }
      m_navlineArr[i].time += (double)offs / cntArr[j];
    }
    */

    timeArr.push_back(-1);
    for (int i = 0, j = 0; i <= m_navlineArr.size(); ++i) {
      if (i > 0 && timeArr[i - 1] != timeArr[i]) {
        navline_t val;
        memset(&val, 0, sizeof(val));
        for (int k = 0; k < cntArr[j]; ++k) {
          val.lat += m_navlineArr[i - 1 - k].lat;
          val.lon += m_navlineArr[i - 1 - k].lon;
          val.h += m_navlineArr[i - 1 - k].h;
          val.heading += m_navlineArr[i - 1 - k].heading;
          val.pitch += m_navlineArr[i - 1 - k].pitch;
          val.speed += m_navlineArr[i - 1 - k].speed;
        }
        val.time = floor(m_navlineArr[i - 1].time);
        val.lat /= cntArr[j];
        val.lon /= cntArr[j];
        val.h /= cntArr[j];
        val.heading /= cntArr[j];
        val.pitch /= cntArr[j];
        val.speed /= cntArr[j];

        //geodeticToCartesian(val.lon, val.lat, val.h, val.x, val.y, val.z);
        geodeticToCartesian(val.lon, val.lat, 0, val.x, val.y, val.z);

        navlineCopy.push_back(val);

        if (i == m_navlineArr.size()) {
          break;
        }

        j++;
      }
    }

    m_navlineArr.swap(navlineCopy);
  }

  onChunkReady();
}

class GeolocatorLS: public GeolocatorBase {
public:
  void onChunkReady() override;

private:
  void calculateFRecv();

private:
  vector<navline_t> m_navlineComb;
};

class GeolocatorLSGSearch: public GeolocatorBase {
public:
  void onChunkReady() override;

private:
  void calculateFRecv();

private:
  vector<navline_t> m_navlineComb;
};

void GeolocatorLS::calculateFRecv()
{
  assert(!m_navlineArr.empty());

  double stime = m_navlineArr.front().time;
  for (int i = 0; i < m_navlineArr.size(); ++i) {
    navline_t &nv = m_navlineArr[i];
    int fromIdx = (nv.time - stime) * NSAMPLES_PER_SEC;
    //int fromIdx = floor(nv.time - stime) * NSAMPLES_PER_SEC;
    int toIdx = NSAMPLES_PER_SEC * IQ_SIZE_SEC;
    if (i != m_navlineArr.size() - 1) {
      toIdx = (m_navlineArr[i + 1].time - stime) * NSAMPLES_PER_SEC;
      //toIdx = (ceil(m_navlineArr[i + 1].time) - stime) * NSAMPLES_PER_SEC;
    }
    int sz = toIdx - fromIdx;
    int sz2 = clp2(sz);
    if (i != m_navlineArr.size() - 1) {
      sz = sz2;
    }
    WrapVecDoub vec(2*sz2);
    for (int j = 0; j < sz2; ++j) {
      if (j < sz) {
        vec.real(j) = (double)m_iqData[2 * (fromIdx + j)];
        vec.imag(j) = (double)m_iqData[2 * (fromIdx + j) + 1];
      } else {
        vec.real(j) = vec.imag(j) = 0.0;
      }
    }
    four1((Doub*)(&vec[0]),sz2,1);
    double maxAmp = 0.0;
    double maxAmpFreq = 0.0; // hz
    double maxIdx = 0;

    //char buf[128];
    //sprintf(buf, "%d.out", i);
    //FILE *fd = fopen(buf, "w");
    //fprintf(fd, "idx,freq\n");
    for (int j = -sz2/2; j <sz2/2; ++j) {
      double amp = log(std::abs(vec[j]));
      if (amp > maxAmp) {
        maxAmp = amp;

        if (j > 0) {
          maxAmpFreq = (double)j / sz2 * NSAMPLES_PER_SEC;
        } else {
          maxAmpFreq = -double(j + sz2/2) / sz2 * NSAMPLES_PER_SEC;
        }

        maxIdx = j;
      }
      //fprintf(fd, "%d,%f\n", j, amp);
    }
    //fclose(fd);

    nv.frecv = -maxAmpFreq;
  }
}

void GeolocatorLSGSearch::calculateFRecv()
{
  assert(!m_navlineArr.empty());

  double stime = m_navlineArr.front().time;
  for (int i = 0; i < m_navlineArr.size(); ++i) {
    navline_t &nv = m_navlineArr[i];
    int fromIdx = (nv.time - stime) * NSAMPLES_PER_SEC;
    //int fromIdx = floor(nv.time - stime) * NSAMPLES_PER_SEC;
    int toIdx = NSAMPLES_PER_SEC * IQ_SIZE_SEC;
    if (i != m_navlineArr.size() - 1) {
      toIdx = (m_navlineArr[i + 1].time - stime) * NSAMPLES_PER_SEC;
      //toIdx = (ceil(m_navlineArr[i + 1].time) - stime) * NSAMPLES_PER_SEC;
    }
    int sz = toIdx - fromIdx;
    int sz2 = clp2(sz);
    if (i != m_navlineArr.size() - 1) {
      sz = sz2;
    }
    WrapVecDoub vec(2*sz2);
    for (int j = 0; j < sz2; ++j) {
      if (j < sz) {
        vec.real(j) = (double)m_iqData[2 * (fromIdx + j)];
        vec.imag(j) = (double)m_iqData[2 * (fromIdx + j) + 1];
      } else {
        vec.real(j) = vec.imag(j) = 0.0;
      }
    }
    four1((Doub*)(&vec[0]),sz2,1);
    double maxAmp = 0.0;
    double maxAmpFreq = 0.0; // hz
    double maxIdx = 0;

    //char buf[128];
    //sprintf(buf, "%d.out", i);
    //FILE *fd = fopen(buf, "w");
    //fprintf(fd, "idx,freq\n");
    for (int j = -sz2/2; j <sz2/2; ++j) {
      double amp = log(std::abs(vec[j]));
      if (amp > maxAmp) {
        maxAmp = amp;

        if (j > 0) {
          maxAmpFreq = (double)j / sz2 * NSAMPLES_PER_SEC;
        } else {
          maxAmpFreq = -double(j + sz2/2) / sz2 * NSAMPLES_PER_SEC;
        }

        maxIdx = j;
      }
      //fprintf(fd, "%d,%f\n", j, amp);
    }
    //fclose(fd);

    nv.frecv = -maxAmpFreq;
  }
}


VecDoub smdFunc(VecDoub_I &arg) { return arg; }

void GeolocatorLSGSearch::onChunkReady()
{
  calculateFRecv();

  std::copy(m_navlineArr.begin(), m_navlineArr.end(), back_inserter(m_navlineComb));

  int sz = m_navlineComb.size();
  assert(sz > 1);

  double minLat = 1e9, maxLat = -1e9, minLon = 1e9, maxLon = -1e9;
  for (auto &v: m_navlineComb) {
    minLat = min(minLat, v.lat);
    maxLat = max(maxLat, v.lat);
    minLon = min(minLon, v.lon);
    maxLon = max(maxLon, v.lon);
  }

  vector<double> rArr(sz), xArr(sz), yArr(sz), bArr(sz), cArr(sz), dArr(sz), kArr(sz);

  double wKm = haversDist(minLat, minLon, minLat, maxLon);
  double hKm = haversDist(minLat, minLon, maxLat, minLon);
  const double stepScale = 0.5;
  double kmPerLat = (maxLat - minLat) / hKm;
  double kmPerLon = (maxLon - minLon) / wKm;
  const double kmWide = 75; //km
  const double latLimit = kmWide * kmPerLat;
  const double lonLimit = kmWide * kmPerLon;
  double latCenter = (maxLat + minLat) / 2;
  double lonCenter = (maxLon + minLon) / 2;
  double lat0 = latCenter - latLimit;
  double lon0 = lonCenter - lonLimit;
  double avgSpeed = 0.0;

  for (int i = 0; i < sz; ++i) {
    auto &v = m_navlineComb[i];
    xArr[i] = (v.lon - lon0) / (2*lonLimit) * kmWide * 2 * 1000;
    yArr[i] = (v.lat - lat0) / (2*latLimit) * kmWide * 2 * 1000;
    avgSpeed += v.speed;
  }

  avgSpeed /= sz;

  MatDoub matA(sz, 3);
  VecDoub vecB(sz, 1.0);
  VecDoub vecSsig(sz, 1.0);

  for (int i = 0; i < sz; ++i) {
    vecB[i] = m_fReceiver + m_navlineComb[i].frecv;
  }

  double bestChiSqr = -1.0, bestGLat = lat0, bestGLon = lon0;
  for (double gLat = latCenter - latLimit; gLat < latCenter + latLimit; gLat += stepScale * kmPerLat) {
    for (double gLon = lonCenter - lonLimit; gLon < lonCenter + lonLimit; gLon += stepScale * kmPerLon) {
      double xcur = (gLon - lon0) / (2*lonLimit) * kmWide * 2 * 1000;
      double ycur = (gLat - lat0) / (2*latLimit) * kmWide * 2 * 1000;

      //if (gLat > 30.670 && gLat < 30.700 && gLon > -85.300 && gLon < -85.100) {
      //  static volatile bool ok = true;
      //  ok = false;
      //}

      for (int i = 0; i < sz; ++i) {
        rArr[i] = sqrt(SQR(xArr[i] - xcur) + SQR(yArr[i] - ycur));
      }

      for (int i = 0; i < sz; ++i) {
        matA[i][0] = (xcur - xArr[i]) / rArr[i] / SPEED_OF_LIGHT;
        matA[i][1] = (ycur - yArr[i]) / rArr[i] / SPEED_OF_LIGHT;
        matA[i][2] = 1;
      }

      Fitsvd lsq(matA, vecB, vecSsig, &smdFunc);
      lsq.fit();

      double femitted = lsq.a[2];
      double xsol = lsq.a[0] / femitted;
      double ysol = lsq.a[1] / femitted;

      const double freqLimit = 500; // hz
      if (abs(femitted - m_fEmitted) > freqLimit) {
        continue;
      }

      double speed = sqrt(SQR(xsol) + SQR(ysol));
      if (abs(speed - avgSpeed) > 10) {
        continue;
      }

      if (bestChiSqr < 0 || bestChiSqr > lsq.chisq) {
        bestChiSqr = lsq.chisq;
        bestGLon = gLon;
        bestGLat = gLat;
      }
    }
  }

  printf("err=%f,lat=%f,lon=%f\n", bestChiSqr / sz, bestGLat, bestGLon);
}

void GeolocatorLS::onChunkReady()
{
  calculateFRecv();

  std::copy(m_navlineArr.begin(), m_navlineArr.end(), back_inserter(m_navlineComb));

  int sz = m_navlineComb.size();
  assert(sz > 1);
  vector<double> alphaArr(sz), distArr(sz), rArr(sz), bArr(sz), cArr(sz), dArr(sz), kArr(sz);

  for (int i = 0; i < sz; ++i) {
    auto &v = m_navlineComb[i];
    bArr[i] = cos(v.heading);
    cArr[i] = sin(v.heading);
    kArr[i] = -(bArr[i] * v.x + cArr[i] * v.y);
    //kArr[i] = -(bArr[i] * v.lon + cArr[i] * v.lat);

    //bArr[i] = cos(v.heading) * cos(v.pitch);
    //cArr[i] = sin(v.heading) * cos(v.pitch);
    //dArr[i] = sin(v.pitch);
    //kArr[i] = -(v.x * bArr[i] + v.y * cArr[i] + v.z * dArr[i]);

    //if (i != sz - 1) {
    //  auto &v1 = m_navlineComb[i+1];
    //  MyVector3d vec(v1.x-v.x, v1.y-v.y, v1.z-v.z);
    //  double norm = MyVectorNormalize(vec);
    //  if (i == 0 || norm > 1e-9) {
    //    bArr[i] = vec.p[0];
    //    cArr[i] = vec.p[1];
    //    dArr[i] = vec.p[2];
    //  } else {
    //    bArr[i] = bArr[i-1];
    //    cArr[i] = cArr[i-1];
    //    dArr[i] = dArr[i-1];
    //  }

    //} else {
    //  bArr[i] = bArr[i-1];
    //  cArr[i] = cArr[i-1];
    //  dArr[i] = dArr[i-1];
    //}

    //kArr[i] = -(v.x * bArr[i] + v.y * cArr[i] + v.z * dArr[i]);
  }

  //double deltaFreq=  m_fEmitted * m_navlineComb.front().speed / SPEED_OF_LIGHT;
  double deltaFreq = 500; // hz
  for (int i = -F0_GRID_STEPS; i <= F0_GRID_STEPS; ++i) {
    double deltaF = ((double)i / F0_GRID_STEPS) * deltaFreq;
    double fEmit = m_fEmitted + deltaF;

    bool bad = false;
    for (int j = 0; j < sz; ++j) {
      auto &v = m_navlineComb[j];
      double v0 = ((m_fReceiver + v.frecv)/ fEmit - 1.0);
      double v1 = SPEED_OF_LIGHT / v.speed * v0;
      if (abs(v1) > 1.0) {
        bad = true;
        break;
      }
      alphaArr[j] = acos(v1);
      distArr[j] = j + 1 == m_navlineComb.size() ? v.speed: getNavDist(m_navlineComb[j], m_navlineComb[j + 1]);
      //distArr[j] = v.speed * 1;
      if (j != 0) {
        distArr[j] += distArr[j - 1];
      }
    }

    if (bad) {
      continue;
    }

    if (alphaArr[0] < alphaArr[1]) {
      rArr[0] = distArr[0] * sin(alphaArr[1]) / sin(alphaArr[1] - alphaArr[0]);
      for (int j = 1; j < sz; ++j) {
        rArr[j] = distArr[j - 1] * sin(alphaArr[0]) / sin(alphaArr[j] - alphaArr[0]);
      }
    } else {
      rArr[0] = distArr[0] * cos(alphaArr[0]) / sin(alphaArr[0] - alphaArr[1]);
      for (int j = 1; j < sz; ++j) {
        rArr[j] = distArr[j - 1] * sin(alphaArr[0]) / sin(alphaArr[0] - alphaArr[j]);
      }
    }
    
    for (int j = 0; j < sz; ++j) {
      rArr[j] = abs(rArr[j]);
    }

    MatDoub matA(sz, 2);
    VecDoub vecB(sz, 1.0);
    VecDoub vecSsig(sz, 1.0);

    for (int j = 0; j < sz; ++j) {
      matA[j][0] = bArr[j];
      matA[j][1] = cArr[j];
      //matA[j][2] = dArr[j];
      //char buf[1024];
      //sprintf(buf, "%f,%f,%f;",bArr[j],cArr[j],dArr[j]);
      
      vecB[j] = rArr[j] * cos(alphaArr[j]) - kArr[j];
      //sprintf(buf, "%f,", vecB[j]);
      //OutputDebugString(buf);
    }

    Fitsvd lsq(matA, vecB, vecSsig, &smdFunc);
    lsq.fit();

    double lon = lsq.a[0];
    double lat = lsq.a[1];
    //double lon, lat, h;
    //cartesianToGeodetic(x, y, z, lon, lat, h);

    printf("%f,%f\n", lon, lat);
  }
}

int main(int argc, char* argv[]){
  //GeolocatorLS geo;
  GeolocatorLSGSearch geo;

#ifdef NMAIN
  Communicator com;
  if(argc > 1) com.message_folder = argv[1];
  if(!com.message_folder.empty() && com.message_folder.back() != '\\' && com.message_folder.back() != '/') com.message_folder += "/";
  string communication_type = argc > 2 ? argv[2] : "files";
  if(communication_type != "files") COMMUNICATION_VIA_FILES = false;
  while(true){
    vector <string> message = split(trim(com.receive_line()));
    if(message.size() < 3){
      cerr << "Incorrect format of the message!" << endl;
      return 1;
    }
    int f_emitter = atoi(message[0].c_str()); 
    string total_lag = message[1];
    string STR = message[2];		
    if(STR == "end"){
      com.send_line("end");
      break;
    }

    geo.onNewChunk(f_emitter, atof(total_lag.c_str()), STR);
    auto latLon = geo.getLatLon();

    string lat = to_string(latLon.first);
    string lon = to_string(latLon.second);
    cerr << lat << "," << lon << endl;
    com.send_line(lat + "," + lon + ",unsolved");  // returning the position of the aircraft at the beginning of the 5 sec. interval
  }
#else

  geo.onNewChunk(134475360, 0, "chunks/chunk0");
  geo.onNewChunk(134475360, 0, "chunks/chunk1");
  geo.onNewChunk(134475360, 0, "chunks/chunk2");
  geo.onNewChunk(134475360, 0, "chunks/chunk3");
  geo.onNewChunk(134475360, 0, "chunks/chunk4");

  //geo.onNewChunk(27325000, 0, "chunk0");
  //geo.onNewChunk(27325000, 1, "chunk1");
  //geo.onNewChunk(27325000, 2, "chunk2");
  //geo.onNewChunk(27325000, 3, "chunk3");
  //geo.onNewChunk(27325000, 4, "chunk4");

#endif


  return 0;
}
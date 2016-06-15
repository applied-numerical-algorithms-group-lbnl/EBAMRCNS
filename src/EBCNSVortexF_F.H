#ifndef _EBCNSVORTEXF_F_H_
#define _EBCNSVORTEXF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCNSVORTEX 
#define GUARDSETCNSVORTEX 
// Prototype for Fortran procedure setcnsvortex ...
//
void FORTRAN_NAME( SETCNSVORTEX ,setcnsvortex )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(ausmallness) );

#define FORT_SETCNSVORTEX FORTRAN_NAME( inlineSETCNSVORTEX, inlineSETCNSVORTEX)
#define FORTNT_SETCNSVORTEX FORTRAN_NAME( SETCNSVORTEX, setcnsvortex)

inline void FORTRAN_NAME(inlineSETCNSVORTEX, inlineSETCNSVORTEX)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(ausmallness) )
{
 CH_TIMELEAF("FORT_SETCNSVORTEX");
 FORTRAN_NAME( SETCNSVORTEX ,setcnsvortex )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REALVECT(aucenter)
      ,CHFt_CONST_REAL(aumach)
      ,CHFt_CONST_REAL(aurnot)
      ,CHFt_CONST_REAL(ausmallness) );
}
#endif  // GUARDSETCNSVORTEX 

#ifndef GUARDCNSVORTEXINIT 
#define GUARDCNSVORTEXINIT 
// Prototype for Fortran procedure cnsvortexinit ...
//
void FORTRAN_NAME( CNSVORTEXINIT ,cnsvortexinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_CNSVORTEXINIT FORTRAN_NAME( inlineCNSVORTEXINIT, inlineCNSVORTEXINIT)
#define FORTNT_CNSVORTEXINIT FORTRAN_NAME( CNSVORTEXINIT, cnsvortexinit)

inline void FORTRAN_NAME(inlineCNSVORTEXINIT, inlineCNSVORTEXINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CNSVORTEXINIT");
 FORTRAN_NAME( CNSVORTEXINIT ,cnsvortexinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCNSVORTEXINIT 

#ifndef GUARDGETMAXVEL 
#define GUARDGETMAXVEL 
// Prototype for Fortran procedure getmaxvel ...
//
void FORTRAN_NAME( GETMAXVEL ,getmaxvel )(
      CHFp_REAL(aumaxvel) );

#define FORT_GETMAXVEL FORTRAN_NAME( inlineGETMAXVEL, inlineGETMAXVEL)
#define FORTNT_GETMAXVEL FORTRAN_NAME( GETMAXVEL, getmaxvel)

inline void FORTRAN_NAME(inlineGETMAXVEL, inlineGETMAXVEL)(
      CHFp_REAL(aumaxvel) )
{
 CH_TIMELEAF("FORT_GETMAXVEL");
 FORTRAN_NAME( GETMAXVEL ,getmaxvel )(
      CHFt_REAL(aumaxvel) );
}
#endif  // GUARDGETMAXVEL 

#ifndef GUARDPOINTCNSVORTEXINIT 
#define GUARDPOINTCNSVORTEXINIT 
// Prototype for Fortran procedure pointcnsvortexinit ...
//
void FORTRAN_NAME( POINTCNSVORTEXINIT ,pointcnsvortexinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTCNSVORTEXINIT FORTRAN_NAME( inlinePOINTCNSVORTEXINIT, inlinePOINTCNSVORTEXINIT)
#define FORTNT_POINTCNSVORTEXINIT FORTRAN_NAME( POINTCNSVORTEXINIT, pointcnsvortexinit)

inline void FORTRAN_NAME(inlinePOINTCNSVORTEXINIT, inlinePOINTCNSVORTEXINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTCNSVORTEXINIT");
 FORTRAN_NAME( POINTCNSVORTEXINIT ,pointcnsvortexinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTCNSVORTEXINIT 

}

#endif

#ifndef _EBEXPLOSIONF_F_H_
#define _EBEXPLOSIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETEXPLOSION 
#define GUARDSETEXPLOSION 
// Prototype for Fortran procedure setexplosion ...
//
void FORTRAN_NAME( SETEXPLOSION ,setexplosion )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(ausize)
      ,CHFp_CONST_REAL(aupnot)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupone)
      ,CHFp_CONST_REAL(aurone)
      ,CHFp_CONST_REAL(ausmallness) );

#define FORT_SETEXPLOSION FORTRAN_NAME( inlineSETEXPLOSION, inlineSETEXPLOSION)
#define FORTNT_SETEXPLOSION FORTRAN_NAME( SETEXPLOSION, setexplosion)

inline void FORTRAN_NAME(inlineSETEXPLOSION, inlineSETEXPLOSION)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(ausize)
      ,CHFp_CONST_REAL(aupnot)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(aupone)
      ,CHFp_CONST_REAL(aurone)
      ,CHFp_CONST_REAL(ausmallness) )
{
 CH_TIMELEAF("FORT_SETEXPLOSION");
 FORTRAN_NAME( SETEXPLOSION ,setexplosion )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REALVECT(aucenter)
      ,CHFt_CONST_REAL(ausize)
      ,CHFt_CONST_REAL(aupnot)
      ,CHFt_CONST_REAL(aurnot)
      ,CHFt_CONST_REAL(aupone)
      ,CHFt_CONST_REAL(aurone)
      ,CHFt_CONST_REAL(ausmallness) );
}
#endif  // GUARDSETEXPLOSION 

#ifndef GUARDEXPLOSIONINIT 
#define GUARDEXPLOSIONINIT 
// Prototype for Fortran procedure explosioninit ...
//
void FORTRAN_NAME( EXPLOSIONINIT ,explosioninit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(oned)
      ,CHFp_CONST_REALVECT(onednorm)
      ,CHFp_CONST_REALVECT(onedorig)
      ,CHFp_BOX(box) );

#define FORT_EXPLOSIONINIT FORTRAN_NAME( inlineEXPLOSIONINIT, inlineEXPLOSIONINIT)
#define FORTNT_EXPLOSIONINIT FORTRAN_NAME( EXPLOSIONINIT, explosioninit)

inline void FORTRAN_NAME(inlineEXPLOSIONINIT, inlineEXPLOSIONINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(oned)
      ,CHFp_CONST_REALVECT(onednorm)
      ,CHFp_CONST_REALVECT(onedorig)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXPLOSIONINIT");
 FORTRAN_NAME( EXPLOSIONINIT ,explosioninit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(oned)
      ,CHFt_CONST_REALVECT(onednorm)
      ,CHFt_CONST_REALVECT(onedorig)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXPLOSIONINIT 

#ifndef GUARDPOINTEXPLOSIONINIT 
#define GUARDPOINTEXPLOSIONINIT 
// Prototype for Fortran procedure pointexplosioninit ...
//
void FORTRAN_NAME( POINTEXPLOSIONINIT ,pointexplosioninit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(oned)
      ,CHFp_CONST_REALVECT(oneDNorm)
      ,CHFp_CONST_REALVECT(oneDOrig) );

#define FORT_POINTEXPLOSIONINIT FORTRAN_NAME( inlinePOINTEXPLOSIONINIT, inlinePOINTEXPLOSIONINIT)
#define FORTNT_POINTEXPLOSIONINIT FORTRAN_NAME( POINTEXPLOSIONINIT, pointexplosioninit)

inline void FORTRAN_NAME(inlinePOINTEXPLOSIONINIT, inlinePOINTEXPLOSIONINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(oned)
      ,CHFp_CONST_REALVECT(oneDNorm)
      ,CHFp_CONST_REALVECT(oneDOrig) )
{
 CH_TIMELEAF("FORT_POINTEXPLOSIONINIT");
 FORTRAN_NAME( POINTEXPLOSIONINIT ,pointexplosioninit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(oned)
      ,CHFt_CONST_REALVECT(oneDNorm)
      ,CHFt_CONST_REALVECT(oneDOrig) );
}
#endif  // GUARDPOINTEXPLOSIONINIT 

}

#endif

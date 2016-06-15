#ifndef _SMOOTHBUMPF_F_H_
#define _SMOOTHBUMPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETSMOOTHBUMP 
#define GUARDSETSMOOTHBUMP 
// Prototype for Fortran procedure setsmoothbump ...
//
void FORTRAN_NAME( SETSMOOTHBUMP ,setsmoothbump )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(ausmallness) );

#define FORT_SETSMOOTHBUMP FORTRAN_NAME( inlineSETSMOOTHBUMP, inlineSETSMOOTHBUMP)
#define FORTNT_SETSMOOTHBUMP FORTRAN_NAME( SETSMOOTHBUMP, setsmoothbump)

inline void FORTRAN_NAME(inlineSETSMOOTHBUMP, inlineSETSMOOTHBUMP)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REALVECT(aucenter)
      ,CHFp_CONST_REAL(aumach)
      ,CHFp_CONST_REAL(aurnot)
      ,CHFp_CONST_REAL(ausmallness) )
{
 CH_TIMELEAF("FORT_SETSMOOTHBUMP");
 FORTRAN_NAME( SETSMOOTHBUMP ,setsmoothbump )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REALVECT(aucenter)
      ,CHFt_CONST_REAL(aumach)
      ,CHFt_CONST_REAL(aurnot)
      ,CHFt_CONST_REAL(ausmallness) );
}
#endif  // GUARDSETSMOOTHBUMP 

#ifndef GUARDSMOOTHBUMPINIT 
#define GUARDSMOOTHBUMPINIT 
// Prototype for Fortran procedure smoothbumpinit ...
//
void FORTRAN_NAME( SMOOTHBUMPINIT ,smoothbumpinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SMOOTHBUMPINIT FORTRAN_NAME( inlineSMOOTHBUMPINIT, inlineSMOOTHBUMPINIT)
#define FORTNT_SMOOTHBUMPINIT FORTRAN_NAME( SMOOTHBUMPINIT, smoothbumpinit)

inline void FORTRAN_NAME(inlineSMOOTHBUMPINIT, inlineSMOOTHBUMPINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SMOOTHBUMPINIT");
 FORTRAN_NAME( SMOOTHBUMPINIT ,smoothbumpinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSMOOTHBUMPINIT 

#ifndef GUARDPOINTSMOOTHBUMPINIT 
#define GUARDPOINTSMOOTHBUMPINIT 
// Prototype for Fortran procedure pointsmoothbumpinit ...
//
void FORTRAN_NAME( POINTSMOOTHBUMPINIT ,pointsmoothbumpinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTSMOOTHBUMPINIT FORTRAN_NAME( inlinePOINTSMOOTHBUMPINIT, inlinePOINTSMOOTHBUMPINIT)
#define FORTNT_POINTSMOOTHBUMPINIT FORTRAN_NAME( POINTSMOOTHBUMPINIT, pointsmoothbumpinit)

inline void FORTRAN_NAME(inlinePOINTSMOOTHBUMPINIT, inlinePOINTSMOOTHBUMPINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTSMOOTHBUMPINIT");
 FORTRAN_NAME( POINTSMOOTHBUMPINIT ,pointsmoothbumpinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTSMOOTHBUMPINIT 

}

#endif

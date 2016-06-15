#ifndef _EBPLANARSHOCKF_F_H_
#define _EBPLANARSHOCKF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETPOSTSHOCKVEL 
#define GUARDGETPOSTSHOCKVEL 
// Prototype for Fortran procedure getpostshockvel ...
//
void FORTRAN_NAME( GETPOSTSHOCKVEL ,getpostshockvel )(
      CHFp_REAL(aupostshockvel) );

#define FORT_GETPOSTSHOCKVEL FORTRAN_NAME( inlineGETPOSTSHOCKVEL, inlineGETPOSTSHOCKVEL)
#define FORTNT_GETPOSTSHOCKVEL FORTRAN_NAME( GETPOSTSHOCKVEL, getpostshockvel)

inline void FORTRAN_NAME(inlineGETPOSTSHOCKVEL, inlineGETPOSTSHOCKVEL)(
      CHFp_REAL(aupostshockvel) )
{
 CH_TIMELEAF("FORT_GETPOSTSHOCKVEL");
 FORTRAN_NAME( GETPOSTSHOCKVEL ,getpostshockvel )(
      CHFt_REAL(aupostshockvel) );
}
#endif  // GUARDGETPOSTSHOCKVEL 

#ifndef GUARDGETPOSTSHOCKTEMP 
#define GUARDGETPOSTSHOCKTEMP 
// Prototype for Fortran procedure getpostshocktemp ...
//
void FORTRAN_NAME( GETPOSTSHOCKTEMP ,getpostshocktemp )(
      CHFp_REAL(aupostshocktemp) );

#define FORT_GETPOSTSHOCKTEMP FORTRAN_NAME( inlineGETPOSTSHOCKTEMP, inlineGETPOSTSHOCKTEMP)
#define FORTNT_GETPOSTSHOCKTEMP FORTRAN_NAME( GETPOSTSHOCKTEMP, getpostshocktemp)

inline void FORTRAN_NAME(inlineGETPOSTSHOCKTEMP, inlineGETPOSTSHOCKTEMP)(
      CHFp_REAL(aupostshocktemp) )
{
 CH_TIMELEAF("FORT_GETPOSTSHOCKTEMP");
 FORTRAN_NAME( GETPOSTSHOCKTEMP ,getpostshocktemp )(
      CHFt_REAL(aupostshocktemp) );
}
#endif  // GUARDGETPOSTSHOCKTEMP 

#ifndef GUARDGETSHOCKLOCX 
#define GUARDGETSHOCKLOCX 
// Prototype for Fortran procedure getshocklocx ...
//
void FORTRAN_NAME( GETSHOCKLOCX ,getshocklocx )(
      CHFp_CONST_REAL(aushockloc)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETSHOCKLOCX FORTRAN_NAME( inlineGETSHOCKLOCX, inlineGETSHOCKLOCX)
#define FORTNT_GETSHOCKLOCX FORTRAN_NAME( GETSHOCKLOCX, getshocklocx)

inline void FORTRAN_NAME(inlineGETSHOCKLOCX, inlineGETSHOCKLOCX)(
      CHFp_CONST_REAL(aushockloc)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETSHOCKLOCX");
 FORTRAN_NAME( GETSHOCKLOCX ,getshocklocx )(
      CHFt_CONST_REAL(aushockloc)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETSHOCKLOCX 

#ifndef GUARDSETPLANARSHOCK 
#define GUARDSETPLANARSHOCK 
// Prototype for Fortran procedure setplanarshock ...
//
void FORTRAN_NAME( SETPLANARSHOCK ,setplanarshock )(
      CHFp_CONST_REAL(aucv)
      ,CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshockdense)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) );

#define FORT_SETPLANARSHOCK FORTRAN_NAME( inlineSETPLANARSHOCK, inlineSETPLANARSHOCK)
#define FORTNT_SETPLANARSHOCK FORTRAN_NAME( SETPLANARSHOCK, setplanarshock)

inline void FORTRAN_NAME(inlineSETPLANARSHOCK, inlineSETPLANARSHOCK)(
      CHFp_CONST_REAL(aucv)
      ,CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(aucenter)
      ,CHFp_CONST_REAL(aupreshockpress)
      ,CHFp_CONST_REAL(aupreshockdense)
      ,CHFp_CONST_INT(aushocknorm)
      ,CHFp_CONST_INT(aushockback) )
{
 CH_TIMELEAF("FORT_SETPLANARSHOCK");
 FORTRAN_NAME( SETPLANARSHOCK ,setplanarshock )(
      CHFt_CONST_REAL(aucv)
      ,CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aums)
      ,CHFt_CONST_REAL(aucenter)
      ,CHFt_CONST_REAL(aupreshockpress)
      ,CHFt_CONST_REAL(aupreshockdense)
      ,CHFt_CONST_INT(aushocknorm)
      ,CHFt_CONST_INT(aushockback) );
}
#endif  // GUARDSETPLANARSHOCK 

#ifndef GUARDPLANARSHOCKINIT 
#define GUARDPLANARSHOCKINIT 
// Prototype for Fortran procedure planarshockinit ...
//
void FORTRAN_NAME( PLANARSHOCKINIT ,planarshockinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_PLANARSHOCKINIT FORTRAN_NAME( inlinePLANARSHOCKINIT, inlinePLANARSHOCKINIT)
#define FORTNT_PLANARSHOCKINIT FORTRAN_NAME( PLANARSHOCKINIT, planarshockinit)

inline void FORTRAN_NAME(inlinePLANARSHOCKINIT, inlinePLANARSHOCKINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PLANARSHOCKINIT");
 FORTRAN_NAME( PLANARSHOCKINIT ,planarshockinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPLANARSHOCKINIT 

#ifndef GUARDPOINTPLANARSHOCKINIT 
#define GUARDPOINTPLANARSHOCKINIT 
// Prototype for Fortran procedure pointplanarshockinit ...
//
void FORTRAN_NAME( POINTPLANARSHOCKINIT ,pointplanarshockinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTPLANARSHOCKINIT FORTRAN_NAME( inlinePOINTPLANARSHOCKINIT, inlinePOINTPLANARSHOCKINIT)
#define FORTNT_POINTPLANARSHOCKINIT FORTRAN_NAME( POINTPLANARSHOCKINIT, pointplanarshockinit)

inline void FORTRAN_NAME(inlinePOINTPLANARSHOCKINIT, inlinePOINTPLANARSHOCKINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTPLANARSHOCKINIT");
 FORTRAN_NAME( POINTPLANARSHOCKINIT ,pointplanarshockinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTPLANARSHOCKINIT 

}

#endif

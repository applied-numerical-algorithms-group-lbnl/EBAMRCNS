#ifndef _EBINFLOWOUTFLOWF_F_H_
#define _EBINFLOWOUTFLOWF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETINFLOWOUTFLOW 
#define GUARDSETINFLOWOUTFLOW 
// Prototype for Fortran procedure setinflowoutflow ...
//
void FORTRAN_NAME( SETINFLOWOUTFLOW ,setinflowoutflow )(
      CHFp_CONST_REAL(aucv)
      ,CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(autempinflow)
      ,CHFp_CONST_REAL(autempinterior)
      ,CHFp_CONST_REAL(aumachinflow)
      ,CHFp_CONST_REAL(aupresinflow)
      ,CHFp_CONST_REAL(aupresinterior) );

#define FORT_SETINFLOWOUTFLOW FORTRAN_NAME( inlineSETINFLOWOUTFLOW, inlineSETINFLOWOUTFLOW)
#define FORTNT_SETINFLOWOUTFLOW FORTRAN_NAME( SETINFLOWOUTFLOW, setinflowoutflow)

inline void FORTRAN_NAME(inlineSETINFLOWOUTFLOW, inlineSETINFLOWOUTFLOW)(
      CHFp_CONST_REAL(aucv)
      ,CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(autempinflow)
      ,CHFp_CONST_REAL(autempinterior)
      ,CHFp_CONST_REAL(aumachinflow)
      ,CHFp_CONST_REAL(aupresinflow)
      ,CHFp_CONST_REAL(aupresinterior) )
{
 CH_TIMELEAF("FORT_SETINFLOWOUTFLOW");
 FORTRAN_NAME( SETINFLOWOUTFLOW ,setinflowoutflow )(
      CHFt_CONST_REAL(aucv)
      ,CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(autempinflow)
      ,CHFt_CONST_REAL(autempinterior)
      ,CHFt_CONST_REAL(aumachinflow)
      ,CHFt_CONST_REAL(aupresinflow)
      ,CHFt_CONST_REAL(aupresinterior) );
}
#endif  // GUARDSETINFLOWOUTFLOW 

#ifndef GUARDEBIOINFLOWOUTFLOWINIT 
#define GUARDEBIOINFLOWOUTFLOWINIT 
// Prototype for Fortran procedure ebioinflowoutflowinit ...
//
void FORTRAN_NAME( EBIOINFLOWOUTFLOWINIT ,ebioinflowoutflowinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_EBIOINFLOWOUTFLOWINIT FORTRAN_NAME( inlineEBIOINFLOWOUTFLOWINIT, inlineEBIOINFLOWOUTFLOWINIT)
#define FORTNT_EBIOINFLOWOUTFLOWINIT FORTRAN_NAME( EBIOINFLOWOUTFLOWINIT, ebioinflowoutflowinit)

inline void FORTRAN_NAME(inlineEBIOINFLOWOUTFLOWINIT, inlineEBIOINFLOWOUTFLOWINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBIOINFLOWOUTFLOWINIT");
 FORTRAN_NAME( EBIOINFLOWOUTFLOWINIT ,ebioinflowoutflowinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBIOINFLOWOUTFLOWINIT 

#ifndef GUARDPOINTEBIOINFLOWOUTFLOWINIT 
#define GUARDPOINTEBIOINFLOWOUTFLOWINIT 
// Prototype for Fortran procedure pointebioinflowoutflowinit ...
//
void FORTRAN_NAME( POINTEBIOINFLOWOUTFLOWINIT ,pointebioinflowoutflowinit )(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTEBIOINFLOWOUTFLOWINIT FORTRAN_NAME( inlinePOINTEBIOINFLOWOUTFLOWINIT, inlinePOINTEBIOINFLOWOUTFLOWINIT)
#define FORTNT_POINTEBIOINFLOWOUTFLOWINIT FORTRAN_NAME( POINTEBIOINFLOWOUTFLOWINIT, pointebioinflowoutflowinit)

inline void FORTRAN_NAME(inlinePOINTEBIOINFLOWOUTFLOWINIT, inlinePOINTEBIOINFLOWOUTFLOWINIT)(
      CHFp_REAL(density)
      ,CHFp_REALVECT(momentum)
      ,CHFp_REAL(energy)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTEBIOINFLOWOUTFLOWINIT");
 FORTRAN_NAME( POINTEBIOINFLOWOUTFLOWINIT ,pointebioinflowoutflowinit )(
      CHFt_REAL(density)
      ,CHFt_REALVECT(momentum)
      ,CHFt_REAL(energy)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTEBIOINFLOWOUTFLOWINIT 

#ifndef GUARDEBIOOUTFLOWBC 
#define GUARDEBIOOUTFLOWBC 
// Prototype for Fortran procedure EBIOOUTFLOWBC ...
//
void FORTRAN_NAME( EBIOOUTFLOWBC ,ebiooutflowbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_EBIOOUTFLOWBC FORTRAN_NAME( inlineEBIOOUTFLOWBC, inlineEBIOOUTFLOWBC)
#define FORTNT_EBIOOUTFLOWBC FORTRAN_NAME( EBIOOUTFLOWBC, ebiooutflowbc)

inline void FORTRAN_NAME(inlineEBIOOUTFLOWBC, inlineEBIOOUTFLOWBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBIOOUTFLOWBC");
 FORTRAN_NAME( EBIOOUTFLOWBC ,ebiooutflowbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBIOOUTFLOWBC 

#ifndef GUARDPOINTEBIOOUTFLOWBC 
#define GUARDPOINTEBIOOUTFLOWBC 
// Prototype for Fortran procedure pointebiooutflowbc ...
//
void FORTRAN_NAME( POINTEBIOOUTFLOWBC ,pointebiooutflowbc )(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_POINTEBIOOUTFLOWBC FORTRAN_NAME( inlinePOINTEBIOOUTFLOWBC, inlinePOINTEBIOOUTFLOWBC)
#define FORTNT_POINTEBIOOUTFLOWBC FORTRAN_NAME( POINTEBIOOUTFLOWBC, pointebiooutflowbc)

inline void FORTRAN_NAME(inlinePOINTEBIOOUTFLOWBC, inlinePOINTEBIOOUTFLOWBC)(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_POINTEBIOOUTFLOWBC");
 FORTRAN_NAME( POINTEBIOOUTFLOWBC ,pointebiooutflowbc )(
      CHFt_VR(flux)
      ,CHFt_VR(prim)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDPOINTEBIOOUTFLOWBC 

#ifndef GUARDEBIOINFLOWBC 
#define GUARDEBIOINFLOWBC 
// Prototype for Fortran procedure EBIOINFLOWBC ...
//
void FORTRAN_NAME( EBIOINFLOWBC ,ebioinflowbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_EBIOINFLOWBC FORTRAN_NAME( inlineEBIOINFLOWBC, inlineEBIOINFLOWBC)
#define FORTNT_EBIOINFLOWBC FORTRAN_NAME( EBIOINFLOWBC, ebioinflowbc)

inline void FORTRAN_NAME(inlineEBIOINFLOWBC, inlineEBIOINFLOWBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBIOINFLOWBC");
 FORTRAN_NAME( EBIOINFLOWBC ,ebioinflowbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBIOINFLOWBC 

#ifndef GUARDPOINTEBIOINFLOWBC 
#define GUARDPOINTEBIOINFLOWBC 
// Prototype for Fortran procedure pointebioinflowbc ...
//
void FORTRAN_NAME( POINTEBIOINFLOWBC ,pointebioinflowbc )(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_POINTEBIOINFLOWBC FORTRAN_NAME( inlinePOINTEBIOINFLOWBC, inlinePOINTEBIOINFLOWBC)
#define FORTNT_POINTEBIOINFLOWBC FORTRAN_NAME( POINTEBIOINFLOWBC, pointebioinflowbc)

inline void FORTRAN_NAME(inlinePOINTEBIOINFLOWBC, inlinePOINTEBIOINFLOWBC)(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_POINTEBIOINFLOWBC");
 FORTRAN_NAME( POINTEBIOINFLOWBC ,pointebioinflowbc )(
      CHFt_VR(flux)
      ,CHFt_VR(prim)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDPOINTEBIOINFLOWBC 

#ifndef GUARDEBIOSLIPWALLBC 
#define GUARDEBIOSLIPWALLBC 
// Prototype for Fortran procedure EBIOSLIPWALLBC ...
//
void FORTRAN_NAME( EBIOSLIPWALLBC ,ebioslipwallbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_EBIOSLIPWALLBC FORTRAN_NAME( inlineEBIOSLIPWALLBC, inlineEBIOSLIPWALLBC)
#define FORTNT_EBIOSLIPWALLBC FORTRAN_NAME( EBIOSLIPWALLBC, ebioslipwallbc)

inline void FORTRAN_NAME(inlineEBIOSLIPWALLBC, inlineEBIOSLIPWALLBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBIOSLIPWALLBC");
 FORTRAN_NAME( EBIOSLIPWALLBC ,ebioslipwallbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBIOSLIPWALLBC 

#ifndef GUARDPOINTEBIOSLIPWALLBC 
#define GUARDPOINTEBIOSLIPWALLBC 
// Prototype for Fortran procedure pointebioslipwallbc ...
//
void FORTRAN_NAME( POINTEBIOSLIPWALLBC ,pointebioslipwallbc )(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_POINTEBIOSLIPWALLBC FORTRAN_NAME( inlinePOINTEBIOSLIPWALLBC, inlinePOINTEBIOSLIPWALLBC)
#define FORTNT_POINTEBIOSLIPWALLBC FORTRAN_NAME( POINTEBIOSLIPWALLBC, pointebioslipwallbc)

inline void FORTRAN_NAME(inlinePOINTEBIOSLIPWALLBC, inlinePOINTEBIOSLIPWALLBC)(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_POINTEBIOSLIPWALLBC");
 FORTRAN_NAME( POINTEBIOSLIPWALLBC ,pointebioslipwallbc )(
      CHFt_VR(flux)
      ,CHFt_VR(prim)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDPOINTEBIOSLIPWALLBC 

}

#endif

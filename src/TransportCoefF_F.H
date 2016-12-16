#ifndef _TRANSPORTCOEFF_F_H_
#define _TRANSPORTCOEFF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETTRANSPORTCOEFFS 
#define GUARDGETTRANSPORTCOEFFS 
// Prototype for Fortran procedure gettransportcoeffs ...
//
void FORTRAN_NAME( GETTRANSPORTCOEFFS ,gettransportcoeffs )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA1(mu)
      ,CHFp_FRA1(lambda)
      ,CHFp_FRA1(bco)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(lambda0)
      ,CHFp_CONST_REAL(bco0)
      ,CHFp_CONST_REAL(temp0) );

#define FORT_GETTRANSPORTCOEFFS FORTRAN_NAME( inlineGETTRANSPORTCOEFFS, inlineGETTRANSPORTCOEFFS)
#define FORTNT_GETTRANSPORTCOEFFS FORTRAN_NAME( GETTRANSPORTCOEFFS, gettransportcoeffs)

inline void FORTRAN_NAME(inlineGETTRANSPORTCOEFFS, inlineGETTRANSPORTCOEFFS)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA1(mu)
      ,CHFp_FRA1(lambda)
      ,CHFp_FRA1(bco)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(lambda0)
      ,CHFp_CONST_REAL(bco0)
      ,CHFp_CONST_REAL(temp0) )
{
 CH_TIMELEAF("FORT_GETTRANSPORTCOEFFS");
 FORTRAN_NAME( GETTRANSPORTCOEFFS ,gettransportcoeffs )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA1(mu)
      ,CHFt_FRA1(lambda)
      ,CHFt_FRA1(bco)
      ,CHFt_CONST_FRA1(temp)
      ,CHFt_CONST_REAL(mu0)
      ,CHFt_CONST_REAL(lambda0)
      ,CHFt_CONST_REAL(bco0)
      ,CHFt_CONST_REAL(temp0) );
}
#endif  // GUARDGETTRANSPORTCOEFFS 

#ifndef GUARDPOINTGETTRANSPORTCOEFFS 
#define GUARDPOINTGETTRANSPORTCOEFFS 
// Prototype for Fortran procedure pointgettransportcoeffs ...
//
void FORTRAN_NAME( POINTGETTRANSPORTCOEFFS ,pointgettransportcoeffs )(
      CHFp_REAL(mu)
      ,CHFp_REAL(lambda)
      ,CHFp_REAL(bco)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(lambda0)
      ,CHFp_CONST_REAL(bco0)
      ,CHFp_CONST_REAL(temp0) );

#define FORT_POINTGETTRANSPORTCOEFFS FORTRAN_NAME( inlinePOINTGETTRANSPORTCOEFFS, inlinePOINTGETTRANSPORTCOEFFS)
#define FORTNT_POINTGETTRANSPORTCOEFFS FORTRAN_NAME( POINTGETTRANSPORTCOEFFS, pointgettransportcoeffs)

inline void FORTRAN_NAME(inlinePOINTGETTRANSPORTCOEFFS, inlinePOINTGETTRANSPORTCOEFFS)(
      CHFp_REAL(mu)
      ,CHFp_REAL(lambda)
      ,CHFp_REAL(bco)
      ,CHFp_CONST_REAL(temp)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(lambda0)
      ,CHFp_CONST_REAL(bco0)
      ,CHFp_CONST_REAL(temp0) )
{
 CH_TIMELEAF("FORT_POINTGETTRANSPORTCOEFFS");
 FORTRAN_NAME( POINTGETTRANSPORTCOEFFS ,pointgettransportcoeffs )(
      CHFt_REAL(mu)
      ,CHFt_REAL(lambda)
      ,CHFt_REAL(bco)
      ,CHFt_CONST_REAL(temp)
      ,CHFt_CONST_REAL(mu0)
      ,CHFt_CONST_REAL(lambda0)
      ,CHFt_CONST_REAL(bco0)
      ,CHFt_CONST_REAL(temp0) );
}
#endif  // GUARDPOINTGETTRANSPORTCOEFFS 

}

#endif

#ifndef _REGTRANSFORMF_F_H_
#define _REGTRANSFORMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCELLCENTERTOFACECENTER 
#define GUARDCELLCENTERTOFACECENTER 
// Prototype for Fortran procedure cellcentertofacecenter ...
//
void FORTRAN_NAME( CELLCENTERTOFACECENTER ,cellcentertofacecenter )(
      CHFp_BOX(dcalcface)
      ,CHFp_FRA1(face)
      ,CHFp_CONST_FRA1(cell)
      ,CHFp_CONST_INT(facedir) );

#define FORT_CELLCENTERTOFACECENTER FORTRAN_NAME( inlineCELLCENTERTOFACECENTER, inlineCELLCENTERTOFACECENTER)
#define FORTNT_CELLCENTERTOFACECENTER FORTRAN_NAME( CELLCENTERTOFACECENTER, cellcentertofacecenter)

inline void FORTRAN_NAME(inlineCELLCENTERTOFACECENTER, inlineCELLCENTERTOFACECENTER)(
      CHFp_BOX(dcalcface)
      ,CHFp_FRA1(face)
      ,CHFp_CONST_FRA1(cell)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_CELLCENTERTOFACECENTER");
 FORTRAN_NAME( CELLCENTERTOFACECENTER ,cellcentertofacecenter )(
      CHFt_BOX(dcalcface)
      ,CHFt_FRA1(face)
      ,CHFt_CONST_FRA1(cell)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDCELLCENTERTOFACECENTER 

}

#endif

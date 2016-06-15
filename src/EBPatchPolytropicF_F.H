#ifndef _EBPATCHPOLYTROPICF_F_H_
#define _EBPATCHPOLYTROPICF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETGAMMAANDSMALL 
#define GUARDSETGAMMAANDSMALL 
// Prototype for Fortran procedure setgammaandsmall ...
//
void FORTRAN_NAME( SETGAMMAANDSMALL ,setgammaandsmall )(
      CHFp_CONST_REAL(augamma) );

#define FORT_SETGAMMAANDSMALL FORTRAN_NAME( inlineSETGAMMAANDSMALL, inlineSETGAMMAANDSMALL)
#define FORTNT_SETGAMMAANDSMALL FORTRAN_NAME( SETGAMMAANDSMALL, setgammaandsmall)

inline void FORTRAN_NAME(inlineSETGAMMAANDSMALL, inlineSETGAMMAANDSMALL)(
      CHFp_CONST_REAL(augamma) )
{
 CH_TIMELEAF("FORT_SETGAMMAANDSMALL");
 FORTRAN_NAME( SETGAMMAANDSMALL ,setgammaandsmall )(
      CHFt_CONST_REAL(augamma) );
}
#endif  // GUARDSETGAMMAANDSMALL 

#ifndef GUARDFLOORSCALAR 
#define GUARDFLOORSCALAR 
// Prototype for Fortran procedure floorscalar ...
//
void FORTRAN_NAME( FLOORSCALAR ,floorscalar )(
      CHFp_FRA1(scal)
      ,CHFp_CONST_REAL(minval)
      ,CHFp_BOX(dcalc) );

#define FORT_FLOORSCALAR FORTRAN_NAME( inlineFLOORSCALAR, inlineFLOORSCALAR)
#define FORTNT_FLOORSCALAR FORTRAN_NAME( FLOORSCALAR, floorscalar)

inline void FORTRAN_NAME(inlineFLOORSCALAR, inlineFLOORSCALAR)(
      CHFp_FRA1(scal)
      ,CHFp_CONST_REAL(minval)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_FLOORSCALAR");
 FORTRAN_NAME( FLOORSCALAR ,floorscalar )(
      CHFt_FRA1(scal)
      ,CHFt_CONST_REAL(minval)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDFLOORSCALAR 

#ifndef GUARDAIRVISCOSITY 
#define GUARDAIRVISCOSITY 
// Prototype for Fortran procedure airviscosity ...
//
void FORTRAN_NAME( AIRVISCOSITY ,airviscosity )(
      CHFp_FRA1(eta)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_BOX(dcalc) );

#define FORT_AIRVISCOSITY FORTRAN_NAME( inlineAIRVISCOSITY, inlineAIRVISCOSITY)
#define FORTNT_AIRVISCOSITY FORTRAN_NAME( AIRVISCOSITY, airviscosity)

inline void FORTRAN_NAME(inlineAIRVISCOSITY, inlineAIRVISCOSITY)(
      CHFp_FRA1(eta)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_AIRVISCOSITY");
 FORTRAN_NAME( AIRVISCOSITY ,airviscosity )(
      CHFt_FRA1(eta)
      ,CHFt_CONST_FRA1(dens)
      ,CHFt_CONST_FRA1(temp)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDAIRVISCOSITY 

#ifndef GUARDAIRTHERMDIFF 
#define GUARDAIRTHERMDIFF 
// Prototype for Fortran procedure airthermdiff ...
//
void FORTRAN_NAME( AIRTHERMDIFF ,airthermdiff )(
      CHFp_FRA1(bco)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_BOX(dcalc) );

#define FORT_AIRTHERMDIFF FORTRAN_NAME( inlineAIRTHERMDIFF, inlineAIRTHERMDIFF)
#define FORTNT_AIRTHERMDIFF FORTRAN_NAME( AIRTHERMDIFF, airthermdiff)

inline void FORTRAN_NAME(inlineAIRTHERMDIFF, inlineAIRTHERMDIFF)(
      CHFp_FRA1(bco)
      ,CHFp_CONST_FRA1(dens)
      ,CHFp_CONST_FRA1(temp)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_AIRTHERMDIFF");
 FORTRAN_NAME( AIRTHERMDIFF ,airthermdiff )(
      CHFt_FRA1(bco)
      ,CHFt_CONST_FRA1(dens)
      ,CHFt_CONST_FRA1(temp)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDAIRTHERMDIFF 

#ifndef GUARDPOINTAIRVISCOSITY 
#define GUARDPOINTAIRVISCOSITY 
// Prototype for Fortran procedure pointairviscosity ...
//
void FORTRAN_NAME( POINTAIRVISCOSITY ,pointairviscosity )(
      CHFp_REAL(eta)
      ,CHFp_CONST_REAL(dens)
      ,CHFp_CONST_REAL(tempin) );

#define FORT_POINTAIRVISCOSITY FORTRAN_NAME( inlinePOINTAIRVISCOSITY, inlinePOINTAIRVISCOSITY)
#define FORTNT_POINTAIRVISCOSITY FORTRAN_NAME( POINTAIRVISCOSITY, pointairviscosity)

inline void FORTRAN_NAME(inlinePOINTAIRVISCOSITY, inlinePOINTAIRVISCOSITY)(
      CHFp_REAL(eta)
      ,CHFp_CONST_REAL(dens)
      ,CHFp_CONST_REAL(tempin) )
{
 CH_TIMELEAF("FORT_POINTAIRVISCOSITY");
 FORTRAN_NAME( POINTAIRVISCOSITY ,pointairviscosity )(
      CHFt_REAL(eta)
      ,CHFt_CONST_REAL(dens)
      ,CHFt_CONST_REAL(tempin) );
}
#endif  // GUARDPOINTAIRVISCOSITY 

#ifndef GUARDPOINTAIRTHERMDIFF 
#define GUARDPOINTAIRTHERMDIFF 
// Prototype for Fortran procedure pointairthermdiff ...
//
void FORTRAN_NAME( POINTAIRTHERMDIFF ,pointairthermdiff )(
      CHFp_REAL(bco)
      ,CHFp_CONST_REAL(dens)
      ,CHFp_CONST_REAL(temp) );

#define FORT_POINTAIRTHERMDIFF FORTRAN_NAME( inlinePOINTAIRTHERMDIFF, inlinePOINTAIRTHERMDIFF)
#define FORTNT_POINTAIRTHERMDIFF FORTRAN_NAME( POINTAIRTHERMDIFF, pointairthermdiff)

inline void FORTRAN_NAME(inlinePOINTAIRTHERMDIFF, inlinePOINTAIRTHERMDIFF)(
      CHFp_REAL(bco)
      ,CHFp_CONST_REAL(dens)
      ,CHFp_CONST_REAL(temp) )
{
 CH_TIMELEAF("FORT_POINTAIRTHERMDIFF");
 FORTRAN_NAME( POINTAIRTHERMDIFF ,pointairthermdiff )(
      CHFt_REAL(bco)
      ,CHFt_CONST_REAL(dens)
      ,CHFt_CONST_REAL(temp) );
}
#endif  // GUARDPOINTAIRTHERMDIFF 

#ifndef GUARDEBPPGRADVEL 
#define GUARDEBPPGRADVEL 
// Prototype for Fortran procedure ebppgradvel ...
//
void FORTRAN_NAME( EBPPGRADVEL ,ebppgradvel )(
      CHFp_FRA(gradvel)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_CONST_INT(gradindex)
      ,CHFp_BOX(dcalc) );

#define FORT_EBPPGRADVEL FORTRAN_NAME( inlineEBPPGRADVEL, inlineEBPPGRADVEL)
#define FORTNT_EBPPGRADVEL FORTRAN_NAME( EBPPGRADVEL, ebppgradvel)

inline void FORTRAN_NAME(inlineEBPPGRADVEL, inlineEBPPGRADVEL)(
      CHFp_FRA(gradvel)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_CONST_INT(gradindex)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_EBPPGRADVEL");
 FORTRAN_NAME( EBPPGRADVEL ,ebppgradvel )(
      CHFt_FRA(gradvel)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(veldir)
      ,CHFt_CONST_INT(gradindex)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDEBPPGRADVEL 

#ifndef GUARDEBPPVELOHALF 
#define GUARDEBPPVELOHALF 
// Prototype for Fortran procedure ebppvelohalf ...
//
void FORTRAN_NAME( EBPPVELOHALF ,ebppvelohalf )(
      CHFp_FRA(gradvel)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_BOX(dcalc) );

#define FORT_EBPPVELOHALF FORTRAN_NAME( inlineEBPPVELOHALF, inlineEBPPVELOHALF)
#define FORTNT_EBPPVELOHALF FORTRAN_NAME( EBPPVELOHALF, ebppvelohalf)

inline void FORTRAN_NAME(inlineEBPPVELOHALF, inlineEBPPVELOHALF)(
      CHFp_FRA(gradvel)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_EBPPVELOHALF");
 FORTRAN_NAME( EBPPVELOHALF ,ebppvelohalf )(
      CHFt_FRA(gradvel)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_INT(veldir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDEBPPVELOHALF 

#ifndef GUARDINCRPDIVU 
#define GUARDINCRPDIVU 
// Prototype for Fortran procedure incrpdivu ...
//
void FORTRAN_NAME( INCRPDIVU ,incrpdivu )(
      CHFp_FRA1(pdivu)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) );

#define FORT_INCRPDIVU FORTRAN_NAME( inlineINCRPDIVU, inlineINCRPDIVU)
#define FORTNT_INCRPDIVU FORTRAN_NAME( INCRPDIVU, incrpdivu)

inline void FORTRAN_NAME(inlineINCRPDIVU, inlineINCRPDIVU)(
      CHFp_FRA1(pdivu)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_INCRPDIVU");
 FORTRAN_NAME( INCRPDIVU ,incrpdivu )(
      CHFt_FRA1(pdivu)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDINCRPDIVU 

#ifndef GUARDINCRRHOUDOTDELE 
#define GUARDINCRRHOUDOTDELE 
// Prototype for Fortran procedure incrrhoudotdele ...
//
void FORTRAN_NAME( INCRRHOUDOTDELE ,incrrhoudotdele )(
      CHFp_FRA1(rhoudotdele)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(specheat)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) );

#define FORT_INCRRHOUDOTDELE FORTRAN_NAME( inlineINCRRHOUDOTDELE, inlineINCRRHOUDOTDELE)
#define FORTNT_INCRRHOUDOTDELE FORTRAN_NAME( INCRRHOUDOTDELE, incrrhoudotdele)

inline void FORTRAN_NAME(inlineINCRRHOUDOTDELE, inlineINCRRHOUDOTDELE)(
      CHFp_FRA1(rhoudotdele)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(specheat)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_INCRRHOUDOTDELE");
 FORTRAN_NAME( INCRRHOUDOTDELE ,incrrhoudotdele )(
      CHFt_FRA1(rhoudotdele)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_REAL(specheat)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDINCRRHOUDOTDELE 

#ifndef GUARDINCRRHSTEMPERATURE 
#define GUARDINCRRHSTEMPERATURE 
// Prototype for Fortran procedure incrrhstemperature ...
//
void FORTRAN_NAME( INCRRHSTEMPERATURE ,incrrhstemperature )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(specheat)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) );

#define FORT_INCRRHSTEMPERATURE FORTRAN_NAME( inlineINCRRHSTEMPERATURE, inlineINCRRHSTEMPERATURE)
#define FORTNT_INCRRHSTEMPERATURE FORTRAN_NAME( INCRRHSTEMPERATURE, incrrhstemperature)

inline void FORTRAN_NAME(inlineINCRRHSTEMPERATURE, inlineINCRRHSTEMPERATURE)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_REAL(specheat)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_INCRRHSTEMPERATURE");
 FORTRAN_NAME( INCRRHSTEMPERATURE ,incrrhstemperature )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_REAL(specheat)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDINCRRHSTEMPERATURE 

#ifndef GUARDFILLPRIMITIVEVECTOR 
#define GUARDFILLPRIMITIVEVECTOR 
// Prototype for Fortran procedure fillprimitivevector ...
//
void FORTRAN_NAME( FILLPRIMITIVEVECTOR ,fillprimitivevector )(
      CHFp_VR(primitive)
      ,CHFp_REAL(dense)
      ,CHFp_REAL(press)
      ,CHFp_REALVECT(vel)
      ,CHFp_REAL(entropy)
      ,CHFp_REAL(internal)
      ,CHFp_REAL(cvtemp)
      ,CHFp_REAL(soundspeed) );

#define FORT_FILLPRIMITIVEVECTOR FORTRAN_NAME( inlineFILLPRIMITIVEVECTOR, inlineFILLPRIMITIVEVECTOR)
#define FORTNT_FILLPRIMITIVEVECTOR FORTRAN_NAME( FILLPRIMITIVEVECTOR, fillprimitivevector)

inline void FORTRAN_NAME(inlineFILLPRIMITIVEVECTOR, inlineFILLPRIMITIVEVECTOR)(
      CHFp_VR(primitive)
      ,CHFp_REAL(dense)
      ,CHFp_REAL(press)
      ,CHFp_REALVECT(vel)
      ,CHFp_REAL(entropy)
      ,CHFp_REAL(internal)
      ,CHFp_REAL(cvtemp)
      ,CHFp_REAL(soundspeed) )
{
 CH_TIMELEAF("FORT_FILLPRIMITIVEVECTOR");
 FORTRAN_NAME( FILLPRIMITIVEVECTOR ,fillprimitivevector )(
      CHFt_VR(primitive)
      ,CHFt_REAL(dense)
      ,CHFt_REAL(press)
      ,CHFt_REALVECT(vel)
      ,CHFt_REAL(entropy)
      ,CHFt_REAL(internal)
      ,CHFt_REAL(cvtemp)
      ,CHFt_REAL(soundspeed) );
}
#endif  // GUARDFILLPRIMITIVEVECTOR 

#ifndef GUARDDUMPPRIMITIVEVECTOR 
#define GUARDDUMPPRIMITIVEVECTOR 
// Prototype for Fortran procedure dumpprimitivevector ...
//
void FORTRAN_NAME( DUMPPRIMITIVEVECTOR ,dumpprimitivevector )(
      CHFp_VR(primitive)
      ,CHFp_REAL(dense)
      ,CHFp_REAL(press)
      ,CHFp_REALVECT(vel)
      ,CHFp_REAL(entropy)
      ,CHFp_REAL(internal)
      ,CHFp_REAL(cvtemp)
      ,CHFp_REAL(soundspeed) );

#define FORT_DUMPPRIMITIVEVECTOR FORTRAN_NAME( inlineDUMPPRIMITIVEVECTOR, inlineDUMPPRIMITIVEVECTOR)
#define FORTNT_DUMPPRIMITIVEVECTOR FORTRAN_NAME( DUMPPRIMITIVEVECTOR, dumpprimitivevector)

inline void FORTRAN_NAME(inlineDUMPPRIMITIVEVECTOR, inlineDUMPPRIMITIVEVECTOR)(
      CHFp_VR(primitive)
      ,CHFp_REAL(dense)
      ,CHFp_REAL(press)
      ,CHFp_REALVECT(vel)
      ,CHFp_REAL(entropy)
      ,CHFp_REAL(internal)
      ,CHFp_REAL(cvtemp)
      ,CHFp_REAL(soundspeed) )
{
 CH_TIMELEAF("FORT_DUMPPRIMITIVEVECTOR");
 FORTRAN_NAME( DUMPPRIMITIVEVECTOR ,dumpprimitivevector )(
      CHFt_VR(primitive)
      ,CHFt_REAL(dense)
      ,CHFt_REAL(press)
      ,CHFt_REALVECT(vel)
      ,CHFt_REAL(entropy)
      ,CHFt_REAL(internal)
      ,CHFt_REAL(cvtemp)
      ,CHFt_REAL(soundspeed) );
}
#endif  // GUARDDUMPPRIMITIVEVECTOR 

#ifndef GUARDFILLCONSERVEDVECTOR 
#define GUARDFILLCONSERVEDVECTOR 
// Prototype for Fortran procedure fillconservedvector ...
//
void FORTRAN_NAME( FILLCONSERVEDVECTOR ,fillconservedvector )(
      CHFp_VR(conserved)
      ,CHFp_REAL(dense)
      ,CHFp_REALVECT(momen)
      ,CHFp_REAL(energ) );

#define FORT_FILLCONSERVEDVECTOR FORTRAN_NAME( inlineFILLCONSERVEDVECTOR, inlineFILLCONSERVEDVECTOR)
#define FORTNT_FILLCONSERVEDVECTOR FORTRAN_NAME( FILLCONSERVEDVECTOR, fillconservedvector)

inline void FORTRAN_NAME(inlineFILLCONSERVEDVECTOR, inlineFILLCONSERVEDVECTOR)(
      CHFp_VR(conserved)
      ,CHFp_REAL(dense)
      ,CHFp_REALVECT(momen)
      ,CHFp_REAL(energ) )
{
 CH_TIMELEAF("FORT_FILLCONSERVEDVECTOR");
 FORTRAN_NAME( FILLCONSERVEDVECTOR ,fillconservedvector )(
      CHFt_VR(conserved)
      ,CHFt_REAL(dense)
      ,CHFt_REALVECT(momen)
      ,CHFt_REAL(energ) );
}
#endif  // GUARDFILLCONSERVEDVECTOR 

#ifndef GUARDDUMPCONSERVEDVECTOR 
#define GUARDDUMPCONSERVEDVECTOR 
// Prototype for Fortran procedure dumpconservedvector ...
//
void FORTRAN_NAME( DUMPCONSERVEDVECTOR ,dumpconservedvector )(
      CHFp_VR(conserved)
      ,CHFp_REAL(dense)
      ,CHFp_REALVECT(momen)
      ,CHFp_REAL(energ) );

#define FORT_DUMPCONSERVEDVECTOR FORTRAN_NAME( inlineDUMPCONSERVEDVECTOR, inlineDUMPCONSERVEDVECTOR)
#define FORTNT_DUMPCONSERVEDVECTOR FORTRAN_NAME( DUMPCONSERVEDVECTOR, dumpconservedvector)

inline void FORTRAN_NAME(inlineDUMPCONSERVEDVECTOR, inlineDUMPCONSERVEDVECTOR)(
      CHFp_VR(conserved)
      ,CHFp_REAL(dense)
      ,CHFp_REALVECT(momen)
      ,CHFp_REAL(energ) )
{
 CH_TIMELEAF("FORT_DUMPCONSERVEDVECTOR");
 FORTRAN_NAME( DUMPCONSERVEDVECTOR ,dumpconservedvector )(
      CHFt_VR(conserved)
      ,CHFt_REAL(dense)
      ,CHFt_REALVECT(momen)
      ,CHFt_REAL(energ) );
}
#endif  // GUARDDUMPCONSERVEDVECTOR 

#ifndef GUARDGETSMALL 
#define GUARDGETSMALL 
// Prototype for Fortran procedure getsmall ...
//
void FORTRAN_NAME( GETSMALL ,getsmall )(
      CHFp_REAL(ausmall)
      ,CHFp_REAL(ausmallp)
      ,CHFp_REAL(ausmallu)
      ,CHFp_REAL(ausmallr) );

#define FORT_GETSMALL FORTRAN_NAME( inlineGETSMALL, inlineGETSMALL)
#define FORTNT_GETSMALL FORTRAN_NAME( GETSMALL, getsmall)

inline void FORTRAN_NAME(inlineGETSMALL, inlineGETSMALL)(
      CHFp_REAL(ausmall)
      ,CHFp_REAL(ausmallp)
      ,CHFp_REAL(ausmallu)
      ,CHFp_REAL(ausmallr) )
{
 CH_TIMELEAF("FORT_GETSMALL");
 FORTRAN_NAME( GETSMALL ,getsmall )(
      CHFt_REAL(ausmall)
      ,CHFt_REAL(ausmallp)
      ,CHFt_REAL(ausmallu)
      ,CHFt_REAL(ausmallr) );
}
#endif  // GUARDGETSMALL 

#ifndef GUARDCONS2PRM 
#define GUARDCONS2PRM 
// Prototype for Fortran procedure cons2prm ...
//
void FORTRAN_NAME( CONS2PRM ,cons2prm )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(u)
      ,CHFp_FRA(q)
      ,CHFp_CONST_INT(logflag)
      ,CHFp_CONST_INT(iverbose) );

#define FORT_CONS2PRM FORTRAN_NAME( inlineCONS2PRM, inlineCONS2PRM)
#define FORTNT_CONS2PRM FORTRAN_NAME( CONS2PRM, cons2prm)

inline void FORTRAN_NAME(inlineCONS2PRM, inlineCONS2PRM)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(u)
      ,CHFp_FRA(q)
      ,CHFp_CONST_INT(logflag)
      ,CHFp_CONST_INT(iverbose) )
{
 CH_TIMELEAF("FORT_CONS2PRM");
 FORTRAN_NAME( CONS2PRM ,cons2prm )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(u)
      ,CHFt_FRA(q)
      ,CHFt_CONST_INT(logflag)
      ,CHFt_CONST_INT(iverbose) );
}
#endif  // GUARDCONS2PRM 

#ifndef GUARDPOINTCONS2PRM 
#define GUARDPOINTCONS2PRM 
// Prototype for Fortran procedure pointcons2prm ...
//
void FORTRAN_NAME( POINTCONS2PRM ,pointcons2prm )(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(logflag) );

#define FORT_POINTCONS2PRM FORTRAN_NAME( inlinePOINTCONS2PRM, inlinePOINTCONS2PRM)
#define FORTNT_POINTCONS2PRM FORTRAN_NAME( POINTCONS2PRM, pointcons2prm)

inline void FORTRAN_NAME(inlinePOINTCONS2PRM, inlinePOINTCONS2PRM)(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(logflag) )
{
 CH_TIMELEAF("FORT_POINTCONS2PRM");
 FORTRAN_NAME( POINTCONS2PRM ,pointcons2prm )(
      CHFt_VR(conserved)
      ,CHFt_VR(primitive)
      ,CHFt_CONST_INT(logflag) );
}
#endif  // GUARDPOINTCONS2PRM 

#ifndef GUARDGETTEMPERATURE 
#define GUARDGETTEMPERATURE 
// Prototype for Fortran procedure gettemperature ...
//
void FORTRAN_NAME( GETTEMPERATURE ,gettemperature )(
      CHFp_FRA1(temperature)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(dcalc)
      ,CHFp_CONST_REAL(specificheat) );

#define FORT_GETTEMPERATURE FORTRAN_NAME( inlineGETTEMPERATURE, inlineGETTEMPERATURE)
#define FORTNT_GETTEMPERATURE FORTRAN_NAME( GETTEMPERATURE, gettemperature)

inline void FORTRAN_NAME(inlineGETTEMPERATURE, inlineGETTEMPERATURE)(
      CHFp_FRA1(temperature)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(dcalc)
      ,CHFp_CONST_REAL(specificheat) )
{
 CH_TIMELEAF("FORT_GETTEMPERATURE");
 FORTRAN_NAME( GETTEMPERATURE ,gettemperature )(
      CHFt_FRA1(temperature)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(dcalc)
      ,CHFt_CONST_REAL(specificheat) );
}
#endif  // GUARDGETTEMPERATURE 

#ifndef GUARDPOINTGETTEMPERATURE 
#define GUARDPOINTGETTEMPERATURE 
// Prototype for Fortran procedure pointgettemperature ...
//
void FORTRAN_NAME( POINTGETTEMPERATURE ,pointgettemperature )(
      CHFp_REAL(temper)
      ,CHFp_CONST_VR(conserved)
      ,CHFp_CONST_REAL(specificheat) );

#define FORT_POINTGETTEMPERATURE FORTRAN_NAME( inlinePOINTGETTEMPERATURE, inlinePOINTGETTEMPERATURE)
#define FORTNT_POINTGETTEMPERATURE FORTRAN_NAME( POINTGETTEMPERATURE, pointgettemperature)

inline void FORTRAN_NAME(inlinePOINTGETTEMPERATURE, inlinePOINTGETTEMPERATURE)(
      CHFp_REAL(temper)
      ,CHFp_CONST_VR(conserved)
      ,CHFp_CONST_REAL(specificheat) )
{
 CH_TIMELEAF("FORT_POINTGETTEMPERATURE");
 FORTRAN_NAME( POINTGETTEMPERATURE ,pointgettemperature )(
      CHFt_REAL(temper)
      ,CHFt_CONST_VR(conserved)
      ,CHFt_CONST_REAL(specificheat) );
}
#endif  // GUARDPOINTGETTEMPERATURE 

#ifndef GUARDPOINTGETTEMPFROMPRIM 
#define GUARDPOINTGETTEMPFROMPRIM 
// Prototype for Fortran procedure pointgettempfromprim ...
//
void FORTRAN_NAME( POINTGETTEMPFROMPRIM ,pointgettempfromprim )(
      CHFp_REAL(temper)
      ,CHFp_CONST_VR(primitive)
      ,CHFp_CONST_REAL(specificheat) );

#define FORT_POINTGETTEMPFROMPRIM FORTRAN_NAME( inlinePOINTGETTEMPFROMPRIM, inlinePOINTGETTEMPFROMPRIM)
#define FORTNT_POINTGETTEMPFROMPRIM FORTRAN_NAME( POINTGETTEMPFROMPRIM, pointgettempfromprim)

inline void FORTRAN_NAME(inlinePOINTGETTEMPFROMPRIM, inlinePOINTGETTEMPFROMPRIM)(
      CHFp_REAL(temper)
      ,CHFp_CONST_VR(primitive)
      ,CHFp_CONST_REAL(specificheat) )
{
 CH_TIMELEAF("FORT_POINTGETTEMPFROMPRIM");
 FORTRAN_NAME( POINTGETTEMPFROMPRIM ,pointgettempfromprim )(
      CHFt_REAL(temper)
      ,CHFt_CONST_VR(primitive)
      ,CHFt_CONST_REAL(specificheat) );
}
#endif  // GUARDPOINTGETTEMPFROMPRIM 

#ifndef GUARDGETFLUX 
#define GUARDGETFLUX 
// Prototype for Fortran procedure getflux ...
//
void FORTRAN_NAME( GETFLUX ,getflux )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(qgdnv)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_REAL(dx) );

#define FORT_GETFLUX FORTRAN_NAME( inlineGETFLUX, inlineGETFLUX)
#define FORTNT_GETFLUX FORTRAN_NAME( GETFLUX, getflux)

inline void FORTRAN_NAME(inlineGETFLUX, inlineGETFLUX)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(qgdnv)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_FRA(flux)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_GETFLUX");
 FORTRAN_NAME( GETFLUX ,getflux )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(qgdnv)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_FRA(flux)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDGETFLUX 

#ifndef GUARDPOINTGETFLUX 
#define GUARDPOINTGETFLUX 
// Prototype for Fortran procedure pointgetflux ...
//
void FORTRAN_NAME( POINTGETFLUX ,pointgetflux )(
      CHFp_VR(flux)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(facedir) );

#define FORT_POINTGETFLUX FORTRAN_NAME( inlinePOINTGETFLUX, inlinePOINTGETFLUX)
#define FORTNT_POINTGETFLUX FORTRAN_NAME( POINTGETFLUX, pointgetflux)

inline void FORTRAN_NAME(inlinePOINTGETFLUX, inlinePOINTGETFLUX)(
      CHFp_VR(flux)
      ,CHFp_VR(primitive)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_POINTGETFLUX");
 FORTRAN_NAME( POINTGETFLUX ,pointgetflux )(
      CHFt_VR(flux)
      ,CHFt_VR(primitive)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDPOINTGETFLUX 

#ifndef GUARDPRED 
#define GUARDPRED 
// Prototype for Fortran procedure pred ...
//
void FORTRAN_NAME( PRED ,pred )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) );

#define FORT_PRED FORTRAN_NAME( inlinePRED, inlinePRED)
#define FORTNT_PRED FORTRAN_NAME( PRED, pred)

inline void FORTRAN_NAME(inlinePRED, inlinePRED)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) )
{
 CH_TIMELEAF("FORT_PRED");
 FORTRAN_NAME( PRED ,pred )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(q)
      ,CHFt_CONST_FRA(dq)
      ,CHFt_FRA(qlo)
      ,CHFt_FRA(qhi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(useflat) );
}
#endif  // GUARDPRED 

#ifndef GUARDPOINTLIMIT 
#define GUARDPOINTLIMIT 
// Prototype for Fortran procedure pointlimit ...
//
void FORTRAN_NAME( POINTLIMIT ,pointlimit )(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REALVECT(veloc)
      ,CHFp_REALVECT(vello)
      ,CHFp_REALVECT(velhi)
      ,CHFp_REAL(press)
      ,CHFp_REAL(prelo)
      ,CHFp_REAL(prehi) );

#define FORT_POINTLIMIT FORTRAN_NAME( inlinePOINTLIMIT, inlinePOINTLIMIT)
#define FORTNT_POINTLIMIT FORTRAN_NAME( POINTLIMIT, pointlimit)

inline void FORTRAN_NAME(inlinePOINTLIMIT, inlinePOINTLIMIT)(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REALVECT(veloc)
      ,CHFp_REALVECT(vello)
      ,CHFp_REALVECT(velhi)
      ,CHFp_REAL(press)
      ,CHFp_REAL(prelo)
      ,CHFp_REAL(prehi) )
{
 CH_TIMELEAF("FORT_POINTLIMIT");
 FORTRAN_NAME( POINTLIMIT ,pointlimit )(
      CHFt_REAL(dense)
      ,CHFt_REAL(denlo)
      ,CHFt_REAL(denhi)
      ,CHFt_REALVECT(veloc)
      ,CHFt_REALVECT(vello)
      ,CHFt_REALVECT(velhi)
      ,CHFt_REAL(press)
      ,CHFt_REAL(prelo)
      ,CHFt_REAL(prehi) );
}
#endif  // GUARDPOINTLIMIT 

#ifndef GUARDPOINTPRED 
#define GUARDPOINTPRED 
// Prototype for Fortran procedure pointpred ...
//
void FORTRAN_NAME( POINTPRED ,pointpred )(
      CHFp_VR(primit)
      ,CHFp_VR(pslope)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) );

#define FORT_POINTPRED FORTRAN_NAME( inlinePOINTPRED, inlinePOINTPRED)
#define FORTNT_POINTPRED FORTRAN_NAME( POINTPRED, pointpred)

inline void FORTRAN_NAME(inlinePOINTPRED, inlinePOINTPRED)(
      CHFp_VR(primit)
      ,CHFp_VR(pslope)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(useflat) )
{
 CH_TIMELEAF("FORT_POINTPRED");
 FORTRAN_NAME( POINTPRED ,pointpred )(
      CHFt_VR(primit)
      ,CHFt_VR(pslope)
      ,CHFt_VR(primlo)
      ,CHFt_VR(primhi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(useflat) );
}
#endif  // GUARDPOINTPRED 

#ifndef GUARDPRM2CONS 
#define GUARDPRM2CONS 
// Prototype for Fortran procedure prm2cons ...
//
void FORTRAN_NAME( PRM2CONS ,prm2cons )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(u)
      ,CHFp_CONST_FRA(q) );

#define FORT_PRM2CONS FORTRAN_NAME( inlinePRM2CONS, inlinePRM2CONS)
#define FORTNT_PRM2CONS FORTRAN_NAME( PRM2CONS, prm2cons)

inline void FORTRAN_NAME(inlinePRM2CONS, inlinePRM2CONS)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(u)
      ,CHFp_CONST_FRA(q) )
{
 CH_TIMELEAF("FORT_PRM2CONS");
 FORTRAN_NAME( PRM2CONS ,prm2cons )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(u)
      ,CHFt_CONST_FRA(q) );
}
#endif  // GUARDPRM2CONS 

#ifndef GUARDPOINTPRM2CONS 
#define GUARDPOINTPRM2CONS 
// Prototype for Fortran procedure pointprm2cons ...
//
void FORTRAN_NAME( POINTPRM2CONS ,pointprm2cons )(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive) );

#define FORT_POINTPRM2CONS FORTRAN_NAME( inlinePOINTPRM2CONS, inlinePOINTPRM2CONS)
#define FORTNT_POINTPRM2CONS FORTRAN_NAME( POINTPRM2CONS, pointprm2cons)

inline void FORTRAN_NAME(inlinePOINTPRM2CONS, inlinePOINTPRM2CONS)(
      CHFp_VR(conserved)
      ,CHFp_VR(primitive) )
{
 CH_TIMELEAF("FORT_POINTPRM2CONS");
 FORTRAN_NAME( POINTPRM2CONS ,pointprm2cons )(
      CHFt_VR(conserved)
      ,CHFt_VR(primitive) );
}
#endif  // GUARDPOINTPRM2CONS 

#ifndef GUARDRIEMANN 
#define GUARDRIEMANN 
// Prototype for Fortran procedure riemann ...
//
void FORTRAN_NAME( RIEMANN ,riemann )(
      CHFp_BOX(d)
      ,CHFp_CONST_FRA(ql)
      ,CHFp_CONST_FRA(qr)
      ,CHFp_FRA(qgod)
      ,CHFp_CONST_INT(facedir) );

#define FORT_RIEMANN FORTRAN_NAME( inlineRIEMANN, inlineRIEMANN)
#define FORTNT_RIEMANN FORTRAN_NAME( RIEMANN, riemann)

inline void FORTRAN_NAME(inlineRIEMANN, inlineRIEMANN)(
      CHFp_BOX(d)
      ,CHFp_CONST_FRA(ql)
      ,CHFp_CONST_FRA(qr)
      ,CHFp_FRA(qgod)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_RIEMANN");
 FORTRAN_NAME( RIEMANN ,riemann )(
      CHFt_BOX(d)
      ,CHFt_CONST_FRA(ql)
      ,CHFt_CONST_FRA(qr)
      ,CHFt_FRA(qgod)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDRIEMANN 

#ifndef GUARDOUTFLOWBC 
#define GUARDOUTFLOWBC 
// Prototype for Fortran procedure OUTFLOWBC ...
//
void FORTRAN_NAME( OUTFLOWBC ,outflowbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_OUTFLOWBC FORTRAN_NAME( inlineOUTFLOWBC, inlineOUTFLOWBC)
#define FORTNT_OUTFLOWBC FORTRAN_NAME( OUTFLOWBC, outflowbc)

inline void FORTRAN_NAME(inlineOUTFLOWBC, inlineOUTFLOWBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OUTFLOWBC");
 FORTRAN_NAME( OUTFLOWBC ,outflowbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOUTFLOWBC 

#ifndef GUARDPOINTOUTFLOWBC 
#define GUARDPOINTOUTFLOWBC 
// Prototype for Fortran procedure pointoutflowbc ...
//
void FORTRAN_NAME( POINTOUTFLOWBC ,pointoutflowbc )(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_POINTOUTFLOWBC FORTRAN_NAME( inlinePOINTOUTFLOWBC, inlinePOINTOUTFLOWBC)
#define FORTNT_POINTOUTFLOWBC FORTRAN_NAME( POINTOUTFLOWBC, pointoutflowbc)

inline void FORTRAN_NAME(inlinePOINTOUTFLOWBC, inlinePOINTOUTFLOWBC)(
      CHFp_VR(flux)
      ,CHFp_VR(prim)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_POINTOUTFLOWBC");
 FORTRAN_NAME( POINTOUTFLOWBC ,pointoutflowbc )(
      CHFt_VR(flux)
      ,CHFt_VR(prim)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDPOINTOUTFLOWBC 

#ifndef GUARDSLIPWALLSOLIDBC 
#define GUARDSLIPWALLSOLIDBC 
// Prototype for Fortran procedure SLIPWALLSOLIDBC ...
//
void FORTRAN_NAME( SLIPWALLSOLIDBC ,slipwallsolidbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) );

#define FORT_SLIPWALLSOLIDBC FORTRAN_NAME( inlineSLIPWALLSOLIDBC, inlineSLIPWALLSOLIDBC)
#define FORTNT_SLIPWALLSOLIDBC FORTRAN_NAME( SLIPWALLSOLIDBC, slipwallsolidbc)

inline void FORTRAN_NAME(inlineSLIPWALLSOLIDBC, inlineSLIPWALLSOLIDBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SLIPWALLSOLIDBC");
 FORTRAN_NAME( SLIPWALLSOLIDBC ,slipwallsolidbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSLIPWALLSOLIDBC 

#ifndef GUARDPOINTRIEMANN 
#define GUARDPOINTRIEMANN 
// Prototype for Fortran procedure pointriemann ...
//
void FORTRAN_NAME( POINTRIEMANN ,pointriemann )(
      CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_VR(primgod)
      ,CHFp_CONST_INT(facedir) );

#define FORT_POINTRIEMANN FORTRAN_NAME( inlinePOINTRIEMANN, inlinePOINTRIEMANN)
#define FORTNT_POINTRIEMANN FORTRAN_NAME( POINTRIEMANN, pointriemann)

inline void FORTRAN_NAME(inlinePOINTRIEMANN, inlinePOINTRIEMANN)(
      CHFp_VR(priml)
      ,CHFp_VR(primr)
      ,CHFp_VR(primgod)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_POINTRIEMANN");
 FORTRAN_NAME( POINTRIEMANN ,pointriemann )(
      CHFt_VR(priml)
      ,CHFt_VR(primr)
      ,CHFt_VR(primgod)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDPOINTRIEMANN 

#ifndef GUARDFLOORCONS 
#define GUARDFLOORCONS 
// Prototype for Fortran procedure floorcons ...
//
void FORTRAN_NAME( FLOORCONS ,floorcons )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) );

#define FORT_FLOORCONS FORTRAN_NAME( inlineFLOORCONS, inlineFLOORCONS)
#define FORTNT_FLOORCONS FORTRAN_NAME( FLOORCONS, floorcons)

inline void FORTRAN_NAME(inlineFLOORCONS, inlineFLOORCONS)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) )
{
 CH_TIMELEAF("FORT_FLOORCONS");
 FORTRAN_NAME( FLOORCONS ,floorcons )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state) );
}
#endif  // GUARDFLOORCONS 

#ifndef GUARDFLOORTEMPERATURE 
#define GUARDFLOORTEMPERATURE 
// Prototype for Fortran procedure floortemperature ...
//
void FORTRAN_NAME( FLOORTEMPERATURE ,floortemperature )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA1(temperat) );

#define FORT_FLOORTEMPERATURE FORTRAN_NAME( inlineFLOORTEMPERATURE, inlineFLOORTEMPERATURE)
#define FORTNT_FLOORTEMPERATURE FORTRAN_NAME( FLOORTEMPERATURE, floortemperature)

inline void FORTRAN_NAME(inlineFLOORTEMPERATURE, inlineFLOORTEMPERATURE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA1(temperat) )
{
 CH_TIMELEAF("FORT_FLOORTEMPERATURE");
 FORTRAN_NAME( FLOORTEMPERATURE ,floortemperature )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA1(temperat) );
}
#endif  // GUARDFLOORTEMPERATURE 

#ifndef GUARDFLOORPRIM 
#define GUARDFLOORPRIM 
// Prototype for Fortran procedure floorprim ...
//
void FORTRAN_NAME( FLOORPRIM ,floorprim )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) );

#define FORT_FLOORPRIM FORTRAN_NAME( inlineFLOORPRIM, inlineFLOORPRIM)
#define FORTNT_FLOORPRIM FORTRAN_NAME( FLOORPRIM, floorprim)

inline void FORTRAN_NAME(inlineFLOORPRIM, inlineFLOORPRIM)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state) )
{
 CH_TIMELEAF("FORT_FLOORPRIM");
 FORTRAN_NAME( FLOORPRIM ,floorprim )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state) );
}
#endif  // GUARDFLOORPRIM 

}

#endif

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDVISCOUSTENSORDOMAINBC_H_
#define _MIXEDVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCFuncEval.H"
#include "NamespaceHeader.H"

///
/**
 */
class MixedViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  ///
  /**
   */
  MixedViscousTensorDomainBC()
  {
    ///set default start length to huge so it defaults to Neumann
    m_startDirichlet = 1.23456789e10;
  }

  ///
  /**
   */
  virtual ~MixedViscousTensorDomainBC()
  {;}


  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);
  
  
  void setStartDirichlet(const Real& a_startDirichlet)
  {
    m_startDirichlet = a_startDirichlet;
  }
private:

  Real m_startDirichlet;

};

#include "NamespaceFooter.H"
#endif

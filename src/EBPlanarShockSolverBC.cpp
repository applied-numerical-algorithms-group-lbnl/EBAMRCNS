#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BoxIterator.H"

#include "EBPlanarShockSolverBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletViscousTensorEBBC.H"
#include "NeumannViscousTensorEBBC.H"
#include "MixedViscousTensorDomainBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "EBPlanarShockF_F.H"
#include "SlipWallViscousTensorDomainBC.H"

#include "NamespaceHeader.H"

bool EBPlanarShockSolverBC::s_printed_stuff = false;

class DVTInflowFunc: public  BaseBCFuncEval
{
public:
  DVTInflowFunc(Real a_inflowvel, int a_inflowdir) 
  {
    m_inflowvel = a_inflowvel;
    m_inflowdir = a_inflowdir;
  }


  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const
  {
    Real retval = 0;
    if(a_comp == m_inflowdir) retval = m_inflowvel;

    return retval;
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const 
  {
    return 0.;
  }
protected:

  Real  m_inflowvel;
  int   m_inflowdir;

};



void 
EBPlanarShockSolverBC::
whereAMI(bool& a_atInflow, 
         bool& a_atOutflow,
         bool& a_atMixed,
         const int&            a_idir, 
         const Side::LoHiSide& a_side)
{

  if(a_idir == m_shockNorm)
    {
      a_atOutflow =  (a_side == Side::Hi);
      
      a_atInflow  =  (a_side == Side::Lo); 
    }
  else
    {
      a_atInflow = false;
      a_atOutflow= false;
      a_atMixed = (a_side == Side::Lo);
    }

}

////
void 
EBPlanarShockSolverBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{

  bool atInflow, atOutflow, atMixed;
  whereAMI(atInflow, atOutflow, atMixed, a_idir, a_side);
  if(atInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      Real value = 0;
      FORT_GETPOSTSHOCKVEL(CHF_REAL(value));
      RefCountedPtr<BaseBCFuncEval> funk(new DVTInflowFunc (value, m_shockNorm));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);

     // NeumannViscousTensorDomainBC neumBC;
     // neumBC.setValue(0.);
     //
     // neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
     // neumBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(atMixed)
    {
      MixedViscousTensorDomainBC mixedBC;
      ParmParse pp;
      if(pp.contains("start_dirichlet"))
        {
          Real startDiri;
          pp.get("start_dirichlet", startDiri);
          mixedBC.setStartDirichlet(startDiri);

          if(!s_printed_stuff)
            {
              pout() << "STARTING NO SLIP CONDITION AT X VALUE = "  << startDiri << endl;
              s_printed_stuff = true;
            }
        }
      mixedBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      mixedBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void 
EBPlanarShockSolverBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{

  bool atInflow, atOutflow, atMixed;
  whereAMI(atInflow, atOutflow, atMixed, a_idir, a_side);

  if(atInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      Real inflowvel;
      FORT_GETPOSTSHOCKVEL(CHF_REAL(inflowvel));

      RefCountedPtr<BaseBCFuncEval> funk(new DVTInflowFunc (inflowvel, m_shockNorm));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
//      NeumannViscousTensorDomainBC neumBC;
//      neumBC.setValue(0.);
//      neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
//      neumBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
//                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(atMixed)
    {
      MixedViscousTensorDomainBC mixedBC;
      ParmParse pp;
      if(pp.contains("start_dirichlet"))
        {
          Real startDiri;
          pp.get("start_dirichlet", startDiri);
          mixedBC.setStartDirichlet(startDiri);

          if(!s_printed_stuff)
            {
              pout() << "STARTING NO SLIP CONDITION AT X VALUE = "  << startDiri << endl;
              s_printed_stuff = true;
            }
        }
      mixedBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      mixedBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                          a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

#include "NamespaceFooter.H"

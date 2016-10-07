#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "MultilevelLinearOp.H"
#include "PoissonUtilities.H"
#include "PolyGeom.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBCellFactory.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"

#include "AllRegularService.H"
#include "GeometryShop.H"
#include "BaseIVFactory.H"
#include "PlaneIF.H"
#include "DEMIF.H"
#include "SphereIF.H"
#include "MultiSphereIF.H"
#include "SphereArrayIF.H"
#include "HelicoilIF.H"
#include "PolarIF.H"
#include "LatheIF.H"
#include "TransformIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "PolynomialIF.H"
#include "TylerChannelIF.H"
#include "SlabService.H"
#include "MultiSlabService.H"
#include "EllipsoidIF.H"
#include "EBSimpleSolver.H"
#include "DataFileIF.H"
#include "CH_Timer.H"
#include "CH_Attach.H"
#include "EBViscousTensorOp.H"
#include "EBViscousTensorOpFactory.H"
#include "NWOEBViscousTensorOp.H"
#include "NWOEBViscousTensorOpFactory.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"

#include "memusage.H"
#include "memtrack.H"


/*****/
/*****/
void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar)
{
  a_boxesCoar.resize(a_boxesFine.size());
  for (int ilev = 0; ilev < a_boxesFine.size(); ilev++)
    {
      a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
      for (int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
        {
          a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
        }
    }
}

void
getSimpleBoxes(Vector<Vector<Box> >&   a_boxes,
               Vector<int>&            a_refRat,
               const Box&              a_domain)
{
  ParmParse pp;
  int maxlev;
  pp.get("max_level", maxlev);
  int amrRef = 2;

  a_refRat.resize(maxlev+1, amrRef);
  a_boxes.resize( maxlev+1);
  Box domlev = a_domain;
  int ishrink = domlev.size(0);
  ishrink /= 2;
  a_boxes[0] = Vector<Box>(1, a_domain);
  for (int ilev = 1 ; ilev <= maxlev; ilev++)
    {
      Box fineBox = refine(domlev, amrRef);
      fineBox.grow(-ishrink);

      a_boxes[ilev] = Vector<Box>(1, fineBox);

      domlev.refine(amrRef);
      ishrink  += 8;
    }
}
void
defineGrids(Vector< DisjointBoxLayout     >& a_grids,
            Vector< EBISLayout            >& a_ebisl,
            const Vector<Vector<Box>      >& a_boxes,
            const Vector<int              >& a_refRat,
            const PoissonParameters        & a_params)
{
  a_grids.resize(a_boxes.size());
  a_ebisl.resize(a_boxes.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
  for (int ilev = 0; ilev < a_boxes.size(); ilev++)
    {
      Vector<int> procs;
      EBEllipticLoadBalance(procs, a_boxes[ilev], domLev);
      a_grids[ilev] = DisjointBoxLayout(a_boxes[ilev], procs);
      int numGhost = 4;
      ebisPtr->fillEBISLayout(a_ebisl[ilev],
                              a_grids[ilev],
                              domLev,
                              numGhost);
      domLev.refine(a_refRat[ilev]);
    }
}

void getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_domBC,
                                RefCountedPtr<BaseEBBCFactory>&                         a_ebBC,
                                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_bcoef,
                                const Vector<DisjointBoxLayout>&                        a_grids,
                                const Vector<EBISLayout>&                               a_ebisl,
                                const PoissonParameters&                                a_params)
{
  CH_TIME("PoissonUtilities::getBCFactories");
  ParmParse pp;
  if (a_params.domBcType == 0)
    {
      pout() << "constant Neumann bcs on domain" << endl;

      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);

      NeumannConductivityDomainBCFactory* domainBCFactory = new NeumannConductivityDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 1)
    {
      pout() << "constant Dirichlet bcs on domain" << endl;
      DirichletConductivityDomainBCFactory* domainBCFactory = new DirichletConductivityDomainBCFactory();
      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 2)
    {
      pout() << "Neumann trig bcs on domain" << endl;
      RealVect     trig = getTrigRV();

      TrigBCFlux* trigFluxPtr = new TrigBCFlux();

      trigFluxPtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigFluxPtr = RefCountedPtr<BaseBCValue>(trigFluxPtr);

      NeumannConductivityDomainBCFactory* domainBCFactory = new NeumannConductivityDomainBCFactory();
      domainBCFactory->setFunction(baseTrigFluxPtr);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 3)
    {
      pout() << "Dirichlet trig bcs on domain" << endl;
      RealVect     trig = getTrigRV();

      TrigBCValue* trigValuePtr = new TrigBCValue();
      trigValuePtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigValuePtr = RefCountedPtr<BaseBCValue>(trigValuePtr);

      DirichletConductivityDomainBCFactory* domainBCFactory = new DirichletConductivityDomainBCFactory();
      domainBCFactory->setFunction(baseTrigValuePtr);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else
    {
      MayDay::Error("Unknown domain BC type");
    }


  if (a_params.ebBcType == 0)
    {
      pout() << "constant Neumann bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      NeumannConductivityEBBCFactory* ebBC = new NeumannConductivityEBBCFactory();
      ebBC->setValue(domBCValue);

      a_ebBC =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (a_params.ebBcType == 1)
    {
      int orderEB = getOrderEB();
      pout() << "constant Dirichlet bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      DirichletConductivityEBBCFactory* ebBC = new DirichletConductivityEBBCFactory();
      ebBC->setOrder(orderEB);
      ebBC->setValue(domBCValue);

      a_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else if (a_params.ebBcType == 2)
    {
      pout() << "Neumann trig bcs on EB" << endl;
      RealVect     trig = getTrigRV();
      TrigBCFlux* trigFluxPtr = new TrigBCFlux();
      trigFluxPtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigFluxPtr =  RefCountedPtr<BaseBCValue>(trigFluxPtr);

      NeumannConductivityEBBCFactory* ebBC = new NeumannConductivityEBBCFactory();
      ebBC->setFunction(baseTrigFluxPtr);

      a_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else if (a_params.ebBcType == 3)
    {
      pout() << "Dirichlet trig bcs on EB" << endl;
      RealVect     trig = getTrigRV();
      RefCountedPtr<BaseBCValue> baseTrigValuePtr;
      DirichletConductivityEBBCFactory* ebBC;
      TrigBCValue* trigValuePtr = new TrigBCValue();
      trigValuePtr->define(trig);
      baseTrigValuePtr =  RefCountedPtr<BaseBCValue>(trigValuePtr);
      ebBC = new DirichletConductivityEBBCFactory();

      int orderEB = getOrderEB();
      ebBC->setOrder(orderEB);
      ebBC->setFunction(baseTrigValuePtr);

      a_ebBC = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else
    {
      MayDay::Error("Unknown EB BC type");
    }
}
/**/
void
defineConductivitySolver( AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          const Vector<EBISLayout>&                    a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params)
{
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    bcoIrreg;
  RefCountedPtr<EBConductivityOpFactory> opFactory;
  defineConductivitySolver(a_solver, a_grids, a_ebisl, a_bottomSolver, a_params,
                           aco, bco, bcoIrreg, opFactory);
}

void
defineConductivitySolver( AMRMultiGrid<LevelData<EBCellFAB> >&                     a_solver,
                          const Vector<DisjointBoxLayout>&                         a_grids,
                          const Vector<EBISLayout>&                                a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&                     a_bottomSolver,
                          const PoissonParameters&                                 a_params,
                          Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           aco,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           bco,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    bcoIrreg,
                          RefCountedPtr<EBConductivityOpFactory>&                  opFactory)
{
  defineConductivityCoef(           aco, bco, bcoIrreg, a_grids, a_ebisl, a_params);
  getConductivityFactory(opFactory, aco, bco, bcoIrreg, a_grids, a_ebisl, a_params);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, *opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ParmParse pp2;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
#ifdef CH_USE_FLOAT
  eps = sqrt(eps);
#endif
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 5;

}


void
defineSlowConductivitySolver( AMRMultiGrid<LevelData<EBCellFAB> >&                     a_solver,
                              const Vector<DisjointBoxLayout>&                         a_grids,
                              const Vector<EBISLayout>&                                a_ebisl,
                              LinearSolver<LevelData<EBCellFAB> >&                     a_bottomSolver,
                              const PoissonParameters&                                 a_params,
                              Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           aco,
                              Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           bco,
                              Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    bcoIrreg,
                              RefCountedPtr<slowEBCOFactory>&                          opFactory)
{
  defineConductivityCoef(           aco, bco, bcoIrreg, a_grids, a_ebisl, a_params);
  getSlowConductivityFactory(opFactory, aco, bco, bcoIrreg, a_grids, a_ebisl, a_params);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, *opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ParmParse pp2;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
#ifdef CH_USE_FLOAT
  eps = sqrt(eps);
#endif
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 5;

}
/**/
Real getDistanceFromAxis(const RealVect& a_xloc,
                         const RealVect& a_axis,
                         const RealVect& a_orig)
{
  RealVect trueNorm, closestPt;
  PolyGeom::pointToLine(closestPt, trueNorm, a_xloc,
                        a_orig, a_axis);
  RealVect distVec = a_xloc - closestPt;
  Real distance;
  PolyGeom::unifyVector(distVec, distance);
  return distance;
}
/**/
void setConductivityCoefs(LevelData<EBCellFAB>          &    a_aco,
                          LevelData<EBFluxFAB>          &    a_bco,
                          LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
                          const Real                    &    a_dx,
                          const  PoissonParameters      &    a_params)
{
  int whichAco, whichBco;
  ParmParse pp;
  pp.get("which_acoef", whichAco);
  pp.get("which_bcoef", whichBco);
  if ((whichAco == 0) && (whichBco == 0))
    {
      Real acoVal, bcoVal;
      pp.get("acoef_value", acoVal);
      pp.get("bcoef_value", bcoVal);
      pout() << "constant acoef = " << acoVal << ", constant bcoef = " << bcoVal << endl;
      for (DataIterator dit = a_aco.dataIterator(); dit.ok(); ++dit)
        {
          a_aco     [dit()].setVal(acoVal);
          a_bco     [dit()].setVal(bcoVal);
          a_bcoIrreg[dit()].setVal(bcoVal);
        }
    }
  else
    {
      MayDay::Error("setConductivityCoef: we have limited acoef and bcoef options available at this time");
    }
  a_aco.exchange(Interval(0,0));
  a_bco.exchange(Interval(0,0));
  a_bcoIrreg.exchange(Interval(0,0));

}
/**/
void
defineConductivityCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                       Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                       Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                       const Vector<DisjointBoxLayout>&                         a_grids,
                       const Vector<EBISLayout>&                                a_ebisl,
                       const PoissonParameters&                                 a_params)
{
  CH_TIME("PoissonUtilities::defineConductivityCoef");
  a_aco.resize(        a_params.numLevels);
  a_bco.resize(        a_params.numLevels);
  a_bcoIrreg.resize(   a_params.numLevels);
  Real dxLev = a_params.coarsestDx[0];
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      LayoutData<IntVectSet> irregSets(a_grids[ilev]);
      int nghost = 1;
      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //has to correspond to number of ghost cells
          Box grownBox = grow(a_grids[ilev].get(dit()), nghost);
          grownBox &= domLev;
          irregSets[dit()] = a_ebisl[ilev][dit()].getIrregIVS(grownBox);
        }
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], irregSets);

      a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebcellfact));
      a_bco[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
      a_bcoIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(a_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));
      setConductivityCoefs(*a_aco[ilev], *a_bco[ilev], *a_bcoIrreg[ilev], dxLev, a_params);
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
  }
}
/**/
void
getConductivityFactory(RefCountedPtr<EBConductivityOpFactory>   &                     a_factory,
                       const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                       const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                       const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                       const Vector<DisjointBoxLayout>&                               a_grids,
                       const Vector<EBISLayout>&                                      a_ebisl,
                       const PoissonParameters&                                       a_params)
{
  CH_TIME("PoissonUtilities::getEBVTOFactory");
  ParmParse pp2;
  Real alpha, beta;
  pp2.get("alpha", alpha);
  pp2.get("beta", beta);

  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  ProblemDomain levDom =  a_params.coarsestDomain;
  ProblemDomain coarDom;
  Vector<EBLevelGrid> eblg(a_grids.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      if (ilev > 0)
        {
          int nref = a_params.refRatio[ilev-1];
          int nvar = 1;
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           coarDom,
                                                                           nref, nvar,
                                                                           *eblg[ilev].getCFIVS()));
          coarDom.refine(a_params.refRatio[ilev]);
        }
      coarDom = levDom;
      levDom.refine(a_params.refRatio[ilev]);
    }

  RefCountedPtr<BaseDomainBCFactory>      domBC;
  RefCountedPtr<BaseEBBCFactory>          ebBC;
  getConductivityBCFactories(domBC, ebBC, a_bco, a_grids, a_ebisl, a_params);

  int relaxType = 1;
  ParmParse pp;
  pp.query("relax_type", relaxType);
  if (relaxType == 1)
    {
      pout() << "using multicolor gauss-seidel relaxation" << endl;
    }
  else if (relaxType == 2)
    {
      pout() << "using gauss-seidel relaxation with red-black ordering" << endl;
    }
  else
    {
      MayDay::Error("invalid relaxation type");
    }
  a_factory = RefCountedPtr<EBConductivityOpFactory>
    (new EBConductivityOpFactory(eblg, quadCFI, alpha, beta, a_aco, a_bco, a_bcoIrreg,
                                 a_params.coarsestDx[0],  a_params.refRatio, domBC, ebBC,
                                 a_params.ghostPhi, a_params.ghostRHS, relaxType));
}

/**/
void
getSlowConductivityFactory(RefCountedPtr<slowEBCOFactory>           &                     a_factory,
                           const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                           const Vector<DisjointBoxLayout>&                               a_grids,
                           const Vector<EBISLayout>&                                      a_ebisl,
                           const PoissonParameters&                                       a_params)
{
  CH_TIME("PoissonUtilities::getEBVTOFactory");
  ParmParse pp2;
  Real alpha, beta;
  pp2.get("alpha", alpha);
  pp2.get("beta", beta);

  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  ProblemDomain levDom =  a_params.coarsestDomain;
  ProblemDomain coarDom;
  Vector<EBLevelGrid> eblg(a_grids.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      if (ilev > 0)
        {
          int nref = a_params.refRatio[ilev-1];
          int nvar = 1;
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           coarDom,
                                                                           nref, nvar,
                                                                           *eblg[ilev].getCFIVS()));
          coarDom.refine(a_params.refRatio[ilev]);
        }
      coarDom = levDom;
      levDom.refine(a_params.refRatio[ilev]);
    }

  RefCountedPtr<BaseDomainBCFactory>      domBC;
  RefCountedPtr<BaseEBBCFactory>          ebBC;
  getConductivityBCFactories(domBC, ebBC, a_bco, a_grids, a_ebisl, a_params);

  int relaxType;
  ParmParse pp;
  pp.get("mg_relax_type", relaxType);
  a_factory = RefCountedPtr<slowEBCOFactory>
    (new slowEBCOFactory(eblg, quadCFI, alpha, beta, a_aco, a_bco, a_bcoIrreg,
                         a_params.coarsestDx[0],  a_params.refRatio, domBC, ebBC,
                         a_params.ghostPhi, a_params.ghostRHS, relaxType));
}
/**/
void
getEBVTOFactory(RefCountedPtr<EBViscousTensorOpFactory>&                       a_factory,
                const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_eta,
                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_lambda,
                const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_etaIrreg,
                const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_lambdaIrreg,
                const Vector<DisjointBoxLayout>&                               a_grids,
                const Vector<EBISLayout>&                                      a_ebisl,
                const PoissonParameters&                                       a_params)
{
  CH_TIME("PoissonUtilities::getEBVTOFactory");
  ParmParse pp2;
  Real alpha, beta;
  pp2.get("alpha", alpha);
  pp2.get("beta", beta);

  Vector<EBLevelGrid> eblg(a_grids.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {

      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }

  RefCountedPtr<BaseDomainBCFactory>      domBC;
  RefCountedPtr<BaseEBBCFactory>          ebBC;
  getViscousBCFactories(domBC, ebBC, a_eta, a_lambda,  a_grids, a_ebisl, a_params);

  a_factory = RefCountedPtr<EBViscousTensorOpFactory>
    (new EBViscousTensorOpFactory(eblg, alpha, beta, a_aco, a_eta, a_lambda, a_etaIrreg, a_lambdaIrreg,
                                  a_params.coarsestDx[0],  a_params.refRatio, domBC, ebBC,
                                  a_params.ghostPhi, a_params.ghostRHS, -1, true));
}

/**/
void
getNWOEBVTOFactory(RefCountedPtr<NWOEBViscousTensorOpFactory>&                    a_factory,
                   const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                   const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_eta,
                   const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_lambda,
                   const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_etaIrreg,
                   const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_lambdaIrreg,
                   const Vector<DisjointBoxLayout>&                               a_grids,
                   const Vector<EBISLayout>&                                      a_ebisl,
                   const PoissonParameters&                                       a_params)
{
  CH_TIME("PoissonUtilities::getEBVTOFactory");
  ParmParse pp2;
  Real alpha, beta;
  pp2.get("alpha", alpha);
  pp2.get("beta", beta);

  Vector<EBLevelGrid> eblg(a_grids.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {

      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }

  RefCountedPtr<BaseDomainBCFactory>      domBC;
  RefCountedPtr<BaseEBBCFactory>          ebBC;
  getViscousBCFactories(domBC, ebBC, a_eta, a_lambda,  a_grids, a_ebisl, a_params);

  a_factory = RefCountedPtr<NWOEBViscousTensorOpFactory>
    (new NWOEBViscousTensorOpFactory(eblg, alpha, beta, a_aco, a_eta, a_lambda, a_etaIrreg, a_lambdaIrreg,
                                     a_params.coarsestDx[0],  a_params.refRatio, domBC, ebBC,
                                     a_params.ghostPhi, a_params.ghostRHS, -1, true));
}
/**/
void
defineViscousTensorCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                        Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_eta,
                        Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_lambda,
                        Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_etaIrreg,
                        Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_lambdaIrreg,
                        const Vector<DisjointBoxLayout>&                         a_grids,
                        const Vector<EBISLayout>&                                a_ebisl,
                        const PoissonParameters&                                 a_params)
{
  CH_TIME("PoissonUtilities::defineViscousTensorCoef");
  a_aco.resize(        a_params.numLevels);
  a_eta.resize(        a_params.numLevels);
  a_lambda.resize(     a_params.numLevels);
  a_etaIrreg.resize(   a_params.numLevels);
  a_lambdaIrreg.resize(a_params.numLevels);
  Real dxLev = a_params.coarsestDx[0];
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      LayoutData<IntVectSet> irregSets(a_grids[ilev]);
      int nghost = 4;
      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //has to correspond to number of ghost cells
          Box grownBox = grow(a_grids[ilev].get(dit()), nghost);
          grownBox &= domLev;
          irregSets[dit()] = a_ebisl[ilev][dit()].getIrregIVS(grownBox);
        }
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], irregSets);

      a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebcellfact));
      a_eta[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
      a_lambda[ilev]      = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(a_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
      a_etaIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(a_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));
      a_lambdaIrreg[ilev] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(a_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));
      setViscousCoefs(*a_aco[ilev], *a_eta[ilev], *a_lambda[ilev], *a_etaIrreg[ilev], *a_lambdaIrreg[ilev], dxLev, a_params);

      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
  }
}
void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          const Vector<EBISLayout>&                    a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params)
{
  RefCountedPtr<EBViscousTensorOpFactory> opFactory;
  defineViscousTensorSolver(a_solver, opFactory, a_grids, a_ebisl, a_bottomSolver, a_params);
 }

void
defineNWOViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                             const Vector<DisjointBoxLayout>&             a_grids,
                             const Vector<EBISLayout>&                    a_ebisl,
                             LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                             const PoissonParameters&                     a_params)
{
  RefCountedPtr<NWOEBViscousTensorOpFactory> opFactory;
  defineNWOViscousTensorSolver(a_solver, opFactory, a_grids, a_ebisl, a_bottomSolver, a_params);
 }

void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                          RefCountedPtr<EBViscousTensorOpFactory>&     a_opFactory,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          const Vector<EBISLayout>&                    a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params)
{
  CH_TIME("PoissonUtilities::defineViscousTensorSolver");
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    lambdaIrreg;

  defineViscousTensorSolver(a_solver,
                            a_grids,
                            a_ebisl,
                            a_bottomSolver,
                            a_params,
                            aco,
                            eta,
                            lambda,
                            etaIrreg,
                            lambdaIrreg,
                            a_opFactory);
}

void
defineNWOViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                             RefCountedPtr<NWOEBViscousTensorOpFactory>&  a_opFactory,
                             const Vector<DisjointBoxLayout>&             a_grids,
                             const Vector<EBISLayout>&                    a_ebisl,
                             LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                             const PoissonParameters&                     a_params)
{
  CH_TIME("PoissonUtilities::defineViscousTensorSolver");
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    lambdaIrreg;

  defineNWOViscousTensorSolver(a_solver,
                               a_grids,
                               a_ebisl,
                               a_bottomSolver,
                               a_params,
                               aco,
                               eta,
                               lambda,
                               etaIrreg,
                               lambdaIrreg,
                               a_opFactory);
}


void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&                    a_solver,
                          const Vector<DisjointBoxLayout>&                        a_grids,
                          const Vector<EBISLayout>&                               a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&                    a_bottomSolver,
                          const PoissonParameters&                                a_params,
                          Vector<RefCountedPtr<LevelData<EBCellFAB> > >       &   aco,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   eta,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   lambda,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   etaIrreg,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   lambdaIrreg,
                          RefCountedPtr<EBViscousTensorOpFactory>&                a_opFactory)


{

  defineViscousTensorCoef(aco, eta, lambda, etaIrreg, lambdaIrreg,  a_grids, a_ebisl, a_params);

  getEBVTOFactory(a_opFactory, aco, eta, lambda, etaIrreg, lambdaIrreg, a_grids, a_ebisl, a_params);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, *a_opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ParmParse pp2;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
#ifdef CH_USE_FLOAT
  eps = sqrt(eps);
#endif
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 5;

}

void
defineNWOViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&                    a_solver,
                             const Vector<DisjointBoxLayout>&                        a_grids,
                             const Vector<EBISLayout>&                               a_ebisl,
                             LinearSolver<LevelData<EBCellFAB> >&                    a_bottomSolver,
                             const PoissonParameters&                                a_params,
                             Vector<RefCountedPtr<LevelData<EBCellFAB> > >       &   aco,
                             Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   eta,
                             Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   lambda,
                             Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   etaIrreg,
                             Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   lambdaIrreg,
                             RefCountedPtr<NWOEBViscousTensorOpFactory>&             a_opFactory)


{

  defineViscousTensorCoef(aco, eta, lambda, etaIrreg, lambdaIrreg,  a_grids, a_ebisl, a_params);

  getNWOEBVTOFactory(a_opFactory, aco, eta, lambda, etaIrreg, lambdaIrreg, a_grids, a_ebisl, a_params);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, *a_opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ParmParse pp2;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
#ifdef CH_USE_FLOAT
  eps = sqrt(eps);
#endif
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 5;

}
void setVelViscous(LevelData<EBCellFAB>&    a_mag,
                   const Real&              a_dx,
                   const PoissonParameters& a_params)
{
  CH_TIME("PoissonUtilities::setVelViscous");
  ParmParse pp;
  bool constant_coef;
  pp.get("use_constant_coef", constant_coef);
  if (constant_coef)
    {
      Real time = 0;
      RealVect vectDx = a_dx*RealVect::Unit;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          setTrigPhi(a_mag, vectDx, a_params, time, idir);
        }
    }
  else
    {
      int whichVel;
      pp.get("which_vel", whichVel);
      //for now just pass through.   might want to do something more than this tho
      setMagResistive(a_mag, a_dx, a_params, whichVel);
    }
}

void setKLVViscous(LevelData<EBCellFAB>&    a_klv,
                   const Real&              a_dx,
                   const PoissonParameters& a_params)
{
  CH_TIME("PoissonUtilities::setKLVViscous");
  int whichVel, whichLambda, whichEta;
  Real eps, alpha, beta;
  ParmParse pp;
  bool constant_coef;
  pp.get("use_constant_coef", constant_coef);
  if (constant_coef)
    {
      Real acoefval, etaval, lambdaval;
      pp.get("acoef_value",  acoefval);
      pp.get("eta_value",    etaval);
      pp.get("lambda_value", lambdaval);
      bool lambda_only;
      pp.get("lambda_only", lambda_only);
      pp.get("alpha",   alpha);
      pp.get("beta",   beta);
      Real alphaeff, betaeff;
      alphaeff =  acoefval*alpha;
      if (lambda_only)
        {
          etaval = 0;
          betaeff = lambdaval*beta;
        }
      else
        {
          lambdaval = -etaval;
          betaeff =    etaval*beta;
        }
      Real time = 0;
      RealVect vectDx = a_dx*RealVect::Unit;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          setTrigKappaLOfPhi(a_klv, vectDx, a_params, alphaeff, betaeff, time, idir);
        }
    }
  else
    {
      pp.get("which_vel", whichVel);
      pp.get("which_eta", whichEta);
      pp.get("which_lambda", whichLambda);
      pp.get("eta_eps",   eps);
      pp.get("alpha",   alpha);
      pp.get("beta",   beta);

      //which_lambda == 2 forces lambda = -factor*eta
      Real factor = 1;
      if (whichLambda == 2)
        {
          pp.get("lambda_factor", factor);
        }
      for (DataIterator dit = a_klv.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<Real>& curFAB = a_klv[dit()].getSingleValuedFAB();
          Box curBox = curFAB.box();

          const RealVect&     trig = getTrigRV();

          RealVect dxv = a_dx*RealVect::Unit;
          for (int comp = 0; comp < a_klv.nComp(); comp++)
            {
              FORT_GETKLVVISCOUS(CHF_FRA1(curFAB,comp),
                                 CHF_CONST_REALVECT(trig),
                                 CHF_CONST_REALVECT(dxv),
                                 CHF_CONST_REALVECT(a_params.probLo),
                                 CHF_CONST_REAL(alpha),
                                 CHF_BOX(curBox),
                                 CHF_INT(comp),
                                 CHF_REAL(eps),
                                 CHF_INT(whichVel),
                                 CHF_INT(whichEta),
                                 CHF_INT(whichLambda),
                                 CHF_CONST_REAL(beta),
                                 CHF_REAL(factor)
                                 );
            }
          const EBISBox& ebisBox = a_klv[dit()].getEBISBox();
          IntVectSet ivsMulti = ebisBox.getIrregIVS(curBox);
          for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              //the divergence of a flux lives at the centroid
              RealVect effProbLo = ebisBox.centroid(vofit());
              effProbLo *= a_dx;
              RealVect xloc = EBArith::getVofLocation(vofit(), dxv, effProbLo);
              Real  volFrac = ebisBox.volFrac(vofit());
              for (int comp = 0; comp < a_klv.nComp(); comp++)
                {
                  Real klvpt;
                  FORT_GETKLVPOINTVISCOUS(CHF_REAL(klvpt),
                                          CHF_CONST_REALVECT(trig),
                                          CHF_CONST_REALVECT(xloc),
                                          CHF_CONST_REAL(alpha),
                                          CHF_INT(comp),
                                          CHF_REAL(eps),
                                          CHF_INT(whichVel),
                                          CHF_INT(whichEta),
                                          CHF_INT(whichLambda),
                                          CHF_CONST_REAL(beta),
                                          CHF_REAL(factor));

                  a_klv[dit()](vofit(), comp) = volFrac*klvpt;
                }
            }
        }
    }
}

void setViscousCoefs(LevelData<EBCellFAB>          &    a_aco,
                     LevelData<EBFluxFAB>          &    a_eta,
                     LevelData<EBFluxFAB>          &    a_lambda,
                     LevelData<BaseIVFAB<Real> >   &    a_etaIrreg,
                     LevelData<BaseIVFAB<Real> >   &    a_lambdaIrreg,
                     const Real                    &    a_dx,
                     const  PoissonParameters      &    a_params)
{
  CH_TIME("PoissonUtilities::setViscousCoefs");
  int whichVel, whichLambda, whichEta;
  bool constant_coef;
  Real eps;
  ParmParse pp;
  pp.get("which_vel", whichVel);
  pp.get("which_eta", whichEta);
  pp.get("which_lambda", whichLambda);
  pp.get("eta_eps",   eps);
  pp.get("use_constant_coef", constant_coef);
  if (constant_coef)
    {
      Real acoefval, etaval, lambdaval;
      pp.get("acoef_value",  acoefval);
      pp.get("eta_value",    etaval);
      pp.get("lambda_value",lambdaval);
      bool lambda_only;
      pp.get("lambda_only", lambda_only);
      if (lambda_only)
        {
          etaval = 0;
        }

      pout() << "setting acoef  to " << acoefval  << " everywhere."<< endl;
      pout() << "setting eta    to " << etaval <<    " everywhere."<< endl;
      pout() << "setting lambda to " << lambdaval << " everywhere."<< endl;
      for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
        {
          a_aco        [dit()].setVal(acoefval);
          a_lambda     [dit()].setVal(lambdaval);
          a_lambdaIrreg[dit()].setVal(lambdaval);
          a_eta        [dit()].setVal(etaval);
          a_etaIrreg   [dit()].setVal(etaval);
        }
    }
  else
    {
      setEtaResistive(a_eta, a_etaIrreg, a_dx, a_params, whichEta);
      pout() << "setting acoef to 1 everywhere --- setViscousCoefs needs to have more options" << endl;
      for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
        {
          a_aco[dit()].setVal(1.0);
        }
      if (whichLambda != 2)
        {
          setEtaResistive(a_lambda, a_lambdaIrreg, a_dx, a_params, whichEta);
        }
      else
        {
          //whichlambda == 2 forces lambda = -factor*eta
          Real factor;
          pp.get("lambda_factor", factor);
          for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
            {
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  a_lambda[dit()][idir].copy(a_eta[dit()][idir]);
                  a_lambda[dit()][idir] *= -factor;
                }
              const IntVectSet& ivs  =  a_lambdaIrreg[dit()].getIVS();
              const EBGraph& ebgraph =  a_lambdaIrreg[dit()].getEBGraph();
              for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
                {
                  a_lambdaIrreg[dit()](vofit(), 0) = -factor*a_etaIrreg[dit()](vofit(),0);
                }
            }
        }
    }
  a_aco        .exchange(Interval(0,0));
  a_eta        .exchange(Interval(0,0));
  a_lambda     .exchange(Interval(0,0));
  a_etaIrreg   .exchange(Interval(0,0));
  a_lambdaIrreg.exchange(Interval(0,0));
}
/********/
void setEtaResistive(LevelData<EBFluxFAB>&        a_eta,
                     LevelData<BaseIVFAB<Real> >& a_etaIrreg,
                     const Real&                  a_dx,
                     const PoissonParameters&     a_params,
                     int a_whichEta)
{
  CH_TIME("PoissonUtilities::setEtaResistive");
  CH_assert(a_eta.nComp() == 1);

  int whichEta;
  ParmParse pp;
  if (a_whichEta < 0)
    {
      pp.get("which_eta", whichEta);
    }
  else
    {
      whichEta = a_whichEta;
    }
  Real eps;
  pp.get("eta_eps", eps);

  for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
    {

      if (whichEta == 0)
        {
          a_etaIrreg[dit()].setVal(1.0);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_eta[dit()][idir].setVal(1.0);
            }
        }
      else
        {
          //          Box cellBox = a_eta[dit()].getRegion();
          Box cellBox = a_eta.disjointBoxLayout().get(dit());
          const EBISBox& ebisBox = a_eta[dit()].getEBISBox();
          IntVectSet ivsIrreg = ebisBox.getIrregIVS(cellBox);
          IntVectSet ivsMulti = ebisBox.getMultiCells(cellBox);
          RealVect dxv = a_dx*RealVect::Unit;
          const RealVect&     trig = getTrigRV();
          Real etaval = 1;
          pp.query("eta_value",    etaval);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              BaseFab<Real>& curFAB = a_eta[dit()][idir].getSingleValuedFAB();
              Box curBox = curFAB.box();

              FORT_GETETARESIST(CHF_FRA1(curFAB, 0),
                                CHF_CONST_REALVECT(trig),
                                CHF_CONST_REALVECT(dxv),
                                CHF_CONST_REALVECT(a_params.probLo),
                                CHF_BOX(curBox),
                                CHF_INT(idir),
                                CHF_REAL(etaval),
                                CHF_REAL(eps),
                                CHF_INT(whichEta));


              for (FaceIterator faceit(ivsMulti, ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
                {
                  RealVect xloc     = EBArith::getFaceLocation(faceit(), dxv, a_params.probLo);
                  for (int comp = 0; comp < a_eta.nComp(); comp++)
                    {
                      Real etaPt;
                      FORT_GETETAPOINTRESIST(CHF_REAL(etaPt),
                                             CHF_CONST_REALVECT(trig),
                                             CHF_CONST_REALVECT(xloc),
                                             CHF_REAL(etaval),
                                             CHF_REAL(eps),
                                             CHF_INT(whichEta));

                      a_eta[dit()][idir](faceit(), comp) = etaPt;
                    }
                }
            }

          for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              RealVect bndryCentroid = ebisBox.bndryCentroid(vofit());
              bndryCentroid *= a_dx;
              bndryCentroid += a_params.probLo;
              RealVect bndryLoc = EBArith::getVofLocation(vofit(), dxv, bndryCentroid);
              for (int comp = 0; comp < a_eta.nComp(); comp++)
                {
                  Real etaBndry;

                  FORT_GETETAPOINTRESIST(CHF_REAL(etaBndry),
                                         CHF_CONST_REALVECT(trig),
                                         CHF_CONST_REALVECT(bndryLoc),
                                         CHF_REAL(etaval),
                                         CHF_REAL(eps),
                                         CHF_INT(whichEta));
                  a_etaIrreg[dit()](vofit(), comp) = etaBndry;
                }
            }
        }
    }
}
/********/
void setMagResistive(LevelData<EBCellFAB>&    a_mag,
  const Real&              a_dx,
  const PoissonParameters& a_params,
  int a_whichMag)
{
  CH_TIME("PoissonUtilities::setMagResistive");
  CH_assert(a_mag.nComp() ==SpaceDim);

  int  whichMag;
  if (a_whichMag < 0)
    {
      ParmParse pp;
      pp.get("which_mag", whichMag);
    }
  else
    {
      whichMag = a_whichMag;
    }
  for (DataIterator dit = a_mag.dataIterator(); dit.ok(); ++dit)
    {
      BaseFab<Real>& curFAB = a_mag[dit()].getSingleValuedFAB();
      Box curBox = curFAB.box();

      const RealVect&     trig = getTrigRV();

      RealVect dxv = a_dx*RealVect::Unit;
      for (int comp = 0; comp < a_mag.nComp(); comp++)
        {
          FORT_GETMAGRESIST(CHF_FRA1(curFAB,comp),
                            CHF_CONST_REALVECT(trig),
                            CHF_CONST_REALVECT(dxv),
                            CHF_CONST_REALVECT(a_params.probLo),
                            CHF_BOX(curBox),
                            CHF_INT(comp),
                            CHF_INT(whichMag));
        }
      const EBISBox& ebisBox = a_mag[dit()].getEBISBox();
      IntVectSet ivsMulti = ebisBox.getMultiCells(curBox);

      for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect xloc = EBArith::getVofLocation(vofit(), dxv, a_params.probLo);
          for (int comp = 0; comp < a_mag.nComp(); comp++)
            {

              Real magpt;
              FORT_GETMAGPOINTRESIST(CHF_REAL(magpt),
                                       CHF_CONST_REALVECT(trig),
                                       CHF_CONST_REALVECT(xloc),
                                       CHF_CONST_INT(comp),
                                       CHF_INT(whichMag));

              a_mag[dit()](vofit(), comp) = magpt;
            }
        }
    }
}
/****************/
void
getEBAMRPFactory(RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&   a_levelOpFactory,
                 const Vector<DisjointBoxLayout>&        a_grids,
                 const Vector< EBISLayout >&             a_ebisl,
                 const PoissonParameters&                a_params,
                 const int&                              a_numPreCondIters,
                 const int&                              a_relaxType,
                 const Real&                             a_time,
                 const Real&                             a_alpha,
                 const Real&                             a_beta)
{
  CH_TIME("PoissonUtilities::getEBAMRPFactory");
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCFactory;
  RefCountedPtr<BaseEBBCFactory>         baseEBBCFactory;

  getBCFactories(baseDomainBCFactory, baseEBBCFactory, a_params);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_params.coarsestDx;
  int numLevels = -1;
  Vector<EBLevelGrid> eblg(a_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  ProblemDomain levDom =  a_params.coarsestDomain;
  ProblemDomain coarDom;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], levDom);
      if (ilev > 0)
        {
          int nref = a_params.refRatio[ilev-1];
          int nvar = 1;
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           coarDom,
                                                                           nref, nvar,
                                                                           *eblg[ilev].getCFIVS()));
          coarDom.refine(a_params.refRatio[ilev]);
        }
      coarDom = levDom;
      levDom.refine(a_params.refRatio[ilev]);
    }
  a_levelOpFactory = RefCountedPtr<EBAMRPoissonOpFactory>(new EBAMRPoissonOpFactory(eblg, a_params.refRatio, quadCFI, //
                                                                                    vectDx, RealVect::Zero,
                                                                                    a_numPreCondIters,a_relaxType,
                                                                                    baseDomainBCFactory, baseEBBCFactory,
                                                                                    a_alpha, a_beta, a_time,a_params.ghostPhi, a_params.ghostRHS,
                                                                                    numLevels));
}


/********/
void
defineSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             const Vector<EBISLayout>&                    a_ebisl,
             LinearSolver<LevelData<EBCellFAB > >&        a_bottomSolver,
             const PoissonParameters&                     a_params,
             const Real&                                  a_time,
             Real                                         a_alpha,
             Real                                         a_beta)
{
    RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > factory;
    defineSolver( a_solver,
                  a_grids,
                  a_ebisl,
                  a_bottomSolver,
                  a_params,
                  a_time,
                  a_alpha,
                  a_beta,
                  factory);

}

void
newAndFillSmoothCellData(Vector<LevelData<EBCellFAB>* >&        a_scalar,
                         const Vector<DisjointBoxLayout>&       a_grids,
                         const Vector<EBISLayout>&              a_ebisl,
                         const PoissonParameters&               a_params,
                         int a_ncomp,
                         bool a_usePhiGhost)
{
  int nlevels = a_grids.size();
  a_scalar.resize(nlevels, NULL);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  IntVect ivghost = a_params.ghostPhi;
  if (!a_usePhiGhost)
    ivghost = a_params.ghostRHS;
  ProblemDomain domLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_scalar[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], a_ncomp, ivghost, ebcellfact);

      EBLevelDataOps::setToZero(*(a_scalar[ilev]));
      for (int icomp= 0; icomp < a_ncomp; icomp++)
        {
          setTrigPhi(*a_scalar[ilev], dxLev, a_params, 0.0, icomp);
        }

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }
}
void
defineSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             const Vector<EBISLayout>&                    a_ebisl,
             LinearSolver<LevelData<EBCellFAB > >&        a_bottomSolver,
             const PoissonParameters&                     a_params,
             const Real&                                  a_time,
             Real                                         a_alpha,
             Real                                         a_beta,
             RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >& operatorFactory)
{
  CH_TIME("PoissonUtilities::defineSolver");
  ParmParse pp;
  int numPreCondIters = 40;
  pp.get("num_pre_cond_iters",numPreCondIters);
  int relaxtype;
  pp.get("mg_relax_type", relaxtype);
  if (relaxtype == 0)
    {
      pout() << "Using levelJacobi" << endl;
    }
  else if (relaxtype == 1)
    {
      pout() << "Using levelMultiColorGS" << endl;
    }
  else if (relaxtype == 2)
    {
      pout() << "Using levelGSRB" << endl;
    }
  else if (relaxtype == 3)
    {
      pout() << "Using slow multicolor relax" << endl;
    }
  else if (relaxtype == 4)
    {
      pout() << "Using levelMultiColorGS with clone" << endl;
    }
  else
    {
      MayDay::Error("bogus relax type in input file");
    }

  int numSmooths;
  pp.get("mg_num_smooths", numSmooths);
  Real hang, eps;
  pp.get("mg_hang", hang);
  pp.get("mg_eps",  eps);
  int numMG, iterMax;
  pp.get("mg_num_cycles", numMG);

  pp.get("mg_iter_max", iterMax);

  getEBAMRPFactory(operatorFactory,
                   a_grids,  a_ebisl,
                   a_params, numPreCondIters, relaxtype,
                   a_time,a_alpha,  a_beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain,
                  *operatorFactory,
                  &a_bottomSolver,
                  a_params.numLevels);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooths, numSmooths, numSmooths,
                               numMG, iterMax, eps, hang, normThresh);
  a_solver.m_verbosity = 4;
  if (pp.contains("multigrid_verbosity"))
    {
      pp.get("multigrid_verbosity", a_solver.m_verbosity);
    }
}
/*****/
/*****/
void
getVectDomainDx(Vector<ProblemDomain>& a_domain, Vector<RealVect>& a_dx, const PoissonParameters& a_params)
{
  CH_TIME("PoissonUtilities::getVectDomainDx");
  a_domain.resize(a_params.numLevels, a_params.coarsestDomain);
  a_dx    .resize(a_params.numLevels, a_params.coarsestDx);
  for (int ilev = 1; ilev < a_dx.size(); ilev++)
    {
      a_dx    [ilev] =        a_dx    [ilev-1]/a_params.refRatio[ilev-1];
      a_domain[ilev] = refine(a_domain[ilev-1],a_params.refRatio[ilev-1]);
    }
}
/********/
void
defineMGBCGSolver(BiCGStabSolver<Vector<LevelData<EBCellFAB>* > >& a_solver,
                  MultilevelLinearOp<EBCellFAB>&                   a_mlOp,
                  const Vector<DisjointBoxLayout>&                 a_grids,
                  const Vector<EBISLayout>&                        a_ebisl,
                  const PoissonParameters&                         a_params,
                  const Real&                                      a_time,
                  Real                                             a_alpha,
                  Real                                             a_beta)
{
  CH_TIME("PoissonUtilities::defineMGBCGSolver");
  Vector<ProblemDomain>  domain;
  Vector<RealVect>       dx;
  getVectDomainDx(domain,dx, a_params);

  ParmParse pp;
  int numPreCondIters = 40;
  pp.get("num_pre_cond_iters",numPreCondIters);
  int relaxtype;
  pp.get("mg_relax_type", relaxtype);
  if (relaxtype == 0)
    {
      pout() << "Using levelJacobi" << endl;
    }
  else if (relaxtype == 1)
    {
      pout() << "Using levelMultiColorGS" << endl;
    }
  else if (relaxtype == 2)
    {
      pout() << "Using levelGSRB" << endl;
    }
  else if (relaxtype == 3)
    {
      pout() << "Using slow multicolor relax" << endl;
    }
  else if (relaxtype == 4)
    {
      pout() << "Using levelMultiColorGS with clone" << endl;
    }
  else
    {
      MayDay::Error("bogus relax type in input file");
    }

  int numSmooths;
  pp.get("mg_num_smooths", numSmooths);
  Real hang, eps;
  pp.get("mg_hang", hang);
  pp.get("mg_eps",  eps);
  int numMG, iterMax;
  pp.get("mg_num_cycles", numMG);
  pp.get("mg_iter_max", iterMax);

  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > opFactory;
  getEBAMRPFactory(opFactory,
                   a_grids,  a_ebisl,
                   a_params, numPreCondIters, relaxtype,
                   a_time,a_alpha,  a_beta);


  //EBAMRPoissonOpFactory* tempPtr = (EBAMRVarCoefPoissonOpFactory *)(&(*opFactory));

  int lBase = 0;

  a_mlOp.m_num_mg_iterations = numMG;
  a_mlOp.m_num_mg_smooth = numSmooths;
  a_mlOp.define(a_grids, a_params.refRatio, domain,
              dx, opFactory, lBase);

  bool homogeneousBC = false;
  a_solver.define(&a_mlOp,homogeneousBC);
  a_solver.m_eps = eps;
  a_solver.m_verbosity = 4;
}
/*****/
int getOrderEB()
{
  CH_TIME("PoissonUtilities::getOrderEB");
  ParmParse pp;
  int orderEB;
  pp.get("order_ebbc", orderEB);
  if (orderEB == 1)
    {
      pout() << "first order EB BC" << endl;
    }
  else if (orderEB == 2)
    {
      pout() << "second order EB BC" << endl;
    }
  else
    {
      MayDay::Error("bogus EBBC order (must be 1 or 2)" );
    }
  return orderEB;
}
/*****/
RealVect getTrigRV()
{
  CH_TIME("PoissonUtilities::getTrigRV");
  RealVect trig;
  ParmParse pp;
  std::vector<Real> trigvec(SpaceDim);
  pp.getarr("trig",trigvec,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      trig[idir] = trigvec[idir];
    }
//    Real pi = 4.*atan(1.0);
//    trig *= pi;
  return trig;
}
/********/
RealVect getSphericalHarmonicRV()
{
  CH_TIME("PoissonUtilities::getSphericalHarmonicRV");
  RealVect lmphase;
  ParmParse pp;
  std::vector<Real> lmphasevec(SpaceDim);
  pp.getarr("lmphase",lmphasevec,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      lmphase[idir] = lmphasevec[idir];
    }
  return lmphase;
}
/********/
void setTrigPhi(LevelData<EBCellFAB>&    a_phi,
                const RealVect&          a_dx,
                const PoissonParameters& a_params,
                const Real&              a_time,
                int a_ivar)
{
  CH_TIME("PoissonUtilities::setTrigPhi");
  CH_assert(a_phi.nComp() > a_ivar);

  int comp = a_ivar;

  int ibox = 0;
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curPhiEBCellFAB = a_phi[dit()];
      Box curPhiBox = curPhiEBCellFAB.getRegion();
      FArrayBox& curPhiFAB = curPhiEBCellFAB.getFArrayBox();
      RealVect     trig = getTrigRV();

      FORT_GETPHI(CHF_FRA1(curPhiFAB,comp),
                  CHF_CONST_REALVECT(trig),
                  CHF_CONST_REALVECT(a_dx),
                  CHF_CONST_REAL(a_time),
                  CHF_CONST_REALVECT(a_params.probLo),
                  CHF_CONST_REALVECT(a_params.probHi),
                  CHF_BOX(curPhiBox));

      const EBISBox& curPhiEBISBox = curPhiEBCellFAB.getEBISBox();
      const EBGraph& curPhiEBGraph = curPhiEBISBox.getEBGraph();

      IntVectSet notRegular = curPhiEBISBox.getIrregIVS(curPhiBox);
      for (VoFIterator vofit(notRegular,curPhiEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5)
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETPHIPOINT(CHF_REAL(onePoint),
                           CHF_CONST_REALVECT(trig),
                           CHF_CONST_REALVECT(x),
                           CHF_CONST_REAL(a_time));

          curPhiEBCellFAB(vof,comp) = onePoint;
        }
      ibox++;
    }
}
/********/
void setSphericalHarmonicPhi(LevelData<EBCellFAB>&    a_phi,
                             const RealVect&          a_dx,
                             const PoissonParameters& a_params,
                             const Real&              a_time)
{
  CH_TIME("PoissonUtilities::setSphericalHarmonicPhi");
  CH_assert(a_phi.nComp() == 1);

  int comp = 0;

  int ibox = 0;
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curPhiEBCellFAB = a_phi[dit()];
      Box curPhiBox = curPhiEBCellFAB.getRegion();
      FArrayBox& curPhiFAB = curPhiEBCellFAB.getFArrayBox();
      RealVect     lmphase = getSphericalHarmonicRV();

      FORT_GETSHPHI(CHF_FRA1(curPhiFAB,comp),
                    CHF_CONST_REALVECT(lmphase),
                    CHF_CONST_REALVECT(a_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_REALVECT(a_params.probLo),
                    CHF_CONST_REALVECT(a_params.probHi),
                    CHF_BOX(curPhiBox));

      const EBISBox& curPhiEBISBox = curPhiEBCellFAB.getEBISBox();
      const EBGraph& curPhiEBGraph = curPhiEBISBox.getEBGraph();

      IntVectSet notRegular = curPhiEBISBox.getIrregIVS(curPhiBox);
      for (VoFIterator vofit(notRegular,curPhiEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5)
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETSHPHIPOINT(CHF_REAL(onePoint),
                             CHF_CONST_REALVECT(lmphase),
                             CHF_CONST_REALVECT(x),
                             CHF_CONST_REAL(a_time));

          curPhiEBCellFAB(vof,comp) = onePoint;
        }
      ibox++;
    }
}
/********/
void setTrigKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                        const RealVect&          a_dx,
                        const PoissonParameters& a_params,
                        const Real&              a_alpha,
                        const Real&              a_beta,
                        const Real&              a_time,
                        int a_ivar)
{
  CH_TIME("PoissonUtilities::setTrigKappaLOfPhi");
  CH_assert(a_kappaLOfPhi.nComp() > a_ivar);

  int comp = a_ivar;
  RealVect     trig = getTrigRV();

  for (DataIterator dit = a_kappaLOfPhi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curKappaLOfPhiEBCellFAB = a_kappaLOfPhi[dit()];
      Box curKappaLOfPhiBox = curKappaLOfPhiEBCellFAB.getRegion();
      FArrayBox& curKappaLOfPhiFAB = curKappaLOfPhiEBCellFAB.getFArrayBox();

      FORT_GETLOFPHI(CHF_FRA1(curKappaLOfPhiFAB,comp),
                     CHF_CONST_REALVECT(trig),
                     CHF_CONST_REALVECT(a_dx),
                     CHF_CONST_REALVECT(a_params.probLo),
                     CHF_CONST_REALVECT(a_params.probHi),
                     CHF_CONST_REAL(a_alpha),
                     CHF_CONST_REAL(a_beta),
                     CHF_CONST_REAL(a_time),
                     CHF_BOX(curKappaLOfPhiBox));

      const EBISBox& curKappaLOfPhiEBISBox = curKappaLOfPhiEBCellFAB.getEBISBox();
      const EBGraph& curKappaLOfPhiEBGraph = curKappaLOfPhiEBISBox.getEBGraph();

      IntVectSet notRegular = curKappaLOfPhiEBISBox.getIrregIVS(curKappaLOfPhiBox);
      for (VoFIterator vofit(notRegular,curKappaLOfPhiEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5 /* + centroid[idir] */)
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETLOFPHIPOINT(CHF_REAL(onePoint),
                              CHF_CONST_REALVECT(trig),
                              CHF_CONST_REALVECT(x),
                              CHF_CONST_REAL(a_alpha),
                              CHF_CONST_REAL(a_beta),
                              CHF_CONST_REAL(a_time));

          Real volfrac = curKappaLOfPhiEBISBox.volFrac(vof);
          curKappaLOfPhiEBCellFAB(vof,comp) = volfrac * onePoint;
        }
    }
}
/********/
void setTrigSource(LevelData<EBCellFAB>&    a_src,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params,
                   const Real&              a_diffConst,
                   const Real&              a_time)
{
  CH_TIME("PoissonUtilities::setTrigSource");
  CH_assert(a_src.nComp() == 1);

  int comp = 0;

  for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curSrcEBCellFAB = a_src[dit()];
      Box curSrcBox = curSrcEBCellFAB.getRegion();
      FArrayBox& curSrcFAB = curSrcEBCellFAB.getFArrayBox();
      RealVect     trig = getTrigRV();

      FORT_GETSRC(CHF_FRA1(curSrcFAB,comp),
                  CHF_CONST_REALVECT(trig),
                  CHF_CONST_REALVECT(a_dx),
                  CHF_CONST_REAL(a_time),
                  CHF_CONST_REAL(a_diffConst),
                  CHF_CONST_REALVECT(a_params.probLo),
                  CHF_CONST_REALVECT(a_params.probHi),
                  CHF_BOX(curSrcBox));

      const EBISBox& curSrcEBISBox = curSrcEBCellFAB.getEBISBox();
      const EBGraph& curSrcEBGraph = curSrcEBISBox.getEBGraph();

      IntVectSet notRegular = curSrcEBISBox.getIrregIVS(curSrcBox);
      for (VoFIterator vofit(notRegular,curSrcEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5)
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETSRCPOINT(CHF_REAL(onePoint),
                           CHF_CONST_REALVECT(trig),
                           CHF_CONST_REALVECT(x),
                           CHF_CONST_REAL(a_time),
                           CHF_CONST_REAL(a_diffConst));

          curSrcEBCellFAB(vof,comp) = onePoint;
        }
    }
}
/********/
void setSphericalHarmonicKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                                     const RealVect&          a_dx,
                                     const PoissonParameters& a_params,
                                     const Real&              a_alpha,
                                     const Real&              a_beta,
                                     const Real&              a_time)
{
  CH_TIME("PoissonUtilities::setSphericalHarmonicKappaLOfPhi");
  CH_assert(a_kappaLOfPhi.nComp() == 1);

  int comp = 0;
  RealVect     lmphase = getSphericalHarmonicRV();

  for (DataIterator dit = a_kappaLOfPhi.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curKappaLOfPhiEBCellFAB = a_kappaLOfPhi[dit()];
      Box curKappaLOfPhiBox = curKappaLOfPhiEBCellFAB.getRegion();
      FArrayBox& curKappaLOfPhiFAB = curKappaLOfPhiEBCellFAB.getFArrayBox();

      FORT_GETLOFSHPHI(CHF_FRA1(curKappaLOfPhiFAB,comp),
                       CHF_CONST_REALVECT(lmphase),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_params.probLo),
                       CHF_CONST_REALVECT(a_params.probHi),
                       CHF_CONST_REAL(a_alpha),
                       CHF_CONST_REAL(a_beta),
                       CHF_CONST_REAL(a_time),
                       CHF_BOX(curKappaLOfPhiBox));

      const EBISBox& curKappaLOfPhiEBISBox = curKappaLOfPhiEBCellFAB.getEBISBox();
      const EBGraph& curKappaLOfPhiEBGraph = curKappaLOfPhiEBISBox.getEBGraph();

      IntVectSet notRegular = curKappaLOfPhiEBISBox.getIrregIVS(curKappaLOfPhiBox);
      for (VoFIterator vofit(notRegular,curKappaLOfPhiEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;
          //the divergence of a flux lives at the centroid
          RealVect centroid = curKappaLOfPhiEBISBox.centroid(vof);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5  + centroid[idir])
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETLOFSHPHIPOINT(CHF_REAL(onePoint),
                                CHF_CONST_REALVECT(lmphase),
                                CHF_CONST_REALVECT(x),
                                CHF_CONST_REAL(a_alpha),
                                CHF_CONST_REAL(a_beta),
                                CHF_CONST_REAL(a_time));

          Real volfrac = curKappaLOfPhiEBISBox.volFrac(vof);
          curKappaLOfPhiEBCellFAB(vof,comp) = volfrac * onePoint;
        }
    }
}
/********/
void setSphericalHarmonicSource(LevelData<EBCellFAB>&    a_src,
                                const RealVect&          a_dx,
                                const PoissonParameters& a_params,
                                const Real&              a_diffConst,
                                const Real&              a_time)
{
  CH_TIME("PoissonUtilities::setSphericalHarmonicSource");
  CH_assert(a_src.nComp() == 1);

  int comp = 0;

  for (DataIterator dit = a_src.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& curSrcEBCellFAB = a_src[dit()];
      Box curSrcBox = curSrcEBCellFAB.getRegion();
      FArrayBox& curSrcFAB = curSrcEBCellFAB.getFArrayBox();
      RealVect     lmphase = getSphericalHarmonicRV();

      FORT_GETSHSRC(CHF_FRA1(curSrcFAB,comp),
                    CHF_CONST_REALVECT(lmphase),
                    CHF_CONST_REALVECT(a_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_REAL(a_diffConst),
                    CHF_CONST_REALVECT(a_params.probLo),
                    CHF_CONST_REALVECT(a_params.probHi),
                    CHF_BOX(curSrcBox));

      const EBISBox& curSrcEBISBox = curSrcEBCellFAB.getEBISBox();
      const EBGraph& curSrcEBGraph = curSrcEBISBox.getEBGraph();

      IntVectSet notRegular = curSrcEBISBox.getIrregIVS(curSrcBox);
      for (VoFIterator vofit(notRegular,curSrcEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect x;

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              x[idir] = a_dx[idir] * (iv[idir] + 0.5)
                + a_params.probLo[idir];
            }

          Real onePoint;

          FORT_GETSHSRCPOINT(CHF_REAL(onePoint),
                             CHF_CONST_REALVECT(lmphase),
                             CHF_CONST_REALVECT(x),
                             CHF_CONST_REAL(a_time),
                             CHF_CONST_REAL(a_diffConst));

          curSrcEBCellFAB(vof,comp) = onePoint;
        }
    }
}
/***************/
void getBCFactories(RefCountedPtr<BaseDomainBCFactory>& a_baseDomainBCFactory,
                    RefCountedPtr<BaseEBBCFactory>&     a_baseEBBCFactory,
                    const PoissonParameters&            a_params)
{
  CH_TIME("PoissonUtilities::getBCFactories");
  ParmParse pp;
  if (a_params.domBcType == 0)
    {
      pout() << "constant Neumann bcs on domain" << endl;

      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);

      NeumannPoissonDomainBCFactory* domainBCFactory = new NeumannPoissonDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 1)
    {
      pout() << "constant Dirichlet bcs on domain" << endl;
      DirichletPoissonDomainBCFactory* domainBCFactory = new DirichletPoissonDomainBCFactory();
      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 2)
    {
      pout() << "Neumann trig bcs on domain" << endl;
      RealVect     trig = getTrigRV();

      TrigBCFlux* trigFluxPtr = new TrigBCFlux();

      trigFluxPtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigFluxPtr = RefCountedPtr<BaseBCValue>(trigFluxPtr);

      NeumannPoissonDomainBCFactory* domainBCFactory = new NeumannPoissonDomainBCFactory();
      domainBCFactory->setFunction(baseTrigFluxPtr);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 3)
    {
      pout() << "Dirichlet trig bcs on domain" << endl;
      RealVect     trig = getTrigRV();

      TrigBCValue* trigValuePtr = new TrigBCValue();
      trigValuePtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigValuePtr = RefCountedPtr<BaseBCValue>(trigValuePtr);

      DirichletPoissonDomainBCFactory* domainBCFactory = new DirichletPoissonDomainBCFactory();
      domainBCFactory->setFunction(baseTrigValuePtr);

      a_baseDomainBCFactory = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else
    {
      MayDay::Error("Unknown domain BC type");
    }


  if (a_params.ebBcType == 0)
    {
      pout() << "constant Neumann bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      NeumannPoissonEBBCFactory* ebBC = new NeumannPoissonEBBCFactory();
      ebBC->setValue(domBCValue);

      a_baseEBBCFactory =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (a_params.ebBcType == 1)
    {
      int orderEB = getOrderEB();
      pout() << "constant Dirichlet bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      DirichletPoissonEBBCFactory* ebBC = new DirichletPoissonEBBCFactory();
      ebBC->setOrder(orderEB);
      ebBC->setValue(domBCValue);

      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else if (a_params.ebBcType == 2)
    {
      pout() << "Neumann trig bcs on EB" << endl;
      RealVect     trig = getTrigRV();
      TrigBCFlux* trigFluxPtr = new TrigBCFlux();
      trigFluxPtr->define(trig);

      RefCountedPtr<BaseBCValue> baseTrigFluxPtr =  RefCountedPtr<BaseBCValue>(trigFluxPtr);

      NeumannPoissonEBBCFactory* ebBC = new NeumannPoissonEBBCFactory();
      ebBC->setFunction(baseTrigFluxPtr);

      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else if (a_params.ebBcType == 3)
    {
      pout() << "Dirichlet trig bcs on EB" << endl;
      RealVect     trig = getTrigRV();
      RefCountedPtr<BaseBCValue> baseTrigValuePtr;
      DirichletPoissonEBBCFactory* ebBC;
      TrigBCValue* trigValuePtr = new TrigBCValue();
      trigValuePtr->define(trig);
      baseTrigValuePtr =  RefCountedPtr<BaseBCValue>(trigValuePtr);
      ebBC = new DirichletPoissonEBBCFactory();

      int orderEB = getOrderEB();
      ebBC->setOrder(orderEB);
      ebBC->setFunction(baseTrigValuePtr);

      a_baseEBBCFactory = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else
    {
      MayDay::Error("Unknown EB BC type");
    }
}
/******/
void compareError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                  const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const Vector< EBISLayout >&              a_ebislFine,
                  const Vector< EBISLayout >&              a_ebislCoar,
                  const PoissonParameters&                 a_paramsFine,
                  const PoissonParameters&                 a_paramsCoar)
{
  CH_TIME("PoissonUtilities::compareError");
  ParmParse pp;
  int testverbosity;
  pp.get("testverbosity", testverbosity);
  const Vector<int> refRat = a_paramsCoar.refRatio;
  Vector<Real> orders;
  EBArith::compareError( orders,
                         a_errorFine,
                         a_errorCoar,
                         a_gridsFine,
                         a_gridsCoar,
                         a_ebislFine,
                         a_ebislCoar,
                         refRat,
                         a_paramsCoar.coarsestDomain.domainBox(),
                         testverbosity);
}
/********/
void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                              Vector< EBISLayout >&            a_ebislCoar,
                              const Vector<DisjointBoxLayout>& a_gridsFine,
                              const PoissonParameters&         a_paramsCoar)
{
  CH_TIME("PoissonUtilities::getCoarseLayoutsFromFine");
  int nlevels = a_paramsCoar.numLevels;
  int numGhost =4;
  const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
  ProblemDomain coarsenedDomainLev = a_paramsCoar.coarsestDomain;
  a_gridsCoar.resize(nlevels);
  a_ebislCoar.resize(nlevels);
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      CH_assert(a_gridsFine[ilev].coarsenable(2));
      coarsen(a_gridsCoar[ilev], a_gridsFine[ilev], 2);

      ebisPtr->fillEBISLayout(a_ebislCoar[ilev],
                              a_gridsCoar[ilev],
                              coarsenedDomainLev,
                              numGhost);
      coarsenedDomainLev.refine(a_paramsCoar.refRatio[ilev]);
    }

}
/********/
void PoissonParameters::coarsen(int a_factor)
{
  CH_TIME("PoissonUtilities::coarsen");
  coarsestDx *= a_factor;
  coarsestDomain.coarsen(a_factor);
}
/********/
void PoissonParameters::refine(int a_factor)
{
  CH_TIME("PoissonUtilities::refine");
  coarsestDx /= a_factor;
  coarsestDomain.refine(a_factor);
}
/********/
void getPoissonParameters(PoissonParameters&  a_params, bool a_forceSingleLevel)
{
  CH_TIME("PoissonUtilities::getPoissonParameters");
  ParmParse pp;

  pp.get("eb_bc_type",a_params.ebBcType);
  pp.get("domain_bc_type",a_params.domBcType);

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
      a_params.ghostPhi[idir] = 4;
      a_params.ghostRHS[idir] = 0;
      if (a_params.ebBcType == 0 || a_params.ebBcType == 2)
        {
          a_params.ghostPhi[idir] = 1;
        }
    }
  if (pp.contains("ghostPhi") && pp.contains("ghostRhs"))
    {
      vector<int> ghost_phi_v(SpaceDim, 4);
      vector<int> ghost_rhs_v(SpaceDim, 4);
      pp.queryarr("ghostPhi",ghost_phi_v,0,SpaceDim);
      pp.queryarr("ghostRhs",ghost_rhs_v,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_params.ghostPhi[idir] = ghost_phi_v[idir];
          a_params.ghostRHS[idir] = ghost_rhs_v[idir];
        }
    }
  if (a_forceSingleLevel)
    {
      a_params.maxLevel  = 0;
      a_params.numLevels = 1;
      a_params.refRatio.resize(1,2);
      a_params.blockFactor = 8;
      a_params.fillRatio = 1;
      a_params.bufferSize = 123;
    }
  else
    {
      pp.get("max_level", a_params.maxLevel);
      a_params.numLevels = a_params.maxLevel + 1;
      pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
      pp.get("block_factor",a_params.blockFactor);
      pp.get("fill_ratio",a_params.fillRatio);
      pp.get("buffer_size",a_params.bufferSize);
    }
  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  std::vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];
  pp.queryarr("is_periodic", is_periodica,0, SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      is_periodic[dim] = (is_periodica[dim] == 1);
      if (is_periodic[dim])
        {
          pout() << "Using Periodic BCs in direction: " << dim << endl;
        }
    }

  a_params.coarsestDomain = ProblemDomain(lo, hi,is_periodic);

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

  pp.get("which_geom",   a_params.whichGeom);
  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.coarsestDx[idir] = a_params.domainLength[idir]/a_params.nCells[idir];
    }
  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;

  a_params.noRefCorners = false;
  pp.query("no_ref_corners", a_params.noRefCorners);

}

RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > newTGASolver(const Vector<DisjointBoxLayout>&       a_grids,
                                                           const Vector<EBISLayout>&              a_ebisl,
                                                           const PoissonParameters&               a_params,
                                                           BiCGStabSolver<LevelData<EBCellFAB> >& a_bottomSolver,
                                                           const int&                             a_coarsestLevel,
                                                           const Real&                            a_diffConst,
                                                           const IntVect&                         a_phiGhost,
                                                           const IntVect&                         a_rhsGhost)
{
  CH_TIME("PoissonUtilities::newTGASolver");
  ParmParse pp;
  int numPreCondIters = 40;
  pp.get("num_pre_cond_iters",numPreCondIters);
  int relaxType;
  pp.get("mg_relax_type", relaxType);
  if (relaxType == 0)
    {
      pout() << "Using levelJacobi" << endl;
    }
  else if (relaxType == 1)
    {
      pout() << "Using levelMultiColorGS" << endl;
    }
  else if (relaxType == 2)
    {
      pout() << "Using levelGSRB" << endl;
    }
  else if (relaxType == 3)
    {
      pout() << "Using slow multicolor relax" << endl;
    }
  else if (relaxType == 4)
    {
      pout() << "Using levelMultiColorGS with clone" << endl;
    }
  else
    {
      MayDay::Error("bogus relax type in input file");
    }
  int numSmooths;
  pp.get("mg_num_smooths", numSmooths);
  Real hang, eps;
  pp.get("mg_hang", hang);
  pp.get("mg_eps",  eps);
  int numMG, iterMax;
  pp.get("mg_num_cycles", numMG);
  pp.get("mg_iter_max", iterMax);
  int odeSolver;
  pp.get("ode_solver_type", odeSolver);
  int verbosity;
  pp.get("verbosity", verbosity);

  ProblemDomain coarsestProbDomain(a_params.coarsestDomain);

  RefCountedPtr<BaseDomainBCFactory> domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>         ebBCFactory;

  getBCFactories(domainBCFactory, ebBCFactory, a_params);

  int numLevels = a_grids.size();
  Vector<EBLevelGrid> eblg(a_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(a_grids.size());
  ProblemDomain levDom =  a_params.coarsestDomain;
  ProblemDomain coarDom;
  for (int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], levDom);
      if (ilev > 0)
        {
          int nref = a_params.refRatio[ilev-1];
          int nvar = 1;
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev],
                                                                           a_ebisl[ilev-1],
                                                                           coarDom,
                                                                           nref, nvar,
                                                                           *eblg[ilev].getCFIVS()));
          coarDom.refine(a_params.refRatio[ilev]);
        }
      coarDom = levDom;
      levDom.refine(a_params.refRatio[ilev]);
    }
  Real thresh = 1.0e-30;

  int numPreCond = 4;
  pp.query("mg_num_precond", numPreCond);

  Real alpha = 1.0;
  Real beta = a_diffConst;
  Real time = 0.;

  EBAMRPoissonOpFactory opFactory(eblg,
                                  a_params.refRatio,
                                  quadCFI,
                                  a_params.coarsestDx,
                                  RealVect::Zero,
                                  numPreCond,
                                  relaxType,
                                  domainBCFactory,
                                  ebBCFactory,
                                  alpha,
                                  beta,
                                  time,
                                  4*IntVect::Unit,
                                  4*IntVect::Unit);

  ProblemDomain level0Dom = eblg[0].getDomain();
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >  solver;

  a_bottomSolver.m_verbosity = 0;

  solver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);
  solver->define(level0Dom,
                 opFactory,
                 &a_bottomSolver,
                 numLevels);

  solver->m_pre        =  numSmooths;
  solver->m_post       =  numSmooths;
  solver->m_bottom     =  numSmooths;
  solver->m_hang       =  hang;
  solver->m_eps        =  eps;
  solver->m_verbosity  =  verbosity;
  solver->m_iterMax    =  iterMax;
  solver->m_normThresh =  thresh;
  solver->setMGCycle(numMG);

  Vector<LevelData<EBCellFAB>* > phi(numLevels);
  Vector<LevelData<EBCellFAB>* > rhs(numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      phi[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1, 4*IntVect::Unit, ebcellfact);
      rhs[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1, 4*IntVect::Unit, ebcellfact);
    }

  int lbase = 0;
  int lmax = numLevels-1;

  solver->init(phi,
               rhs,
               lmax,
               lbase);

  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > retval = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> >(solver,
                                       opFactory,
                                       level0Dom,
                                       a_params.refRatio,
                                       numLevels,
                                       verbosity));

  return retval;


}
/********/
void definePoissonGeometry(const PoissonParameters&  a_params)
{
  CH_TIME("PoissonUtilities::definePoissonGeometry");
  int max_level = a_params.maxLevel;

  ProblemDomain finestDomain = a_params.coarsestDomain;
  for (int ilev = 0; ilev <  max_level; ilev++)
    {
      finestDomain.refine(a_params.refRatio[ilev]);
    }

  RealVect origin = RealVect::Zero;

  RealVect fineDx = a_params.coarsestDx;
  int ebMaxCoarsen = -1;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= a_params.refRatio[ilev];
    }
  int whichgeom = a_params.whichGeom;
  int ebMaxSize = a_params.maxGridSize;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  ParmParse pp;
  if (!pp.contains("ebis_file"))
    {
      int verbosity = 0;
      if (whichgeom == 0)
        {
          //allregular
          pout() << "all regular geometry" << endl;

          AllRegularService regserv;
          ebisPtr->define(finestDomain, origin, fineDx[0], regserv, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 1)
        {
          pout() << "ramp geometry" << endl;
          RealVect rampNormal;
          vector<Real>  rampNormalVect(SpaceDim);
          pp.getarr("ramp_normal",rampNormalVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              rampNormal[idir] = rampNormalVect[idir];
            }

          Real rampAlpha;
          pp.get("ramp_alpha",rampAlpha);

          RealVect rampPoint = RealVect::Zero;
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (rampNormal[idir] != 0.0)
            {
              rampPoint[idir] = rampAlpha / rampNormal[idir];
              break;
            }
          }

          bool inside = true;

          PlaneIF ramp(rampNormal,rampPoint,inside);

          GeometryShop workshop(ramp,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 2)
        {
          pout() << "slab geometry" << endl;
          vector<int> slab_lo(SpaceDim);
          pp.getarr("slab_lo",slab_lo,0,SpaceDim);
          vector<int> slab_hi(SpaceDim);
          pp.getarr("slab_hi",slab_hi,0,SpaceDim);
          IntVect lo, hi;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              lo[idir] = slab_lo[idir];
              hi[idir] = slab_hi[idir];
            }
          Box coveredBox(lo,hi);
          SlabService slab(coveredBox);
          //this generates the new EBIS
          RealVect origin = RealVect::Zero;
          ebisPtr->define(finestDomain, origin, fineDx[0], slab, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 3)
        {
          pout() << "multi slab geometry" << endl;
          int numSlabs;
          pp.get("num_slabs",numSlabs);
          Vector<Box> coveredBoxes(numSlabs);
          for (int islab=0; islab<numSlabs; islab++)
            {
              char slabLoString[80];
              char slabHiString[80];
              sprintf(slabLoString, "slab_lo_%d", islab);
              sprintf(slabHiString, "slab_hi_%d", islab);
              vector<int> slab_lo(SpaceDim);
              pp.getarr(slabLoString,slab_lo,0,SpaceDim);
              vector<int> slab_hi(SpaceDim);
              pp.getarr(slabHiString,slab_hi,0,SpaceDim);
              IntVect lo, hi;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  lo[idir] = slab_lo[idir];
                  hi[idir] = slab_hi[idir];
                }
              Box coveredBox(lo,hi);
              coveredBoxes.push_back(coveredBox);
            }
          MultiSlabService multiSlab(coveredBoxes);
          //this generates the new EBIS
          RealVect origin = RealVect::Zero;
          ebisPtr->define(finestDomain, origin, fineDx[0], multiSlab, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 4)
        {
          RealVect cylinderAxis, cylinderCenterPt;
          vector<Real>  cylinderAxisVect(SpaceDim);
          vector<Real>  cylinderCentVect(SpaceDim,0.0);
          pp.getarr("cylinder_axis",     cylinderAxisVect, 0, SpaceDim);
          pp.queryarr("cylinder_center_pt",cylinderCentVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylinderAxis[idir]     = cylinderAxisVect[idir];
              cylinderCenterPt[idir] = cylinderCentVect[idir];
            }
          Real sum;
          PolyGeom::unifyVector(cylinderAxis, sum);

          Real cylinderRadius;
          pp.get("cylinder_radius",cylinderRadius);

          bool cylinderInside = true;
          pp.query("cylinder_inside",cylinderInside);

          pout() << "using a tilted cylinder implicit function" << endl;

          TiltedCylinderIF tunnel(cylinderRadius, cylinderAxis, cylinderCenterPt, cylinderInside);
          GeometryShop workshop(tunnel,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 5)
        {
          pout() << "sphere geometry" << endl;
          vector<Real> sphere_center(SpaceDim);
          pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
          RealVect sphereCenter;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              sphereCenter[idir] = sphere_center[idir];
            }
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);

          bool insideRegular = false;
          pp.query("inside",insideRegular);

          SphereIF implicit(sphereRadius,sphereCenter,insideRegular);

          GeometryShop workshop(implicit,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 6)
        {
          pout() << "multisphere geometry" << endl;

          int numSpheres;
          pp.get("num_spheres", numSpheres);

          Vector<Real>     radius(numSpheres);
          Vector<RealVect> center(numSpheres);

          for (int isphere = 0; isphere < numSpheres; isphere++)
            {
              char radiusString[80];
              char centerString[80];
              sprintf(radiusString, "sphere_radius_%d", isphere);
              sprintf(centerString, "sphere_center_%d", isphere);

              Real sphereRadius;
              pp.get(radiusString, sphereRadius);

              vector<Real> sphere_center(SpaceDim);
              pp.getarr(centerString,sphere_center, 0, SpaceDim);

              RealVect sphereCenter;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  sphereCenter[idir] = sphere_center[idir];
                }

              center[isphere] = sphereCenter;
              radius[isphere] = sphereRadius;
            }

          bool insideRegular = false;
          pp.query("inside",insideRegular);

          MultiSphereIF impMultiSphere(radius,center,insideRegular);

          GeometryShop workshop(impMultiSphere,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 8)
        {
          pout() << "parabolic mirror geometry" << endl;
          Real amplitude;
          RealVect center;
          vector<Real> centervec;
          int updir;
          pp.get("mirror_updir", updir);
          pp.get("mirror_amplitude", amplitude);
          pp.getarr("mirror_center",centervec, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }

          Vector<PolyTerm> poly;

          PolyTerm mono;
          Real coef;
          IntVect powers;

          if (updir != 0)
          {
            // x^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[0] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // x term
            coef = -2.0*amplitude*center[0];
            powers = IntVect::Zero;
            powers[0] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // y or z term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[0]*center[0] + center[updir];
            powers = IntVect::Zero;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }
          else
          {
            // y^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[1] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // y term
            coef = -2.0*amplitude*center[1];
            powers = IntVect::Zero;
            powers[1] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // x term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[1]*center[1] + center[updir];
            powers = IntVect::Zero;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }

          bool inside = (amplitude >= 0);

          PolynomialIF mirror(poly,inside);

          GeometryShop workshop(mirror,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 9)
        {
          pout() << "ovoid geometry" << endl;
          Real radius, axialdist;
          RealVect center;
          vector<Real> centervec;
          pp.get("ovoid_radius", radius);
          pp.get("ovoid_axial_dist", axialdist);
          int axis;
          pp.get("ovoid_axis", axis);
          pp.getarr("ovoid_center_lo",centervec, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }

          RealVect centerLo = center;

          SphereIF sphereLo(radius,centerLo,true);

          RealVect centerHi = center;
          centerHi[axis] += axialdist;

          SphereIF sphereHi(radius,centerHi,true);

          UnionIF ends(sphereLo,sphereHi);

          RealVect cylinderAxis = RealVect::Zero;
          cylinderAxis[axis] = 1.0;

          TiltedCylinderIF tube(radius,cylinderAxis,centerLo,true);
          PlaneIF loBound(cylinderAxis,centerLo,true);
          PlaneIF hiBound(cylinderAxis,centerHi,false);

          IntersectionIF loHi(loBound,hiBound);
          IntersectionIF finiteTube(tube,loHi);

          UnionIF capsule(finiteTube,ends);

          ComplementIF ovoid(capsule,true);

          GeometryShop workshop(ovoid,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 11)
        {
          pout() << "ellipsoid geometry" << endl;
          vector<Real> ellipsoid_center(SpaceDim);
          vector<Real> ellipsoid_radii(SpaceDim);
          vector<Real> ellipsoid_xaxis(SpaceDim);
          pp.getarr("ellipsoid_center",ellipsoid_center, 0, SpaceDim);
          pp.getarr("ellipsoid_radii", ellipsoid_radii, 0, SpaceDim);
          pp.getarr("ellipsoid_xaxis", ellipsoid_xaxis, 0, SpaceDim);
          int fluid_inside;
          pp.get("ellipsoid_fluid_inside",fluid_inside);
          bool insideCalc = (fluid_inside==1);
          RealVect ellipsoidCenter;
          RealVect ellipsoidRadii;
          RealVect ellipsoidXAxis;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              ellipsoidCenter[idir] = ellipsoid_center[idir];
              ellipsoidRadii[idir]  = ellipsoid_radii[idir];
              ellipsoidXAxis[idir]  = ellipsoid_xaxis[idir];
            }
          Real sum;
          PolyGeom::unifyVector(ellipsoidXAxis, sum);
          EllipsoidIF ellipsoid(ellipsoidRadii, ellipsoidCenter, insideCalc);
          RealVect origxaxis = BASISREALV(0);
          TransformIF rotazoid(ellipsoid);
          rotazoid.rotate(origxaxis, ellipsoidXAxis, ellipsoidCenter);
          GeometryShop workshop(rotazoid,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 13)
        {
          pout() << "tyler channel geometry" << endl;
          Real yDomainLength = a_params.domainLength[1];
          Real x1, x2, y1, y2;
          pp.get("tyler_x1", x1);
          pp.get("tyler_x2", x2);
          pp.get("tyler_y1", y1);
          pp.get("tyler_y2", y2);

          TylerChannelIF channel(x1, x2, y1, y2, yDomainLength);
          GeometryShop workshop(channel,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 16)
        {
          pout() << "sphere array" << endl;
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);

          vector<Real> first_sphere_center(SpaceDim);
          vector<Real> sphere_spacing(SpaceDim);
          pp.getarr("first_sphere_center",first_sphere_center, 0, SpaceDim);
          pp.getarr("sphere_spacing",sphere_spacing, 0, SpaceDim);
          RealVect firstCenter, sphereSpacing;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              firstCenter[idir] = first_sphere_center[idir];
              sphereSpacing[idir] = sphere_spacing[idir];
            }
          SphereArrayIF implicit(sphereRadius, firstCenter, sphereSpacing);
          GeometryShop workshop(implicit,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 19)
        {
          pout() << "wedge geometry" << endl;
          RealVect rampNormal1,rampNormal2;
          vector<Real>  rampNormalVect1(SpaceDim),rampNormalVect2(SpaceDim);
          pp.getarr("ramp_normal1",rampNormalVect1, 0, SpaceDim);
          pp.getarr("ramp_normal2",rampNormalVect2, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              rampNormal1[idir] = rampNormalVect1[idir];
              rampNormal2[idir] = rampNormalVect2[idir];
            }

          RealVect rampPoint1,rampPoint2;
          vector<Real>  rampPointVect1(SpaceDim),rampPointVect2(SpaceDim);
          pp.getarr("ramp_point1",rampPointVect1, 0, SpaceDim);
          pp.getarr("ramp_point2",rampPointVect2, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            rampPoint1[idir] = rampPointVect1[idir];
            rampPoint2[idir] = rampPointVect2[idir];
          }

          bool inside = true;

          PlaneIF ramp1(rampNormal1,rampPoint1,inside);
          PlaneIF ramp2(rampNormal2,rampPoint2,inside);

          IntersectionIF wedge(ramp1,ramp2);
          ComplementIF outside(wedge,true);

          GeometryShop workshop(outside,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 88)
        {
          pout() << "parabolic mirror geometry with sphere" << endl;
          Real amplitude;
          RealVect center;
          vector<Real> centervec;
          int updir;
          pp.get("mirror_updir", updir);
          pp.get("mirror_amplitude", amplitude);
          pp.getarr("mirror_center",centervec, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }

          Vector<PolyTerm> poly;

          PolyTerm mono;
          Real coef;
          IntVect powers;
          if (updir != 0)
          {
            // x^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[0] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // x term
            coef = -2.0*amplitude*center[0];
            powers = IntVect::Zero;
            powers[0] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

#if CH_SPACEDIM==3
            // z^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[2] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // z term
            coef = -2.0*amplitude*center[2];
            powers = IntVect::Zero;
            powers[2] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
#endif
            // y or z term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[0]*center[0] + center[updir];
            powers = IntVect::Zero;
#if CH_SPACEDIM==3
            coef += amplitude*center[2]*center[2];
#endif

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }
          else
          {
            // y^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[1] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // y term
            coef = -2.0*amplitude*center[1];
            powers = IntVect::Zero;
            powers[1] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

#if CH_SPACEDIM==3
            // z^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[2] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // z term
            coef = -2.0*amplitude*center[2];
            powers = IntVect::Zero;
            powers[2] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
#endif
            // x term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[1]*center[1] + center[updir];
#if CH_SPACEDIM==3
            coef += amplitude*center[2]*center[2];
#endif
            powers = IntVect::Zero;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }


          
          bool inside = (amplitude >= 0);

          PolynomialIF mirror(poly,inside);
          vector<Real> sphere_center(SpaceDim);
          pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
          RealVect sphereCenter;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              sphereCenter[idir] = sphere_center[idir];
            }
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);

          bool insideRegular = false;
          pp.query("inside",insideRegular);

          SphereIF sphere(sphereRadius,sphereCenter,insideRegular);

          IntersectionIF ends(sphere, mirror);
          GeometryShop workshop(ends,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else
        {
          //bogus which_geom
          pout() << " bogus which_geom input = "
                 << whichgeom << endl;
          MayDay::Error();
        }

    }
  else
    {
      std::string ebis_file;
      pp.get("ebis_file",ebis_file);
      pout() << " recreating  geometry from file " << ebis_file << endl;
      //define ebis anew from file input
#ifdef CH_USE_HDF5
      HDF5Handle handleIn(ebis_file, HDF5Handle::OPEN_RDONLY);
      ebisPtr->define(handleIn);
      handleIn.close();
#else
      MayDay::Error("No hdf5 available to process input ebis");
#endif
    }
  //  pout() << "Low swirl burner geometry done" << endl;
}
/***/
/***/
void getAllIrregRefinedLayouts(Vector<DisjointBoxLayout>& a_grids,
                               Vector<EBISLayout>&        a_ebisl,
                               const PoissonParameters&   a_params)
{
  CH_TIME("PoissonUtilities::getAllIrregRefinedLayouts");
  ParmParse pp;
  if (a_params.whichGeom == 0)
    {
      pout() << "all regular geometry" << endl;
      Vector<int> proc(1, 0);
      Vector<Box> coarBoxes(1, a_params.coarsestDomain.domainBox());
      a_grids.resize(a_params.numLevels);
      if (a_params.numLevels == 1)
        {
          domainSplit(a_params.coarsestDomain, coarBoxes,
                      a_params.maxGridSize, a_params.blockFactor);
          mortonOrdering(coarBoxes);
          EBEllipticLoadBalance(proc, coarBoxes, a_params.coarsestDomain);
        }
      a_grids[0] = DisjointBoxLayout(coarBoxes, proc, a_params.coarsestDomain);

      ProblemDomain coarProbDom = a_params.coarsestDomain;
      Box coarBox = coarProbDom.domainBox();
      for (int ilev = 1; ilev < a_params.numLevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if (iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(    a_params.refRatio[ilev-1]);
          coarProbDom.refine(a_params.refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          a_grids[ilev] = DisjointBoxLayout(refBoxes, proc,coarProbDom);
        }

      //fill the ebislayouts
      a_ebisl.resize(a_params.numLevels);
      int numGhost = 4;
      ProblemDomain domainLev = a_params.coarsestDomain;
      const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          //generate ebislayout
          ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],domainLev, numGhost);
          domainLev.refine(a_params.refRatio[ilev]);
        }

    }
  else
    {
      //split up coarsest domain by max box size and
      //make a dbl at coarsest level
      Vector<Box> boxesCoarsest;
      Vector<int> procsCoarsest;

      domainSplit(a_params.coarsestDomain, boxesCoarsest,
                  a_params.maxGridSize, a_params.blockFactor);

      mortonOrdering(boxesCoarsest);
      EBEllipticLoadBalance(procsCoarsest, boxesCoarsest, a_params.coarsestDomain);
      DisjointBoxLayout dblCoarsest(boxesCoarsest, procsCoarsest);

      //make a ebislayout at coarsest level.
      EBISLayout ebislCoarsest;
      const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
      ebisPtr->fillEBISLayout(ebislCoarsest,dblCoarsest,
                              a_params.coarsestDomain, 0);

      //query the input file for "which_tags" to use
      int which_tags = -99;
      int queryAns = pp.query("which_tags",which_tags);
      if (queryAns==0)
        {//file doesn't have "which_tags", tag EB
          which_tags = 0;
        }
      //make the tags
      IntVectSet tagsCoarsestLocal;
      if (which_tags == 0)
        {
          pout() << "tag all irregular cells" << endl;
          //tag all irregular coarse iv's
          IntVect loiv = a_params.coarsestDomain.domainBox().sideEnd(Side::Lo);
          IntVect hiiv = a_params.coarsestDomain.domainBox().sideEnd(Side::Hi);
          IntVect miiv = hiiv + IntVect::Unit;
          miiv /= 2;
          Box loBox(loiv, miiv);
          Box hiBox(miiv, hiiv);

          for (DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
            {
              const EBISBox& ebisBox = ebislCoarsest[dit()];
              const Box&       box = dblCoarsest.get(dit());
              tagsCoarsestLocal |= ebisBox.getIrregIVS(box);
              if (a_params.noRefCorners)
                {
                  tagsCoarsestLocal -= loBox;
                  tagsCoarsestLocal -= hiBox;
                }
            }
          // buffer the irregular cells
          int ebTagBufferSize = 0;
          pp.query("eb_tag_buffer", ebTagBufferSize);
          tagsCoarsestLocal.grow(ebTagBufferSize);
        }
      else if (which_tags == 1)
        {//make some arbitrary boxes
          pout() << "arbitrary tag set 1" << endl;
          //set up coarse tags
          int nc = a_params.coarsestDomain.size(0);
          int nmi = nc/2;//16
          int nqu = nc/4;//8
          int ntf = (nc*3)/4;  //24
          int nte = (nc*3)/8; //12
          int nfe = (nc*5)/8; //20
#if (CH_SPACEDIM ==2)
          Box boxf1(IntVect(0,  nqu), IntVect(nmi-1,ntf-1));
          Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
          Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
          Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
          Box boxf1(IntVect(0,  nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
          Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
          Box boxf3(IntVect(nqu,  0,  0), IntVect(nfe-1,nqu-1,nqu-1));
          Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
          tagsCoarsestLocal |= boxf1;
          tagsCoarsestLocal |= boxf2;
          tagsCoarsestLocal |= boxf3;
          tagsCoarsestLocal |= boxf4;
        }
      else if (which_tags == 2)
        {//make an arbitrary box
          pout() << "arbitrary tag set 2" << endl;
          //set up coarse tags
          int nc = a_params.coarsestDomain.size(0);
          int nmi = nc/2;//16
          //          int nqu = nc/4;//8
          int ntf = (nc*3)/4;  //24
          int nte = (nc*3)/8; //12
          int nfe = (nc*5)/8; //20
#if (CH_SPACEDIM ==2)
          //Box boxf1(IntVect(0,  nqu), IntVect(nmi-1,ntf-1));
          Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
          //Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
          //Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
          //Box boxf1(IntVect(0,  nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
          Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
          //Box boxf3(IntVect(nqu,  0,  0), IntVect(nfe-1,nqu-1,nqu-1));
          //Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
          //tagsCoarsestLocal |= boxf1;
          tagsCoarsestLocal |= boxf2;
          //tagsCoarsestLocal |= boxf3;
          //tagsCoarsestLocal |= boxf4;
        }
      else if (which_tags == 3)
        {
          pout() << "tagging lower half of domain" << endl;
          //set up coarse tags
          int nc = a_params.coarsestDomain.size(0);

          Box halfBox = a_params.coarsestDomain.domainBox();
          halfBox.growHi(0, -nc/2);
          tagsCoarsestLocal = IntVectSet(halfBox);
        }
      else
        {
          MayDay::Error("PoissonUtilities:: unexpected which_tags case");
        }

      //generate vector of grids.
//      print_memory_line("before BRMeshRefine");
//      UnfreedMemory();
      BRMeshRefine gridder(a_params.coarsestDomain, a_params.refRatio,
                           a_params.fillRatio,      a_params.blockFactor,
                           a_params.bufferSize,     a_params.maxGridSize);
//      print_memory_line("after BRMeshRefine");
//      UnfreedMemory();

      Vector<Vector<Box> > newMeshes(a_params.numLevels);
      Vector<Vector<Box> > oldMeshes(a_params.numLevels);
      //      oldMeshes[0]= Vector<Box>(1, a_params.coarsestDomain.domainBox());
      oldMeshes[0]= boxesCoarsest;
      for (int ilev = 1; ilev < a_params.numLevels; ilev++)
        {
          oldMeshes[ilev] = Vector<Box>(1, refine(oldMeshes[ilev-1][0], a_params.refRatio[ilev-1]));
        }
      int baseLevel = 0;
//      print_memory_line("before regrid");
//      UnfreedMemory();
      gridder.regrid(newMeshes, tagsCoarsestLocal, baseLevel, a_params.maxLevel, oldMeshes);
//      print_memory_line("after regrid");
//      UnfreedMemory();

      /*
        If this is  a single level calculation, then we can make the layout more efficient by
        omitting completely covered boxes in the layout.  This is perilous in an AMR calculation,
        however because then one might violate proper nesting constraints.  There be dragons.
      */

      bool useMortonOrdering =  true;
      bool useEllipticLoadBalance = true;
      pp.query("use_morton_ordering", useMortonOrdering);
      pp.query("use_timed_load_balance", useEllipticLoadBalance);
      Vector<Vector<int> > newProcs(a_params.numLevels);
      ProblemDomain domLevel = a_params.coarsestDomain;
      if(useMortonOrdering)
        {
          pout() << "using morton ordering for boxes" << endl;
        }
      else
        {
          pout() << "NOT using morton ordering for boxes" << endl;
        }
      if(useEllipticLoadBalance)
        {
          pout() << "using timed load balance" << endl;
        }
      else
        {
          pout() << "using standard load balance" << endl;
        }

      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          if(useMortonOrdering)
            {
              mortonOrdering(newMeshes[ilev]);
            }
          if(useEllipticLoadBalance)
            {
              EBEllipticLoadBalance(newProcs[ilev], newMeshes[ilev], domLevel);
            }
          else
            {
              LoadBalance(newProcs[ilev], newMeshes[ilev]);
            }
          domLevel.refine(a_params.refRatio[ilev]);
        }

      a_grids.resize(a_params.numLevels);
      a_ebisl.resize(a_params.numLevels);
      int numGhost = 4;
      ProblemDomain domainLev = a_params.coarsestDomain;
//      print_memory_line("before fillEBISLayout");
//      UnfreedMemory();
      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          a_grids[ilev] = DisjointBoxLayout(newMeshes[ilev], newProcs[ilev],domainLev);
          //generate ebislayout
          ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],
                                  domainLev, numGhost);
          domainLev.refine(a_params.refRatio[ilev]);
        }
 //     print_memory_line("after fillEBISLayout");
 //     UnfreedMemory();
    }
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = a_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += a_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += a_grids[ilev].size();
      pout() << "getAllIrregRefineLayouts:level[" << ilev
             << "], number of boxes = " << a_grids[ilev].size()
             << ", number of points = " << pointsThisLevel << endl;
    }
  pout() << "getAllIrregRefineLayouts:"
         <<  "   total boxes = " << totalBoxes
         <<  ", total points = " << totalPoints <<  endl;
}

/***/

void getViscousBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_domBC,
                           RefCountedPtr<BaseEBBCFactory>&                         a_ebbc,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_eta,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_lambda,
                           const Vector<DisjointBoxLayout>&                        a_grids,
                           const Vector<EBISLayout>&                               a_ebisl,
                           const PoissonParameters&                                a_params)
{
  CH_TIME("PoissonUtilities::getViscousBCFactories");
  ParmParse pp;
  if (a_params.domBcType == 0)
    {
      pout() << "constant Neumann bcs on domain" << endl;

      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);

      NeumannViscousTensorDomainBCFactory* domainBCFactory = new NeumannViscousTensorDomainBCFactory();
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else if (a_params.domBcType == 1)
    {
      pout() << "constant Dirichlet bcs on domain" << endl;
      DirichletViscousTensorDomainBCFactory* domainBCFactory = new DirichletViscousTensorDomainBCFactory();
      Real domBCValue;
      pp.get("domain_bc_value", domBCValue);
      domainBCFactory->setValue(domBCValue);

      a_domBC = RefCountedPtr<BaseDomainBCFactory>(domainBCFactory);
    }
  else
    {
      MayDay::Error("Unknown domain BC type");
    }


  if (a_params.ebBcType == 0)
    {
      pout() << "constant Neumann bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      NeumannViscousTensorEBBCFactory* ebBC = new NeumannViscousTensorEBBCFactory();
      ebBC->setValue(domBCValue);

      a_ebbc =RefCountedPtr<BaseEBBCFactory>( ebBC);
    }
  else if (a_params.ebBcType == 1)
    {
      pout() << "constant Dirichlet bcs on EB" << endl;
      Real domBCValue;
      pp.get("eb_bc_value", domBCValue);
      DirichletViscousTensorEBBCFactory* ebBC = new DirichletViscousTensorEBBCFactory();
      ebBC->setValue(domBCValue);

      a_ebbc = RefCountedPtr<BaseEBBCFactory>(ebBC);
    }
  else
    {
      MayDay::Error("Unknown EB BC type");
    }
}

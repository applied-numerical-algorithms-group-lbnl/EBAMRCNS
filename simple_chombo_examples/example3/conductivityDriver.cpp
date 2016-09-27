#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBConductivityOp.H"
#include "EBConductivityOpFactory.H"
#include "EBLevelDataOps.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "EBSimpleSolver.H"
#include "BiCGStabSolver.H"
#include "EBEllipticLoadBalance.H"
#include "EBLoadBalance.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "SphereIF.H"
#include "EBAMRDataOps.H"
#include "GeometryShop.H"
#include "BaseIVFactory.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include "memusage.H"
#include "memtrack.H"

class PoissonParameters
{
public:
  IntVect       nCells;
  int           maxGridSize;
  int           blockFactor;
  int           bufferSize;
  Real          fillRatio;
  int           maxLevel;
  int           numLevels;
  Vector<int>   refRatio;
  ProblemDomain coarsestDomain;
  RealVect      coarsestDx;
  RealVect      domainLength;
  RealVect      probLo;
  RealVect      probHi;
  void coarsen(int a_factor);
  void  refine(int a_factor);

};

/********/
void getPoissonParameters(PoissonParameters&  a_params)
{
  CH_TIME("PoissonUtilities::getPoissonParameters");
  ParmParse pp;


  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
  pp.get("block_factor",a_params.blockFactor);
  pp.get("fill_ratio",a_params.fillRatio);
  pp.get("buffer_size",a_params.bufferSize);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  a_params.coarsestDomain = ProblemDomain(Box(lo, hi));

  a_params.domainLength = RealVect::Unit;

  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.coarsestDx[idir] = a_params.domainLength[idir]/a_params.nCells[idir];
    }
  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;
}
/***/
void getAllIrregRefinedLayouts(Vector<DisjointBoxLayout>& a_grids,
                               Vector<EBISLayout>&        a_ebisl,
                               const PoissonParameters &  a_params)
{
  CH_TIME("PoissonUtilities::getAllIrregRefinedLayouts");

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

  //make the tags
  IntVectSet tagsCoarsestLocal;
  pout() << "tag all irregular cells" << endl;

  for (DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = ebislCoarsest[dit()];
      const Box&       box = dblCoarsest.get(dit());
      tagsCoarsestLocal |= ebisBox.getIrregIVS(box);
    }

  //generate vector of grids.
  BRMeshRefine gridder(a_params.coarsestDomain, a_params.refRatio,
                       a_params.fillRatio,      a_params.blockFactor,
                       a_params.bufferSize,     a_params.maxGridSize);

  Vector<Vector<Box> > newMeshes(a_params.numLevels);
  Vector<Vector<Box> > oldMeshes(a_params.numLevels);

  oldMeshes[0]= boxesCoarsest;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      oldMeshes[ilev] = Vector<Box>(1, refine(oldMeshes[ilev-1][0], a_params.refRatio[ilev-1]));
    }
  int baseLevel = 0;
  gridder.regrid(newMeshes, tagsCoarsestLocal, baseLevel, a_params.maxLevel, oldMeshes);

  Vector<Vector<int> > newProcs(a_params.numLevels);
  ProblemDomain domLevel = a_params.coarsestDomain;
  pout() << "using morton ordering for boxes" << endl;
  pout() << "using timed load balance" << endl;

  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
          
      mortonOrdering(newMeshes[ilev]);
      EBEllipticLoadBalance(newProcs[ilev], newMeshes[ilev], domLevel);
      domLevel.refine(a_params.refRatio[ilev]);
    }

  a_grids.resize(a_params.numLevels);
  a_ebisl.resize(a_params.numLevels);
  int numGhost = 4;
  ProblemDomain domainLev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      a_grids[ilev] = DisjointBoxLayout(newMeshes[ilev], newProcs[ilev],domainLev);
      //generate ebislayout
      ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],
                              domainLev, numGhost);
      domainLev.refine(a_params.refRatio[ilev]);
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

  int ebMaxSize = a_params.maxGridSize;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();


  pout() << "sphere geometry" << endl;
  RealVect sphereCenter = 0.5*RealVect::Unit;
  Real sphereRadius = 0.1;
  bool insideRegular = false;
  SphereIF implicit(sphereRadius,sphereCenter,insideRegular);

  GeometryShop workshop(implicit,0,fineDx);
  //this generates the new EBIS
  ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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

  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      EBFluxFactory        ebfluxfact(a_ebisl[ilev]);
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev]);

      a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(a_grids[ilev], 1, 4*IntVect::Unit, ebcellfact));
      a_bco[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(a_grids[ilev], 1, 4*IntVect::Unit, ebfluxfact));
      a_bcoIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(a_grids[ilev], 1, 4*IntVect::Unit, baseivfact));

      ParmParse pp;
      Real acoVal, bcoVal;
      pp.get("acoef_value", acoVal);
      pp.get("bcoef_value", bcoVal);

      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_aco     [ilev])[dit()].setVal(acoVal);
          (*a_bco     [ilev])[dit()].setVal(bcoVal);
          (*a_bcoIrreg[ilev])[dit()].setVal(bcoVal);
        }
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
  pout() << "homobeneous Dirichlet bcs on domain and eb" << endl;
  DirichletConductivityDomainBCFactory* domBCRaw = new DirichletConductivityDomainBCFactory();
  domBCRaw->setValue(0.);
  DirichletConductivityEBBCFactory* ebBCRaw = new DirichletConductivityEBBCFactory();
  //2 if you want johansen, 1 if you want least squares
  ebBCRaw->setOrder(1);
  ebBCRaw->setValue(0.);

  RefCountedPtr<BaseDomainBCFactory>      domBC(domBCRaw);
  RefCountedPtr<BaseEBBCFactory>          ebBC(ebBCRaw);

  int relaxType = 1;
  pout() << "using multicolor gauss-seidel relaxation" << endl;
  IntVect ghost = 4*IntVect::Unit;
  a_factory = RefCountedPtr<EBConductivityOpFactory>
    (new EBConductivityOpFactory(eblg, quadCFI, alpha, beta, a_aco, a_bco, a_bcoIrreg,
                                 a_params.coarsestDx[0],  a_params.refRatio, domBC, ebBC,
                                 ghost, ghost, relaxType));
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
/********/
void solve(const PoissonParameters&  a_params)
{
  int nvar = 1;
  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;
  getAllIrregRefinedLayouts(grids, ebisl, a_params);


  //define  data
  Vector<LevelData<EBCellFAB>* > phi(a_params.numLevels);
  Vector<LevelData<EBCellFAB>* > rhs(a_params.numLevels);
  ParmParse pp;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      EBCellFactory factory(ebisl[ilev]);
      phi[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, 4*IntVect::Unit, factory);
      rhs[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, 4*IntVect::Unit, factory);

      //for now just set phi to zero and rhs to -1.
      EBLevelDataOps::setVal(*phi[ilev], 0.0);
      EBLevelDataOps::setVal(*rhs[ilev], 1.0);
    }

  //create the solver
  AMRMultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  solver" << endl;

  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  bicgstab.m_verbosity = 0;
  EBSimpleSolver simp;
  simp.setNumSmooths(100);

  LinearSolver<LevelData<EBCellFAB> >* botSolve = &bicgstab;

  defineConductivitySolver(solver, grids, ebisl, *botSolve, a_params);

  pout() << "solving " << endl;
  //solve the equation
  solver.init(phi, rhs, a_params.maxLevel, 0);

  int numSolves;

  pp.get("number_of_solves", numSolves);
  for(int isolve = 0; isolve < numSolves; isolve++)
    {
      pout() << "starting solve number " << isolve << endl;
      EBAMRDataOps::setVal(phi, 0.0);
      solver.solveNoInit(phi, rhs, a_params.maxLevel, 0);
    }


#ifdef CH_USE_HDF5
  bool fileOut;
  pp.get("do_file_output", fileOut);
  if (fileOut)
    {
      pout() << "outputting the answer to file" << endl;

      //output the answer
      char charstr[100];
      sprintf(charstr, "phi.%dd.hdf5", SpaceDim);
      writeEBAMRname(&phi, charstr);

      sprintf(charstr, "rhs.%dd.hdf5", SpaceDim);
      writeEBAMRname(&rhs, charstr);
    }
#endif
  //clean up memory
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      delete phi[ilev] ;
      delete rhs[ilev] ;
    }
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters params;

    //read params from file
    getPoissonParameters(params);

    //define geometry from given params
    definePoissonGeometry(params);



    //solve the stinking problem and output everything
    solve(params);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

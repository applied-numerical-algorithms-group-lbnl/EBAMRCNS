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


#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBAMRIO.H"
#include "VoFIterator.H"
#include "SphereIF.H"
#include "EBArith.H"
#include "AggStencil.H"
#include "RegFortranF_F.H"
#include "GeometryShop.H"


void
getFluxStencil(VoFStencil      &   a_fluxStencil,
               const FaceIndex &   a_face,
               const EBISBox   &   a_ebisBox,
               const Real      &   a_dx)
{

  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights

  //everything is more complicated near a coarse-fine interface
  IntVectSet cfivs;

  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     cfivs,
                                                     a_ebisBox,
                                                     a_ebisBox.getDomain());
  a_fluxStencil.clear();
  for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      if (!a_face.isBoundary())
        {
          //face centered flux is just the centered difference gradient.
          faceCentSten.add(face.getVoF(Side::Hi),  1.0/a_dx, 0);
          faceCentSten.add(face.getVoF(Side::Lo), -1.0/a_dx, 0);
        }
      else
        {
          //the boundary condition handles this one.
        }

      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
}
/***/
void
getKappaDivFStencil(VoFStencil       &   a_vofStencil,
                    const VolIndex   &   a_vof,
                    const EBISBox    &   a_ebisBox,
                    const Real       &   a_dx)
{
  a_vofStencil.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, faces[iface],
                             a_ebisBox, a_dx);

              Real areaFrac = a_ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/a_dx;
              a_vofStencil += fluxStencil;
            }
        }
    }
}
///
void
getStencilAtPoint(RefCountedPtr<BaseStencil> & a_stencil, 
                  RefCountedPtr<BaseIndex>   & a_destVoF,
                  const VolIndex             & a_vof, 
                  const EBISBox              & a_ebis, 
                  const Real                 & a_dx)
{
  VoFStencil vofsten;
  getKappaDivFStencil(vofsten, a_vof, a_ebis, a_dx);
  //someone more clever than I about casting could probably avoid these.
  a_destVoF = RefCountedPtr<BaseIndex  >(new VolIndex(a_vof));
  a_stencil = RefCountedPtr<BaseStencil>(new VoFStencil(vofsten));
}
/****/
void  
getAggStencils(LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB> >  >&  a_stencils, 
               const LevelData<EBCellFAB>                                    &  a_klp, 
               const LevelData<EBCellFAB>                                    &  a_phi, 
               const DisjointBoxLayout                                       &  a_grids, 
               const EBISLayout                                              &  a_ebisl, 
               const Real                                                    &  a_dx)
{
  //allocate one aggstencil per box
  a_stencils.define(a_grids);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box    & grid = a_grids[dit()];
      const EBISBox& ebis = a_ebisl[dit()];
      //only need special stencils on cut cells
      IntVectSet ivsIrreg = ebis.getIrregIVS(grid);
      VoFIterator vofit(ivsIrreg, ebis.getEBGraph());
      Vector<VolIndex> vofs = vofit.getVector();
      Vector<RefCountedPtr< BaseIndex>   > destVoFs(vofs.size());
      Vector<RefCountedPtr< BaseStencil> > stencils(vofs.size());
      for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          //get the stencil for this vof and put it into the vector
          getStencilAtPoint(stencils[ivof], destVoFs[ivof],
                            vofs[ivof], ebis, a_dx);
        }
      //send all those slow stencils to make one fast one
      a_stencils[dit()] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB> >
        (new AggStencil<EBCellFAB, EBCellFAB>(destVoFs, stencils, a_phi[dit()], a_klp[dit()]));
    }
}

//evaluate kappa*Laplacian(phi)
void
getKappaLaplacian(LevelData<EBCellFAB>         &   a_klp, 
                  LevelData<EBCellFAB>         &   a_phi, 
                  const DisjointBoxLayout      &   a_grids, 
                  const EBISLayout             &   a_ebisl, 
                  const Real                   &   a_dx)
{
  //this would be member data of some class
  LayoutData<RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > > stencils;
  getAggStencils(stencils, a_klp, a_phi, a_grids, a_ebisl, a_dx);
  //fill ghost cells (no flux through domain or EB here)
  //you can archive the ExchangeCopier here to improve performance 
  a_phi.exchange();
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      
      //because i did this incrementally
      a_klp[dit()].setVal(0.);
      const Box    & grid = a_grids[dit()];
      const EBISBox& ebis = a_ebisl[dit()];
      //evaluate everywhere as if there is there is no EB

      BaseFab<Real>& regKlp = a_klp[dit()].getSingleValuedFAB();
      BaseFab<Real>& regPhi = a_phi[dit()].getSingleValuedFAB();
      
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          Box loBox, hiBox, centerBox;
          int hasLo, hasHi;
          EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox,
                              ebis.getDomain(), grid, idir);
          FORT_INCR1DLAPLACIAN(CHF_FRA1(regKlp, 0),
                               CHF_FRA1(regPhi, 0),
                               CHF_REAL(a_dx),
                               CHF_INT(idir),
                               CHF_BOX(loBox),
                               CHF_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_INT(hasHi),
                               CHF_BOX(centerBox));
        }

      //fix up at irregular cells using aggregated stencil
      //incrOnly=false is telling the stencil to set the answer to zero first
      bool incrOnly = false;
      int varDest = 0;
      stencils[dit()]->apply(a_klp[dit()], a_phi[dit()], varDest, incrOnly);
    }
}

//fill phi = (some smooth function)
void
fillPhiData(LevelData<EBCellFAB>    &   a_phi, 
            const DisjointBoxLayout &   a_grids, 
            const EBISLayout        &   a_ebisl, 
            const Real              &   a_dx)
{
  static const Real g_pi = 4.*atan(1.);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box    & grid = a_grids[dit()];
      const EBISBox& ebis = a_ebisl[dit()];
      IntVectSet ivs(grid);
      for(VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect vofloc = EBArith::getVoFLocation(vofit(), a_dx, RealVect::Zero);
          Real value = 1;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              value *= cos(g_pi*vofloc[idir]);
              a_phi[dit()](vofit(), 0) = value;
            }
        }
    }
}
//make data and evaluate kappa*Laplacian(data)
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Scoping trick
  {

    //make grids and geometry
    DisjointBoxLayout grids;
    int nx = 64;
    int maxGrid = 32;
    int blockFactor = 16;
    Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);

    Vector<Box> boxes;
    domainSplit(domain,  boxes,
                maxGrid, blockFactor);

    Vector<int> procs;
    LoadBalance(procs, boxes);
    grids.define(boxes, procs);
    
    int nghost = 4;

    //this is the implicit function
    Real radius = 0.1;
    RealVect center = 0.5*RealVect::Unit;
    
    //the false means the fluid is outside the sphere
    SphereIF sphere(radius, center, false);
    RealVect origin = RealVect::Zero;
    Real dx = 1.0/nx;
    
    GeometryShop shop(sphere, 0, dx*RealVect::Unit);
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    //this defines and the geometric moments at all refinements
    ebisPtr->define(domain, origin, dx, shop, maxGrid, -1);

    EBISLayout  ebisl;

    ebisPtr->fillEBISLayout(ebisl,
                            grids,
                            domain,
                            nghost);

    //define the data
    EBCellFactory fact(ebisl);
    LevelData<EBCellFAB> phi(grids, 1, nghost*IntVect::Unit, fact);
    LevelData<EBCellFAB> lph(grids, 1, nghost*IntVect::Unit, fact);


    //fill phi = (some smooth function)
    pout() << "filling phi with some smooth data" << endl;
    fillPhiData(phi, grids, ebisl, dx);

    //evaluate kappa*Laplacian(phi)
    pout() << "evaluating kappa*Laplacian(phi) with homogeneous Neumann BCs everywhere" << endl;
    pout() << "the answer will look odd near the sphere" << endl;
    getKappaLaplacian(lph, phi, grids, ebisl, dx);


#ifdef CH_USE_HDF5
    pout() << "writing out the answers"<< endl;
    writeEBLevelname(&(phi),"pltPhi.hdf5");
    writeEBLevelname(&(lph),"pltKappaLapl.hdf5");
#endif
  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

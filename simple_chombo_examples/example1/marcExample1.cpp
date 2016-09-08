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

#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBAMRIO.H"
#include "VoFIterator.H"
#include "SphereIF.H"
#include "GeometryShop.H"

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Scoping trick
  {

    DisjointBoxLayout gridsFine,  gridsCoar;
    int nxFine = 64;
    int maxGridFine = 32;
    int blockFactorFine = 16;
    Box domainFine(IntVect::Zero, (nxFine-1)*IntVect::Unit);
    Box domainCoar = domainFine;
    domainCoar.coarsen(2);

    Vector<Box> boxes;
    domainSplit(domainFine,  boxes,
                maxGridFine, blockFactorFine);

    Vector<int> procs;
    LoadBalance(procs, boxes);
    gridsFine.define(boxes, procs);
    coarsen(gridsCoar, gridsFine, 2);
    
    int nghost = 4;

    //this is the implicit function
    Real radius = 0.45;
    RealVect center = 0.5*RealVect::Unit;
    
    SphereIF sphere(radius, center, true);
    RealVect origin = RealVect::Zero;
    Real fineDx = 1.0/nxFine;
    
    GeometryShop shop(sphere, 0, fineDx*RealVect::Unit);
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    //this defines the geometric moments
    ebisPtr->define(domainFine, origin, fineDx, shop, maxGridFine, -1);

    EBISLayout  ebislFine, ebislCoar;
    ebisPtr->fillEBISLayout(ebislCoar,
                            gridsCoar,
                            domainCoar,
                            nghost);

    ebisPtr->fillEBISLayout(ebislFine,
                            gridsFine,
                            domainFine,
                            nghost);

    EBCellFactory factFine(ebislFine);
    EBCellFactory factCoar(ebislCoar);
    LevelData<EBCellFAB> dataFine(gridsFine, 1, nghost*IntVect::Unit, factFine);
    LevelData<EBCellFAB> dataCoar(gridsCoar, 1, nghost*IntVect::Unit, factCoar);

    //they are over the same proc layout so you can use the same iterator
    for(DataIterator dit = gridsFine.dataIterator(); dit.ok(); ++dit)
      {
        Box        bCoar = gridsCoar[dit()];
        EBISBox ebisCoar = ebislCoar[dit()];
        IntVectSet ivsCoar(bCoar);
        for(VoFIterator vofit(ivsCoar, ebisCoar.getEBGraph()); vofit.ok(); ++vofit)
          {
            IntVect iv = vofit().gridIndex();
            dataCoar[dit()](vofit(), 0) =  (Real)(iv[0]+ iv[1]);
          }
        Box bFine        = gridsFine[dit()];
        EBISBox ebisFine = ebislFine[dit()];
        IntVectSet ivsFine(bFine);
        for(VoFIterator vofit(ivsFine, ebisFine.getEBGraph()); vofit.ok(); ++vofit)
          {
            IntVect iv = vofit().gridIndex();
            dataFine[dit()](vofit(), 0) = (Real)(iv[0]+ iv[1]);
          }
      }

#ifdef CH_USE_HDF5
    writeEBLevelname(&(dataFine),"pltFine.hdf5");
    writeEBLevelname(&(dataCoar),"pltCoar.hdf5");
#endif
  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

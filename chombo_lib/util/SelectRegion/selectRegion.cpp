#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code reads a plotfile and writes another.
// The output file contains a subset of the input file defined by a index
// range on the coarsest level.

#include <iostream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"

// One more function for MPI
void dumpmemoryatexit();

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif

  // ------------------------------------------------
  // parse command line and input file
  // ------------------------------------------------
  // infile must be first
  if (argc < 2)
    {
      cerr << "  need inputs file" << endl;
      abort();
    }

  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, in_file);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  // declare variable to store hierarchy
  string                         inFileName;
  Vector<DisjointBoxLayout>      inGrids;
  Vector<LevelData<FArrayBox>* > inData;
  Vector<string>                 inVars;
  Box                            inDomain;
  Real                           inDx;
  Real                           inDt;
  Real                           inTime;
  Vector<int>                    inRefRatio;
  int                            inNumLevels;

  pp.get("infile" ,inFileName);

  string outFileName;
  pp.get("outfile",outFileName);

  Vector<int> loEndVect;
  pp.getarr("lo_end",loEndVect,0,SpaceDim);

  Vector<int> hiEndVect;
  pp.getarr("hi_end",hiEndVect,0,SpaceDim);

  IntVect loEnd;
  IntVect hiEnd;

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    loEnd[idir] = loEndVect[idir];
    hiEnd[idir] = hiEndVect[idir];
  }

  Box boxOfInterest(loEnd,hiEnd);

  ReadAMRHierarchyHDF5(inFileName,
                       inGrids,
                       inData,
                       inVars,
                       inDomain,
                       inDx,
                       inDt,
                       inTime,
                       inRefRatio,
                       inNumLevels);

  Vector<LevelData<FArrayBox>* > outData(inNumLevels);
  Vector<DisjointBoxLayout>      outGrids(inNumLevels);

  // allocate output data -- same domain as input
  for (int level = 0; level < inNumLevels; level++)
    {
      const DisjointBoxLayout curDBL = inGrids[level];

      Vector<Box> outBoxes;
      Vector<int> outProcs;

      for (LayoutIterator lit = curDBL.layoutIterator(); lit.ok(); ++lit)
      {
        Box curBox = curDBL.get(lit());

        if (curBox.intersects(boxOfInterest))
        {
          outBoxes.push_back(curBox & boxOfInterest);
        }
      }

      LoadBalance(outProcs,outBoxes);

      outGrids[level].define(outBoxes,outProcs);

      outData[level] = new LevelData<FArrayBox>(outGrids[level],
                                                inVars.size(),
                                                inData[level]->ghostVect());
      // copy data for this level
      Interval inInterval(0,inVars.size());
      inData[level]->copyTo(*outData[level]);

      boxOfInterest.refine(inRefRatio[level]);
    }

  WriteAMRHierarchyHDF5(outFileName,
                        inGrids,
                        outData,
                        inVars,
                        inDomain,
                        inDx,
                        inDt,
                        inTime,
                        inRefRatio,
                        inNumLevels);

  // clean up memory
  for (int level = 0; level < inNumLevels; level++)
    {
      delete inData[level];
      delete outData[level];
    } 

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main

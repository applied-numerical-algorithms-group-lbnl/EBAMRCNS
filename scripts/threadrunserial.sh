#!/bin/csh -f
if ($#argv != 5) then
  echo "Usage: threadrun minthread maxthread timedir exec inputs "
  exit 1
endif

set minProcs = $1
set maxProcs = $2

set timedir = $3
set exec    = $4
set inputs  = $5

set procs = $minProcs
echo "executable name = $exec"
echo "input file name = $inputs"
echo "time directory  = $timedir"

if (! -e $timedir) then
    mkdir $timedir
endif

while ($procs <= $maxProcs)
 echo "number of threads = $procs"
 set command = "setenv OMP_NUM_THREADS $procs"
 echo $command
 $command
 set outfile = "screen.out"
if (-e $outfile) then
    rm $outfile
endif
 set command = "$exec $inputs > $outfile"
 #set command = "cp $inputs $outfile"
 echo $command 
 $command
 set pprocs = `echo $procs | awk '{printf("%06d",$1);}'`    
set outdir  = $timedir/amrg.time.$pprocs.thread
if (! -e $outdir) then
   mkdir $outdir
endif 

 set command = "cp time.table $outfile  $outdir"
 echo $command
 $command

#  echo $dir

  set procs = `expr 1 +  $procs`
end

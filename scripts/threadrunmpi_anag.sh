#!/bin/csh -f
if ($#argv != 7) then
  echo "Usage: threadrun minthread maxthread timedir exec inputs minproc maxproc"
  exit 1
endif

set minThread = $1
set maxThread = $2

set timedir = $3
set exec    = $4
set inputs  = $5
set minProc  = $6
set maxProc  = $7

set procs   = $minProc
echo "executable name = $exec"
echo "input file name = $inputs"
echo "time directory  = $timedir"
echo "minProcs = $minProc"
echo "maxProcs = $maxProc"

if (! -e $timedir) then
    set command = "mkdir $timedir"
    $command
endif

while ($procs   <= $maxProc)
 echo "number of procs = $procs"
 set threads = $minThread

while ($threads <= $maxThread)
 echo "number of threads = $threads"
 set command =     "setenv OMP_NUM_THREADS $threads"
 echo $command
 $command

 set command = "mpirun -np $procs -machinefile $machfile $exec $inputs"
 echo $command
 $command

 set procfile = `echo $procs | awk '{printf("%06d",$1);}'`    
 
 set command =  "mv time.table.0 $timedir/amrg.time.$threads.threads.$procs.procs"
 echo $command
 $command

 set threads = `expr 1 +  $threads`
end
 set procs = `expr 1 +  $procs`

end

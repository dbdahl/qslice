export OMP_NUM_THREADS=1

rnd=2
dte=240626

num_of_trials=$(Rscript functions/num_of_trials_all.R $rnd)
echo $num_of_trials

for ii in $( seq 1 $num_of_trials )
do
  echo "Rscript --no-restore 2_run_trials_all.R $rnd $ii $dte &> logs/all_rnd$rnd\_job$ii\_dte$dte.txt"
done | parallel --jobs 10
wait

echo "finished"

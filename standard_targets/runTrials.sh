export OMP_NUM_THREADS=1

target="normal"
rnd=1
dte=240308

num_of_trials=$(Rscript functions/num_of_trials.R $target $rnd)
echo $num_of_trials

for ii in $( seq 1 $num_of_trials )
do
  echo "Rscript --no-restore 2_run_trial.R $target $rnd $ii $dte &> logs/$target\_rnd$rnd\_job$ii\_dte$dte.txt"
done | parallel --jobs 10
wait

echo "finished"

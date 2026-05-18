export OMP_NUM_THREADS=1

dte=260518
dat="db40"

num_of_trials=$(Rscript functions/num_of_trials.R $dat)
echo $num_of_trials

for ii in $( seq 1 $num_of_trials )
do
  echo "Rscript --no-restore 2_run_trials_all.R $dat $ii $dte &> logs/$dat\_job$ii\_dte$dte.txt"
done | parallel --jobs 50
wait

echo "finished"

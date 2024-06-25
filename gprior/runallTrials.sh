export OMP_NUM_THREADS=1

dte=240520

num_of_trials=$(Rscript functions/num_of_trials.R)
echo $num_of_trials

for ii in $( seq 1 $num_of_trials )
do
  echo "Rscript --no-restore 2_run_trials_all.R $ii $dte &> logs/all_job$ii\_dte$dte.txt"
done | parallel --jobs 10
wait

echo "finished"

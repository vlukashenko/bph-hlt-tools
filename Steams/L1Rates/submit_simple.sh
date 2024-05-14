Executable   = step3_rate_estimation.sh
arguments    = 5000000
#Requirements = MY.WantOS == "el7"
Output       = condor_outputs/$(Cluster)_$(Process).stdout
Error        = condor_outputs/$(Cluster)_$(Process).stderr
Log          = condor_outputs/$(Cluster)_$(Process).log
+MaxRuntime  = 60*60*6
notification = Always
#notify_user = your.email
queue
#queue maxE from runs.txt

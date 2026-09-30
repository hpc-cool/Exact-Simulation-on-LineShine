cd EXP00

time mpirun -np 8 --allow-run-as-root -map-by ppr:1:NUMA:pe=32 -x LD_LIBRARY_PATH --bind-to core -x UCX_TSL=self,sm -x UCX_RC_VERBS_TX_MIN_SGE=2 -x UCX_UD_VERBS_TX_MIN_SGE=1 ./nemo
tail -n 1 run.stat

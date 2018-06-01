# NB this command is specific to the configuration on lxplus and is not gaurenteed elsewhere
#outdir="/afs/cern.ch/work/s/sethzenz/ws/" # can't set absolute path on lsf because we're expecting to stage
queue="1nd"
useAAA=0
version="1030c"
# doHTXS is false but rapidity filter will be applied
fggRunJobs.py --load bbhth_sig_jobs_2017.json -d bbhth_sig_jobs_$version -x cmsRun workspaceStd.py maxEvents=-1 -n 500 -q $queue -D -P useAAA=$useAAA doHTXS=False doFiducial=False tthTagsOnly=True puTarget=2.588e+05,1.084e+06,2.012e+06,3.779e+06,4.056e+06,5.878e+06,6.452e+06,6.833e+06,9.253e+06,2.182e+07,4.37e+07,8.265e+07,1.315e+08,1.886e+08,2.665e+08,3.753e+08,5.235e+08,6.978e+08,8.712e+08,1.031e+09,1.172e+09,1.285e+09,1.372e+09,1.44e+09,1.498e+09,1.553e+09,1.6e+09,1.629e+09,1.635e+09,1.615e+09,1.571e+09,1.509e+09,1.431e+09,1.339e+09,1.239e+09,1.137e+09,1.038e+09,9.472e+08,8.657e+08,7.967e+08,7.438e+08,7.105e+08,6.982e+08,7.061e+08,7.299e+08,7.616e+08,7.907e+08,8.053e+08,7.953e+08,7.554e+08,6.86e+08,5.944e+08,4.911e+08,3.873e+08,2.924e+08,2.122e+08,1.486e+08,1.01e+08,6.693e+07,4.343e+07,2.772e+07,1.748e+07,1.092e+07,6.786e+06,4.208e+06,2.614e+06,1.631e+06,1.025e+06,6.507e+05,4.176e+05,2.709e+05,1.776e+05,1.174e+05,7.814e+04,5.224e+04,3.5e+04,2.346e+04,1.571e+04,1.05e+04,6987,4631,3053,2002,1304,843.6,541.8,345.3,218.3,136.8,85.01,52.34,31.94,19.3,11.55,6.85,4.021,2.337,1.344,0.7655,0.4315

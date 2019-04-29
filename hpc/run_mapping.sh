#!/bin/bash
export PYTHONPATH=/home/syang29/.conda/envs/my_env/lib/python2.7/site-packages:$PYTHONPATH
cd /work/syang29/Monsif_pbs/
python mapping_serial.py -inputCSV $1 -out output/$2.out 2>&1

#!/bin/bash

export PATH=/usr/local/cuda-10.0/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-10.0/lib64:$LD_LIBRARY_PATH

conda deactivate

#python3 hybrid-RNN-evaluation_ec2.py
#python3 grid_search_hybrid_RNN_ec2.py
python3 hybrid-RNN-evaluation_TSS_ec2.py

#sudo init 0

import os
import subprocess
import cellbender
#import scanpy as sc
#import matplotlib.pyplot as plt
#import numpy as np
#import numba

os.chdir("/work/ABG/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X")
sample_dirs = [d for d in os.listdir() if d.startswith('Sample_') and os.path.isdir(d)]
for sample_dir in sample_dirs:
    input_file = os.path.join(sample_dir, f"{sample_dir}_output", 'outs', 'raw_feature_bc_matrix.h5')
    if os.path.isfile(input_file):
        #os.chdir()
        #os.makedirs(f"{sample_dir}_cb_10k",exist_ok =True)
        output_file = os.path.join(sample_dir,f"{sample_dir}_cb_10k", f"{sample_dir}_cb_10k.h5")
        # Construct and run the cellbender command
        command = f"cellbender remove-background --input {input_file} --output {output_file} --total-droplets-included 10000 --cuda"
        subprocess.run(command, shell=True)
    else:
        print(f"raw_feature_bc_matrix.h5 not found in {sample_dir}")


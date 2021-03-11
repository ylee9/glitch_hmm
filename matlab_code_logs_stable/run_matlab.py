#!/usr/bin/env python3
import glob
import subprocess
import os
from mpi4py import MPI
import sys

def chunk(list, n):
    return [list[i::n] for i in range(n)]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    #data_dir = os.getenv('DATA_DIR')
    data_dir = sys.argv[1]
    print(data_dir)
    z_files = glob.glob(data_dir + '/data/point_*_z.dat')
    uuids = [x.split("_")[-2].strip() for x in z_files]
    if os.path.isdir(data_dir + '/results_data'):
        uuids_to_remove = []
        for uuid in uuids:
            if os.path.isfile(data_dir + '/results_data/point_' + uuid + '_path.dat'):
                uuids_to_remove.append(uuid)

        for uuid in uuids_to_remove:
            uuids.remove(uuid)

    uuids = chunk(uuids, comm.Get_size())
else:
    uuids = None

uuids = comm.scatter(uuids, root=0)

print("Got " + str(len(uuids)) + " UUIDs")
print(uuids)
for uuid in uuids:
    #matlab_code_dir = os.getenv('MATLAB_CODE_DIR')
    matlab_code_dir = sys.argv[2]
    print(os.getenv('PATH'))
    cmd = "export UUID=" + uuid + ";matlab -singleCompThread -nodesktop -r \"run('" + matlab_code_dir + "/process_single_uuid.m');exit\""
    print(subprocess.check_output(cmd, shell=True, executable='/bin/zsh'))



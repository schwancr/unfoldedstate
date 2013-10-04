#!/usr/bin/env python

"""
A simple script that computes the crysol SAXS profile for each
snapshot in an MSMBuilder Trajectories/ dir.

-- TJL Oct 2013
"""


import sys
import glob
import os
import os.path.join as pjoin
import argparse

import tables
import numpy as np
import multiprocessing as mp

from mdtraj import io
from msmbuilder import Trajectory


crysol_params = {
                'lm' : 15,    # Maximum order of harmonics (default 15) 
                'fb' : 17,    # Order of Fibonacci grid (default 17)
                'sm' : 0.5,   # Maximum scattering vector (default 0.5)
                'ns' : 51,    # num points in computed curve (default 51)  
                'un' : 1,     # Angular units (default 1)                
                'dns': 0.334, # Solvent density (default 0.334 e/A**3)   
                'dro': 0.03,  # Contrast of hydration shell (0.03 e/A**3)
                }


DEVNULL = open(os.devnull,'w')

def setup_tmp_dir():

    if not os.path.exists('/scratch/tjlane/crysol'):
        os.system('mkdir -p /scratch/tjlane/crysol')

    global TEMPDIR
    TEMPDIR = '/scratch/tjlane/crysol'

    return


def compute_crysol(trajectory, save_to=None):
    """
    Compute crysol for all the snapshots in an msmbuilder trajectory.
   
    Parameters
    ----------
    trajectory : msmbulder.Trajectory.trajectory
        The trajectory to compute SAXS for

    save_to : str or None
        If this is a string, will save to an h5 file of that name.    

    Returns
    -------
    q_values : np.ndarray
        The q_values at which the scattering was computed, in ()

    scattering_pred : np.ndarray
        The estimated integrated intensity for each `q_value`
    """

    if type(trajectory) == str:
        trajectory = Trajectory.load_from_file(trajectory)

    for i in range(len(trajectory)):

        frame = trajectory[i]

        pdbfn = '%s/tmp4crysol.pdb' % TEMPDIR
        frame.save_to_pdb(pdbfn)

        # run crysol comand line
        args = ['/%s %s' % kv for kv in crysol_params.items()]
        cmd = ['crysol', pdbfn] + args + ['2>&1 /dev/null']
        subprocess.check_call(cmd, shell=True)

        # parse the output
        intensities_output = 'tmp4crysol00.int'
        if not os.path.exists(intensities_output):
             raise IOError('crysol output not found')

        d = genfromtxt(intensities_output, skip_header=1)
        q_values = d[:,0]

        # initialize output space
        if not scattering_pred:
            scattering_pred = np.zeros((len(trajectory), d.shape[0]))

        scattering_pred[i,:] = d[:,3]

        os.remove(pdbfn)
        os.remove(intensities_output)
        os.remove('tmp4crysol00.alm')
        os.remove('tmp4crysol00.log')


    if save_to:
        io.saveh(save_to, q_values=q_values, saxs=scattering_pred)
        print "Saved: %s" % save_to

    return q_values, scattering_pred


def main(in_dir, extension, out_dir):

    pool = mp.Pool(mp.cpu_count())

    print "Searching for files matching:", in_dir + '*' + extension
    trajectory_list = glob(in_dir + '*' + extension)
    print "Found: %d trajectory files" % len(trajectory_list)

    output_names = [ pjoin(out_dir, os.path.basename(t_fn).split('.')[:-1] + '-crysol.h5') for t_fn in trajectory_list ]

    result = pool.map_async(compute_crysol, zip(trajectory_list, output_names))
    result.wait()

    DEVNULL.close()

    return


if __name__ == '__main__':

    argparse.ArgumentParser(description='Map crysol onto trajectory data')
    parser.add_argument('in_dir', help='Directory containing trajectories to compute SAXS for')
    parser.add_argument('-e', '--extension', help='Trajectory extension to look for. Default: lh5', default='lh5') 
    parser.add_argument('-o', '--out_dir', help='Directory to save output in. Default: crysol.', default='crysol')
    args = parser.parse_args()

    main(args.in_dir, args.extension, args.out_dir)


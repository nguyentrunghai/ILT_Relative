
import argparse
import os

from _process_yank_outputs import read_reduced_energy
from _mbar import mbar_weights_in_holo

parser = argparse.ArgumentParser()
parser.add_argument('--yank_dir', type=str )
parser.add_argument('--repeat', type=int, default=3)
parser.add_argument('--nc_file', type=str, default='complex-implicit.nc')
parser.add_argument('--ligand_name', type=str, default=None)

args = parser.parse_args()

yank_nc_files = []
for repeat in range(args.repeat):
    yank_nc_files.append( os.path.join(args.yank_dir, "Repeat%d" %repeat, "output", args.nc_file)  )

print yank_nc_files

u_kln , N_k = read_reduced_energy(yank_nc_files)
weights = mbar_weights_in_holo(u_kln, N_k)

out_file = open(args.ligand_name + '.weig', 'w')
out_file.write( "# mbar weights of %s\n" %args.ligand_name )
for i in range( len(weights) ):
    out_file.write("%10d %20.10e\n" %( i, weights[i] ) )
out_file.close()



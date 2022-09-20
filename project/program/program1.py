import sys
import numpy as np
from utils import *
import subprocess
from subprocess import Popen, PIPE, STDOUT
from tqdm import tqdm

if(len(sys.argv)!=2):
    assert print('usage: python program_1.py options.dat')

#parsing options.inp
with open(sys.argv[1],'r') as file:
    buffer=file.read().split('\n')[:-1]

sel_dihs_inp_path=buffer[0].split(';')[0].split()[0]
n_atoms=int(buffer[1].split(';')[0].split()[0])
n_molecules=int(buffer[2].split(';')[0].split()[0])
trj_paths=buffer[3].split(';')[0].split()[0]
input_itp=buffer[4].split(';')[0].split()[0]
gmx_dir=buffer[5].split(';')[0].split()[0]

#reading what dihedral angels should be read
with open(sel_dihs_inp_path,'r') as file:
    sel_dihs=np.array([ [int(el2) for el2 in el.split()]for el in file.read().split('\n')[:-1]])


manage_dir('result')

#parsing itp file (contains topology of molecule)
dihedrals, atom_types, atom_type_names=itp_parser(input_itp,gmx_dir)

#parsing paths of trj
with open(trj_paths,'r') as file:
    trj_paths_names={line.split(":")[0].split()[0] : line.split(":")[1].split()[0] for line in file.read().split('\n')[:-1]}


#for each of sel_dih gets infroamtion what types/how many of dihedrals are depended on this degree of freedom
#degree of freedom = axis desribed by 2nd and 3rd atom of dihedral
tors_info=get_types_from_tors(sel_dihs, atom_type_names,dihedrals)


#calculating populations
for i_sel, sel_dih in tqdm(enumerate(sel_dihs), total=len(sel_dihs)):
    xy={}
    for trj_path, trj_name in zip(trj_paths_names.keys(), trj_paths_names.values()):
        manage_dir('run')
        #create ndx
        with open('run/ndx.ndx','w') as file:
            file.write('[ dih1 ]\n')
            for i_mol in range(n_molecules):
                file.write('%d %d %d %d\n' % (sel_dih[0]+i_mol*n_atoms,
                                            sel_dih[1]+i_mol*n_atoms,
                                            sel_dih[2]+i_mol*n_atoms,
                                            sel_dih[3]+i_mol*n_atoms,))
        #run calc
        cmd=('gmx angle -type dihedral -n run/ndx.ndx -f %s -od run/result.xvg' % trj_path)
        subprocess.run(cmd.split(), capture_output=True)
        angs=read_angs('run/result.xvg')
        xy[trj_name]=angs

    #generate plot
    generate_plot(xy, tors_info[i_sel], "-".join([str(el)for el in sel_dih]))

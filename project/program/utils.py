import os
from subprocess import Popen, PIPE, STDOUT
import numpy as np
from matplotlib import pyplot as plt

def read_angs(path):
    with open(path, 'r') as file:
        angs={ int(line.split()[0]):float(line.split()[1]) for line in file.read().split('\n')[:-1] if line[0] != '#' and line[0] != '@'}
    return angs

def generate_plot(xys, tor_info, title):
    text_kwargs = dict(ha='center', va='center', fontsize=20, color='black')
    text_kwargs2 = dict(ha='center', va='center', fontsize=10, color='blue')
    text_kwargs3 = dict(ha='center', va='center', fontsize=10, color='red')
    fig, ax1 = plt.subplots(figsize=(10,7))

    M=np.array([np.array([el[key]for key in el.keys()]).min() for el in xys.values()])
    min_point=np.min(M)
    M=np.array([np.array([el[key]for key in el.keys()]).max() for el in xys.values()])
    max_point=np.max(M)

    for label in xys.keys():
        xy_dict=xys[label]
        ax1.plot(list(xy_dict.keys()), list(xy_dict.values()),label=label)

    # print(tor_info)
    for i,dih in enumerate(tor_info):
        dih_raw=dih.split(':')[0].split('-')
        if(dih_raw[0][0]=="H" or dih_raw[-1][0]=="H"):
            ax1.text(0,max_point-((max_point-min_point)/100)-i*0.002,dih, text_kwargs2)
        else:
            ax1.text(0,max_point-((max_point-min_point)/100)-i*0.002,dih, text_kwargs3)

    ax1.set(xlabel='dih_ang [*]', ylabel='Population',
           title=title)

    ax1.grid()
    ax1.legend()
    fig.savefig("result/"+title+".png")

def get_types_from_tors(sel_dihs, atom_type_names,dihedrals):
    result=[]
    for sel_dih in sel_dihs:
        types_in_sel_dih=[]
        for dihedral in dihedrals:
            if([sel_dih[1],sel_dih[2]]==dihedral[1:3] or dihedral[1:3]==[sel_dih[2],sel_dih[1]]):
                dih_is=np.array(dihedral)-1
                dih_type=[atom_type_names[i] for i in dih_is]
                types_in_sel_dih.append(dih_type)
        current_dict={}
        for dih in types_in_sel_dih:
            dihrev=[dih[3],dih[2],dih[1],dih[0]]
            keys_list=[el.split('-') for el in list(current_dict.keys())]
            if dih not in keys_list and dihrev not in keys_list:
                current_dict["-".join(dih)]=1
            else:
                if(dih in keys_list):
                    current_dict["-".join(dih)]+=1
                else:
                    current_dict["-".join(dihrev)]+=1
        result.append([key+":"+str(current_dict[key]) for key in current_dict.keys()])
    return result

def manage_dir(dir):
    if(os.path.isdir(dir)):
        os.system('rm -r '+dir)
        os.mkdir(dir)
    else:
        os.mkdir(dir)

def itp_parser(input_itp, gmx_dir):
    #parsing itp file for dihedrals
    with open(input_itp, 'r') as file:
        buffer=file.read().split("\n")

    dihedrals=[]
    for i1 in range(len(buffer)):
        line=buffer[i1]
        if("[ dihedrals ]" in line):
            for i2 in range(len(buffer)-i1):
                line = buffer[i1+i2+1].split(";")[0]
                if(len(line.split())==5):
                    dihedrals.append([int(_) for _ in line.split()[:-1]])
                else:
                    break
            break

    atom_types=[]
    for i1 in range(len(buffer)):
        line=buffer[i1]
        if("[ atoms ]" in line):
            for i2 in range(len(buffer)-i1):
                line = buffer[i1+i2+2].split(";")[0]
                if(len(line.split())==8):
                    atom_types.append(line.split()[1])
                else:
                    break
            break

    atom_type_names = []
    for i in range(len(atom_types)):
        phrase=("%s" % atom_types[i])
        cmd = ['grep',phrase,gmx_dir+'/share/gromacs/top/oplsaa.ff/ffnonbonded.itp']
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        atom_type_names.append(stdout.decode('utf-8').split()[1])

    return dihedrals, atom_types, atom_type_names

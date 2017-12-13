#! /usr/bin/python
#
# Simple parser to extract contacts
# Chain ids or TER not required
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from ForceField import VdwParamset
from ResLib import  ResiduesDataLib
from StructureWrapper import *
from Bio.PDB.ResidueDepth import *
from Bio.PDB.NACCESS import *

import sys
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
COVLNK = 2.0
HBLNK  = 3.5

all_polars = [
    'N', 'ND1', 'ND2', 'NE',  'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
    'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',  'OG1', 'OH',
    'S', 'SD',  'SG'
]
backbone_polars =  ['N','O']
waternames = ['WAT','HOH']

def main():

    parser = argparse.ArgumentParser(
                prog='polarContacts',
                description='Polar contacts detector'
            )

    parser.add_argument(
        '--backonly',
        action='store_true',
        dest='backonly',
        help='Restrict to backbone'
    )
    parser.add_argument( #Argument to plot the different number of polar contacts that each residue has
        '--plot',
        action='store_true',
        dest='plotMode',
        help='Restrict to sidechain'
    )
    parser.add_argument(
        '--nowats',
        action='store_true',
        dest='nowats',
        help='Exclude water molecules'
    )
    
    parser.add_argument(
        '--diel',
        type= float,
        action='store',
        dest='diel',
        default = 1.0,
        help='Relative dielectric constant'
    )
    
    parser.add_argument(
        '--vdw',
        action='store',
        dest='vdwprm',
        help='VDW Paramters file'
    )
    
    parser.add_argument(
        '--rlib',
        action='store',
        dest='reslib',
        help='AminoAcid library'
    )

    parser.add_argument( 
        '--index',
        action ='store',
        dest= 'index',
        help= 'Only select the molecules with more stable contact by an index'
    )

    parser.add_argument(
        '--surf',
        action='store_true',
        dest='surf',
        help='Usa ASA',
    )

    parser.add_argument('pdb_path')

    args = parser.parse_args()

    print ("Settings")
    print ("--------")
    for k,v in vars(args).items():
        print ('{:10}:'.format(k),v)

    backonly = args.backonly
    nowats =args.nowats
    plotMode= args.plotMode
    pdb_path = args.pdb_path
    vdwprm = args.vdwprm
    reslib = args.reslib
    index= args.index
    surf= args.surf
    diel = args.diel
    
# Load VDW parameters
    vdwParams = VdwParamset(vdwprm)
    print ("{} atom types loaded".format(vdwParams.ntypes))

# Load AA Library
    aaLib = ResiduesDataLib(reslib)
    print ("{} amino acid atoms loaded".format(aaLib.nres))
    
    if not pdb_path:
        parser.print_help()
        sys.exit(2)

    parser = PDBParser(PERMISSIVE=1)

    try:
        st = parser.get_structure('st', pdb_path)
    except OSError:
        print ("#ERROR: loading PDB")
        sys.exit(2)

# Checking for model
    if len(st) > 1:
        print ("#WARNING: Several Models found, using only first")

# Using Model 0 any way
    st = st[0]
    # Getting surfaces
    if surf:
        res_surfaces = NACCESS(st,naccess_binary='/Users/daniel/Downloads/BioPhysics-energies0/NACCESS/naccess')
        at_surfaces = NACCESS_atomic(st,naccess_binary='/Users/daniel/Downloads/BioPhysics-energies0/NACCESS/naccess')
        print ("Surfaces obtained from NACCESS")

# Making a list of polar atoms
    polats = []
    if backonly:
        selected_atoms = backbone_polars

    else:
        selected_atoms = all_polars

    for at in st.get_atoms():
        if at.id in selected_atoms:
            polats.append(at)
#Searching for contacts under HNLNK on diferent residues
    nbsearch = NeighborSearch(polats)
    hblist = []
    for at1, at2 in nbsearch.search_all(HBLNK):
        if at1.get_parent() == at2.get_parent():
            continue
 #Discard covalents and neighbours
        if (at1-at2) < COVLNK:
            continue
        if abs(at2.get_parent().id[1] - at1.get_parent().id[1]) == 1:
            continue
# remove waters
        if nowats:
            if at1.get_parent().get_resname() in waternames \
                or at2.get_parent().get_resname() in waternames:
                continue
                
        atom1 = Atom(at1,1,aaLib,vdwParams)
        atom2 = Atom(at2,1,aaLib,vdwParams)        
        if at1.get_serial_number() < at2.get_serial_number():
            hblist.append([at1, at2])
        else:
            hblist.append([at2, at1])
       
    print ()
    print ("Polar contacts")
    print ('{:13} {:13} {:6} '.format(
            'Atom1','Atom2','Dist (A)')
    )
    if plotMode: #Dictionary to save as a key the residues with polar contacts and as values the number of contacts on sidechain or mainchain
        sidecontacts={}
        maincontacts= {}
    res= set() #A set to store the residues involved on polar contacts. stored as a set to avoid repeated residues
    
    for hb in sorted (hblist,key=lambda i: i[0].get_serial_number()):
        r1 = hb[0].get_parent()
        r2 = hb[1].get_parent()
        res.add(str(r1.id[1])+ r1.get_resname())
        print ('{:14} {:14} {:6.3f} '.format(
            r1.get_resname()+' '+str(r1.id[1])+hb[0].id,
            r2.get_resname()+' '+str(r2.id[1])+hb[1].id,
            hb[0] - hb[1]
            )
        )
        if surf:
            print ('{:14} ({:6.3f}) {:14} ({:6.3f}) {:6.3f} '.format(
                hb[0].atid(), float(hb[0].at.xtra['EXP_NACCESS']),
                hb[1].atid(), float(hb[1].at.xtra['EXP_NACCESS']),
                hb[0].at - hb[1].at
                )
            )
        if plotMode:
            if hb[0].id in backbone_polars and hb[1].id in backbone_polars:
                if r1.id[1] not in maincontacts:
                    maincontacts[r1.id[1]] = 1
                    if r1.id[1] not in sidecontacts:
                        sidecontacts[r1.id[1]] = 0
                else:
                    maincontacts[r1.id[1]] += 1
           
            else:
                if r1.id[1] not in sidecontacts:
                    
                    if r1.id[1] not in maincontacts:
                        maincontacts[r1.id[1]] = 0

                    sidecontacts[r1.id[1]] = 1
                
                else:
                    sidecontacts[r1.id[1]] += 1
    if plotMode: #Plot the number of polar contacts of each residue on a bar graph.
                    #Contacts on mainchain are blue bars, contacts on sidechain red bars

        N= len(res)
        nindex= np.arange(N)
        main_con= maincontacts.values()
        side_con= sidecontacts.values()
        fig,ax = plt.subplots()
        bar_width= 0.35

        rec1=ax.bar(nindex, main_con, bar_width, color='b',label='Mainchain contacts')
        rec2=ax.bar(nindex + bar_width, side_con, bar_width, color='r',label='Sidechain contacts')
        ax.set_xlabel('AminoAcid Number')
        ax.set_ylabel('Number of polar contacts')
        ax.set_xticks(nindex+bar_width /2, res)
        ax.set_xticklabels(res)
        ax.legend()
        fig.tight_layout()
        plt.show()





        
    print ()
    print ("Residue interactions")

# Making list or residue pairs to avoid repeated pairs
    respairs = []
    for hb in hblist:
        r1 = Residue(hb[0].get_parent(), 1, aaLib, vdwParams)
        r2 = Residue(hb[1].get_parent(), 1, aaLib, vdwParams)
        if [r1,r2] not in respairs:
            respairs.append([r1,r2])
    l=[]
    for rpair in sorted(respairs, key=lambda i: i[0].resNum()):            
        eint = rpair[0].elecInt(rpair[1],diel)
        evdw = rpair[0].vdwInt(rpair[1])
        print (
            '{:10} {:10} {: 8.4f} {: 8.4f} {: 8.4f}'.format(
                rpair[0].resid(), 
                rpair[1].resid(),
                eint,
                evdw,
                eint+evdw)
        )
        l.append([rpair[0].resid(),rpair[0],rpair[1].resid(),rpair[1],eint,evdw, eint+evdw])
    
    if index is not None: #Select only the residues with less energy (the most stables ones)
        for e, element in enumerate(sorted(l, key=lambda i: i[6])):
            if e < int(index):
                print(element)
    if surf:
        srfr1 = float(rpair[0].residue.xtra['EXP_NACCESS']['all_polar_rel'])
        srfr2 = float(rpair[1].residue.xtra['EXP_NACCESS']['all_polar_rel'])
        # Define 30% threshold for buried
        if (srfr1>30.) or (srfr2>30.):
            diel0=80.0
        else:
            diel0=diel
        eint = rpair[0].elecInt(rpair[1],diel0)
        print (
            '{:10} ({:>8.3f}) {:10}  ({:>8.3f}) (e: {:4.1f}) {:>8.4f} {:>8.4f} {:>8.4f}'.format(
            rpair[0].resid(), srfr1,
            rpair[1].resid(), srfr2,
            diel0, eint,
            evdw,
            eint+evdw)
        )

        
if __name__ == "__main__":
    main()

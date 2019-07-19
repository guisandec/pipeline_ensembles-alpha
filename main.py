# Archivo principal del pipeline_ensembles
# Autor: GUISANDE DONADIO, CE

import os, glob, shutil, sys, csv
from pprint import pprint as pp
from Bio.PDB import PDBParser, PDBIO, PPBuilder
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_protein# FUNCIONES Y DATOS GLOBALES

path = ""
verbose = True
data = {}

def cargar_txt_en_lista(file_path):
    return_list = []
    with open(file_path,"r") as openfile:
        for lines in openfile:
            return_list.append(lines.replace("\n",""))
            
    return return_list

def get_uniprot_from_pdb(pdb_id,chain):
    with open("archivos_importantes/pdb_chain_uniprot.csv","r") as openfile:
        temp = csv.reader(openfile,delimiter=",")
        for item in temp:
            if (item[0].upper()== pdb_id) and  (item[1]==chain):
                return item[2]
                break
            else:
                continue
    return None

def get_seq_from_uniprot(uniprot):
    seqfile = "fasta/seq_"+uniprot+".fasta"
    cmd = "wget https://www.uniprot.org/uniprot/"+uniprot+".fasta -O "+seqfile
    os.system(cmd)
    if os.path.isfile:
        for entry in SeqIO.parse(seqfile,"fasta"):
            seq = str(entry.seq)
    return seq

def bajar_estructura(pdb_id):
    os.chdir(path+"/ent_files/")
    filename ="pdb"+pdb_id.lower()+".ent"
    url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"+pdb_id.lower()[1:3]+"/pdb"+pdb_id.lower()+".ent.gz"
    if not os.path.isfile(filename):
        #os.system("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"+pdb_id.lower()+".ent.gz")
        wget_code = os.system("wget "+url)
        if wget_code == 2048:
            os.chdir(path)
            return ("OBSOLETE")
        #Uncompress
        os.system("gunzip pdb"+pdb_id.lower()+".ent.gz")
    else:
        if verbose: print (">> File "+filename+" already exists")
    os.chdir(path)
    return filename


def split_pdb_by_chain(pdb_id):
    if not os.path.isdir("pdb_chains/"+pdb_id.upper()):
            os.mkdir("pdb_chains/"+pdb_id.upper())
    actual_pdbfile = PDBParser().get_structure(pdb_id,"ent_files/pdb"+pdb_id.lower()+".ent")
    return_dict = dict()
    for model in actual_pdbfile:
        for chain in model:
            outfilename = pdb_id.upper() + "-" + str(model.get_id()+1) +  "_" + str(chain.get_id()) + ".pdb"
            if not os.path.isfile("pdb_chains/"+pdb_id.upper()+"/"+outfilename):
                io = PDBIO()
                io.set_structure(chain)
                io.save("pdb_chains/"+pdb_id.upper()+"/"+outfilename)
            ppb = PPBuilder().build_peptides(chain)
            this_seq = Seq("",generic_protein)
            for pp in ppb:
                 this_seq += pp.get_sequence()
            return_dict[outfilename]=this_seq
    return return_dict


def read_expdata(path_to_pdbfile):
    with open(path_to_pdbfile,"r") as openfile:
        expdata = ""
        nummodel = 1
        for row in openfile:
            if "EXPDTA" in row:
                expdata +=  (row[9:-1].strip())
                if "X-RAY" in expdata:
                    break
                else:
                    continue
            elif "NUMMDL" in row:
                nummodel = int((row[9:-1].strip()))
                break
            elif "REMARK" in row:
                break
            elif "ATOM" in row:
                break
    return expdata,nummodel



def read_atom_full(pdb_file):
    return_lst = []
    with open(pdb_file,"r") as openfile:
        for row in openfile:
            if row[0:5] == "ATOM ":
                tmp_lst = []
                tmp_lst.append(int(row[6:10+1])) #Integer serial Atom serial number.
                tmp_lst.append(row[12:15+1].strip()) #Atom name Atom name.
                tmp_lst.append(row[16+1]) #Character altLoc Alternate location indicator.
                tmp_lst.append(row[17:19+1]) #Residue name resName Residue name.
                tmp_lst.append(row[21+1]) #Character chainID Chain identifier.
                tmp_lst.append(int(row[22:25+1])) #Integer resSeq Residue sequence number.
                #tmp_lst.append(row[26+1]) #AChar iCode Code for insertion of residues.
                #tmp_lst.append(row[30:37+1]) #Real(8.3) x Orthogonal coordinates for X in Angstroms.
                #tmp_lst.append(row[38:45+1]) #Real(8.3) y Orthogonal coordinates for Y in Angstroms.
                #tmp_lst.append(row[46:53+1]) #Real(8.3) z Orthogonal coordinates for Z in Angstroms.
                #tmp_lst.append(row[54:59+1]) #Real(6.2) occupancy Occupancy.
                #tmp_lst.append(row[60:65+1]) #Real(6.2) tempFactor Temperature factor.
                #tmp_lst.append(row[76:77+1]) #LString(2) element Element symbol, right-justified.
                #tmp_lst.append(row[78:89+1].replace("\n","")) #LString(2) charge Charge on the atom.
                return_lst.append(tmp_lst)
    return return_lst

def print_fixed(query_str,fix=100):
    for index,char in enumerate(query_str):
        print (char,end="")
        if (index+1)%fix ==0:
            print ("")
    print ("")
    return

def save_alidic_to_fasta(ali_dic,outpath):
    with open (outpath,"w") as openfile:
        for itemsin in ali_dic:
            openfile.write (">"+item+"\n")
            openfile.write (ali_dic[item]+"\n")


d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X'}

import sys
# import subprocess
import os
import getopt
from numpy import mean
# import re

#declaring the input variables
args=sys.argv[1:]                                                                   # get arguments from usr's input; filename is sys.arg[0], so args start from [1:]
last_200_pdb = ''
distance_xvg = ''
delta_distance = 0.2
ligand_name = ''
output_name = ''
max_delta_distance = 1


#optarg for the input example-c 2 -i PLDXPAL -t done.trg -p 50 -N 6 -n 80000
try:
   opts, args = getopt.getopt(args,"h:f:d:t:n:o:",["help","trajectory_file =",          # getopot.getopt(sys.arg, short_option‘-h,-i,-t,-p,etc’, long_option'--help,--input_seq,--receptor)
                                    "distance_xvg =",                                   # with usr's input e.g -i PLDXPAL -c2, the 'getopt' function can grab them and save them seprately into 'opts' and 'args'
                                    "delta_distance =",
                                    "ligand_name =",
                                    "output_name ="])
except getopt.GetoptError:
   print ('path_detect.py -f <trajectoryfile> -d <distance.xvg> -n <ligand_name> -t <distance_interval> -o <output_name>')
   sys.exit(2)                                                                      # Exiting the program raises a SystemExit exception, 0 means normal exit, and others are abnormal exits.
 
for opt, arg in opts:                                                               # Generate several pairs of value, e.g: opr,arg = -i,PLDXPAL
   if opt == '-h':
      print ('path_detect.py -f <trajectoryfile> -d <distance.xvg> -n <ligand_name> -t <distance_interval> -o <output_name>')
      sys.exit()
   elif opt in ("-f", "--trajectory_file"):
      last_200_pdb = arg
   elif opt in ("-d", "--distance_xvg"):
      distance_xvg = arg
   elif opt in ("-t", "--delta_distance"):
      delta_distance = float(arg)
   elif opt in ("-n", "--ligand_name"):
       ligand_name = arg
   elif opt in ("-o", "--output_name"):
       output_name = arg
      
# os.system('gmx')
print("Hello World!")

class tauRAMD_PATH():
    time = []
    distance = []
    picked_time = []
    picked_distance = []
    picked_frame_index = []
    target_ligands = []                                                             # 包括所有目标ligand的第一个原子的序号
    ave_x_target_ligands= []
    ave_y_target_ligands= []
    ave_z_target_ligands= []


    atomic_index = []                                                               # each empty list here is to used as the charactors for object.
    atomic_name = []
    residue_name = []
    chain_name = "A"
    residue_index = []
    X_peratom = []
    Y_peratom = []
    Z_peratom = []
    bfactor_per_factor = []
    charge_per_factor = []
    Atomtype_per_atom = []
    
    def __init__(self,distance_xvg,delta_distance,last_200_pdb,ligand_name,output_name,max_delta_distance):
        
        self.distanceReader(distance_xvg)
        self.find_path(delta_distance,max_delta_distance)
        self.PDBreader(last_200_pdb)
        self.ligand_finder(ligand_name)
        self.PDBwriter(output_name,ligand_name)
    
    def distanceReader(self,distance_file):
        self.time.clear()
        self.distance.clear()
        f = open(distance_file, "r")                                                # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
        for line in f:                                                              # iterate each line in file "f"                     
                if(line.split()[0] != "@" and line.split()[0] != "@TYPE" and line.split()[0] != "#"):               # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
                    self.time.append(float(line.split()[0]))                        # The second column is the atomic number
                    self.distance.append(float(line.split()[1]))                          # The 3rd column is the atom name C CA CD1 CD2 and so on
        f.close()
        
    def find_path(self,interval,max_interval):
        self.picked_distance.clear()
        self.picked_time.clear()
        self.picked_frame_index.clear()
        self.picked_time.append(self.time[0])
        self.picked_distance.append(self.distance[0])
        self.picked_frame_index.append(0)
        for i in range (0 ,len(self.time)):
            if(self.distance[i] - self.picked_distance[-1] >= float(interval) and self.distance[i] - self.picked_distance[-1] <= float(max_interval)):
                self.picked_time.append(self.time[i])
                self.picked_distance.append(self.distance[i])
                self.picked_frame_index.append(i)
                
    def PDBreader(self,filename):
        # aa = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
        self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors
        self.atomic_index.clear()
        self.atomic_name.clear()
        self.residue_name.clear()
        # self.chain_name.clear()
        self.residue_index.clear()
        self.X_peratom.clear()
        self.Y_peratom.clear()
        self.Z_peratom.clear()
        self.bfactor_per_factor.clear()
        self.charge_per_factor.clear()
        self.Atomtype_per_atom.clear()
        f = open(filename, "r")                                                     # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
        # for line in f:                                                              # iterate each line in file "f"                     
        #         if(line.split()[0] in["ATOM","HETATM"] and line.split()[2] not in ["CL","NA","K"]):       # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
        #             self.atomic_index.append(float(line.split()[1]))                # The second column is the atomic number
        #             self.atomic_name.append(line.split()[2])                        # The 3rd column is the atom name C CA CD1 CD2 and so on
        #             self.residue_name.append(line.split()[3])                       # Column 4 is the residue name TYR ALA etc.
        #             # self.chain_name.append(line.split()[4])                         # The 5th column is the name of the chain it is on
        #             self.residue_index.append(float(line.split()[4]))               # The sixth column is the residue number
        #             #self.X_peratom.append(float(re.split('-| |', string)))
        #             self.X_peratom.append(float(line.split()[5]))                   # Column 7 is the x-coordinate of the atom
        #             self.Y_peratom.append(float(line.split()[6]))                   # The 8th column is the Y-coordinate of the atom
        #             self.Z_peratom.append(float(line.split()[7]))                   # The ninth column is the Z-coordinate of the atom
        #             self.bfactor_per_factor.append(float(line.split()[8]))          # The 10th column is the B-factor of the atom, which is used to judge the activity level
        #             self.charge_per_factor.append(float(line.split()[9]))          # Column 11 is the charge of the residue
        #             try:
        #                 self.Atomtype_per_atom.append(line.split()[10])                 # Column 12 is the atomic type of the atom C, H, O, N, S, CL, etc.
        #             except:
        #                 self.Atomtype_per_atom.append("l")
        #             #print(line)               
                   
        for line in f:                                                              # iterate each line in file "f"                     
                if(line.split()[0] in["ATOM","HETATM"] and line.split()[2] not in ["CL","NA","K"]):       # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
                    self.atomic_index.append(float(line.split()[1]))                # The second column is the atomic number
                    self.atomic_name.append(line.split()[2])                        # The 3rd column is the atom name C CA CD1 CD2 and so on
                    self.residue_name.append(line.split()[3])                       # Column 4 is the residue name TYR ALA etc.
                    # self.chain_name.append(line.split()[4])                         # The 5th column is the name of the chain it is on
                    self.residue_index.append(float(line.split()[4]))               # The sixth column is the residue number
                    # #self.X_peratom.append(float(re.split('-| |', string)))
                    # self.X_peratom.append(float(line.split()[5]))                   # Column 7 is the x-coordinate of the atom
                    # self.Y_peratom.append(float(line.split()[6]))                   # The 8th column is the Y-coordinate of the atom
                    # self.Z_peratom.append(float(line.split()[7]))                   # The ninth column is the Z-coordinate of the atom
                    # self.bfactor_per_factor.append(float(line.split()[8]))          # The 10th column is the B-factor of the atom, which is used to judge the activity level
                    # self.charge_per_factor.append(float(line.split()[9]))          # Column 11 is the charge of the residue
                    # record_type = line[0:6].strip()  # 记录类型（ATOM/HETATM）
                    # self.atomic_index = int(line[6:11].strip())  # 原子序号
                    # self.atomic_name = line[12:16].strip()  # 原子名称
                    # alt_loc = line[16].strip()  # 替代位置指示符
                    # residue_name = line[17:20].strip()  # 残基名称
                    # chain_id = line[21].strip()  # 链标识符
                    # residue_num = int(line[22:26].strip())  # 残基序号
                    # insert_code = line[26].strip()  # 插入代码
                    x = float(line[30:38].strip())  # X坐标
                    self.X_peratom.append(x)
                    y = float(line[38:46].strip())  # Y坐标
                    self.Y_peratom.append(y)
                    z = float(line[46:54].strip())  # Z坐标
                    self.Z_peratom.append(z)
                    occupancy = float(line[54:60].strip()) if line[54:60].strip() else 0.0  # 占用率
                    self.bfactor_per_factor.append(occupancy)
                    temp_factor = float(line[60:66].strip()) if line[60:66].strip() else 0.0  # 温度因子
                    self.charge_per_factor.append(temp_factor)
                    try:
                        element = line[76:78].strip()  # 元素符号
                        self.Atomtype_per_atom.append(element)
                    except:
                        self.Atomtype_per_atom.append("l")
                    # charge = line[78:80].strip()  # 电荷
                    # try:
                    #     self.Atomtype_per_atom.append(line.split()[10])                 # Column 12 is the atomic type of the atom C, H, O, N, S, CL, etc.
                    # except:
                    #     self.Atomtype_per_atom.append("l")
                    #print(line)
                    
    def ligand_finder(self,ligandName):
        number_of_atoms_in_ligand = int(self.residue_name.count(ligandName)/len(self.time))
        self.target_ligands.clear()
        c = 0
        for i in range(0, len(self.atomic_index)):
            if(self.residue_name[i] == ligandName and self.residue_name[i-1] != ligandName and c in self.picked_frame_index):
                self.target_ligands.append(i)
                b = i + number_of_atoms_in_ligand
                self.ave_x_target_ligands.append(mean(self.X_peratom[i:b]))
                self.ave_y_target_ligands.append(mean(self.Y_peratom[i:b]))
                self.ave_z_target_ligands.append(mean(self.Z_peratom[i:b]))
                c += 1
            elif(self.residue_name[i] == ligandName and self.residue_name[i-1] != ligandName):
                c += 1
                
    def PDBwriter(self,filename,ligandName):
        number_of_atoms_in_protein = int((len(self.atomic_index)-self.residue_name.count(ligandName))/len(self.time))
        f = open(filename, "w")                                                             # e.g: f = linesplit[0]+"_PO3.pdb"
        for i in range (0 ,len(self.ave_x_target_ligands)):                                         # Create a loop, i is a sequence starting from 0, and the number of atoms is the length  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # Formatted output, %4s, right-aligned, the output occupies 4 columns in total. If the length is less than 4 columns, the left end will be filled with spaces. If it is greater than 4 columns, the actual length will be output as a string
                                             i+1,                          # %7d, right-aligned, the output occupies a total of 7 columns, if the length is less than 7 columns, the left end is filled with spaces, signed decimal certificate integer
                                             "O1",                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                             ligandName,                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                             self.chain_name,                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                             i+1,                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                             self.ave_x_target_ligands[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.ave_y_target_ligands[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.ave_z_target_ligands[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             0,                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             0,
                                             "O"), file = f )                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal


        for i in range(0,number_of_atoms_in_protein):
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM",
                                                                           self.atomic_index[i]+len(self.ave_x_target_ligands),
                                                                           self.atomic_name[i],                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                                                           self.residue_name[i],                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                                                           self.chain_name,                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                                                           self.residue_index[i]+len(self.ave_x_target_ligands),                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                                                           self.X_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                                                           self.Y_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                                                           self.Z_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                                                           self.bfactor_per_factor[i],                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                                                           self.charge_per_factor[i],
                                                                           self.Atomtype_per_atom[i]), file = f )
                                                                           
        print("END", file = f)
        f.close()        
            
                
x = tauRAMD_PATH(distance_xvg, delta_distance, last_200_pdb, ligand_name, output_name, max_delta_distance)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    

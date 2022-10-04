printf "c\nc\n" | gmx trjcat -f 1_md.xtc 2_md_ramd.xtc -settime
printf "1\n 17\n" | gmx trjconv -s 0_md.tpr -f 3_trajout.xtc -center -pbc nojump -o 4_md_noPBC.pdb			# The md.tpr is grab from previous 100ns conventional MD simulation
printf "1\n 0\n" | gmx trjconv -s 4_md_noPBC.pdb -f 4_md_noPBC.pdb -fit rot+trans -o 5_md_noPBC_fit.pdb
c=`more md_noPBC_fit.pdb | grep TITLE | tail -200 | awk '{print $6}' | head -1`; printf "0\n" | gmx trjconv -s 5_md_noPBC_fit.pdb -f 5_md_noPBC_fit.pdb -b $c -o 6_last200_skip2.pdb -skip 2
gmx distance -s 6_last200_skip2.pdb -f 6_last200_skip2.pdb -select 'com of group 1 plus com of group 16' -oav 7_distance_pro_mol.xvg
python $HOME/workspace/bash_script/tauRAMD_path_finder/path_detect.py -f 6_last200_skip2.pdb -d 7_distance_pro_mol.xvg  -n MOL -t 0.2 -o path.pdb 

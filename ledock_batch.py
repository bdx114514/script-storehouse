import os,sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser

op_path = r"."
list_dir = os.listdir(op_path)
def updateFile(file,old_str,new_str):
    file_data = ""
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if old_str in line:
                line = line.replace(old_str,new_str)
            file_data += line
    with open(file,"w",encoding="utf-8") as f:
        f.write(file_data)
os.system("touch affinity_from_ledock.txt")
list_need = []
for i in range(len(list_dir)):
  n_path = op_path + "/" + str(list_dir[i])
  if os.path.isdir(n_path):
    os.chdir(n_path)          #判断是否为文件夹，是则进入#
    os.system("lepro "+ str(list_dir[i]) +"_protein.pdb") #执行lepro预处理#
    p = PDBParser(PERMISSIVE=1)
    id=str(list_dir[i])+"_pocket"
    file=str(list_dir[i])+"_pocket.pdb"
    s = p.get_structure(id,file)
    list_x=[]
    list_y=[]
    list_z=[]
    for chains in s:
     for chain in chains:
      for residue in chain:
       for atom in residue:
          x = atom.get_vector()[0]
          list_x.append(x)
          y=atom.get_vector()[1]
          list_y.append(y)
          z=atom.get_vector()[2]
          list_z.append(z)
    xmax = str(max(list_x))
    xmin = str(min(list_x))
    ymax = str(max(list_y))
    ymin = str(min(list_y))
    zmax = str(max(list_z))
    zmin = str(min(list_z))            #通过结合位点的pdb坐标获得对接box的坐标#
    updateFile("dock.in","20","50")
    updateFile("dock.in","xmax",xmax)
    updateFile("dock.in","xmin",xmin)
    updateFile("dock.in","ymax",ymax)
    updateFile("dock.in","ymin",ymin)
    updateFile("dock.in","zmax",zmax)
    updateFile("dock.in","zmin",zmin)  #将对接box坐标更新到参数文件里,并修改生成poses数量为50#
    os.system("ls *.mol2 > ligands")
    os.system("ledock dock.in")    
    with open(str(list_dir[i])+"_ligand.dok", "r", encoding="utf-8") as x:
       line = x.readline()
       counts =1
       while line:
        if counts>=2:
         break
        line=x.readline()
        score = line[line.rfind('Score'):]
        counts+=1
        score_t = str(list_dir[i]) + "|" + score
        list_need.append(score_t)        #每组对接生成的poses会根据affinity从大到小排序，将每组对接的第一个亲和力保存在列表里#
    print("redocking of "+str(list_dir[i])+" has completed")
    print(score_t)
    os.chdir(r"..") #返回上级菜单#
aaa = np.array(list_need)   
np.savetxt('./affinity_from_ledock.txt',list_need,fmt='%s') #输出结果#  
  
  

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:17:21 2021

@author: Mohan
"""

'''
pass genealogy file as first argument, all_nodes_haplotypes as second argument, and founder ID as third argument
    E.g.:
>python ./IBD_traceback.py ./new index genealogy.csv ./all_nodes_haplotypes.txt ./proband_haplotypes.txt 100
'''

import sys
import csv


def genealogy_dict(gen_file):
    gen={}
    with open(gen_file, mode='r') as genealogy:
        csv_gen = csv.reader(genealogy)    
        header=next(csv_gen)
        
        if header[1]=="mother" and header[2]=="father":
            mother=1
            father=2
        if header[2]=="mother" and header[1]=="father":
            mother=2
            father=1    
            
        for row in csv_gen:
            gen[row[0]]=(row[father],row[mother])
    return(gen)

# class hap:
#     def __init__(self, nRec, XO_start, RecPos):
#         self.nRec      = int(nRec)
#         self.XO_start  = int(XO_start)
#         self.RecPos    = [int(x) for x in RecPos]

def create_haps(an_file, num_nodes):
    hap_dict = {}
    for i in range(num_nodes):
        line = next(an_file)
        line  = line.rstrip().split("}{")

        temp  = line[0].split(";")
        ID    = temp[1]
        temp1 = temp[2].split(",")
        temp2 = temp[3].split(",")  

        if int(temp1[0]) == 0:
            hap_dict[ID+".1"] = hap(temp1[0],temp1[1],[])
        else:
            hap_dict[ID+".1"] = hap(temp1[0],temp1[1],temp1[2:])
        if int(temp2[0]) == 0:
            hap_dict[ID+".2"] = hap(temp2[0],temp2[1],[])
        else:
            hap_dict[ID+".2"] = hap(temp2[0],temp2[1],temp2[2:])
    return(hap_dict)

def traceback(gen, haps, start, position):
    multi = False
    path = [] 
    ah = start #active haplotype name

    if ah[-2:]   == ".1":
        next_node = gen[ah[:-2]][0]
    if ah[-2:]   == ".2":
        next_node = gen[ah[:-2]][1]

    while next_node != '335':

        a_hap = haps[ah]

        if ah[-2:]   == ".1":
            path.append(next_node)
            next_node = gen[ah[:-2]][0]
            
            if a_hap.nRec == 0:
                if a_hap.XO_start == 0:
                    ah=next_node+".1"
                elif a_hap.XO_start == 1:
                    ah=next_node+".2"

            elif a_hap.nRec > 0:
                crossover_count = 0
                for recpos in a_hap.RecPos:
                    if recpos <= position[0]:
                        crossover_count += 1
                    elif recpos > position[0] and recpos<position[1]:
                        multi=True
                        # position = (position[0],recpos)
                        break
                # print(crossover_count)
                if a_hap.XO_start == 0:
                    if crossover_count%2==0:
                        ah=next_node+".1"
                    elif crossover_count%2==1:
                        ah=next_node+".2"                            
                elif a_hap.XO_start == 1:
                    if crossover_count%2==0:
                        ah=next_node+".2"
                    elif crossover_count%2==1:
                        ah=next_node+".1" 

        elif ah[-2:]   == ".2":
            path.append(next_node)
            next_node = gen[ah[:-2]][1]
            
            if a_hap.nRec == 0:
                if a_hap.XO_start == 0:
                    ah=next_node+".1"
                elif a_hap.XO_start == 1:
                    ah=next_node+".2"

            elif a_hap.nRec > 0:
                crossover_count = 0

                for recpos in a_hap.RecPos:
                    if recpos <= position[0]:
                        crossover_count += 1
                    elif recpos > position[0] and recpos<position[1]:
                        multi=True
                        # position = (position[0],recpos)
                        break
                # print(crossover_count)
                if a_hap.XO_start == 0:
                    if crossover_count%2==0:
                        ah=next_node+".1"
                    elif crossover_count%2==1:
                        ah=next_node+".2"                            
                elif a_hap.XO_start == 1:
                    if crossover_count%2==0:
                        ah=next_node+".2"
                    elif crossover_count%2==1:
                        ah=next_node+".1"                 
    
    path = tuple(path)
    return(path, multi)

def check_adj(list_seg): 
    right_boundaries = [x[1] for x in list_seg]
    left_boundaries  = [x[0] for x in list_seg]
    adjacent = [x for x in left_boundaries if x in right_boundaries]

    flag = False
    if len(adjacent) > 0 :
        flag = True

    while len(adjacent) > 0:
        replace = []
        position = adjacent[0]
        # print(position)
        index1  = right_boundaries.index(position)
        index2  =  left_boundaries.index(position)
        seg1 = list_seg[index1]
        seg2 = list_seg[index2]
        replace.append( (seg1[0],seg2[1]) )
        # print(list_seg)

        list_seg.remove(seg1)
        list_seg.remove(seg2)
        list_seg.extend(replace)

        right_boundaries = [x[1] for x in list_seg]
        left_boundaries  = [x[0] for x in list_seg]
        adjacent = [x for x in left_boundaries if x in right_boundaries]

    return(list_seg, flag)

def read_pro_file(numPro, pfile):
    line = next(pfile)
    line = line.rstrip()
    temp = line.split("}{")
    ID=temp[0].split(";")[1]  
    hap1=temp[1].lstrip("{").split(";")
    hap2=temp[2].rstrip("}").split(";")  
    origins_ch1 = hap1[1::2]
    origins_ch2 = hap2[1::2]

    pos = {"ch1":[], "ch2":[]}
    if "335.1" in origins_ch1 or "335.2" in origins_ch1:
        for i,x in enumerate(origins_ch1):
            if x=="335.1" or x=="335.2":
                pos["ch1"].append( (int(hap1[i*2]), int(hap1[2 + i*2])) )

    if "335.1" in origins_ch2 or "335.2" in origins_ch2:
        for i,x in enumerate(origins_ch2):
            if x=="335.1" or x=="335.2":
                pos["ch2"].append( (int(hap2[i*2]), int(hap2[2 + i*2])) )

    return(pos["ch1"],pos["ch2"])

if __name__ == "__main__" :
    genealogy_file  = sys.argv[1]
    all_nodes_haplo = sys.argv[2]
    proband_haplo   = sys.argv[3]
    proID           = sys.argv[4]
    founderID       = sys.argv[5]
    gen = genealogy_dict(genealogy_file)

    pfile   = open(proband_haplo)
    afile   = open(all_nodes_haplo)
    line = next(pfile)
    line = line.split(";")
    simulNo = int(line[0])
    numPro  = int(line[1])
    line    = next(afile)
    numInd  = int(line.split(";")[2])

    rfile   = open("traceback.csv", "w")
    csv_results = csv.writer(rfile)
    csv_results.writerow(["simulNo","path", "Lpos","Rpos"])

    pathlist = []
    for i in range(simulNo):
        print(i)
        hap_dict = create_haps(afile, numInd)
        t1,t2    = read_pro_file(numPro, pfile)
        if len(t1) > 0:
            for t_seg in t1:
                path, multi = traceback(gen,hap_dict,"222.1",t_seg)
                if multi:
                    csv_results.writerow([i+1,0,t_seg[0],t_seg[1]])
                else:
                    if not path in pathlist:
                        pathlist.append(path)
                        csv_results.writerow([i+1,len(pathlist),t_seg[0],t_seg[1]])
                    else:
                        csv_results.writerow([i+1,pathlist.index(path)+1,t_seg[0],t_seg[1]])
        if len(t2) > 0 :
            for t_seg in t2:
                path, multi = traceback(gen,hap_dict,"222.2",t_seg)
                if multi:
                    csv_results.writerow([i+1,0,t_seg[0],t_seg[1]])
                else:
                    if not path in pathlist:
                        pathlist.append(path)
                        csv_results.writerow([i+1,len(pathlist),t_seg[0],t_seg[1]])
                    else:
                        csv_results.writerow([i+1,pathlist.index(path)+1,t_seg[0],t_seg[1]])

    pfile.close()
    afile.close()
    rfile.close()
    print(pathlist)
        
def check_overlaps(interval1, interval2):
    left_pos1  = interval1[0]
    right_pos1 = interval1[1]
    left_pos2  = interval2[0]
    right_pos2 = interval2[1]
    
    #find if the segments overlap, and if they do identify the bounds of the overlapping region
    if left_pos1 <= left_pos2 and right_pos1 > left_pos2:
        if right_pos1 >= right_pos2:
            return(left_pos2, right_pos2)
        else:
            return(left_pos2, right_pos1)
    
    elif left_pos2 <= left_pos1 and right_pos2 > left_pos1:
        if right_pos1 >= right_pos2:
            return(left_pos1, right_pos2)
        else:
            return(left_pos1, right_pos1)   

    else:
        return False

#we store the segment positions in lists of tuples, we wan't pairwise comparison between each element of the two lists of positions and return all overlaps
def check_overlap_list(list1, list2): 
    return_list = [] #list of tuples (left bound, right bound) that store all the overlapping positions 
    for i in list1:
        for j in list2:
            overlap = check_overlaps(i,j)
            
            if overlap:
                return_list.append(overlap)
                   
    return(return_list)

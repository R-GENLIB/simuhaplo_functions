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

def traceback(gen, haps, founderID, start, position):
    multi = False
    path = [] 
    ah = start #active haplotype name

    if ah[-2:]   == ".1":
        next_node = gen[ah[:-2]][0]
    if ah[-2:]   == ".2":
        next_node = gen[ah[:-2]][1]

    while next_node != founderID:

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

def read_pro_file(numPro, pfile, proID, founderID):
    for i in range(numPro):

        line = next(pfile)
        line = line.rstrip()
        temp = line.split("}{")
        ID=temp[0].split(";")[1]  
        if ID = proID:
            hap1=temp[1].lstrip("{").split(";")
            hap2=temp[2].rstrip("}").split(";")  
            origins_ch1 = hap1[1::2]
            origins_ch2 = hap2[1::2]

            pos = {"ch1":[], "ch2":[]}
            if founderID+".1" in origins_ch1 or founderID+".2" in origins_ch1:
                for i,x in enumerate(origins_ch1):
                    if x==founderID+".1" or x==founderID+".2":
                        pos["ch1"].append( (int(hap1[i*2]), int(hap1[2 + i*2])) )

            if founderID+".1" in origins_ch2 or founderID+".2" in origins_ch2:
                for i,x in enumerate(origins_ch2):
                    if x==founderID+".1" or x==founderID+".2":
                        pos["ch2"].append( (int(hap2[i*2]), int(hap2[2 + i*2])) )

    return(pos["ch1"],pos["ch2"])

if __name__ == "__main__" :
    genealogy_file  = sys.argv[1]
    all_nodes_haplo = sys.argv[2]
    proband_haplo   = sys.argv[3]
    proID           = str(sys.argv[4])
    founderID       = str(sys.argv[5])
    gen = genealogy_dict(genealogy_file)

    pfile   = open(proband_haplo)
    afile   = open(all_nodes_haplo)
    line = next(pfile)
    line = line.split(";")
    simulNo = int(line[0])
    numPro  = int(line[1])
    line    = next(afile)
    numInd  = int(line.split(";")[2])

    rfile   = open("traceback.csv", "w", newline='')
    csv_results = csv.writer(rfile)
    csv_results.writerow(["simulNo","path", "Lpos","Rpos"])

    pathlist = []
    for i in range(simulNo):
        print(i)
        hap_dict = create_haps(afile, numInd)
        t1,t2    = read_pro_file(numPro, pfile, proID, founderID)
        if len(t1) > 0:
            for t_seg in t1:
                path, multi = traceback(gen,hap_dict,founderID, proID+".1",t_seg)
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
                path, multi = traceback(gen,hap_dict,founderID, proID+".2",t_seg)
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
    
    o_pathlist = open("pathlist.txt","w")
    for path in pathlist:
        o_pathlist.write(",".join(path))
    o_pathlist.close()

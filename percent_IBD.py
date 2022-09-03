import sys
import csv
import os


def get_probandID(fileobj):
    first_line = next(fileobj)
    first_line = first_line.strip()
    first_line = first_line.split(";")
    
    simulNo = int(first_line[0])
    probandNo = int(first_line[1])
    # Morgan_Len = float(first_line[2])
    
    return simulNo, probandNo

#this function reads 1 simulation worth of results from the proband_haplotype file, and turns it into the dictionary format
#The key;s of the dict are proband ID, the value is another dict with keys 1,2, which are the 2 chromosomes for each individual
#finally, 1 and 2, are a dictionary for both chromosomes which have a key for every segment (segment ID), the corresponding value is a list of tuples. Each tuple is the positions for the segment. List of tuples (left boundary position, right boundary position) since a segment ID can appear in multiple discontinous segments.
# 3 nested dictionaries: [probandID][chromosome][segmentID]
def read_haplos(fileobj, probandNo):
    haplo_dict = {}
    for j in range(probandNo):
        line = next(fileobj)
        line = line.rstrip()
        temp = line.split("}{")
        ID=temp[0].split(";")[1]  
        hap1=temp[1].lstrip("{").split(";")
        hap2=temp[2].rstrip("}").split(";")     
        
        #this list contains the founder of origin for each segment in the probands haplotypes. Both chromosomes are stored in a list
        #theyre in order so we can use the index to get positions
        origins_ch1 = hap1[1::2]
        origins_ch2 = hap2[1::2]
        haplo_dict[ID] = {1:{},2:{}} #store each chromosome seperately

        for index, founderID in enumerate(origins_ch1):
            if founderID not in haplo_dict[ID][1]: #check if key already exists
                haplo_dict[ID][1][founderID] = [(int(hap1[index*2]), int(hap1[2 + index*2]))]
            else: #If another segment later in the chromosome is also from the same founder we just add the position to the list of positions
                haplo_dict[ID][1][founderID].append(  (int(hap1[index*2]), int(hap1[2 + index*2]))  )

        for index, founderID in enumerate(origins_ch2):
            if founderID not in haplo_dict[ID][2]: #check if key already exists
                haplo_dict[ID][2][founderID] = [(float(hap2[index*2]), float(hap2[2 + index*2]))]
            else:
                haplo_dict[ID][2][founderID].append(  (float(hap2[index*2]), float(hap2[2 + index*2]))  )

    return haplo_dict

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

#fore recursive use
def check_overlap_list2(list_pos): #this one returns the positions of the original segments (if they overlap) not the positions of the overlap itself
    return_list = []  
    for i in range(len(list_pos)):
        for j in range(i):
            overlap = check_overlaps(list_pos[i],list_pos[j])
            
            if overlap:
                return_list.append(list_pos[i])
                return_list.append(list_pos[j])
                return(return_list)
    if len(return_list) > 0:
        return(list(set(return_list)))
    else:
        return([])

#if two IBD segments are completely adjacent (right boundary of first segment = left boundary of second) then we merge them into one IBD segment
#if both probands inherit a longer stretch made of multiple contiguous founder segments
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


def combine(list_seg1, list_seg2):
    x = list_seg1 + list_seg2
    # print(list_seg1)
    # print(list_seg2)
    
    overlaps = check_overlap_list2(x)
    while len(overlaps) > 0:
        seg1 = overlaps[0]
        seg2 = overlaps[1]
        x.remove(seg1)
        x.remove(seg2)
        x.append( ( min(seg1[0],seg2[0]), max(seg1[1],seg2[1]) ) )
        overlaps = check_overlap_list2(x)
    return(x)

#checks both chromosomes of a pair of individuals 
def check_IBD_diploid(haplo_pro1, haplo_pro2):
    #this function will check what proportion of proband1's haplotype is IBD to proband2. 

    total_len_IBD =0
    numSegs = 0
    longest_seg = 0

    IBD1_1, _ = check_IBD_single(haplo_pro1[1], haplo_pro2[1])
    IBD1_2, _ = check_IBD_single(haplo_pro1[1], haplo_pro2[2])
    IBD2_1, _ = check_IBD_single(haplo_pro1[2], haplo_pro2[1])
    IBD2_2, _ = check_IBD_single(haplo_pro1[2], haplo_pro2[2])

    IBDp1h1 = combine(IBD1_1, IBD1_2)
    IBDp1h2 = combine(IBD2_1, IBD2_2)
    IBDp2h1 = combine(IBD1_1, IBD2_1)
    IBDp2h2 = combine(IBD1_2, IBD2_2)

    all_seg = IBDp1h1 + IBDp1h2 + IBDp2h1 + IBDp2h2
    if len(all_seg) > 0 :
        all_seg_len = [x[1] - x[0] for x in all_seg]
        numSegs = len(all_seg)
        total_len_IBD = sum(all_seg_len)/4
        longest_seg = max(all_seg_len)
        
    return(total_len_IBD, numSegs, longest_seg)
#compares haploid pair

def check_IBD_single(hap1, hap2):
    segIDs_1 = set(hap1.keys())
    segIDs_2 = set(hap2.keys())
    
    intersection = segIDs_1.intersection(segIDs_2)
    IBD_segs = []
    for segID in intersection:
        IBD_segs.extend(check_overlap_list(hap1[segID], hap2[segID]))

    IBD_segs, flag = check_adj(IBD_segs)
    return(IBD_segs, flag)

def all_pairs_IBD(pfile, results_file):
    csv_results = csv.writer(results_file)
    csv_results.writerow(["simulNo", "pro1", "pro2", "sharedIBD","numSeg","longestSeg"])
    
    simulNo, probandNo = get_probandID(pfile)
    for k in range(simulNo):
        haplo_dict = read_haplos(pfile, probandNo)
        proband_list = list(haplo_dict.keys())
        # print(len(proband_list))
        #iterate through all pairs of probands
        for n in range(len(proband_list)):
            for o in range (n):
                pro1 = proband_list[n]
                pro2 = proband_list[o]
                # print(k,pro1,pro2)
                # print("--------")
            
                shared_IBD, numSeg, longest_seg = check_IBD_diploid(haplo_dict[pro1],haplo_dict[pro2])
                if shared_IBD > 0:
                    csv_results.writerow([k+1, pro1, pro2, shared_IBD, numSeg, longest_seg])

if __name__ == "__main__":

    pfile = sys.argv[1]
    results_file = open("HBD.csv", mode='w', newline='')
    all_pairs_IBD(pfile, results_file)               
    pfile.close()
    results_file.close()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def ecDNAMetrics(gfile, cfile):
    cycle_frac = 0
    t_n_ratio = 0
    bfb_score = 0
    
    s,c  = readConvertedCycles(cfile, gfile)
    frac_table = cycleFractionTable(c,s)
    
    circular_frac = frac_table[frac_table['Circular']]
    if circular_frac.shape[0] > 0:
        largest_cyc_size = circular_frac['Size'].iloc[0]
        cycle_frac = circular_frac['CN Fraction'].iloc[0] #move this chunk to another function
    
    n_trans, cns, all_cns = count_transitions(gfile)
    n_cn_states = len(deduce_states2(all_cns))
    
    t_n_ratio = n_trans/n_cn_states
    
    n_chrs = len(spanningChrs(gfile))
    
    return (cycle_frac, t_n_ratio, bfb_score, n_chrs), n_trans, n_cn_states

def ecDNAContext(metrics, n_trans, n_cn, t_n_cutoff = 4, cycle_cutoff = 0.2, bfb_cutoff = 1):
    if metrics[0] >= cycle_cutoff and (metrics[1] < t_n_cutoff or n_cn <= 1):
        return "Episome"
    
    elif metrics[0] >= cycle_cutoff and metrics[1] >= t_n_cutoff and n_cn > 1:
        return "Both; Check for BFB"
    
    elif metrics[1] >= t_n_cutoff and n_cn > 1:
        n_chrs = metrics[3]
        if n_chrs > 1:
            return "Chromoplexy"
        else:
            return "Chromothripsis"
    else:
        return "Neither; Check for BFB"
        
def spanningChrs(graphf):
    #given an amplicon graph file, get the chromosomes involved in the amplicon
    data = open(graphf).readlines()
    unique_chrs = []
    for line in data:
        if line.startswith("sequence"):
            fields = line.split()
            start_chr = fields[1].split(":")[0]
            end_chr = fields[2].split(":")[0]
            if start_chr not in unique_chrs:
                unique_chrs.append(start_chr)
                
            if end_chr not in unique_chrs:
                unique_chrs.append(end_chr)
    return unique_chrs


def seg_str_int(segs):
    res = []
    for seg in segs:
        i = int(seg[:len(seg)-1])
        if seg[-1] == '-': i *= -1
        res.append(i)
    return res
    
def readCyclef(file):
    pref = file.split("/")[-1].split("amplicon")[0]
    intervalFiles = []
    segments = {}
    cycles = {}
    
    data = open(file).readlines()
    for line in data:
        if line.startswith("Interval"):
            linesplit = line.split()
            intervalFile = pref + linesplit[2] + "_" + linesplit[3] + "_" + linesplit[4] + "_cnseg.txt"
            intervalFiles.append(intervalFile)
        if line.startswith("Segment"):
            ls = line.split()
            segments[int(ls[1])] =  [ls[2], int(ls[3]), int(ls[4])]
        if line.startswith("Cycle"):
            ls = line.split(";")
            cycle_num = int(ls[0].split("=")[1])
            copy_num = float(ls[1].split("=")[1])
            cycle_segs = ls[2].split("=")[1].strip().split(",")
            cycle_segs = seg_str_int(cycle_segs)
            cycles[cycle_num] = [copy_num, cycle_segs]
    return intervalFiles, segments, cycles

def readConvertedCycles(cycles_file, graph_file):
    segments = {}
    cycles = {}
    graph_data = open(graph_file)
    sequences = []
    for line in graph_data:
        if line.startswith("sequence"):
            sequences.append(line.split())
    
    data = open(cycles_file).readlines()
    for line in data:
        if line.startswith("Segment"):
            ls = line.split()
            seg_cn = float(sequences[int(ls[1]) - 1][3])
            segments[int(ls[1])] =  [ls[2], int(ls[3]), int(ls[4]), seg_cn]
        if line.startswith("Cycle"):
            ls = line.split(";")
            cycle_num = int(ls[0].split("=")[1])
            copy_num = float(ls[1].split("=")[1])
            cycle_length = ls[2].split("=")[1]
            cycle_segs = ls[3].split("=")[1].strip().split(";")[0].split(",")
            cycle_segs = seg_str_int(cycle_segs)
            cycles[cycle_num] = [copy_num, cycle_length, cycle_segs]
    return segments, cycles


def get_plausible_frac(cyclef, graphf):
    graph_segments = [] #store size and copy numbers
    graph_data = open(graphf)
    for line in graph_data:
        if line.startswith("sequence"):
            fields = line.split()
            seg_size = int(fields[5])
            seg_cn = float(fields[3])
            graph_segments.append((seg_size, seg_cn))

    max_prop = 0
    cycle_data = open(cyclef)
    max_segs = []
    cycle_cn = 0
    cycle_frac = 0
    for line in cycle_data:
        if line.startswith("Cycle"):
            ls = line.split(";")
            prop = float(ls[2].split("=")[1])
            cycle_cn = float(ls[1].split("=")[1])
            segs = ls[3].rstrip().split("=")[1].split(",")
            if '0+' in segs:
                prop = 0
            
            if prop > max_prop:
                max_prop = prop
                max_segs = segs
    total_weighted_cn = 0
    for seg in graph_segments:
        total_weighted_cn += seg[0]*seg[1]
    if max_prop > 0:
        #get fraction of copy number contributed by segments in this cycle
        #scale factor is given by copy count
        print(max_segs)
        max_segs = seg_str_int(max_segs)
        print(max_segs)
        cycle_size = 0
        for seg in max_segs:
            cycle_size += graph_segments[seg - 1][0]
        
        cycle_frac = cycle_cn*cycle_size/total_weighted_cn
    return max_prop, cycle_frac, cycle_cn
    
def cycleContribution2(cycles, segments):
    #find copy number contribution of overlapping cycles for each region in each interval
    seg_contribs = []
   
    for segID in segments:
        chrom, start, end, seg_cn = segments[segID]
        row = [segID, chrom, start, end, seg_cn]
        total = 0
        circular_total = 0
        for cID in cycles:
            cycle_cn = cycles[cID][0]
            cycle_segs = cycles[cID][1]
            cur_seg_count = 0
            circular = False if cycle_segs[0] == 0 else True
            for cseg in cycle_segs:
                if cseg == segID or cseg == -segID:
                    cur_seg_count += 1
            cycle_share = cur_seg_count*100*cycle_cn/seg_cn
            total += cycle_share
            if circular: circular_total += cycle_share
            row.append(cycle_share)
        row.append(end - start)
        row.append(total)
        row.append(circular_total)
        seg_contribs.append(row)
    return pd.DataFrame(seg_contribs, columns=['Segment', 'Chr','Start','End','Segment_CN'] + [cID for cID in cycles] + ['Size','Total','Circular Total'])
                    
    
def cycleFractionTable(cycles, segments):
    res = []
    total_cn = 0 #get total weighted copy number for the amplicon
    # segments are not overlapping
    for segID in segments:
        chrom, start, end, seg_cn = segments[segID]
        total_cn += seg_cn*(abs(end - start))

    for cycID in cycles:
        cycle_cn, cycle_size, cycle_segs = cycles[cycID]
        circular = False if 0 in cycle_segs else True
        cycle_size_num = int(cycle_size.split("Kbp")[0])*1000
        cycle_frac = cycle_size_num*cycle_cn/total_cn
        row = [cycID, cycle_size,cycle_cn, circular, cycle_frac]
        res.append(row)
    return pd.DataFrame(res, columns=['Cycle','Size','CN','Circular','CN Fraction']).sort_values(by='CN Fraction', ascending=False)
        
def isEpisome(frac_df):
    frac_df_circular = frac_df[frac_df['Circular'] == True]
    
    if (frac_df_circular.shape[0] > 0):
        if frac_df_circular['CN Fraction'].max() > 0.7:
            return True
    return False

def pos_str_int(pos):
    #get -i or +i int from i- or i+
    res = int(pos[:len(pos)-1])
    if pos[-1] == '-':
        res *= -1
    return res

def count_transitions(graphf):
    prev_cn = None
    prev_spos = None
    prev_epos = None
    prev_schr, prev_echr = None, None
    lines = open(graphf).readlines()
    all_cns = []
    cns = []
    n_trans = 0
    for line in lines:
        if line.startswith("sequence"):
            line = line.split()
            schr, spos = line[1].split(":")
            echr, epos = line[2].split(":")
            spos = pos_str_int(spos)
            epos = pos_str_int(epos)
            cn = float(line[3])
            all_cns.append(cn)
            
            if prev_cn:
                if prev_schr == schr and prev_echr == echr and abs(abs(prev_epos) - abs(spos)) <= 2000:
                    #if round(cn) == round(prev_cn):
                    if abs(cn - prev_cn) < 1:
                        prev_epos = epos
                        continue
                    else:
                        n_trans += 1
                       
                        cns.append(round(cn))
                
            prev_cn, prev_spos, prev_epos, prev_schr, prev_echr = cn, spos, epos, schr, echr
    return n_trans, cns, all_cns

def deduce_states2(all_cns):
    all_round = sorted(list(map(lambda x: round(x), all_cns)))
    cur_cns = []
    count_dict = {}
    for i in range(len(all_round)):
        if len(cur_cns) > 0:
            if abs(cur_cns[0] - all_round[i]) > 2:
                count_dict[sum(cur_cns)/len(cur_cns)] = len(cur_cns)
                cur_cns = []
        cur_cns.append(all_round[i])
    count_dict[sum(cur_cns)/len(cur_cns)] = len(cur_cns)
    return count_dict
        
#oscillating copy numbers signify chromothripsis
#write a function that takes in CN file or graph file as input and outputs chains of segments with oscillating copy numbers. 
def oscillating_cn(graphf, show_plot=True):
    prev_cn = None
    prev_spos = None
    prev_epos = None
    prev_schr, prev_echr = None, None
    cn_diff = []
    all_cns = []
    cns = []
    lines = open(graphf).readlines()
    n_seg = 0
    n_seq = 0
    n_seq_arr = []
    if show_plot:
        plt.figure(figsize = (15,7))
        #plt.scatter(range(1,len(cn_diff) + 1, 1),cn_diff)
        print("Oscillations with length > 3:")
        print("#Segs\tMax CN\tMin CN")
        
    for line in lines:
        if line.startswith("sequence"):
            n_seq += 1
            line = line.split()
            schr, spos = line[1].split(":")
            echr, epos = line[2].split(":")
            spos = pos_str_int(spos)
            epos = pos_str_int(epos)
            cn = float(line[3])
            all_cns.append(cn)
   
            if prev_cn:
                
                if prev_schr == schr and prev_echr == echr and abs(abs(prev_epos) - abs(spos)) <= 2000:
                    if round(cn) == round(prev_cn):
               
                        prev_epos = epos
                        continue
                    else:
                        n_seg += 1
                        cn_diff.append(cn - prev_cn)
                        if show_plot:
                            plt.scatter(n_seq, cn - prev_cn, color='blue')
                        cns.append(cn)
                        n_seq_arr.append(n_seq)
                else:
                    if show_plot:
                        plt.axvline(n_seq,color='black' )
           
              
                        
                
            prev_cn, prev_spos, prev_epos, prev_schr, prev_echr = cn, spos, epos, schr, echr
            
    if show_plot:
        plt.xlim(0,n_seq + 1)
        plt.ylim(max(-40, min(cn_diff, default=0) - 2), min(40, max(cn_diff, default=0) + 2))
        
    l = cn_diff
    n_osc = 0
    n_seg_osc = 0
    cur_length = 0
    last_i = 0
    max_cont = 0
    cl_arr = []
    s1 = []
    s2 = []
    for i in range(len(l) - 2):
        if l[i]*l[i+1] < 0 and l[i+1]*l[i+2] < 0 and min(abs(l[i]), abs(l[i+1]), abs(l[i+2])) >= 1:
            if abs(1 - abs(l[i+1])/abs(l[i])) < 0.1  and abs(1 - abs(l[i+2])/abs(l[i+1])) < 0.1:
                #if abs(1 - abs(l[i+1])/abs(l[i])) > 0.1 or abs(1 - abs(l[i+2])/abs(l[i+1])) > 0.1:
                    #print("check this")
                    #print(l[i],l[i+1], l[i+2]) 
                if last_i > 0  and last_i == i + 1:
                    cur_length += 1
                    last_i = i+2
                    n_seg_osc += 1
                else:
                
                    if show_plot:
                        if cur_length > 3:
                            print(str(cur_length) + "\t" + str(round(max(cns[i],cns[i+1],cns[i+2]),2))\
                                  + "\t" + str(round(min(cns[i],cns[i+1],cns[i+2]),2)))
                    cl_arr.append(cur_length)
                    s1.append(max(cns[i],cns[i+1],cns[i+2]))
                    s2.append(min(cns[i],cns[i+1],cns[i+2]))
                    if cur_length > max_cont:
                        max_cont = cur_length
                    n_seg_osc += 3
                    cur_length = 3
                    last_i = i+2
                
                n_osc += 1
                if show_plot:
                    plt.scatter(n_seq_arr[i], l[i],color='red',marker='x')
                    plt.xlabel("Segment Number")
                    plt.ylabel("Copy Number Change")
                    plt.title("CN changes and Oscillations")
                    
    
  

    def deduce_states():
        cn_avg_states = []
        for i in range(len(s1)):
           
            if cl_arr[i] > 3:
                if len(cn_avg_states) == 0:
                    if s1[i] > 2 and s1[i] - s2[i] >= 1:
                        cn_avg_states.append(s1[i])
                        cn_avg_states.append(s2[i])
                else:
                    if s1[i] > 2 and s1[i] - s2[i] >= 1:
                        for s in [s1[i],s2[i]]:
                            found = False
                            for j in range(len(cn_avg_states)):
                                if abs(cn_avg_states[j] - s) < 2:
                                    found = True
                            if found:
                                min_j = None
                                cur_min = 1000
                                for j in range(len(cn_avg_states)):
                                    if abs(cn_avg_states[j] - s) < cur_min:
                                        cur_min = abs(cn_avg_states[j] - s)
                                        min_j = j
                                        
                                cn_avg_states[min_j] += s
                                cn_avg_states[min_j]/= 2
                            if not found:
                               
                                cn_avg_states.append(s)
        return cn_avg_states
    
    states = sorted(deduce_states(),reverse=True)
    
    state_dict = {}
    
    for cn in all_cns:
        
        for s in states:
            if abs(cn - s) < 1:
                state_dict.setdefault(s,0)
                state_dict[s] += 1
            
                
    state_dict = {k: v for k, v in sorted(state_dict.items(), key=lambda item: item[1], reverse=True)}
    
    if show_plot:
        plt.show()
                            
    return max_cont, n_osc, n_seg, n_seg_osc, state_dict, len(all_cns)


def findOverlapSegs(region, segments):
    #given segments as dictionary, find segmentIDs that
    #have some overlap with the region given by chri:start->end
    region = region.split(":")
    chromosome = region[0]
    regStart = int(region[1].split("->")[0])
    regEnd = int(region[1].split("->")[1])
    
    overlapSegs = []
    for seg in segments:
        if segments[seg][0] == chromosome:
            seg_start = segments[seg][1]
            seg_end = segments[seg][2]
            if (seg_start >= regStart and seg_start <= regEnd) or (seg_end >= regStart and seg_end <= regEnd)\
                or (seg_start <= regStart and seg_end >= regEnd):
                overlapSegs.append(seg)
    return overlapSegs

def readIntervalFile(file):
    regions = {}
    data = open(file).readlines()[1:]
    for line in data:
        ls = line.split()
        reg = ls[0] + ":" + ls[1] + "->" + ls[2]
        cpcount = float(ls[3])
        regions[reg] = cpcount
    return regions

def getOverlap(segment, region):
    seg_start = segment[1]
    seg_end = segment[2]
    region_chr, region = region.split(":")[0], region.split(":")[1]
    if segment[0] != region_chr: 
        return 0
    frac = 0
    region = region.split("->")
    reg_start, reg_end = int(region[0]), int(region[1])
    if seg_end <= reg_end and seg_end >= reg_start:
        frac = (seg_end - reg_start)/(reg_end - reg_start)
    elif seg_start >= reg_start and seg_start <= reg_end:
        frac = (reg_end - seg_start)/(reg_end - reg_start)
    elif seg_end <= reg_end and seg_start >= reg_start:
        frac = (seg_end - seg_start)/(reg_end - reg_start)
    elif seg_end >= reg_end and seg_start <= reg_start:
        frac = 1
    return frac

def findCycleContribution(cycle_file, dir_path = "CCLE_AA_outputs/"):
    #find copy number contribution of overlapping cycles for each region in each interval
    intervalFiles, segments, cycles = readCyclef(cycle_file)
    reg_contribs = []
   
    for intFile in intervalFiles:
        f = dir_path + intFile
        regions = readIntervalFile(f)
        
        for reg in regions:
            reg_cn = regions[reg]
            #get list of segments that have some overlap with reg
            overlap_segs = findOverlapSegs(reg, segments)
            print(overlap_segs)
            row = [reg, reg_cn]
            total = 0
            for cycleID in cycles:
                cycle_cn = cycles[cycleID][0]
                cycle_segs = cycles[cycleID][1]
                #print(cycle_segs)
                share = sum([getOverlap(segments[abs(seg)], reg) for seg in cycle_segs if seg != 0])
                #print(share)
                cycle_cn_share = 100*share*cycle_cn/reg_cn
                total += cycle_cn_share
                row.append(cycle_cn_share)
            row.append(total)
            reg_contribs.append(row)
    return pd.DataFrame(reg_contribs, columns = ["Region", "Region_CN"] + [cycleID for cycleID in cycles] + ["Total"])


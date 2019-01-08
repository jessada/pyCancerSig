#!/usr/bin/env python
import sys
import argparse
import os
import numpy
from operator import itemgetter, attrgetter
import readVCF
from collections import defaultdict

# define constant
event_list = []
EVENT_BND = "BND"
event_list.append(EVENT_BND)
EVENT_DEL = "DEL"
event_list.append(EVENT_DEL)
EVENT_DUP = "DUP"
event_list.append(EVENT_DUP)
EVENT_INV = "INV"
event_list.append(EVENT_INV)
log10_list = []
LOG10_2_3 = "log10_2_3"
log10_list.append(LOG10_2_3)
LOG10_3_4 = "log10_3_4"
log10_list.append(LOG10_3_4)
LOG10_4_5 = "log10_4_5"
log10_list.append(LOG10_4_5)
LOG10_5_6 = "log10_5_6"
log10_list.append(LOG10_5_6)
LOG10_6_7 = "log10_6_7"
log10_list.append(LOG10_6_7)
LOG10_7_8 = "log10_7_8"
log10_list.append(LOG10_7_8)
LOG10_8_9 = "log10_8_9"
log10_list.append(LOG10_8_9)
LOG10_9up = "log10_9up"
log10_list.append(LOG10_9up)

def EFF(effects):
    snp_dictionary={}
    for effect in effects:
        if "HIGH" in effect or "MODERATE" in effect and "protein_coding" in effect:
            variant=effect.split("(")[0]
            if "[" in variant:
                variant=variant.split("[")[0]
            gene=effect.split("|")[5]
            feature=effect.split("|")[3]
            #for each snp, each unique effect of each gene shoudl only be reported once, ie not multiple intron variant GENEX for snp z
            if not gene in snp_dictionary:
                snp_dictionary[gene]={variant:feature}
            else:
                snp_dictionary[gene].update({variant:feature})
        #if sequence_feature is not the only entry of a gene, 
    return(snp_dictionary)
    
def ANN(effects):
    snp_dictionary={}
    for effect in effects:
        if "HIGH" in effect or "MODERATE" in effect or "MODIFIER":
            variant=effect.split("|")[1]
            if "[" in variant:
                variant=variant.split("[")[0]
            gene=effect.split("|")[3]
            feature=effect.split("|")[10]
            #for each snp, each unique effect of each gene shoudl only be reported once, ie not multiple intron variant GENEX for snp z
            if not gene in snp_dictionary:
                snp_dictionary[gene]={}
            snp_dictionary[gene][variant]=feature
            #if sequence_feature is not the only entry of a gene, 
    return(snp_dictionary)
    
parser = argparse.ArgumentParser("""extract structural variant profile from FindSV vcf""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--frequency',type=float,default = 0.01,help="frequency threshold, more common variants are not included in the profile")
parser.add_argument('--out',type=str,required=True,help="output file name")
args, unknown = parser.parse_known_args()

# initiate event counter
event_counter = {}
for event in event_list:
    event_counter[event] = defaultdict(int)

# count events
for line in open(args.vcf):
    if not "#" == line[0]:
        chrA, posA, chrB, posB,event_type,INFO,format = readVCF.readVCFLine(line)
        content=line.strip().split("\t")

        if chrB not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]:
            continue

        frq=0
        if "FRQ" in INFO:
            frq=float(INFO["FRQ"])

        if frq > args.frequency:
            continue

        # set default signal strength to be very high to keep result sfrom CNVnator
        signalPE=1000000

        if "PE" in format:
            signalPE =int(format["PE"][0])
        elif  "DV" in format:
            signalPE = int(format["DV"][0])

        if signalPE < 10:
            continue

        length = float("inf")
        if chrA == chrB:
            length = int(posB)-int(posA)

        log_len = numpy.log10(length)
        if length == float("inf"):
            event_counter[event_type][LOG10_9up] += 1
            continue
        if log_len > 9:
            event_counter[event_type][LOG10_9up] += 1
            continue
        if log_len > 8:
            event_counter[event_type][LOG10_8_9] += 1
            continue
        if log_len > 7:
            event_counter[event_type][LOG10_7_8] += 1
            continue
        if log_len > 6:
            event_counter[event_type][LOG10_6_7] += 1
            continue
        if log_len > 5:
            event_counter[event_type][LOG10_5_6] += 1
            continue
        if log_len > 4:
            event_counter[event_type][LOG10_4_5] += 1
            continue
        if log_len > 3:
            event_counter[event_type][LOG10_3_4] += 1
            continue
        if log_len > 2:
            event_counter[event_type][LOG10_2_3] += 1
            continue

# write to output
header = "profile_id"
stat_txt = os.path.splitext(os.path.basename(args.vcf))[0]
if stat_txt.endswith("_FindSV"):
    stat_txt = stat_txt[:-7]
#print(stat_txt, file=sys.stderr)
for event in event_list:
    for log_len in log10_list:
        header += "\t" + event + "_" + log_len
        event_count = event_counter[event][log_len]
        if event_count == 0:
            event_count = 1
        stat_txt += "\t" + str(event_count)


with open(args.out, "w") as f:
    f.write(header+"\n")
    f.write(stat_txt+"\n")
    f.close()

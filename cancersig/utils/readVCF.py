"""
Author: Jesper Eisfeldt, https://github.com/J35P312
"""
import sys,re

def readVCFLine(line):
    if line[0] == "#":
        return(None)

    variation   = line.rstrip().split("\t")
    event_type=""
    chrA=variation[0].replace("chr","").replace("Chr","").replace("CHR","");
    posA=int(variation[1]);
    posB=0;

    description ={}
    INFO=variation[7].split(";");
    for tag in INFO:
        tag=tag.split("=")
        if(len(tag) > 1):
            description[tag[0]]=tag[1];

    format={}
    format_keys={}
    if len(variation) > 8:
        format_string=variation[8].split(":")

        i = 0;
        for key in format_string:
            format_keys[i]=key
            format[key]=[]
            i += 1
        format_fields=variation[9:]
        for sample in format_fields:
            i=0
            format_string= sample.split(":")
            for i in range(0,len(format_keys)):
                format[format_keys[i]].append(format_string[i])
                i += 1
           
    #Delly translocations
    if("TRA" in variation[4]):
        event_type="BND"
        chrB=description["CHR2"]
        posB=int(description["END"]);
        if chrA > chrB:
            chrT= chrA
            chrA = chrB
            
            chrB= chrT
            tmpPos=posB
            posB=posA
            PosA=tmpPos

    #intrachromosomal variant
    elif(not  "]" in variation[4] and not "[" in variation[4]):
        
        chrB=chrA;
        posB=int(description["END"]);
        #sometimes the fermikit intra chromosomal events are inverted i.e the end pos is a lower position than the start pos
        if(posB < posA):
            tmp=posB
            posB=posA
            posA=tmp
        event_type=variation[4].strip("<").rstrip(">");

        if "<" in variation[4] and ">" in variation[4]:
            if "DUP" in event_type:
                event_type="DUP"
            if "INS" in event_type:
                event_type="BND"
        else:
            if "SVTYPE" in description:
                event_type=description["SVTYPE"];
                if "INS" in event_type:
                    event_type = "BND"
    #if the variant is given as a breakpoint, it is stored as a precise variant in the db
    else:
        B=variation[4];

        B=re.split("[],[]",B);
        for string in B:
            if string.count(":"):
                lst=string.split(":");
                chrB=lst[0].replace("chr","").replace("Chr","").replace("CHR","")
                posB=int(lst[len(lst)-1]);
                if chrA > chrB:
                    chrT = chrA
                    chrA = chrB
                    chrB = chrT
                        
                    tmpPos=posB
                    posB=posA
                    posA=tmpPos
                elif chrA == chrB and posA > posB:
                    tmpPos=posB
                    posB=posA
                    posA=tmpPos                   
                
        event_type="BND"
    return( chrA, posA, chrB, posB,event_type,description,format);

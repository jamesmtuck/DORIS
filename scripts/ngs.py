'''
Author: Kevin Volkel
filename: DORIS_NGS.py
Description: Strand analaysis for the DORIS project
'''

###GLOBALS###

reverse_primer={
    "type1":"CAGGTACGCAGTTAGCACTC" #constant across all file-sets
}

payload={ #dictionary used for payloads 
    "type1": "CGTACTGCTCGATGAGTACTCTGCTCGACGAGATGAGACGAGTCTCTCGTAGACGAGAGCAGACTCAGTCATCGCGCTAGAGAGCA", #conanical payload1 5-bases (used in file sets 1-3), G at end is from promoter
    "type2": "ACTGCTCGATGAGTACTCTGCTCGACGAGATGAGACGAGTCTCTCGTAGACGAGAGCAGACTCAGTCATCGCGCTAGAGAGCATAGAGTCGTG", #conanical payload 3 bases (used in file sets 1-3)
    "type3": "CGTACTGCTCGATGAGTACTCTGCTCGACGAGATGAGACGAGTCTCTCGTAGACGAGAGCA" #payload for file-set 4
    }


payload_to_forward={#sequences between the payload and forward sequences, used for file-set 1
    "type1":"GCGCGCTATAGTGAGTCGTATTANNNNN",
    "type2":"GCGCGCTATAGTGAGTCGTATTA",
    "type3":""#empty--> not used
    }


forward_primer={#dictionary used for forward primers
    "type1": "TCCGTAGTCATATTGCCACG", #conanical forward primer, set-1 files
    "type2": "GCGCGCAAAAAAAAAAAAAA", #forward primer for set-2 files
    "type3": "GCGCGCAGTCAGATCAGCTATTACTG", #forward primer for set-3 files
    "type4": "GACTCAGTCATCGCGCTAG" #forward primer for set-4 files
}




'''
complete strand = reverse_primer+*5-barcode*+*payload*+*3-barcode*+payload_to_forward+forward_primer

--barcode will depend on what sequence occurs after the reverse_primer. Test if there is an indiviual base G: indicates 5 base barcode. 
Or test for a short sequence after the primer: ACTGCT indicates a 3 base barcode.

--there should be two payload versions for each set corresponding to each barcode

'''



file_strand_archs={ #holds characteristics for each file studied
    "dsLi": {
        "reverse_primer":"type1",
        "forward_primer":"type1",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type1",
        "payload_to_forward_2":"type2",
        "set": 1
    },
    "hyLi": {
        "reverse_primer":"type1",
        "forward_primer":"type1",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type1",
        "payload_to_forward_2":"type2",
        "set": 1
    },
    "Pool_Stock": {
        "reverse_primer":"type1",
        "forward_primer":"type1",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type1",
        "payload_to_forward_2":"type2",
        "set": 1
    },
    "dsLi_cDNA_Tailed": {
        "reverse_primer":"type1",
        "forward_primer":"type2",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 2
    },
    "hyLi_cDNA_Tailed": {
        "reverse_primer":"type1",
        "forward_primer":"type2",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 2
    },
    "hyLi_cDNA_GC": {
        "reverse_primer":"type1",
        "forward_primer":"type3",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 3
    },
    "dsLi_cDNA_GC": {
        "reverse_primer":"type1",
        "forward_primer":"type3",
        "payload_1":"type1",
        "payload_2":"type2",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 3
    },
    "hyLi_cDNA_26": {
        "reverse_primer":"type1",
        "forward_primer":"type4",
        "payload_1":"type3",
        "payload_2":"type3",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 4
    },
    "dsLi_cDNA_26": {
        "reverse_primer":"type1",
        "forward_primer":"type4",
        "payload_1":"type3",
        "payload_2":"type3",
        "payload_to_forward_1":"type3",
        "payload_to_forward_2":"type3",
        "set": 4
    }
}
    
    
##############END GLOBALS#############################



reverse_table={'A':'T','T':'A','G':'C','C':'G','N':'N', '\x00':'N'}
def calculate_reverse_compliment(strand):
    reverse=strand[::-1]
    #print strand
    return ''.join([reverse_table[nuc] for nuc in reverse])


def convertBaseHelper(base,dec,s):
    bases = ['A', 'C', 'G', 'T']
    m = int(dec % base)
    q = int(dec / base)
    #print(s)
    s = s + bases[m]
    #print(s)
    if q > 0:
        return convertBaseHelper(base,q,s)
    else:
        return s
    
def convertBase(base,dec,length):
    bases=['A','C','G','T']
    s = convertBaseHelper(base,dec,'')
    s = s.ljust(length,bases[0])
    return s
    
def get_experiment_base_name(experiment_file):
    exper_base=experiment_file.split("_rep_")[0]
    assert not exper_base == ""
    if "rep" in exper_base:
        exper_base=exper_base.split("-rep-")[0]
        assert not exper_base == ""
    return exper_base



barcode_to_index_dict={}
def init_barcode_to_index_dict():
    #for however many barcodes we have (1088) generate the corresponding sequence
    #mapping is simply a base 4 mapping
    global barcode_to_index_dict
    index_array=range(0,1088)
    base_3_array=[]
    base_5_array=[]
    for i in range(0,64):
        barcode_to_index_dict[convertBase(4,i,3)[::-1]]=i #reverse it so most sig character is from left to right
        base_3_array.append(convertBase(4,i,3)[::-1])
    for i in range(64,1088):
        barcode_to_index_dict[convertBase(4,i-64,5)[::-1]]=i
        base_5_array.append(convertBase(4,i-64,5)[::-1])
        #print convertBase(4,i-64,5)[::-1]
    assert len(barcode_to_index_dict)==1088
    assert len(base_3_array)==64
    assert len(base_5_array)==1024
    return base_3_array+base_5_array
    
    
#mapping: barcode --> index
def barcode_to_index(barcode):
    global barcode_to_index_dict
    #print barcode_to_index_dict
    return barcode_to_index_dict[barcode]



def analyze_strands(file_base_name,file_key_name,strand_array,data_dictionary):
    global reverse_primer
    global file_strand_archs
    global payload
    #initialize keys used for data collection, key to data dictionary should acount for file name
    count_key=file_key_name+"_"+"count"
    error_rate_key=fq_file+"_"+"error_rate"
    for strand in strand_array:
        if len(strand)<100:
            continue #remove largely short strands
        reverse_p=reverse_primer[file_strand_archs[file_base_name]["reverse_primer"]]
        reverse_p_RC=calculate_reverse_compliment(reverse_p)

        find_FW=strand.find(reverse_p)
        find_RC=strand.find(reverse_p_RC)
        _strand=""
        if not find_FW == -1:
            #found forward strand
            start_point=find_FW
            end_point=min(find_FW+len(reverse_p),len(strand))
            _strand=strand
            #print "forward way start {} end {}".format(start_point,end_point)
            assert end_point>=0
            assert start_point>=0
        elif not find_RC == -1:
            #found the reverse complement
            _strand=calculate_reverse_compliment(strand)
            end_point=min(len(strand)-find_RC,len(_strand)-1)
            start_point=end_point-len(reverse_p)
            #print "reverse way start {} end {}".format(start_point,end_point)
            assert end_point>=0
            assert start_point>=0
        else:
            continue #did not find either

        assert not _strand == ""
        #know the bound points, end_points are non-inclusive, process _strand now
        if _strand[end_point]=='G':
            N_5_payload=payload[file_strand_archs[file_base_name]["payload_1"]]
            N_5_start=min(end_point+6,len(_strand)-1)
            N_5_end=min(end_point+6+len(N_5_payload),len(_strand))
            lv_N_5_payload=lv.distance(N_5_payload,_strand[N_5_start:N_5_end])
            payload_ld=lv.distance(payload[file_strand_archs[file_base_name]["payload_1"]],payload[file_strand_archs[file_base_name]["payload_2"]])
            if lv_N_5_payload < payload_ld/2 or (payload_ld==0 and lv_N_5_payload<7):
                #we have a NNNNN-5 base barcode strand
                barcode_start=end_point+1
                barcode_end=barcode_start+5
                barcode=_strand[barcode_start:barcode_end]
                assert len(barcode)==5
                data_dictionary[count_key][barcode_to_index(barcode)]+=1 #count barcode occurance
            if lv_N_5_payload < payload_ld/2 or (payload_ld==0 and lv_N_5_payload<7):
                edit_operations=lv.editops(N_5_payload,_strand[end_point+6:end_point+6+len(N_5_payload)])
                for op in edit_operations:
                    if data_dictionary[error_rate_key]["NNNNN"][op[0]][op[1]]==-1:
                        data_dictionary[error_rate_key]["NNNNN"][op[0]][op[1]]=1
                    else:
                        data_dictionary[error_rate_key]["NNNNN"][op[0]][op[1]]+=1
                data_dictionary[error_rate_key]["NNNNN"]["total"]+=1
        else:
            #assume we have a NNN-3 base barcode strand
            #use part of promoter region to be reference closest to the 3-base barcode, common across all the files studied
            reference="GCGCGC"
            #print "NNN 3 barcode"
            reference_start=_strand.find(reference)
            if reference_start==-1: continue
            #reference occurs after the barcode
            NNN_3_payload=payload[file_strand_archs[file_base_name]["payload_2"]]
            lv_NNN_3_payload=lv.distance(NNN_3_payload,_strand[end_point:end_point+len(NNN_3_payload)])
            pre_barcode_start=reference_start-4
            payload_ld=lv.distance(payload[file_strand_archs[file_base_name]["payload_1"]],payload[file_strand_archs[file_base_name]["payload_2"]])
            if _strand[pre_barcode_start]=='A' and  lv_NNN_3_payload<payload_ld/2: 
                #high confidence we found a barcode
                barcode_start=pre_barcode_start+1
                barcode=_strand[barcode_start:reference_start]
                assert len(barcode)==3
                data_dictionary[count_key][barcode_to_index(barcode)]+=1
            if lv_NNN_3_payload<payload_ld/2:
                edit_operations=lv.editops(NNN_3_payload,_strand[end_point:end_point+len(NNN_3_payload)])
                for op in edit_operations:
                    if data_dictionary[error_rate_key]["NNN"][op[0]][op[1]]==-1:
                        data_dictionary[error_rate_key]["NNN"][op[0]][op[1]]=1
                    else:
                        data_dictionary[error_rate_key]["NNN"][op[0]][op[1]]+=1
                data_dictionary[error_rate_key]["NNN"]["total"]+=1
if __name__=="__main__":
    import argparse
    import os
    import pandas as pd #going to dump things out using data frames
    import pickle as pi
    import numpy as np
    import Levenshtein as lv
    barcode_array=init_barcode_to_index_dict() #this array will be used as row labels in the data frame 
    
    parser = argparse.ArgumentParser(description="Analyze strands for Doris")
    parser.add_argument('--range', dest="fq_range", action="store", default="1-10", help="Range for fastq files")
    parser.add_argument('--fastq_directory',dest="fq_dir", action="store", default=None, help="Directory for stripped fastq files")

    args = parser.parse_args()

    lower,upper=args.fq_range.split("-")
    lower_int = int(lower)
    upper_int = int(upper)

    sample_files=[]

    dir_sorted=os.listdir(args.fq_dir)
    dir_sorted.sort()

    data_dictionary={}
    
    for _file in dir_sorted[lower_int:upper_int+1]: 
        if os.path.isfile(args.fq_dir+'/'+_file):
            sample_files.append(_file)

    if not os.path.exists("DORIS_DATA"):
        os.mkdir("DORIS_DATA")

    for fq_file in sample_files:
        print (fq_file)
        file_path=args.fq_dir+'/'+fq_file
        file_base_name=get_experiment_base_name(_file)
        sequence_strands=[line.rstrip('\n') for line in open(file_path)]
        #generate result dictionary keys 
        count_key=fq_file+"_"+"count"
        total_reads_key=fq_file+"_"+"total_reads"
        error_rate_key=fq_file+"_"+"error_rate"
        if count_key not in data_dictionary:
            #initialize an array of counters for counting each barcode
            data_dictionary[count_key]=[0]*1088
        if total_reads_key not in data_dictionary:
            data_dictionary[total_reads_key]=['--']*1088
            data_dictionary[total_reads_key][0]=len(sequence_strands)
        if error_rate_key not in data_dictionary:
            data_dictionary[error_rate_key]={}
            data_dictionary[error_rate_key]["NNN"]={}
            data_dictionary[error_rate_key]["NNNNN"]={}
            data_dictionary[error_rate_key]["NNN"]["replace"]=[-1]*200
            data_dictionary[error_rate_key]["NNN"]["delete"]=[-1]*200
            data_dictionary[error_rate_key]["NNN"]["insert"]=[-1]*200
            data_dictionary[error_rate_key]["NNN"]["total"]=0
            data_dictionary[error_rate_key]["NNNNN"]["replace"]=[-1]*200
            data_dictionary[error_rate_key]["NNNNN"]["delete"]=[-1]*200
            data_dictionary[error_rate_key]["NNNNN"]["insert"]=[-1]*200
            data_dictionary[error_rate_key]["NNNNN"]["total"]=0
        analyze_strands(file_base_name,fq_file,sequence_strands,data_dictionary)
        
    #make a data frame using the collected data
    count_dump_path="DORIS_DATA/"+get_experiment_base_name(sample_files[0])+'_count.csv'
    error_rate_dump_path="DORIS_DATA/"+get_experiment_base_name(sample_files[0])+'_error_rate.csv'
    count_file=open(count_dump_path,'w+')
    error_rate_file=open(error_rate_dump_path,'w+')
    
    count_dict={}
    error_rate_dict={}
    for key in sorted(data_dictionary.keys()):
        if "count" in key or "total_reads" in key:
            count_dict[key]=data_dictionary[key]
        elif "error_rate" in key:
            for strand_type in data_dictionary[key]:
                for error_type in data_dictionary[key][strand_type]:
                    if "total" in error_type: continue
                    error_rate_dict_key=key.split('error_rate')[0]+error_type+"_"+strand_type
                    error_rate_dict[error_rate_dict_key]=data_dictionary[key][strand_type][error_type]
                    for count_index, count in enumerate(error_rate_dict[error_rate_dict_key]):
                        if count==-1:
                            error_rate_dict[error_rate_dict_key][count_index]=0
                        else:
                            error_rate_dict[error_rate_dict_key][count_index]=float(count)/float(data_dictionary[key][strand_type]["total"])

    error_rate_frame=pd.DataFrame(error_rate_dict)
    count_frame=pd.DataFrame(count_dict,index=barcode_array)
    count_frame.to_csv(count_dump_path)
    error_rate_frame.to_csv(error_rate_dump_path)
    #seriialize the results dictionary
    #data_dict_pickle_path="DORIS_DATA/"+get_experiment_base_name(sample_files[0])+'_pickled.dict'
    
            

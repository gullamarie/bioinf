#!/usr/bin/python

"""
Author: M. Gulla

This script reads merged SNV and INDEL VCF files genereated by SNV caller 
STRELKA2 that are subsequently merged with mergeVCF Picard tool 

Optional argument that must be set: 
    - both_criteria (bool) : Filter on both PASS quality marking or 
        just by n.o. read counts 
        False: Only filter by read count values in relation to threshold 
        True: Filter both by read count threshold and by quality=PASS
    - write_summary (bool) : Wheter or not to write summary file with 
        filtering statistics for each processed sample 
"""
# Edit these config variables to preferred settings !!!  
both_criteria = 2 #1
write_summary = True

#import csv
#import shlex
import re
#import math
import datetime
#import sys
import operator
import collections
import numpy as np
#from collections import defaultdict
#from operator import itemgetter
#import os.path

samples=["1-T-1-N", "2-T-1-N", "3-T-3-N", "4-T-4-N", "5-T-5-N", "6-T-6-N"]


# "downstream"
# "exonic"
# "exonic;splicing"
# "intergenic"
# "intronic"
# "ncRNA_exonic"
# "ncRNA_intronic"
# "ncRNA_splicing"
# "ncRNA_UTR3"
# "ncRNA_UTR5"
# "splicing"
# "upstream"
# "upstream;downstream"
# "UTR3"
# "UTR5"
# "UTR5;UTR3"


filter_percentage_alt_reads = 5.0
filter_t_alt_reads = 5
filter_n_alt_reads = 1 # 1: Accept 0 or 1 alt reads in normal
filter_t_reads = 15
filter_n_reads = 10

filter_func = ["exonic", "exonic;splicing", "splicing"]#, "intronic", "ncRNA_exonic"]
# , "UTR3", "UTR5", "exonic;splicing", "splicing"]
# , "downstream", "intergenic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
#"ncRNA_splicing", "ncRNA_UTR3", "ncRNA_UTR5", "splicing", "upstream", 
#"upstream;downstream", "UTR5;UTR3"]

snv_algorithm = ["strelka"]
indel_algorithm = ["strelka"]
# filter_mutect_judgement = ["KEEP", "REJECT"]
# filter_strelka_judgement = ["PASS"]

def passes_read_filter_with_msg(t_ref_count, t_alt_count, n_ref_count, n_alt_count): 
    #set filtering thresholds
    filter_percentage_alt_reads = 5.0
    filter_t_alt_reads = 5
    filter_n_alt_reads = 1 # 1: Accept 0 or 1 alt reads in normal
    filter_t_reads = 15 #15
    filter_n_reads = 10 #10
    
    ret = False
    
    if not (t_alt_count >= filter_t_alt_reads):
        #print("t_alt_count < filter_t_alt_reads")
        #few_t_alt_c += 1
        pass
         
    else: 
        #t_alt_count >= filter_t_alt_reads passed
        if not (n_alt_count <= filter_n_alt_reads): 
            #print("n_alt_count > filter_n_alt_reads")
            pass 
            #many_n_alt_c += 1
        else: 
            #n_alt_count <= filter_n_alt_reads passed
            if not ((t_ref_count + t_alt_count) >= filter_t_reads):
                #print("(t_ref_count + t_alt_count) < filter_t_reads") 
                #few_t_total_c += 1
                pass
            else: 
                 #(t_ref_count + t_alt_count) >= filter_t_reads
                if not ((n_ref_count + n_alt_count) >= filter_n_reads):
                    #print("(n_ref_count + n_alt_count) < filter_n_reads")
                    #few_n_total_c += 1 
                    pass
                else:
                     if not (100.0*t_alt_count/(t_alt_count + t_ref_count) >= filter_percentage_alt_reads ):
                          #print("100.0*t_alt_count/(t_alt_count + t_ref_count) < filter_percentage_alt_reads")
                          #low_perc += 1
                          pass
                     else: 
                          ret = True

    return ret               

def getOldKey(dictionary, new_key):
     pos = new_key[2]
     new_ref = new_key[4]
     new_alt = new_key[5]
     keys = list(dictionary.keys())
     positions = [pos - 1, pos, pos + 1, pos + 2]
     #find key in dict with matching position
     
     return_key = False
     #loop through all possible positions 
     for pp in positions:
          old_key = (new_key[0:2]) + (pp,) + (new_key[3],)
          
          for key in keys: 
               if not len(np.intersect1d(old_key, key)) == 4: 
                    continue
               else: 
                    #old_key in key: use key 
                    #print(key)
                    old_ref = key[4]
                    old_alt = key[5]
                    if new_ref == '-': 
                         if new_alt in old_alt:
                                             
                              #found matching key, return this oldkey
                              return_key = key
                              #continue
                                   
                    elif new_alt == '-': 
                         if new_ref in old_ref: 
                              return_key = key
                              #continue
                    

     if return_key: 
          return return_key
     else: 
          return False
          
          
bases = {
    'A' : 1,
    'C' : 2, 
    'G' : 3, 
    'T' : 4
    }

if write_summary: 
    #write qc file for all samples 
    summary_file = 'annotation_file_combo_raw.txt'
    summary_w = open(summary_file, 'w')
    header_line = 'sample_name\tno_variants_prefilter\tno_variants_postfilter\tno_snvs\tno_indels\tno_poor-qual_filter-passed\tpercentage_lowEVS\tpercentage_lowDepth\n'
    summary_w.write(header_line)
    summary_w.close()



all_snvs = {}
snvs_dict = {}
indels_dict = {}
all_muts = {}

total_unfiltered_snvs = 0
total_filtered_snvs = 0


for sample in samples:
     print(sample)
     tumor_id = sample.split("_")[0]

     snvs = {} 
     
     #clarification: name of variable snv dict contains both indels and snvs
     indel_count = 0
     snv_count = 0
     
     qual_lowevs = 0
     qual_lowdepth = 0
     passed_poor_qual = 0

     # STRELKA combo vcf file:
     #snv_file = "/merged_output/%s-somatic.combined.vcf" % (tumor_id)
     snv_file = "%s-somatic.combined.vcf" % (tumor_id)
     with open(snv_file) as f:
        
          for line in f:
               #do something with data
               cnt = 0
               is_snv = False
               row = line.split('\t')
        	
               cnt = cnt + 1
            
               if row[0] == "#" or row[0] == "contig":
                    header = row
               elif row[0] == '#CHROM':
                    n_idx = row.index('NORMAL')
                    t_idx = row.index('TUMOR\n')
               
               elif (row[0] != "chrM" and row[0][:3] == "chr" and len(row[0]) <= 5) :
                    #print len(row)
                    total_unfiltered_snvs += 1
                    
                    chr = row[0]
                    pos = int(row[1])
                    ref = row[3]
                    alt = row[4]
                    qual = row[6]
                    
                    if (len(ref) == 1 and len(alt) == 1) :
                         #snv
                         is_snv = True
                         var_type = 'snv'
                         snv_count += 1
                         
                         n_ref_count = int(row[n_idx].split(':')[bases[ref]+4].split(',')[0])
                         n_alt_count = int(row[n_idx].split(':')[bases[alt]+4].split(',')[0])
                         t_ref_count = int(row[t_idx].split(':')[bases[ref]+4].split(',')[0])
                         t_alt_count= int(row[t_idx].split(':')[bases[alt]+4].split(',')[0])
                         
                         
                    else: 
                         #indel 
                         var_type = 'indel'
                         is_snv = False
                         indel_count += 1 
                         
                         n_ref_count = int(row[n_idx].split(':')[3].split(',')[0])
                         n_alt_count = int(row[n_idx].split(':')[4].split(',')[0])
                         t_ref_count = int(row[t_idx].split(':')[3].split(',')[0])
                         t_alt_count= int(row[t_idx].split(':')[4].split(',')[0])
                    
                    #find zygosity
                    if row[t_idx].split(':')[0] == '0/1':
                         zyg = 'het'
                    if row[t_idx].split(':')[0] == '0/0':
                         zyg = 'hom'
                         
                    
              
                    if both_criteria: 
                         #run PASS and read count filtering 
                         if (qual == 'PASS' and passes_read_filter_with_msg(t_ref_count, t_alt_count, n_ref_count, n_alt_count) ):
                              snvs[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              total_filtered_snvs += 1
                              
                              #sort into snvs or indels dict 
                              if is_snv: 
                                  snvs_dict[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              else: 
                                  indels_dict[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              all_muts[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                         else:
                              pass
                    else: 
                         #only run read count filering 
                         if passes_read_filter_with_msg(t_ref_count, t_alt_count, n_ref_count, n_alt_count):
                              snvs[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              total_filtered_snvs += 1
                              
                              #sort into snvs or indels dict 
                              if is_snv: 
                                  snvs_dict[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              else: 
                                  indels_dict[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                              all_muts[(sample, chr, pos, var_type, ref, alt)] = [t_ref_count, t_alt_count, n_ref_count, n_alt_count]
                         else: 
                              pass
                    
                    #check ref alt ratio in ref to zyg
                    ratio = (t_ref_count/(t_alt_count+0.1))
                    if ratio < 0.75 and zyg == 'hom':
                         print(t_ref_count, t_alt_count, ratio)
                              
                    if (not qual == 'PASS' and passes_read_filter_with_msg(t_ref_count, t_alt_count, n_ref_count, n_alt_count) ):
                         #print("Passes read filter but poor quality: ", chr, pos)
                         passed_poor_qual += 1 
                    
                    if not (qual == 'PASS'):
                         if ',' in qual: # > one quality mark 
                              qual = qual.split(',')
                              if 'LowEVS' in qual[0] or 'LowEVS' in qual[1]: 
                                   qual_lowevs += 1
                              if 'LowDepth' in qual[0] or 'LowDepth' in qual[1]:
                                   qual_lowdepth += 1
                         else: # one quality mark 
                              if qual == 'LowEVS': 
                                   qual_lowevs += 1
                              elif qual == 'LowDepth': 
                                   qual_lowdepth += 1
                                 
     
     for k, v in snvs.items():
          if len(v) < 4: 
               print(k, v)
               break


     print("Total snvs: ", total_unfiltered_snvs)
     print("Total snvs after filtering: ", total_filtered_snvs)

	
     #Check that all v values are in list 
     for k, v in snvs.items():
          if not type(v) == list: 
               print("type (v) not list")
               print(k, v)
     
     print("Reading anno", sample)

     anno_file = "%s-merged-strelka2.vcf.annovar.summary.hg38_multianno.csv" % (tumor_id)


     with open(anno_file) as f:
        
          for line in f:    
               
               cnt = 0
               row = line.split(",")
    
               cnt = cnt + 1
    
               if row[0] == "Func":
                    #header = row
                    pass
               #check this one 
               elif row[0] != "chrM" and row[0][:3] == "chr" and len(row[0]) <= 5 :
                    chr = row[0]
                    pos = int(row[1])
                    ref = row[3]
                    alt = row[4]
                    

                    if (len(ref) > 1 or len(alt) > 1):
                         if '-' in ref or '-' in alt: 
                              var_type = 'indel'
                    elif (len(ref) == 1 and len(alt) == 1):
                         if '-' in ref or '-' in alt:
                              var_type = 'indel'
                         else: 
                              var_type = 'snv'
                    
                    #find func
                    func = row[5].strip('"')
                    if func in filter_func and len(func) < 5:
                         print(chr, pos, ref, alt, func)

                    gene = row[6].strip('"')
                    syn = row[8].strip('"')
                    AAChange = row[9].strip('"')
                    if '"' in AAChange:
                         AAChange = AAChange.strip('"')
                    
                    dbsnp144 = row[17]
                    if '"' in dbsnp144:
                         dbsnp144 = dbsnp144.strip('"')
                    
                    
                    new_key = (sample, chr, pos, var_type, ref, alt)
                    
                    
                    if var_type == 'snv':
                         orig_key = (sample, chr, pos, var_type, ref, alt)
                         if orig_key in snvs:
                              insert_list = [dbsnp144, gene, func, syn, AAChange]
                                                         

                              values = snvs[orig_key] + insert_list
                              snvs[orig_key] = values 
                              

                              
                    elif var_type == 'indel':
                        
                         
                         old_key = getOldKey(snvs, new_key)
                         
                         if old_key: #found key in snvs: 
                              try: 
                                   #replace old key with new
                                   values = list(snvs[old_key])
                                   snvs[new_key] = snvs.pop(old_key) + [func, gene, syn, AAChange, dbsnp144]
                                              
                              except KeyError: 
                                   print('insuccsessful, key error in snv dict', new_key)
                                   pass

     c = 0              
     for k, v in snvs.items():
          if not len(v) > 4: 
               print(k, v)
               c += 1
     print(c, " instances of variations without annotation info" )
               
               

     print("Done reading anno.", sample)
     print("Length:", len(snvs), sample)

     all = len(snvs)
     cnt = 0
    
     #Double check filtering criteria in snv-dict 

     for k, v in snvs.items():
          #print(type(v[0]))

          reads = v[0:4]
          try:
               fun = v[4]
          
          except KeyError: 
               print('no annotation info ', k)
               continue
          

          if fun not in filter_func: 
               if int(reads[1]) < filter_t_alt_reads: 
                    if int(reads[0]) + int(reads[1]) < filter_t_reads:
                         if int(reads[2]) + int(reads[3]) < filter_n_reads: 
                              if 100.0*float(reads[1])/(reads[0]+reads[1]) < filter_percentage_alt_reads or int(reads[3]) > filter_n_alt_reads:
                                   print("Deleting from dict: ")   #(k, v)
                                   #print(k, v)
                                   del snvs[k]
    
     print("Done deletion.")

     filtered = len(snvs)
     print(sample + ": Preprocessed " + str(total_unfiltered_snvs) + " SNV variants. Kept " + str(filtered) + ".")
     
     all_snvs.update(snvs)
     total_filtered_snvs += filtered
     

     #write to anno file here 
     if write_summary:
          #sep = '\t'
          outline = [str(i) for i in [sample, total_unfiltered_snvs, total_filtered_snvs, snv_count, indel_count, passed_poor_qual, qual_lowevs, qual_lowdepth]]
          outline = '\t'.join(outline) + '\n'
          #read in lines in file, append new
          summary_o = open(summary_file, 'r')
          summary_lines = summary_o.readlines()
          summary_lines.append(outline)
          summary_o.close()
          #write out lines
          summary_w = open(summary_file, 'w')
          for line in summary_lines:
               summary_w.write(line)
          summary_w.close()

     
#end samples loop


#Sort all found variants 
all_muts = dict(all_snvs.items())
all_muts = sorted(all_snvs.items(), key=operator.itemgetter(0))

all_muts = collections.OrderedDict(sorted(all_snvs.items()))

d = datetime.datetime.now().strftime("%y-%m-%d-%H-%M")

out_file = "MUTATIONS-SNVS-strelka-%s.csv" % (d)
#out_file = "MUTATIONS-STRELKA2-COMBINED-anno-format-%s.csv" % (d)

mywriter = open(out_file, 'w')
samples = [str(i) for i in samples]
filter_func = [str(i) for i in filter_func]
smpls = ",".join(samples)
fncs = ",".join(filter_func)

mywriter.write("# Samples included: %s" % (smpls) + '\n')
mywriter.write("# Funcs: %s" % (fncs) + '\n')
mywriter.write("# filter_n_reads: %s" % (str(filter_n_reads)) + '\n')
mywriter.write("# filter_t_reads: %s" % (str(filter_t_reads)) + '\n')
mywriter.write("# filter_percentage_alt_reads: %s" % (str(filter_percentage_alt_reads)) + '\n')
header_line = ','.join(["sample","chromosome","position","type","ref","alt","novel/dbsnp","covered","t_ref_count","t_alt_count","n_ref_count","n_alt_count","dbsnp138NonFlagged","gene","func","syn","Otherinfo","Obs","AAChange","judgement","filtered"])
mywriter.write(header_line)

#write out all info from dict

for k, v in all_muts.items():
     outlist = list(k)
     v = list(v)
     for i in range(0,len(v)): 
          el = v[i]
          if (i == 0):
               l = str(el)
               l = re.sub('[()]', '', l) 
               l = l.split(',')
               l = [j.strip() for j in l]
               outlist.extend(l)

          elif (i == 1) and ('"' in str(el)):    
               outlist.append(el.strip('"')) 
          else: 
               outlist.append(el)
     
     outlist = [str(i) for i in outlist]     
     outline = ','.join(outlist) + '\n'
     mywriter.write(outline)     
     
mywriter.close()

print("Created file: ", out_file)




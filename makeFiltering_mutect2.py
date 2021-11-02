#!/usr/bin/python

import csv
import shlex
import re
import math
import datetime
import sys
from operator import itemgetter
import os.path


samples_t = ["T-7", "T-8"]
samples_n = ["N-3", "N-4"]

filter_percentage_alt_reads = 5.0
filter_percentage_alt_reads_normal = 1.0
filter_t_alt_reads = 5
filter_n_alt_reads = 1 # 1: Accept 0 or 1 alt reads in normal
filter_t_reads = 15
filter_n_reads = 15
filter_func = ["exonic", "exonic;splicing", "splicing"]


all_snvs = {}
all_mut = {}
idx = 0

for sample in samples_t:
	sample_n = samples_n[idx]
	idx += 1

	print sample, "running."
	snvs = {}

	# MUTECT2:
	snv_file = "/MUTECT-2/%s_%s_MUTECT2_filtered.vcf" % (sample, sample_n)
	snvso = open(snv_file, "r")
	snvs_lines = snvso.readlines()

	i = 0
	for row in csv.reader(snvs_lines, delimiter="\t"):
		i += 1
		if i % 100000 == 0:
			print i
		if row[0][0] == "#":
			header = row
		elif row[0] != "chrM" and row[0][:3] == "chr" and len(row[0]) <= 5:
			if row[0] == "chrX":
				chr = 23
			elif row[0] == "chrY":
				chr = 24
			else:
			chr = int(row[0][3:])
			pos = int(row[1])
			ref = row[3]
			alt = row[4]
			if ',' in alt:
				alts = alt.split(',')

			else 
				if len(ref) > len(alt): # deletion
					if "," in alt:
						alt1 = "-"
						reflist = ref.split(",")
						ref1 = reflist[1:]
					else:
						ref = ref[1:]
						alt = "-"
						pos += 1
					if "," in alt:
						alt1 = 
				elif len(ref) < len(alt): # insertion
					ref = "-"
					alt = alt[1:]
				elif "," in alt:
					alt 
					alt1 = alt[0]
					ref1 = ref[0]
				covered = "."
			t_ref_count = int(row[9].split(":")[1].split(",")[0])
			t_alt_count = int(row[9].split(":")[1].split(",")[1])
			n_ref_count = int(row[10].split(":")[1].split(",")[0])
			n_alt_count = int(row[10].split(":")[1].split(",")[1])
			filtered = "."
			keep = "."

			snvs[(sample, chr, pos, ref, alt)] = (".", covered, t_ref_count, t_alt_count, n_ref_count, n_alt_count, keep, filtered)

	snvso.close

	print "Done Mutect"
	
	#  MUTECT annotering:

	anno_file = "/ANNOVAR/MUTECT2/%s/%s-mutect2.vcf.annovar.summary.hg38_multianno.csv" % (sample, sample)

	annos = open(anno_file, "r")
	anno_lines = annos.readlines()

	i = 0
	for row in csv.reader(anno_lines, delimiter=","):
		i += 1
		if i % 100000 == 0:
			print i
		if row[0] == "Chr":
               	        header = row
		elif row[0] != "chrM" and row[0][:3] == "chr" and len(row[0]) <= 5: # Avoid all extra chromosomes
			if row[0] == "chrX":
				chr = 23
			elif row[0] == "chrY":
				chr = 24
			else:
				chr = int(row[0][3:])
			pos = int(row[1])
			ref = row[3]
			alt = row[4]
			func = row[5]
			gene = row[6].replace(",",";")
			syn = row[8]
			AAChange = row[9]
			dbsnp144 = row[17]

			if (sample, chr, pos, ref, alt) in snvs:
				snvs[(sample, chr, pos, ref, alt)] = (snvs[(sample, chr, pos, ref, alt)], dbsnp144, gene, func, syn, AAChange)
			else:
				print "Not found: ", sample, chr, pos, ref, alt, gene, "\n"
	annos.close

	print "Done Anno Mutect"

        all = len(snvs)
	print "All SNVs:", all
	i = 0
       	for k, v in snvs.items():
               	if v[3] not in filter_func or int(v[0][3]) < int(filter_t_alt_reads) or int(v[0][2]) + int(v[0][3]) < int(filter_t_reads) or int(v[0][4]) + int(v[0][5]) < int(filter_n_reads) or 100.0*float(v[0][3])/(v[0][2]+v[0][3]) < filter_percentage_alt_reads or int(v[0][5]) > filter_n_alt_reads or 100.0*float(v[0][5])/(v[0][4] + v[0][5]) > filter_percentage_alt_reads_normal:
                       	del snvs[k]

	print "Done SNV deletion"

        filtered = len(snvs)
       	print sample + ": Preprocessed " + str(all) + " candidate SNVs. Kept " + str(filtered) + "."
       	all_snvs = dict(all_snvs.items() + snvs.items())


all_muts = dict(all_snvs.items())

all_muts_sorted = sorted(all_muts.keys(), key=itemgetter(0, 1, 2))

d = datetime.datetime.now().strftime("%y-%m-%d-%H-%M")

out_file = open("MUTATIONS-%s.csv" % (d), "wb")

mywriter = csv.writer(out_file)

smpls = "; ".join(samples_t)
fncs = "; ".join(filter_func)

mywriter.writerow(["# Samples included: %s" % (smpls)])
mywriter.writerow(["# Funcs: %s" % (fncs)])
mywriter.writerow(["# filter_t_reads: %s" % (str(filter_t_reads))])
mywriter.writerow(["# filter_n_reads: %s" % (str(filter_n_reads))])
mywriter.writerow(["# filter_t_alt_reads: %s" % (str(filter_t_alt_reads))])
mywriter.writerow(["# filter_n_alt_reads: %s" % (str(filter_n_alt_reads))])
mywriter.writerow(["# filter_percentage_alt_reads: %s" % (str(filter_percentage_alt_reads))])
mywriter.writerow(["# filter_percentage_alt_reads_normal: %s" % (str(filter_percentage_alt_reads_normal))])


mywriter.writerow(["sample","chromosome","position","ref","alt","novel/dbsnp","covered","t_ref_count","t_alt_count","n_ref_count","n_alt_count","dbsnp144","gene","function","syn/nonsyn","Detailed_info","judgement","failure_reasons"])

for k in all_muts_sorted:
	if int(all_muts[k][0][2]) + int(all_muts[k][0][3]) >= int(filter_t_reads) and int(all_muts[k][0][4]) + int(all_muts[k][0][5]) >= int(filter_n_reads) and 100.0*float(all_muts[k][0][3])/(float(all_muts[k][0][2]) + float(all_muts[k][0][3])) >= filter_percentage_alt_reads and int(all_muts[k][0][5]) <= filter_n_alt_reads and all_muts[k][3] in filter_func:
		if k[1] == 23:
			mywriter.writerow([k[0], "chrX", k[2], k[3], k[4], all_muts[k][0][0], all_muts[k][0][1], all_muts[k][0][2], all_muts[k][0][3], all_muts[k][0][4], all_muts[k][0][5], all_muts[k][1], all_muts[k][2], all_muts[k][3], all_muts[k][4], all_muts[k][5], all_muts[k][0][6], all_muts[k][0][7]])
                elif k[1] == 24:
                        mywriter.writerow([k[0], "chrY", k[2], k[3], k[4], all_muts[k][0][0], all_muts[k][0][1], all_muts[k][0][2], all_muts[k][0][3], all_muts[k][0][4], all_muts[k][0][5], all_muts[k][1], all_muts[k][2], all_muts[k][3], all_muts[k][4], all_muts[k][5], all_muts[k][0][6], all_muts[k][0][7]])
                else:
			print k
			print all_muts[k]
                        mywriter.writerow([k[0], "chr" + str(k[1]), k[2], k[3], k[4], all_muts[k][0][0], all_muts[k][0][1], all_muts[k][0][2], all_muts[k][0][3], all_muts[k][0][4], all_muts[k][0][5], all_muts[k][1], all_muts[k][2], all_muts[k][3], all_muts[k][4], all_muts[k][5], all_muts[k][0][6], all_muts[k][0][7]])


out_file.close

print "Created file: MUTATIONS-%s.csv" % (d)

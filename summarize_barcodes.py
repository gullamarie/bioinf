"""
Reads in all bam-files in subfolders and makes a 
table with names of bam-files with matching barcode. 


Matches tumor and normal TCGA files

"""

import os


# Read in all bam files
bam_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(".") for f in filenames if os.path.splitext(f)[1] == '.bam']

# Read in match-file
matchfile = open('match.txt', 'r')
# store info in dict
match = {}
for line in matchfile.readlines():
        line = line.split()
        if len(line) > 1:
                match[line[0]] = line[1]

#print(match)




def find_bam_pairs(bam_files):
        #
        barcodes=[]
        ids={}
        complete={}
        pairs=[]
        out_tumor = open('missing_normal.txt', 'w')
        out_normal = open('missing_tumor.txt', 'w')

        for entry in bam_files:
                entry_split = entry.split("/")
                #print(entry_split)

                filepath = entry_split[2]

                #if not barcode, get from matchfile
                if not ('TCGA-') in filepath:
                        barcode_wh = match[filepath]
                        barcode = barcode[0:12]
                        #print(barcode_wh, barcode)
                else:
                        barcode = entry_split[2][:12]
                        barcode_wh = entry_split[2][:16]
                #print(barcode_wh, filepath)
                #print('\n')

                ids[barcode_wh] = [barcode_wh, filepath]

                #check if already in barcodes
                if barcode in barcodes:
                        pair_barcode_wh = complete[barcode]
                        if (pair_barcode_wh != barcode_wh) & (pair_barcode_wh[0:12] == barcode_wh):
                                pairs.append((pair_barcode_wh, barcode_wh))
                else:
                        barcodes.append(barcode)
                        complete[barcode]=barcode_wh

        print(len(barcodes))
        print(len(ids.keys()))

        #find missing ones
        for entry in ids.keys():
                #print(entry)
                barcode_wh = entry #str(complete[entry])
                #print(barcode_wh, type(barcode_wh))
                if not barcode_wh in str(pairs):

                        #test for tumor
                        if '01A' in entry:
                                outline = barcode_wh + '\t' + str(ids[barcode_wh]) + '\n'
                                out_tumor.write(outline)
                        #test for normal:
                        if ('11A' in entry) | ( '10A' in entry):
                                outline = barcode_wh + '\t' + str(ids[barcode_wh]) + '\n'
                                out_normal.write(outline)


        out_tumor.close()
        out_normal.close()
        return pairs, ids


def make_bam_file(pairs, ids):
        outfile = open('decodefile.txt', 'w')

        for pair in pairs:

                #classify T vs N
                if '01A' in pair[0]:
                        tumor_barcode = pair[0]
                        normal_barcode = pair[1]
                elif ('11A' in pair[0]) | ( '10A' in pair[0]):
                        tumor_barcode = pair[1]
                        normal_barcode = pair[0]

                #get info
                tumor_list = ids[tumor_barcode]
                tumor_info = '\t'.join(tumor_list)
                normal_list = ids[normal_barcode]
                normal_info = '\t'.join(normal_list)
                outline = tumor_info + '\t' + normal_info + '\t' + '\n'
                outfile.write(outline)

        outfile.close()



def main():

        bam_files=[os.path.join(dp, f) for dp, dn, filenames in os.walk(".") for f in filenames if os.path.splitext(f)[1] == '.bam']



        bam_pairs, bam_ids = find_bam_pairs(bam_files)

        for elem in bam_pairs:
                print(elem)

        print("******")

        make_bam_file(bam_pairs, bam_ids)



if __name__ == '__main__':
        main()





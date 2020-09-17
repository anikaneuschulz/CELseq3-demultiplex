import pysam
import sys
import os
import pdb
import subprocess

# to prepare index_file:
# unzip fastq that has barcode read
# cat {input} | paste -d '\t' - - - - | sort -k1,1 > {output}
# cut -f1 {output} > {output_readnames}
#  cut -f2 {output} | cut -c7-14 > {output_barcodes}
# paste {output_readnames} {output_barcodes} > {output_readnames_barcodes}


reads = pysam.AlignmentFile(sys.argv[1], "rb")
#this has readnames and barcodes
index_file = sys.argv[2]
#this has all barcodes that were used in the sample
whitelist = sys.argv[3]

corrected_keys = {}
barcode_dict = {}
final_readname_dict = {}

#group readnames by barcodes
print("gathering raw barcodes...")
index = open(index_file)
for line in index:
    l = line.split("\t")
    readname = l[0]
    barcode = l[1]
    if not barcode in barcode_dict:
        barcode_dict[barcode] = []
        barcode_dict[barcode].append(readname.split(" ")[0])
    else:
        barcode_dict[barcode].append(readname.split(" ")[0])
index.close()


#max_BC = ""
#maxentry = 0
#for entry in barcode_dict:
#    leng = len(barcode_dict[entry])
#    if leng > maxentry:
#        maxentry = leng
#        max_BC = entry

#print(max_BC)
#print(maxentry)

print("loading whitelist...")
#load whitelist of Barcodes
wlist = {"####"}
with open(whitelist, "r") as w:
    for line in w:
        wlist.add(line)
wlist.remove("####")

print("correcting barcodes...")
# correct barcodes for one allowed mismatch
for key in barcode_dict:
    for wkey in wlist:
        diff = 0
        for ch1, ch2 in zip(key, wkey):
            if ch1 != ch2:
                diff+= 1
                if diff == 2:
                    break
                        #save some time here
        if diff <= 1:
            if not wkey in corrected_keys:
                corrected_keys[wkey] = []
            corrected_keys[wkey].append(key)
            break

print("gathering readnames...")
# prepare final lists of readnames using corrected barcodes
for key in corrected_keys:
    for entry in corrected_keys[key]:
        if not key in final_readname_dict:
            final_readname_dict[key] = []
        final_readname_dict[key]+=barcode_dict[entry]

# convert lists to sets so membership tests don't take literally forever
for key in final_readname_dict:
    final_readname_dict[key] = set(final_readname_dict[key])




print("preparing reads for demultiplexing...")
# prepare dictionary if ready attibuted to barcodes
final_reads_dict = {}
reads_processed = 0
sys.stdout.write('0 reads processed')
for read in reads.fetch():
    reads_processed+= 1
    sys.stdout.write('\r' + str(reads_processed) + ' reads processed.')
    readname = "@" + read.query_name
    for key in final_readname_dict:
        if readname in final_readname_dict[key]:
            if not key in final_reads_dict:
                final_reads_dict[key] = []
            final_reads_dict[key].append(read)
            break


# writing reads...
for barcode in final_reads_dict:
    print("\nwriting barcode " + barcode.rstrip())
    ofn = sys.argv[1][:-4] + "_" + barcode.rstrip() + ".sam"
    outfile = pysam.AlignmentFile(ofn, "w", template=reads)
    for read in final_reads_dict[barcode]:
        outfile.write(read)



#pdb.set_trace()


outfile.close()



#ofn = sys.argv[1][:-4] + "_cluster_" + str(key) + ".sam"
#outfile = pysam.AlignmentFile(ofn, "w", template=bamfile)

#pdb.set_trace()

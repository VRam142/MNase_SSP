import os,sys,re
import pysam
import random
import pickle
from collections import Counter
import gzip
from bx.intervals.intersection import Intersecter, Interval
import random
import numpy as np
def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def readIterator(filenames,chrom):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom):
        yield read
      input_file.close()

def calculatePerBase(filenames, length_hist, chrom):
  for read in readIterator(filenames,chrom):	
    if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
    if isSoftClipped(read.cigar): continue
    if read.is_paired:
      if read.mapq < 5: continue
      if read.mate_is_unmapped: continue
      if read.rnext != read.tid: continue
      if read.is_read1:
        if read.isize == 0: continue
        rstart = min(read.pos,read.pnext)+1 # 1-based
        rend = rstart+abs(read.isize)-1 # end included
        rlength = rend-rstart+1
	if rlength > 300: continue 
	length_hist[rlength] += 1


def main():
  genome_file = open(sys.argv[1])
  filename = sys.argv[2]
  sim_type = sys.argv[3]
  filenames = [filename]
  chroms = []
  length_hist = Counter()
  for line in genome_file:
    split = line.split()
    if sim_type == "True":
      chrom = split[0].split("chr")[1]
      try: int(chrom) + 1
      except:
        if chrom == "X" or chrom == "Y": chrom = chrom
        else: continue
    else:
      chrom = split[0]
    calculatePerBase(filenames, length_hist, chrom)
  for i in range(20,300):
        print "%s\t%s\t%s\t%s" % (i, float(length_hist[i]) / np.sum(length_hist.values()), length_hist[i], np.sum(length_hist.values()))
  genome_file.close()

if __name__ == "__main__":
  main()

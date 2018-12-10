import os,sys,re
import pysam
import random
import numpy as np
from collections import Counter
import gzip
import random
from intervaltree import Interval, IntervalTree

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

def readIterator(filenames,chrom, start, end):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom, start, end):
        yield read
      input_file.close()

def calculatePerBase(filenames, tss, sim_type, valid_chroms, minimum, maximum):
  all_pos = {}
  tree = IntervalTree()
  vals = []
  intervals = []
  for i in range(minimum, maximum, 10):
	vals.append(i)
  vals.append(i + 10)
  for i in range(len(vals) - 1):
        tree[vals[i]:vals[i + 1]  - 1] = "%s.%s" % (vals[i], vals[i + 1] - 1)
	all_pos["%s.%s" % (vals[i], (vals[i + 1] - 1))] = []
	intervals.append("%s.%s" % (vals[i], vals[i + 1] - 1))
#  if mode == "long": 
#    minLength = long_min
#    maxLength = long_max
#  elif mode == "short":
#    minLength = short_min
#    maxLength = short_max
  for line in tss:
    posRange = {}
    split = line.rstrip('\n').split('\t')
    if split[0] not in valid_chroms: 
      continue
    if sim_type:
      chrom = split[0].split('chr')[1]
    else:
      chrom = split[0]
    if len(split) >= 3:
      t_start, strand = int(split[1]), split[2]
    else:
      t_start = int(split[1])
      strand = "+"
    start = t_start - 1000
    if start < 0: start = 0
    end = t_start + 1000
    for interval in intervals:
      posRange[interval] = [0] * 2001
    for read in readIterator(filenames, chrom, start, end):	
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
	  rmid = (int(rend) + int(rstart)) / 2
          rlength = rend-rstart+1
          if tree.overlaps(rlength):
            interval = sorted(tree[rlength])[0].data
            i = rmid
            if start <= i <= end:
              if strand == "+":
                posRange[interval][i - start] += 1
              if strand == "-":
                posRange[interval][-i + end] += 1
    for interval in intervals:
      all_pos[interval].append(posRange[interval])
  for i in all_pos:
    all_pos[i] = np.array(all_pos[i])
  return all_pos

def main():
  matrices = {}
  valid = {}
  tfile = open(sys.argv[1])
  valid_chroms = open(sys.argv[2])
  for line in valid_chroms:
    split = line.split()
    valid[split[0]] = True
  tss = tfile.readlines()
  filename = sys.argv[3]
  label = sys.argv[4]
  sim_type = sys.argv[5]
  #10 bp chunks from min to max
  minimum = int(sys.argv[6])
  maximum = int(sys.argv[7])
  if sim_type == "True":
        sim_type = True
  else:
        sim_type = False
  filenames = [filename]
  matrices = calculatePerBase(filenames, tss, sim_type, valid, minimum, maximum)
  np.savez_compressed("Hmaps_%s" % label, **matrices) 
  tfile.close()
  valid_chroms.close()

if __name__ == "__main__":
  main()

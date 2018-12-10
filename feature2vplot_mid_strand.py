import os,sys,re
import pysam
import random
import numpy as np
from collections import Counter
import gzip
from bx.intervals.intersection import Intersecter, Interval
import random

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

def calculatePerBase(filenames, tss, mode, sim_type, valid_chroms):
  all_pos = {"+":{}, "-":{}}  
  for i in range(30,301):
    all_pos["+"][i] = [0] * 2001  	
    all_pos["-"][i] = [0] * 2001  	
#  print len(all_pos)
#  print len(all_pos[0])
  for line in tss:
    split = line.rstrip('\n').split('\t')
    if split[0] not in valid_chroms: continue
    if sim_type:
      chrom = split[0].split('chr')[1]
      try: 
        int(chrom)
      except:
	if chrom == "X" or chrom =="Y":
	  chrom = chrom
	else:
          continue
    else:
      chrom = split[0]
    if len(split) >= 3:
      t_start, strand = int(split[1]), split[2]
    else: 
      t_start = int(split[1])
      strand = "+"
    start = t_start - 1000
    if start <= 0: continue
    end = t_start + 1000
    maxLength = 300
    minLength = 30
    for read in readIterator(filenames, chrom, start, end):
      if read.mapq < 5: continue	
      if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
      if isSoftClipped(read.cigar): continue
      if read.is_paired:
        if read.mate_is_unmapped: continue
        if read.rnext != read.tid: continue
        if read.is_read1:
          if read.isize == 0: continue
	  if read.is_reverse:
	    stranded = "+"
	  else:
	    stranded = "-"
          rstart = min(read.pos,read.pnext)+1 # 1-based
          rend = rstart+abs(read.isize)-1 # end included
	  rmid = (int(rend) + int(rstart)) / 2
          rlength = rend-rstart+1
#          print rlength
          if minLength <= rlength <= maxLength:
            if rmid >= start and rmid <= end:
              if strand == "+":
#               print rmid - start
                all_pos[stranded][rlength][rmid - start] += 1
              if strand == "-":
 #              print -rmid + end
                all_pos[stranded][rlength][-rmid + end] += 1
  tot_pos = {"+":[], "-":[]}
  for i in reversed(range(30,301)):
    for stranded in tot_pos:
      tot_pos[stranded].append(all_pos[stranded][i])    
  return [np.array(tot_pos["+"]), np.array(tot_pos['-'])]

def main():
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
  if sim_type == "True":
	sim_type = True
  else:
	sim_type = False
# sim_type = True
  filenames = [filename]
  vplots = calculatePerBase(filenames, tss, "long", sim_type, valid)
  np.save("Vplot.plus.mid.%s" % label, vplots[0])
  np.save("Vplot.minus.mid.%s" % label, vplots[1])
  tfile.close()
  valid_chroms.close()

if __name__ == "__main__":
  main()

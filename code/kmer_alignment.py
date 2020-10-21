#!/usr/bin/env python3
import subprocess
import os
import tempfile
import shutil

def overlapping_substring_length(query, subject, min_k=1):
  # get the length of maximum common substring
  # min_k provides a lower bound to check
  
  # check trivial case
  if query == subject:
    return len(query)
    
  # search every substring of the query
  for k in range(len(query)-1, min_k-1, -1):
    for sub_start in range(0, len(query)-k+1):
      if query[sub_start:(sub_start+k)] in subject:
        return k
        
  # no match
  return 0

def group_kmers(kmers):
  # identify K and verify all seqs are length K
  K = len(kmers[0])
  if not all([len(seq) == K for seq in kmers]):
    raise ValueError("All elements of kmers must be the same length")
  
  # initial assortment 
  kmer_groups = [ [kmers[0]] ]  # initial 'founder' group with one member
  for kk in kmers[1:]:
    # try to find a group for this kmer (kk)
    found = False
    for gid in range(len(kmer_groups)):
      # iterate through each kmer in this group trying to find
      # a match with one bp offset
      for gk in kmer_groups[gid]:
        if overlapping_substring_length(kk, gk, K-1) >= K-1:
          kmer_groups[gid].append(kk)
          found = True
          break
      
      # no need to keep searching if group already found
      if found:
        break
        
    # didn't find an existing group, create a new one
    if not found:
      kmer_groups.append([kk])
  
  # pass back through, trying to merge smaller groups with larger ones
  # first, reorder the groups in increasing order of size
  small2large_order = [x[0] for x in sorted(enumerate([len(g) for g in kmer_groups]), key=lambda x: x[1])]
  kmer_groups = [kmer_groups[i] for i in small2large_order]
  for gid in range(len(kmer_groups)-1):
    found = False
    for kk in kmer_groups[gid]:
      for nid in range(len(kmer_groups)-1, gid, -1):
        for nk in kmer_groups[nid]:
          if overlapping_substring_length(kk, nk, K-1) >= K-1:
            kmer_groups[nid].extend(kmer_groups[gid])
            kmer_groups[gid].clear()
            found = True
            break
        if found:
          break
      if found:
        break
  
  # reorder groups in decreasing order of size
  large2small_order = [x[0] for x in sorted(enumerate([len(g) for g in kmer_groups]), key=lambda x: -x[1])]
  return [kmer_groups[i] for i in large2small_order if len(kmer_groups[i]) > 0  ]


def get_abyss_alignment(seqs, abyss_path=None):
  # if abyss path not defined, try to get it from system
  if abyss_path is None:
    abyss_path = shutil.which("abyss-align")
  
  # verify abyss installed and executable
  if abyss_path is None:
    raise FileNotFoundError("Could not find abyss-align on the system path.")
  elif not os.path.exists(abyss_path):
    raise FileNotFoundError("abyss_path ({}) does not exist.".format(abyss_path))
  if not os.access(abyss_path, os.X_OK):
    raise PermissionError("abyss_path ({}) is not executable.".format(abyss_path))
  
  # write sequences to temp FastA file
  tmp_fa = tempfile.NamedTemporaryFile()
  with open(tmp_fa.name, 'w') as f:
    for i, s in enumerate(seqs):
      f.write(">{}\n{}\n".format(i+1, s))
  
  # run abyss-align on the kmers
  try:
    abyss_proc = subprocess.run([abyss_path, tmp_fa.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  finally:
    # delete temp FastA
    tmp_fa.close()
  
  if abyss_proc.returncode != 0:
    raise RuntimeError("abyss failed: {}".format(abyss_proc.stderr.decode("utf-8")))
  
  # extract sequence from output
  aligned_seq = abyss_proc.stdout.decode("utf-8").split('\n')[-4]

  return aligned_seq
  
if __name__ == "__main__":
  import argparse
  import sys
  
  # get arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--abyss_path", type=str, default=None, help="Path to the abyss-align executable")
  parser.add_argument("-f", "--kmer_file", type=str, help="File containing list of kmers (one per line)")
  args = parser.parse_args()
  
  # load kmers
  if args.kmer_file:
    with open(args.kmer_file, 'r') as f:
      kmers = [s.strip() for s in f.readlines()]
  else:
    kmers = [s.strip() for s in sys.stdin.readlines()]
  
  # group kmers and get extended alignments
  kmer_groups = group_kmers(kmers)
  kmer_extended_seqs = [get_abyss_alignment(kmers, args.abyss_path).upper() for kmers in kmer_groups]
  
  # print out aligned seqs in decreasing order of length
  for seq in sorted(kmer_extended_seqs, key=lambda x: -len(x)):
    print(seq)

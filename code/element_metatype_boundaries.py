#!/usr/bin/python3

import TransposableElements as TE

def load_meta_types(filenm):
	# construct dictionary of repName -> meta type
	meta_types = {}
	with open(filenm, 'r') as mf:
		for line in mf:
			(repnm, metanm) = line.rstrip('\n').split('\t')
			meta_types[repnm] = metanm
	return meta_types

def load_elements(filenm, meta_types):
	# read in elements
	elements = {}
	with open(filenm, 'r') as ef:
		for line in ef:
			entry_elem = TE.RMEntryAsTE(line)
			entry_elem.sub[0].type.Meta = meta_types[entry_elem.sub[0].type.Name]
	
			if entry_elem.id in elements:
				elements[entry_elem.id].merge(entry_elem)
			else:
				elements[entry_elem.id] = entry_elem
	return elements

def format_meta_entry(gpos, eid, sub_num, meta, meta_num):
	return "{}\t{}\t{}\t{}\t{}\t{}:{}:{}".format(gpos.chrom, gpos.start, gpos.end, gpos.strand, eid, sub_num, meta, meta_num)
	
def save_element_metatype_boundaries(elements, f):
	if type(f) == str:
		outf = open(f, 'wt')
	
	def write_entry(*args):
		outf.write(format_meta_entry(*args) + "\n")	
	
	for elem_id in sorted(elements.keys()):
		e = elements[elem_id]
		s = e.sub
		curr_sub_meta = s[0].type.Meta
		curr_meta_5 = s[0].pos.five_prime()
		curr_meta_N = 1
		meta_counts = {curr_sub_meta: 1}
		for sid in range(1, len(s)):
			if s[sid].type.Meta != curr_sub_meta:
				# write out previous meta-subelement
				meta_pos = TE.GenomicPosition(e.chrom, curr_meta_5, s[sid-1].pos.three_prime(), e.strand, strict=False)
				write_entry(meta_pos, elem_id, curr_meta_N, curr_sub_meta, meta_counts[curr_sub_meta])
				
				# set up next meta-subelement
				curr_sub_meta = s[sid].type.Meta
				curr_meta_N += 1
				curr_meta_5 = s[sid].pos.five_prime()
	
				if curr_sub_meta not in meta_counts:
					meta_counts[curr_sub_meta] = 0
	
				meta_counts[curr_sub_meta] += 1
	
		meta_pos = TE.GenomicPosition(e.chrom, curr_meta_5, s[-1].pos.three_prime(), e.strand, strict=False)
		write_entry(meta_pos, elem_id, curr_meta_N, curr_sub_meta, meta_counts[curr_sub_meta])
		
	if type(f) == str:
		outf.close()

if __name__ == "__main__":
	import argparse
	import sys
	
	parser = argparse.ArgumentParser()
	parser.add_argument("element_tsv", type=str, help="TSV file of elements, one subelement per line")
	parser.add_argument("meta_tsv", type=str, help="TSV file with two columns: repName, meta-type")
	parser.add_argument("-o", "--outfile", type=str, help="Output file (ignore or use '-' for stdout)")
	args = parser.parse_args()
	
	if not args.outfile or args.outfile == "-":
		args.outfile = sys.stdout
		
	meta_types = load_meta_types(args.meta_tsv)
	elements = load_elements(args.element_tsv, meta_types)
	
	save_element_metatype_boundaries(elements, args.outfile)

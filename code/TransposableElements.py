#!/usr/bin/python3

from itertools import groupby
import pyranges as pr

def RMEntryAsTE(rm_line):
	entry = rm_line.rstrip('\n').split('\t')
	etype = ElementType(entry[5], entry[6], entry[4])
	epos = GenomicPosition(entry[0], int(entry[1]), int(entry[2]), entry[3])
	ereppos = [int(entry[7]), int(entry[8])] if len(entry) >= 9 else None
	eid = int(entry[9]) if len(entry) >= 10 else None
	
	return TransposableElement(Subelement(etype, epos, ereppos), eid)

class ElementType:
	def __init__(self, rep_class, rep_family, rep_name, meta_type=None):
		self.Class = rep_class
		self.Family = rep_family
		self.Name = rep_name
		self.Meta = meta_type

	def __eq__(self, other):
		if isinstance(other, ElementType):
			return self.Class == other.Class and \
				   self.Family == other.Family and \
				   self.Name == other.Name
		return False

	def same_family_as(self, other):
		if isinstance(other, ElementType):
			return self.Class == other.Class and self.Family == other.Family
		return False

	def meta_matches(self, other):
		return self.same_family_as(other) and (self.Meta == other.Meta or self.Meta is None or other.Meta is None)

class GenomicPosition:
	def __init__(self, chrom, start, end, strand="*", strict=True):
		if not isinstance(start, int):
			raise ValueError("start must be an int")
		if not isinstance(end, int):
			raise ValueError("end must be an int")
		if strand not in ["+", "-", "*"]:
			raise ValueError("strand must be +, -, or *")
		if end < start:
			if strict:
				raise ValueError("start must be >= end")
			else:
				start, end = end, start

		self.chrom = chrom
		self.start = start
		self.end = end
		self.strand = strand

	def five_prime(self):
		if self.strand == "-":
			return self.end
		return self.start

	def three_prime(self):
		if self.strand == "-":
			return self.start
		return self.end

	def compatible_with(self, other, exact_strand=True):
		if isinstance(other, GenomicPosition):
			if exact_strand:
				return self.chrom == other.chrom and self.strand == other.strand
			else:
				return self.chrom == other.chrom and (self.strand == other.strand or self.strand == "*" or other.strand == "*")
		return False
	
	def flanked(self, flank):
		return GenomicPosition(self.chrom, self.start - flank, self.end + flank, self.strand)

	def flanked_3(self, flank):
		if self.strand == "-":
			return GenomicPosition(self.chrom, self.start - flank, self.end, self.strand)
		return GenomicPosition(self.chrom, self.start, self.end + flank, self.strand)

	def flanked_5(self, flank):
		if self.strand == "-":
			return GenomicPosition(self.chrom, self.start, self.end + flank, self.strand)
		return GenomicPosition(self.chrom, self.start - flank, self.end, self.strand)

	def __str__(self):
		return str.format("{}\t{}\t{}\t{}", self.chrom, self.start, self.end, self.strand)

	def pr(self):
		return pr.PyRanges(chromosomes=[self.chrom], starts=[self.start], ends=[self.end], strands=[self.strand])	

class Subelement:
	"""
	A subelement of a transposable element.
	Has an ElementType, a GenomicPosition, and optionally a position denoting
	the location of the subelement within the canonical reference for its ElementType.
	"""
	def __init__(self, elem_type, pos, rep_pos):
		self.type = elem_type
		self.pos = pos
		self.rep_pos = rep_pos

	#def __init__(self, rep_class, rep_family, rep_name, \
	#			 chrom, start, end, strand, \
	#			 rep_start, rep_end):
	#	self.type = ElementType(rep_class, rep_family, rep_name)
	#	self.pos = GenomicPosition(chrom, start, end, strand)
	#	self.rep_pos = [rep_start, rep_end]

	def compatible_with(self, other, check_meta=True):
		if isinstance(other, Subelement) and self.pos.compatible_with(other.pos) \
		   and (not check_meta or self.type.meta_matches(other.type)):
			return True
		return False
			
	def fuse(self, other):
		"""
		Fuse two subelements together.
		
		Args:
			other: Another Subelement.   
		"""
		if not isinstance(other, Subelement):
			raise TypeError("Attempted to fuse with non-Subelement")
		if not self.compatible_with(other):
			raise ValueError("Attempted to merge incompatible Subelements")

		self.pos.start = min(self.pos.start, other.pos.start)
		self.pos.end   = min(self.pos.end,   other.pos.end)
		
		if self.rep_pos is not None and other.rep_pos is not None and self.type.Name == other.type.Name:
			self.rep_pos = [min(self.rep_pos[0], other.rep_pos[0]), max(self.rep_pos[1], other.rep_pos[1])]
		else:
			self.rep_pos = None
	

def reduce_runs(l):
	return [x[0] for x in list(groupby(l))]		
	
# TODO take care of order 
# enforce order of subelements by pos.start
# backwards if - strand
class TransposableElement:
	def __init__(self, subelements, elem_id=None):
		if not isinstance(subelements, list):
			if isinstance(subelements, Subelement):
				subelements = [subelements]
			else:
				raise TypeError("subelements must either be a list of Subelements or a single Subelement")
		else:
			if not all([isinstance(s, Subelement) for s in subelements]):
				raise TypeError("subelements is a list containing at least one non-Subelement")

		if len(subelements) > 1 and not all([subelements[0].compatible_with(s, check_meta=False) for s in subelements[1:]]):
			raise ValueError("Attempted to create TransposableElement with incompatible Subelements")

		self.strand = subelements[0].pos.strand
		self.chrom = subelements[0].pos.chrom
		self.sub = TransposableElement.arrange_subelements(subelements, self.strand)
		self.id = elem_id

	@staticmethod
	def arrange_subelements(subs, strand=None):
		if strand is None:
			strand = s[0].pos.strand
		return sorted(subs, key=lambda s: s.pos.five_prime(), reverse=(strand == "-"))

	def compatible_with(self, other):
		if isinstance(other, TransposableElement):
			return self.sub[0].compatible_with(other.sub[0], check_meta=False)
		return False

	def names(self):
		return [s.type.Name for s in self.sub]
	
	def meta(self):
		return reduce_runs([s.type.Meta for s in self.sub])

	def meta_str(self):
		return "-".join(self.meta())

	def merge(self, other):
		if not isinstance(other, TransposableElement):
			raise TypeError("Attempting to merge a non-TransposableElement")
		if not self.compatible_with(other):
			raise ValueError("Attempted to merge incompatible transposable elements")

		self.sub = TransposableElement.arrange_subelements(self.sub + other.sub, strand=self.strand)
	
	def first(self):
		return self.sub[0]

	def last(self):
		return self.sub[-1]

	def span(self):
		return GenomicPosition(self.chrom, self.first().pos.five_prime(), self.last().pos.three_prime(), self.strand, strict=False)

class ERV(TransposableElement):
	def missing_5prime_LTR(self):
		return self.first().type.Meta != "LTR"
	
	def missing_3prime_LTR(self):
		return self.last().type.Meta != "LTR"

	def is_solo_LTR(self):
		return self.meta() == ["LTR"]
	
	def is_complete(self):
		return not (self.missing_5prime_LTR() or self.missing_3prime_LTR())

	def is_fully_structured(self):
		return self.is_complete() and "internal" in self.meta()

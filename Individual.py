#!/usr/bin/env python
# encoding: utf-8
"""
Individual.py

Created by Joshua Shapiro on 2008-08-11.
"""
from __future__ import division
import operator
import itertools

from numpy import *
from Chromosomes import *

def newChromosomes(parent, n = None, cM = 200, chrNames=None, interference = "absent"):
    """Generates a list of n chromosomes, each with length cM (possibly a list, in which case n is optional)"""
    if parent == None:
        raise ValueError, "Must specify parent identifier"
    if n == None and type(cM) == type(1):
        raise ValueError, "Must specify number of chromosomes"
    if type(cM) == type(1):
        cM = [cM for c in range(n)]
    elif len(cM) != n  and n != None:
        raise ValueError, "List of cM lengths must equal n, if specified"
    if chrNames == None:
        chrNames = [x+1 for x in range( len(cM) )]
    elif len(chrNames) != len(cM):
        raise ValueError, "chrNames must equal number of chrs requested"
    
    return [Chromosome(cM = cM, name = name, newParent = parent, interference = interference) for name, cM in itertools.izip(chrNames, cM)]

def alleleCompare(a,b, P1):
  if a == b :
    if a == P1:
      genotype = 2
    else:
      genotype = 0
  else:
    genotype = 1
  return genotype

class Diploid(object):
    """ Diploid individual, monoecious"""
    def __init__(self,  name = None, chromosome_set = None, newChr = None, cM = 200, chrNames = None):
        self.name = name
        if (chromosome_set == None and newChr == None) or (chromosome_set != None and newChr != None):
              raise ValueError, "Must specify only one of either the number of new Chromosomes to add or sets of the chromosomes themselves."
        elif chromosome_set != None:
            self.chromosome_set = chromosome_set
        else:
            self.chromosome_set = (newChromosomes(parent = self.name, n = newChr, cM = cM, chrNames = chrNames),
                                   newChromosomes(parent = self.name, n = newChr, cM = cM, chrNames = chrNames) )
        for chr_list in self.chromosome_set:
            chr_list.sort(key= operator.attrgetter("name"))
       
    def getNChr(self):
        return len(self.chromosome_set[0])
    nChr = property(fget= getNChr, doc = "Number of chromosomes")
       
    def make_gamete(self):
        gamete = []
        for maternal, paternal in itertools.izip (self.chromosome_set[0], self.chromosome_set[1]):
            gamete.append(maternal.recombine(paternal)[0])
        return gamete

    def mate(self, partner, nOffspring = 1):
        offspring = [Diploid(chromosome_set = (self.make_gamete(), partner.make_gamete()))   
                     for _ in xrange(nOffspring) ]
        return offspring
    
    def getAllGenos(self, interval= 1, cM = True, reference = None):
      """Gets the genotypes for evenly spaced intervals.  Spacing is either by cM or if cM=False, by fractions of the genome.
         The return is the locations, with the genotype coded as 0,1,2 as the number of alleles equal to the reference (arbitrary if not specified)."""
      if not reference:
        reference = self.chromosome_set[0][0].getParentAtLocation(0)
      genotype = []
      if cM:
        for chrom in range(self.nChr):
          chrLength = self.chromosome_set[0][chrom].cM
          nPos = int(chrLength//interval)
          positions = [i * interval for i in xrange(nPos)] + [chrLength]
          allelesA = self.chromosome_set[0][chrom].getParentAtMapLocs(positions)
          allelesB = self.chromosome_set[1][chrom].getParentAtMapLocs(positions)
          loci = [ (chrom, pos, alleleCompare(a,b,reference) )
                  for pos, a, b 
                  in itertools.izip(positions, allelesA, allelesB)]
          genotype += loci
      else:
        nPos = int(1//interval)
        positions = [i * interval for i in xrange(nPos)] + [1.0]
        for chrom in range(self.nChr):
          allelesA = self.chromosome_set[0][chrom].getParentAtLocations(positions)
          allelesB = self.chromosome_set[1][chrom].getParentAtLocations(positions)
          loci = [ (chrom, pos, alleleCompare(a,b,reference) )
                  for pos, a, b 
                  in itertools.izip(positions, allelesA, allelesB)]
          genotype += loci
      return genotype      
                  
                
        
class Haploid(object):
    """Haploid individual, monoecious"""
    def __init__(self, name = None, chromosomes = None, newChr = None, cM = 200, chrNames = None):
        self.name = name
        if (chromosomes == None and newChr == None) or (chromosomes != None and newChr != None):
            raise ValueError, "Must specify only one of either the number of new Chromosomes to add or a list of the chromosomes themselves."
        elif chromosomes != None:
            self.chromosomes = chromosomes
        elif newChr != None:
            self.chromosomes = newChromosomes(parent = self.name, n = newChr, cM = cM, chrNames = chrNames)
        self.chromosomes.sort(key= operator.attrgetter("name"))
    
    def getNChr(self):
        return len(self.chromosomes)
    nChr = property(fget= getNChr, doc = "Number of chromosomes")
    
    def mate(self, mate, nOffspring = 1):
        if type(mate) != type (self):
            raise TypeError, "Haploid individuals can only mate with other haploids"
        if self.nChr != mate.nChr:
            raise ValueError, "Individuals to mate must have the same number of chromosomes"
        
        if nOffspring == 1:
            offspring = Haploid( chromosomes = [ selfChr.recombine(mateChr)[0] for selfChr, mateChr in itertools.izip(self.chromosomes, mate.chromosomes) ] )
        else:
            offspring = [Haploid( chromosomes = [ selfChr.recombine(mateChr)[0] for selfChr, mateChr in itertools.izip(self.chromosomes, mate.chromosomes) ] )
                        for __ in range(nOffspring)]
        return offspring
    
        
if __name__ == '__main__':
    Aparent = Diploid(name = "A", newChr = 10)
    Bparent = Diploid(name = "B", newChr = 10)
    F1s = Aparent.mate(Bparent, nOffspring = 4)
    F2s = F1s[0].mate(F1s[1], nOffspring = 4)
    
    print "Offspring Chromosomes:"
    print "F1 Ind 1: %s" % F1s[0].name 
    for chr in F1s[0].chromosome_set[0]:
        print "Chr %s: %s" % (chr.name, chr.segments)
    print "F2 Ind 2"
    for chr in F1s[1].chromosome_set[1]:
        print "Chr %s: %s" % (chr.name, chr.segments)
        
    for chr in F2s[1].chromosome_set[1]:
        print "Chr %s: %s" % (chr.name, chr.segments)
        
    print F2s[0].getAllGenos(interval = 50, reference = "A")

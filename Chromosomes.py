#!/usr/bin/env python
# encoding: utf-8
"""
Chromosomes.py

Created by Joshua Shapiro on 2008-08-11.
"""
from __future__ import division 
import itertools
import operator


from numpy import *


def generateBreaksPoisson(cM = 200):
    breaks = random.uniform(size = random.poisson(cM/100.0))
    breaks.sort()
    return breaks


class Chromosome(object):
    """Chromosome object which contains information on parentage of segments"""
    def __init__(self,  cM=200, name=None, segments=None, newParent=None, interference = "absent"):
        super(Chromosome, self).__init__()
        self.name = name
        self.cM = cM
        self.segments = segments
        #segments are lists of tuples where the first position is the start location, the second is the parent of origin
        self.interference = interference
        
        if newParent != None:
            if self.segments != None:
                raise ValueError, "Can't set segment information in a new parent"
            self.segments = [(0,newParent)]
        
        if not (self.interference in ["complete", "absent"]):
          raise ValueError, "Interference must be one of 'complete' or 'absent'."
        
    def __eq__(self, other):
      if  isinstance(other, Chromosome):
        if self.name != other.name:
          return False
        if len(self.segments) != len(other.segments):
          return False
        for s1, s2 in itertools.izip(self.segments, other.segments):
          if s1 != s2: return False
        return True  
      else: return NotImplemented
    
    def __ne__(self, other):
      equal_result = self.__eq__(other)
      if (equal_result is not NotImplemented):
          return not equal_result
      return NotImplemented
    
    def getParentAtLocation(self, loc):
        """gets the Parental Identity for a chromosomal location"""
        if loc < 0 or loc > 1:
            raise ValueError, "Location must be in range [0,1]"
            
        i = 0
        while i < len(self.segments) and self.segments[i][0] < loc :
            i+=1
        return self.segments[i-1][1]
    
    def getParentAtLocations(self, locs):
        """gets the Parental Identity for a list of chromosomal locations"""
        parents = [''] * len(locs)
        order = [i for _,i in sorted(itertools.izip( locs, range(len(locs)) ))]
        if locs[order[0]] < 0 or locs[order[-1]] > 1 :
            raise ValueError, "Locations must be in range [0,1]"
        i = 0
        for n in order:
          while i < len(self.segments) and self.segments[i][0] < locs[n] :
            i+=1
          parents[n] = self.segments[i-1][1]
        return parents
    
    
    def getParentAtMapLoc(self, mapLoc):
      """gets the Parental identity for a given cM position"""
      if mapLoc < 0 or mapLoc > self.cM:
        raise ValueError, "Map location must be withing the range of the chromosome."
      loc = cM/float(self.cM)
      return self.getParentAtLocation(loc)
      
    def getParentAtMapLocs(self, mapLocs):
        """gets the Parental Identity for a list of chromosomal locations"""
        parents = [''] * len(mapLocs)
        locs = [ml/float(self.cM) for ml in mapLocs]
        order = [i for _,i in sorted(itertools.izip( locs, range(len(locs)) ))]
        if locs[order[0]] < 0 or locs[order[-1]] > 1 :
            raise ValueError, "Locations must be in range [0,1]"
        i = 0
        for n in order:
          while i < len(self.segments) and self.segments[i][0] < locs[n] :
            i+=1
          parents[n] = self.segments[i-1][1]
        return parents
    
    def recombine(self, mate, interference = None):
        if self.name != mate.name:
            raise ValueError, "Chromosome names are not the same; can't recombine between them." 
        if self == mate:
          #shortcut: any recombinants would be identical anyway
          return (Chromosome(name = self.name, cM = self.cM,  segments = self.segments), 
                  Chromosome(name = self.name, cM = self.cM,  segments = self.segments))
        
        if method == None:
          method = self.recombMethod
        
        segments1 = list(self.segments)
        segments2 = list(mate.segments)
        brokenSegments1 = list()
        brokenSegments2 = list()
        
        if interference == "absent":
          crossOvers = generateBreaksPoisson(self.cM)
        elif interference == "complete":
          crossOvers = random.uniform(size=1)
        else:
          raise ValueError, "Interference setting must be one of 'absent' or 'complete'."
        
        for crossOver in crossOvers:
            tempSeg = list()
            #move segments to temp list until you get past the crossover
            while len(segments1) > 0 and segments1[0][0] < crossOver:
                tempSeg.append(segments1.pop(0))
            #add in new start segment
            segments1.insert(0, (crossOver, tempSeg[-1][1]))
            brokenSegments1.append(tempSeg)
            
            tempSeg=[]
            while len(segments2) > 0 and segments2[0][0] < crossOver:
                tempSeg.append(segments2.pop(0))
            segments2.insert(0, (crossOver, tempSeg[-1][1]))
            brokenSegments2.append(tempSeg)
        #all breakpoints done, now just add the remainder of the segments to our lists
        brokenSegments1.append(segments1)
        brokenSegments2.append(segments2)
        
        #combine breakpoint lists, alternating parental lists
        chr1 = list()
        chr2 = list()
        for index, bs1 in enumerate(brokenSegments1):
            if index % 2 == 0 : 
                chr1 += bs1 
                chr2 += brokenSegments2[index]
            else:
                chr1 += brokenSegments2[index]
                chr2 += bs1
        
        #remove redundant segments
        chr1 = [x for i,x in enumerate(chr1) if (i == 0 or chr1[i][1] != chr1[i-1][1])]
        chr2 = [x for i,x in enumerate(chr2) if (i == 0 or chr2[i][1] != chr2[i-1][1])]
        
        
        if random.binomial(1,0.5): #randomly order xover products
            return (Chromosome(name = self.name, cM = self.cM,  segments = chr1), 
                    Chromosome(name = self.name, cM = self.cM,  segments = chr2))
        else:
            return (Chromosome(name = self.name, cM = self.cM,  segments = chr2), 
                    Chromosome(name = self.name, cM = self.cM,  segments = chr1))            


        
if __name__ == '__main__':
    a = Chromosome(newParent = "Blue")
    b = Chromosome(segments = [(0,1)])
    newChrs = a.recombine(b)
    print newChrs[0].segments
    newChrs = newChrs[0].recombine(newChrs[1])
    print newChrs[0].segments
    print newChrs[0].getParentAtLocation(0.5)
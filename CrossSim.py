#!/usr/bin/env python
# encoding: utf-8
"""
CrossSim.py

Simulates crosses among haploid individuals.  8 parent crosses with various crossing schemes.

Created by Joshua Shapiro on 2008-10-06.
"""
from __future__ import division
import sys
import os
import optparse


from numpy import *
import Crosses
import Individual

def getOptions():
  """Get command line options"""
  parser = optparse.OptionParser()
  parser.add_option("-r", "--replicates", 
                    action = "store", dest = "reps",
                    type = "int", default = 1,
                    help = "Number of simulations to run")
  parser.add_option("-s","--scheme",
                    action = "store", dest = "scheme", default = 'collab',
                    help = "Crossing scheme. One of: 'collab' (collaborative cross), 'round' (round robin crosses), 'random' (random crossing)")
  parser.add_option("-g", "--generations",
                    action = "store", dest = "generations", 
                    type = "int", default = 0,
                    help = "Number of generations fo simulate. For Collaborative cross, this is the number of additional random mating generations")
  parser.add_option("-n", "--nInd", 
                    action = "store", dest = "nInd", 
                    type = "int", default = None,
                    help = "Number of individuals to simulate for each generation.")

  parser.add_option("-f", "--finalInd", 
                    action = "store", dest = "fInd", 
                    type = "int", default = 315,
                    help = "Final number of individuals to simulate.  If it is greater than the number of individuals present in the last generation, an additional round of crosses may be added")
  (options, args) = parser.parse_args()
  return options

class SimStats(object):
  """structure for holding simulation statistics, filled dynamically"""
  def __init__(self):
    pass
    

def main():
  options = getOptions()
  parentIDs = [0,1,2,3,4,5,6,7]
  parents = [Individual.Haploid(name = p, newChr = 1) for p in parentIDs ] 
  allStats = list()
   
  for i in range(0, int(options.reps)):
    stats = SimStats()
    if options.scheme == 'collab':
      population = Crosses.collabCross(parents)
      if options.generations > 0:
        for i in range(0, options.generations - 1):
          population =  Crosses.randomCross(population, nOffspring = options.nInd)
        population = Crosses.randomCross(population, nOffspring = options.fInd)
    
    elif options.scheme == 'random':
      
      if options.generations == 0:
        gens = 1
      else:
        gens = options.generations
      population = Crosses.abaCross(parents)
      for i in range(0, gens - 1):
        population = Crosses.randomCross(population, nOffspring = options.nInd)
      population = Crosses.randomCross(population, nOffspring = options.fInd)
    
    elif options.scheme == 'bigRand':
      if options.generations == 0:
        gens = 1
      else:
        gens = options.generations
      population = Crosses.randInfCross(parents, options.fInd, generations = gens)
    #crosses done, data in population list
    #calculate stats for each cross
    
    #longest unrecombined segment
    breaks = list()
    for ind in population:
      breaks += [s[0] for s in ind.chromosomes[0].segments]
    breaks.sort()
    breaks = [x for i,x in enumerate(breaks) if (i == 0 or breaks[i] != breaks[i-1])]
    breaks.append(1.0)
    lastbp = 0
    maxSeg = 0
    segs = list()
    for bp in breaks:
      segLen = float(bp - lastbp)
      lastbp = bp
      segs.append(segLen * 200)
    
    stats.medSeg = median(segs)
    stats.meanSeg = mean(segs)
    stats.varSeg = var(segs)
    stats.maxSeg = max(segs)
    
    #proportion of genome covered for each parent
    sites = range(0,201)
    freqs = [[0 for _ in parents] for _ in sites]
    for i in sites:
      for ind in population:
        #increment count at that parental location
        freqs[i][ind.chromosomes[0].getParentAtLocation(float(i)/200)] += 1
    
    missingSegment = 0
    anyMissing = 0
    lowSegment = 0
    anyLow = 0
    for site in freqs:
      hasMissing = 0
      hasLow = 0
      for parent in site:  
        if parent == 0:
          missingSegment += 1
          hasMissing = 1
        if parent < 20:
          lowSegment += 1
          hasLow = 1
      anyMissing += hasMissing
      anyLow += hasLow
        
    stats.missing = missingSegment/len(freqs) * len(parents)
    stats.low = lowSegment/len(freqs) * len (parents)
    stats.anyMissing = anyMissing/len(freqs)
    stats.anyLow = anyLow/len(freqs)
    
    allStats.append(stats)
  
  print "median\tmean\tvar\tmax\tFracMissing\tFracLowFreq\ttotFracMissing\ttotFracLow" 
  for stat in allStats:
    print "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (stat.medSeg, stat.meanSeg, stat.varSeg,  stat.maxSeg, stat.anyMissing, stat.anyLow, stat.missing, stat.low)     

if __name__ == '__main__':
  main()


#!/usr/bin/env python
# encoding: utf-8
"""
Crosses.py

Created by Joshua Shapiro on 2008-08-11.
Copyright (c) 2008 Princeton University. All rights reserved.
"""

import Individual
from numpy import *
import itertools

def randomCross(population, nOffspring = None):
    if nOffspring == None:
        nOffspring = len(population)
    r1 = random.randint(0, len(population), size = nOffspring)
    r2 = random.randint(0, len(population), size = nOffspring)
    newPop = [ population[a].mate(population[b]) for a,b in itertools.izip(r1, r2) ]
    return newPop


def collabCross(population, reciprocal = False):
    """Collaborative cross scheme.  Optimized for 2^n individuals"""
    #first cross
    C1 = list()
    for ind1 in population:
        for ind2 in population:
            if (reciprocal and ind1 != ind2) or ind1.name < ind2.name:
                newInd = ind1.mate(ind2)
                newInd.name = len(C1)
                newInd.parents = (ind1.name, ind2.name)
                C1.append(newInd)
    #remaining crosses
    while len(C1[0].parents) <= len(population)/2 :
        C2 = list()
        print len(C1)
        for ind1 in C1:
            for ind2 in C1:
                if (reciprocal and ind1 != ind2) or ind1.name < ind2.name:
                    #check for shared parents
                    share = False
                    for par in ind1.parents:
                        if par in ind2.parents:
                            share = True
                            break
                    if not share and (reciprocal or ind1.name < ind2.name):
                        newInd = ind1.mate(ind2)
                        newInd.name = len(C2)
                        newInd.parents = ind1.parents + ind2.parents
                        C2.append(newInd)
        C1 = list(C2)

    return(C2)

def rrCross(population, reciprocal = False):
    """Round robin cross for a population of individuals"""
    shiftedPop = list(population)
    shiftedPop.insert(0, shiftedPop.pop())
    offspring = [a.mate(b) for a,b in itertools.izip(population, shiftedPop)]
    if reciprocal:
        offspring += [b.mate(a) for a,b in itertools.izip(population, shiftedPop)]
    return offspring

def abaCross(population, reciprocal = False):
    """All by All crossing scheme for a population"""
    if reciprocal:
        offspring = [a.mate(b) for a in population for b in population if a != b]
    else:
        offspring = [a.mate(b) for a in population for b in population if a.name < b.name]
    return offspring

def randInfCross(population, nOffspring, generations):
  """Random cross for near infinite population sizes: all segregants are independent at every generation"""
  gen1Size = nOffspring * (2 ** (generations -1))
  offspring = randomCross(population, nOffspring = gen1Size)
  for _ in range(0, generations - 1):
    offspring = [offspring[i].mate(offspring[i + 1]) for i in range(0, len(offspring), 2)]
  return offspring

if __name__ == '__main__':
    myPop = [Individual.Haploid(name = x, newChr = 2) for x in ["a","b","c","d","e","f","g", "h"] ] 
    nextGen = randomCross(myPop, nOffspring = 10)
    
    print "Offspring:"
    for i, ind in enumerate(nextGen):
        print "Ind %s:" %i
        for chr in ind.chromosomes:
            print "Chr %s: %s" % (chr.name, chr.segments)
    
    myPop = [Individual.Haploid(name = x, newChr = 0) for x in range(16) ] 
    ccr = collabCross(myPop)
    print "Collaborative Cross: NInd = %d" % len(ccr)
    for i in range(10):
        print ccr[i].parents
    for chr in ccr[0].chromosomes:
        print "Chr %s: %s" % (chr.name, chr.segments)  
	
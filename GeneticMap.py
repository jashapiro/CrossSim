#!/usr/bin/env python
# encoding: utf-8
"""
GeneticMap.py

Created by Joshua Shapiro on 2008-10-02.
"""
from __future__ import division
import sys
import os

class GeneticMap(object):
  """Defines a relationship between the physical and genetic maps for each marker position"""
  def __init__(self, markers = [], mapLength = None, physLength = None):
    self.markers = list(markers)
    if len(self.markers) > 0:
      chrom = self.markers[0].chrom
    for m in self.markers:
      if m.chrom != chrom:
        raise ValueError, "All markers in a chromosome map must be on the same chromosome"
    self.markers.sort(key = operator.attrgetter("cM"))
    if mapLength:
      self.mapLength = mapLength
    elif len(self.markers) > 0:
      mapLength = self.markers[-1].cM
    else:
      pass
    if physLength:
      self.physLength = physLength:
    elif len(self.markers) > 0:
      physLength = self.markers[-1].bp
    else:
      pass
  
  nMarkers = property(fget = lambda self: len(self.markers), doc= "Number of markers in the map")
  markerNames = property (fget = lambda self: [m.name for m in self.markers], doc = "List of the marker names") 
  
  def addMarkers(self, markers):
    """Adds a marker to the chromosomal map"""
    for marker in markers:
      if marker.chrom != self.markers[0].chrom:
        raise ValueError, "All markers in a chromosome map must be on the same chromosome"
    self.markers += markers
    self.markers.sort(key = operator.attrgetter("cM"))
  
  def getMarkerMapPosition(self, marker):
    """Get the map postion in centimorgans of the named marker"""
    if marker in self.markerNames:
      return self.markers[self.markerNames.index(marker)].cM
    else:
      return None
  
  def getMarkerPhysPosition(self, marker):
    """Get the map postion in centimorgans of the named marker"""
    if marker in self.markerNames:
      return self.markers[self.markerNames.index(marker)].bp
    else:
      return None
  
  def getMarkersByPhysicalRange(self, start, end):
    """returns a list of markers in a sequence range"""
    return [m for m in self.markers if (m.bp >= start and m.bp <= end)]
  def getMarkersByMapRange(self, start, end):
    """returns a list of markers in a map distance (cM) range""" 
    return [m for m in self.markers if (m.cM >= start and m.cM <= end)]
        

class Marker(object):
  """Defines the location of a marker both by name, chromosome and genetic/physical positions"""
  def __init__(self, name, chrom, cM, bp = None):
    self.name = name
    self.chrom = chrom
    self.cM = cM
    self.bp = bp

class AlleleMap(object):
  """Mapping of allele values to marker names for a single individual"""
  def __init__(self, alleleMap):
    super(AlleleMap, self).__init__()
    self.map = alleleMap
    if self.map != None and  type(self.map) != type (dict()):
      except(TypeError, "alleleMap must be a dictionary of dictionaries")
  
  def addMarker(self, marker, markerValue):
    if type(marker) == type(Marker()):
      markerName = marker.name
    else:
      markerName = marker
    if markerName in self.map and self.map[markerName] != markerValue:
      raise ValueError, "Marker '%s' is already assigned to a different value in the allele map" % (markerName)
    else:
      self.map[markerName] = markerValue
  
  def getMarker(self, marker):
    if type(marker) == type(Marker()):
      markerName = marker.name
    else:
      markerName = marker
    return self.map[markerName]


      
    
    

if __name__ == '__main__':
  pass
#!/usr/bin/python
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of DrNUM.                                          +
# +                                                                      +
# + Copyright 2013 numrax GmbH, enGits GmbH                              +
# +                                                                      +
# + DrNUM is free software: you can redistribute it and/or modify        +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + DrNUM is distributed in the hope that it will be useful,             +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Patch:

  def __init__(self):
    self.x1 = 0.
    self.y1 = 0.
    self.z1 = 0.
    self.x2 = 0.
    self.y2 = 0.
    self.z2 = 0.
    self.Ni = 0.
    self.Nj = 0.
    self.Nk = 0.
    self.hi = 1.e99
    self.hj = 1.e99
    self.hk = 1.e99
    self.SeekLayersI1 = 0
    self.SeekLayersI2 = 0
    self.SeekLayersJ1 = 0
    self.SeekLayersJ2 = 0
    self.SeekLayersK1 = 0
    self.SeekLayersK2 = 0
    self.index = -1
    self.neighI1 = set([])
    self.neighI2 = set([])
    self.neighJ1 = set([])
    self.neighJ2 = set([])
    self.neighK1 = set([])
    self.neighK2 = set([])
    self.overlap_factor = 1.0;
    self.neigh_overlap_factor = 0.0;
    self.overlap = 2
    self.name = "N/A"
    self.i1_flux = "far"
    self.i2_flux = "far"
    self.j1_flux = "far"
    self.j2_flux = "far"
    self.k1_flux = "far"
    self.k2_flux = "far"
    self.dim_set = False

  def setOverlap(self, overlap):
    self.overlap = overlap

  def setName(self, name):
    self.name = name

  def Li(self):
    return self.x2 - self.x1

  def Lj(self):
    return self.y2 - self.y1

  def Lk(self):
    return self.z2 - self.z1

  def addNeighI1(self, neigh):
    self.neighI1.add(neigh.index)
    self.SeekLayersI1 = self.overlap
    neigh.neighI2.add(self.index)
    neigh.SeekLayersI2 = neigh.overlap

  def addNeighI2(self, neigh):
    self.neighI2.add(neigh.index)
    self.SeekLayersI2 = self.overlap
    neigh.neighI1.add(self.index)
    neigh.SeekLayersI1 = neigh.overlap

  def addNeighJ1(self, neigh):
    self.neighJ1.add(neigh.index)
    self.SeekLayersJ1 = self.overlap
    neigh.neighJ2.add(self.index)
    neigh.SeekLayersJ2 = neigh.overlap

  def addNeighJ2(self, neigh):
    self.neighJ2.add(neigh.index)
    self.SeekLayersJ2 = self.overlap
    neigh.neighJ1.add(self.index)
    neigh.SeekLayersJ1 = neigh.overlap

  def addNeighK1(self, neigh):
    self.neighK1.add(neigh.index)
    self.SeekLayersK1 = self.overlap
    neigh.neighK2.add(self.index)
    neigh.SeekLayersK2 = neigh.overlap

  def addNeighK2(self, neigh):
    self.neighK2.add(neigh.index)
    self.SeekLayersK2 = self.overlap
    neigh.neighK1.add(self.index)
    neigh.SeekLayersK1 = neigh.overlap

  def addTopNeigh(self, neigh):
    self.addNeighK2(neigh)

  def addBottomNeigh(self, neigh):
    self.addNeighK1(neigh)

  def addLeftNeigh(self, neigh):
    self.addNeighI1(neigh)

  def addRightNeigh(self, neigh):
    self.addNeighI2(neigh)

  def addFrontNeigh(self, neigh):
    self.addNeighJ1(neigh)

  def addBackNeigh(self, neigh):
    self.addNeighJ2(neigh)

  def setRes(self, h):
    h = float(h)
    self.hi = h
    self.hj = h
    self.hk = h

  def setResI(self, h):
    self.hi = float(h)

  def setResJ(self, h):
    self.hj = float(h)

  def setResK(self, h):
    self.hk = float(h)

  def setBox(self, x1, y1, z1, x2, y2, z2):
    self.x1 = float(x1)
    self.y1 = float(y1)
    self.z1 = float(z1)
    self.x2 = float(x2)
    self.y2 = float(y2)
    self.z2 = float(z2)

  def updateResolution(self):
    self.Ni = int((self.x2 - self.x1)/self.hi)
    self.Nj = int((self.y2 - self.y1)/self.hj)
    self.Nk = int((self.z2 - self.z1)/self.hk)
    while self.Ni*self.hi < self.x2 - self.x1:
      self.Ni += 1
    while self.Nj*self.hj < self.y2 - self.y1:
      self.Nj += 1
    while self.Nk*self.hk < self.z2 - self.z1:
      self.Nk += 1

  def setDim(self, ni, nj, nk):
    self.Ni = ni
    self.Nj = nj
    self.Nk = nk
    self.hi = (self.x2 - self.x1)/float(self.Ni)
    self.hj = (self.y2 - self.y1)/float(self.Nj)
    self.hk = (self.z2 - self.z1)/float(self.Nk)
    self.dim_set = True

  def inflateI1(self, hn):
    #self.x1 -= self.overlap_factor*(0.5*hn + 2*self.hi)
    self.x1 -= self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hi

  def inflateI2(self, hn):
    #self.x2 += self.overlap_factor*(0.5*hn + 2*self.hi)
    self.x2 += self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hi

  def inflateJ1(self, hn):
    #self.y1 -= self.overlap_factor*(0.5*hn + 2*self.hj)
    self.y1 -= self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hj

  def inflateJ2(self, hn):
    #self.y2 += self.overlap_factor*(0.5*hn + 2*self.hj)
    self.y2 += self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hj

  def inflateK1(self, hn):
    #self.z1 -= self.overlap_factor*(0.5*hn + 2*self.hk)
    self.z1 -= self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hk

  def inflateK2(self, hn):
    #self.z2 += self.overlap_factor*(0.5*hn + 2*self.hk)
    self.z2 += self.neigh_overlap_factor*hn + self.overlap*self.overlap_factor*self.hk

  def toString(self):
    txt = "1001 // index=%d name='%s'\n{\n" % (self.index, self.name)
    txt += "  " + str(self.x1) + " " + str(self.y1) + " " + str(self.z1) + "\n"
    txt += "  1 0 0\n"
    txt += "  0 1 0\n"
    txt += "  1\n"
    txt += "  " + str(self.Ni) + " " + str(self.Nj) + " " + str(self.Nk) + "\n"
    txt += " "
    txt += " " + str(self.SeekLayersI1);
    txt += " " + str(self.SeekLayersI2);
    txt += " " + str(self.SeekLayersJ1);
    txt += " " + str(self.SeekLayersJ2);
    txt += " " + str(self.SeekLayersK1);
    txt += " " + str(self.SeekLayersK2);
    txt += "\n"
    txt += "  " + str(self.x2 - self.x1) + " " + str(self.y2 - self.y1) + " " + str(self.z2 - self.z1) + "\n"
    txt += "  fx fy fz\n"
    txt += " "
    if len(self.neighI1) == 0:
      txt += " " + self.i1_flux
    else:
      txt += " 0"
    if len(self.neighI2) == 0:
      txt += " " + self.i2_flux
    else:
      txt += " 0"
    if len(self.neighJ1) == 0:
      txt += " " + self.j1_flux
    else:
      txt += " 0"
    if len(self.neighJ2) == 0:
      txt += " " + self.j2_flux
    else:
      txt += " 0"
    if len(self.neighK1) == 0:
      txt += " " + self.k1_flux
    else:
      txt += " 0"
    if len(self.neighK2) == 0:
      txt += " " + self.k2_flux
    else:
      txt += " 0"
    txt += "\n"
    txt += "  0\n"
    txt += "}\n"
    return txt

  def numCells(self):
    return self.Ni*self.Nj*self.Nk

  def numOverlapCells(self):
    N = 0
    N += self.SeekLayersI1*self.Nj*self.Nk
    N += self.SeekLayersI2*self.Nj*self.Nk
    N += self.SeekLayersJ1*self.Ni*self.Nk
    N += self.SeekLayersJ2*self.Ni*self.Nk
    N += self.SeekLayersK1*self.Ni*self.Nj
    N += self.SeekLayersK2*self.Ni*self.Nj
    return N

  def setI1Flux(self, flux):
    self.i1_flux = flux

  def setI2Flux(self, flux):
    self.i2_flux = flux

  def setJ1Flux(self, flux):
    self.j1_flux = flux

  def setJ2Flux(self, flux):
    self.j2_flux = flux

  def setK1Flux(self, flux):
    self.k1_flux = flux

  def setK2Flux(self, flux):
    self.k2_flux = flux

  def pointInside(self, x, y, z):
    point_inside = True
    if x < self.x1 or x > self.x2:
      point_inside = False
    if y < self.y1 or y > self.y2:
      point_inside = False
    if z < self.z1 or z > self.z2:
      point_inside = False
    return point_inside



class Mesh:

  def __init__(self):
    self.patches = []

  def findPatch(self, x, y, z):
   for patch in self.patches:
     if patch.pointInside(x, y, z):
       return patch.index
   return -1

  def printPatches(self):
    for patch in self.patches:
      print "idx", patch.index, "x1", patch.x1, "x2", patch.x2, "y1", patch.y1, "y2", patch.y2, "z1",patch.z1, "z2", patch.z2

  def getPatch(self, index):
    return patches[index]

  def createPatch(self):
    new_patch = Patch()
    new_patch.index = len(self.patches)
    self.patches.append(new_patch)
    return new_patch

  def save(self, file_name):
    f = open(file_name, "w")
    for patch in self.patches:
      f.write(patch.toString() + "\n")
    f.write("0\n")
    f.close

  def toString(self):
    txt = ""
    for patch in self.patches:
      txt += patch.toString()
      txt += "\n"
    txt += "0\n"
    return txt

  def update(self):
    for i in range(0, len(self.patches)):
      self.patches[i].updateResolution()

    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighI1) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighI1:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateI1(hn)
      if len(self.patches[i].neighI2) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighI2:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateI2(hn)
      if len(self.patches[i].neighJ1) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighJ1:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateJ1(hn)
      if len(self.patches[i].neighJ2) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighJ2:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateJ2(hn)
      if len(self.patches[i].neighK1) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighK1:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateK1(hn)
      if len(self.patches[i].neighK2) > 0:
        hn = 0.0
        for neigh in self.patches[i].neighK2:
          hn = max(hn, self.patches[neigh].hi)
        self.patches[i].inflateK2(hn)

    for i in range(0, len(self.patches)):
      self.patches[i].updateResolution()

  def numCoreCells(self):
    N = 0
    for patch in self.patches:
      N += patch.numCells() - patch.numOverlapCells()
    return N

  def numOverlapCells(self):
    N = 0
    for patch in self.patches:
      N += patch.numOverlapCells()
    return N

  def minCoreSize(self):
    N = self.numCoreCells()
    for patch in self.patches:
      N = min(N, patch.numCells() - patch.numOverlapCells())
    return N

  def averageCoreSize(self):
    N = self.numCoreCells()
    for patch in self.patches:
      N = N + patch.numCells() - patch.numOverlapCells()
    return N/len(self.patches)

  def maxCoreSize(self):
    N = 0
    for patch in self.patches:
      N = max(N, patch.numCells() - patch.numOverlapCells())
    return N

  def numCells(self):
    return self.numCoreCells() + self.numOverlapCells()

  def setMinDim(self, N_min):
    for i in range(0, len(self.patches)):
      self.patches[i].hi = min(self.patches[i].hi, self.patches[i].Li()/N_min)
      self.patches[i].hj = min(self.patches[i].hj, self.patches[i].Lj()/N_min)
      self.patches[i].hk = min(self.patches[i].hk, self.patches[i].Lk()/N_min)

  def setOverlapFactor(self, of):
    for i in range(0, len(self.patches)):
      self.patches[i].overlap_factor = of

  def setNeighbourOverlapFactor(self, of):
    for i in range(0, len(self.patches)):
      self.patches[i].neigh_overlap_factor = of

  def printInfo(self):
    print "\n"
    print "number of patches       : %d" % len(self.patches)
    print "total number of cells   : %.2f*10^6" % (1e-6*self.numCells())
    print "                        : %d" % self.numCells()
    print "number of core cells    : %.2f*10^6" % (1e-6*self.numCoreCells())
    print "                        : %d" % self.numCoreCells()
    print "number of overlap cells : %.2f*10^6" % (1e-6*self.numOverlapCells())
    print "                        : %d" % self.numOverlapCells()
    print "overlap ratio in %%      : %.2f" % (100*float(self.numOverlapCells())/float(self.numCoreCells()))
    print "smallest core size      : %d" % self.minCoreSize()
    print "largest core size       : %d" % self.maxCoreSize()
    print "average core size       : %d" % self.averageCoreSize()

  def setMaxRes(self, max_h):
    for i in range(0, len(self.patches)):
      self.patches[i].hi = min(max_h, self.patches[i].hi)
      self.patches[i].hj = min(max_h, self.patches[i].hj)
      self.patches[i].hk = min(max_h, self.patches[i].hk)

  def setGrading1(self, growth_factor):
    done = False
    while not done:
      done = True
      for i in range(0, len(self.patches)):
        h = self.patches[i].hi
        h = min(h, self.patches[i].hj)
        h = min(h, self.patches[i].hk)
        h = h*growth_factor
        neighbours  = []
        neighbours += self.patches[i].neighI1
        neighbours += self.patches[i].neighI2
        neighbours += self.patches[i].neighJ1
        neighbours += self.patches[i].neighJ2
        neighbours += self.patches[i].neighK1
        neighbours += self.patches[i].neighK2
        for j in neighbours:
          if h < self.patches[j].hi:
            self.patches[j].hi = h
            done = False
          if h < self.patches[j].hj:
            self.patches[j].hj = h
            done = False
          if h < self.patches[j].hk:
            self.patches[j].hk = h
            done = False

  def setGrading2(self, growth_factor):
    done = False
    while not done:
      done = True
      for i in range(len(self.patches)):
        h = self.patches[i].hi
        h = min(h, self.patches[i].hj)
        h = min(h, self.patches[i].hk)
        h = h*growth_factor
        neighbours1 = set([])
        neighbours1 = neighbours1.union(self.patches[i].neighI1)
        neighbours1 = neighbours1.union(self.patches[i].neighI2)
        neighbours1 = neighbours1.union(self.patches[i].neighJ1)
        neighbours1 = neighbours1.union(self.patches[i].neighJ2)
        neighbours1 = neighbours1.union(self.patches[i].neighK1)
        neighbours1 = neighbours1.union(self.patches[i].neighK2)

        neighbours2 = set([])
        for j in neighbours1:
          neighbours2 = neighbours2.union(self.patches[j].neighI1)
          neighbours2 = neighbours2.union(self.patches[j].neighI2)
          neighbours2 = neighbours2.union(self.patches[j].neighJ1)
          neighbours2 = neighbours2.union(self.patches[j].neighJ2)
          neighbours2 = neighbours2.union(self.patches[j].neighK1)
          neighbours2 = neighbours2.union(self.patches[j].neighK2)

        neighbours2 = neighbours2.difference(neighbours1)

        neighbours = neighbours1
        for j in neighbours2:
          neighbours3 = set([])
          neighbours3 = neighbours3.union(self.patches[j].neighI1.intersection(neighbours1))
          neighbours3 = neighbours3.union(self.patches[j].neighI2.intersection(neighbours1))
          neighbours3 = neighbours3.union(self.patches[j].neighJ1.intersection(neighbours1))
          neighbours3 = neighbours3.union(self.patches[j].neighJ2.intersection(neighbours1))
          neighbours3 = neighbours3.union(self.patches[j].neighK1.intersection(neighbours1))
          neighbours3 = neighbours3.union(self.patches[j].neighK2.intersection(neighbours1))
          if len(neighbours3) > 1:
            neighbours.add(j)

        for j in neighbours:
          if h < self.patches[j].hi:
            self.patches[j].hi = h
            done = False
          if h < self.patches[j].hj:
            self.patches[j].hj = h
            done = False
          if h < self.patches[j].hk:
            self.patches[j].hk = h
            done = False

  def createRectGrid(self, x, y, z, h_snap = -1):
    if h_snap > 0:
      for i in range(len(x)):
        x[i] = int(x[i]/h_snap)*h_snap
      for i in range(len(y)):
        y[i] = int(y[i]/h_snap)*h_snap
      for i in range(len(y)):
        z[i] = int(z[i]/h_snap)*h_snap
    patch = []
    for i in range(len(x) - 1):
      patch.append([])
      for j in range(len(y) - 1):
        patch[i].append([])
        for k in range(len(z) - 1):
          patch[i][j].append(self.createPatch())
          patch[i][j][k].setName(str(i) + ',' + str(j) + ',' + str(k))
          patch[i][j][k].setBox(x[i], y[j], z[k], x[i+1], y[j+1], z[k+1])
          if i > 0:
            patch[i][j][k].addNeighI1(patch[i-1][j][k])
          if j > 0:
            patch[i][j][k].addNeighJ1(patch[i][j-1][k])
          if k > 0:
            patch[i][j][k].addNeighK1(patch[i][j][k-1])
    return patch

  def setI1Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighI1) == 0:
        self.patches[i].setI1Flux(flux)

  def setI2Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighI2) == 0:
        self.patches[i].setI2Flux(flux)

  def setJ1Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighJ1) == 0:
        self.patches[i].setJ1Flux(flux)

  def setJ2Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighJ2) == 0:
        self.patches[i].setJ2Flux(flux)

  def setK1Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighK1) == 0:
        self.patches[i].setK1Flux(flux)

  def setK2Flux(self, flux):
    for i in range(0, len(self.patches)):
      if len(self.patches[i].neighK2) == 0:
        self.patches[i].setK2Flux(flux)

  def setOverlap(self, overlap):
    for i in range(0, len(self.patches)):
      self.patches[i].setOverlap(overlap)

  def setAllFluxes(self, flux):
    self.setI1Flux(flux)
    self.setI2Flux(flux)
    self.setJ1Flux(flux)
    self.setJ2Flux(flux)
    self.setK1Flux(flux)
    self.setK2Flux(flux)

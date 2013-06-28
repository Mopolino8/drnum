#!/usr/bin/python

class Patch:
  
  def __init__(self):
    self.x1 = 0
    self.y1 = 0
    self.z1 = 0
    self.x2 = 0
    self.y2 = 0
    self.z2 = 0
    self.Ni = 0
    self.Nj = 0
    self.Nk = 0
    self.hi = 0
    self.hj = 0
    self.hk = 0
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
    self.overlap = 2
    
  def setNeighI1(self, neigh):
    self.neighI1.add(neigh.index)
    self.SeekLayersI1 = self.overlap
    neigh.neighI2.add(self.index)
    neigh.SeekLayersI2 = neigh.overlap

  def setNeighI2(self, neigh):
    self.neighI2.add(neigh.index)
    self.SeekLayersI2 = self.overlap
    neigh.neighI1.add(self.index)
    neigh.SeekLayersI1 = neigh.overlap
    
  def setNeighJ1(self, neigh):
    self.neighJ1.add(neigh.index)
    self.SeekLayersJ1 = self.overlap
    neigh.neighJ2.add(self.index)
    neigh.SeekLayersJ2 = neigh.overlap

  def setNeighJ2(self, neigh):
    self.neighJ2.add(neigh.index)
    self.SeekLayersJ2 = self.overlap
    neigh.neighJ1.add(self.index)
    neigh.SeekLayersJ1 = neigh.overlap

  def setNeighK1(self, neigh):
    self.neighK1.add(neigh.index)
    self.SeekLayersK1 = self.overlap
    neigh.neighK2.add(self.index)
    neigh.SeekLayersK2 = neigh.overlap

  def setNeighK2(self, neigh):
    self.neighK2.add(neigh.index)
    self.SeekLayersK2 = self.overlap
    neigh.neighK1.add(self.index)
    neigh.SeekLayersK1 = neigh.overlap
    
  def setTopNeigh(self, neigh):
    self.setNeighK2(neigh)
    
  def setBottomNeigh(self, neigh):
    self.setNeighK1(neigh)
    
  def setLeftNeigh(self, neigh):
    self.setNeighI1(neigh)
      
  def setRightNeigh(self, neigh):
    self.setNeighI2(neigh)
    
  def setFrontNeigh(self, neigh):
    self.setNeighJ1(neigh)
    
  def setBackNeigh(self, neigh):
    self.setNeighJ2(neigh)
    
  def setRes(self, h):
    self.hi = h
    self.hj = h
    self.hk = h
    
  def setResI(self, h):
    self.hi = h

  def setResJ(self, h):
    self.hj = h

  def setResK(self, h):
    self.hk = h

  def setBox(self, x1, y1, z1, x2, y2, z2):
    self.x1 = x1
    self.y1 = y1
    self.z1 = z1
    self.x2 = x2
    self.y2 = y2
    self.z2 = z2
    
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
    
  def inflateI1(self, hn):
    self.x1 -= 0.5*hn + 2*self.hi
    
  def inflateI2(self, hn):
    self.x2 += 0.5*hn + 2*self.hi
    
  def inflateJ1(self, hn):
    self.y1 -= 0.5*hn + 2*self.hj
    
  def inflateJ2(self, hn):
    self.y2 += 0.5*hn + 2*self.hj
    
  def inflateK1(self, hn):
    self.z1 -= 0.5*hn + 2*self.hk
    
  def inflateK2(self, hn):
    self.z2 += 0.5*hn + 2*self.hk
    
  def toString(self):
    txt = "1001\n{\n"
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
      txt += " far"
    else:
      txt += " 0"
    if len(self.neighI2) == 0:
      txt += " far"
    else:
      txt += " 0"
    if len(self.neighJ1) == 0:
      txt += " far"
    else:
      txt += " 0"
    if len(self.neighJ2) == 0:
      txt += " far"
    else:
      txt += " 0"
    if len(self.neighK1) == 0:
      txt += " far"
    else:
      txt += " 0"
    if len(self.neighK2) == 0:
      txt += " far"
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

    
class Mesh:
  
  def __init__(self):
    self.patches = []
    
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
      
  def numCells(self):
    return self.numCoreCells() + self.numOverlapCells()
    
  def printInfo(self):
    print "\n"
    print "number of patches       : %d" % len(self.patches)
    print "total number of cells   : %d" % self.numCells()
    print "number of core cells    : %d" % self.numCoreCells()
    print "number of overlap cells : %d" % self.numOverlapCells()
    print "overlap ratio in %%      : %.2f" % (100*float(self.numOverlapCells())/float(self.numCoreCells()))




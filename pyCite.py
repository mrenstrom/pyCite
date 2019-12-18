
import os
import gzip
from collections import defaultdict
from sys import argv
import Levenshtein as lv
#
# constants
#
antibodies = {
'GCATTCTGTCACCTA':'CD47'   ,
'CTTTGTACGAGCAAA':'CD52'   ,
'TCCTTTCCTGATAGG':'CD56'   ,
'TCCCTTGCGATTTAC':'CD45'   ,
'TATCCCTTGGGATGG':'CD3'    ,
'CTGGGCAATTACTCG':'CD19'   ,
'CAATCAGACCTATGA':'CD14'   ,
'TAACTCAGGGCCTAT':'CD33'   ,
'GCAGAAATCTCCCTT':'CD34'   ,
'CAGCCCGATTAAGGT':'Beta2 microglobulin',
'GCATTGTACGATTCA':'CD90'   ,
'AGACTAATAGCTGAC':'CD117'  ,
'CAGCCATTCATTAGG':'CD10'   ,
'TCAATCCTTCCGCTT':'CD45RA' ,
'CTTCACTCTGTCAGG':'CD123'  ,
'TGGATTCCCGGACTT':'CD7'    ,
'GTTTCCTTGACCAAG':'CD201 EPCR',
'TTCCGAGGATGATCT':'CD49f'  ,
'TGTTCCCGCTCAACT':'CD4'    ,
'TGGCTTCAGGTCCTA':'CD44'   ,
'GCTGCGCTTTCCATT':'CD8a'   ,
'TCTCAGACCTCCGTA':'CD14'   ,
'AAGTTCACTCTTTGC':'CD16'   ,
'CTCCGAATCATGTTG':'CD45RO' ,
'TTCTGGGTCCCTAGA':'CD20'   ,
'TACGCCTATAACTTG':'CD11c'  ,
'TGGTAACGACAGTCC':'CD133'  ,
'CATTAACGGGATGCC':'CD5'    ,
'GTCTCTTGGCTTAAA':'CD69'   ,
'AATAGCGAGCAAGTA':'HLA-DR' ,
'GACAAGTGATCTGCA':'CD11b'  ,
'AAGTATGCCCTACGA':'CD64'   ,
'CATTAACGGGATGCC':'CD235ab',
'TTTCAACGCCCTTTC':'CD13'   ,
'TGTACCCGCTTGTGA':'CD38'   ,
'TGCAATTACCCGGAT':'CD45'   ,
'TCACCAGTACCTAGT':'CD15'   ,
'CCGTGTTCCTCATTA':'CD71'}

hashtags = {
'GTCAACTCTTTAGCG':'Hashtag1',
'TGATGGCCTATTGGG':'Hashtag2',
'TTCCGCCTCTCTTTG':'Hashtag3'}


class Cell:
  def __init__(self,bc):
    self._adtDict = {}
    self._bc = bc

  def addADT(self,adt,umi):
    try:
      # this adt already added
      umiDict = self._adtDict[adt]
      umiDict[umi] += 1
    except KeyError:
      # a new adt
      umiDict = defaultdict(int)
      umiDict[umi] += 1
      self._adtDict[adt] = umiDict
  def numUMI(self):
    count = 0
    for adt,umiDict in self._adtDict.items():
      count += len(umiDict)
    return(count)
  #
  #
  #
  def numHit(self):
    count = 0
    for adt,umiDict in self._adtDict.items():
      for k,v in umiDict.items():
        count += v
    return(count)
  #
  #
  #
  def merge(self,cellObj):
    for adt,umiD in cellObj._adtDict.items():
      try:
        # this obj already has adt
        targetUMID = self._adtDict[adt]
        for k,v in umiD.items():
          targetUMID[k] += v
      # add num dict for this new adt
      except KeyError:
        targetUMID = defaultdict(int)
        for k,v in umiD.items():
          targetUMID[k] += v
        self._adtDict[adt] = targetUMID
  # 
  # number of unqiue UMIs for adt
  #
  def getADTCount(self,adtSeq):
    if adtSeq in self._adtDict.keys():
      return(len(self._adtDict[adtSeq]))
    return(0)
  #
  #
  #
  def errorCorrectUMI(self):
    for ab,umiDict in self._adtDict.items():
      #sorted small to large
      sd = sorted(umiDict.items(),key=lambda x:umiSortFunc(x[0],x[1]))
      for i in range(0,len(sd)):
        um_i = sd[i][0]
        c_i  = umiDict.get(um_i,0)
        for j in range(len(sd)-1,0,-1):
          if i != j:
           um_j = sd[j][0]
           c_j  = umiDict.get(um_j,0)
           # has entry already been merged
           if c_j == 0: continue
           d = lv.distance(um_i,um_j)
           if d <= 2:
             #print("merge")
             #print("  {}   {}  ->   {}   {}".format(um_i,c_i,um_j,c_j))
             #print("  confirn count i = {}    count j = {}".format(umiDict[um_i],umiDict[um_j]))
             umiDict[um_j] += umiDict[um_i]
             #print("  final count = {}".format(umiDict[um_j]))
             del umiDict[um_i]
             break
  #
  # remove all items from adt dict - makes object empyt
  #
  def clear(self):
    self._adtDict = {}

#----------------------------------------------------------------------------------------------------
#
# simulataneously read the R1 and R2 files for ADT/HTO
# R1 contains 16 bit cell barcode then 10 bit UMI
# R2 contains 15 bit ADT code
#
#
# the rd dictionary contains a key for each cell barcode and a UMI dictionary for each unique UMI
#
#----------------------------------------------------------------------------------------------------
def qCheck(seq,q,rd):
  for i,ch in enumerate(q):
    qScore = ord(ch) - 33
    if qScore < 20:
       seq = seq[0:i] + 'N' + seq[i+1:]
       
  badBases = 0
  for ch in seq:
    if ch == 'N': badBases += 1
 
  #if badBases > 3:
  #  print("read {}  seq {}".format(rd,q))
  #  print( [ord(n)-33 for n in q])
      

  return((badBases,seq))
#----------------------------------------------------------------------------------------------------
#
# errorCorrectADT
#
# look for exact match of adt code in antibody list
#
# if not found look for close hamming match
#
#----------------------------------------------------------------------------------------------------
def errorCorrectADT(adt):
  #
  # exact match
  #
  if adt in adtEntries.keys():
    return(adt)
  #
  # error-correct
  #
  for ab in adtEntries.keys():
    d = lv.distance(ab,adt)
    if d <= 3:
      return(ab)
  return("")

#----------------------------------------------------------------------------------------------------
#
# readInputFiles
#
# read R1 and R2 sequence files
#   find cell barcode and UMI from R1
#   find antibode code from R2
#   add to cell barcode dict  [barcode]
#     each cell barcode entry contains adt dict of UMIs
#
#----------------------------------------------------------------------------------------------------
def readInputFiles(readDir,dataFiles):
  count   = 0
  qFail   = 0
  adtFail = 0
  rd = {}
  for df in dataFiles:
    print("Process {}".format(df))
    dataFile1 = readDir + df + "_R1_001.fastq.gz"
    dataFile2 = readDir + df + "_R2_001.fastq.gz"

    with gzip.open(dataFile1,mode='rt') as fh1:
      with gzip.open (dataFile2,mode='rt') as fh2:
        while 1:  
          count += 1
          id1= fh1.readline().strip()
          if id1 == '':break
          id1 = id1.split(' ')[0]
          seq1 = fh1.readline().strip()
          e1 = fh1.readline().strip()
          q1 = fh1.readline().strip()
          
          id2 = fh2.readline().strip()
          if id2 == '':break
          id2 = id2.split(' ')[0]
          seq2 = fh2.readline().strip()[0:15]
          e2 = fh2.readline().strip()
          q2 = fh2.readline().strip()[0:15]
          #
          # Q check
          #
          bad1,seq1 = qCheck(seq1,q1,1)
          bad2,seq2 = qCheck(seq2,q2,2)
  
          if (bad1 > 2) or (bad2 > 3):
            qFail += 1
            continue
  
          bc = seq1[:16]
          umi = seq1[16:26]
          adt = errorCorrectADT(seq2)
          if adt == "":
            # counld not error-correct
            adtFail += 1
            continue
          #
          # read 1 and read 2 IDs must match
          # 
          if id1 != id2:
              print("Fatal Error, R1 and R2 sequences not aligned {} {}".format(id1,id2))
              quit(-1)
          #
          # add to dict
          #
          try:
              adtObj = rd[bc]
              adtObj.addADT(adt,umi)
          except KeyError:
              adtObj = Cell(bc)
              rd[bc] = adtObj
              adtObj.addADT(adt,umi)
  
          if count % 100000 == 0:
              print(count)
    print("Total lines read = {}".format(count))

  return(rd,qFail,adtFail)
#
# buildKMER
#
# build sequence-matching structure
#
def buildKMER(rd):
  print("build kmer struct for {} barcodes".format(len(rd)))
  kmer = [{},{},{},{}]
  #
  #
  #
  for bc,adtObj in rd.items():
    if adtObj.numUMI() > 1:
      #print(bc)
      kl = []
      for i in range(4):
          kl.append(bc[4*i:4*i+4])
        
      for i,k in enumerate(kl):
          #print(i,k)
          try:
              kmer[i][k].append(bc)
          except KeyError:
              kmer[i][k] = [bc]
      #print('len kmer = {}'.format(len(kmer)))
    
  print('KMER Done')
  return(kmer)
#
# errorCorrectCellBarcodes
#
# sd - sorted list of cell barcode, adt dict
#    - sorted by how many umi found
#
def errorCorrectCellBarcodes(sd,rd):
  print("Error-correct {} cell barcodes".format(len(sd)))
  bcCount = 0
  match = 0
  for bc,cellObj in sd:
    numUMI = cellObj.numUMI()
    if numUMI > 500:
      print("break for numUMI = {}".format(numUMI))
      break
    #
    # cell barcodes are 16 bases, break down into 4 groups of 4 bases
    #
    kl = []
    for i in range(0,4):
        kl.append(bc[4*i:4*i+4])
    #
    # look in kmer lib for each group of 4 bases at the proper location, gather list
    # of barcodes that match
    #
    rl = defaultdict(int)
    for i,k in enumerate(kl):
        kmerLib = kmer[i]
        if k in kmerLib.keys():
            bcList = kmerLib[k]
            for cellBC in bcList:
                rl[cellBC] += 1
    #print("length of mer match list = {}".format(len(rl)))
    #
    # filter for at least 2 kmer matches and numUMI > 0
    #
    matL = []
    for mBC,merMatchCount in rl.items():
      if merMatchCount >=2:
        
        numUMI = rd[mBC].numUMI()

        #print(merMatchCount,numUMI)
        if numUMI > 0:
          matL.append((mBC,numUMI))
    #
    # for the filtered list see if any within hamming 2
    #
    #print("Length numADT list = {}".format(len(matL)))

    if len(matL) > 0:
      #
      # sort targets based on num UMI
      # 
      rs = sorted(matL,key =lambda x:x[1],reverse=True)
      #
      # see if match really within hamming 2
      #
      for bestBC,numUMI in rs:
        d =  lv.distance(bc,bestBC)
        #print("{}  ->  {}  numUMI = {}, d = {}".format(bc,bestBC,numUMI,d))
        if (d > 0) and (d <= 2):
            #print("{}  {} mer match count = {}, levD = {}".format(bc,bestBC,numUMI,d))
            match += 1
            # add to good...it might match master list and not be in 
            targetCellObj = rd[bestBC]
            targetCellObj.merge(cellObj)

            #
            # clear merged bacrode-cellObj
            #
            rd[bc].clear()
            break

    bcCount += 1
    if bcCount % 1000 == 0:
        print("loop count = {}, match = {}".format(bcCount,match))
  # 
  # clear deleted sequences
  #
  retD = {}
  for cellBC,cellObj in rd.items():
    if  cellObj.numUMI() > 10:
      retD[cellBC] = cellObj


  return(retD)
#
# errorCorrectUMI
#
# results stored in barcode dict  [barcode] -> adtDICT
#
# adtDict is a umiDict for each antibody   [adt] -> umiDICT
#
# umiDict for each antibody has seq,count
#
# each umi is good for one final count for antibody abundance, no matter how
# many are detected. UMIs may have erorrs so merge UMI withing 2
#
def umiSortFunc(umiK,umiV):
  bad = len([x for x in umiK if x == 'N'])
  if bad > 0:
    return(-bad)
  return(umiV)

def errorCorrectUMI(rd):
  for bc,cellObj in rd.items():
    cellObj.errorCorrectUMI()


#-----------------------------------------------------------------------------------------------
#
#
# main script
#
#
#-----------------------------------------------------------------------------------------------

adtType,readDir = argv[1:3]
readFiles = argv[3:]
readFile = readFiles[0] 
#
# must be adt or hto
#
if adtType == "adt":
  adtEntries = antibodies
elif adtType == "hto":
  adtEntries = hashtags
else:
  print("Error, invalid adt type: {}".format(adtType))
  print("must be adt or hto")
  quit(-1)
#
#
# read in R1 and R@ source files
#  R1 = cell barcode + UMI
#  R2 = ADT or HTO   -  antibody code
#
print("Read fastq data files")
print(readFiles)

rd,qFail,adtFail = readInputFiles(readDir,readFiles)

print("Read {} unique cell barcodes".format(len(rd)))
print("Q1 fail           = {}".format(qFail))
print("ADT sequence fail = {}".format(adtFail))

#
#
# logging
#
with open(readFile + "_pre_hitcount.txt",'w') as fh:
  ld = sorted(rd.items(), key=lambda x:x[1].numHit(),reverse=True)
  for bc,cellObj in ld:
    if cellObj.numHit() > 1:
      fh.write("{}\n".format(cellObj.numHit()))
#
#
# build kmer object for fast barcode error correction
#
#
kmer = buildKMER(rd)
#
# attempt to error correct cell barcodes
#
#   sort cell barcodes by umi count
#
# sd = list(barcode, umiList)
#
sd = sorted(rd.items(), key=lambda x:x[1].numUMI())
#
#
#
rd = errorCorrectCellBarcodes(sd,rd)
#
# for each cell barcode, error-correct included UMI for each ADT/HTO 
#
errorCorrectUMI(rd)
#
#
# write out debug info
#
#
with open(readFile+"_debug.txt",'w') as fh:
  sd = sorted(rd.items(),key = lambda x:x[1].numUMI(),reverse=True)
  for i,(bc,cellObj) in enumerate(sd):
    fh.write("{}\n".format(bc))
    for adt,umiDict in cellObj._adtDict.items():
      fh.write("  {}   {}     {}\n".format(adt,adtEntries[adt],len(umiDict)))
      for umi,count in umiDict.items():
        fh.write("    {}  {}\n".format(umi,count))
    if i > 10:break


#
#
# logging
#
with open(readFile + "_post_UMIcount.txt",'w') as fh:
  ld = sorted(rd.items(), key=lambda x:x[1].numUMI(),reverse=True)
  for bc,cellObj in ld:
    if cellObj.numUMI() > 1:
      fh.write("{}\n".format(cellObj.numUMI()))
#
#
# write cell barcode matrix
#
#
with open(readFile+"_CBM.txt",'w') as fh:
  #
  # prepare
  #
  sAB = sorted(adtEntries.items(),key=lambda x:x[1])
  # sort by num unique UMI
  sd = sorted(rd.items(),key = lambda x:x[1].numUMI(),reverse=True)
  #
  # header
  #
  fh.write(" \t")
  for bc,cellObj in sd:
    fh.write("{}\t".format(bc))
  fh.write("\n")
  #
  # data
  #
  for antiSeq,antiName in sAB:
    fh.write("{}\t".format(antiName))
    cellCount = 0
    for bc,cellObj in sd:
      try:
        count = cellObj.getADTCount(antiSeq)
      except KeyError:
        count = 0
      fh.write("{}\t".format(count))
      cellCount += 1
      if cellCount > 4000:break
    fh.write("\n")
    


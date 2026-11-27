# -*- coding: utf-8 -*-
""" imRNAquant/f_imRNAquant.py
imRNAquant v1.0.0 by Ozkurede (2026)
Functions used by imRNAquant.py.
u_qb
"""
import datetime, os, copy, gzip

def udate17 ():    
   return datetime.datetime.now().strftime ("(%Y%m%d_%H%M%S)")
def udate19 ():    
   return datetime.datetime.now().strftime ("%Y %m %d %H:%M.%S")

def p (*x_, n = None, sep = "", end = ""): 
   output = ""    
   for x in x_:
      p = str (x)
      output = output + p + sep
   # remove the last separator if exists.
   output = output [0 : len (output) - len (sep)]
   print (output, end = end, flush = True)

def pl (*x_, n = None, sep = "", end = "\n"): # p | end = "\n"   
   p (*x_, n = n, sep = sep, end = end)

def pll (*x_, n = None, sep = "", end = "\n"): # p | pre = "\n", end = "\n"
   pl ()
   p (*x_, n = n, sep = sep, end = end)
    
def ps (*x_, n = None, sep = " ", end = ""): # p | sep = " ". 
   p (*x_, n = n, sep = sep, end = end)

def dec (value, n = 3): # show n digits after decimal point
   return round (value, n)

class pwClass:
   
   def __init__ (self):
      self.out = ''
   
   def makeString (self, *s_):
      string = ""
      for s in s_: 
         string = string + str (s)
      return string
   
   def prints (self, pstr):
      print (pstr, end = "") 
       
   def s (self, *x_):
      pstr = self.makeString (*x_) 
      self.prints (pstr)
      self.out += pstr

   def l (self, *x_):
      pstr = self.makeString (*x_) + "\n" 
      self.prints (pstr)
      self.out += pstr

   def ll (self, *x_):
      self.l ()
      self.l (*x_)

def comp (s):
   s = s.upper ()
   sComp = ""
   comp_table = {"A":"T", "T":"A", "G":"C", "C":"G","N":"N"}
   for nt in s:
      if nt in comp_table:
         ntComp = comp_table [nt]
      else: 
         ntComp = "X"; 
         print ("Warning: Unknown nt found and replaced with X")
      sComp = ntComp + sComp
   return sComp

def load_mtDseq (mtDseqFile):
   mtDseq = ""
   with open (mtDseqFile, "r") as file:
      for line in file:
         if line.startswith (">") == 0:
            mtDseq = mtDseq + line.rstrip ()
   return mtDseq

def read (seq, a, b):

   def modNoZeroIndex (a, l):
      if (a < 0):
         if (a < (0 - l)): # if more than 1 times passed the 0, rewind to the first pass.
             a = a % l - l 
         # passed 0 once, map to positive while skipping 0 index.
         a = a + l + 1
      if (a == 0): # if a is 0, set it to -1 (last nt on the period)
         a = a + l 
      if (a > l): # if higher than l, map it to positive. 
         a = a % l
      return (a)

   l = len (seq)
   a = modNoZeroIndex (a, l)
   b = modNoZeroIndex (b, l)

   if (a <= b):
      return (seq [a-1 : b])
   else:
      return (seq [a-1 : l] + seq [0 : b])

def procRoot (read_file, frag__, CL): # the processor; 
   
   read_process_limit = CL ['read_process_limit'] # D_qbh remove before release.
   
   # read_process_limit = 300000000 # process upto 30M reads per sample.

   root = read_file ["id"]
   pl ('›››---Function run started: procRoot (' + root + ')')
   pl ('›››---  Module, ', __name__, '; parent, ', os.getppid(),'; self, ', os.getpid ())
   n__ = {}
   n__ [root] = copy.deepcopy (frag__)
   Rfile__ = {}
   Rfile__ [root] = {}

   i_Rpair = 0
   for Rtype in ['reverse', 'forward']:
      if read_file [Rtype] == None: continue # uD_qb8 test this.

      i_Rpair += 1 # Number of Rfiles loaded so far for this root. if 2, perform calculations in this loop.
      fn = read_file [Rtype]

      R = {'fileId': root + '_' + Rtype,
           'fileName': fn, 
           'file': CL ['readFileFolder'] + "/" + fn, 
           'fileType': Rtype}

      pll ('Loading reads from ', R ['fileId'], ': ', R ['file'])

      with gzip.open (R ['file'], 'rt') as fileIn:
         content = [line.strip () for i, line in enumerate (fileIn) if i % 4 == 1 and i < 4 * read_process_limit]
      n_reads = len (content)

      pl ("  Loaded ", n_reads , " reads from ", R ['fileId'])
      pl ()
      
      R ['n_reads']  = n_reads
      Rfile__ [root] [Rtype] = copy.deepcopy (R)

      for fragId in n__[root]: 
         if 'count' not in n__ [root][fragId]: # Create in R1 loop; do not overwrite in R2 loop.
            n__ [root][fragId]['count'] = {}
         seqSelf = n__ [root][fragId]['seq']
         seqComp = n__ [root][fragId]['seqComp']

         if Rtype == 'reverse':
            countComp_R1 = sum (seqComp in line for line in content)
            count_label = fragId + '_comp.' + R ['fileId']
            n__ [root][fragId]['count']['comp.R1'] = countComp_R1 
            p (count_label, ", ", countComp_R1, ". ")
            del countComp_R1, count_label
         
         if Rtype == 'forward':
            countSelf_R2 = sum (seqSelf in line for line in content)
            count_label = fragId + '_self.' + R ['fileId']
            n__ [root][fragId]['count']['self.R2'] = countSelf_R2
            p (count_label, ", ", countSelf_R2, ", ")
            del countSelf_R2, count_label

         if read_file ['reverse'] != None and read_file ['forward'] != None:
            if i_Rpair == 2:
               pl ("Paired-end reads, this is the second one. Calculating countLeqv.")
               n__ [root][fragId]['count']['Leqv'] = n__ [root][fragId]['count']['comp.R1'] + n__ [root][fragId]['count']['self.R2']
            else: 
               pl ("Paired-end reads; this is the first one.")

         if read_file ['reverse'] == None or read_file ['forward'] == None:
            if i_Rpair == 1:
               pl ("Single-end reads. Assigning countLeqv.")
               count_key = 'self.R2' if Rtype == 'forward' else 'comp.R1'
               n__ [root][fragId]['count']['Leqv'] = n__[root][fragId]['count'][count_key]
               
            else:
               pl ('ERROR: Not paired-end but i_Rpair is not 1')

   pl ("End of ", R ['fileId'], " count.")
   return root, Rfile__ [root], n__ [root]

def titleRow (read_files_):
   title_RfileList =''
   for read_file in read_files_: 
      title_RfileList += read_file ['id'] + ','
   title = 'siteId,fragType,fragSubType,relToRef,' + title_RfileList + 'contAvg,,siteId,' + title_RfileList
   return title

def calc_contMean (count__, cont_, absolute = None):
   if len (count__) == 0: return 'noData'
   if len (cont_) == 0: return 'noCont'
   contSum = 0
   for R, value in count__.items():
      if R in cont_: 
         # if any control value is below 10, return -1 (noncalculated).
         if value == 'Ref0' or value < 10: 
            
            if absolute == 'absolute': return -1
         if value == 'Ref0':
            contSum += 0
         else:
            contSum += value
   contMean = dec (contSum / len (cont_))
   return contMean

def calc_foldContMean (count__, contMean):
   foldContMean__ = {}
   for sampId, value in count__.items():
      if contMean == 'Unavail': foldContMean__[sampId] = 'Unavail'
      elif contMean == 0: foldContMean__ [sampId] = 'Avg0'
      elif value == 'Ref0': foldContMean__ [sampId] = 'vRef0'
      elif type (value) == str: foldContMean__ [sampId] = 'str' #value is other string
      elif contMean == -1: foldContMean__ [sampId] = '!LowCt'
      else: foldContMean__ [sampId] = dec (value / contMean)
   return foldContMean__

# End.
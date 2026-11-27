# -*- coding: utf-8 -*-
""" imRNAquant/imRNAquant.py
imRNAquant v1.0.0 by Ozkurede (2026)
u_qb
―――
Function: 
   Quantify mitochondrial imRNAs using parallel processing on multiple cores.

Input:
   mtDNA sequence file
      Selected in CONFIG.py
      Multiple reference genome sequences are provided under imRNAquant/resources/
   mtDNA map file (gene annotations)
      Selected in CONFIG.py
      Multiple map files, based on the sequence file, are provided under imRNAquant/resources/
   mtDNA UJI file (uncleaved junction identifier definitions)
      Selected in CONFIG.py
      Multiple UJI sets, based on the sequence file, are provided under imRNAquant/resources/
   fastq.gz read files 
      Folder and file names are set in CONFIG.py
      The files are provided by the user
Output:
   output/log.txt, process log
   output/report.csv, imRNA quantification report table
"""

import os, json, copy
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from f_imRNAquant import (udate17, udate19, dec, pl, pwClass, comp, load_mtDseq, read, 
                        procRoot, titleRow, calc_contMean, calc_foldContMean)
from CONFIG import load_config

pl ('RUNNING ', __name__, ' ', os.getpid ()) # informative; not needed.

if __name__ == '__main__':
   pw = pwClass()
   
   ### LOAD Configuration from config.py #######################################################
   CL = load_config ()

   #i Assign condId_ and testId_ based on "CL['samples_']" list in CONFIG.py
   contId_ = []; testId_ = []
   for read_file in CL ['samples_']:
      if read_file ['group'] == 'Control': contId_.append (read_file ['id'])
      if read_file ['group'] == 'Test': testId_.append (read_file ['id'])
 
   if len (contId_) < 3:
      raise SystemExit (
      "imRNAquant requires at least 3 control samples but " +
      str (len (contId_)) + " given." +
      " Modify CONFIG.py to have at least 3 control samples")
      
   CL ['jsidelenPost'] = 30
   CL ['jsidelen'] = 15
   refFrag_extension = {'refPre': 'tail', 'refPost': 'head'}

   # Output files
   log_file = "output/imRNAquant_" + udate17 () + "_log.txt"
   report_file = "output/imRNAquant_" + udate17 () + "_report.csv"

   # print main settings on screen and add to log
   pw = pwClass()
   pw.ll ("imRNAquant started, ", udate19 ())
   pw.l ("$mtDseq_file, ", CL ['mtDseq_file'] )
   pw.l ("$jsidelenPost, ", CL ['jsidelenPost'])
   pw.l ("$jrefinlen, ", CL ['jrefinlen'])
   pw.l ("$readFileFolder, ", CL ['readFileFolder'])
   pw.l ("$samples_, ", CL ['samples_'])
   pw.l ("$log_file, ", log_file)
   pw.ll ("$report_file, ", report_file)

   input ("Press ENTER to continue (Ctrl+C to stop) >")

   #i Load mtDseq
   pw.ll ("Loading mtDNA sequence from ", CL ['mtDseq_file'])
   mtDseq = load_mtDseq (CL ['mtDseq_file'])
   pw.l ("Loaded, ", len (mtDseq), " nts, $mtDseq")

   #i Load site__ data using the map (annotation) file.
   pw.ll ("Loading mtDmap from the file, ", CL ['mtDmap_file'])
   with open (CL ['mtDmap_file'], "r") as f: 
      site__ = json.load (f)
   pw.l ("Loaded, ", len (site__) - 2, " sites, $site__")

   #i Load UJIjson__ data using the UJI file.
   pw.ll ("Loading junctions from the file, ", CL ['mtD_UJI_file'])
   with open (CL ['mtD_UJI_file'], "r") as f: 
      UJIjson__ = json.load (f)
   pw.l ("Loaded, ", len (UJIjson__) - 2, " junctions, $UJIjson__")

   #i Add, to the frag__ dict, junctions and associated reference segments based on UJI definitions 
   pw.ll ('Generating frag__ dictionary (UJIs and reference segment sequences) using UJI definitions.')
   frag__ = {}

   for key, j in UJIjson__.items():
      if key [0] == '#': continue # skip noncontent keys
      seqstart = j ['UJIpos'] 
      seqend = j ['UJIpos'] + CL ['jsidelenPost'] - 1

      t_seq = read (mtDseq, seqstart, seqend)
      frag__ [key] = {
        'fragId': key,
        'fragType': 'junction',
        'jType': j ['type'],
        'jPos': j ['pos'],
        'pos': j ['pos'],
        'UJIpos': j ['UJIpos'],
        'seq': t_seq,
        'seqComp': comp (t_seq),
        'refPre': j ['refPre'],
        'refPost': j ['refPost']}

      pw.l ("Added junction, ", key )

      # for param, value in frag__ [key].items(): pw.l ('   ', param, ': ', value) # display detailed info
      #i Add reference segments for this UJI to the frag list. refPre → $refName + '_tail' ; refPost → $refName + '_head'

      frag__ [key]['refSet'] = {}

      ### BEGIN Add upsteam and downstream reference segments ###
      for refType in {'refPre', 'refPost'}:
         t_refTypeSub = refType + '_' + refFrag_extension [refType]
         t_ref = j [refType]
         t_refSubId = t_ref + '_' + refFrag_extension [refType]
         if refType == 'refPre': t_refSub_pos = site__ [t_ref]['end'] - CL ['jrefinlen']
         if refType == 'refPost': t_refSub_pos = site__ [t_ref]['start'] + CL ['jrefinlen']
         t_refSub_seq = read (mtDseq, t_refSub_pos - CL ['jsidelen'], t_refSub_pos + CL ['jsidelen'] - 1)

         #i Add to refSet (reference set) independent of its presence as a frag in frag__
         frag__ [key]['refSet'] [t_refTypeSub] = t_refSubId
         if t_refSubId not in frag__:
            frag__ [t_refSubId] = {
               'fragId': t_refSubId,
               'fragType': 'ref',
               'refSubType': refFrag_extension [refType],
               'refSubPos': t_refSub_pos,
               'pos': t_refSub_pos - CL ['jsidelen'], # seq strat pos
               'seq': t_refSub_seq,
               'seqComp': comp (t_refSub_seq),
               'refForh': key }

            pw.l ('  for ', key, ' Added ', t_refSubId )
            # for param, value in frag__[t_refSubId].items(): pw.l ('   ', param, ': ', value) # display detailed info
         else:
            pw.l ('  for ', key, ', tried adding ', t_refSubId, '; already added for ', frag__ [t_refSubId] ['refForh'])

         ### BEGIN Add refMid reference segment ###      
         #i Add, if not added already, refMid segments of refPre or refPost genes. refMid → $refName + '_mid'
         t_refMidId = t_ref + '_' + 'mid'
         if CL ['refMidFrag'] == True:
            t_refTypeMid = refType + '_mid'
            frag__ [key]['refSet'] [t_refTypeMid] = t_refMidId

            #i Add to refSet (reference set) independent of its presence as a frag in frag__
            frag__ [key]['refSet'] [t_refTypeSub] = t_refSubId
            if  t_refMidId not in frag__:
               t_refMid_pos = (site__ [t_ref] ["start"] + site__ [t_ref] ['end']) // 2
               t_refMid_seq =  read (mtDseq, t_refMid_pos - CL ['jsidelen'], t_refMid_pos + CL ['jsidelen'] - 1)
               frag__ [t_refMidId] = {
                   'fragId': t_refMidId,
                   'fragType': 'ref',
                   'refSubType': 'mid',
                   'refSubPos': t_refMid_pos,
                   'pos': t_refMid_pos - CL ['jsidelen'],
                   'seq': t_refMid_seq,
                   'seqComp': comp (t_refMid_seq),
                   'refForh': key }
               pw.l ('  for ', key, ' Added ', t_refMidId)
               # for param, value in frag__[t_refMidId].items(): pw.l ('   ', param, ': ', value) # display detailed info
            else:
               pw.l ('  for ', key, ', tried adding ', t_refMidId, '; already added for ', frag__ [t_refMidId] ['refForh'])
         ### END Add refMid reference segment ###
      ### END Add upsteam and downstream reference segments ###
   pw.l ("Added junctions and associated reference segments to frag__ dictionary.")
   pw.l ("Done: ", len (frag__), " items in the count list.")  
   pw.l ("A total of ", len (frag__), " segments will be searched in ", len (CL ['samples_']), " samples.")

### BEGIN Process read files and perform imRNA calculations ###################
   Rfile__ = {}
   n__ = {}

   ### BEGIN Parallel processing of read files ###
   with ProcessPoolExecutor (max_workers = CL ['MC_max']) as executor:
      # result_ = [
      #     procRoot(sample, frag__, CL)
      #     for sample in CL["samples_"]
      # ]

      result_ = list (executor.map (procRoot, CL ['samples_'], repeat (frag__), repeat (CL)))
   ### END Parallel processing of read files ###   

   ### BEGIN Assign results from parallel processes to variables ###  
   for root, Rfile_root, n_root in result_:
      Rfile__[root] = copy.deepcopy (Rfile_root)
      n__[root] = copy.deepcopy (n_root)
   ### END Assign results from parallel processes to variables ###
   pw.l ("Completed counting.", udate19 ())

   ### BEGIN Calculate relative counts ###
   pw.ll ("Calculating relative counts.")
   for R in n__:
      for fragId in n__ [R]:
         if n__ [R][fragId]['fragType'] != 'junction':
            continue
         n__ [R][fragId]['countRelTo'] = {}
         frag = n__ [R][fragId]

         for refSetLabel in frag ['refSet']:
            refSubId = frag ['refSet'] [refSetLabel]
            n__ [R] [fragId] ['countRelTo'] [refSetLabel] = {}
            n__ [R] [fragId] ['countRelTo'] [refSetLabel] ['refSubId'] = refSubId
            n__ [R] [fragId] ['countRelTo'] [refSetLabel] ['refSubCount'] = n__ [R][refSubId]['count']
            n__ [R] [fragId] ['countRelTo'] [refSetLabel] ['ratio'] = {}

            for countLabel in frag ['count']:
               if n__ [R] [refSubId] ['count'][countLabel] == 0: 
                  n__ [R] [fragId] ['countRelTo'] [refSetLabel] ['ratio'][countLabel] = "Ref0"
               else:
                  value = 1000 * frag ['count'][countLabel] / n__ [R][refSubId] ['count'][countLabel]
                  n__ [R][fragId] ['countRelTo'][refSetLabel] ['ratio'][countLabel] = dec (value)
   ### END Calculate relative counts ###

   pw.ll ("Calculated relative counts. ", udate19 (),"\n  Writing output to file: ", log_file)
   with open (log_file, "w", encoding = "utf-8-sig") as file: 
      file.write (pw.out)
   pw.l ('  Saved output file. ' + udate19 ())
### END Process read files and perform imRNA calculations #####################

### BEGIN Generate report file in comma-separated values (csv) format #########
   pw.out = '' # Reset pw.out, start generating output for the report file

   def table_fragCount (rep_fragType = None, rep_fragSubType = None):
      rep_fragType_text = rep_fragType if rep_fragType != None else ""
      rep_fragSubType_text = rep_fragSubType if rep_fragSubType != None else ""
      
      pw.ll ('Reporting absolute counts (heavy-strand transcripts). ', 
             rep_fragType_text, ' ', rep_fragSubType_text)
      
      # Print the title row
      pw.l (titleRow (CL ['samples_']))

      # Calculate and print data rows
      for fragId, frag in frag__.items ():
         # Skip unrequested items
         if rep_fragType != None and frag ['fragType'] != rep_fragType: continue
         if rep_fragSubType != None:
            fragSubKey = {'junction':'jType', 'ref':'refSubType'}[rep_fragType]
            if rep_fragSubType != frag [fragSubKey]: continue

         # Set fragSubKey key based on fragType string
         fragSubKey = {'junction':'jType', 'ref':'refSubType'} [frag ['fragType']]
         
         # Print frag info; keep printing data on this row.
         pw.s (fragId, ',', frag ['fragType'], ',' , frag [fragSubKey], ',none,')
         
         # Get count data for each sample and print
         count__ = {}
         for read_file in CL ['samples_']:
            R = read_file ['id']
            count__ [R] = n__ [R][fragId]['count']['Leqv']
            pw.s (count__ [R], ',')
         
         # Calculate and print contMean; then Print fragId; 
         contMean = calc_contMean (count__, contId_, 'absolute') 
         pw.s (contMean, ',,', fragId, ',')
         
         # Calculate and print foldContMean values for each sample
         foldContMean__ = calc_foldContMean (count__, contMean)
         for read_file in CL ['samples_']: 
            sId = read_file ['id']
            pw.s (foldContMean__ [sId], ',')
         
         pw.l () # move to the next row.

   def table_relCount (report_jType, report_refTypeSub):
      
      pw.ll ('Reporting relative count for jType ', report_jType, ' normalized to ', report_refTypeSub)
      pw.l (titleRow (CL ['samples_']))
      
      for fragId, frag in frag__.items ():
         if frag ['fragType'] != 'junction': continue
         if frag ['jType'] != report_jType: continue
         
         #i the cols are fragId, fragType, fragSubType; then values  
         pw.s (fragId, ',', frag ['fragType'], ',' , frag ['jType'], ',', frag ['refSet'] [report_refTypeSub], ',')

         # Load absolute counts, calc contMean to assess consistent detection
         count__ = {}
         for read_file in CL ['samples_']:
            sId = read_file ['id']
            count__ [sId] = n__ [sId][fragId]['count']['Leqv']
         contMean = calc_contMean (count__, contId_, 'absolute') 
         
         relCount__ = {}
         for read_file in CL ['samples_']:
            sId = read_file ['id']
            if contMean == -1: # Do not normalize if not consistently detected
               relCount__ [sId] = '!LowCt'
            else:
               relCount__ [sId] = n__[sId][fragId]['countRelTo'][report_refTypeSub]['ratio']['Leqv']         
            
            pw.s (relCount__ [sId], ',')
               

         # Do not calculate relContMean and foldContMean__ if controls include ref0 or contMean == -1 
         controls_include_ref0 = any (
            n__ [sId][fragId]['countRelTo'][report_refTypeSub]['ratio']['Leqv'] == 'Ref0' 
            for sId in contId_)
         
         if controls_include_ref0:
            relContMean = "Unavail"
         else:
            relContMean = -1 if contMean == -1 else calc_contMean (relCount__, contId_)
            
         foldContMean__ = calc_foldContMean (relCount__, relContMean)
         
         pw.s (relContMean, ',,', fragId, ',')
         for read_file in CL ['samples_']:
            sId = read_file ['id']
            pw.s (foldContMean__ [sId],',')
         pw.l ()

   # all frag counts, including UJIs and reference segments.
   table_fragCount ()
   # tR–mR UJI counts
   table_fragCount ('junction', 'tR-mR')
   # mR–tR UJI counts
   table_fragCount ('junction', 'mR-tR')
   
   # Calculate values normalized to mid-frags if this option is selected
   if CL ['refMidFrag'] == True:
      table_relCount ('tR-mR', 'refPost_mid')
      table_relCount ('mR-tR', 'refPre_mid')
      
   # Normalized counts for UJIs with 5' tRNA, relative to nearby segment in the mRNA gene.
   table_relCount ('tR-mR', 'refPost_head')
   # Normalized counts fo UJI with 3' tRNA, relative to nearby segment in the mRNA gene.
   table_relCount ('mR-tR', 'refPre_tail')

   # Optional. Display the report table during run on screen.
   print (pw.out.replace (',', '\t'))

   pl ("  Writing to file: ", report_file)
   with open (report_file, "w", encoding = "utf-8-sig") as file:
      file.write (pw.out)
   print ("  Saved output file. " + udate19 ())
### END Generate report file in comma-separated values (csv) format ###########
   print ("■ Run complete.")
# End.
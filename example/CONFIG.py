# -*- coding: utf-8 -*-
""" imRNAquant/example/CONFIG.py 
imRNAquant v1.0.0 by Ozkurede (2026)
u_qb
―――
Function: Configure imRNAquant settings based on toy example RNA-seq dataset
"""

### SET CONGIGURATION VALUES BELOW ###

def load_config ():
   
   CL = {} # Stores configuration values.
   
   ############################################################################
   # imRNA QUANTIFICATION SETTINGS
   ############################################################################
   
   # Maximum number of reads to process per sample.
   # The default is 50 million, which is more than most datasets provide.
   #   Increase this limit if you want to process more
   #   reads per sample, for example by setting it to 100 million or 1 billion.
   # Default, 50 million
   CL ['read_process_limit'] = 50000000 
   
   # Distance in nucleotides from the reference gene boundary to the reference
   #   segment. Range: 100–200 nt. 
   # No need to change unless mtDNA mutations are found on reference sites.
   # Default, 100 nt.
   CL ['jrefinlen'] = 100
   
   # If True, also report values normalized to segments from the middle of each
   #   reference gene. Use these values to check consistency with the results
   #   obtained using the standard reference segments.
   # Default: True.
   CL ['refMidFrag'] = True

   ############################################################################
   # PARALLEL PROCESSING
   ############################################################################

   # Maximum number of worker processes. 
   # Use no more workers than the number of CPU cores available to the program.
   # Set this value to 1 or 2 if memory is limited.
   # Default, 4
   CL ['MC_max'] = 4

   ############################################################################
   # REFERENCE files for mtDNA sequence, map, and UJIs
   ############################################################################
   
   # Use a matching set of sequence, annotation, and UJI files for the species
   #   and strain being analyzed.
   # After imRNAquant finishes, review the log to confirm the settings before
   #   interpreting the results.
   
   # mtDNA sequence file (FASTA)
   # Default value for mouse samples, "resources/mouse_mtDNA_seq.fa"
   # Default value for human samples, "resources/human_mtDNA_seq.fa"
   CL ['mtDseq_file'] = "resources/mouse_mtDNA_seq.fa"
   
   # mtDNA map (annotation) file (JSON)
   # Default value for mouse samples, "resources/mouse_mtDNA_map.json"
   # Default value for human samples, "resources/human_mtDNA_map.json"
   CL ['mtDmap_file'] = "resources/mouse_mtDNA_map.json" 
   
   # mtDNA UJI description file (JSON)
   # Default value for mouse samples, "resources/mouse_mtDNA_UJI.json"
   # Default value for human samples, "resources/human_mtDNA_UJI.json"
   CL ['mtD_UJI_file'] = "resources/mouse_mtDNA_UJI.json"
   
   
   ############################################################################
   # SAMPLES and associted READ FILES
   ############################################################################
   
   # Edit the RNAseq fastq.gz folder and sample list below to describe samples and associated files
   
   # The example below contains four control samples and four test samples
   #   with paired-end reads. It uses the toy fastq.gz files in example/fastq folder.
   
   
   # Absolute path to the folder containing the read files (.fastq.gz)
   # Example for local RNAseq files:
   # CL ['readFileFolder'] = 'c:/RNAseq/experimentName'
   CL ['readFileFolder'] = 'example/fastq'
   
   # Sample List
   # 
   # Add one entry per sample and assign each sample to "Control" or "Test".
   # Include at least three samples with group set to "Control".
   #
   # For SINGLE-END data, enter the fastq.gz filename under "forward" or "reverse"
   # according to the orientation of the reads relative to the RNA transcript.
   # Set the other field to None.
   #
   # Forward-oriented reads:
   #   "forward": "sample1_R1.fastq.gz", 
   #   "reverse": None
   #
   # Reverse-oriented reads:
   #   "forward": None, 
   #   "reverse": "sample1_R1.fastq.gz"
   #
   # An "_R1" or "_1" filename suffix identifies the read number.
   # Use your library preparation protocol to determine the read orientation.
   #
   # For PAIRED-END data in which read 1 is reverse-oriented, enter the read 1
   #   filename under "reverse" and the read 2 filename under "forward". Common
   #   suffixes are "_1" or "_R1" for read 1 and "_2" or "_R2" for read 2.
   #
   #
   # The example below contains four control samples and four test samples
   #   with paired-end reads.
   #
   CL ['samples_'] = [
      {  # sample 1 
         "id": "Sample1", # Sample name, no spaces.
         "group": "Control",
         "reverse": "s1_1.fastq.gz",
         "forward": "s1_2.fastq.gz",
      },
      {  # sample 2
         "id": "Sample2", # Sample name, no spaces.
         "group": "Control",   
         "reverse": "s2_1.fastq.gz",
         "forward": "s2_2.fastq.gz",
      },
      {  # sample 3 
         "id": "Sample3", # Sample name, no spaces.
         "group": "Control",   
         "reverse": "s3_1.fastq.gz",
         "forward": "s3_2.fastq.gz",
      },
      {  # sample 4
         "id": "Sample4", # Sample name, no spaces.
         "group": "Control",
         "reverse": "s4_1.fastq.gz",
         "forward": "s4_2.fastq.gz",
      },
      {  # sample 5
         "id": "Sample5", # Sample name, no spaces.
         "group": "Test",
         "reverse": "s5_1.fastq.gz",
         "forward": "s5_2.fastq.gz",
      },
      {  # sample 6
         "id": "Sample6", # Sample name, no spaces.
         "group": "Test",
         "reverse": "s6_1.fastq.gz",
         "forward": "s6_2.fastq.gz",
      },
      {  # sample 7
         "id": "Sample7", # Sample name, no spaces.
         "group": "Test",
         "reverse": "s7_1.fastq.gz",
         "forward": "s7_2.fastq.gz",
      },
      {  # sample 8
         "id": "Sample8", # Sample name, no spaces.
         "group": "Test",
         "reverse": "s8_1.fastq.gz",
         "forward": "s8_2.fastq.gz",
      },
   ]

  
   return CL
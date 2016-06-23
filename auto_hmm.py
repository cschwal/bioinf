#! /usr/bin/env python

'''
auto_hmm.py
  takes an list of GI numbers, downloads them, aligns them using MUSCLE,
  then builds an HMM using HMMER and optionally presses the HMM
  
  contact: jonathan.tietz@gmail.com
'''

# import modules
import csv, re, os, math, argparse, subprocess
from Bio import Entrez, Seq, SeqIO	# BioPython

# input options
parser = argparse.ArgumentParser(description='auto_hmm -- Generates HMM from a protein GI list')
parser.add_argument('-i', help='Input GI list', required=True)                        # input filename -- GI list
parser.add_argument('-e', help='Email address', required=True)                                          # email address is required for Entrez identification
parser.add_argument('-o', help='Output file base (no suffix)', default='autohmm_out') 
parser.add_argument('-c', help='Column containing GI numbers', type=int, default=0)
parser.add_argument('--verbose', help='Give verbose output (for debugging)', action='store_true', default=False) 
parser.add_argument('-p', help='Press HMM afterwards', action='store_true', default=False)
parser.add_argument('-n', help='HMM name')
args = parser.parse_args()

# initialize

Entrez.email = args.e
gi_list_file = args.i
fasta_file = args.o + '_sequences.fasta'
align_file = args.o + '_alignment.clw'
hmm_file = args.o + '_hmm.hmm'
press = args.p
v = args.verbose
if args.n:
  hmm_name = args.n
else:
  hmm_name = args.o
  
# functions

def csv_to_list(filename, column, verbose, comma_delimited):
  '''
                filename (str) - RODEO architecture csv file to open
                  column (int) - column of csv file containing protein GIs
                verbose (bool) - output errors/status
        comma_delimited (bool) - True to output as a comma-delimited list; False to return as standard list
  '''
  # open file
  try:
    f = open(filename, 'rU')
  except IOError:
    print "ERROR .. Could not read file in function rodeo_csv_to_list .. file = ", filename
    sys.exit()
  protein_list = []
  final_list = []
  with f:
    reader = csv.reader(f)
    for row in reader:
      protein_list.append(row[column].strip())	# just in case there's extra whitespace, strip()
  f.close()
  # condense list to unique values
  if verbose and len(protein_list) > len(set(protein_list)):
    print "Removing %d redundant members" % (len(protein_list) - len(set(protein_list)))
  protein_list = set(protein_list)
  # validate the format of each one
  for element in protein_list:
    match = re.search(r'^\d+$', element)
    if not match:
      if verbose: print "Removing bad protein GI: %s" % element
    else:
      final_list.append(element)
  # recapitulate list
  if verbose:
    print("List: ", final_list)
  if comma_delimited:
    return ",".join(final_list)
  else:
    return final_list

def gi_list_fetch_to_fasta(list_of_gis, fasta_filename):
  # post GI list to Entrez
  search_results = Entrez.read(Entrez.epost("protein", id=",".join(list_of_gis)))
  webenv = search_results["WebEnv"]
  query_key = search_results["QueryKey"]
  
  # perform EFetch on GI list and download results, saving to file
  fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", webenv=webenv, query_key=query_key)
  data = fetch_handle.read()
  fetch_handle.close()
  if v: print data
  with open(fasta_filename, 'w+') as f:
    f.write(data)

def program_exists(executable_name):
  successful = True
  try:
    subprocess.call(['which', executable_name])
  except:
    successful = False
  return successful

def align_proteins(fasta_filename, alignment_filename):
  os.system("muscle -in %s -out %s -clw" % (fasta_filename, alignment_filename))
  return True

def make_hmm(alignment_filename, hmm_filename, press_option, hmm_named):
  os.system("hmmbuild -n %s %s %s" % (hmm_named, hmm_filename, alignment_filename))
  if press_option:
    os.system("hmmpress %s" % hmm_filename)
  if v: print ' ... Pressed hmm'
  
# main

def main():

  # convert CSV to GI list
  if v: print ' ... Reading GI list',
  gi_list = csv_to_list(gi_list_file,args.c,v,False)
  if v: print ' - DONE!'

  # download proteins & save as a FASTA
  if v: print ' ... Fetching proteins via Entrez & saving as FASTA',
  gi_list_fetch_to_fasta(gi_list, fasta_file)
  if v: print ' - DONE!'
  
  # build alignment of proteins
  if v: print ' ... Building alignment with MUSCLE & saving in CLW format',
  did_alignment = False
  if program_exists('muscle'):
    did_alignment = align_proteins(fasta_file, align_file)
  else:
    print ' ... Critical error: MUSCLE not available'
  if v: print ' - DONE!'

  # build HMM from alignment
  if v: print ' ... Building HMM from alignment using HMMER',
  if not did_alignment:
    print ' ... Critical error: no alignment was generated'
  elif program_exists('hmmbuild'):
    make_hmm(align_file, hmm_file, press, hmm_name)
  else:
    print ' ... Critical error: HMMER not available'
  if v: print ' - DONE!'

  return

if __name__ == '__main__':
  main()
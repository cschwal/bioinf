#! /usr/bin/env python

'''
phmmerall.py
  takes an input FASTA file and performs phmmer exhaustively
'''

# import modules
import csv, re, os, math, argparse
from Bio import Entrez, Seq, SeqIO	# BioPython

# input options
parser = argparse.ArgumentParser(description='phmmerall -- Generates all v all phmmer for sequence similarity network generation')
parser.add_argument('-i', help='Input protein GI-containing CSV', required=True)                        # input filename -- either GI list or CSV containing GI list
parser.add_argument('-e', help='Email address', required=True)                                          # email address is required for Entrez identification
parser.add_argument('-o', help='Output file base (no suffix)', default='phmmerall_out') 
parser.add_argument('--verbose', help='Give verbose output (for debugging)', action='store_true') 
parser.add_argument('-c', help='Column containing GI numbers', type=int, default=0)
args = parser.parse_args()

# initialize

Entrez.email = args.e
arch_file = args.i
fasta_file = args.o + '.fasta'
hmmer_file = args.o + '.txt'
csv_file = args.o + '.csv'
v = args.verbose

# functions

def rodeo_csv_to_list(filename, column, verbose, comma_delimited):
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

def phmmerall(input_filename, output_filename):
  os.system("phmmer -o tmp.tmp --tblout %s --noali %s %s" % (output_filename, input_filename, input_filename))
  os.system("rm tmp.tmp")

def parse_phmmer(phmmer_filename, verbose):
  '''
      converts phmmer space-delimited file to list of tuples for alignment scores and E-values
  '''
  # open file
  try:
    f = open(phmmer_filename, 'rU')
  except IOError:
    print "ERROR .. Could not read file in function parse_phmmer .. file = ", phmmer_filename
    sys.exit()
  result_list = []
  interactions = []
  with f:
    reader = f.readlines()
    for row in reader:
      comments = re.match(r'^\#.*$',row)
      if not comments:
        values = re.match(r'^(\S+?)\s+?\S+?\s+?(\S+?)\s+?\S+?\s+?(\S+?)\s+?(\S+?)\s+?',row)
        if not values and verbose: # If didn't match expected regex, warn user
          print("Encountered bad pHMMER row ...",row)
        if not (values.group(1) == values.group(2)):    # only add these if it's not a self-reference
          if float(values.group(3)) == 0:
            log_eval = -10000
          else:
            log_eval = math.log(float(values.group(3)))
          quad = (values.group(1), values.group(2), values.group(3), log_eval, values.group(4))
          # make sure not already duplicated
          if (values.group(2), values.group(1)) not in interactions:
            result_list.append(quad)
            interactions.append((values.group(1),values.group(2)))
  f.close()
  return result_list

def phmmer_tuples_to_csv(phmmer_tuples, csv_filename):
  ''' 
      takes output of parse_phmmer and saves to csv
  '''
  with open(csv_filename, 'w+') as f:
    fwriter = csv.writer(f, delimiter=',')
    fwriter.writerow(('Node 1', 'Node 2', 'E-value', 'log10 E-value', 'score'))
    for row in phmmer_tuples:
      fwriter.writerow(row)
  
  # main

def main():
  # convert CSV to GI list
  gi_list = rodeo_csv_to_list(arch_file,args.c,v,False)
  
  # post GI list to Entrez
  search_results = Entrez.read(Entrez.epost("protein", id=",".join(gi_list)))
  webenv = search_results["WebEnv"]
  query_key = search_results["QueryKey"]
  
  # perform EFetch on GI list and download results, saving to file
  fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", webenv=webenv, query_key=query_key)
  data = fetch_handle.read()
  fetch_handle.close()
  if v: print data
  with open(fasta_file, 'w+') as f:
    f.write(data)
  phmmerall(fasta_file, hmmer_file)
  phmmer_tuples_to_csv(parse_phmmer(hmmer_file, v),csv_file)

if __name__ == '__main__':
  main()
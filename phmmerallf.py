#! /usr/bin/env python

'''
phmmerallf.py
  takes an input FASTA file and performs phmmer exhaustively
  then compares all sequences to an input pHMM database
  and then optionally gives a sequence-similarity network (SSN)
  in XGMML format suitable for Cytoscape
  
  contact: jonathan.tietz@gmail.com
'''

# import modules
import csv, re, os, math, argparse
from Bio import Entrez, Seq, SeqIO	# BioPython

# input options
parser = argparse.ArgumentParser(description='phmmerall -- Generates all v all phmmer for sequence similarity network generation')
parser.add_argument('-i', help='Input protein FASTA file', required=True)                        # input filename -- FASTA file
parser.add_argument('-e', help='Email address', required=True)                                          # email address is required for Entrez identification
parser.add_argument('-o', help='Output file base (no suffix)', default='phmmerall_out') 
parser.add_argument('--verbose', help='Give verbose output (for debugging)', action='store_true', default=False) 
parser.add_argument('-c', help='Column containing GI numbers', type=int, default=0)
parser.add_argument('-d', help='pHMM database', required=True)
parser.add_argument('--xgmml', help='Output .xgmml file for Cytoscape', action='store_true')
threshold = parser.add_mutually_exclusive_group()
threshold.add_argument('--eval', help='E-value threshold for phmmer', type=float, default=10)
threshold.add_argument('--score', help='Score threshold for phmmer', type=float)
args = parser.parse_args()

# initialize

Entrez.email = args.e
arch_file = args.i
fasta_file = args.o + '.fasta'
hmmer_file = args.o + '_phmmer.txt'
hmmscan_file = args.o + '_hmmscan.txt'
phmmer_csv_file = args.o + '_phmmer.csv'
hmmscan_csv_file = args.o + '_hmmscan.csv'
hmmscan_top_csv_file = args.o + '_top_hmmscan.csv'
network_file = args.o + '_network.xgmml'
if args.eval:
  eval_threshold = " -E " + str(args.eval)
if args.score:
  eval_threshold = " -T " + str(args.score)
v = args.verbose
hmm_db = args.d

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
  os.system("phmmer -o tmp2.tmp --tblout %s --noali %s %s %s" % (output_filename, eval_threshold, input_filename, input_filename))
  os.system("rm tmp2.tmp")
  return

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
    fwriter.writerow(('Node 1', 'Node 2', 'E-value', 'Log10 E-value', 'Score'))
    for row in phmmer_tuples:
      fwriter.writerow(row)
 
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

def hmmscanall(input_filename, hmm_database, output_filename):
  os.system("hmmscan -o tmp3.tmp --tblout %s --noali %s %s" % (output_filename, hmm_database, input_filename))
  os.system("rm tmp3.tmp")
  return

def parse_hmmscan(hmmscan_filename, verbose):
  '''
        parses hmmscan outputs as a list of tuples for each identifier
  ''' 
   # open file
  try:
    f = open(hmmscan_filename, 'rU')
  except IOError:
    print "ERROR .. Could not read file in function parse_hmmscan .. file = ", hmmscan_filename
    sys.exit()
  result_list = []
  with f:
    reader = f.readlines()
    for row in reader:
      comments = re.match(r'^\#.*$',row)
      if not comments:
        values = re.match(r'(\S+?)\s+?(\S+?)\s+?(\S+?)\s+?\S+?\s+?(\S+?)\s+?(\S+?)\s+?\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+\S+?\s+(.+)\n',row)
        if not values: # If didn't match expected regex, warn user
          if verbose: print("Encountered bad pHMMER row ...",row)
        else:
          pent = (values.group(3), values.group(1), values.group(2), values.group(4), values.group(5), values.group(6))
          result_list.append(pent)
  f.close()
  return result_list  

def hmmscan_tuples_to_csv(hmmscan_tuples, csv_filename):
  ''' 
      takes output of parse_hmmscan and saves to csv
  '''
  with open(csv_filename, 'w+') as f:
    fwriter = csv.writer(f, delimiter=',')
    fwriter.writerow(('GI Number', 'HMM Name', 'HMM ID', 'E-value', 'Score', 'HMM Description'))
    for row in hmmscan_tuples:
      fwriter.writerow(row)

def parse_fasta_to_identifier_list(fasta_filename):
  '''
       gets a list of identifiers from a standard formatted FASTA file
  '''
  # open file
  try:
    f = open(fasta_filename, 'rU')
  except IOError:
    print "ERROR .. Could not read file in function parse_fasta_to_identifier_list .. file = ", fasta_filename
    sys.exit()
  result_list = []
  with f:
    reader = f.readlines()
    for row in reader:
      if re.match(r'>',row):
        values = re.match(r'>(.+?)\s(.+?)\n',row)
        if values:
          result_list.append((values.group(1), values.group(2)))
        else:
          values = re.match(r'>(.+?)\n',row)
          result_list.append((values.group(1), values.group(1)))
  f.close()
  return result_list
   
def sortkey(item): return float(item[3])    
      
def clean_for_xgmml(hmmscan_list, edge_list, fasta_file):
  '''
       clean up a list with multiple pHMM hit nodes to give only the top one
  '''
  # generate list of unique nodes from hmmscan_list 
  nodes = []
  for item in hmmscan_list:
    nodes.append(item[0])
  nodes_with_hmms = set(nodes)

  # check if any are not in hmmscan_list
  identifiers = parse_fasta_to_identifier_list(fasta_file)
  identifers = set(identifiers)
  all_nodes = []
  nodes_without_hmms = []
  for item in identifiers:
    all_nodes.append(item[0])
  for t in all_nodes:
    if t not in nodes_with_hmms:  
      nodes_without_hmms.append(t)
  
  # find the highest-scoring pHMM for each node (as sorted by e-value)
  node_best_hmms = []
  for item in nodes_with_hmms:
    hmm_list = sorted([t for t in hmmscan_list if t[0] == item], key=sortkey)
    node_best_hmms.append(hmm_list[0])
  
  ''' 
            will still need to see about singletons ...
  '''
  
  # add arbitrary values for the pHMM-less
  for item in nodes_without_hmms:
    node_best_hmms.append((item,"none","none","none","none","none"))
  
  return node_best_hmms

def output_to_xgmml(node_list, edge_list, hmm_list, xgmml_file):
  with open(xgmml_file, "w+") as f:
    # write XGMML header
    f.write('<!-- Database: 20151013 -->\n')
    f.write('<graph label="pHMMer all vs all network" xmlns="http://www.cs.rpi.edu/XGMML">\n')
    
    # write nodes
    counter = 0             #initialize counter
    sorted_hmms = sorted(hmm_list)
    cross_references = {}
    sorted_nodes = sorted(node_list)
    for node in sorted_nodes:
      if node[0] != sorted_hmms[counter][0]:
        print 'MISMATCH ERROR!'
        sys.exit()
      f.write('<node id="%s" label="%s">\n' % ('id' + str(counter), 'id' + str(counter)))
      f.write('  <att name="Identifier" type="string" value="%s" />\n' % node[0])
      f.write('  <att name="Description" type="string" value="%s" />\n' % node[1])
      f.write('  <att name="Check ID" type="string" value="%s" />\n' % sorted_hmms[counter][0])
      if sorted_hmms[counter][1] == 'none':
        f.write('  <att name="HMM Name" type="string" value="No HMM match" />\n')
      else:
        f.write('  <att name="HMM Name" type="string" value="%s" />\n' % sorted_hmms[counter][1])
        f.write('  <att name="HMM ID" type="string" value="%s" />\n' % sorted_hmms[counter][2])
        f.write('  <att name="E-value" type="real" value="%s" />\n' % sorted_hmms[counter][3])
        f.write('  <att name="Score" type="real" value="%s" />\n' % sorted_hmms[counter][4])
        f.write('  <att name="HMM Description" type="string" value="%s" />\n' % sorted_hmms[counter][5])
      f.write('</node>\n')
      cross_references[node[0]] = 'id' + str(counter)
      counter += 1
    
    # write edges
    for edge in edge_list:
      first_id = cross_references[edge[0]]
      second_id = cross_references[edge[1]]
      f.write('<edge id="%s,%s" label="%s,%s" source="%s" target="%s">\n' % (first_id, second_id, first_id, second_id, first_id, second_id))
      f.write('  <att name="E-value" type="real" value="%s" />\n' % edge[2])
      f.write('  <att name="Log10 E-value" type="real" value="%s" />\n' % edge[3])
      f.write('  <att name="Score" type="real" value="%s" />\n' % edge[4])
      f.write('</edge>\n')
    
    # write footer
    f.write('</graph>')

  return



# main
  

def main():
  '''# convert CSV to GI list
  if v: print ' ... Converting CSV to GI list',
  gi_list = rodeo_csv_to_list(arch_file,args.c,v,False)
  if v: print ' - DONE!'

  # download proteins & save as a FASTA
  if v: print ' ... Fetching proteins via Entrez & saving as FASTA',
  gi_list_fetch_to_fasta(gi_list, fasta_file)
  if v: print ' - DONE!'
  '''
  
  fasta_file = args.i
  # perform all-vs-all phmmer and save unique results to a CSV file
  if v: print ' ... Performing all-vs-all pHMMER ...',
  '''phmmerall(fasta_file, hmmer_file)'''
  if v: print 'and saving results to CSV',
  phmmer_results = parse_phmmer(hmmer_file, v)
  phmmer_tuples_to_csv(phmmer_results, phmmer_csv_file)
  if v: print ' - DONE!'
  
  # perform hmmscan against pHMM database for all FASTAs and save results to a CSV file
  if v: print ' ... Performing search against pHMM database: %s ...' % hmm_db,
  '''hmmscanall(fasta_file, hmm_db, hmmscan_file)'''
  hmmscan_results = parse_hmmscan(hmmscan_file, v)
  if v: print 'and saving results to CSV',
  hmmscan_tuples_to_csv(hmmscan_results, hmmscan_csv_file)
  if v: print ' - DONE!'
  
  # because we want a single top pHMM hit, we'll process that before sending to XGMML
  if args.xgmml:
    if v: print ' ... Generating XGMML for Cytoscape ...',
    xgmml_nodes = clean_for_xgmml(hmmscan_results, phmmer_results, fasta_file)
    output_to_xgmml(parse_fasta_to_identifier_list(fasta_file), phmmer_results, xgmml_nodes, network_file)
    if v: print ' - DONE!'
  
if __name__ == '__main__':
  main()
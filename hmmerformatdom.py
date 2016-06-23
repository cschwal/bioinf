#! /usr/bin/env python

'''
hmmerformat.py
  takes input hmmscan output files (if using --domtblout option in hmmscan) and converts them to readable HTML
  
  contact: jonathan.tietz@gmail.com
'''

# import modules
import csv, re, os, math, argparse, operator
from Bio import Entrez, Seq, SeqIO	# BioPython

# input options
parser = argparse.ArgumentParser(description='hmmerformat - Produces an HTML readout of hmmscan results')
parser.add_argument('-i', help='Input hmmscan --domtblout output file', required=True)                        # input filename -- either GI list or CSV containing GI list
parser.add_argument('-o', help='Output file base (no suffix)', default='hmmerformat_out') 
parser.add_argument('-f', help='FASTA file containing protein sequences') 
parser.add_argument('-e', help='E-value significance threshold (default 1e-5)', default='1e-5', type=float)

args = parser.parse_args()

# initialize

input_file = args.i
output_html_file = args.o + '.html'
eval_threshold = args.e

# functions

def parse_hmmscan(hmmscan_filename):
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
        values = re.match(r'(\S+?)\s+?(\S+?)\s+?\S+?\s+?(\S+?)\s+?\S+?\s+?\S+?\s+?(\S+?)\s+?\S+?\s+?\S+?\s+?\S+?\s+?\S+?\s+?(\S+?)\s+?(\S+?)\s+?(\S+?)\s+?\S+?\s+?\S+?\s+?\S+?\s+?\S+?\s+?\S+?\s+?(\S+?)\s+?(\S+?)\s+?\S+?\s+?(.+)\n',row)
        if not values: # If didn't match expected regex, warn user
          print("Encountered bad pHMMER row ...",row)
        else:
          '''
            regex capture groups are:
              (1) Target name    (2) Target accession    (3) Query name    (4) E-value (5) c-Value (6) i-Value
              (7) Score     (8) Coord start    (9) Coord end     (10) Description
          '''
          props = (values.group(3), values.group(1), values.group(2), values.group(4), values.group(5), values.group(6), values.group(7), values.group(8), values.group(9), values.group(10))
          result_list.append(props)
  f.close()
  return result_list  

def parse_fasta(fasta_filename):
  '''
       gets a list of identifiers and sequences from a standard formatted FASTA file
  '''
  result_list = {}
  order_counter = 0
  for seq_record in SeqIO.parse(args.f, "fasta"):
    result_list[seq_record.id] = [seq_record.seq, len(seq_record), order_counter]
    order_counter += 1
  return result_list

def sortkey(item):
  return float(item[3])
  
def write_html(gene_list, hmmscan, html_file):
  with open(html_file, 'w') as f:
    # print header info
    f.write('<!DOCTYPE html>\n'+'<html lang="en">\n'+'<head>\n'
    +'<meta charset="utf-8">\n'
    +'<meta http-equiv="X-UA-Compatible" content="IE=edge">\n'
    +'<meta name="robots" content="noindex, nofollow" />\n'
    +'<meta name="googlebot" content="noindex" />\n'
    +'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
    +'<meta name="description" content="">\n'
    +'<meta name="author" content="">\n'
    +'<link rel="icon" href="../../favicon.ico">\n'
    +'<title>hmmerformat</title>\n'
    +'<!-- Latest compiled and minified CSS -->\n'
    +'<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">\n'
    +'<!-- Optional theme -->\n'
    +'<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css" integrity="sha384-fLW2N01lMqjakBkx3l/M9EahuwpSfeNvV63J5ezn3uZzapT0u7EYsXMjQV+0En5r" crossorigin="anonymous">\n'
    +'</head>\n'
    +'<body>\n'
    +'  <div class="container">\n'
    +'    <div class="row">\n'
    +'        <div class="col-md-12">\n'
    +'            <h1>Content</h1>\n')
    # obtain sorted gene list
    unsorted_gene_list = []
    for gene in gene_list:
      unsorted_gene_list.append((gene, gene_list[gene][2]))
    unsorted_gene_list = set(unsorted_gene_list)
    sorted_gene_list = sorted(unsorted_gene_list, key=operator.itemgetter(1))
    # print genes
    for gene in sorted_gene_list:
      # get hmms
      hmmscan_specific = [x for x in hmmscan if x[0] == gene[0]]
      hmmscan_specific = sorted(hmmscan_specific, key=sortkey)
      how_many_hmms = len(hmmscan_specific)
      hmmscan_threshold = [x for x in hmmscan_specific if float(x[3]) < eval_threshold]
      how_many_significant_hmms = len(hmmscan_threshold)
      # print sequence information
      f.write( "<h2>%s</h2>" % gene[0])
      f.write('                  <p><b>%d</b> pHMM matches, <b>%d</b> of which meet %g threshold</p>\n' % (how_many_hmms, how_many_significant_hmms, eval_threshold))
      # print progress bar of significant hits
      if how_many_significant_hmms > 0: f.write('<div class="panel panel-default"><div class="panel-heading">Domain locations for significant pHMM hits</div><div class="panel-body">')
      for row in hmmscan_threshold:
        # calculate start and end coordinates
        if gene[0] == 0: 
          start_coord = 0
        else: 
          start_coord = 100 * int(row[7]) / int(gene_list[gene[0]][1])
        if gene[0] == 0:
          end_coord = 100
        else: 
          end_coord = 100 * int(row[8]) / int(gene_list[gene[0]][1])
        accession = row[2]
        f.write('\n      <div class="progress">\n'
        +'<div class="progress-bar" style="background-image:none;background-color:#f5f5f5;-webkit-box-shadow:inset 0 1px 2px rgba(0,0,0,.1);box-shadow:inset 0 1px 2px rgba(0,0,0,.1);width: %f%%">&nbsp;</div>\n' % start_coord
        +'<div class="progress-bar" style="width: %f%%">%s</div></div>' % (float(end_coord - start_coord), accession))
      if how_many_significant_hmms > 0: f.write('</div></div>')
      f.write( '<div class="panel panel-default">')
      f.write( '<div class="panel-heading">Sequence (%d aa)</div>' % gene_list[gene[0]][1])
      f.write( '<div class="panel-body" style="word-wrap:break-word;"><p><samp>%s</samp></p><p><a href="http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&LINK_LOC=protein&PAGE_TYPE=BlastSearch">BLASTP against nr database</a></p></div>' % (gene_list[gene[0]][0], gene_list[gene[0]][0]))
      f.write( '</div>\n')
      # print hmmscan information
      f.write('                <div class="panel panel-default">\n'
      +'                  <div class="panel-heading">pHMM Matches</div>\n'
      +'                    <table class="table table-striped">\n'
      +'                        <tr>\n'
      +'                            <th>Target Name</th>\n'
      +'                            <th>Accession</th>\n'
      +'                            <th>E-value</th>\n'
      +'                            <th>Individual E-value<br />(i-Value)</th>\n'
      +'                            <th>Condensed E-value<br />(c-Value)</th>\n'
      +'                            <th>Score</th>\n'
      +'                            <th>Domain <br />Start</th>\n'
      +'                            <th>Domain <br />End</th>\n'
      +'                            <th>Description</th>\n'
      +'                        </tr>')
      # for each gene, print pHMM information
      for row in hmmscan_specific:
        if float(row[3]) > eval_threshold:
          muted = ' class="text-muted"'
        else:
          muted = ''
        f.write('                        <tr%s>\n' % muted
        +'                            <td>%s</td>\n' % row[1]
        +'                            <td>%s</td>\n' % row[2]
        +'                            <td>%s</td>\n' % row[3]
        +'                            <td>%s</td>\n' % row[5]
        +'                            <td>%s</td>\n' % row[4]
        +'                            <td>%s</td>\n' % row[6]
        +'                            <td>%s</td>\n' % row[7]
        +'                            <td>%s</td>\n' % row[8]
        +'                            <td>%s</td>\n' % row[9]
        +'                        </tr>')
      if how_many_hmms == 0:
        f.write('                        <tr><td colspan="9">no pHMMs matched</td></tr>')
      f.write('                    </table>\n'
      +'                </div>')
    # print footer
    f.write('            </div>\n'
    +'        </div>\n'
    +'    </div>\n'
    +'  </body>\n'
    +'  <!-- /container -->\n'
    +'    <!-- Bootstrap core JavaScript\n'
    +'    ================================================== -->\n'
    +'    <!-- Placed at the end of the document so the pages load faster -->\n'
    +'    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>\n'
    +'    <script>window.jQuery || document.write(\'<script src="../../assets/js/vendor/jquery.min.js"><\/script>\')</script>\n'
    +'    <script src="../../dist/js/bootstrap.min.js"></script>\n'
    +'    <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->\n'
    +'    <script src="../../assets/js/ie10-viewport-bug-workaround.js"></script>\n'
    +'    </html>\n'
    +'<!-- Latest compiled and minified JavaScript -->\n'
    +'<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>')
  return
   
# main
  

def main():
  # convert FASTA to identifier list (if applicable)
  gene_list = parse_fasta(args.f)
  
  # fetch hmmscan results for each gene
  hmmscan_results = parse_hmmscan(input_file)
  
  # print HTML for each
  write_html(gene_list, hmmscan_results, output_html_file)

  
if __name__ == '__main__':
  main()
#!/usr/bin/python

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO, Entrez
import argparse

# import options
parser = argparse.ArgumentParser(description='Basic ORF Diagram Drawer')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', help='Input nucleotide record')
group.add_argument('-f', help='Input GenBank-format file')
parser.add_argument('-e', help='Email address (required for Entrez)', required=True)
parser.add_argument('--start', help='Start coordinate', default=0)
parser.add_argument('--end', help='End coordinate')
parser.add_argument('-o', help='Output filename base (no extension!)', required=True)
parser.add_argument('--pdf', help='Make PDF', action='store_true')
parser.add_argument('--png', help='Make PNG', action='store_true')
parser.add_argument('--svg', help='Make SVG', action='store_true')
parser.add_argument('--eps', help='Make EPS', action='store_true')
parser.add_argument('-a', help='Feature to annotate', choices=['CDS', 'gene'], default='CDS')
parser.add_argument('--height', help='Height in cm', type=float, default=2)
parser.add_argument('--width', help='Width in cm', type=float, default=10)
parser.add_argument('--shape', help='Shape', choices=['circular', 'linear'], default='linear')
args = parser.parse_args()

Entrez.email = args.e

if args.i:
	handle = Entrez.efetch(db="nuccore", id=args.i, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()
if args.f:
	record = SeqIO.read(args.f, "genbank")

gd_diagram = GenomeDiagram.Diagram(record.description)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

# alternate feature colors
for feature in record.features:
	if feature.type != args.a:
		continue
	if len(gd_feature_set) % 2 == 0:
		color = colors.blue
	else:
		color = colors.lightblue
	gd_feature_set.add_feature(feature, color=color, label=True, label_angle=0, label_size=2, label_position="middle", sigil="BIGARROW")

# set start and end from input coordinates, or default to entire record if not given
if args.end:
	endDiag = args.end
else:
	endDiag = len(record)

# generate diagram format	
gd_diagram.draw(format=args.shape, pagesize=(args.height*cm,args.width*cm), fragments=1, start=args.start, end=endDiag)

# save files -- defaults to a PDF output if no formats are specified, but otherwise only gives the ones desired
if args.pdf:
	pdfname = args.o + ".pdf"
	gd_diagram.write(pdfname, "PDF")
if args.png:
	pngname = args.o + ".png"
	gd_diagram.write(pngname, "PNG")
if args.eps:
	epsname = args.o + ".eps"
	gd_diagram.write(epsname, "EPS")
if args.svg:
	svgname = args.o + ".svg"
	gd_diagram.write(svgname, "SVG")
if not args.pdf and not args.png and not args.eps and not args.svg:
	pdfname = args.o + ".pdf"
	gd_diagram.write(pdfname, "PDF")
#!/usr/bin/env python
import mimeo
import argparse
import os

def mainArgs():
	parser = argparse.ArgumentParser(description='Find all high-identity segments shared between genomes.',
									 prog='mimeo-map')
	# Input options
	parser.add_argument('--adir',type=str,default='A_genome_split',help='Name of directory containing sequences from A genome.')
	parser.add_argument('--bdir',type=str,default='B_genome_split',help='Name of directory containing sequences from B genome.')
	parser.add_argument('--afasta',type=str,default=None,help='A genome as multifasta.')
	parser.add_argument('--bfasta',type=str,default=None,help='B genome as multifasta.')
	parser.add_argument('-r','--recycle',action="store_true",help='Use existing alignment "--outfile" if found.')
	# Output options
	parser.add_argument('-d', '--outdir',type=str,default=None,help='Write output files to this directory. (Default: cwd)')
	parser.add_argument('--gffout',type=str,default="mimeo_B_in_A.gff3",help='Name of GFF3 annotation file.')
	parser.add_argument('--outfile',type=str,default="mimeo_alignment.tab",help='Name of alignment result file.')
	parser.add_argument('--verbose',action="store_true",default=False,help='If set report LASTZ progress.')
	parser.add_argument('--label',type=str,default="BHits",help='Set annotation TYPE field in gff.')
	parser.add_argument('--prefix',type=str,default="BHit",help='ID prefix for B-genome hits annotated in A-genome.')
	# LASTZ options
	parser.add_argument('--lzpath',type=str,default="lastz",help='Custom path to LASTZ executable if not in $PATH.')
	parser.add_argument('--minIdt',type=int,default=60,help='Minimum alignment identity to report.')
	parser.add_argument('--minLen',type=int,default=100,help='Minimum alignment length to report.')
	parser.add_argument('--hspthresh',type=int,default=3000,help='Set HSP min score threshold for LASTZ.')
	args = parser.parse_args()
	return args

def main():
	# Get cmd line args
	args = mainArgs()
	# Check for required programs.
	tools = [args.lzpath]
	missing_tools = []
	for tool in tools:
		missing_tools += mimeo.missing_tool(tool)
	if missing_tools:
		print('WARNING: Some tools required by mimeo could not be found: ' +
			  ', '.join(missing_tools))
		print('You may need to install them to use all features.')
	# Set output paths
	adir_path,bdir_path,outdir,outtab,gffout = mimeo.set_paths(adir=args.adir,bdir=args.bdir,afasta=args.afasta,bfasta=args.bfasta,outdir=args.outfile,outtab=args.outtab,gffout=args.gffout)
	# Get all B to A alignment pairs
	pairs =  mimeo.get_all_pairs(Adir=adir_path,Bdir=bdir_path)
	# Get chrm lens for GFF header
	chrLens = mimeo.chromlens(seqDir=adir_path)
	# Do not realign if outtab exists AND recycle mode is set
	if not args.recycle or not os.path.isfile(outtab):
		# Compose alignment commands
		cmds = mimeo.map_LZ_cmds(lzpath=args.lzpath,pairs=pairs,minIdt=args.minIdt,minLen=args.minLen,hspthresh=args.hspthresh,outfile=outtab,verbose=args.verbose,reuseTab=args.recycle)
		# Run alignments
		mimeo.run_cmd(cmds,verbose=args.verbose)
	#Import alignment as df
	alignments = mimeo.import_Align(infile=outtab,prefix=args.prefix,minLen=100,minIdt=95)
	# Write to GFF3
	with open(gffout, 'w') as f:
		for x in mimeo.writeGFFlines(alnDF=alignments,chrlens=chrLens,ftype=args.label):
			f.write(x)
	print("Finished!")



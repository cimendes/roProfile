'''
Generation of pan-genome profile files using Roary output.

input:
	-roary output path
	-gff path

'''

import csv, argparse, time, sys, os
from Bio import SeqIO
import operator
import pandas



class Logger(object):
	def __init__(self, out_directory):
		self.logfile = os.path.join(out_directory, "run.log")
		if os.path.isfile(self.logfile):
			print "Logfile already exists! It will be overwritten..." + "\n"
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")
	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()
	def flush(self):
		pass

def paserGenePA(filename):
	#parser for Roary's gene_presence_absence.csv file

	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)

		isolates=reader.next()[14:] #save all filenames for the isolates

		IsolateGeneList={} #dictionary contaning all genes in the pan-genome per isolate
		geneList=[] # for control purposes! should be the same size as the pan-genome obtained in Roary
		coreGenes=[] #genes present in ALL isolates

		for row in reader:
			geneList.append(row[0])
			genes=row[14:]

			#check if core (present in ALL isolates)
			if '' not in genes:
				coreGenes.append(row[0])

			for i in range(0,len(genes)):
				genename=genes[i]
				if '\t' in genename: # multiple genes per locus
					genename=genename.split('\t')
				else:
					genename=[genename]
				isolate= isolates[i]

				if isolate not in IsolateGeneList.keys():
					gene={}
					gene[row[0]]=[genename]
					IsolateGeneList[isolate]=gene
				else:
					value=IsolateGeneList[isolate]
					value[row[0]]=[genename]
					IsolateGeneList[isolate]=value

	print "pangenome size: " + str(len(geneList)) + " genes"

	print "core size: " + str(len(coreGenes)) + " genes"

	return IsolateGeneList, coreGenes

def parserGFF(gffDir, profileDict):

	print "Reading... "
	for item in os.listdir(gffDir):
		
		print item
		filename=item.replace('.gff', '')

		#cleaning temp files if they exist
		if os.path.isfile('temp_genes_gff.txt'):
			os.remove('temp_genes_gff.txt')
		if os.path.isfile('temp_contigs.fasta'):
			os.remove('temp_contigs.fasta')

		#separating the gff into 2 different files: one with the annotations and another with the conting sequences
		with open(gffDir+item, 'r') as in_handle, open('temp_genes_gff.txt', 'a') as temp_genes, open('temp_contigs.fasta', 'a') as temp_contigs:
			for line in in_handle: 
				if not line.startswith('##'):
					if '\t' in line:
						temp_genes.write(line)
					else:
						temp_contigs.write(line)
		
		#parsing the annotation file into a dictionary
		gffFiles={}
		with open('temp_genes_gff.txt', 'r') as temp_genes:
			for line in temp_genes:
				line=line.split('\t')
				ID=line[-1].split(';')
				locusID=str(ID[0].split('=')[1])
				contig=line[0]
				begining=int(line[3])
				end=int(line[4])
				location=[contig, begining, end]
				gffFiles[locusID]=location
		
		#parsing the sequence file into a SeqIO dictionary. one contig per entry
		handle = open("temp_contigs.fasta", "rU")
		records_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		handle.close()

		#removing temp files
		os.remove('temp_genes_gff.txt')
		os.remove('temp_contigs.fasta')

		#assigning the gene sequence to each gene in the data scrutute
		for genegroup, gene in profileDict[filename].items():
			if len(gene[0])>1: # multitple genes in one locus 
				multipleGenes=[]
				for item in gene[0]:
					geneInfo=gffFiles[item]
					contigSeq=records_dict[geneInfo[0]].seq
					geneseq=str(contigSeq[geneInfo[1]:geneInfo[2]])
					multipleGenes.append([item,geneseq])
				profileDict[filename][genegroup]=multipleGenes
			else:
				if gene[0][0] != '': # make sure the gene exists in this isolate
					geneInfo= gffFiles[gene[0][0]]
					contigSeq=records_dict[geneInfo[0]].seq
					geneseq=str(contigSeq[geneInfo[1]:geneInfo[2]])
					gene[0].append(geneseq)
	
	return profileDict

def doProfileSequences(sequenceDict):
	profileSeqDict={}
	profileFile={}
	header=[]

	#profile algorythm algorithm
	for isolate, genes in sequenceDict.items():
		print "file " + str(isolate)
		#print len(genes) #control - should be same size as pan-genome
		if isolate not in profileFile.keys():
			profileFile[isolate]={}
		for geneGroup, geneInfo in genes.items():
			if geneGroup not in header:
				header.append(geneGroup)
			if geneGroup not in profileSeqDict.keys():
				temp={}
				temp['LNF']=0
				profileSeqDict[geneGroup]=temp #initiate new entry with the profile number for Locus Not Found
			if len(geneInfo)>1: #multiple genes in one loci
				sequence=''
				for item in geneInfo:
					sequence+=item[1]
				if sequence not in profileSeqDict[geneGroup].keys():
					number=len(profileSeqDict[geneGroup])
					profileSeqDict[geneGroup][sequence]=number
					profileFile[isolate][geneGroup]=profileSeqDict[geneGroup][sequence]

				else: 
					profileNum=profileSeqDict[geneGroup][sequence]
					profileFile[isolate][geneGroup]=profileNum #!!!!
			else:
				if geneInfo[0][0] != '': # make sure the gene exists in this isolate
					if geneInfo[0][1] not in profileSeqDict[geneGroup].keys():
						number=len(profileSeqDict[geneGroup])
						profileSeqDict[geneGroup][geneInfo[0][1]]=number
						profileFile[isolate][geneGroup]=profileSeqDict[geneGroup][geneInfo[0][1]]

					else:
						profileNum=profileSeqDict[geneGroup][geneInfo[0][1]]
						profileFile[isolate][geneGroup]=profileNum 
				else:
					geneSeq='LNF'
					profileNum=profileSeqDict[geneGroup][geneSeq]
					profileFile[isolate][geneGroup]=profileNum 
	
	return profileSeqDict, profileFile, header

def writeFiles(profileSeqDict, profileFile, header, coregenes, corearg):
	
	if not corearg: 

		print "...creating pan-genome profile file..."

		with open ("pan-profile.tsv", 'w') as profileOutFile:
			profileOutFile.write('Isolate\t'+"\t".join(header)+'\n') #header of the file
			for isolate,profile in profileFile.items():
				toWrite=[]
				toWrite.append(isolate)
				for item in header:
					allele=str(profile[item])
					toWrite.append(allele)
				profileOutFile.write('\t'.join(toWrite)+'\n')
	else:

		print "... creating core-genome profile file..."

		with open("core-profile.tsv", 'w') as coreProfileOutFile:
			coreProfileOutFile.write('Isolate\t'+"\t".join(coregenes)+'\n') #header of the file, only with core genes
			for isolate,profile in profileFile.items():
				toWrite=[]
				toWrite.append(isolate)
				for item in coregenes:
					allele=str(profile[item])
					toWrite.append(allele)
				coreProfileOutFile.write('\t'.join(toWrite)+'\n')

	print "...creating sequence files..."

	sequenceDir =str(os.getcwd() + '/sequences/')
	if not os.path.exists(sequenceDir):
		os.makedirs(sequenceDir)
	else:
		print "sequence directory already exists! emptying directory..."
		for file in os.listdir(sequenceDir):
			os.remove(sequenceDir+file)

	for geneGroup, profile in profileSeqDict.items():
		if geneGroup in coregenes:
			with open("%score.%s.fasta" % (sequenceDir,geneGroup), 'w') as fasta:
				sorted_info=sorted(profile.items(), key=operator.itemgetter(1))
				for item in sorted_info:
					if item[1]==0: #don't write the Locus Not Found
						pass
					else:
						fasta.write('>%s\n' % item[1]) #profile number
						fasta.write(item[0]+'\n') #sequence
		else:
			with open("%saccessory.%s.fasta" % (sequenceDir,geneGroup), 'w') as fasta:
				sorted_info=sorted(profile.items(), key=operator.itemgetter(1))
				for item in sorted_info:
					if item[1]==0: #don't write the Locus Not Found
						pass
					else:
						fasta.write('>%s\n' % item[1]) #profile number
						fasta.write(item[0]+'\n') #sequence

def transposeMatrix(filename):
	df=pandas.read_csv(filename, sep='\t', header=0, index_col=0)
	transpose=df.transpose()
	transpose.to_csv(path_or_buf='gene_presence_absence_profile.tsv', sep='\t')

def main():

	version=1.0

	parser = argparse.ArgumentParser(description='Generation of pan-genome profile files using Roary output. By default, it will generate a profile for the full pan-genome, with Locus Not Fund represented as 0 ', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-r', '--roary', help='Path to directory containing all output files from Roary (https:/sanger-pathogens.github.io/Roary)')
	parser.add_argument('-d', '--gffdir', help='Path to directory containing all gff files used in the Roary analysis.')
	parser.add_argument('-cg','--core', help='Generate profile file for the core-genome only', required= False, default=False, action='store_true')
	parser.add_argument('-t', '--transpose', help= 'transpose the gene presence absence rtab file from roary to be used as profile', required=False, default=False, action='store_true')
	parser.add_argument('--version', help='Display version, and exit.', default=False, action='store_true')

	args = parser.parse_args()

	#version
	if args.version:
		print sys.stdout, "Current version: %s" %(version)
		sys.exit(0)

	#check if something is missing
	if  args.roary is None or args.gffdir is None:
		parser.print_usage()
		if args.roary is None:
			print "error: argument -r/--roary is required"
		if args.gffdir is None:
			print "error: argument -d/--gffdir is required"
		sys.exit(1)

	start_time = time.time()

	sys.stdout = Logger("./")


	#START

	if args.transpose:
		print 'transposing gene presence matrix...'
		transposeMatrix(args.roary+'gene_presence_absence.Rtab')

	print 'parsing gene presence and absence file...'
	profileDict, coregenes =paserGenePA(args.roary + 'gene_presence_absence.csv')

	print 'parsing gff files...'
	sequenceDict=parserGFF(args.gffdir, profileDict)

	print "creating profiles..."
	profileSeqDict, profileFileDict, header=doProfileSequences(sequenceDict)

	print "writing files..."
	writeFiles(profileSeqDict, profileFileDict, header, coregenes, args.core)

	
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken,3600)
	minutes, seconds = divmod(rest, 60)
	print "Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

	print "Finished"
	sys.exit(0)


if __name__ == "__main__":
    main()
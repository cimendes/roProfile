'''
Generation of pan-genome profile files using Roary output. 

input:
	-roary output path
	-gff path

'''

import csv, argparse, time, sys, os
from Bio import SeqIO
import operator
import mpld3
from mpld3 import utils, plugins 
import matplotlib.pyplot as plt
plt.switch_backend('agg') 

class Logger(object):
	#class for run.log
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

class ClickInfo2(plugins.PluginBase):
    """Plugin for getting info on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clickinfo2", ClickInfo2);
    ClickInfo2.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo2.prototype.constructor = ClickInfo2;
    ClickInfo2.prototype.requiredProps = ["id"];
    ClickInfo2.prototype.defaultProps = {labels:null}   function ClickInfo2(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    ClickInfo2.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        labels = this.props.labels;
        obj.elements().on("mousedown",
                          function(d, i){ 
                            window.open(labels[i], '_blank')});
    }
    """
    def __init__(self, points, labels):
        self.points = points
        self.labels = labels
        """if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None"""
        self.dict_ = {"type": "clickinfo2","id": utils.get_id(points, "pts"),"labels": labels}

def paserGenePA(filename):
	#parser for Roary's gene_presence_absence.csv file

	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)

		isolates=reader.next()[14:] #save all filenames for the isolates

		IsolateGeneList={} #dictionary contaning all genes in the pan-genome per isolate
		geneList=[] # for control purposes! should be the same size as the pan-genome obtained in Roary
		coreGenes=[] #genes present in ALL isolates
		toRemove=[] #genes with multiple alleles

		for row in reader:
			geneList.append(row[0])
			genes=row[14:]

			#check if core (present in ALL isolates)
			if '' not in genes:
				coreGenes.append(row[0])

			for i in range(0,len(genes)):
				genename=genes[i]
				if '\t' in genename: # multiple genes per locus
					toRemove.append(row[0])
					break
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

	#removing genes with multiple alleles
	for isolate, geneGroupList in IsolateGeneList.items():
		for geneGroup in geneGroupList.keys():
			if geneGroup in toRemove:
				del geneGroupList[geneGroup]


	print "\tpangenome size: " + str(len(geneList)) + " genes"

	print "\tcore size: " + str(len(coreGenes)) + " genes"

	print "\t\t%s loci to be removed due to multiple alleles.. " % (str(len(toRemove)))
	
	#removal of loci with multiple alleles
	for item in coreGenes:
		if item in toRemove:
			coreGenes.remove(item)

	return IsolateGeneList, coreGenes, toRemove

def parserGFF(gffDir, profileDict, threshold, coregenes):
	#parser for GFF3 files, retrieving the geneIDs and sequence coords

	print "Reading... "
	for item in os.listdir(gffDir):
		
		print '\t' + str(item)
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
		
		#parsing the feature file into a dictionary
		gffFiles={}
		with open('temp_genes_gff.txt', 'r') as temp_genes:
			for line in temp_genes:
				line=line.split('\t')
				ID=line[-1].split(';')
				locusID=str(ID[0].split('=')[1])
				contig=line[0]
				begining=int(line[3])-1 #to get the full sequence
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
			if gene[0][0] != '': # make sure the gene exists in this isolate
				geneInfo= gffFiles[gene[0][0]]
				contigSeq=records_dict[geneInfo[0]].seq
				geneseq=str(contigSeq[geneInfo[1]:geneInfo[2]])
				gene[0].append(geneseq)

	#filter geneGroups by size (0.2+-mode?)
	geneGroupLen={}
	for filename,geneGroupDict in profileDict.items():
		lengthSeq=[]
		for geneGroup, alleles in geneGroupDict.items():
			if  len(alleles[0]) >1: #has sequence
				if geneGroup not in geneGroupLen.keys():
					geneGroupLen[geneGroup]=[len(alleles[0][1])]
				else:
					geneGroupLen[geneGroup].append(len(alleles[0][1]))

	#calculating the mode for each loci
	print 'Calculating Gene Size Mode...'
	for geneGroup, lens in geneGroupLen.items():
		moda=max(set(lens), key=lens.count)
		geneGroupLen[geneGroup]=moda

	print "\tUsing a threshold of " + str(threshold)
	groupsToRemove=[]
	for filename,geneGroupDict in profileDict.items():
		for geneGroup, alleles in geneGroupDict.items():
			if  len(alleles[0]) >1: #has sequence
				if len(alleles[0][1]) > geneGroupLen[geneGroup]+(geneGroupLen[geneGroup]*threshold) or len(alleles[0][1]) < geneGroupLen[geneGroup]-(geneGroupLen[geneGroup]*threshold):
					groupsToRemove.append(geneGroup)

	print "Genes to be removed due to allele size variation: " + str(len(set(groupsToRemove)))
	count=0
	for item in set(groupsToRemove):
		if item in coregenes:
			count+=1
			coregenes.remove(item)
	print "\t... of which %s are in core genome." % (count)

	print "removing genes..."
	for filename,geneGroupDict in profileDict.items():
		for geneGroup in geneGroupDict.keys():
			if geneGroup in set(groupsToRemove):
				del geneGroupDict[geneGroup]
				
	return profileDict, coregenes, groupsToRemove

def doProfileSequences(sequenceDict):
	#function for allele number atribution

	profileSeqDict={}
	profileFile={}
	header=[]

	#profile algorythm algorithm
	for isolate, genes in sequenceDict.items():
		print "\tfile " + str(isolate)
		#print len(genes) #control - should be same size as pan-genome
		if isolate not in profileFile.keys():
			profileFile[isolate]={}
		for geneGroup, geneInfo in genes.items():
			#print geneInfo
			if geneGroup not in header:
				header.append(geneGroup)
			if geneGroup not in profileSeqDict.keys():
				temp={}
				temp['LNF']=0
				profileSeqDict[geneGroup]=temp #initiate new entry with the profile number for Locus Not Found
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
	#function to write profile and sequence files

	print "\t...creating pan-genome profile file..."

	with open ("pan-profile.tsv", 'w') as profileOutFile:
		profileOutFile.write('Isolate\t'+"\t".join(header)+'\n') #header of the file
		for isolate,profile in profileFile.items():
			toWrite=[]
			toWrite.append(isolate)
			for item in header:
				allele=str(profile[item])
				toWrite.append(allele)
			profileOutFile.write('\t'.join(toWrite)+'\n')

	print "\t\tpan-genome profile size: " + str(len(header)) + " genes"

	if corearg:
		print "\t... creating core-genome profile file..."

		with open("core-profile.tsv", 'w') as coreProfileOutFile:
			coreProfileOutFile.write('Isolate\t'+"\t".join(coregenes)+'\n') #header of the file, only with core genes
			for isolate,profile in profileFile.items():
				toWrite=[]
				toWrite.append(isolate)
				for item in coregenes:
					allele=str(profile[item])
					toWrite.append(allele)
				coreProfileOutFile.write('\t'.join(toWrite)+'\n')

		print "\t\tcore-genome profile size: " + str(len(coregenes)) + " genes"


	print "\t...creating sequence files..."

	sequenceDir =str(os.getcwd() + '/sequences/')
	if not os.path.exists(sequenceDir):
		os.makedirs(sequenceDir)
	else:
		print "\t\tsequence directory already exists! emptying directory..."
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
	#function to obtain the gene presence and absence profile
	import pandas

	df=pandas.read_csv(filename, sep='\t', header=0, index_col=0)
	transpose=df.transpose()
	transpose.to_csv(path_or_buf='gene_presence_absence_profile.tsv', sep='\t')


def distributionGraph(filename):
	#function to obtain the interactive gene frequency graph for the pangenome, in html
	genes=[]
	coreSize=0

	with open(filename, 'r') as fh:
		reader = csv.reader(fh, delimiter='\t')
		isolates=reader.next()[1:] 

		for row in reader:
			geneName=row[0]
			numbers = [ int(x) for x in row[1:]]
			freq=sum(numbers)
			genes.append((geneName, freq))
			if freq == len(isolates):
				coreSize+=1
	x=[]
	y=[]
	xLabels=[]

	for item in genes:
		y.append(item[1])
		x.append(genes.index(item))
		xLabels.append(item[0])

	fig,ax = plt.subplots(figsize=(12, 9))  

	line=ax.plot(x, y, '.', color="#3F5D7D")
	plt.title("Pan-Genome Frequency")
	plt.ylabel('Frequency')
	plt.xlabel('Gene')
	plt.text(1.1, 0.9,  s=" Pan-Genome: %s genes | Core Genome: %s genes | Acessory Genome: %s genes" % (str(len(genes)), str(coreSize), str(len(genes)-coreSize)),fontsize=10, horizontalalignment='right',verticalalignment='center', transform = ax.transAxes) 
	ax.yaxis.labelpad = 40

	mpld3.plugins.connect(fig, plugins.PointLabelTooltip(line[0],labels=xLabels))

	mpld3.save_html(fig,'panGenome.html')

def makeLogFiles(removedMultiple, removedSize):
	#function to write the log file for the removed loci

	with open ('removedLoci.log', 'w') as logFile:
		logFile.write('- Removed loci due to multiple alleles:\n')
		logFile.write('\t'.join(str(s) for s in set(removedMultiple)) + '\n')
		logFile.write('- Removed loci due to size variation:\n')
		logFile.write('\t'.join(str(s) for s in set(removedSize)) + '\n')
		logFile.write('- Total loci removed: %s' % (len(set(removedMultiple))+len(set(removedSize))))

def newGenePA(filename, removedMultiple, removedSize):
	with open(filename, 'r') as csvfile:
		with open('clean_gene_presence_absence.csv', 'w') as outFile:
			for line in csvfile:
				items=line.split(',')
				geneGroup=items[0][1:-1]
				if geneGroup in removedMultiple or geneGroup in removedSize:
					pass
				else:
					outFile.write(line)

def main():

	version='1.3.0'

	parser = argparse.ArgumentParser(description='Generation of pan-genome profile files using Roary output (https:/sanger-pathogens.github.io/Roary). By default, it will generate a profile for the full pan-genome, with Locus Not Fund represented as 0.', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-r', '--roary', help='Path to directory containing all output files from Roary.')
	parser.add_argument('-d', '--gffdir', help='Path to directory containing all gff files used in the Roary analysis.')
	parser.add_argument('-c','--core', help='Generate profile file for the core-genome only, with genes present in all isolates.', required= False, default=False, action='store_true')
	parser.add_argument('-t', '--transpose', help= 'Transpose the gene presence absence rtab file from roary to be used as profile.', required=False, default=False, action='store_true')
	parser.add_argument('-f', '--frequency', help= 'Generate pan-genome frequency plot.', required=False, default=False, action='store_true')
	parser.add_argument('-th', '--threshold', nargs='?', type=float, help='Threshold for the allele size (default=0.2).', required=False, default=0.2)
	parser.add_argument('-g', '--genefile', help="Obtain a roary's gene presence and absence csv file without the removed loci.", required=False, default=False, action='store_true')
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

	print 'parsing gene presence and absence file...'
	profileDict, coregenes, removedMultiple = paserGenePA(args.roary + 'gene_presence_absence.csv')

	if args.transpose:
		print 'transposing gene presence matrix...'
		transposeMatrix(args.roary+'gene_presence_absence.Rtab')

	if args.frequency:
		print 'generating pan-genome frequency plot..'
		distributionGraph(args.roary+'gene_presence_absence.Rtab')

	print 'parsing gff files...'
	sequenceDict, coregenes, removedSize=parserGFF(args.gffdir, profileDict, args.threshold, coregenes)

	print "creating profiles..."
	profileSeqDict, profileFileDict, header=doProfileSequences(sequenceDict)

	print "writing files..."
	writeFiles(profileSeqDict, profileFileDict, header, coregenes, args.core)
	makeLogFiles(removedMultiple, removedSize)

	if args.genefile:
		print 'generating new gene presence and absence file...'
		newGenePA(args.roary + 'gene_presence_absence.csv', removedMultiple, removedSize)


	
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken,3600)
	minutes, seconds = divmod(rest, 60)
	print "\nRuntime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

	print "Finished"
	sys.exit(0)


if __name__ == "__main__":
    main()
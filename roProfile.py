'''
Generation of pan-genome profile files using Roary output.

input:
	-roary output path
	-gff path

'''

import csv, argparse, time, sys, os
from BCBio import GFF



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
		geneList=[] # TODO remove! for control!

		for row in reader:
			geneList.append(row[0])
			genes=row[14:]

			for i in range(0,len(genes)):
				genename=genes[i]
				isolate= isolates[i]

				#gene[row[0]]=[genename]
				if isolate not in IsolateGeneList.keys():
					gene={}
					gene[row[0]]=[genename]
					#print gene
					IsolateGeneList[isolate]=gene
					#print IsolateGeneList[isolate]
				else:
					#print IsolateGeneList[isolate]
					value=IsolateGeneList[isolate]
					value[row[0]]=[genename]
					IsolateGeneList[isolate]=value
		#print len(geneList)
		#print len(IsolateGeneList['SDSE_ERR109325_02062016'])
		#print IsolateGeneList
	return IsolateGeneList

def parserGFF(gffDir, profileDict):
	print "Reading... "
	for item in os.listdir(gffDir):

		in_handle = open(gffDir+item)
		
		filename=item.replace('.gff', '')
		print filename

		#for faster parsing
		listOFgenes=[]
		for gene in profileDict[filename].keys():
			listOFgenes.append((gene,profileDict[filename][gene][0]))

		for line in in_handle: 
			if not line.startswith('##'):
				if line.startswith('gnl') or line.startswith('NZ') or line.startswith('gi'):
					line=line.split('\t')
					ID=line[-1].split(';')
					locusID=str(ID[0].split('=')[1])
					#print locusID
					#print profileDict[filename]
					for gene in listOFgenes:
						if locusID== gene[1]:
							#print "found!"
							contig=line[0]
							begining=int(line[3])
							end=int(line[4])
							sequenceInfo=[contig, begining, end]
							profileDict[filename][gene[0]].append(sequenceInfo)
							#print profileDict[filename][gene[0]]
		in_handle.close()

		parseGGF4sequence=open(gffDir+item)
		for rec in GFF.parse(parseGGF4sequence):
			#print rec.id
			for geneID, genes in profileDict[filename].items():
				#print genes
				if len(genes)>1:
					if genes[1][0] == rec.id:
						geneSeq=rec.seq[genes[1][1]:genes[1][2]]
						genes.append(str(geneSeq))
						#print genes

	for key, value in profileDict.items():
		print key
		for k,v in value.items():
			print v




def main():

	version=0.01

	parser = argparse.ArgumentParser(description='Generation of pan-genome profile files using Roary output.', epilog='by C I Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-r', '--roary', help='Path to directory containing all output files from Roary (https:/sanger-pathogens.github.io/Roary)')
	parser.add_argument('-d', '--gffdir', help='Path to directory containing all gff files used in the Roary analysis.')
	parser.add_argument('--version', help='Display version, and exit.', default=False,action='store_true')

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


	#STRART
	print 'parsing gene presence and absence file...'
	profileDict=paserGenePA(args.roary + 'gene_presence_absence.csv')

	print 'parsing gff files...'
	sequenceDict=parserGFF(args.gffdir, profileDict)


	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken,3600)
	minutes, seconds = divmod(rest, 60)
	print "Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

	print "Finished"
	sys.exit(0)


if __name__ == "__main__":
    main()
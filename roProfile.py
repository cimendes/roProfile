'''
Generation of pan-genome profile files using Roary output.

input:
	-roary output path
	-gff path

'''

import csv, argparse, time, sys, os
from gff3 import Gff3 #https://pypi.python.org/pypi/gff3/0.2.0
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

		indexDict={}
		filenames=reader.next()[14:] #save all filenames
		for name in filenames:
			indexDict[filenames.index(name)]=name

		IsolateGeneList={}
		geneList=[] # TODO remove! for control!

		for row in reader:
			geneList.append(row[0])
			genes=row[14:]

			for i in range(0,len(genes)):
				genename=genes[i]
				filename= indexDict[i]

				if filename not in IsolateGeneList.keys():
					IsolateGeneList[filename]=[(row[0],genename)]
				else:
					IsolateGeneList[filename].append((row[0],genename))
	return IsolateGeneList

def parserGFF(gffDir, profileDict):

	for item in os.listdir(gffDir):
		filename=item.replace('.gff', '')
		
		listOFgenes=[]
		for gene in profileDict[filename]:
			listOFgenes.append(gene[1])
		#print listOFgenes

		in_handle = open(gffDir+item)

		realIDsList={}
		for line in in_handle: 
			if not line.startswith('##'):
					if line.startswith('gnl') or line.startswith('NZ') or line.startswith('gi'):
						line=line.split('\t')
						lala=line[-1].split(';')
						locusID=str(lala[0].split('=')[1])
						start=int(line[3])
						end=int(line[4])
						#print locusID
						if locusID in listOFgenes:
							#print locusID
							realID=line[0]
							#print realID
							if realID not in realIDsList.keys():
								realIDsList[realID] = [[locusID, start, end]]
							else:
								realIDsList[realID].append([locusID, start, end])



		#print realIDsList
		in_handle.close()

		parseGGF4sequence=open(gffDir+item)
		print 'Got Here!'
		for rec in GFF.parse(parseGGF4sequence):
			#print rec.id
			if rec.id in realIDsList.keys():
   				#print len(rec.seq)
   				#print rec.id
   				contigSeq=rec.seq

   				for value in realIDsList[rec.id]:
   					geneSeq=contigSeq[value[1]:value[2]]
   					#print len(geneSeq)
   					value.append(geneSeq)
		parseGGF4sequence.close()
		print realIDsList



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
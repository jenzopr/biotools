#-----------------------------------------------------------------------
#
# build_annotation_table.py version 1.0
# @author Jens Preußner
# @contact jens.preussner@mpi-bn.mpg.de
#
#-----------------------------------------------------------------------
import argparse, json, os, re, sys
from Bio import SeqIO
from Bio.Blast import NCBIXML

version = '1.0.1'
authors = 'Jens Preussner'

#-----------------------------------------------------------------------
# Subroutines 
#-----------------------------------------------------------------------
def guessFileFormat(f):
	ext = os.path.splitext(f)[-1].lower()
	if ext == '.xml' or ext == '.txt' or ext == '.tab':
		return(ext)
	else:
		return(None)

def saveTrinotateCDS(query, CDSid, description):
	global annotation
	if query in annotation:
		# Create a pattern for type, len, start, end and strand of the cds
		p = re.compile('type:([a-z0-9_]+).*len:([0-9]+).*'+query+':([0-9]+)-([0-9]+)\((.)\)')
		m = p.search(description)
		if m is not None:
			CDStype = m.group(1).lower()
			CDSlength = m.group(2)
			CDSstart = m.group(3)
			CDSend = m.group(4)
			CDSstrand = m.group(5)
			if 'CDS' in annotation[query]:
				cds = annotation[query]['CDS']
			else:
				cds = {}
			if CDSid not in cds:
				cds[CDSid] = {'start':CDSstart, 'end':CDSend, 'strand':CDSstrand, 'length':CDSlength, 'type':CDStype}
			else:
				cds[CDSid]['start'] = CDSstart
				cds[CDSid]['end'] = CDSend
				cds[CDSid]['strand'] = CDSstrand
				cds[CDSid]['length'] = CDSlength
				cds[CDSid]['type'] = CDStype
			annotation[query]['CDS'] = cds 
	elif args.verbose:
		print("Can't find query "+query+" in fasta/annotation file. Skipping this query..")	

def processSignalP(l):
	global annotation
	for f in l:
		for entry in f:
			fields = entry.split('\t')
			recordPattern = re.compile('cds\.(.*)\|(m\.[0-9]+)',re.IGNORECASE)
			m = recordPattern.search(fields[0])
			if m is not None :
				recordId = m.group(1).lower()
				CDSid = m.group(2).lower()
				if recordId in annotation:
					if 'CDS' in annotation[recordId]:
               		                 	cds = annotation[recordId]['CDS']
					else:
						cds = {}
					if CDSid not in cds:
						cds[CDSid] = {'signalp':[{'start':fields[3], 'end':fields[4], 'score':fields[5]}]}
					else:
						if 'signalp' in cds[CDSid]:
							cds[CDSid]['signalp'].append({'start':fields[3], 'end':fields[4], 'score':fields[5]})
						else:
							cds[CDSid]['signalp'] = [{'start':fields[3], 'end':fields[4], 'score':fields[5]}]
					annotation[recordId]['CDS'] = cds
				elif args.verbose:
					print("Can't find query "+query+" in fasta/annotation file. Skipping this query..")
			
def processTMHMM(l):
	global annotation
	for f in l:
		for entry in f:
			fields = entry.split('\t')
			recordPattern = re.compile('cds\.(.*)\|(m\.[0-9]+)',re.IGNORECASE)
			m = recordPattern.search(fields[0])
			if m is not None:
				recordId = m.group(1).lower()
				CDSid = m.group(2).lower()
				if recordId in annotation:
					if 'CDS' in annotation[recordId]:
						cds = annotation[recordId]['CDS']
					else:
						cds = {}
					helices = fields[4].split('=')[1]
					topology = fields[5].split('=')[1].rstrip('\n')
					if CDSid not in cds:
						cds[CDSid] = {'tmhmm':[{'topology':topology, 'helices':helices}]}
					else:
						if 'tmhmm' in cds[CDSid]:
							cds[CDSid]['tmhmm'].append({'topology':topology, 'helices':helices})
						else:
							cds[CDSid]['tmhmm'] = [{'topology':topology, 'helices':helices}]
					annotation[recordId]['CDS'] = cds
			elif args.verbose:
				print("Can't find query "+query+" in fasta/annotation file. Skipping this query..")

def processPfam(l):
	global annotation
	for f in l:
		for entry in f:
			fields = entry.split()
			recordPattern = re.compile('(.*)\|(m\.[0-9]+)',re.IGNORECASE)
			m = recordPattern.search(fields[3])
			if m is not None:
				recordId = m.group(1).lower()
				CDSid = m.group(2).lower()
				if recordId in annotation:
					if 'CDS' in annotation[recordId]:
						cds = annotation[recordId]['CDS']
					else:
						cds = {}
					if CDSid not in cds:
						cds[CDSid] = {'domains':[{'name':fields[0], 'acc':fields[1], 'evalue':fields[12], 'start':fields[19], 'end':fields[20]}]}
					else:
						if 'domains' in cds[CDSid]:
							cds[CDSid]['domains'].append({'name':fields[0], 'acc':fields[1], 'evalue':fields[12], 'start':fields[19], 'end':fields[20]})
						else:
							cds[CDSid]['domains'] = [{'name':fields[0], 'acc':fields[1], 'evalue':fields[12], 'start':fields[19], 'end':fields[20]}]
					annotation[recordId]['CDS'] = cds
			elif args.verbos:
				print("Can't find query "+query+" in fasta/annotation file. Skipping this query..")

def parse(file,filetype,flavour,pattern,cutoff='1e-15'):
	global annotation
	global islowqual

	def formatAnnotation(identifier, evalue, description, bitscore, program, version, database):
		return {'id':identifier, 'eval':evalue, 'description':description, 'bitscore':bitscore, 'program':program, 'version':version, 'database':database}
	
	def updateAnnotation(query, accession, description, evalue, bitscore, program, version, database, lowqual):
		global annotation
		splitAccession = accession.split('|',1)
		# If there is no annotation at all, just add the current annotation
		if len(annotation[query]) == 1:
			annotation[query][splitAccession[0]] = formatAnnotation(splitAccession[1].rstrip('|'),evalue, description, bitscore, program, version, database)
			if getGeneAlias(description) is not None:
				annotation[query]['alias'].insert(0,getGeneAlias(description))
			islowqual[query] = lowqual
		# If there is already a annotation and it is of the same quality as the current, update the annotation, if the score is better
		elif lowqual == islowqual[query]:
			if splitAccession[0] not in annotation[query]:
				annotation[query][splitAccession[0]] = formatAnnotation(splitAccession[1].rstrip('|'),evalue, description, bitscore, program, version, database)
				if getGeneAlias(description) is not None:
					annotation[query]['alias'].insert(0,getGeneAlias(description))
			elif float(evalue) < float(annotation[query][splitAccession[0]]['eval']):
				annotation[query][splitAccession[0]] = formatAnnotation(splitAccession[1].rstrip('|'),evalue, description, bitscore, program, version, database)
				if getGeneAlias(description) is not None:
					annotation[query]['alias'].insert(0,getGeneAlias(description))
		# If the existing annotation is of low quality but the current is of high quality, delete all annotations and save the high quality annotation
		elif lowqual == False and islowqual[query] == True:
			annotation[query] = {'alias':[]}
			annotation[query][splitAccession[0]] = formatAnnotation(splitAccession[1].rstrip('|'),evalue, description, bitscore, program, version, database)
			if getGeneAlias(description) is not None:
				annotation[query]['alias'].insert(0,getGeneAlias(description))
			islowqual[query] = False
	
	def cleanQuery(x):
		return x.split(' ',1)[0]
	
	def getGeneAlias(x):
		p = re.compile('GN=[a-z0-9]+',re.IGNORECASE)
		m = p.search(x)
		# Try to find previous alias of that query, splitting a potetntial defline into accession and description
		q = query.split(' ',1)[0]
		if m is not None :
			alias = m.group()[3:].lower()
			return(alias)
		else:
			return(None)
	
	if filetype == '.xml':
		from Bio.Blast import NCBIXML
		fileHandle = NCBIXML.parse(file)
		for blastRecord in fileHandle:
			for alignment in blastRecord.alignments:
				lowqual = False
				if pattern.search(alignment.title) is not None:
					lowqual = True
				if float(alignment.hsps[0].expect) <= float(cutoff):
					query = cleanQuery(blastRecord.query)
					try:
						# Remove unwanted internal BLAST identifiers starting with gln - use identifier from hit_def instead (XML specific)
						if alignment.hit_id.startswith('gnl'):
							alignment.hit_id = alignment.hit_def.split(' ',1)[0]
						updateAnnotation(query,alignment.hit_id,alignment.hit_def,alignment.hsps[0].expect,alignment.hsps[0].bits,blastRecord.application,blastRecord.version,blastRecord.database,lowqual)
					except (Exception, KeyError) as e:
						sys.stderr.write(e)
						break
				else:
					break
	if filetype == '.tab':
		for entry in file:
			fields = entry.split('\t')
			lowqual = False
			if pattern.search(fields[4]) is not None:
				lowqual = True
			if float(fields[3]) <= float(cutoff):
				query = cleanQuery(fields[0])
				try:
					updateAnnotation(query,fields[1],fields[4],fields[3],fields[2],flavour,'undef','undef', lowqual)
				except (Exception, KeyError) as e:
					sys.stderr.write(e)
					break

def processFilesInList(l,e,program):
	for f in l:
		ftype = guessFileFormat(f.name)
		if args.verbose:
			print('Processing file: '+f.name+' of type '+str(type(f)))
			if ftype is None:
				print('Could not automatically determine file type, I will assume XML format!')
		if ftype is None:
			ftype = '.xml'
		if ftype == '.txt' or ftype == '.xml' or ftype == '.tab':
			parse(f,ftype,program,pattern,e)
#-----------------------------------------------------------------------
# Add arguments and program discription
#-----------------------------------------------------------------------
metavars=('file1','file2')
parser = argparse.ArgumentParser(description='Builds an annotation table from a fasta file and different BLAST results.',epilog='Naming conventions: BLAST results should end in .xml for XML format, .tab for tabular format and .txt for plain BLAST output.')
group_in = parser.add_mutually_exclusive_group(required=True)
group_transdecoder = parser.add_argument_group(title=None, description='Results from the Trinotate annotation pipeline can be included using:')
group_in.add_argument('--fasta','-f', type=argparse.FileType('r'), help='A file with sequences in fasta format. Sequence identifiers have to match those in the BLAST results.',metavar='sequences.fasta')
group_in.add_argument('--annotation', '-a', type=argparse.FileType('r'), help='An existing json annotation file from previous runs of this tool.',metavar='annotation.json')
parser.add_argument('--out','-o',required=True, type=argparse.FileType('w'), help='A writeable file to which output is directed.',metavar='output.tabular')
group_transdecoder.add_argument('--cds', type=argparse.FileType('r'), help='Transdecoders CDS file (required for use with Trinotate results).', metavar='transdecoder.cds', default=sys.stdin)
group_transdecoder.add_argument('--signalp', nargs='*', type=argparse.FileType('r'), help='Trinotates output for SignalP in tabular format.', metavar=metavars, default=sys.stdin)
group_transdecoder.add_argument('--tmhmm', nargs='*', type=argparse.FileType('r'), help='Trinotates output for TMHMM in tabular format.', metavar=metavars, default=sys.stdin)
group_transdecoder.add_argument('--pfam', nargs='*', type=argparse.FileType('r'), help='Trinotates output for Pfam search in tabular format.', metavar=metavars, default=sys.stdin)
parser.add_argument('--blastp', nargs='*', type=argparse.FileType('r'), help='BLAST results from BLASTP searches in plain, tabular or XML format (default).', metavar=metavars, default=sys.stdin)
parser.add_argument('--blastn', nargs='*', type=argparse.FileType('r'), help='BLAST results from BLASTN searches in plain, tabular or XML format (default).', metavar=metavars, default=sys.stdin)
parser.add_argument('--blastx', nargs='*', type=argparse.FileType('r'), help='BLAST results from BLASTX searches in plain, tabular or XML format (default).', metavar=metavars, default=sys.stdin)
parser.add_argument('--tblastn', nargs='*', type=argparse.FileType('r'), help='BLAST results from TBLASTN searches in plain, tabular or XML format (default).', metavar=metavars, default=sys.stdin)
parser.add_argument('--tblastx', nargs='*', type=argparse.FileType('r'), help='BLAST results from TBLASTX searches in plain, tabular or XML format (default).', metavar=metavars, default=sys.stdin)
parser.add_argument('--eblastp', type=float, help='E-value treshold for BLASTP results (default 1e-15).',default=1e-15)
parser.add_argument('--eblastn', type=float, help='E-value treshold for BLASTN results (default 1e-15).',default=1e-15)
parser.add_argument('--eblastx', type=float, help='E-value treshold for BLASTX results (default 1e-15).',default=1e-15)
parser.add_argument('--etblastn', type=float, help='E-value treshold for TBLASTN results (default 1e-15).',default=1e-15)
parser.add_argument('--etblastx', type=float, help='E-value treshold for TBLASTX results (default 1e-15).',default=1e-15)
parser.add_argument('--lowqual', '-l', help='A comma-separated list of keywords that mark a annotation as low quality if found in the hit description. Low quality annotations are discarded, if a high quality annotation is found. Defaults to predicted, cdna, clone, riken, hypothetical, uncharacterized.', default='predicted,cdna,clone,riken,hypothetical,uncharacterized')
#parser.add_argument('--noexclude','-n', help='Disables exclusion of a hit containing a keyword given by the --exclude/-e option.', action='store_true')
parser.add_argument('--verbose','-v', help='Increase output verbosity.', action='store_true')
parser.add_argument('--json','-j',type=argparse.FileType('w'), help='Output annotation in JSON format.')
parser.add_argument('--version', action='version', version='%(prog)s version '+version+' by '+authors)
args = parser.parse_args()

#-----------------------------------------------------------------------
# Needed variables for annotation table output
#----------------------------------------------------------------------
annotation = {}
islowqual = {}
#----------------------------------------------------------------------
# Read input
# a) FASTA file for a new annotation table
# b) An existing annotation table for updating annotations
# c) An (optional) file with CDS sequences
#----------------------------------------------------------------------
if args.fasta:
	for record in SeqIO.parse(args.fasta, 'fasta'):
		annotation[record.id] = {}
		annotation[record.id]['alias'] = []
		annotation[record.id]['cds'] = []
		islowqual[record.id] = True

if args.annotation:
	annotation = json.load(args.annotation)

if args.cds:
	recordPattern = re.compile('cds\.(.*)\|(m\.[0-9]+)',re.IGNORECASE)
	for record in SeqIO.parse(args.cds, 'fasta'):
		m = recordPattern.search(record.id)
		if m is not None :
			recordId = m.group(1).lower()
			CDSid = m.group(2)
			saveTrinotateCDS(recordId, CDSid, record.description)

#----------------------------------------------------------------------
# Parse exclude list and generate a regular expression pattern
#----------------------------------------------------------------------
pattern = re.compile('|'.join(args.lowqual.split(',')),re.IGNORECASE)
if args.verbose:
	print('The low quality pattern is: '+str(pattern))

#----------------------------------------------------------------------
# Read BLAST results and store them in the annotation
# variable if threshholds like evalue, bitscore... are met
#----------------------------------------------------------------------
if type(args.blastp) is list:
	processFilesInList(args.blastp,args.eblastp,'blastp')
if type(args.blastn) is list:
	processFilesInList(args.blastn,args.eblastn,'blastn')
if type(args.tblastn) is list:
	processFilesInList(args.tblastn,args.etblastn,'tblastn')
if type(args.tblastx) is list:
	processFilesInList(args.tblastx,args.etblastx,'tblastx')
if type(args.blastx) is list:
	processFilesInList(args.blastx,args.eblastx,'blastx')
if type(args.signalp) is list:
	processSignalP(args.signalp)
if type(args.tmhmm) is list:
	processTMHMM(args.tmhmm)
if type(args.pfam) is list:
	processPfam(args.pfam)
#----------------------------------------------------------------------
# Print results
#----------------------------------------------------------------------

if args.json is not None:
	json.dump(annotation,args.json,indent=1)

for query in sorted(annotation):
	if not annotation[query]:
		continue
	args.out.write(query+'\t')
	if annotation[query]['alias'] is not None:
		args.out.write(','.join(list(set(annotation[query]['alias']))))
	args.out.write('\t')
	annotationList = []
	filtered = dict((key,value) for key,value in annotation[query].items() if key != 'alias')
	for type in sorted(filtered, key=lambda x: str(filtered.get(x)['eval'])):
		annotationList.append(type+'|'+annotation[query].get(type)['id'])
	args.out.write(','.join(annotationList))
	args.out.write('\n')


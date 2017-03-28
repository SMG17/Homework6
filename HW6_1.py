import gffutils

#db = gffutils.create_db('saccharomyces_cerevisiae.gff', dbfn='yeast.db', merge_strategy = "merge", force = True, verbose = True)

db = gffutils.FeatureDB('yeast.db') #calling the database

#for type in db.featuretypes():
#	print type

#for mRNA in db.features_of_type("mRNA"):
	#print mRNA.start, mRNA.stop	
	#print mRNA["Name"]
	#for CDS in db.children(mRNA, featuretype = "CDS"):
		#print CDS.start, CDS.stop
	#for CDS in db.features_of_type("CDS"):
		#for parent in db.parents(CDS):
			#print parent.featuretype
mRNA_total = 0.0
mRNA_with_intron = []
for mRNA in db.features_of_type("mRNA"):
	mRNA_total += 1
	for intron in db.children(mRNA, featuretype = 'intron'):
		mRNA_with_intron.append(mRNA["Name"])

print len(mRNA_with_intron)/mRNA_total

output = open("gene_lengths.txt","w")
for mRNA in db.features_of_type("mRNA"):
	mRNA_length = 0.0	
	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDS_length = CDS.stop - CDS.start +1
		mRNA_length += CDS_length
	output.write("%s\t%d\n" % (mRNA["Name"],mRNA_length))

output.close()
	

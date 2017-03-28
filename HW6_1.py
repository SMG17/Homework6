import gffutils

#db = gffutils.create_db('saccharomyces_cerevisiae.gff', dbfn='yeast.db', merge_strategy = "merge", force = True, verbose = True)

db = gffutils.FeatureDB('yeast.db') #calling the database

for type in db.featuretypes():
	print type

#for mRNA in db.features_of_type("mRNA"):
	#print mRNA.start, mRNA.stop	
	#print mRNA["Name"]
	#for CDS in db.children(mRNA, featuretype = "CDS"):
		#print CDS.start, CDS.stop
	#for CDS in db.features_of_type("CDS"):
		#for parent in db.parents(CDS):
			#print parent.featuretype

for gene in db.features_of_type("gene"):
	for i in db.children(gene, featuretype = 'intron', order_by = 'start'):
		print i

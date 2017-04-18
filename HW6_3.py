import gffutils
import pysam

db = gffutils.FeatureDB('yeast.db')
bamFile = pysam.AlignmentFile("output.sorted.bam","rb")
output = open("FPKM.txt","w")

ReadCountPerGene = []
GeneNames = []
TotalReadCount = 0.0
LengthPerGene = []

for mRNA in db.features_of_type("mRNA"):
	if mRNA.chrom == "chrmt":
		break
	GeneNames.append(mRNA["Name"])
	count = bamFile.count(mRNA.chrom,mRNA.start,mRNA.stop)
	ReadCountPerGene.append(count)
	TotalReadCount += count
	mRNA_length = 0.0
	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDS_length = CDS.stop - CDS.start + 1
		mRNA_length += CDS_length
	LengthPerGene.append(mRNA_length)

for i in range(len(GeneNames)):
	FPKM = float(ReadCountPerGene[i])/(LengthPerGene[i]*TotalReadCount)*10**9
	output.write("%s\t%d\n" % (GeneNames[i],FPKM))

output.close()

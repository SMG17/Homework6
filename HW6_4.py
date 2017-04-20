import gffutils
import pysam
import numpy as np
import scipy.stats as st

db = gffutils.FeatureDB('yeast.db')
bamFileBY = pysam.AlignmentFile("outputBY.sorted.bam","rb")
bamFileRM = pysam.AlignmentFile("outputRM.sorted.bam","rb")
output = open("BY_RM_counts.txt","w")

GeneNames = []
k_BY = []
k_RM = []
N_BY = 0.0
N_RM = 0.0
LengthPerGene = []

for mRNA in db.features_of_type("mRNA"):
	if mRNA.chrom == "chrmt":
		continue
	GeneNames.append(mRNA["Name"])
	countBY = bamFileBY.count(mRNA.chrom,mRNA.start,mRNA.stop) #k11 or countBY
	countRM = bamFileRM.count(mRNA.chrom,mRNA.start,mRNA.stop) #k21 or countRM
	k_BY.append(countBY) 
	k_RM.append(countBY)
	N_BY += countBY #N1 or TotalReadCountBY
	N_RM += countRM #N2 or TotalReadCountRM

for i in range(len(GeneNames)):	
	lambda_mRNA = (float(k_BY[i])/N_BY+float(k_RM[i])/N_RM)/2
	#compute lambdas for each experiment separately
	lambda_BY = float(k_BY[i])/N_BY
	lambda_RM = float(k_RM[i])/N_RM
	#compute probability of the data uner the null
	prob_data_null = st.poisson.pmf(k_BY[i],N_BY*lambda_mRNA)*st.poisson.pmf(k_RM[i],N_RM*lambda_mRNA)
	#compute probability of the date under the alt
	prob_data_alt = st.poisson.pmf(k_BY[i],N_BY*lambda_BY)*st.poisson.pmf(k_RM[i],N_RM*lambda_RM)
	#compute the LRT
	LRT = 2*(np.log(prob_data_alt) - np.log(prob_data_null))
	#compute the p-value
	p_value = st.chi2.sf(LRT,df=1)
	output.write("%s\t%d\t%d\t%f\n" % (GeneNames[i],k_BY[i],k_RM[i],p_value))

output.close()

#@author narumeena
#@description  fast mapping of chrom cords 


import csv
import pandas as pd
import multiprocessing 
from joblib import Parallel, delayed

import pandas as pd
import boto3
import io
from itertools import product
from functools import partial



bucket = 'com-rosettahub-default-narendra.meena'
s3 = boto3.client('s3',aws_access_key_id='AKIAJWADVB5BCJNZSYMQ',aws_secret_access_key='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')
#obj = s3.get_object(Bucket=bucket, Key='prabakaranLab/cosmic/noncoding/CosmicNonCodingVariants-14-12-17-SP.vcf')


#vcf file, list of mutation 

vcfMutation                 = pd.read_csv(io.BytesIO(s3.get_object(Bucket=bucket, Key='prabakaranLab/cosmic/noncoding/CosmicNonCodingVariants-14-12-17-SP.vcf')['Body'].read()),delimiter='\t', dtype='unicode',skiprows=13)
vcfMutation["#CHROM"]       = 'chr' + vcfMutation["#CHROM"].astype(str)
chromList                   = vcfMutation["#CHROM"].unique()

#print(vcfMutation)

#mutationChrom = pd.DataFrame()
#sORfsChrom    = pd.DataFrame()
#print(vcfMutation.loc[vcfMutation['#CHROM'] == 'chr1'])



#sORFs file, list of sORFS 
sORFs = pd.read_csv(io.BytesIO(s3.get_object(Bucket=bucket, Key='prabakaranLab/arka_replication/human_PLsorf_database.txt')['Body'].read()),delimiter='\t', dtype='unicode')

#print(sORFs)

#print(sORFs.Chromosome.unique())

num_cores = multiprocessing.cpu_count()
    
#results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)


def mappingMutation(mut,sORfsChrom, mutationChrom ):
    #for orf in sORfsChrom.index.tolist():
    for orf in range(len(sORfsChrom.Chromosome)):
        #print(sORfsChrom)
        with open('/home/ngs/narumeena/Documents/arkarplication/sorfs/noncoding/CosmicNonCodingVariantsMappedToHumanPLsORFs.csv', 'a') as gh:
            writer = csv.writer(gh)
	    #print(sORfsChrom.Chromosome[orf])
	    #print(mutationChrom['#CHROM'][mut])
            if(mutationChrom['#CHROM'][mut]==sORfsChrom.Chromosome[orf]):
                #print("chrome match")
		#print(float( sORfsChrom.Start_position[orf]))
		#print(float( mutationChrom.POS[mut]))
                #print(float( sORfsChrom.End_position[orf]))
                if(float( sORfsChrom.Start_position[orf]) <= float( mutationChrom.POS[mut]) <= float( sORfsChrom.End_position[orf])):
                    print(mutationChrom.iloc[mut].astype(str).tolist() + sORfsChrom.iloc[orf].astype(str).tolist())
                    writer.writerow(mutationChrom.iloc[mut].astype(str).tolist() + sORfsChrom.iloc[orf].astype(str).tolist())


        gh.close()

pool = multiprocessing.Pool(100)
#pool.map(partial(mappingMutation, sORfsChrom=sORFs, mutationChrom=vcfMutation),range(len(vcfMutation['#CHROM'])))
#Loop for chromosomes
for chrom in chromList:
    print(chrom)
    mutationChrom   = vcfMutation.loc[vcfMutation['#CHROM'] == chrom]
    sORfsChrom      = sORFs.loc[sORFs['Chromosome'] == chrom]
    sORfsChrom.index = range(len(sORfsChrom.index))
    #print( sORfsChrom)
    mutationChrom.index = range(len(mutationChrom.index))
    pool.map(partial(mappingMutation, sORfsChrom=sORfsChrom, mutationChrom=mutationChrom),range(len(mutationChrom['#CHROM'])))
    #print(sORfsChrom.index.tolist())
    #print(sORfsChrom.loc[61001,'Chromosome'])S
    #print(sORfsChrom.Chromosome )
    #for mut in range(len(mutationChrom['#CHROM'])):
    #print(mutationChrom.index.tolist())
    #zip(*pool.map(mappingMutation, iter(mutationChrom.index.tolist())))

    #pool = multiprocessing.Pool(1000)
    #mut = mutationChrom.index.tolist()
    #print(len(mut))
    #if len(mut)!=0:
    	#pool.map(partial(mappingMutation, sORfsChrom=sORfsChrom, mutationChrom=mutationChrom),mut)
    #Parallel(n_jobs=1000)(delayed(mappingMutation)(mut,sORfsChrom, mutationChrom ) for mut in iter(mutationChrom.index.tolist()))
	#close the pool and wait for the work to finish
    #pool.close()
    #pool.join()

    #print(chrom)


#result = sORFs[(sORFs.Start_position<=vcfMutation.POS )&(sORFs.End_position>=vcfMutation.POS)&(sORFs.Chromosome==vcfMutation['#CHROM'])]


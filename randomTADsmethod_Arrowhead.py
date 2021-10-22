import pandas as pd 
import random
import numpy
from pybedtools import BedTool
import pybedtools
import sys


dataset_name=sys.argv[1]
TADcaller=sys.argv[2]
bin_size=sys.argv[5]
#Set species to define correct chromsome lists
species=sys.argv[3]

#Upload TAD data and sorted gene lists
TADdata_string=str(str(dataset_name)+"_"+str(TADcaller)+str(bin_size)+"kb.txt")
TADdata = pd.read_csv(TADdata_string,sep='\t')
TADdata['length'] = TADdata['end'] - TADdata['start']

#Define chromsomes
if species=="Human":
	chromnumarg=sys.argv[4]
	if int(chromnumarg) <= 22:
		chrom=str("chr")+chromnumarg
	elif int(chromnumarg)==23:
		chrom=str("chrX")
	elif int(chromnumarg)==24:
		chrom=str("chrY")           
	chrom_size=pd.read_csv("hg19.chrom.sizes",sep='\t',header=None, names=["chr","size"])
	genes=pd.read_csv("hg19_proteincodinggenes_biomart_sorted.bed",sep='\t',header=None, names=["chr","start", "end","strand","gene"])

elif species=="Mouse":
	chromnumarg=sys.argv[4]
	if int(chromnumarg) <= 19:
		chrom=str("chr")+chromnumarg
	elif int(chromnumarg)==20:
		chrom=str("chrX")           
	elif int(chromnumarg)==21:
		chrom=str("chrY")

	chrom_size=pd.read_csv("mm10.chrom.sizes",sep='\t',header=None, names=["chr","size"])
	genes=pd.read_csv("mm10_proteincodinggenes_biomart_sorted.bed",sep='\t',header=None, names=["chr","start", "end","strand","gene"])
	

#Open outfile
outfile_string=str(str(dataset_name)+"_"+str(TADcaller)+"_"+str(chrom)+"_"+str(bin_size)+"kb_RANDOM.txt")
outfile=open(outfile_string,"w")
headings_outfile="TADID"+"\t"+"TAD_chr"+"\t"+"new_TAD_start"+"\t"+"new_TAD_end"+"\t"+"no_genes_in_TAD"+"\t"+"gene"+"\t"+"gene_chr"+"\t"+"gene_strand"+"\t"+"gene_start"+"\t"+"gene_end"+"\t"+"region_size"+"\n"
outfile.write(headings_outfile)


#Print the current chromosome
print(chrom)

#Extract data for current chromsome
chrTADdata=TADdata.loc[TADdata['chr'] == chrom]
chrgenes=genes.loc[genes['chr']==chrom]
chrchrom_size=chrom_size.loc[chrom_size['chr']==chrom]
chrchrom_size_val=int(chrchrom_size["size"])

#Shuffle TADs
if len(chrTADdata)>0:
	chrTADdata = chrTADdata.sample(frac=1)
else:
	print("No TADs on chromosome")

#Sort genes by their end position
chrgenes=chrgenes.sort_values(by=['start'])

#Intalise empty dataframe to keep track of new TADs
random6TADs=pd.DataFrame(columns=['chr', 'start', 'end'])


#Loop over TADs on current chromsome
for row in chrTADdata.itertuples():
	no_genes=row._2
	length=row.length
	TADID=row._1

	TADstart=row.start
	TADend=row.end

	#Calculate a random number between 1 and 10 to be used to decide whether current TAD can form non-nested overlaps (10% probability)
	tenpercent=random.randint(1,10)

	#Define variables to allow proposed new TAD position to be checked
	RanTADnogenes="NA"
	nonacceptedoverlap="NA"
	whilecounter=0
	broken=False

	#Test random positions for new TAD location
	while RanTADnogenes!=no_genes or nonacceptedoverlap!=False:
		whilecounter+=1

		#If new location can not be found in 10000 attempts break
		if whilecounter==10000:
			broken=True
			break

		#Generate a random number between 1 and the total length of current chrom minus the size of the current TAD
		randomstart=random.randint(1, chrchrom_size_val-length)
		randomend=randomstart+length

		#Make genes and random TAD into bed format
		newTAD=BedTool(str(chrom+"\t"+str(randomstart)+"\t"+str(randomend)), from_string=True)
		bedchrgenes=BedTool.from_dataframe(chrgenes)

		if len(random6TADs.index)!=0:

			#Checkoverlap with exisiting TADs
			bedrandom6TADs=BedTool.from_dataframe(random6TADs)
			checkoverlap=newTAD.intersect(bedrandom6TADs,wao=True)
			pandascheckoverlap=checkoverlap.to_dataframe(names=['newTADchr', 'newTADstart', 'newTADend','existingTADchr','existingTADstart','existingTADend','existingTADID','overlap'])
			pandascheckoverlap=pandascheckoverlap.loc[pandascheckoverlap['existingTADchr']!="."]

			if len(pandascheckoverlap.index)==0:
				nonacceptedoverlap=False
			else:
				overlapingTADs=pandascheckoverlap['existingTADID'].tolist()

				#Check if TAD is nested
				#New TAD inside existing TADs
				TADintersect_newinside=bedrandom6TADs.intersect(newTAD, wao=True, F=1)
				TADintersect_newinside=TADintersect_newinside.to_dataframe(names=['existingTADchr', 'existingTADstart', 'existingTADend','existingTADID', 'newTADchr', 'newTADstart', 'newTADend',  'overlap'])
				TADintersect_newinside=TADintersect_newinside.loc[TADintersect_newinside['newTADchr']!="."]
				TADintersect_newinsideTADs=TADintersect_newinside['existingTADID'].tolist()

				#Existing TAD inside new TAD
				TADintersect_insidenew=newTAD.intersect(bedrandom6TADs, wao=True, F=1)
				TADintersect_insidenew=TADintersect_insidenew.to_dataframe(names=['newTADchr', 'newTADstart', 'newTADend', 'existingTADchr', 'existingTADstart', 'existingTADend','existingTADID',  'overlap'])
				TADintersect_insidenew=TADintersect_insidenew.loc[TADintersect_insidenew['existingTADchr']!="."]
				TADintersect_insidenew_TADs=TADintersect_insidenew['existingTADID'].tolist()


				nestedTADs=set(TADintersect_newinsideTADs+TADintersect_insidenew_TADs)
				overlapingTADs=set(overlapingTADs)
				if nestedTADs==overlapingTADs:
					nonacceptedoverlap=False
				else:
					if tenpercent==1:
						nonacceptedoverlap=False
					else:
						nonacceptedoverlap=True
		else:
			nonacceptedoverlap=False


		if nonacceptedoverlap==False:
			#Full intersect with genes
			geneintersect=newTAD.intersect(bedchrgenes, wao=True, F=1)
			#Convert back to pandas dataframe
			pandasgeneintersect=geneintersect.to_dataframe(names=['TADchr', 'TADstart', 'TADend', 'genechrom', 'genestart', 'geneend', 'genestrand', 'gene', 'overlap'])
			#Filter out row created when no genes overlap
			pandasgeneintersect=pandasgeneintersect.loc[pandasgeneintersect['genechrom']!="."]

			RanTADnogenes=len(pandasgeneintersect.index)


		else:
			RanTADnogenes="NA"
			continue


	#Write random TAD to file and add accepted random TAD position to dataframe
	if broken==False:
		random6TADs=random6TADs.append({'chr': chrom, 'start': randomstart, 'end': randomend, 'id': TADID}, ignore_index=True)


		if RanTADnogenes==0:
			newline=str(TADID)+"\t"+str(chrom)+"\t"+str(randomstart)+"\t"+str(randomend)+"\t"+str(RanTADnogenes)+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str(int(randomend)-int(randomstart))+"\n"
			outfile.write(newline)		
			
		else:
			for generow in pandasgeneintersect.itertuples():
				newline=str(TADID)+"\t"+str(chrom)+"\t"+str(randomstart)+"\t"+str(randomend)+"\t"+str(RanTADnogenes)+"\t"+str(generow.gene)+"\t"+str(generow.genechrom)+"\t"+str(generow.genestrand)+"\t"+str(generow.genestart)+"\t"+str(generow.geneend)+"\t"+str(int(randomend)-int(randomstart))+"\n"
				outfile.write(newline)
	
	else:
		newline=str(TADID)+"\t"+str(chrom)+"\t"+str("NA")+"\t"+str("NA")+"\t"+str(no_genes)+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\t"+str("NA")+"\n"
		outfile.write(newline)

		
	pybedtools.cleanup(remove_all=True)

#snakemake met marijn en roel
rule all:
#uitlezen van id's uit "RNA-Seq-counts.txt"
rule read:
	input:
		"data/RNA-Seq-counts.txt"
	output:
		"data/geneids.txt"	
	run:
	
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			for i in readfile:
				z = i.split("\t",1)
				id = z[0]
				outputfile.write(id+"\n")
		

		
#Herkennen van het type gene id en converteren naar handzame gene identifiers
rule converter:
	input:
		"data/geneids.txt"
	output:
		"data/newgeneids.txt"
	run:
		
		from Bio import Entrez
		Entrez.email = 'marijnspr@gmail.com'
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			for line in readfile:
				lines = line.split('\n')
				id = lines[0]
				handle = Entrez.esearch(db="protein", term=id)
				record = Entrez.read(handle)
				handle.close()
				
				searchid = (record['IdList'])
				searchid = str(searchid[0])
				outputfile.write(searchid+'\n')
		


#ophalen van de sequentie
rule seqfinder:
	input:
		"data/newgeneids.txt"
	output:
		"data/sequenties.fasta"
	run:
		from Bio import Entrez
		Entrez.email = 'marijnspr@gmail.com'
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			for line in readfile:
				lines = line.split('\n')
				id = lines[0]
				handle = Entrez.efetch(db="nuccore", id=id, rettype="fasta", retmode="text")
				record = handle.read()
				handle.close()
				outputfile.write(record+'\n')
		
	
#ophalen van de functies
rule functionfinder:
	input:
		"data/newgeneids.txt"
	output:
		"data/functions.txt"
	run:
		from Bio import Entrez
		import re
		Entrez.email = 'marijnspr@gmail.com'
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			for line in readfile:
				lines = line.split('\n')
				id = lines[0]
				handle = Entrez.efetch(db="protein", id=id, rettype="ft", retmode="text")
				record = handle.read()
				handle.close()
				function = re.search('note.*',record).group(0)
				function = function.strip('note\t')
				outputfile.write(str(function)+'\n')

				
#ophalen van pubMedIDs 1/2
rule recordfinder:
	input:
		"data/sequenties.fasta",
		"data/newgeneids.txt"
	output:
		"data/records.txt"
	run:
		from Bio import Entrez
		Entrez.email = 'marijnspr@gmail.com'
		recordList = []
		idlist = []
		counter = 0
		with open(input[0], "r") as fasta, open(input[1], "r") as readfile, open(output[0], 'w') as outputfile:
			for line in readfile:
				lines = line.split('\n')
				id = lines[0]
				idlist.append(id)
			for line in fasta:
				if ">" in line:
					if "=" in line:
						lines = line.split("=")
						id = lines[1]
						handle = Entrez.esearch(db="pubmed", term=id, rettype="uilist", retmode="text")
						record = handle.read()
						handle.close()
					else:
						record = ""
					recordList.append(record)
					
			for i in idlist:
				outputfile.write("Gene: "+idlist[counter]+'\n'+recordList[counter]+'\n'+"-")
				counter+=1
		

# ophalen van pubMedIDs 2/2
rule pubmedfinder:
	input:
		"data/records.txt"
	output:
		"data/pubmedids.txt"
	run:
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			for line in readfile:
				lines = line.split('\n')
				if "<" in line:
					if "<Id>" in line:
						line = line.strip("<Id>")
						line = line[:-6]
						outputfile.write(line+'\n')
				else:
					outputfile.write(line)
	
	
#sortern op aantal pubMedIDs
rule pubmedsorter:
	input:
		"data/pubmedids.txt"
	output:
		"data/sortedpubmedids.txt"
	run:
		counter = 0
		with open(input[0], "r") as readfile, open(output[0], 'w') as outputfile:
			content = readfile.read()
			publist = content.split("-")
			publist.sort(key = len, reverse=True)
			for i in publist:	
				outputfile.write(publist[counter])
				counter+=1
	

#zoeken welk gen in welke pathway voorkomt aan de hand van KEGG_id
rule pathwayfinder:
	input:
		"data/geneids.txt"
	output:
		"data/pathways.txt"
	run:
		from bioservices import KEGG
		s = KEGG()
		idlist = []
		counter = 0
		
		with open(input[0], "r") as readfile, open(output[0],"w") as outputfile:
			for line in readfile:
				lines = line.split("\n")
				id = lines[0]
				
				idlist.append(id)
				pathway = s.get_pathway_by_gene(id,"lpl")
				if pathway == None:
					pathway = ""
				outputfile.write(idlist[counter]+'\t'+str(pathway)+'\n')
				counter+=1
			
rule resultaten:
	input:
		"data/geneids.txt",
		"data/newgeneids.txt",
		"data/sequenties.fasta",
		"data/functions.txt",
		"data/pubmedids.txt",
		"data/pathways.txt",
		"data/sortedpubmedids.txt"
	output:
		"results.txt"
	run:
		counter = 0
		counter2 = 0
		geneidlijst = []
		convertedgeneidlijst = []
		seqlijst = []
		functionlijst = []
		pubmedidlijst = []
		pathwayLijst = []
		sortedpubmedidlijst = []
		countpubmedidlijst = []
		
		with open(input[0], "r") as readfile1, open(input[1], "r") as readfile2, open(input[2], "r") as readfile3, open(input[3], "r") as readfile4, open(input[4], "r") as readfile5, open(input [5], "r") as readfile6, open(input[6], "r") as readfile7, open(output[0], 'w') as outputfile:
			outputfile.write("Results of the workflow"+'\n'+'\n')
			
			outputfile.write("Genen in volgorde van hoeveelheid voorkomen in artikelen (dalend van veel naar weinig):"+'\n')
			for line in readfile7:
				if "Gene:" in line:
					countpubmedidlijst = []
					lines = line.split('\n')
					sortedPubMedID = lines[0]
					sortedpubmedidlijst.append(sortedPubMedID)
					outputfile.write(sortedpubmedidlijst[counter2]+'\n')
					counter2+=1
			
			outputfile.write('\n'+"Rapport per Gen:"+'\n'+'\n')
			for line in readfile1:
				lines = line.split('\n')
				geneID = lines[0]
				geneidlijst.append(geneID)
			for line in readfile2:
				lines = line.split('\n')
				convertedGeneID = lines[0]
				convertedgeneidlijst.append(convertedGeneID)
			content = readfile3.read()
			seqlijst = content.split(">")
			for line in readfile4:
				lines = line.split('\n')
				function = lines[0]
				functionlijst.append(function)
			for line in readfile5:
				lines = line.split('\n')
				pubMedID = lines[0]
				pubmedidlijst.append(pubMedID)
			for line in readfile6:
				lines = line.split('\n')
				pathway = lines[0]
				pathwayLijst.append(pathway)
				
				outputfile.write("Gene ID:"+'\t'+geneidlijst[counter]+'\n'+"Converted Gene ID:"+'\t'+convertedgeneidlijst[counter]+'\n'+"Sequence:"+'\t'+seqlijst[counter][:-6]+'\n'+"Function:"+'\t'+functionlijst[counter]+'\n'+"PubMed entries:"+'\t'+pubmedidlijst[counter]+'\n'+"KEGG Pathways:"+'\t'+pathwayLijst[counter]+'\n'+'\n')
				counter+=1

rule visualize:
	output:
		"dag.pdf"
	shell:
		"snakemake --forceall --dag resultaten | dot -Tpdf > dag.pdf"
		
		

	
	

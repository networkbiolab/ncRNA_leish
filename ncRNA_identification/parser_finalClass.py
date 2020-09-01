import sys
name = sys.argv[1]
file = open(name,"r")
lista1=[]
for line in file:
		aux= line[0:-1].split("\t")
		aux2= aux[3].split(";")
		
		if "CD-box" in aux2 and "snRNA" in aux2:
			define="CD-box"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and "antisense" in aux2: 
			define=""
			if "miRNA" in aux2 and "sRNA" in aux2: 
				define="Unclassified"
			else:
				define="Cis-reg"
				print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and "miRNA" in aux2: 
			define=""
			if "piRNA" in aux2:
				define= "miRNA"
			else:
				if "snRNA" in aux2:
					define="Unclassified"
				else: 
					define="Cis-reg"
					print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and  "piRNA"  in aux2:
			define=""
			if "sRNA" in aux2: 
				define="Cis-reg"
			if "rRNA" in aux2: 
				define="rRNA"		 
			else:
				define= "Cis-reg"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" and "sRNA" in aux2: 
			define= "Cis-reg"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "HACA-box" in aux2 and "snRNA" in aux2: 
			define="HACA-box"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "piRNA" in aux2: 
			define=""
			if"rRNA" in aux2: 
				define="rRNA" 
			else: 
				define="miRNA"
				print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "rRNA" in aux2: 
			define="rRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "snRNA" in aux2: 
			define= "snRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "snoRNA" in aux2: 
			define="snoRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "piRNA" in aux2 and "snoRNA" in aux2:
			define="snoRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "snRNA" in aux2 and "snoRNA" in aux2: 
			define="snoRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Unclassified"in aux2 and "miRNA" in aux2:
			define="miRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Gene" in aux2 and "SRPRNA" in aux2:
			define="SRPRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "lncRNA" in aux2 and "miRNA" in aux2: 
			define="lncRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "lncRNA" in aux2 and "piRNA" in aux2: 
			define="lncRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "lncRNA" in aux2 and "siRNA" in aux2: 
			define="lncRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "siRNA" in aux2: 
			define: "miRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "piRNA,miRNA" in aux2:
			define="miRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "piRNA" in aux2 and "sisRNA" in aux2: 
			define="sisRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "rRNA" in aux2 and "siRNA" in aux2:
			define="rRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "rRNA" in aux2 and "snRNA" in aux2: 
			define="rRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif  "tRNA" in aux2 and "tmRNA" in aux2: 
			define="tRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Unclassified" in aux2 and "snoRNA"  in aux2: 
			define="snoRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and "rRNA" in aux2: 
			define="rRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and "Unclassified" in aux2: 
			define="Cis-reg"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "miRNA" in aux2 and "tRNA" in aux2: 
			define="tRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "RUF21" in aux2 and "miRNA" in aux2:
			define="RUF21"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "gRNA" in aux2 and "rRNA" in aux2: 
			define="gRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "SRPRNA" in aux2 and "piRNA" in aux2:
			define="SRPRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "SRPRNA" in aux2 and "Unclassified" in aux2: 
			define="SRPRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "lncRNA" in aux2 and "rRNA" in aux2:
			define="rRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "lncRNA" in aux2 and "snRNA" in aux2: 
			define="lncRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "piRNA" in aux2 and "snRNA" in aux2: 
			define="piRNA"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		elif "Cis-reg" in aux2 and "miRNA" in aux2 and "piRNA" in aux2 and "snRNA" in aux2: 
			define="Unclassified"
			print(aux[0]+"\t"+aux[1]+"\t"+aux[2]+"\t"+define)
		
		
		
		
		else:
			print(line[:-1])
		

import sys
name = sys.argv[1]
file = open(name,"r")
lista1=[]
for line in file:
		aux= line[0:-1].split("\t")
		
		if "mRNA" in aux[7]:
			print(line,end="")

import sys
name = sys.argv[1]
file = open(name,"r")
for i,line in enumerate(file):
	if i == 0: print(line, end="")
	else:
		line= line.replace("\n","")
		aux= line.split("\t")
		
		cont = 0
		for j in aux[1:]:
			if j <= "0.5": cont+=1
		if cont >= 1: 
			print(line+"\n", end="")

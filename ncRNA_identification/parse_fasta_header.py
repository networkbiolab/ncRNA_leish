import sys
name = sys.argv[1]
file = open(name,"r")
count=0
for line in file:
	
	if(line[0] == ">"):
		aux=line[0:-1].split(":")
		aux1=aux[0].replace(">","")
		aux2=aux[2].split("_")
		count+=1
			
		print(">"+aux2[0]+"|"+aux1+"|"+str(count))		
	else: 
		print(line,end="")

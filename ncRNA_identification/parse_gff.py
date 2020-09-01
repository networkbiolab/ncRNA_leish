import re 
import sys
name = sys.argv[1]
#file1 = open(name,"r")
name2 = sys.argv[2]

diccionario={}
f= open(name2,"r")


for line in f:
	aux = line[:-1].split("\t")
	diccionario[aux[0]]=aux[1]
	
f.close()
#print(diccionario)

with open(name, "r") as file1:
	for i in file1:
		for key in diccionario.keys():
			if key in i: 
				i = i.replace(key, diccionario[key])
				print(i,end="")




# use : python3 parse_gff.py file.gff dicctionary.file

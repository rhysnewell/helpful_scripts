#!/usr/bin/env python3
import sys

def main():
	try:
	    faa = sys.argv[1]
	    fna = sys.argv[2]
	except IndexError:
	    print("Usage: <Fasta Amino Acids> <Fasta Nucleotides> \n", 
	    "Fasta Amino Acids: Output from EnrichM clutering annotations \n",
	    "Fasta Nucleotides: Output from prodigal")
	    return
	
	for line in open(fna):
	    if line.startswith('>'):
	        dna_id = line.split(' ')[0].replace('>', '')
	        for prot in open(faa):
	            if prot.startswith('>'):
	                prot_id = prot.split(' ')[0].replace('>', '')
	                if dna_id == prot_id:
	                    print(prot.strip())
	    else:
	        print(line.strip())

if __name__ == "__main__":
	main()
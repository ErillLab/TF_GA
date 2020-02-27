# Main execution

from Bio import SeqIO
import time
from objects.organismFactory import OrganismFactory
import random


POPULATION_LENGTH = 10
POSITIVE_FILENAME = '61383_ref_mLynCan4_v1.p_chrA1.fa'
organismPopulation = []
positiveDataset = []
negativeDataset = []
threshold = 5


def main():
    
    # TEST READ  DATASET
    # mySeq = readFastaFile(POSITIVE_FILENAME)
    # print(len(mySeq))

    # Generate initial population
    organismFactory = OrganismFactory()

    for i in range(POPULATION_LENGTH):
        organismPopulation.append(organismFactory.getOrganism())
    print("Total population: "+str(len(organismPopulation)))
    
    algo = 10

    

    while algo > threshold:
        
        # Shuffle population
        
        random.shuffle(organismPopulation)
        for i in range(0, len(organismPopulation) - 1, 2):
            print(i)
            org1 = organismPopulation[i]
            org2 = organismPopulation[i+1]
            id1 = org1.getID()            
            id2 = org2.getID()            
            print("Pareja: {} {}".format(id1, id2))
            
        algo -= 2


def combineOrganisms(organism1, organism2):
    






    
# Reads a fasta files and return an array of DNA sequences
def readFastaFile(filename):
    dataset = []
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')    

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    return dataset


if __name__ == '__main__':
    initial = time.time()
    main()
    print("Time: {}".format(time.time()-initial))




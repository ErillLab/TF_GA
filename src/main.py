# Main execution

from Bio import SeqIO
import time
from objects.organismFactory import OrganismFactory
import random
import copy

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
    # Main loop, it iterates untill organisms do not get a significant change.
    while algo > threshold:
        
        # Shuffle population
        random.shuffle(organismPopulation)

        # Iterate over pairs of organisms
        for i in range(0, len(organismPopulation) - 1, 2):
            print(i)
            org1 = organismPopulation[i]
            org2 = organismPopulation[i+1]
            id1 = org1.getID()            
            id2 = org2.getID()            
            print("Pareja: {} {}".format(id1, id2))
            children = combineOrganisms(org1, org2, organismFactory)
            org1.print()
            org2.print()
            for child in children:
                child["child"].print()

            print("---- objects ---")
            print(children)

            break
        # Update organisms changes
        break
        algo -= 4

# Gets 2 organisms, and returns 2 children with format (child, similarity to parent 1, similarity to parent 2)
def combineOrganisms(organism1, organism2, organismFactory):
    # Save the number of nodes from the parents
    nNodesOrg1 = organism1.countNodes()
    nNodesOrg2 = organism2.countNodes()
    
    # Create the 2 childs and assign new IDs
    child1 = copy.deepcopy(organism1)
    child2 = copy.deepcopy(organism2)

    child1.setID(organismFactory.getID())
    child2.setID(organismFactory.getID())

    # Select random nodes from each child
    randomNode1 = random.randint(0,nNodesOrg1 - 1)
    randomNode2 = random.randint(0,nNodesOrg2 - 1)
    node1 = child1.getNode(randomNode1)
    node2 = child2.getNode(randomNode2)

    # Save the number of nodes taken from each  child
    nNodesFromOrg1 = node1.countNodes()
    nNodesFromOrg2 = node2.countNodes()

    # parentNodeX has an specific scructure:
    # parentNodeX[0] = isRootNode
    # parentNodeX[1] = parentNode
    # parentNodeX[2] = isLeftSideNode
    # 
    # Search is done by ID, so this should be checked bc 
    # you could find duplicated IDs!!

    parentNode1 = child1.getParent(node1.ID)
    parentNode2 = child2.getParent(node2.ID)
    
    #Set nodes in 
    if parentNode1[0]:
        # Its the root node of child 1
        child1.setRootNode(node2)
    else:
        if parentNode1[2]:
            # Child on left side
            parentNode1[1].setNode1(node2)
        else:
            # Child on left side
            parentNode1[1].setNode2(node2)

    if parentNode2[0]:
        # Its the root node of child 2
        child2.setRootNode(node1)
    else: 
        if parentNode2[2]:
            # Child on left side
            parentNode2[1].setNode1(node1)
        else:
            # Child on left side
            parentNode2[1].setNode2(node1)


    # Swap nodes
    #print("NODE1.ID = {}\nNODE2.ID = {}".format(node1.ID, node2.ID)) 
    #child1.setNode(node2, node1.ID)
    #print("NODE1.ID = {}\nNODE2.ID = {}".format(node1.ID, node2.ID)) 
    #child2.setNode(node1, node2.ID)
    #print("NODE1.ID = {}\nNODE2.ID = {}".format(node1.ID, node2.ID)) 
    
    #tmp = node1
    #node1 = node2
    #node2 = tmp
    

    nNodesChild1 = child1.countNodes()
    nNodesChild2 = child2.countNodes()
    
    # dictionary with an organism and similarities to each parent
    child1Similarities = {
            "simOrg1": (nNodesChild1 - nNodesFromOrg2) / nNodesChild1,
            "simOrg2": nNodesFromOrg2 / nNodesChild1,
            "child": child1
            }
    

    child2Similarities = {
            "simOrg1": nNodesFromOrg1 / nNodesChild2,
            "simOrg2": (nNodesChild2 - nNodesFromOrg1) / nNodesChild2,
            "child": child2
            }
    return [child1Similarities, child2Similarities]


    
# Reads a fasta files and return an array of DNA sequences (strings)
def readFastaFile(filename):
    dataset = []
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')    

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    return dataset

# Entry point to app execution
# It calculates the time, but could include other app stats
if __name__ == '__main__':
    initial = time.time()
    main()
    print("\n")
    print("-"*50)
    print(" "*20+"STATS")
    print("-"*50)
    print("Time: {}".format(time.time()-initial))
    print("-"*50)




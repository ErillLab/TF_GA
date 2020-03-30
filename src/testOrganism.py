# Tests an organism Fitness

from searchOrganisms import readFastaFile, readJsonFile, exportOrganism
from objects.organismFactory import OrganismFactory

import random
CONFIG_FILE = "config.json"

def main():


    config = readJsonFile(CONFIG_FILE)
    positivePath = config["main"]["DATASET_BASE_PATH_DIR"] + config["main"]["POSITIVE_FILENAME"]
    negativePath = config["main"]["DATASET_BASE_PATH_DIR"] + config["main"]["NEGATIVE_FILENAME"]
    inputOrganismsPath = "inputOrganisms.json"

    positiveDataset = readFastaFile(positivePath)
    negativeDataset = readFastaFile(negativePath)
    print("{} {}".format(len(positiveDataset), len(negativeDataset)))

    organismFactory = OrganismFactory(config["organism"], config["organismFactory"], config["connector"], config["pssm"])

    aOrganisms = organismFactory.importOrganisms("inputOrganisms.json")
    random.shuffle(negativeDataset)

    for org in aOrganisms:

        org.print()
        
        p1 = org.getScore(positiveDataset)
        n1 = org.getScore(negativeDataset[:50])
        print("ORG {} Positive: {} Negative: {}\n".format(org.ID, p1, n1))

        exportOrganism(org, positiveDataset, "{}positive_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID))
        exportOrganism(org, negativeDataset[:50], "{}negative_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID))
if __name__ == '__main__':
    main()

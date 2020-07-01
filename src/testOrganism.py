# Tests an organism Fitness

from searchOrganisms import readFastaFile, readJsonFile, exportOrganism
from objects.organismFactory import OrganismFactory

CONFIG_FILE = "config.json"


def main():

    config = readJsonFile(CONFIG_FILE)
    positivePath = (
        config["main"]["DATASET_BASE_PATH_DIR"] + config["main"]["POSITIVE_FILENAME"]
    )
    negativePath = (
        config["main"]["DATASET_BASE_PATH_DIR"] + config["main"]["NEGATIVE_FILENAME"]
    )
    complexityFactor = config["main"]["COMPLEXITY_FACTOR"]
    MAX_SEQUENCES_TO_FIT_POS = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    MAX_SEQUENCES_TO_FIT_NEG = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]

    inputOrganismsPath = config["main"]["INPUT_FILENAME"]
    mean_nodes = 3.0
    mean_fitness = 150
    positiveDataset = readFastaFile(positivePath)
    positiveDataset.sort()
    negativeDataset = readFastaFile(negativePath)
    print("{} {}".format(len(positiveDataset), len(negativeDataset)))

    organismFactory = OrganismFactory(
        config["organism"],
        config["organismFactory"],
        config["connector"],
        config["pssm"],
    )

    aOrganisms = organismFactory.importOrganisms(inputOrganismsPath)
    # random.shuffle(negativeDataset)

    for org in aOrganisms:

        # org.print()
        nodes = org.countNodes()

        p1 = org.getSeqSetFitness(positiveDataset[:MAX_SEQUENCES_TO_FIT_POS])
        n1 = org.getSeqSetFitness(negativeDataset[:MAX_SEQUENCES_TO_FIT_NEG])
        # p1 = 20
        # n1 = org.getSeqSetFitness(negativeDataset[31:32])
        c1 = org.getComplexity(mean_nodes, mean_fitness)
        
        # Score
        fitness = p1 - n1
        effectiveFitness = fitness - complexityFactor * c1
        print(
            "ORG {} N: {:.2f} P: {:.2f} N: {:.2f} C: {:.2f} F: {:.2f} EF: {:.2f}\n".format(
                org.ID, nodes, p1, n1, c1, fitness, effectiveFitness
            )
        )

        exportOrganism(
            org,
            positiveDataset,
            "{}positive_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID),
            organismFactory,
        )
        # exportOrganism(
        #     org,
        #     negativeDataset[31:32],
        #     "{}negative_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID),
        #     organismFactory,
        # )
        
        exportOrganism(
            org,
            negativeDataset[:50],
            "{}negative_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID),
            organismFactory,
        )


if __name__ == "__main__":
    main()

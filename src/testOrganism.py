"""Tests an organism Fitness
"""

from searchOrganisms import read_fasta_file, read_json_file, export_organism
from objects.organismFactory import OrganismFactory

CONFIG_FILE = "config.json"


def main():
    """Main execution for the test organisms

    """

    config = read_json_file(CONFIG_FILE)
    posititve_path = (
        config["main"]["DATASET_BASE_PATH_DIR"]
        + config["main"]["POSITIVE_FILENAME"]
    )
    negative_path = (
        config["main"]["DATASET_BASE_PATH_DIR"]
        + config["main"]["NEGATIVE_FILENAME"]
    )
    complexity_factor = config["main"]["COMPLEXITY_FACTOR"]
    max_sequences_to_fit_pos = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    max_sequences_to_fit_neg = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]

    input_organisms_path = config["main"]["INPUT_FILENAME"]
    mean_nodes = 3.0
    mean_fitness = 150
    positive_dataset = read_fasta_file(posititve_path)
    positive_dataset.sort()
    negative_dataset = read_fasta_file(negative_path)
    print("{} {}".format(len(positive_dataset), len(negative_dataset)))

    organism_factory = OrganismFactory(
        config["organism"],
        config["organismFactory"],
        config["connector"],
        config["pssm"],
    )

    a_organisms = organism_factory.import_organisms(input_organisms_path)
    # random.shuffle(negativeDataset)

    for org in a_organisms:

        # org.print()
        nodes = org.count_nodes()

        p_1 = org.get_seq_set_fitness(
            positive_dataset[:max_sequences_to_fit_pos]
        )
        n_1 = org.get_seq_set_fitness(
            negative_dataset[:max_sequences_to_fit_neg]
        )
        # p1 = 20
        # n1 = org.getSeqSetFitness(negativeDataset[31:32])
        c_1 = org.get_complexity(mean_nodes, mean_fitness)

        # Score
        fitness = p_1 - n_1
        effective_fitness = fitness - complexity_factor * c_1
        print(
            (
                "ORG {} N: {:.2f} P: {:.2f} N: {:.2f} C: {:.2f} F: {:.2f}"
                + " EF: {:.2f}\n"
            ).format(org._id, nodes, p_1, n_1, c_1, fitness, effective_fitness)
        )

        export_organism(
            org,
            positive_dataset,
            "{}positive_{}".format(
                config["main"]["RESULT_TEST_BASE_PATH_DIR"], org._id
            ),
            organism_factory,
        )
        # exportOrganism(
        #     org,
        #     negativeDataset[31:32],
        #     "{}negative_{}".format(config["main"]["RESULT_TEST_BASE_PATH_DIR"], org.ID),
        #     organismFactory,
        # )

        export_organism(
            org,
            negative_dataset[:50],
            "{}negative_{}".format(
                config["main"]["RESULT_TEST_BASE_PATH_DIR"], org._id
            ),
            organism_factory,
        )


if __name__ == "__main__":

    # TODO: Add the profiler here to improve the fitness calculation performance
    main()

"""Tests an organism Fitness
"""

from search_organisms import read_fasta_file, read_json_file, export_organism
from objects.organism_factory import OrganismFactory

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
    max_sequences_to_fit_pos = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    max_sequences_to_fit_neg = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]

    input_organisms_path = config["main"]["INPUT_FILENAME"]
    positive_dataset = read_fasta_file(posititve_path)
    positive_dataset.sort()
    negative_dataset = read_fasta_file(negative_path)
    #print("{} {}".format(len(positive_dataset), len(negative_dataset)))
    
    genome_length = config["main"]["GENOME_LENGTH"]
    
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
        
        # Boltzmannian fitness
        performance1 = org.get_boltz_fitness(positive_dataset[:max_sequences_to_fit_pos],
                                             negative_dataset[:max_sequences_to_fit_neg],
                                             genome_length)
        boltz_fitness = performance1["score"]
        
        # Gini coefficient
        gini_coeff = performance1["avg_gini"]
        
        # Discriminative fitness
        P = org.get_discriminative_fitness(positive_dataset[:max_sequences_to_fit_pos])["score"]
        
        N = org.get_discriminative_fitness(negative_dataset[:max_sequences_to_fit_neg])["score"]
        
        discr_fitness =  P - N
        
        
        print(
            (
                "Org {} Nodes: {:.2f} GiniPSSMs: {:.2f} P: {:.2f} N: {:.2f}"
                + " DiscrF: {:.2f} BoltzF: {:.2f}\n"
            ).format(
                org._id,  # "Org"
                nodes,  # "Nodes"
                gini_coeff, # GiniPSSMs
                P,  # "P"
                N,  # "N"
                discr_fitness,  # "DiscrF"
                boltz_fitness,  # "BoltzF"
                )
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


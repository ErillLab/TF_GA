"""
Oragnism object
It allocates the full data structure
"""

import random
import numpy as np


class OrganismObject:
    """Organism object
    """

    def __init__(self, _id: int, conf: dict, root_node=None) -> None:
        """Organism constructor

        Args:
            _id: Organism identifier
            conf: Configuration from JSON file
            root_node (Node): Top-level node for the tree data structure
        """
        self._id = _id
        self.root_node = root_node
        self.num_nodes = 0

        self.cumulative_fit_method = conf["CUMULATIVE_FIT_METHOD"]
        self.mutate_probability_substitute_pssm = conf[
            "MUTATE_PROBABILITY_SUBSTITUTE_PSSM"
        ]
        self.mutate_probability_rise_child = conf[
            "MUTATE_PROBABILITY_RISE_CHILD"
        ]
        self.mutate_probability_sunk_child = conf[
            "MUTATE_PROBABILITY_SUNK_CHILD"
        ]
        self.mutate_probability_node_mutation = conf[
            "MUTATE_PROBABILITY_NODE_MUTATION"
        ]
        self.min_nodes = conf["MIN_NODES"]
        self.max_nodes = conf["MAX_NODES"]
        self.is_tracked = False

    # Setters an getters
    def set_root_node(self, root_node) -> None:
        """Setter root_node

        Args:
            root_node (Node): Top-level node for the tree data structure
        """
        self.root_node = root_node

    def get_id(self) -> int:
        """Getter _id

        Returns:
            _id of the organism
        """
        return self._id

    def set_id(self, _id: int) -> None:
        """Setter _id

        Args:
            _id: ID to to set in the organism
        """
        self._id = _id

    def set_is_tracked(self, new_tracked: bool):
        """Setter is_tracked

        Args:
            new_tracked: True if the organism should be tracked.
                         False otherwise.
        """
        self.is_tracked = new_tracked

    def mutate(self, org_factory) -> None:
        """Mutates an organism based on JSON configured probabilities

        Args:
            org_factory (OrganismFactory): Factory of organisms and node
                                           components
        """

        # Substitute a random node by a random PSSM
        if random.random() < self.mutate_probability_substitute_pssm:

            n_nodes = self.count_nodes()
            random_node = random.randint(0, n_nodes - 1)
            substituted_node = self.get_node(random_node)
            parent_node = self.get_parent(substituted_node._id)
            # TODO: Set the length in the PSSM config
            new_node = org_factory.createPSSM(org_factory.PWM_LENGTH)

            if parent_node["is_root_node"]:

                self.root_node = new_node
            else:

                if parent_node["is_left_side"]:

                    parent_node["self"].set_node1(new_node)
                else:

                    parent_node["self"].set_node2(new_node)

            self.reset_ids()

        # Rise the level of a child
        if random.random() < self.mutate_probability_rise_child:

            n_nodes = self.count_nodes()
            random_node = random.randint(0, n_nodes - 1)
            rised_node = self.get_node(random_node)
            parent_node = self.get_parent(rised_node._id)

            #
            if not rised_node._id == self.root_node._id:

                parent1 = self.get_parent(rised_node._id)
                parent2 = self.get_parent(parent1["self"]._id)

                if parent2["is_root_node"]:
                    self.root_node = rised_node
                else:
                    if parent2["is_left_side"]:
                        parent2["self"].set_node1(rised_node)
                    else:
                        parent2["self"].set_node2(rised_node)
                self.reset_ids()

        # Add complexity to the organism
        if random.random() < self.mutate_probability_sunk_child:
            n_nodes = self.count_nodes()
            random_node = random.randint(0, n_nodes - 1)
            suken_node = self.get_node(random_node)
            parent_node = self.get_parent(suken_node._id)

            new_node = org_factory.create_connection(0)

            if parent_node["is_root_node"]:
                self.root_node = new_node
            else:
                if parent_node["is_left_side"]:

                    parent_node["self"].set_node1(new_node)
                else:
                    parent_node["self"].set_node2(new_node)

            if random.random() < 0.5:
                new_node.set_node1(suken_node)
            else:
                new_node.set_node2(suken_node)
            self.reset_ids()

        # Mutate a random node
        if random.random() < self.mutate_probability_node_mutation:

            n_nodes = self.count_nodes()
            random_node = random.randint(0, n_nodes - 1)
            mutated_node = self.get_node(random_node)
            mutated_node.mutate(org_factory)

    def get_complexity(self, mean_nodes: float, mean_fitness: float) -> float:
        """Returns the implicit complexity assiciated to the  current organism

        Args:
            mean_nodes: Average number of nodes of the population
            mean_fitness: Average fitness of the population

        Returns:
            Complexity penalty assiciated to the organism

        """
        # Complexity is calculed as:
        # meanFitnessScore * # nodes / mean_nodes

        # Check complexity of the organism
        # If its over/under organism MAX/MIN apply an extra complexity penalty
        extra_penalty = 0.0
        extra_penalty_factor = 300
        base_penalty = 0.0
        nodes = self.count_nodes()

        base_penalty = mean_fitness * self.num_nodes / mean_nodes

        # introduce bound for number of nodes
        # organisms containing more nodes get an additional penalty factor
        if nodes < self.min_nodes:
            extra_penalty = (self.min_nodes - nodes) * extra_penalty_factor
        if nodes > self.max_nodes:
            extra_penalty = (nodes - self.max_nodes) * extra_penalty_factor

        return base_penalty + extra_penalty

    def get_seq_fitness(self, s_dna: str) -> dict:
        """Return the fitness of the organism for a given DNA sequence

        Args:
            s_dna: DNA sequence to analize

        Returns:
           score, blocked and blockers
        """

        # call recursively to get the total fitness of the organism
        node_root = self.root_node.get_placement_2(s_dna, len(s_dna))

        if len(node_root) < 1:
            print("Too few placement options")
            return {
                "energy": -1000,
                "position": 0,
                "lock_vector": []
                    }

        # return score, blocks and blokcers in that sequence
        return node_root[0]

    def get_seq_set_fitness(self, a_dna: list) -> float:
        """Return the total Fitness for an array of DNA sequences and the
        fitness method

        Args:
            a_dna: list of dna sequences

        Returns:
            score assigned to the organism
        """

        score = 0

        # sum method returns the sum of all fitness to DNA sequences
        if self.cumulative_fit_method == "sum":

            for s_dna in a_dna:
                sfit = self.get_seq_fitness(s_dna)
                score += sfit["energy"]

        # mean method returns the mean of all fitness to SNA sequences
        elif self.cumulative_fit_method == "mean":

            scores = []
            for s_dna in a_dna:
                sfit = self.get_seq_fitness(s_dna)
                scores.append(sfit["energy"])
            score = np.mean(scores)

        return score

    def get_node(self, objective: int):
        """Returns a node Number N based on in-order search. [0-N)

        Args:
            objective: _id of the objective

        Returns:
            Node with specified _id if found. None otherwise
        """
        node_count = 0
        return self.root_node.get_node(objective, node_count)

    def set_node(self, node, _id: int) -> None:
        """Set the node in the _id node position

        Args:
            node (Node): Node to set in the tree
            _id: ID of the position to set the node
        """

        print("node._id = {} ID to change {}".format(node._id, _id))
        if self.root_node._id == _id:
            self.root_node = node
        else:
            self.root_node.setNode(node, _id)

    def get_parent(self, _id: int) -> dict:
        """Get the parent node of a given _id and if it is the left child

        Args:
            _id: ID of the child

        Returns:
            dictionary with the keys:
            "is_root_node": True if its root of the organism. False ortherwise.
            "self":  Parent node of the _id child
            "is_left_side": True if ID is from left side of the parent
        """

        if self.root_node._id == _id:
            return {"is_root_node": True}

        return self.root_node.get_parent(_id)

    def count_nodes(self) -> int:
        """Returns the number of nodes of the organism

        Returns:
            Number of nodes of hte organism
        """

        self.num_nodes = self.root_node.count_nodes()

        return self.num_nodes

    def reset_ids(self) -> None:
        """Reset _ids of the full organism
        """
        first_id = 0
        self.root_node.reset_id(first_id)

    def print(self) -> None:
        """Prints the whole tree data structure
        """
        print()
        print("***** Organism {} *****".format(self._id))
        self.root_node.print(1)

    def export(self, filename: str) -> None:
        """Exports the whole tree data structure

        Args:
            filename: Name of the file to export the organism
        """
        organism_file = open(filename, "w+")
        organism_file.write("***** Organism {} *****".format(self._id))

        self.root_node.export(organism_file, 0)
        organism_file.write("\n")
        organism_file.close()

    def export_results(self, a_dna: list, filename: str) -> None:
        """Exports all DNA sequences organism binding to a file

        Args:
            filename: Name of the file to export sequences
            a_dna: list fo sequences to export

        """

        # Sort the array, so its always shown in the same order
        # sorting is done by sequence, so first sequences start with "AAA.."
        a_dna.sort()

        results_file = open(filename, "w+")

        # for every DNA sequence
        for s_dna in a_dna:

            # call fitness evaluation for sequence
            sfit = self.get_seq_fitness(s_dna.lower())

            # write out the sequence
            results_file.write("\n{}\n".format(s_dna))

            # create an empy positions map
            map_positions = "-" * len(s_dna)

            # positions for PSSMs are in blocks and blocked lists, returned by
            # getSeqFitness. we zip them and then iterate over the zip to
            # print the PSSMs in their locations respective to the sequence
            positions = sfit["lock_vector"]
            for pssm in positions:
                # print _id, capped to the length of PSSM
                _p = round(pssm["position"])
                pssm_str = (str(pssm["id"]) * pssm["length"])[:pssm["length"]]

                # fill up map at correct positions
                map_positions = (
                    map_positions[:_p] + pssm_str + map_positions[_p + pssm["length"]:]
                )

            # write map to file for this sequence
            results_file.write(map_positions + "\n")
            # resultsFile.write(str(stuff) + "\n")

        results_file.close()

    def print_result(self, s_dna: str) -> str:
        """Prints the results of s_dna binding sites to stdout

        Args:
            s_dna: DNA sequence to export

        Returns:
            DNA sequence and binding sites of the organisms recognizer
        """

        s_dna = s_dna.lower()

        # call fitness evaluation for sequence
        sfit = self.get_seq_fitness(s_dna.lower())

        # create an empy positions map
        map_positions = "-" * len(s_dna)
        
        # positions for PSSMs are in blocked and blocked lists, returned by
        # getSeqFitness. we zip them and then iterate over the zip to
        # print the PSSMs in their locations respective to the sequence
        positions = sfit["lock_vector"]
        for pssm in positions:
            # print _id, capped to the length of PSSM
            _p = round(pssm["position"])
            pssm_str = (str(pssm["id"]) * pssm["length"])[:pssm["length"]]

            # fill up map at correct positions
            map_positions = (
                map_positions[:_p] + pssm_str + map_positions[_p + pssm["length"]:]
            )
            # handle two-digit positions, by alterning between digits

        # return map for this sequence
        return "{}\n{}".format(s_dna, map_positions)

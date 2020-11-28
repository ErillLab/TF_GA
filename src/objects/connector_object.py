"""C object
Connects two nodes at a specific distance

"""

# pylint: disable=E0402
# type: ignore
import random
from .node_object import Node
import numpy as np
import math


def norm_pdf(x, mu, sigma):
    if sigma != 0:
        var = float(sigma)**2
        denom = (2*math.pi*var)**.5
        num = math.exp(-(float(x)-float(mu))**2/(2*var))
        p = num/denom
    else:
        # when sigma is 0
        if x == mu:
            p = 1
        else:
            p = 0
    return p


def norm_cdf(x, mu, sigma):
    # Cumulative distribution function for the normal distribution
    z = (x-mu)/abs(sigma)
    return (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0



# pylint: disable=R0902
class ConnectorObject(Node):
    """Connector Object is a node that connects two nodes
    """

    # pylint: disable=R0913
    def __init__(
            self,
            _mu: int,
            _sigma: int,
            config: dict,
            node_1: Node = Node(),
            node_2: Node = Node(),
    ):
        """Connector constructor gets mu, sigma and can get the two initial nodes.

        Args:
            _mu: Mean distance between node1 and node2
            _sigma: Variance between node 1 and node2
            config: Configurations loadad from config.json
            node_1: Conceptual left node
            node_2: Conceptual right node

        """
        super().__init__()
        self._mu = _mu  # Mean discance
        self._sigma = _sigma  # Variance between elements
        self._id = 0
        self.mutate_probability_sigma = config["MUTATE_PROBABILITY_SIGMA"]
        self.mutate_probability_mu = config["MUTATE_PROBABILITY_MU"]
        self.mutate_probability_swap = config["MUTATE_PROBABILITY_SWAP"]
        self.tau = config["TAU"]
        self.mutate_variance_sigma = config["MUTATE_VARIANCE_SIGMA"]
        self.mutate_variance_mu = config["MUTATE_VARIANCE_MU"]
        self.placement_options = config["PLACEMENT_OPTIONS"]

        self.node1 = node_1
        self.node2 = node_2

    # pylint: enable=R0913
    # Setters
    def set_mu(self, _mu: int) -> None:
        """Set mu variable

        Args:
            _mu: Mean distance between node1 and node2
        """
        self._mu = _mu

    def set_sigma(self, sigma: int) -> None:
        """Set sigma variable

        Args:
            sigma: Variance between node 1 and node2
        """
        self._sigma = sigma

    def set_node1(self, node1: Node) -> None:
        """Set left node

        Args:
            node1: Conceptual left node
        """
        self.node1 = node1

    def set_node2(self, node2: Node) -> None:
        """Set right node

        Args:
            node2: Conceptual right node
        """
        self.node2 = node2

    def count_nodes(self) -> int:
        """Counts the number of nodes below the node (including itself)

        Returns:
            Number of nodes below the current connector
        """
        return self.node1.count_nodes() + self.node2.count_nodes() + 1

    def get_node(self, objective: int, node_count: int) -> Node:
        """Get a specific node based on a count and the objective node

        Args:
            objective: ID of the objectivo node to return
            node_count: Number of nodes analized

        Returns:
            Node with the objective. None otherwise
        """
        left_nodes = self.node1.count_nodes()
        node_number = left_nodes + node_count
        returned_node = None

        if node_number == objective:
            # We are on the node we are searching
            returned_node = self
        elif node_number < objective:
            # The node is on the right side
            returned_node = self.node2.get_node(objective, node_number + 1)
        elif node_number > objective:
            # The node is on the left side
            returned_node = self.node1.get_node(objective, node_count)

        return returned_node

    def get_parent(self, _id: int) -> dict:
        """Get the parent node of a given ID and if it is the left child

        Args:
            _id: Identificator of the child node

        Returns:
            dictionary with the parent info if exists:
            "isRootNode": True if the _id if from the root node and has.
                          no parent
            "self":       Node object with the parent of the _id child.
            "isLeftSide": True if _id child is on the left side.
                          False otherwise.

            if the parent does not exist, returns None
        """

        if self.node1._id == _id:
            # Return itself and specify child is on left side
            return {"is_root_node": False, "self": self, "is_left_side": True}

        if self.node2._id == _id:
            # Return itself and specify child is on right side
            return {"is_root_node": False, "self": self, "is_left_side": False}

        check_node = None
        check_node = self.node1.get_parent(_id)
        if check_node is None:
            check_node = self.node2.get_parent(_id)
        return check_node

    def get_all_pssm(self) -> list:
        """returns an array of all pssm objects of the organism

        Returns:
            List with all the pssm
        """
        return self.node1.getAllPssm() + self.node2.getAllPssm()


    # pylint: enable=R1702
    # pylint: enable=R0915
    def get_placement(
                self,
            s_dna: str,
            s_dna_len: int,
            automatic_placement_options: int,
            is_automatic_placement_options: bool
    ) -> list:
        """Compute the best option to connect its nodes.

        Args:
            s_dna: DNA sequence
            s_dna_len: length of the DNA sequence
            automatic_placement_options: Computed automatic placement options
            is_automatic_placement_options: true if automatic placement options should be plaed

        Returns:
            list of the best placed nodes with this connector

        Description:
            This function implements the placement behavior for connectors.
            The placement problem is defined as who to best position an
            organism on a sequence (i.e. how to maximize its fitness given
            the sequence).
            The implementation in this function follows the recursive 
            formulation of the organism. 
            
            Connectors receive from lower nodes (either connectors or PSSMs)
            a ranked, truncated list of possible placements, which includes the
            positions that have been blocked for PSSMs under those placements.
            
            The connector then iterates over all possible combinations of
            these placements, computing the overall energy (i.e. including its
            daugther nodes energies and its own energy term), and returns a
            ranked list of placements (and their blocked positions) to the 
            upper level connector.
        """

        possible_candidates = []

        possibilities_node_1 = self.node1.get_placement(
            s_dna,
            s_dna_len,
            automatic_placement_options,
            is_automatic_placement_options
                )
        possibilities_node_2 = self.node2.get_placement(
            s_dna,
            s_dna_len,
            automatic_placement_options,
            is_automatic_placement_options
                )
        

        # for all possible daughter placement combinations
        for possibility_1 in possibilities_node_1:
            for possibility_2 in possibilities_node_2:

                # Check that the placement of PSSMs for one daugther (i.e. its
                # lock vector) does not conflict with the other one
                is_overlapping = False
                for pssm_1 in possibility_1["lock_vector"]:
                    # if is_overlapping:
                    #     continue
                    for pssm_2 in possibility_2["lock_vector"]:
                        # if is_overlapping:
                        #     continue

                        # It can be done in one sentence but in multiple lines
                        # it looks more understandable.
                        # check for daughter node conflict:
                        # check for overlap of the form:
                        #   2222
                        # 1111
                        if (pssm_2["position"] < pssm_1["position"] + pssm_1["length"]
                                and pssm_2["position"] >= pssm_1["position"]):
                            is_overlapping = True

                        # check for overlap of the form:
                        # 2222
                        #  1111
                        if (pssm_1["position"] < pssm_2["position"] + pssm_2["length"]
                                and pssm_1["position"] >= pssm_2["position"]):
                            is_overlapping = True

                # if there is overlap, disregard combination, go for next
                # daughter node 2 placement option
                if is_overlapping:
                    continue

                # if there is no overlap, compute the overall energy of the
                # arrangement and add it to the list of possible placements
                
                # Distance d
                d = possibility_2["position"] - possibility_1["position"]
                
                # Numerator                
                numerator = norm_pdf(d, self._mu, self._sigma)
                
                # Normalize by AUC within the range of observable d values
                
                max_d = s_dna_len - 1  # Maximum d observable
                min_d = -1 * max_d  # Minimum d observable
                
                if self._sigma == 0:
                    auc = 1.0  # all the gaussian is within the [-(L-1), +(L-1)] range
                else:
                    auc = norm_cdf(max_d, self._mu, self._sigma) - norm_cdf(min_d, self._mu, self._sigma)
                
                
                # avoid zero-division error
                # This will never happen, unless an organism evolves a really extreme sigma
                if auc < 1e-100:
                    auc = 1e-100 
                    print("AUC was 0 with mu =", self._mu, "and sigma =", self._sigma)                
                                
                # avoid log(0) error when computing e_connector
                if numerator < 1e-100:
                    numerator = 1e-100
                
                # Normalize
                numerator = numerator / auc
                
                # Denominator
                # p(d) according to null model
            	# The probabity to observe d depends on the number of ways of getting d,
            	# given two randomly selected positions within a sequence of lentgth L.
            	# There are L - abs(d) ways of getting two random position at
            	# distance d within a sequence of length L.
            	# The total number of equally likely couples (-> distances) is L*(L-1).
            	# Therefore  p(d|null_model) = (L-abs(d)) / (L*(L-1))
                denominator = (s_dna_len - abs(d)) / (s_dna_len * (s_dna_len-1))
                
                # compute additive connector energy term
                e_connector = np.log2(numerator / denominator)
                
                # compute overall placement energy (connector + children)
                energy = possibility_1["energy"] + possibility_2["energy"] + e_connector

                # add placement to list of possible placements
                # The candidate placement includes:
                # the position (average of daughter nodes for connector)
                # the energy (connector + daughter nodes)
                # the blocked PSSM positions inherited from the daugther nodes
                # the list of all the scores of the recognizers
                possible_candidates.append({
                    "position": (possibility_1["position"] + possibility_2["position"]) / 2,
                    "energy": energy,
                    "lock_vector": possibility_1["lock_vector"] + possibility_2["lock_vector"],
                    "recognizers_scores": possibility_1["recognizers_scores"] + possibility_2["recognizers_scores"]
                    })

        # reverse sort the list of candidate placements based on energy
        possible_candidates.sort(key=lambda c: c["energy"], reverse=True)
        options = automatic_placement_options if is_automatic_placement_options else self.placement_options

        # return the truncated (at self.placement_options) list of best
        # possible placements
        return possible_candidates[:options]

    def set_node(self, node, _id) -> None:
        """Sets the node on a given ID

        Args:
            node: Node to insert in the tree
            _id: id to insert the node

        """

        if self.node1._id == _id:
            self.node1 = node
        elif self.node2._id == _id:
            self.node2 = node
        else:
            self.node1.set_node(node, _id)
            self.node2.set_node(node, _id)

    def reset_id(self, new_id: int) -> int:
        """Resets the id parameter based on a given_ID

        Args:
            new_id: Base id to set the subtree IDs

        Returns:
            id of the last id assigned in the tree

        """
        new_id = self.node1.reset_id(new_id)
        self._id = new_id
        new_id = self.node2.reset_id(new_id + 1)
        return new_id

    # pylint: disable=W0613
    def mutate(self, org_factory) -> None:
        """mutation for a connector

        Args:
            org_factory(organism_factory): Organism Facory
        """
        # print("Mutating Connector..." + str(self.ID))
        if random.random() < self.mutate_probability_sigma:
            # Alter sigma
            self._sigma += random.uniform(
                -self.mutate_variance_sigma, self.mutate_variance_sigma
            )

        if random.random() < self.mutate_probability_mu:
            # Alter mu
            self._mu += random.uniform(
                -self.mutate_variance_mu, self.mutate_variance_mu
            )

        if random.random() < self.mutate_probability_swap:
            # Swap connectors
            tmp_node = self.node1
            self.node1 = self.node2
            self.node2 = tmp_node

    # pylint: enable=W0613

    def print(self, distance: int) -> None:
        """It prints the connector mu, sigma values and its children values in
           tree structure

        Args:
            distance: Depth in the tree
        """
        print(
            "   |" * distance
            + " - C"
            + str(self.ID)
            + " m: {} s: {}".format(self.mu, self.sigma)
        )
        self.node1.print(distance + 1)
        self.node2.print(distance + 1)

    def export(self, export_file, level: int) -> None:
        """Exports Connector data to the given file

        Args:
            export_file (file): File to export the conector
            level: Depth in the tree

        """
        export_file.write(
            "\n"
            + "   |" * level
            + " - C"
            + str(self._id)
            + " m: {} s: {}".format(self._mu, self._sigma)
        )
        self.node1.export(export_file, level + 1)
        self.node2.export(export_file, level + 1)

    # pylint: disable=R0201
    def is_connector(self) -> bool:
        """node is connector

        Returns:
            True because is a connector
        """
        return True

    def is_pssm(self) -> bool:
        """node is not a pssm

        Returns:
            False because is a connector
        """
        return False

"""C object
Connects two nodes at a specific distance

"""

# pylint: disable=E0402
# type: ignore
import random
from .node_object import Node
import numpy as np


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

    # pylint: disable=R1702
    # pylint: disable=R0915
    def get_placement(
        self, sDNA: str, sDNAlen: int, blocks: list, blockers: list
    ) -> dict:
        """DEPRECATED

        This is a position/score propagation method, defined for connector
        objects.
        It is invoked by the placement method in the organism, for the root
        connector object, and calls itself recursively.
        The function calls itself until reaching terminal connectors, which
        call onto PSSM objects.
        At that point, the call to the getPlacement function in PSSM nodes
        leads to the evaluation of the PSSM node across all the sequence,
        and it returns the score/position pairs, sorted by descending score.

        The connector function then propagates this up, taking the middle
        position between both PSSMs and adding the connector energy to the
        energies provided by the PSSMs.

        The connector determines (i.e. freezes) the PSSM locations, adding
        them to the block list that is passed as a parameter.

        Further connector objects proceed in the same manner, computing
        middle distance and adding their energy contribution, until the
        energy reaches the root node, and is returned as the fitness for
        the organism.

        The energy contribution of each connector is: EN1 + EN2 + Tau * EC,
        where EN1 is the energy of its daugher element 1, and EN2 that of
        daughter element 2. The EC connector energy component is an
        exponential function controlled by the difference in the observed
        distance of the elements  of the connector with respect to an ideal
        mean distance (mu), and modulated by a dispersion parameter (sigma).
        Tau controls the "weight" of the connector contribution to energy.

        Args:
            sDNA: DNA sequence to compute energy
            sDNAlen: length of the DNA sequence
            blocks: TODO
            blockers: ids of the blocked nodes

        Returns:
            dictionary with the info about the energy
            "pspairs": TODO
            "blocked": TODO
            "blocker": TODO

        """

        # This tau shows how much value we give to the connector fit
        tau = self.tau

        # ask daughter nodes what their placement is
        # placement call will return a vector of positions (for PSSMs) and
        # their scores (sorted descending by score) , as well as an udpated
        # block/blocker vector
        node1 = self.node1.get_placement(sDNA, sDNAlen, blocks, blockers)
        node2 = self.node2.get_placement(sDNA, sDNAlen, blocks, blockers)

        # precompute connector energy term (not dependent on PSSM placement)
        logterm = np.log10(10 + self._sigma ** 2)

        maxenergy = -np.inf
        maxposition = 0
        max1 = 0
        max2 = 0
        placecnt = 0
        placeopt = self.placement_options
        confset = False
        # iterate over all possible configurations of sub-node placements
        # and determine the optimal one. this goes on for a _minimum_ number of
        # iterations (user-defined), but continues if a satisfactory
        # configuration has not been found (for instance, because
        # all tried configurations generated conflicts with blocked positions)
        # the two for loops iterate over possible subset of configurations
        # if a valid solution is not found, they iterate over next subset
        # this will take place when at least one of the nodes is a PSSM and has
        # returned more than one possible position (although in theory
        # connectors could also return more than one position)

        while not confset:
            # set range according to availability of positions
            if len(node1["pspairs"]) > placeopt:
                range1st = placecnt
                range1nd = placecnt + placeopt
            else:
                range1st = 0
                range1nd = 1
            if len(node2["pspairs"]) > placeopt:
                range2st = placecnt
                range2nd = placecnt + placeopt
            else:
                range2st = 0
                range2nd = 1

            for n1count in range(range1st, range1nd):
                for n2count in range(range2st, range2nd):
                    # compute connector energy terms that depend on PSSM
                    # placement
                    numerator = (
                        self._mu
                        - (
                            node2["pspairs"][n2count]["pos"]
                            - node1["pspairs"][n1count]["pos"]
                        )
                    ) ** 2
                    exponent = -1.0 * numerator / (1 + 2 * (self._sigma ** 2))
                    expterm = np.exp(exponent)
                    # compute additive connector energy term
                    e_connector = (tau / logterm) * expterm
                    # submodel energy: additive
                    try:
                        energy = (
                            node1["pspairs"][n1count]["energy"]
                            + node2["pspairs"][n2count]["energy"]
                        ) + e_connector
                        # submodel position: average of daughter positions
                        position = (
                            node1["pspairs"][n1count]["pos"]
                            + node2["pspairs"][n2count]["pos"]
                        ) / 2

                    except Exception as e:
                        print(e)
                        print("Values of\nn1count: {}\nn2count: {}\n".format(n1count, n2count))
                        print("{} {}".format(len(node2["pspairs"]), placeopt))
                        print("{} {}".format(len(node1["pspairs"]), placeopt))

                    # if BOTH nodes are pssms, avoid pssm overlapping
                    if self.node1.is_pssm() and self.node2.is_pssm():
                        # determine largest PSSM
                        p_len = (
                            self.node1.get_length()
                            if self.node1.get_length()
                            > self.node2.get_length()
                            else self.node2.get_length()
                        )
                        # determine if overlap exists, skip combo if there is
                        if (
                                abs(
                                    round(node2["pspairs"][n2count]["pos"])
                                    - round(node1["pspairs"][n1count]["pos"])
                                ) < p_len
                        ):
                            continue
                            #pass

                    # if ONE of the nodes is a PSSM, make sure its position
                    # does not overlap with blocked positions
                    if self.node1.is_pssm():
                        p_len = self.node1.get_length()
                        startpos1 = round(
                            node1["pspairs"][n1count]["pos"]
                        ) - round(p_len / 2)
                        blocked = False
                        for jpos in range(startpos1, startpos1 + p_len):
                            if jpos in blocks:
                                blocked = True
                                break
                        # if position has been blocked, skip this combo
                        if blocked:
                            continue
                    if self.node2.is_pssm():
                        p_len = self.node2.get_length()
                        startpos2 = round(
                            node2["pspairs"][n2count]["pos"]
                        ) - round(p_len / 2)
                        blocked = False
                        for jpos in range(startpos2, startpos2 + p_len):
                            if jpos in blocks:
                                blocked = True
                                break
                        # if position has been blocked, skip this combo
                        if blocked:
                            continue
                    # if this energy is the best so far, annotate it
                    # set confset to true, since we have obtained a valid
                    # configuration
                    if energy > maxenergy:
                        maxenergy = energy
                        maxposition = position
                        max1 = n1count
                        max2 = n2count
                        confset = True

                # increase base configuration counter
                placecnt = placecnt + 1

            # print('Max1: ',node1['pspairs'][max1]['pos'])
            # print('Max2: ',node2['pspairs'][max2]['pos'])

        # block the positions of this connector's PSSMs
        if self.node1.is_pssm():
            n1length = self.node1.get_length()
            blockstartpos1 = round(node1["pspairs"][max1]["pos"]) - round(
                n1length / 2
            )
            for blockade in range(n1length):
                blocks.append(blockstartpos1 + blockade)
                blockers.append(self.node1._id)
        if self.node2.is_pssm():
            n2length = self.node2.get_length()
            blockstartpos2 = round(node2["pspairs"][max2]["pos"]) - round(
                n2length / 2
            )
            for blockade in range(n2length):
                blocks.append(blockstartpos2 + blockade)
                blockers.append(self.node2._id)

        pair = {"pos": maxposition, "energy": maxenergy}

        return {"pspairs": [pair], "blocked": blocks, "blocker": blockers}

    # pylint: enable=R1702
    # pylint: enable=R0915
    def get_placement_2(
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

        possibilities_node_1 = self.node1.get_placement_2(
            s_dna,
            s_dna_len,
            automatic_placement_options,
            is_automatic_placement_options
                )
        possibilities_node_2 = self.node2.get_placement_2(
            s_dna,
            s_dna_len,
            automatic_placement_options,
            is_automatic_placement_options
                )

        # precompute the log term of the connector energy component, which is
        # not dependent on the placement of daugther nodes
        logterm = np.log10(10 + self._sigma ** 2)

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
                
                # compute the rest of the connector energy term
                numerator = (
                    self._mu
                    - (
                        possibility_2["position"] -
                        possibility_1["position"]
                    )
                ) ** 2
                exponent = -1.0 * numerator / (1 + 2 * (self._sigma ** 2))
                expterm = np.exp(exponent)
                # compute additive connector energy term
                e_connector = (self.tau / logterm) * expterm

                # compute overall placement energy (connector + children)
                energy = possibility_1["energy"] + possibility_2["energy"] + e_connector

                # add placement to list of possible placements
                # include the position (average of daughter nodes for connector)
                # the energy (connector + daughter nodes)
                # and the blocked PSSM positions inherited from the daugther nodes
                possible_candidates.append({
                    "position": (possibility_1["position"] + possibility_2["position"]) / 2,
                    "energy": energy,
                    "lock_vector": possibility_1["lock_vector"] + possibility_2["lock_vector"]
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
            self._sigma += random.randint(
                -self.mutate_variance_sigma, self.mutate_variance_sigma
            )

        if random.random() < self.mutate_probability_mu:
            # Alter mu
            self._mu += random.randint(
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

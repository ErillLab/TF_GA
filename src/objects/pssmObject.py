"""P object
Saves al specific PSSM data structure
"""
import random
from objects.nodeObject import Node
import numpy as np


class PssmObject(Node):
    """PSSM object
    """

    def __init__(self, pwm, config: dict) -> None:
        """PSSM constructor

        Args:
            pwm (numpy.array): PWM
            config: configuration from JSON file
        """

        super().__init__()
        self._id = 0
        self.length = len(pwm)  # length of the numpy array
        self.pwm = pwm  # numpy array of dictionaries
        self.pssm = None
        self.optimal_combination: list = []
        self.mutate_probability_random_col = config[
            "MUTATE_PROBABILITY_RANDOM_COL"
        ]
        self.mutate_probability_flip_cols = config[
            "MUTATE_PROBABILITY_FLIP_COL"
        ]
        self.mutate_probability_flip_rows = config[
            "MUTATE_PROBABILITY_FLIP_ROW"
        ]
        self.mutate_probability_shift_left = config[
            "MUTATE_PROBABILITY_SHIFT_LEFT"
        ]
        self.mutate_probability_shift_right = config[
            "MUTATE_PROBABILITY_SHIFT_RIGHT"
        ]
        self.pseudo_count = config["PSEUDO_COUNT"]
        self.upper_print_probability = config["UPPER_PRINT_PROBABILITY"]
        self.scan_reverse_complement = config["SCAN_REVERSE_COMPLEMENT"]
        # It first calculates PSSM Matrix based on  pwm
        self.recalculate_pssm()

    def count_nodes(self) -> int:
        """Counts the number of nodes. Leaves always return 1 as it is a single
        node

        Returns:
           The number of nodes(1)
        """
        return 1

    def get_node(self, objective: int, node_count: int) -> Node:
        """Get an specific node

        Args:
            objective: ID of the objective node to find
            node_count: number of nodes before the current node

        Returns:
            itself if its the node, otherwise return None

        """
        return self if objective == node_count else None

    def get_parent(self, _id: int) -> Node:
        """pssm objects cannot be parents

        Args:
            _id: Id of the child node

        Returns:
            Node with ID _id. PSSM recognizers cannot be parents.
        """
        return None

    def get_length(self) -> int:
        """Length of the pssm recognizer

        Returns:
            Length of the pwm
        """
        return self.length

    def mutate(self, org_factory) -> None:
        """Mutation operators associated to the PSSM recognizer

        Args:
            org_factory (OrganismFactory): Cretes organisms and Node components
        """

        if random.random() < self.mutate_probability_random_col:

            # Randomize PSSM column
            new_col = org_factory.get_pwm_column()
            # Select a random col in self.pwm
            column_to_update = random.randint(0, self.length - 1)
            # Insert it in that position
            self.pwm[column_to_update] = new_col

        if random.random() < self.mutate_probability_flip_cols:
            # Flip two columns

            col1, col2 = random.sample(range(self.length), 2)
            # Select two random columns and swap it
            tmp_col = self.pwm[col1]
            self.pwm[col1] = self.pwm[col2]
            self.pwm[col2] = tmp_col

        if random.random() < self.mutate_probability_flip_rows:
            # Flip two rows
            # Save values of two rows and swap it (its gonna be one by one)
            bases = ["a", "c", "g", "t"]
            random.shuffle(bases)
            base1, base2 = bases[:2]

            # Swap rows
            for i in range(self.length):
                tmp_base = self.pwm[i][base1]
                self.pwm[i][base1] = self.pwm[i][base2]
                self.pwm[i][base2] = tmp_base

        if random.random() < self.mutate_probability_shift_left:
            # Shift to right/left
            self.pwm = np.roll(self.pwm, -1)

        if random.random() < self.mutate_probability_shift_right:
            # Shift to right/left
            self.pwm = np.roll(self.pwm, 1)

        self.recalculate_pssm()

    # Calculate self.pssm based on self.pwm
    def recalculate_pssm(self) -> None:
        """ Calculates the PSSM based on the pwm values
        """
        tmp_pssm = []
        for column in self.pwm:
            # From pwm to pssm
            # log2(base/0.25) = log2(4.0*base)
            decimals = 2
            tmp_bases = []
            # cast to float so round function does not become crazy
            tmp_bases.append(
                float(np.log2(4.0 * column["c"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["t"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["g"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["a"] + self.pseudo_count))
            )

            tmp_pssm.append(
                {
                    "c": round(tmp_bases[0], decimals),
                    "t": round(tmp_bases[1], decimals),
                    "g": round(tmp_bases[2], decimals),
                    "a": round(tmp_bases[3], decimals),
                }
            )
        self.pssm = np.array(tmp_pssm)
        # Also calculate the optimal pssm combinations
        self.optimal_combination = [""]
        for position in tmp_pssm:
            max_bases = []
            max_base_score = float("-inf")
            for base in position:
                if position[base] > max_base_score:
                    max_bases = [base]
                    max_base_score = position[base]
                elif position[base] == max_base_score:
                    max_bases.append(base)

            tmp_optimal = []
            for base in max_bases:
                for comb in self.optimal_combination:
                    tmp_optimal.append(comb + base)

            self.optimal_combination = tmp_optimal
        # print(self.optimal_combination)

    def get_placement(
            self, s_dna: str, s_dna_len: int, blocks: list, blockers: list
    ) -> dict:
        """Searchs sequence with PSSM ojbect returns position and score

        Args:
            s_dna: DNA sequence
            s_dna_len: length of the DNA sequence
            blocks: TODO
            blockers: IDs of the blocked nodes

        Returns:
            dictionary with info about the energy
            "pspairs": TODO
            "blocked": TODO
            "blocker": TODO

        """

        scores_on_seq = []
        pssm_length = self.length

        # Check every position possible on the sequence
        for pos in range(s_dna_len - pssm_length):
            mypos = float(pos) + pssm_length / 2
            blocked = False
            for jpos in range(pos, pos + pssm_length):
                if jpos in blocks:
                    blocked = True
                    break
            # if the position has not been blocked by another PSSM
            if not blocked:
                # Compute the score of PSSM at position in sequence
                # append them to list
                scores_on_seq.append(
                    {
                        "pos": mypos,
                        "energy": self.get_score(
                            s_dna[pos: pos + pssm_length]
                        ),
                    }
                )
            # otherwise just add large negative
            else:
                scores_on_seq.append({"pos": mypos, "energy": -10000.0})

        # print (blocks)
        # print(scoresonseq,'\n')
        # sort list and return it
        scores_on_seq.sort(key=lambda k: k["energy"], reverse=True)

        return {
            "pspairs": scores_on_seq,
            "blocked": blocks,
            "blocker": blockers,
        }

    def get_all_pssm(self) -> list:
        """Adds himself as a pssm recognizer

        Returns:
            list of pssms
        """
        return [self]

    def get_score(self, s_dna: str) -> float:
        """a score to that DNA secuence

        Args:
            s_dna: dna partial sequence (length of the pssm)

        Returns:
            score assigned to s_dna. If reverse sequence is better, reverse
            score is returned

        """

        complement = {"a": "t", "t": "a", "g": "c", "c": "g"}
        # gets a score from pssm
        score = 0
        score_reverse = 0
        str_length = len(s_dna)
        for i in range(str_length):

            score += self.pssm[i][s_dna[i]]
            score_reverse += self.pssm[str_length - i - 1][
                complement[s_dna[str_length - i - 1]]
            ]
        # Returns the max binding score
        return (
            score
            if score > score_reverse or not self.scan_reverse_complement
            else score_reverse
        )

    def set_node(self, node: Node, _id: int) -> None:
        """Nodes cannot be setted from recognizer objects

        Args:
            node: Node to insert in the tree
            _id: id to insert the node
        """

    def reset_id(self, new_id: int) -> int:
        """Resets the id on the whole tree structure

        Args:
            new_id: id to assign to the current recognizer

        Returns:
            The next id to assign
        """

        self._id = new_id
        return new_id + 1

    def print(self, distance: int) -> None:
        """Print PSSM object (similar to Logo format)

        Args:
            distance: Depth in the tree
        """

        recognized = ""

        for position in self.pwm:
            base = "a"
            # Find max base
            for new_base in position.keys():
                if position[new_base] > position[base]:
                    base = new_base

            if position[base] >= self.upper_print_probability:
                base = base.upper()
            recognized += base

        print("   |" * distance + " - " + recognized)

    def export(self, export_file, level: int) -> None:
        """Exports pssm to a file

        Args:
            export_file: File to write the output
            level: Depth in the tree
        """
        recognized = ""

        for position in self.pwm:
            base = "a"
            # Find max base
            for new_base in position.keys():
                if position[new_base] > position[base]:
                    base = new_base
            # Change to uppercase based on probability
            if position[base] >= self.upper_print_probability:
                base = base.upper()
            recognized += base

        export_file.write("\n" + "   |" * level + " - " + recognized)
        # exportFile.write("\n")

    # pylint: disable=R0201
    def is_connector(self) -> bool:
        """node is not a connector

        Returns:
            False because is a pssm recognizer

        """
        return False

    def is_pssm(self):
        """node is a pssm recognizer

        Returns:
            True because is a pssm recognizer

        """
        return True

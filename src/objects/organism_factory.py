"""Organism Factory creates organisms, connectors and pssms
"""

import random
import json
import numpy
from .organism_object import OrganismObject
from .connector_object import ConnectorObject
from .pssm_object import PssmObject


class OrganismFactory:
    """Factory
    """

    def __init__(self, conf_org, conf_org_fac, conf_con, conf_pssm) -> None:

        self._id = 0
        self.initial_connector_probability = conf_org_fac[
            "INITIAL_CONNECTOR_PROBABILITY"
        ]
        self.reducer_probability_factor = conf_org_fac[
            "REDUCER_PROBABILITY_FACTOR"
        ]

        self.min_mu = conf_org_fac["MIN_MU"]
        self.max_mu = conf_org_fac["MAX_MU"]

        self.min_sigma = conf_org_fac["MIN_SIGMA"]
        self.max_sigma = conf_org_fac["MAX_SIGMA"]
        # The number of positions, PSSM object can recognize

        self.pwm_length = conf_org_fac["PWM_LENGTH"]

        self.pwm_probability_step = conf_org_fac[
            "PWM_PROBABILITY_STEP"
        ]  # It should be a BASE_PROBABILITY divisor Ex: 1, 2, 4, 5, 10, 25...
        self.pwm_probability_base = conf_org_fac["PWM_PROBABILITY_BASE"]
        self.pwm_probability_decimals = conf_org_fac[
            "PWM_PROBABILITY_DECIMALS"
        ]

        self.conf_org = conf_org
        self.conf_con = conf_con
        self.conf_pssm = conf_pssm

    def get_id(self) -> int:
        """Gives a new ID for an organism
        TODO: This should be a function so all the count of IDs, including
        assigned outside the class, keep consistency between all organisms

        Returns:
           a new non-repeated ID
        """
        self._id += 1
        return self._id

    def get_organism(self) -> OrganismObject:
        """It creates and returns a full organism datastructure

        Returns:
            A new organism based on JSON config file
        """

        new_organism = OrganismObject(self.get_id(), self.conf_org)
        root_node = None

        # Based on a random probability, we assign a connector or a PSSM object
        # to the root Node
        if random.random() < self.initial_connector_probability:
            root_node = self.create_connection(
                self.initial_connector_probability
                * self.reducer_probability_factor
            )
        else:
            root_node = self.create_pssm(self.pwm_length)

        new_organism.set_root_node(root_node)
        new_organism.reset_ids()

        return new_organism

    def create_connection(
            self, connection_probability: float
    ) -> ConnectorObject:
        """It returns a connector object with its nodes assigned depending on
        connectionProbability

        Args:
            connection_probability: Probability to generate a connector instead
                                    of a recognizer

        Returns:
            A new Connection with recognizers included
        """

        # Assign a random value to mu and sigma
        _mu = random.randint(self.min_mu, self.max_mu)
        _sigma = random.randint(self.min_sigma, self.max_sigma)

        # Create the new connection
        new_connection = ConnectorObject(_mu, _sigma, self.conf_con)

        # Set the connection node to a connector or PSSM object, based on a
        # random probability. If connector object is selected, probability of
        # getting another connector is reduced by PROBABILITY_REDUCED_FACTOR
        node1 = None
        if random.random() < connection_probability:
            node1 = self.create_connection(
                connection_probability * self.reducer_probability_factor
            )
        else:
            node1 = self.create_pssm(self.pwm_length)

        new_connection.set_node1(node1)

        node2 = None
        if random.random() < connection_probability:
            node2 = self.create_connection(
                connection_probability * self.reducer_probability_factor
            )
        else:
            node2 = self.create_pssm(self.pwm_length)

        new_connection.set_node2(node2)

        return new_connection

    def create_pssm(self, length: int) -> PssmObject:
        """It return a PSSM object with a specific length

        Args:
            min_length: minimum positions a recognizer can recognize
            max_length: maximum positions a recognizer can recognize

        Returns:
            A pssm with an initializated PWM
        """
        pwm = []
        # Generate as much as needed
        for _ in range(length):
            pwm.append(self.get_pwm_column())

        return PssmObject(numpy.array(pwm), self.conf_pssm)

    def get_pwm_column(self) -> dict:
        """Generates a single column of the pwm

        Returns:
            a random probability for each base [a, c, g, t]
        """

        initial_probability = (
            self.pwm_probability_base / self.pwm_probability_step
        )
        probabilities: list = []
        # Left probability is
        left_probability = initial_probability
        # Minimum and maximum number of probabilities to be generated
        min_probability = 0
        max_probability = 4
        # Number of decimals on the probability
        decimals = 2
        # Generate 4 random probabilities out of initial_probability, one for
        # each base

        # Add a probability while we have less than 3 and and total probability
        # is not 1
        while (
                left_probability > min_probability
                and len(probabilities) < max_probability - 1
        ):
            new_probability = random.randint(0, left_probability)
            probabilities.append(float(new_probability))
            left_probability -= new_probability
        # Add the last probability or fill with 0 probability
        if left_probability > 0:
            probabilities.append(initial_probability - sum(probabilities))
        else:
            while len(probabilities) < max_probability:
                probabilities.append(0.0)

        # Shuffle the array is needed so high probability is not always on
        # first positions
        random.shuffle(probabilities)

        # Transform probabilities array from integer
        # [0-(BASE_PROBABILITY / STEP)] to complementary float
        # probabilities [0.0-1.0]
        numpy_probabilities = (
            numpy.array(probabilities)
            * self.pwm_probability_step
            * (1 / self.pwm_probability_base)
        )
        probabilities = numpy_probabilities.tolist()

        # Return object with "decimals" decimals probability to each base
        return {
            "a": round(probabilities[0], decimals),
            "g": round(probabilities[1], decimals),
            "c": round(probabilities[2], decimals),
            "t": round(probabilities[3], decimals),
        }

    def import_organisms(self, file_name: str) -> list:
        """Import Organisms from file

        Args:
            file_name: Name of the file with the organisms to read as an input

        Returns:
            a list of organisms objects read from the file
        """
        organism_list = []
        organism_json = {}
        with open(file_name) as json_file:
            organism_json = json.load(json_file)

        for organism in organism_json:

            new_organism = OrganismObject(self.get_id(), self.conf_org)
            root_node = None

            if organism["rootNode"]["objectType"] == "pssm":
                root_node = self.import_pssm(organism["rootNode"])
            else:
                root_node = self.import_connector(organism["rootNode"])

            new_organism.set_root_node(root_node)
            new_organism.reset_ids()

            if "isTracked" in organism.keys():
                new_organism.set_is_tracked(organism["isTracked"])

            organism_list.append(new_organism)

        return organism_list

    def import_connector(self, connector: dict) -> ConnectorObject:
        """Import Connector from JSON object

        Args:
            connector: connector in dictionary format

        Returns:
            Connector object from given connector dictionary
        """
        new_connector = ConnectorObject(
            connector["mu"], connector["sigma"], self.conf_con
        )

        node1 = None
        if connector["node1"]["objectType"] == "pssm":
            node1 = self.import_pssm(connector["node1"])
        else:
            node1 = self.import_connector(connector["node1"])

        new_connector.set_node1(node1)

        node2 = None
        if connector["node2"]["objectType"] == "pssm":
            node2 = self.import_pssm(connector["node2"])
        else:
            node2 = self.import_connector(connector["node2"])

        new_connector.set_node2(node2)

        return new_connector

    def import_pssm(self, pssm: dict) -> PssmObject:
        """Import PSSM from JSON object

        Args:
            pssm: pssm recognizer in dictionary format

        Returns:
            PSSM Object from given  pssm dictionary

        """
        return PssmObject(numpy.array(pssm["pwm"]), self.conf_pssm)

    def export_organisms(self, a_organisms: list, filename: str) -> None:
        """Export a list of organisms to JSON format

        Args:
            a_organisms: list of organisms to export
            filename: name of the file to export all the organisms
        """
        list_json_organisms = []
        for o_organism in a_organisms:
            organism = {}
            if o_organism.root_node.is_connector():
                organism["rootNode"] = self.export_connector(
                    o_organism.root_node
                )
            else:
                organism["rootNode"] = self.export_pssm(o_organism.root_node)
            list_json_organisms.append(organism)

        with open(filename, "w+") as json_file:
            json.dump(list_json_organisms, json_file, indent=2)

    def export_connector(self, o_connector: ConnectorObject) -> dict:
        """Export connector object

        Args:
            o_connector: Connector to export

        Returns:
            Connector in dictionary format
        """
        connector = {}
        connector["objectType"] = "connector"
        connector["mu"] = o_connector._mu
        connector["sigma"] = o_connector._sigma

        # Check if its pssm
        if o_connector.node1.is_connector():
            connector["node1"] = self.export_connector(o_connector.node1)
        else:
            connector["node1"] = self.export_pssm(o_connector.node1)

        if o_connector.node2.is_connector():
            connector["node2"] = self.export_connector(o_connector.node2)
        else:
            connector["node2"] = self.export_pssm(o_connector.node2)

        return connector

    def export_pssm(self, o_pssm: PssmObject) -> dict:
        """Export PSSM object

        Args:
            o_pssm: PSSM object to export

        Returns:
            pssm in dictionary format

        """
        pssm = {}
        pssm["objectType"] = "pssm"
        pssm["pwm"] = o_pssm.pwm.tolist()
        return pssm

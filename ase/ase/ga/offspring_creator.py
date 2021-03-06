"""Base module for all operators that create offspring."""
import numpy as np
from random import random

from ase import Atoms
from ase.ga.atoms_attach import enable_raw_score_methods


class OffspringCreator(object):
    """Base class for all procreation operators

    Parameters:

    verbose: Be verbose and print some stuff

    """
    def __init__(self, verbose=False):
        self.descriptor = 'OffspringCreator'
        self.verbose = verbose

    def get_new_individual(self, parents):
        """Function that returns a new individual.
        Overwrite in subclass."""
        raise NotImplementedError

    def finalize_individual(self, indi):
        """Call this function just before returning the new individual"""
        indi.info['key_value_pairs']['origin'] = self.descriptor

        enable_raw_score_methods(indi)

        return indi

    @classmethod
    def initialize_individual(cls, parent, indi=None):
        """Initializes a new individual that inherits some parameters
        from the parent, and initializes the info dictionary.
        If the new individual already has more structure it can be
        supplied in the parameter indi."""
        if indi is None:
            indi = Atoms(pbc=parent.get_pbc(), cell=parent.get_cell())
        indi.info['keywords'] = []
        ## key_value_pairs for numbers and strings
        indi.info['key_value_pairs'] = {}
        ## data for lists and the like
        indi.info['data'] = {}

        return indi


class OperationSelector(object):
    """Class used to randomly select a procreation operation
    from a list of operations.

    Parameters:

    probabilities: A list of probabilities with which the different
        mutations should be selected. The norm of this list
        does not need to be 1.

    oplist: The list of operations to select from.
    """
    def __init__(self, probabilities, oplist):
        assert len(probabilities) == len(oplist)
        self.oplist = oplist
        self.rho = np.cumsum(probabilities)

    def __get_index__(self):
        v = random() * self.rho[-1]
        for i in range(len(self.rho)):
            if self.rho[i] > v:
                return i

    def get_new_individual(self, candidate_list):
        """Choose operator and use it on the candidate. """
        to_use = self.__get_index__()
        return self.oplist[to_use].get_new_individual(candidate_list)

    def get_operator(self):
        """Choose operator and return it."""
        to_use = self.__get_index__()
        return self.oplist[to_use]

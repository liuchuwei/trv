from collections import defaultdict
from typing import *

from models.gtf_model import GtfModel
from models.segment_model import SegmentModel


def dictionary_intersect(d1: DefaultDict[Any, List], d2: DefaultDict[Any, List]) -> DefaultDict[Any, List]:
    """Set intersection (&) operation on default dicitonary

    Arguments
    ---------
    d1: defaultdict
        First default dict
    d2: defaultdict
        Second default dict

    Returns
    -------
    A dictionary with the key the set intersection of the keys.
    If same key is present the entry will be combined using __add__
    """
    keys_set = set(d1) & set(d2)
    return defaultdict(list, ((k, d1[k] + d2[k]) for k in keys_set))


class MoleModel:
    def __init__(self) -> None:
        self.mappings_record: DefaultDict[GtfModel, List[SegmentModel]] = None
        # self.final_report: Tuple[int, int, int, int, int, int] = None

    def add_mappings_record(self, mappings_record: DefaultDict[GtfModel, List[SegmentModel]]) -> None:
        if self.mappings_record is None:
            self.mappings_record = mappings_record
        else:
            self.mappings_record = dictionary_intersect(self.mappings_record, mappings_record)

from typing import *

from models.read_model import Read


class ExonFeature:

    def __init__(self, start: int, end: int, kind: int, exin_no: str, transcript_model: Any=None) -> None:
        self.start = start
        self.end = end
        self.transcript_model = transcript_model
        self.kind = kind  # it should be ord("e"), ord("i"), ord("m"), ....
        self.exin_no = int(exin_no)
        self.is_validated = False

    def __lt__(self, other: Any) -> bool:
        if self.start == other.start:
            return self.end < other.end
        return self.start < other.start

    def ends_upstream_of(self, read: Read) -> bool:
        """The following situation happens
                                                            Read
                                               *|||segment|||-?-||segment|||????????
                ???????|||||Ivl|||||||||*

        """
        return self.end < read.pos  # NOTE: pos is diffetent from start, consider chagning

    def doesnt_start_after(self, segment: Tuple[int, int]) -> bool:
        """One of the following situation happens

                            *||||||segment|||||????????
            *||||Ivl|||||*
                *|||||||||||||Ivl||||||||||????????????
                                    *|||||||||||||Ivl||||||||||????????????
                                              *|||||||||||||Ivl||||||||||????????????

        """
        return self.start < segment[-1]

    def contains(self, segment: Tuple[int, int], minimum_flanking: int = 5) -> bool:
        """One of following situation happens

            *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

                  *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

                      *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking
        """
        return (segment[0] + minimum_flanking >= self.start) and (segment[-1] - minimum_flanking <= self.end) and (
                    (segment[-1] - segment[0]) > minimum_flanking)

    def start_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int = 5) -> bool:
        """The following situation happens

          *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking

        """
        return (segment[0] + minimum_flanking < self.start) and (segment[-1] - minimum_flanking > self.start)

    def end_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int = 5) -> bool:
        """The following situation happens

                                      *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking

        """
        return (segment[0] + minimum_flanking < self.end) and (segment[-1] - minimum_flanking > self.end)

    def intersects(self, segment: Tuple[int, int], minimum_flanking: int=5) -> bool:
        return (segment[-1] - minimum_flanking > self.start) and\
               (segment[0] + minimum_flanking < self.end)  # and ((segment[-1] - segment[0]) > minimum_flanking)

    def get_upstream_exon(self) -> Any:
        """To use only for introns. Returns the vcy.Feature corresponding to the neighbour exon downstream

        Note
        ----
        In a 15 exons transcript model:
        Upstream to intron10 is exon9 or the interval with inxex `18` if strand "+".
        Upstream to intron10 is exon11 or the interval with inxex `8` if strand "-"
        """
        if self.transcript_model.chromstrand[-1] == "+":
            ix = (self.exin_no * 2) - 2
        else:
            # in the case on strand -
            ix = len(self.transcript_model.list_features) - 2 * self.exin_no - 1
        return self.transcript_model.list_features[ix]

    # if self.chromstrand[-1] == "+":
    #             intron_number = self.list_features[-1].exin_no
    #         else:
    #             intron_number = self.list_features[-1].exin_no - 1

    def get_downstream_exon(self) -> Any:
        """To use only for introns. Returns the vcy.Feature corresponding to the neighbour exon downstream

        Note
        ----
        In a 15 exons transcript model:
        Downstream to intron10 is exon11 or the interval with index `20` if strand "+".
        Downtream to intron10 is exon10 or the interval with index `10` if strand "-"
        """
        if self.transcript_model.chromstrand[-1] == "+":
            ix = self.exin_no * 2
        else:
            # in the case on strand -
            ix = len(self.transcript_model.list_features) - 2 * self.exin_no + 1
        return self.transcript_model.list_features[ix]

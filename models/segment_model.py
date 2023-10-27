from typing import Tuple

from models.exon_model import ExonFeature


class SegmentModel:

    def __init__(self, segment: Tuple[int, int], feature: ExonFeature, is_spliced: bool = False) -> None:
        self.segment = segment
        self.feature = feature
        self.is_spliced = is_spliced  # this is really BAM_CREF_SKIP

    @property
    def skip_makes_sense(self) -> bool:
        """If the SKIP in the segment matches some extremity of the feature and therefore can be interpreted as a splice event
        """
        if not self.is_spliced:
            return True  # NOTE: maybe here I should raise an error because the property is not supposed to be called
        else:
            if abs(self.feature.start - self.segment[0]) <= 6 or abs(self.feature.end - self.segment[1]) <= 6:
                return True
            else:
                return False

    @property
    def maps_to_exon(self) -> bool:
        return self.feature.kind == 101  # ord("e")

    @property
    def maps_to_intron(self) -> bool:
        return self.feature.kind == 105  # ord("i")

from typing import *

class Read:
    """ Container for reads from sam alignment file """

    def __init__(self, umi: str, chrom: str, strand: str, pos: int, segments: List, clip5: Any, clip3: Any, ref_skipped: bool) -> None:
         self.umi, self.chrom, self.strand, self.pos, self.segments, self.clip5, self.clip3, self.ref_skipped = \
            umi, chrom, strand, pos, segments, clip5, clip3, ref_skipped

    @property
    def is_spliced(self) -> bool:
        return self.ref_skipped  # len(self.segments) > 1

    @property
    def start(self) -> int:
        return self.segments[0][0]

    @property
    def end(self) -> int:
        return self.segments[-1][1]

    @property
    def span(self) -> int:
        return self.end - self.start + 1

    def __lt__(self, other: Any) -> bool:
        if self.chrom == other.chrom:
            if self.start == other.start:
                return self.end < other.end
            return self.start < other.start
        return self.chrom < other.chrom

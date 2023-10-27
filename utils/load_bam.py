import logging
import random
import string
from typing import List, Tuple

import pysam

from models.read_model import Read

def parse_cigar_tuple(cigartuples: List[Tuple], pos: int) -> Tuple[List[Tuple[int, int]], bool, int, int]:
    segments = []
    hole_to_remove = set()
    ref_skip = False
    clip5 = clip3 = 0
    p = pos
    for i, (operation_id, length) in enumerate(cigartuples):
        if operation_id == 0:  # vcy.CIGAR[operation_id] == "BAM_CMATCH"
            segments.append((p, p + length - 1))
            p += length
        elif operation_id == 3:  # A splice || vcy.CIGAR[operation_id] == 'BAM_CREF_SKIP'
            ref_skip = True
            p += length
        elif operation_id == 2:  # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
            if length <= 3:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
            p += length
        elif operation_id == 4:  # bases at 5' or 3' are NOT part of the alignment || vcy.CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
            if p == pos:
                clip5 = length  # At start of alignment
            else:
                clip3 = length  # Must be at end of alignment vcy.CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
            p += length
        elif operation_id == 1:  # An insertion BAM_CINS
            if length <= 3:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
            # else do nothing
            # NOTE: maybe we should make so that the reads get discarded
        elif operation_id == 5:  # BAM_CHARD_CLIP
            logging.warn("Hard clip was encountered! All mapping are assumed soft clipped")

    # Merge segments separated by small insertions and deletions
    for a, b in enumerate(sorted(hole_to_remove)):  # NOTE maybe sorted is not required realy
        segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])

    return segments, ref_skip, clip5, clip3


def iter_alignment(bamfile):
    fin = pysam.AlignmentFile(bamfile)  # type: pysam.AlignmentFile
    for i, read in enumerate(fin):
        if i % 10000000 == 0:
            logging.debug(f"Read first {i // 1000000} million reads")  # count
        if read.is_unmapped:
            continue
        if read.get_tag("NH") != 1:
            continue
        umi = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in
                      range(12))  # read.get_tag(self.umibarcode_str)

        strand = '-' if read.is_reverse else '+'
        chrom = fin.get_reference_name(read.rname)

        if chrom.startswith('chr'):
            if "_" in chrom:
                chrom = chrom.split("_")[1]
            else:
                chrom = chrom[3:]
                if chrom == "M":
                    chrom = "MT"

        pos = read.reference_start + 1  # reads in pysam are always 0-based, but 1-based is more convenient to wor with in bioinformatics
        segments, ref_skipped, clip5, clip3 = parse_cigar_tuple(read.cigartuples, pos)

        read_object = Read(umi, chrom, strand, pos, segments, clip5, clip3, ref_skipped)

        if read_object.span > 3000000:  # Longest locus existing
            logging.warn(f"Trashing read, too long span\n{read.tostring(fin)}")
        else:
            yield read_object

    fin.close()
    yield None
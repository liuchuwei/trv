from collections import defaultdict
from typing import Tuple, List

from models.exon_model import ExonFeature
def LoadMaskGTF(mask_gtf, tolerance = 5):
    gtf_lines = [line for line in open(mask_gtf) if not line.startswith('#')]

    def sorting_key(entry: str) -> Tuple[str, bool, int, str]:
        """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
        x = entry.split("\t")
        return (
        x[0], x[6] == "+", int(x[3]), entry)  # The last element of the touple corresponds to the `last resort comparison`


    gtf_lines = sorted(gtf_lines, key=sorting_key)



    line = gtf_lines.pop(0)
    fields = line.rstrip().split('\t')
    chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
    if chrom[:3].lower() == "chr":
        chrom = chrom[3:]
    start = int(start_str)
    end = int(end_str)
    chromstrand: str = chrom + strand

    # Set this tu the current entry
    curr_chrom: str = chrom
    curr_feature_class: str = feature_class
    curr_feature_type: str = feature_type
    curr_start: int = start
    curr_end: int = end
    curr_n: int = 1
    curr_strand: str = strand
    curr_tags: str = tags
    curr_chromstrand: str = chromstrand

    mask_chromstrand = defaultdict(list)  # type: Dict[str, List]
    repeat_ivls_list: List[ExonFeature] = []

    for line in gtf_lines:

        fields = line.rstrip().split('\t')
        chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields

        if chrom[:3].lower() == "chr":
            chrom = chrom[3:]
        start = int(start_str)
        end = int(end_str)
        chromstrand = chrom + strand

        if chromstrand != curr_chromstrand:
            mask_chromstrand[curr_chromstrand] = repeat_ivls_list
            repeat_ivls_list = []
            curr_chrom = chrom
            curr_strand = strand
            curr_chromstrand = curr_chrom + curr_strand

        if start > curr_end + tolerance:
            repeat_ivls_list.append(ExonFeature(start=curr_start, end=curr_end, kind=ord("r"), exin_no=curr_n))

            curr_start = start
            curr_end = end
            curr_n = 1  # number of original repeat regions included
            curr_tags = tags

        else:
            curr_end = end
            curr_n += 1
            gap = start - curr_end
            curr_tags = f"{curr_tags} gap {gap}; {tags}" if gap > 0 else curr_tags + tags

    n = 0
    for chromstrand, feature_list in mask_chromstrand.items():
        feature_list.sort()  # relies on the __lt__ method of vcy.Feature
        n += len(feature_list)

    return mask_chromstrand
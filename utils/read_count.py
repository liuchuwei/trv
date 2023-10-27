from typing import *
import numpy as np
from models.molecular_model import MoleModel
from collections import defaultdict

def reverse(strand: str) -> str:
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"
    else:
        raise ValueError(f"Unknown strand {strand}")


def count(molitem, dict_layers_columns, geneid2ix):

    spliced = dict_layers_columns["spliced"]
    unspliced = dict_layers_columns["unspliced"]
    ambiguous = dict_layers_columns["ambiguous"]
    spanning = dict_layers_columns["spanning"]

    if len(set(i.geneid for i in molitem.mappings_record.keys())) == 1:
        gene_check: Set[str] = set()

        has_onlyintron_model = 0
        has_only_span_exin_model = 1
        has_onlyexo_model = 0
        has_mixed_model = 0
        multi_gene = 0
        for transcript_model, segments_list in molitem.mappings_record.items():
            gene_check.add(transcript_model.geneid)
            if len(gene_check) > 1:
                multi_gene = 1
            has_introns = 0
            has_exons = 0
            has_exseg_with_spliced_flag = 0
            has_exin_intron_span = 0
            for segment_match in segments_list:
                if segment_match.maps_to_intron:
                    has_introns = 1
                    if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                        downstream_exon = segment_match.feature.get_downstream_exon()
                        if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                            has_exin_intron_span = 1
                    if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                        upstream_exon = segment_match.feature.get_upstream_exon()
                        if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                            has_exin_intron_span = 1
                elif segment_match.maps_to_exon:
                    has_exons = 1
                    if segment_match.is_spliced:
                        has_exseg_with_spliced_flag = 1
            if has_introns and not has_exons:
                has_onlyintron_model = 1
            if has_exons and not has_introns:
                has_onlyexo_model = 1
            if has_exons and has_introns and not has_exin_intron_span:
                has_valid_mixed_model = 1
                has_mixed_model = 1
            if not has_exin_intron_span:
                has_only_span_exin_model = 0

        if multi_gene:
            # Many genes are compatible with the observation, do not count
            return
        else:
            if not len(molitem.mappings_record):
                # NOTE it does not happen for Smartseq2
                # No gene is compatible with the observation, do not count
                return
            else:
                if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                    # More common situation, normal exonic read, count as spliced
                    gene_ix = geneid2ix[transcript_model.geneid]
                    spliced[gene_ix, 0] += 1
                    return
                if has_only_span_exin_model:
                    # NOTE This is what I want to count as spanning
                    # All the compatible transcript models have spanning exon-intron boundaries, count unspliced
                    gene_ix = geneid2ix[transcript_model.geneid]
                    spanning[gene_ix, 0] += 1
                    return
                if has_onlyintron_model and not has_mixed_model and not has_onlyexo_model:
                    gene_ix = geneid2ix[transcript_model.geneid]
                    unspliced[gene_ix, 0] += 1
                    return
                if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                    # Ambiguity among the transcript models compatible with the mapping, most common case! Count ambiguous
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, 0] += 1
                    return
                if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                    # NOTE has_mixed model is used only here in this logic
                    # Ambiguity among the transcript models compatible with the mapping. Count ambiguous
                    gene_ix = geneid2ix[transcript_model.geneid]
                    ambiguous[gene_ix, 0] += 1
                    return

def ReadCount(Read, feature_indexes, mask_indexes, geneid2ix):
    reads_to_count = []

    for r in Read:

        if r is None:
            molitems: DefaultDict[str, MoleModel] = defaultdict(MoleModel)
            reads_to_count.sort()

            for read in reads_to_count:
                # Consider the correct strand
                ii = feature_indexes[f"{read.chrom}{read.strand}"]
                iir = feature_indexes[f"{read.chrom}{reverse(read.strand)}"]
                iim = mask_indexes[f"{read.chrom}{read.strand}"]
                iimr = mask_indexes[f"{read.chrom}{reverse(read.strand)}"]

                # Check if read is fully inside a masked region, in that case skip it
                if iim.has_ivls_enclosing(read) or iimr.has_ivls_enclosing(read):
                    continue

                # Look for overlap between the intervals and the read
                mappings_record = ii.find_overlapping_ivls(read)

                if len(mappings_record):
                    umi = f"{read.umi}"
                    molitems[umi].add_mappings_record(mappings_record)

                mappings_record_r = iir.find_overlapping_ivls(read)
                if len(mappings_record_r):
                    umi = f"{read.umi}"
                    molitems[umi].add_mappings_record(mappings_record_r)

            dict_layers_columns: Dict[str, np.ndarray] = {}
            shape = (len(geneid2ix), 1)

            for layer_name in ['spliced', 'unspliced', 'ambiguous', 'spanning']:
                dict_layers_columns[layer_name] = np.zeros(shape, dtype='uint16', order="C")


            for umi, molitem in molitems.items():

                count(molitem, dict_layers_columns, geneid2ix)

        if r is not None:
            reads_to_count.append(r)

    return dict_layers_columns
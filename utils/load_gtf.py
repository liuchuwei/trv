import re
from collections import OrderedDict
from typing import Tuple, Dict, List

from models.gene_model import GeneInfo
from models.gtf_model import GtfModel
from models.exon_model import ExonFeature

def LoadGTF(gtf_file):

    gtf_lines = [line for line in open(gtf_file) if not line.startswith('#')]

    def sorting_key(entry: str) -> Tuple[str, bool, int, str]:
        """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
        x = entry.split("\t")
        return (x[0], x[6] == "+", int(x[3]),
                entry)  # The last element of the touple corresponds to the `last resort comparison`

    gtf_lines = sorted(gtf_lines, key=sorting_key)

    chrm_strand : Dict[str, Dict[str, GtfModel]] = {}
    curr_chrom = None
    gtf_models : List[GtfModel] = []

    for nth_line, line in enumerate(gtf_lines):
        fields = line.rstrip().split('\t')
        chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields

        # removing the chr from the name and removing part after the . if present
        if "chr" in chrom[:4]:
            chrom = chrom[3:]  # NOTE before it was chrom[3:].split(".")[0]
        else:
            pass

        # extract exon information
        ## Define the regexes for the parsing
        regex_trid = re.compile('transcript_id "([^"]+)"')
        regex_trname = re.compile('transcript_name "([^"]+)"')
        regex_geneid = re.compile('gene_id "([^"]+)"')
        regex_genename = re.compile('gene_name "([^"]+)"')
        regex_exonno = re.compile('exon_number "*?([\w]+)')  # re.compile('exon_number "([^"]+)"')

        if chrom + strand != curr_chrom:

            if curr_chrom is not None:  # Every time with exception with first and the last chromosome
                chrm_strand[curr_chrom] = features
                gtf_models.append(features)

            features = OrderedDict()
            curr_chrom = chrom + strand

        if feature_type in ("exon"):
            trid = regex_trid.search(tags).group(1)
            _trname_search = regex_trname.search(tags)

            if _trname_search is None:
                trname = trid
            else:
                trname = _trname_search.group(1)

            geneid = regex_geneid.search(tags).group(1)
            _genename_search = regex_genename.search(tags)

            if _genename_search is None:
                genename = geneid
            else:
                genename = _genename_search.group(1)

            exonno = regex_exonno.search(tags).group(1)

            start = int(start_str)
            end = int(end_str)

            chromstrand = chrom + strand

            try:
                features[trid].append_exon(ExonFeature(start=start, end=end, kind=ord("e"), exin_no=exonno))
            except KeyError:
                features[trid] = GtfModel(trid=trid, trname=trname, geneid=geneid, genename=genename,
                                                     chromstrand=chromstrand)
                features[trid].append_exon(ExonFeature(start=start, end=end, kind=ord("e"), exin_no=exonno))


    chrm_strand[curr_chrom] = features
    gtf_models.append(features)
    geneid2ix: Dict[str, int] = {}
    genes: Dict[str, GeneInfo] = {}

    for gtf_model in gtf_models:

        for name, trmodel in gtf_model.items():

            if trmodel.geneid in geneid2ix:
                if genes[trmodel.geneid].start > trmodel.start:
                    genes[trmodel.geneid].start = trmodel.start
                if genes[trmodel.geneid].end < trmodel.end:
                    genes[trmodel.geneid].end = trmodel.end
            else:
                geneid2ix[trmodel.geneid] = len(geneid2ix)
                genes[trmodel.geneid] = GeneInfo(trmodel.genename, trmodel.geneid, trmodel.chromstrand,
                                                          trmodel.start, trmodel.end)

    return chrm_strand, geneid2ix, genes
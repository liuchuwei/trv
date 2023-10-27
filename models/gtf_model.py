from typing import List
from models.exon_model import ExonFeature

class GtfModel:

    def __init__(self, trid: str, trname: str, geneid: str, genename: str, chromstrand: str) -> None:
        self.trid = trid
        self.trname = trname
        self.geneid = geneid
        self.genename = genename
        self.chromstrand = chromstrand
        self.list_features: List[ExonFeature] = []

    def __iter__(self) -> ExonFeature:
        for i in self.list_features:
            yield i
    @property
    def start(self) -> int:
        """ NOTE: This should be accessed only after the creation of the transcript model is finished
        (i.e.) after append_exon has been called to add all the exons/introns
        """
        return self.list_features[0].start

    @property
    def end(self) -> int:
        """NOTE: This should be accessed only after the creation of the transcript model is finished
        (i.e.) after append_exon has been called to add all the exons/introns
        """
        return self.list_features[-1].end

    def append_exon(self, exon_feature: ExonFeature) -> None:
        """Append an exon and create an intron when needed

        Arguments
        ---------
        exon_feature: ExonFeature
            A feature object represneting an exon to add to the gtf model.
        """
        exon_feature.transcript_model = self
        if len(self.list_features) == 0:
            # first/last exon
            self.list_features.append(exon_feature)
        else:
            # Some exon already existed
            if self.chromstrand[-1] == "+":
                intron_number = self.list_features[-1].exin_no
            else:
                intron_number = self.list_features[-1].exin_no - 1
            self.list_features.append(ExonFeature(start=self.list_features[-1].end + 1,
                                                  end=exon_feature.start - 1,
                                                  kind=ord("i"),
                                                  exin_no=intron_number,
                                                  transcript_model=self))
            self.list_features.append(exon_feature)
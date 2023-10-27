from collections import defaultdict
from typing import List, Dict

from models.exon_model import ExonFeature
from models.gtf_model import GtfModel
from models.read_model import Read
from models.segment_model import SegmentModel


class FeatureIndex:
    """ Search help class used to find the intervals that a read is spanning """
    def __init__(self, ivls: List[ExonFeature]=[]) -> None:
        self.ivls = ivls
        self.ivls.sort()  # NOTE: maybe I am sorting twice check what I do upon creation
        self.iidx = 0  # index of the current interval
        self.maxiidx = len(ivls) - 1
        # NOTE needs to be changed to consider situations of identical intervals but different objects

    @ property
    def last_interval_not_reached(self) -> bool:
        return self.iidx < self.maxiidx


    def has_ivls_enclosing(self, read: Read) -> bool:
        """Finds out if there are intervals that are fully containing all the read segments

        Args
        ----
        read: vcy.Read
            the read object to be analyzed

        Returns
        -------
        respones: bool
            if one has been found

        """
        if len(self.ivls) == 0:
            return False

        ivl = self.ivls[self.iidx]  # current interval

        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and ivl.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            ivl = self.ivls[self.iidx]

        for segment in read.segments:
            segment_matchtype = 0
            # Local search for each segment move a little forward (this just moves a coupple of intervals)
            i = self.iidx
            ivl = self.ivls[self.iidx]
            while i < self.maxiidx and ivl.doesnt_start_after(segment):
                matchtype = 0  # No match
                if ivl.contains(segment):
                    matchtype = 1
                if ivl.start_overlaps_with_part_of(
                        segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= 2
                if ivl.end_overlaps_with_part_of(segment):  # NOTE: should this be elif or it makes sense to allow both?
                    matchtype |= 2

                segment_matchtype |= matchtype
                # move to the next interval
                i += 1
                ivl = self.ivls[i]

            # If one of the segments does not match inside a repeat we return false
            if segment_matchtype ^ 1:
                return False
        # If I arrive at this point of the code all the segments matched inside
        return True

    def find_overlapping_ivls(self, read: Read) -> Dict[GtfModel, List[SegmentModel]]:
        """Finds the possible overlaps between Read and Features and return a 1 read derived mapping record

        Arguments
        ---------
        read: vcy.Read
            the read object to be analyzed

        Returns
        -------
        mapping_record: Dict[vcy.TranscriptModel, List[vcy.SegmentMatch]]
            A record of the mappings by transcript model.
            Every entry contains a list of segment matches that in turn contains information on the segment and the feature

        Note
        ----
        - It is possible that a segment overalps at the same time an exon and an intron (spanning segment)
        - It is not possible that a segment overalps at the same time two exons. In that case the read is splitted
        into two segments and  the Read attribute `is_spliced == True`.
        - Notice that the name of the function might be confousing. if there is a non valid overallapping an empty mappign record will be return
        - Also notice that returning an empty mapping record will cause the suppression of the counting of the molecule

        """

        mapping_record: Dict[GtfModel, List[SegmentModel]] = defaultdict(list)

        if len(self.ivls) == 0:
            return mapping_record

        feature: ExonFeature = self.ivls[self.iidx]  # current interval
        # Move forward until we find the position we will never search left of again (because the reads are ordered)
        while self.last_interval_not_reached and feature.ends_upstream_of(read):
            # move to the next interval
            self.iidx += 1
            feature = self.ivls[self.iidx]

        # Loop trough the mapping segments of a read (e.g. just one of an internal exon, generally 2 for a splice. for intron???)
        for seg_n, segment in enumerate(read.segments):
            # Local search for each segment move a little forward (this just moves a couple of intervals)
            i = self.iidx
            feature = self.ivls[i]
            while i < self.maxiidx and feature.doesnt_start_after(segment):
                # NOTE in this way the checks will be repeated for every clone of an interval in each transcript model
                # I might want to cache the check results of the previous feature
                # it was  if feature.contains(segment) or feature.start_overlaps_with_part_of(segment) or feature.end_overlaps_with_part_of(segment)
                # but I changed to the more simple
                if feature.intersects(segment) and (segment[-1] - segment[0]) > 5:
                    mapping_record[feature.transcript_model].append(SegmentModel(segment, feature, read.is_spliced))
                #  move to the next interval
                i += 1
                feature = self.ivls[i]

        # NOTE: Removing first the one with less match and then requiring a splicing matching the transcript model is very stringent
        # It could be that for short ~10bp SKIP sequences the alligner has made a mistake and this might kill the whole molecule
        # NOTE: the code below is not very efficient
        if len(mapping_record) != 0:
            # Remove transcript models that are suboptimal match
            # in alternative one could use len(read.segment)
            max_n_segments = len(max(mapping_record.values(), key=len))
            for tm, segmatch_list in list(mapping_record.items()):
                if len(segmatch_list) < max_n_segments:
                    del mapping_record[tm]

        # NOTE: the code below is not very efficient, would be nice to avoid for loops
        # NOTE: potentailly bad effects: it could kill a lot of molecules if transcript models are not annotated correctly or missing
        if len(mapping_record) != 0:
            # A SKIP mapping needs to be explainaible by some kind of exon-exon exon-intron junction!
            # So if it falls internally, the TM needs to be removed from the mapping record
            for tm, segmatch_list in list(mapping_record.items()):
                for sm in segmatch_list:
                    if not sm.skip_makes_sense:
                        del mapping_record[tm]
                        break

        return mapping_record

# -*- coding: utf-8 -*-

import argparse

import os
import sys

from tookit import Tookit
from collections import defaultdict
from itertools import chain
from typing import DefaultDict, Dict, Set


from models.feature_index import FeatureIndex
from utils.get_output import GetOutput
from utils.load_bam import iter_alignment
from utils.load_gtf import LoadGTF
from utils.load_mask_gtf import LoadMaskGTF
from utils.read_count import ReadCount

def BuildingIndex(tools, FLAGS):

    cmd = "%s \
    --runMode genomeGenerate \
    --genomeDir %s \
    --genomeFastaFiles %s \
    --sjdbOverhang 100 \
    --sjdbGTFfile %s \
    --runThreadN 8 \
    --limitGenomeGenerateRAM=100000000000" % (tools.STAR, FLAGS.reference, FLAGS.fa, FLAGS.gtf)

    os.system(cmd)

def Align1st(tools, FLAGS):

    cmd = "%s \
    --genomeDir %s \
    --readFilesIn %s \
    --runThreadN 16 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --readFilesCommand cat \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMtype None \
    --outSAMmode None \
    --outFileNamePrefix %s" % (tools.STAR, FLAGS.reference, " ".join(FLAGS.fastq), FLAGS.prefix)

    os.system(cmd)


def InterIndex(tools, FLAGS):

    fldir = FLAGS.out + "/" + FLAGS.prefix
    outtab = fldir + "/" + FLAGS.prefix + "SJ.out.tab"

    cmd = "%s \
        --runMode genomeGenerate \
        --genomeDir %s \
        --genomeFastaFiles %s \
        --sjdbOverhang 100 \
        --runThreadN 16 \
        --sjdbFileChrStartEnd %s" % (tools.STAR, fldir, FLAGS.fa, outtab)

    os.system(cmd)


def Align2nd(tools, FLAGS):

    fldir = FLAGS.out + "/" + FLAGS.prefix

    cmd = "%s \
    --genomeDir %s \
    --readFilesIn %s \
    --runThreadN 16 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 0 \
    --readFilesCommand cat \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMheaderHD @HD VN:1.4 \
    --outSAMattrRGline ID:571 \
    --outFileNamePrefix %s" % (tools.STAR, fldir, " ".join(FLAGS.fastq), FLAGS.prefix)

    os.system(cmd)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Mapping and splicing calculate.')
    parser.add_argument('-gtf', '--gtf', required=True, help="Path of gtf file")
    parser.add_argument('-mask_gtf', '--mask_gtf', required=True, help="Path of mask gtf file")
    parser.add_argument('-fa', '--fa', required=False, help="Path of fasta file")
    parser.add_argument('-fastq', '--fastq', required=True, dest="fastq", nargs='+', help="Path of fastq files")
    parser.add_argument('-reference', '--reference', required=True, help="Directory of reference")
    parser.add_argument('-prefix', '--prefix', required=True, help="prefix of out files")
    parser.add_argument('-out', '--out', required=True, help="Output directory")

    args = parser.parse_args(sys.argv[1:])

    global FLAGS
    FLAGS = args

    tools = Tookit()

    '1.BuildingIndex'
    if not os.path.exists(FLAGS.reference + "/" + "SAindex"):
        BuildingIndex(tools, FLAGS)

    '2.Alignment 1st'
    fldir = FLAGS.out + "/" + FLAGS.prefix
    if not os.path.exists(fldir):
        os.mkdir(fldir)

    Align1st(tools, FLAGS)

    Align1out = [item for item in os.listdir("./") if item.startswith(FLAGS.prefix)]

    for item in Align1out:
        cmd = "mv %s %s" % (item, fldir)
        os.system(cmd)

    '3.Intermediate Index Generation'
    InterIndex(tools, FLAGS)

    '4.Alignment 2nd'
    Align2nd(tools, FLAGS)

    Align2out = [item for item in os.listdir("./") if item.startswith(FLAGS.prefix)]
    fldir = FLAGS.out + "/" + FLAGS.prefix
    for item in Align2out:
        cmd = "mv %s %s" % (item, fldir)
        os.system(cmd)

    '5.Tidy'
    fldir = FLAGS.out + "/" + FLAGS.prefix
    fls = os.listdir(fldir)
    fls.remove(FLAGS.prefix + "Aligned.sortedByCoord.out.bam")

    for item in fls:
        fl = fldir + "/" + item
        cmd = "rm %s" % (fl)
        os.system(cmd)

    '6.Sort bam files'
    fldir = FLAGS.out + "/" + FLAGS.prefix
    bamfls = fldir + "/" + FLAGS.prefix + "Aligned.sortedByCoord.out.bam"
    sortfls = fldir + "/" + FLAGS.prefix + "sort.bam"
    cmd = '%s sort -l 7 -m 1M -t NOTAG -O BAM -@ 1 -o %s %s' % (tools.samtools, sortfls, bamfls)
    # os.system(cmd)

    '7.Splicing calculate'
    # read gtf files
    print('loading gtf files...')
    gtf_file = FLAGS.gtf
    chrm_strand, geneid2ix, genes = LoadGTF(gtf_file)

    feature_indexes: DefaultDict[str, FeatureIndex] = defaultdict(FeatureIndex)
    for chromstrand_key, annotions_ordered_dict in chrm_strand.items():
        feature_indexes[chromstrand_key] = FeatureIndex(
            sorted(chain.from_iterable(annotions_ordered_dict.values())))

    # read mask gtf files
    print('loading mask gtf files...')
    mask_gtf = FLAGS.mask_gtf
    mask_chromstrand = LoadMaskGTF(mask_gtf)

    mask_indexes: DefaultDict[str, FeatureIndex] = defaultdict(FeatureIndex)
    for chromstrand_key, annotions_list in mask_chromstrand.items():
        mask_indexes[chromstrand_key] = FeatureIndex(annotions_list)  # This suould be sorted

    # read bamfiles
    bamfile = sortfls
    Read = iter_alignment(bamfile)

    # counting
    print('Counting...')
    readcount = ReadCount(Read, feature_indexes, mask_indexes, geneid2ix)

    # output
    print('Finish!...')

    out = GetOutput(genes, readcount)

    fldir = FLAGS.out + "/" + FLAGS.prefix
    outfl = fldir + "/" + FLAGS.prefix + "splice.tsv"
    out.to_csv(outfl, index=False)
# trv

## Dependence
[![](https://img.shields.io/badge/STAR-2.7-green)](https://github.com/alexdobin/STAR) 
[![](https://img.shields.io/badge/samtools-v1.1.7-green)](https://github.com/samtools/samtools)  

## Genome
[![](https://img.shields.io/badge/hg38-fasta-orange)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz)
[![](https://img.shields.io/badge/genecodeV44-gtf-orange)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz)


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#Intallation">Intallation</a>
    </li>
    <li><a href="#Usage">Usage</a></li>
    <li><a href="#References">References</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#Contact">Contact</a></li>
  </ol>
</details>


## Intallation
1.Clone the project
   ```sh
   git clone https://github.com/liuchuwei/caNano.git
   ```
2.Install conda environment
   ```sh
   conda env create -f caNano.yaml
   ```
3.prepare tookit: 

check and modify the tookit.py file (in 'utils' directory).

## Usage

   ```sh
   python caNano.py -gtf <gencode.gtf> -mask_gtf <rmsk.gtf> -fastq <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,... -reference <reference directory> -prefix <prefix> -out <out directory>
   ```

## License
Distributed under the GPL-2.0 License License. See LICENSE for more information.

## Contact
liuchw3@mail2.sysu.edu.cn

## Reference

    
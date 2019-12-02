# Tools for the quality control of duplex sequencing data

This tools allow a deeper insight into duplex sequencing data. They can be used within the [Galaxy platform](http://usegalaxy.org) under the section Du Novo, on your local Galaxy installation from the [toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/duplex_family_size_distribution) and on the command line. 

## Dependencies
This tools were build with Python 2.7.

## Usage
A detailed description of all tools can be found on [Galaxy](http://usegalaxy.org) about its parameters, input and output files.

### TD: Tag distance analysis of duplex tags
Tags used in Duplex Sequencing (DS) are randomized barcodes, e.g 12 base pairs long. Since each DNA fragment is labeled by two tags at each end there are theoretically 4 to the power of (12+12) unique combinations. However, the input DNA in a typical DS experiment contains only ~1,000,000 molecules creating a large tag-to-input excess (4^24 ≫ 1,000,000). Because of such excess it is highly unlikely to tag distinct input DNA molecules with highly similar barcodes.

This tool calculates the number of nucleotide differences among tags (tag distance), also known as [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance). In this context the Hamming distance is simply the number of differences between two tags. The tool compares in a randomly selected subset of tags (default n=1000), the difference between each tag of the subset with the tags of the complete dataset. Each tag will differ by a certain number of nucleotides with the other tags; yet the tool uses the smallest difference observed with any other tag.

The tool can be used this tool via the command line with its default settings as the following:

`$ python2 td.py --inputFile tag_file.tabular --inputName1 tag_file.tabular --sample_size 1000 --subset_tag 0 --nproc 8 --rel_freq --minFS 1 --maxFS 0 --nr_above_bars --output_pdf out_file.pdf --output_tabular out_file.tabular --output_chimeras out_file_chimeras.tabular`

### FSD: Family Size Distribution of duplex sequencing tags
This tool provides a computationally very fast insight into the distribution of the family sizes of ALL tags from a Duplex Sequencing experiment (DS) and gives a first assessment of the distribution of PE-reads in families with 1 member up to >20 members. This information is very useful in early decision steps of the analysis parameters, such as the minimum number of PE-reads to build the single stranded consensus sequence (SSCS). Moreover, this tool can compare several datasets or different steps in the analysis pipeline to monitor data loss or gain (e.g families re-united with barcode correction tool from the Du Novo Analysis Pipeline). In an extension of this tool, each family is stratified into SSCS (ab/ba) and DSC and visualizes the allocation of DCSs respective to SSCS-ab and SSCS-ba. This is quite handy to better understand the relationship of SSCS to DCS per family and identify sources of bias (e.g. more SSCS to DCS in a particular family size, or more forward ab than reverse ba reads).

The tool can be used this tool via the command line with its default settings as the following:

`$ python2 fsd.py --inputFile1 tag_file.tabular --inputName1 tag_file.tabular --inputFile2 tag_file2.tabular --inputName2 tag_file2.tabular --inputFile3 tag_file3.tabular --inputName3 tag_file3.tabular --inputFile4 tag_file4.tabular --inputName4 tag_file4.tabular --log_axis --rel_freq --output_pdf out_file.pdf --output_tabular out_file.tabular`

### FSD regions: Family Size Distribution of user-specified regions in the reference genome
This tool provides a computationally very fast insight into the distribution of the family sizes of ALL tags from a Duplex Sequencing (DS) experiment that were aligned to different regions targeted in the reference genome.

The tool can be used this tool via the command line with its default settings as the following:

`$ python2 fsd_regions.py --inputFile tag_file.tabular --inputName1 tag_file.tabular --bamFile DCS.bam --rangesFile regions.bed --output_pdf out_file.pdf --output_tabular out_file.tabular`

### FSD Before/After: Family Size Distribution of duplex sequencing tags during Du Novo analysis
This tool will create a distribution of family sizes from tags of various steps of the [Du Novo Analysis Pipeline](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1039-4).

The tool can be used this tool via the command line with its default settings as the following:

`$ python2 fsd_beforevsafter.py --inputFile_SSCS tag_file.tabular --inputName1 tag_file.tabular --makeDCS DCS.fasta --afterTrimming DCS_trimmed.fasta --bamFile DCS.bam --output_pdf out_file.pdf --output_tabular out_file.tabular`



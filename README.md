# Tools for the quality control of duplex sequencing data

This tools allow a deeper insight into duplex sequencing data. They can be used within the [Galaxy platform](http://usegalaxy.org) under the section Du Novo, on your local Galaxy installation from the [toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/duplex_family_size_distribution) and on the command line. 

## Dependencies
This tools were build with Python 2.7.

## Usage
### TD: Tag distance analysis of duplex tags
Tags used in Duplex Sequencing (DS) are randomized barcodes, e.g 12 base pairs long. Since each DNA fragment is labeled by two tags at each end there are theoretically 4 to the power of (12+12) unique combinations. However, the input DNA in a typical DS experiment contains only ~1,000,000 molecules creating a large tag-to-input excess (4^24 ≫ 1,000,000). Because of such excess it is highly unlikely to tag distinct input DNA molecules with highly similar barcodes.

This tool calculates the number of nucleotide differences among tags (tag distance), also known as [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance). In this context the Hamming distance is simply the number of differences between two tags. The tool compares in a randomly selected subset of tags (default n=1000), the difference between each tag of the subset with the tags of the complete dataset. Each tag will differ by a certain number of nucleotides with the other tags; yet the tool uses the smallest difference observed with any other tag.

This tools expects a tabular file with the tags of all families, the family sizes and information about forward (ab) and reverse (ba) strands, which is created by the `Make Families` or `Correct Barcodes` of the Du Novo pipeline described in [Stoler *et al.* 2016]. A detailed describtion on how to produce this file in [Galaxy](http://usegalaxy.org). 

The tool can be used this tool via the command line with its default settings as the following:
$ python2 td.py --inputFile tag_file.tabular --inputName1 tag_file.tabular --sample_size 0 --subset_tag 0 --nproc 8 --rel_freq --minFS 1 --maxFS 0 --nr_above_bars --output_pdf out_file.pdf --output_tabular out_file.tabular --output_chimeras out_file_chimeras.tabular


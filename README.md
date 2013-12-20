# bcbio.variation.recall

Parallel merging, squaring off and ensemble calling for genomic variants.
Provide a general framework meant to combine multiple variant calls, either from
single individuals, batched family calls, or multiple approaches on the same
sample. Splits inputs based on shared genomic regions without variants, allowing
independent processing of smaller regions with variant calls. Handles:

- Merging multiple samples, called independently, into a single final VCF file.
- Squaring off multiple samples, called independently, by recalling at all
  identified genomic positions.
- Ensemble calling on single samples using inputs from multiple variant callers.

This is a work in progress.

## Usage

The executable `bcbio-variation-recall` bash script contains a ready to run jar
file. Pre-built distributions will be available. To create a development version
run `make` and the executable will be available in the `bin` directory. This
requires [leiningen].

[leiningen]: http://leiningen.org/

### Merge

    Merge multiple VCF files together, running in parallel over genomic regions.

    Usage: bcbio.variation.recall merge [options] out-file ref-file vcf-files

      out-file:  VCF (or bgzipped VCF) file to write merged output to
      ref-file:  FASTA format genome reference file
      vcf-files: VCF files to merge. Can be specified on the command line
                 or as a text file containing paths to files for processing

    Options:
      -c, --cores  Number of cores to use
      -h, --help

## Thank you

External software provides the underlying algorithms. This tool is a framework
for pulling them together. The following command line programs need to be on
your path:

- [freebayes][freebayes]
- [glia][glia]
- [vcflib][vcflib]
- [bedtools][bedtools]
- [bcftools (0.20+, with htslib)][bcftools]
- [samtools][samtools]

The [bcbio-nextgen][bcbio-nextgen] pipeline installs all this software automatically.

[bcbio-nextgen]: https://github.com/chapmanb/bcbio-nextgen
[bedtools]: http://bedtools.readthedocs.org/en/latest/
[vcflib]: https://github.com/ekg/vcflib
[bcftools]: https://github.com/samtools/bcftools
[freebayes]: https://github.com/ekg/freebayes
[glia]: https://github.com/ekg/glia
[samtools]: http://samtools.sourceforge.net/

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html

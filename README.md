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

External command line programs provide the underlying algorithms. This tool
is a framework for pulling them together. The following programs need to
be on your path:

- [bedtools][bedtools]
- [vcflib][vcflib]
- [bcftools (0.20+, with htslib)][bcftools]
- [freebayes][freebayes]
- [glia][glia]
- [samtools][samtools]

The [[bcbio-nextgen][bcbio-nextgen]] pipeline installs all this software automatically.

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

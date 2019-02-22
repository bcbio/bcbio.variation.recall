## 0.2.4 (22 Feb 2019)

- Fix platypus calling to correctly add contig headers to avoid GATK4 GatherVcf
  errors.

## 0.2.3 (12 Feb 2019)

- Fix wrapper script to avoid java in `/bin/java` without `JAVA_HOME` set.

## 0.2.2 (1 Feb 2019)

- Handle input to bcftools concat when all original files have no variants.
- Remove problematic INFO and FORMAT fields for FreeBayes. In complex decomposed, it
  retains more AD/AO/QA/GL fields then the headers allow, causing bcftools failures.
- Fix logging to correctly raise errors on failures.

## 0.2.1 (17 Jan 2019)

- Clarify filters on depth (DP) which fail on recent bcftools when present
  in both INFO and FORMAT key/values.

## 0.1.9 (3 August 2018)

- Remove usage of gatk-framework, which is no longer included in bcbio. Replace
  GATK usage with bcftools merge and concat.

## 0.1.8 (18 May 2018)

- Fix tabix index error for variants present at position 1 of contigs.

## 0.1.7 (15 August 2016)

- Avoid FreeBayes error when running `--variant-input` without an empty VCF.

## 0.1.6 (25 April 2016)

- Handle larger merge sizes with additional memory specifications up to 3000 or
  more samples.
- Use local temporary directories for GATK calls, avoiding filling up temporary
  space. Thanks to Roman Valls Guimerà.
- Correctly parse sample names from large numbers of inputs without causing
  overflow error.
- Do not silently exit if running into Java memory errors.
- Automatically convert ensemble output files to bgzipped output.

## 0.1.5 (15 April 2016)

- Add genotype qualities (GQ) to output of FreeBayes recalls.
- Use samtools instead of scramble for CRAM integration.

## 0.1.4 (15 October 2015)

- Annotate ensemble variant calls with names of callers supporting an ensemble variant.

## 0.1.3 (13 October 2015)

- Support latest GATK and htsjdk to handle issues with gVCFs on GATK 3.4 and
  different command line options in gatk-framework.

## 0.1.2 (25 April 2015)

- Correctly handle copies of single files in transactional directories.

## 0.1.1 (15 April 2015)

- Correctly intersect files when one of the inputs is empty. Thanks to Daryl
  Waggott.

## 0.1.0 (7 April 2015)

- Fix command line sorting issues with X/Y chromosomes. Thanks to Lorena
  Pantano.
- Separate command line help library into bcbio.run for use outside of
  bcbio.variation.recall.

## 0.0.9 (2 April 2015)

- FreeBayes: update with options for latest validated version 0.9.21-7

## 0.0.8 (20 March 2015)

- Avoid copy errors when merging recalls with only a single region.

## 0.0.7 (23 February 2015)

- Fix errors encountered on large runs, adding additional temporary directory
  usage during copy and spit commands to prevent partial files.

## 0.0.6 (30 January 2015)

- Correctly handle identically named input files for recalling.
- Avoid filtering low quality reference calls when recalling with Platypus.
- Correctly support vcffixup stdin in recent vcflib.
- Support new sambamba index changes.

## 0.0.5 (1 December 2014)

- Ensemble calling: handle variant inputs with multiple identical calls.

## 0.0.4 (28 October 2014)

- FreeBayes recalling fixes to support GATK compatibility: remove
  duplicate alternative alleles and ceil very low FreeBayes quality scores.

## 0.0.3 (26 October 2014)

- Fix FreeBayes recalling/squaring off to better distinguish reference and no
  call. Avoid filtering low-quality reference positions and assigning as no
  call while still removing problem variants.
- Fix ensemble calling for cases with closely spaced variants containing
  identical ref/alts.

## 0.0.2 (29 September 2014)

- Initial implementation of intersection based ensemble calling. Handle
  germline and tumor/normal cases.
- Use new picard and htsjdk libraries for variant manipulation.

## 0.0.1 (20 September 2014)

- Initial release with support for parallel genome merging and squaring
  off/joint calling with FreeBayes, Platypus and samtools.

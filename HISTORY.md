## 0.0.6 (in progress)

- Correctly handle identically named input files for recalling.
- Avoid filtering low quality reference calls when recalling with Platypus.

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

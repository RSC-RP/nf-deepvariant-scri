include { SOMALIER_EXTRACT } from '../modules/nf-core/somalier/extract'
include { SOMALIER_RELATE } from '../modules/nf-core/somalier/relate'

workflow somalier {
    take:
    input_ch
    bams_ch
    fasta
    fai
    sites

    main:
    SOMALIER_EXTRACT(bams_ch, fasta, fai, sites)
    // Need sites file as param
    // Need Python module to build the pedigree file from the input CSV
    // Need to collect extract files
    // SOMALIER_RELATE()

    emit:
    extract = SOMALIER_EXTRACT.out.extract
}
include { SOMALIER_EXTRACT } from '../modules/nf-core/somalier/extract'
include { SOMALIER_RELATE } from '../modules/nf-core/somalier/relate'
include { MAKE_PEDIGREE } from '../modules/local/python/make_pedigree'

workflow somalier {
    take:
    input_ch
    bams_ch
    fasta
    fai
    sites
    cohort_name

    main:
    SOMALIER_EXTRACT(bams_ch, fasta, fai, sites)
    // Get sample name from BAM for pedigree
    sample_lookup = SOMALIER_EXTRACT.out.extract
        .map{ meta, extract ->
            meta + [ sample_name: extract.baseName ]
        }
        .collectFile(name: 'sample_lookup.txt'){ meta ->
            "${meta.id}\t${meta.sample_name}\n"
        }
    // Build the pedigree file from the input CSV
    MAKE_PEDIGREE(input_ch, sample_lookup, cohort_name)
    ped_ch = MAKE_PEDIGREE.out
        .map{ ped -> [[id: cohort_name], ped]}
    // Collect extract files, join to pedigree
    relate_input = SOMALIER_EXTRACT.out.extract
        .map{ _meta, extract -> extract }
        .collect()
        .map{ extract_list -> [[id: cohort_name], extract_list]}
        .join(ped_ch)
    
    SOMALIER_RELATE(relate_input, [])

    emit:
    extract = SOMALIER_EXTRACT.out.extract
    lookup = sample_lookup
    pedigree = MAKE_PEDIGREE.out
    html = SOMALIER_RELATE.out.html
    pairs_tsv = SOMALIER_RELATE.out.pairs_tsv
    samples_tsv = SOMALIER_RELATE.out.samples_tsv
}

include { MAKE_EXAMPLES_TRIO } from '../modules/local/deepvariant/make_examples_trio/main'
include { CALL_VARIANTS_TRIO } from '../modules/local/deepvariant/call_variants_trio/main'
include { MAKE_EXAMPLES_SINGLE } from '../modules/local/deepvariant/make_examples_single/main'
include { CALL_VARIANTS_SINGLE } from '../modules/local/deepvariant/call_variants_single/main'

// Generalized subworkflows to run make_examples and call_variants
// Any individual
workflow deepvariant {
    take:
        bam_ch // tuple with meta, bam, and index. Meta has proband_sex, proband_id, mother_id, father_id, and id
        fasta
        fai
        par_bed
        x_only
        y_only
        test_bams
        genome_ver
        chromnames
        deepvar_model
        nshards
    main:
        MAKE_EXAMPLES_SINGLE(bam_ch, fasta, fai, par_bed, x_only, y_only, test_bams, genome_ver, chromnames)
        MAKE_EXAMPLES_SINGLE.out.proband_tfrecord
            .join(MAKE_EXAMPLES_SINGLE.out.example_info)
            .set{ single_tfrecords }
        CALL_VARIANTS_SINGLE(single_tfrecords, deepvar_model, nshards)
        CALL_VARIANTS_SINGLE.out
    emit:
        CALL_VARIANTS_SINGLE.out
}

// Any family
workflow deeptrio {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        test_bams
        genome_ver
        chromnames
        deepvar_model
        nshards
    main:
        MAKE_EXAMPLES_TRIO(bam_ch, fasta, fai, test_bams, genome_ver, chromnames, deepvar_model)
        MAKE_EXAMPLES_TRIO.out.proband_tfrecord
            .map{ meta, me, gvcf -> [[proband_id: meta.proband_id], meta, me, gvcf] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ _proband, meta, me, gvcf, ei -> [meta, me, gvcf, ei] }
            .set{ all_proband_tfrecords }
        MAKE_EXAMPLES_TRIO.out.father_tfrecord
            .map{ meta, me, gvcf -> [[proband_id: meta.proband_id], meta, me, gvcf] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ _proband, meta, me, gvcf, ei -> [meta, me, gvcf, ei] }
            .set{ all_father_tfrecords }
        MAKE_EXAMPLES_TRIO.out.mother_tfrecord
            .map{ meta, me, gvcf -> [[proband_id: meta.proband_id], meta, me, gvcf] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ _proband, meta, me, gvcf, ei -> [meta, me, gvcf, ei] }
            .set{ all_mother_tfrecords }
        all_proband_tfrecords
            .concat(all_father_tfrecords, all_mother_tfrecords)
            .filter{ meta, _me, _gvcf, _ei -> meta.id != "" }
            .set{ all_me_tfrecords }
        CALL_VARIANTS_TRIO(all_me_tfrecords, deepvar_model, nshards)
        CALL_VARIANTS_TRIO.out
    emit:
        CALL_VARIANTS_TRIO.out
}

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
    main:
        MAKE_EXAMPLES_SINGLE(bam_ch, fasta, fai, par_bed, x_only, y_only)
        MAKE_EXAMPLES_SINGLE.out.proband_tfrecord
            .join(MAKE_EXAMPLES_SINGLE.out.example_info)
            .set{ single_tfrecords }
        CALL_VARIANTS_SINGLE(single_tfrecords)
        CALL_VARIANTS_SINGLE.out
            .set{ cv_out }
    emit:
        cv_out
}

// Any family
workflow deeptrio {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
    main:
        MAKE_EXAMPLES_TRIO(bam_ch, fasta, fai)
        MAKE_EXAMPLES_TRIO.out.proband_tfrecord
            .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ [it[1], it[2], it[3], it[4]] }
            .set{ all_proband_tfrecords }
        MAKE_EXAMPLES_TRIO.out.father_tfrecord
            .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ [it[1], it[2], it[3], it[4]] }
            .set{ all_father_tfrecords }
        MAKE_EXAMPLES_TRIO.out.mother_tfrecord
            .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
            .join(MAKE_EXAMPLES_TRIO.out.example_info)
            .map{ [it[1], it[2], it[3], it[4]] }
            .set{ all_mother_tfrecords }
        all_proband_tfrecords
            .concat(all_father_tfrecords, all_mother_tfrecords)
            .filter{ it[0].id != "" }
            .set{ all_me_tfrecords }
        CALL_VARIANTS_TRIO(all_me_tfrecords)
        CALL_VARIANTS_TRIO.out
            .set{ cv_out }
    emit:
        cv_out
}

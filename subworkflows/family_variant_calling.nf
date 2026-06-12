include { deeptrio } from "./make_examples_call_variants"
include { deepvariant } from "./make_examples_call_variants"
include { deepvariant as deepvariant1 } from "./make_examples_call_variants"

// subworkflows for different family structures
// Singletons
workflow mecv_single {
    take:
        bam_ch // tuple with meta, bam, and index. Meta has proband_sex, proband_id, mother_id="", father_id=""
        fasta
        fai
        par_bed
    main:
        bam_ch
            .map{ [it[0] + [id: it[0].proband_id], it[1], it[2]] }
            .set{ bam_ch }
        deepvariant(bam_ch, fasta, fai, par_bed, false, false)
        deepvariant.out
            .set{ cv_out }
    emit:
        cv_out
}

// Male trio; DeepTrio misses father's X chromosome
workflow mecv_maletrio {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        par_bed
    main:
        // family
        deeptrio(bam_ch, fasta, fai)
        // just the father
        bam_ch
            .map{ [it[0] + [id: it[0].father_id], it[1][1], it[2][1]] }
            .set{ bam_dad }
        deepvariant(bam_dad, fasta, fai, par_bed, true, false) // last two booleans specify X chrom
        deepvariant.out
           .map{ [[id: it[0].id, proband_id: it[0].proband_id, role: "parent", sex: it[0].sex], it[1], it[2]] }
           .set { cv_dad }
        deeptrio.out
           .join(cv_dad, remainder: true)
           .map{
             if(it.size == 5){
                [it[0], it[1] + it[3], it[2] + it[4]] // combines output from deeptrio and deepvariant
             }
             else{
                [it[0], it[1], it[2]]
             }
           }
           .set{ cv_out }
    emit:
        cv_out
}

// Female trio, or female plus father; DeepTrio misses Y chromosome.
workflow mecv_femaletrio_dadduo {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        par_bed
    main:
        // family
        deeptrio(bam_ch, fasta, fai)
        // just the father
        bam_ch
            .map{ [it[0] + [id: it[0].father_id], it[1][1], it[2][1]] }
            .set{ bam_dad }
        deepvariant(bam_dad, fasta, fai, par_bed, false, true) // last two booleans specify Y chrom
        deepvariant.out
           .map{ [[id: it[0].id, proband_id: it[0].proband_id, role: "parent", sex: it[0].sex], it[1], it[2]] }
           .set { cv_dad }
        deeptrio.out
           .join(cv_dad, remainder: true)
           .map{
             if(it.size == 5){
                [it[0], it[1] + it[3], it[2] + it[4]] // combines output from deeptrio and deepvariant
             }
             else{
                [it[0], it[1], it[2]]
             }
           }
           .set{ cv_out }
    emit:
        cv_out
}

// Duo of male proband and mother; DeepTrio misses proband Y chromosome.
workflow mecv_malemomduo {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        par_bed
    main:
        // family
        deeptrio(bam_ch, fasta, fai)
        // just the proband
        bam_ch
            .map{ [it[0] + [id: it[0].proband_id], it[1][0], it[2][0]] }
            .set{ bam_pro }
        deepvariant(bam_pro, fasta, fai, par_bed, false, true) // last two booleans specify Y chrom
        deepvariant.out
           .map{ [[id: it[0].id, proband_id: it[0].proband_id, role: "child", sex: it[0].sex], it[1], it[2]] }
           .set { cv_pro }
        deeptrio.out
           .join(cv_pro, remainder: true)
           .map{
             if(it.size == 5){
                [it[0], it[1] + it[3], it[2] + it[4]] // combines output from deeptrio and deepvariant
             }
             else{
                [it[0], it[1], it[2]]
             }
           }
           .set{ cv_out }
    emit:
        cv_out
}

// Duo of female proband and mother; nothing missed by DeepTrio.
workflow mecv_femalemomduo {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        par_bed
    main:
        // family
        deeptrio(bam_ch, fasta, fai)
           .set{ cv_out }
    emit:
        cv_out
}

// Male duo with father; DeepTrio misses X chromosome for both.
workflow mecv_maledadduo {
    take:
        bam_ch // tuple with meta, bams, and indices. Meta has proband_sex, proband_id, mother_id, father_id
        fasta
        fai
        par_bed
    main:
        // family
        deeptrio(bam_ch, fasta, fai)
        // just the father
        bam_ch
            .map{ [it[0] + [id: it[0].father_id], it[1][1], it[2][1]] }
            .set{ bam_dad }
        deepvariant(bam_dad, fasta, fai, par_bed, true, false) // last two booleans specify X chrom
        deepvariant.out
           .map{ [[id: it[0].id, proband_id: it[0].proband_id, role: "parent", sex: it[0].sex], it[1], it[2]] }
           .set { cv_dad }
        // just the proband
        bam_ch
            .map{ [it[0] + [id: it[0].proband_id], it[1][0], it[2][0]] }
            .set{ bam_pro }
        deepvariant1(bam_pro, fasta, fai, par_bed, true, false) // last two booleans specify X chrom
        deepvariant1.out
           .map{ [[id: it[0].id, proband_id: it[0].proband_id, role: "child", sex: it[0].sex], it[1], it[2]] }
           .set { cv_pro }
        deeptrio.out
           .join(cv_dad.concat(cv_pro), remainder: true)
           .map{
             if(it.size == 5){
                [it[0], it[1] + it[3], it[2] + it[4]] // combines output from deeptrio and deepvariant
             }
             else{
                [it[0], it[1], it[2]]
             }
           }
           .set{ cv_out }
    emit:
        cv_out
}

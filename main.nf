include { mecv_single } from "../subworkflows/family_variant_calling"
include { mecv_maletrio } from "../subworkflows/family_variant_calling"
include { mecv_femaletrio_dadduo } from "../subworkflows/family_variant_calling"
include { mecv_malemomduo } from "../subworkflows/family_variant_calling"
include { mecv_femalemomduo } from "../subworkflows/family_variant_calling"
include { mecv_maledadduo } from "../subworkflows/family_variant_calling"
include { POSTPROCESS_VARIANTS } from '../modules/local/deepvariant/postprocess_variants/main'
include { GLNEXUS } from '../modules/local/glnexus/main'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'

// workflow for variant calling on trios, duos, or singletons
workflow CALL_TRIOS {
    // Reference sequence for alignment and genotype calling, which may be different from Ensembl.
    Channel.fromPath(file(params.fasta_bams, checkIfExists: true))
        .collect() // allows the reference to be used with multiple input VCFs
        .set{ fasta_bams }
    fai_file2 = file("${params.fasta_bams}.fai")
    if( fai_file2.exists() ){
        Channel.fromPath(fai_file2)
            .collect()
            .set{ fai_bams }
    }
    else {
        // Run Samtools Faidx
        SAMTOOLS_FAIDX(fasta_bams)
        SAMTOOLS_FAIDX.out.fai
            .collect()
            .set{ fai_bams }
    }

    // Read in sample list
    Channel.fromPath(file(params.sample_bams, checkIfExists: true))
        .splitCsv(header: true, sep: ',')
        .map{
            row ->
                [ [proband_sex: row.proband_sex, proband_id: row.proband_id, father_id: row.father_id, mother_id: row.mother_id], // meta for family
                  [[id: row.proband_id, proband_id: row.proband_id, ord: 0],
                   [id: row.father_id, proband_id: row.proband_id, ord: 1],
                   [id: row.mother_id, proband_id: row.proband_id, ord: 2]], // meta for individuals
                  [row.proband_bam, row.father_bam, row.mother_bam],
                  [row.proband_index, row.father_index, row.mother_index]
                ]
            }
        .transpose()
        .filter{ it[1].id != "" & it[2] != "" }
        .map{ [it[0], it[1], file(it[2], checkIfExists: true), it[3]] }
        .set{ input_ch }
    
    // check sex specification
    input_ch
        .map{ assert ["Male", "male", "M", "Female", "female", "F"].contains(it[0].proband_sex) }

    // Get a channel just for looking up metadata again
    input_ch
        .map{ [it[1], it[0]] }
        .set{ meta_lookup }

    // Split off into samples input as bams vs fastqs
    input_ch
        .branch{ it ->
            fastq: it[2].name.endsWith(".fq.gz") | it[2].name.endsWith(".fastq.gz") | it[2].name.endsWith(".fq") | it[2].name.endsWith(".fastq")
            bam: true // get everything else
        }
        .set{ input_ch }

    // Align FASTQ files
    input_ch.fastq
        .map{ [it[1], it[2]] } // get metadata and fastq
        .set{ fastq_ch }
    fasta_bams
        .map{ [[id: it.simpleName], it]}
        .collect()
        .set{ fasta_for_bwa }
    Channel.fromPath(params.bwa_index)
        .map{ [[id: it.simpleName], it]}
        .collect()
        .set{ bwa_index }
    BWA_MEM(fastq_ch, bwa_index, fasta_for_bwa, true)

    // Combine new aligned bams with any existing bams
    BWA_MEM.out.bam
        .join(meta_lookup)
        .map{ [it[2], it[0], it[1], ""]}
        .concat(input_ch.bam)
        .set{ bam_ch }
    
    // Index any BAMS with missing index
    bam_ch
        .branch{ it ->
            needs_index: it[3] == ""
            indexed: true
        }
        .set{ bam_ch }
    bam_ch.needs_index
        .map{ [it[1], it[2]] }
        .set{ bams_to_index }
    SAMTOOLS_INDEX(bams_to_index)

    // Join back in with indexed BAMs
    bam_ch.indexed
        .map{ [it[0], it[1], it[2], file(it[3], checkIfExists: true)] }
        .set{ bams_indexed }
    bams_to_index
       .join(SAMTOOLS_INDEX.out.bai)
       .join(meta_lookup)
       .map{ [it[3], it[0], it[1], it[2]] }
       .concat(bams_indexed)
       .set{ bam_ch }

    // Get BAMS back into proband-father-mother order and split families from singletons
    bam_ch
        .map{ [ it[0], it[1].plus([bam: it[2], index: it[3]])] }
        .groupTuple(sort: { it.ord } )
        .map{ [ it[0], it[1].bam, it[1].index ] }
        .branch{ it ->
            single: it[0].father_id == "" & it[0].mother_id == ""
            maletrio: ["Male", "male", "M"].contains(it[0].proband_sex) & it[0].father_id != "" & it[0].mother_id != ""
            femaleWdad: ["Female", "female", "F"].contains(it[0].proband_sex) & it[0].father_id != ""
            malemomduo: ["Male", "male", "M"].contains(it[0].proband_sex) & it[0].father_id == "" & it[0].mother_id != ""
            femalemomduo: ["Female", "female", "F"].contains(it[0].proband_sex) & it[0].father_id == "" & it[0].mother_id != ""
            maledadduo: ["Male", "male", "M"].contains(it[0].proband_sex) & it[0].father_id != "" & it[0].mother_id == ""
        }
        .set{bam_ch}
    
    // Variant calling on families
    Channel.fromPath(file(params.par_bed, checkIfExists: true))
        .collect()
        .set{ par_bed }
    mecv_maletrio(bam_ch.maletrio, fasta_bams, fai_bams, par_bed)
    mecv_femaletrio_dadduo(bam_ch.femaleWdad, fasta_bams, fai_bams, par_bed)
    mecv_malemomduo(bam_ch.malemomduo, fasta_bams, fai_bams, par_bed)
    mecv_femalemomduo(bam_ch.femalemomduo, fasta_bams, fai_bams, par_bed)
    mecv_maledadduo(bam_ch.maledadduo, fasta_bams, fai_bams, par_bed)
    // Variant calling on singletons
    mecv_single(bam_ch.single, fasta_bams, fai_bams, par_bed)

    // Join together families and singletons
    mecv_single.out
        .concat(mecv_maletrio.out,
                mecv_femaletrio_dadduo.out,
                mecv_malemomduo.out,
                mecv_femalemomduo.out,
                mecv_maledadduo.out)
        .set{ call_variants_out }

    // Postprocess both together
    POSTPROCESS_VARIANTS(call_variants_out, fasta_bams, fai_bams, par_bed)
    POSTPROCESS_VARIANTS.out
        .map{ [it[1], it[2]] }
        .collect()
        .set{ all_ind_vcfs }

    // gl_nexus for joint VCF
    GLNEXUS(all_ind_vcfs, params.cohort_name) 
}

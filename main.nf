include { make_bams } from './subworkflows/make_bams'
include { mecv_single } from "./subworkflows/family_variant_calling"
include { mecv_maletrio } from "./subworkflows/family_variant_calling"
include { mecv_femaletrio_dadduo } from "./subworkflows/family_variant_calling"
include { mecv_malemomduo } from "./subworkflows/family_variant_calling"
include { mecv_femalemomduo } from "./subworkflows/family_variant_calling"
include { mecv_maledadduo } from "./subworkflows/family_variant_calling"
include { POSTPROCESS_VARIANTS } from './modules/local/deepvariant/postprocess_variants/main'
include { GLNEXUS } from './modules/local/glnexus/main'
include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'

params {
    sample_bams: Path // Input CSV listing BAM paths by trio/duo/singleton
    deepvar_model: String // "WGS" or "WES"
    genome_ver: String // "hg19" or "hg38"
    chromnames: String // "g1k", "ensembl", or "ucsc"
    cohort_name: String // Name used in prefix of output files
    make_examples_nshards: Integer = 32 // Parallelization for DeepVariant "make examples" step
    test_bams: Boolean = false // If true, only runs a tiny portion of the genome
    fasta_bams: Path // Reference genome that was used to align BAMs. Check BAM header if unsure. If starting from FASTQ, align to this reference.
    bwa_index: Path // folder containing the BWA index corresponding to fasta_bams
    par_bed: Path // haploid regions for DeepVariant
    glnexus_filter: Boolean = true // Whether to use DeepVariant_unfiltered (false) or DeepVariantWES or DeepVariantWGS (true)
}

// workflow for variant calling on trios, duos, or singletons
workflow {
    main:
    assert params.deepvar_model == "WGS" | params.deepvar_model == "WES"
    assert params.genome_ver == "hg19" | params.genome_ver == "hg38"
    assert params.chromnames == "g1k" | params.chromnames == "ensembl" | params.chromnames == "ucsc"

    // Reference sequence for alignment and genotype calling, which may be different from Ensembl.
    channel.fromPath(file(params.fasta_bams, checkIfExists: true))
        .map{ fasta -> [[id: fasta.simpleName], fasta]}
        .collect() // allows the reference to be used with multiple input VCFs
        .set{ fasta_bams }
    fai_file2 = file("${params.fasta_bams}.fai")
    if( fai_file2.exists() ){
        channel.fromPath(fai_file2)
            .map{ fai -> [[id: fai.simpleName], fai]}
            .collect()
            .set{ fai_bams }
    }
    else {
        // Run Samtools Faidx
        SAMTOOLS_FAIDX(fasta_bams, false)
        SAMTOOLS_FAIDX.out.fai
            .collect()
            .set{ fai_bams }
    }

    channel.fromPath(params.bwa_index)
        .map{ index -> [[id: index.simpleName], index]}
        .collect()
        .set{ bwa_index }

    // Read in sample list
    channel.fromPath(file(params.sample_bams, checkIfExists: true))
        .set{ input_ch }
    
    // Align and index if necessary
    make_bams(input_ch, fasta_bams, bwa_index)
    bam_ch = make_bams.out
    
    // Variant calling on families
    channel.fromPath(file(params.par_bed, checkIfExists: true))
        .collect()
        .set{ par_bed }
    mecv_maletrio(bam_ch.maletrio, fasta_bams, fai_bams, par_bed, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)
    mecv_femaletrio_dadduo(bam_ch.femaleWdad, fasta_bams, fai_bams, par_bed, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)
    mecv_malemomduo(bam_ch.malemomduo, fasta_bams, fai_bams, par_bed, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)
    mecv_femalemomduo(bam_ch.femalemomduo, fasta_bams, fai_bams, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)
    mecv_maledadduo(bam_ch.maledadduo, fasta_bams, fai_bams, par_bed, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)
    // Variant calling on singletons
    mecv_single(bam_ch.single, fasta_bams, fai_bams, par_bed, params.test_bams, params.genome_ver, params.chromnames, params.deepvar_model, params.make_examples_nshards)

    // Join together families and singletons
    mecv_single.out
        .concat(mecv_maletrio.out,
                mecv_femaletrio_dadduo.out,
                mecv_malemomduo.out,
                mecv_femalemomduo.out,
                mecv_maledadduo.out)
        .set{ call_variants_out }

    // Postprocess both together
    POSTPROCESS_VARIANTS(call_variants_out, fasta_bams, fai_bams, par_bed, params.chromnames)
    POSTPROCESS_VARIANTS.out
        .map{ _meta, gvcf, tbi -> [gvcf, tbi] }
        .collect()
        .set{ all_ind_vcfs }

    // gl_nexus for joint VCF
    GLNEXUS(all_ind_vcfs, params.cohort_name, params.glnexus_filter, params.deepvar_model)

    // OUTPUT
    publish:
    gvcfs = POSTPROCESS_VARIANTS.out
    bcf = GLNEXUS.out
}

output {
    gvcfs: Channel<Path> {
        path "gvcfs"
        index {
            path 'gvcfs.json'
        }
    }
    bcf: Channel<Path> {
        path '.'
    }
}

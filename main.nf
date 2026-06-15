include { mecv_single } from "./subworkflows/family_variant_calling"
include { mecv_maletrio } from "./subworkflows/family_variant_calling"
include { mecv_femaletrio_dadduo } from "./subworkflows/family_variant_calling"
include { mecv_malemomduo } from "./subworkflows/family_variant_calling"
include { mecv_femalemomduo } from "./subworkflows/family_variant_calling"
include { mecv_maledadduo } from "./subworkflows/family_variant_calling"
include { POSTPROCESS_VARIANTS } from './modules/local/deepvariant/postprocess_variants/main'
include { GLNEXUS } from './modules/local/glnexus/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
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

    // Read in sample list
    channel.fromPath(file(params.sample_bams, checkIfExists: true))
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
        .filter{ _meta_family, meta_individual, bam, _bai -> meta_individual.id != "" & bam != "" }
        .map{ meta_family, meta_individual, bam, bai -> [meta_family, meta_individual, file(bam, checkIfExists: true), bai] }
        .set{ input_ch }
    
    // check sex specification
    input_ch
        .map{ _meta_family, meta_individual, _bam, _bai ->
        assert ["Male", "male", "M", "Female", "female", "F"].contains( meta_individual.proband_sex)
        }

    // Get a channel just for looking up metadata again
    input_ch
        .map{ meta_family, meta_individual, _bam, _bai -> [meta_individual, meta_family] }
        .set{ meta_lookup }

    // Split off into samples input as bams vs fastqs
    input_ch
        .branch{ _meta_family, _meta_individual, bam, _bai ->
            fastq: bam.name.endsWith(".fq.gz") | bam.name.endsWith(".fastq.gz") | bam.name.endsWith(".fq") | bam.name.endsWith(".fastq")
            bam: true // get everything else
        }
        .set{ input_ch }

    // Align FASTQ files
    input_ch.fastq
        .map{ _meta_family, meta_individual, fastq, _empty -> [meta_individual, fastq] } // get metadata and fastq
        .set{ fastq_ch }
    channel.fromPath(params.bwa_index)
        .map{ index -> [[id: index.simpleName], index]}
        .collect()
        .set{ bwa_index }
    BWA_MEM(fastq_ch, bwa_index, fasta_bams, true)

    // Combine new aligned bams with any existing bams
    BWA_MEM.out.bam
        .join(meta_lookup)
        .map{ meta_individual, bam, meta_family -> [meta_family, meta_individual, bam, ""]}
        .concat(input_ch.bam)
        .set{ bam_ch }
    
    // Index any BAMS with missing index
    bam_ch
        .branch{ _meta_family, _meta_individual, _bam, bai ->
            needs_index: bai == ""
            indexed: true
        }
        .set{ bam_ch }
    bam_ch.needs_index
        .map{ _meta_family, meta_individual, bam, _bai -> [meta_individual, bam] }
        .set{ bams_to_index }
    SAMTOOLS_INDEX(bams_to_index)

    // Join back in with indexed BAMs
    bam_ch.indexed
        .map{ meta_family, meta_individual, bam, bai -> [meta_family, meta_individual, bam, file(bai, checkIfExists: true)] }
        .set{ bams_indexed }
    bams_to_index
       .join(SAMTOOLS_INDEX.out.bai)
       .join(meta_lookup)
       .map{ meta_individual, bam, bai, meta_family -> [meta_family, meta_individual, bam, bai] }
       .concat(bams_indexed)
       .set{ bam_ch }

    // Get BAMS back into proband-father-mother order and split families from singletons
    bam_ch
        .map{ meta_family, meta_individual, bam, bai ->
            [ meta_family, meta_individual.plus([bam: bam, index: bai])]
        }
        .groupBy()
        .map{ meta_family, meta_ind_bams ->
            def mylist = meta_ind_bams.toList().sort{ it -> it.ord }
            [
                meta_family,
                mylist.collect{ it -> it.bam },
                mylist.collect{ it -> it.index }
            ]
        }
        .branch{ meta_family, _bams, _bais ->
            single: meta_family.father_id == "" & meta_family.mother_id == ""
            maletrio: ["Male", "male", "M"].contains(meta_family.proband_sex) & meta_family.father_id != "" & meta_family.mother_id != ""
            femaleWdad: ["Female", "female", "F"].contains(meta_family.proband_sex) & meta_family.father_id != ""
            malemomduo: ["Male", "male", "M"].contains(meta_family.proband_sex) & meta_family.father_id == "" & meta_family.mother_id != ""
            femalemomduo: ["Female", "female", "F"].contains(meta_family.proband_sex) & meta_family.father_id == "" & meta_family.mother_id != ""
            maledadduo: ["Male", "male", "M"].contains(meta_family.proband_sex) & meta_family.father_id != "" & meta_family.mother_id == ""
        }
        .set{bam_ch}
    
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
}

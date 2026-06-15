include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

// Subworkflow to set up BAMs channel from sample sheet
workflow make_bams {
    take:
        input_ch // sample sheet, unprocessed
        fasta // reference genome fasta
        bwa_index // BWA index

    main:
    input_ch = input_ch
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
    
    // check sex specification
    input_ch
        .map{ _meta_family, meta_individual, _bam, _bai ->
        assert ["Male", "male", "M", "Female", "female", "F"].contains( meta_individual.proband_sex)
        }

    // Get a channel just for looking up metadata again
    meta_lookup = input_ch
        .map{ meta_family, meta_individual, _bam, _bai -> [meta_individual, meta_family] }

    // Split off into samples input as bams vs fastqs
    input_ch = input_ch
        .branch{ _meta_family, _meta_individual, bam, _bai ->
            fastq: bam.name.endsWith(".fq.gz") | bam.name.endsWith(".fastq.gz") | bam.name.endsWith(".fq") | bam.name.endsWith(".fastq")
            bam: true // get everything else
        }
    
    // Align FASTQ files
    fastq_ch = input_ch.fastq
        .map{ _meta_family, meta_individual, fastq, _empty -> [meta_individual, fastq] } // get metadata and fastq

    BWA_MEM(fastq_ch, bwa_index, fasta, true)

    // Combine new aligned bams with any existing bams
    bam_ch = BWA_MEM.out.bam
        .join(meta_lookup)
        .map{ meta_individual, bam, meta_family -> [meta_family, meta_individual, bam, ""]}
        .concat(input_ch.bam)
        .branch{ _meta_family, _meta_individual, _bam, bai ->
            needs_index: bai == ""
            indexed: true
        }
    
    // Index any BAMs that need it
    bams_to_index = bam_ch.needs_index
        .map{ _meta_family, meta_individual, bam, _bai -> [meta_individual, bam] }
    SAMTOOLS_INDEX(bams_to_index)

    // Join back in with indexed BAMs
    bams_indexed = bam_ch.indexed
        .map{ meta_family, meta_individual, bam, bai -> [meta_family, meta_individual, bam, file(bai, checkIfExists: true)] }
    bam_ch = bams_to_index
       .join(SAMTOOLS_INDEX.out.bai)
       .join(meta_lookup)
       .map{ meta_individual, bam, bai, meta_family -> [meta_family, meta_individual, bam, bai] }
       .concat(bams_indexed)

    emit:
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
}
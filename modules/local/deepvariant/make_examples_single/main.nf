process MAKE_EXAMPLES_SINGLE {
    tag "$meta.id"
    container "docker://google/deepvariant:1.8.0-gpu"
    label "make_examples"

    input:
    tuple val(meta), path(bam), path(bai) // meta has proband_sex, proband_id, father_id, mother_id, id
    path(fasta_bams)
    path(fai_bams)
    path(par_bed)
    val(x_only) // boolean
    val(y_only) // boolean

    output:
    tuple val(meta2), path("make_examples*.tfrecord*.gz"), path("gvcf*.tfrecord*.gz"), emit: proband_tfrecord
    tuple val(meta2), path("make_examples*.tfrecord*.example_info.json"), emit: example_info

    script:
    def args = task.ext.args ?: ''

    assert !(x_only & y_only)

    is_male = meta.id == meta.father_id || meta.proband_sex == "Male" || meta.proband_sex == "male" || meta.proband_sex == "M"
    if(!is_male){
        assert meta.id == meta.mother_id || meta.proband_sex == "Female" || meta.proband_sex == "female" || meta.proband_sex == "F"
    }
    if(params.test_bams){
        assert meta.proband_id.startsWith("HG00") && bam.size() < 50000000
    }
    def proband_id = meta.id

    meta2 = meta + [sex: is_male ? "Male" : "Female"]

    // Set start and end points for X and Y
    if(params.annovar_buildver == "hg19"){
        startX = 2734540
        endX = 154997473
        startY = 2649521
        endY = 59034049
    }
    if(params.annovar_buildver == "hg38"){
        startX = 2781480
        endX = 155701382
        startY = 2781480
        endY = 56887902
    }

    // Chromosome prefix
    if(params.chromnames == "g1k" || params.chromnames == "ensembl"){
        pr = ""
    }
    if(params.chromnames == "ucsc"){
        pr = "chr"
    }

    // Set up code that will be the same for every run of make_examples
    mecmd = ["time seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples",
             "--mode calling --ref ${fasta_bams}",
             "--channel_list=read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size",
             "--pileup_image_height 100 --task {} --reads=${bam} --sample_name ${proband_id} $args"].join(' ').trim()
    if(is_male){
        mecmd = mecmd + " --haploid_contigs=\"${pr}X,${pr}Y\" --par_regions_bed=\"${par_bed}\""
    }

    if(x_only){
        if(params.test_bams){
            regions = "--regions ${pr}X:122530000-122550000"
        }
        else{
            regions = "--regions ${pr}X:${startX}-${endX}"
        }
        chunk = '2' // See make_examples_trio, male trio example for which region each chunk covers.
    }
    else if(y_only){
        regions = "--regions ${pr}Y:${startY}-${endY}"
        chunk = '4'
    }
    else if(params.test_bams){
        // Tiny example region of the genome for demo
        regions = "--regions \"${pr}2:179300000-179600000 ${pr}16:450000-470000 ${pr}X:122530000-122550000\""
        chunk = '1'
    }
    else if(is_male){
        regions = ["--regions",
                   "\"${pr}1 ${pr}2 ${pr}3 ${pr}4 ${pr}5 ${pr}6 ${pr}7 ${pr}8 ${pr}9 ${pr}10",
                   "${pr}11 ${pr}12 ${pr}13 ${pr}14 ${pr}15 ${pr}16 ${pr}17 ${pr}18 ${pr}19",
                   "${pr}20 ${pr}21 ${pr}22 ${pr}X ${pr}Y\""].join(' ')
        chunk = '1'
    }
    else{
        regions = ["--regions",
                   "\"${pr}1 ${pr}2 ${pr}3 ${pr}4 ${pr}5 ${pr}6 ${pr}7 ${pr}8 ${pr}9 ${pr}10",
                   "${pr}11 ${pr}12 ${pr}13 ${pr}14 ${pr}15 ${pr}16 ${pr}17 ${pr}18 ${pr}19",
                   "${pr}20 ${pr}21 ${pr}22 ${pr}X\""].join(' ')
        chunk = '1'
    }

    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    $regions \
    --examples make_examples${chunk}.tfrecord@${task.cpus}.gz \
    --gvcf gvcf${chunk}.tfrecord@${task.cpus}.gz
    """

    stub:
    """
    mkdir tmp2
    touch make_examples1.tfrecord-00000-of-00002.gz
    touch make_examples1.tfrecord-00001-of-00002.gz

    touch make_examples1.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples1.tfrecord-00001-of-00002.gz.example_info.json

    touch gvcf1.tfrecord-00000-of-00002.gz
    touch gvcf1.tfrecord-00001-of-00002.gz
    """
}

process MAKE_EXAMPLES_TRIO {
    tag "$meta.proband_id"
    container "docker://google/deepvariant:deeptrio-1.8.0-gpu"
    label "make_examples"

    input:
    tuple val(meta), path(bams), path(bais) // meta has proband_sex, proband_id, father_id, mother_id
    path(fasta_bams)
    path(fai_bams)

    output:
    tuple val(meta2), path("make_examples*_child.tfrecord*.gz"), path("gvcf*_child.tfrecord*.gz"), emit: proband_tfrecord
    tuple val(meta3), path("make_examples*.tfrecord*.example_info.json"), emit: example_info
    tuple val(meta4), path("make_examples*_parent1.tfrecord*.gz"), path("gvcf*_parent1.tfrecord*.gz"), emit: father_tfrecord, optional: true
    tuple val(meta5), path("make_examples*_parent2.tfrecord*.gz"), path("gvcf*_parent2.tfrecord*.gz"), emit: mother_tfrecord, optional: true

    script:
    def args = task.ext.args ?: ''

    is_male = meta.proband_sex == "Male" || meta.proband_sex == "male" || meta.proband_sex == "M"
    if(!is_male){
        assert meta.proband_sex == "Female" || meta.proband_sex == "female" || meta.proband_sex == "F"
    }
    if(params.test_bams){
        assert is_male && meta.proband_id.startsWith("HG00") && bams[0].size() < 50000000
    }
    def proband_id = meta.proband_id
    def father_id = meta.father_id
    def mother_id = meta.mother_id
    meta2 = [id: proband_id, proband_id: proband_id, role: "child", sex: meta.proband_sex]
    meta3 = [proband_id: proband_id]
    meta4 = [id: father_id, proband_id: proband_id, role: "parent", sex: "Male"]
    meta5 = [id: mother_id, proband_id: proband_id, role: "parent", sex: "Female"]
    if(params.annovar_buildver == "hg19"){
        par1endX = 2734539
        startX = 2734540
        endX = 154997473
        par2start = 154997474
        par2end = 155270560
        startY = 2649521
        endY = 59034049
    }
    if(params.annovar_buildver == "hg38"){
        par1endX = 2781479
        startX = 2781480
        endX = 155701382
        par2start = 155701383
        par2end = 156030895
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
    autosomes = "${pr}1 ${pr}2 ${pr}3 ${pr}4 ${pr}5 ${pr}6 ${pr}7 ${pr}8 ${pr}9 ${pr}10 ${pr}11 ${pr}12 ${pr}13 ${pr}14 ${pr}15 ${pr}16 ${pr}17 ${pr}18 ${pr}19 ${pr}20 ${pr}21 ${pr}22"

    // Set up code that will be the same for every run of make_examples
    mecmd = ["time seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples",
             "--mode calling --ref ${fasta_bams}",
             "--channel_list=read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size",
             "--pileup_image_height_child 100 --pileup_image_height_parent 100 --task {}",
             "--reads=${bams[0]} --sample_name ${proband_id} $args"].join(' ').trim()
    
    // Tiny example region of the genome for demo
    if(params.test_bams)
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --reads_parent2=${bams[2]} \
    --sample_name_parent1 ${father_id} \
    --sample_name_parent2 ${mother_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${pr}2:179300000-179600000 ${pr}16:450000-470000"

    $mecmd \
    --reads_parent1=${bams[2]} \
    --sample_name_parent1 ${mother_id} \
    --examples make_examples2.tfrecord@${task.cpus}.gz \
    --gvcf gvcf2.tfrecord@${task.cpus}.gz \
    --regions ${pr}X:122530000-122550000

    # Rename mother output to parent2
    for x in make_examples2_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done
    for x in gvcf2_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done
    """

    // Full genome, with sex chromosome handled based on proband sex and trio/duo status
    else if(is_male && father_id != "" && mother_id != "") // male, both parents
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --reads_parent2=${bams[2]} \
    --sample_name_parent1 ${father_id} \
    --sample_name_parent2 ${mother_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${autosomes} ${pr}X:1-${par1endX}"

    $mecmd \
    --reads_parent1=${bams[2]} \
    --sample_name_parent1 ${mother_id} \
    --examples make_examples2.tfrecord@${task.cpus}.gz \
    --gvcf gvcf2.tfrecord@${task.cpus}.gz \
    --regions ${pr}X:${startX}-${endX}

    # Rename mother output to parent2
    for x in make_examples2_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done
    for x in gvcf2_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done

    $mecmd \
    --reads_parent1=${bams[1]} \
    --reads_parent2=${bams[2]} \
    --sample_name_parent1 ${father_id} \
    --sample_name_parent2 ${mother_id} \
    --examples make_examples3.tfrecord@${task.cpus}.gz \
    --gvcf gvcf3.tfrecord@${task.cpus}.gz \
    --regions ${pr}X:${par2start}-${par2end}

    $mecmd \
    --reads_parent1=${bams[1]} \
    --sample_name_parent1 ${father_id} \
    --examples make_examples4.tfrecord@${task.cpus}.gz \
    --gvcf gvcf4.tfrecord@${task.cpus}.gz \
    --regions ${pr}Y:${startY}-${endY}

    rm gvcf4_parent2.tfrecord*.gz
    rm make_examples4_parent2.tfrecord*.gz
    """

    else if(!is_male && father_id != "" && mother_id != "") // female, both parents
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --reads_parent2=${bams[2]} \
    --sample_name_parent1 ${father_id} \
    --sample_name_parent2 ${mother_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${autosomes} ${pr}X"
    """

    else if(father_id == "" && mother_id != "") // male or female, mother only
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --sample_name_parent1 ${mother_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${autosomes} ${pr}X"

    # Rename mother output to parent2
    for x in make_examples1_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done
    for x in gvcf1_parent1.tfrecord*.gz; do
        mv -- "\$x" "\${x//parent1/parent2}"
    done
    """

    else if(is_male && father_id != "" && mother_id == "") // male, father only
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --sample_name_parent1 ${father_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${autosomes} ${pr}X:1-${par1endX}"

    $mecmd \
    --reads_parent1=${bams[1]} \
    --sample_name_parent1 ${father_id} \
    --examples make_examples3.tfrecord@${task.cpus}.gz \
    --gvcf gvcf3.tfrecord@${task.cpus}.gz \
    --regions "${pr}X:${par2start}-${par2end} ${pr}Y"

    # males with only father provided may need singleton variant calling on X chromosome.
    """

    //if(!is_male && father_id != "" && mother_id == "")
    else // female, father only
    """
    mkdir tmp2
    export TMPDIR=tmp2

    $mecmd \
    --reads_parent1=${bams[1]} \
    --sample_name_parent1 ${father_id} \
    --examples make_examples1.tfrecord@${task.cpus}.gz \
    --gvcf gvcf1.tfrecord@${task.cpus}.gz \
    --regions "${autosomes} ${pr}X"
    """

    stub:
    if(father_id != "" && mother_id != "")
    """
    mkdir tmp2
    touch make_examples1_child.tfrecord-00000-of-00002.gz
    touch make_examples1_child.tfrecord-00001-of-00002.gz
    touch make_examples2_child.tfrecord-00000-of-00002.gz
    touch make_examples2_child.tfrecord-00001-of-00002.gz
    touch make_examples1_parent1.tfrecord-00000-of-00002.gz
    touch make_examples1_parent1.tfrecord-00001-of-00002.gz
    touch make_examples2_parent1.tfrecord-00000-of-00002.gz
    touch make_examples2_parent1.tfrecord-00001-of-00002.gz
    touch make_examples1_parent2.tfrecord-00000-of-00002.gz
    touch make_examples1_parent2.tfrecord-00001-of-00002.gz
    touch make_examples2_parent2.tfrecord-00000-of-00002.gz
    touch make_examples2_parent2.tfrecord-00001-of-00002.gz

    touch make_examples1.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples1.tfrecord-00001-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00001-of-00002.gz.example_info.json

    touch gvcf1_child.tfrecord-00000-of-00002.gz
    touch gvcf1_child.tfrecord-00001-of-00002.gz
    touch gvcf2_child.tfrecord-00000-of-00002.gz
    touch gvcf2_child.tfrecord-00001-of-00002.gz
    touch gvcf1_parent1.tfrecord-00000-of-00002.gz
    touch gvcf1_parent1.tfrecord-00001-of-00002.gz
    touch gvcf2_parent1.tfrecord-00000-of-00002.gz
    touch gvcf2_parent1.tfrecord-00001-of-00002.gz
    touch gvcf1_parent2.tfrecord-00000-of-00002.gz
    touch gvcf1_parent2.tfrecord-00001-of-00002.gz
    touch gvcf2_parent2.tfrecord-00000-of-00002.gz
    touch gvcf2_parent2.tfrecord-00001-of-00002.gz
    """

    else if(father_id == "" && mother_id != "")
    """
    mkdir tmp2
    touch make_examples1_child.tfrecord-00000-of-00002.gz
    touch make_examples1_child.tfrecord-00001-of-00002.gz
    touch make_examples2_child.tfrecord-00000-of-00002.gz
    touch make_examples2_child.tfrecord-00001-of-00002.gz
    touch make_examples1_parent2.tfrecord-00000-of-00002.gz
    touch make_examples1_parent2.tfrecord-00001-of-00002.gz
    touch make_examples2_parent2.tfrecord-00000-of-00002.gz
    touch make_examples2_parent2.tfrecord-00001-of-00002.gz

    touch make_examples1.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples1.tfrecord-00001-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00001-of-00002.gz.example_info.json

    touch gvcf1_child.tfrecord-00000-of-00002.gz
    touch gvcf1_child.tfrecord-00001-of-00002.gz
    touch gvcf2_child.tfrecord-00000-of-00002.gz
    touch gvcf2_child.tfrecord-00001-of-00002.gz
    touch gvcf1_parent2.tfrecord-00000-of-00002.gz
    touch gvcf1_parent2.tfrecord-00001-of-00002.gz
    touch gvcf2_parent2.tfrecord-00000-of-00002.gz
    touch gvcf2_parent2.tfrecord-00001-of-00002.gz
    """

    //if(father_id != "" && mother_id == "")
    else
    """
    mkdir tmp2
    touch make_examples1_child.tfrecord-00000-of-00002.gz
    touch make_examples1_child.tfrecord-00001-of-00002.gz
    touch make_examples2_child.tfrecord-00000-of-00002.gz
    touch make_examples2_child.tfrecord-00001-of-00002.gz
    touch make_examples1_parent1.tfrecord-00000-of-00002.gz
    touch make_examples1_parent1.tfrecord-00001-of-00002.gz
    touch make_examples2_parent1.tfrecord-00000-of-00002.gz
    touch make_examples2_parent1.tfrecord-00001-of-00002.gz

    touch make_examples1.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples1.tfrecord-00001-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00000-of-00002.gz.example_info.json
    touch make_examples2.tfrecord-00001-of-00002.gz.example_info.json

    touch gvcf1_child.tfrecord-00000-of-00002.gz
    touch gvcf1_child.tfrecord-00001-of-00002.gz
    touch gvcf2_child.tfrecord-00000-of-00002.gz
    touch gvcf2_child.tfrecord-00001-of-00002.gz
    touch gvcf1_parent1.tfrecord-00000-of-00002.gz
    touch gvcf1_parent1.tfrecord-00001-of-00002.gz
    touch gvcf2_parent1.tfrecord-00000-of-00002.gz
    touch gvcf2_parent1.tfrecord-00001-of-00002.gz
    """
}

process GLNEXUS {
    //container = 'https://depot.galaxyproject.org/singularity/glnexus:1.4.1--h671cb6e_1' // no jemalloc; slow
    container 'docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1'
    publishDir "$params.outdir", mode: 'copy', overwrite: true

    input:
        path(gvcfs)
        val(cohort_name)
        val(glnexus_filter)
        val(deepvar_model)
    output:
        path(outbcf)
    script:
        outbcf = "${cohort_name}.glnexus.bcf"
        if(glnexus_filter){
            mod = deepvar_model
        }
        else {
            mod = '_unfiltered'
        }
    """
    glnexus_cli \
        -t ${task.cpus} \
        --config DeepVariant${mod} \
        *.gvcf.gz > ${outbcf}
    """

    stub:
    """
    touch "${cohort_name}.glnexus.bcf"
    """
}

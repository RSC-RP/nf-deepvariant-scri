process CALL_VARIANTS_SINGLE {
    container "docker://google/deepvariant:1.8.0-gpu"
    tag "$meta.id"
    label "call_variants"

    // me_tfrecord is output from make_examples.  role is "child" or "parent".
    input:
    tuple val(meta), path(me_tfrecord), path(gvcf_tfrecord), path(example_info) // meta has id, proband_id, and role

    output:
    tuple val(meta), path("call_variants*.tfrecord.gz"), path(gvcf_tfrecord)

    script:

    sample_id = meta.id

    // Build argument for --examples and --outfile (list of call sets to loop through)
    tfr_first = me_tfrecord.findAll{ it.name ==~ /.*tfrecord-00000-of.*/ }
    tfr_prefix = tfr_first.collect{ it.simpleName - "make_examples" }.sort()
    tfr_string = tfr_prefix.join(' ')

    """
    mkdir tmp
    export TMPDIR=tmp
    export TF_GPU_ALLOCATOR=cuda_malloc_async
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib/python3.10/dist-packages/tensorrt_libs

    for value in $tfr_string
    do
        call_variants \
        --batch_size 384 \
        --outfile "call_variants\${value}_${sample_id}.tfrecord.gz" \
        --examples "make_examples\${value}.tfrecord@${params.make_examples_nshards}.gz" \
        --checkpoint /opt/models/${params.deepvar_model.toLowerCase()}
    done
    """

    stub:
    """
    mkdir tmp
    touch call_variants1_${sample_id}.tfrecord.gz
    touch call_variants2_${sample_id}.tfrecord.gz
    """
}

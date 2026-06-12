process POSTPROCESS_VARIANTS {
    container "docker://google/deepvariant:deeptrio-1.8.0-gpu"
    publishDir "${params.outdir}/gvcfs/", mode: 'copy', overwrite: true
    tag "$meta.id"

    input:
    tuple val(meta), path(cv_tfrecord), path(gvcf_tfrecord)
    path(fasta_ensembl)
    path(fai_ensembl)
    path(par_bed)

    output:
    tuple val(meta), path(out_gvcf), path("${out_gvcf}.tbi")

    script:
    sample_id = meta.id
    out_gvcf = "${sample_id}.gvcf.gz"
    out_vcf = "${sample_id}.vcf.gz"

    // Build argument for --nonvariant_site_tfrecord_path and --infile
    tfr_prefix = gvcf_tfrecord.collect{ it.simpleName.replaceFirst(~/^gvcf/, "") }.unique( false ).sort() // 1_parent1
    tfr_string = tfr_prefix.join(' ')
    cv_shards = cv_tfrecord[0].simpleName.replaceFirst(~/.*\-of\-/, "").toInteger()
    gvcf_shards = gvcf_tfrecord[0].name.replaceFirst(".gz", "").replaceFirst(~/.*\-of\-/, "").toInteger()

    // --cpus needs to be set to 0 if not using parallel processing
    pvcpu = task.cpus == 1 ? 0 : task.cpus

    // haploid contigs for males
    is_male = meta.sex == "Male" || meta.sex == "male" || meta.sex == "M"
    if(is_male){
        pr = params.chromnames == "ucsc" ? "chr" : ""
        hapstring = "--haploid_contigs=\"${pr}X,${pr}Y\" --par_regions_bed=\"${par_bed}\""
    }
    else{
        assert meta.sex == "Female" || meta.sex == "female" || meta.sex == "F"
        hapstring = ""
    }

    """
    mkdir tmp
    export TMPDIR=tmp

    for value in $tfr_string
    do
        postprocess_variants \
        --ref $fasta_ensembl \
        --infile "call_variants\${value}_${sample_id}@${cv_shards}.tfrecord.gz" \
        --outfile temp.vcf.gz \
        --nonvariant_site_tfrecord_path "gvcf\${value}.tfrecord@${gvcf_shards}.gz" \
        --gvcf_outfile temp.gvcf.gz \
        $hapstring \
        --cpus $pvcpu

        if test -f "$out_gvcf"; then
            mv $out_gvcf temp1.gvcf.gz
            bcftools concat -O z -o $out_gvcf --threads ${task.cpus} temp1.gvcf.gz temp.gvcf.gz
        else
            mv temp.gvcf.gz $out_gvcf
        fi
    done

    bcftools index -t --threads ${task.cpus} $out_gvcf
    """

    stub:
    """
    mkdir tmp
    touch "${sample_id}.gvcf.gz"
    touch "${sample_id}.gvcf.gz.tbi"
    touch temp.gvcf.gz
    touch temp.gvcf.gz.tbi
    touch temp1.gvcf.gz
    touch temp.vcf.gz
    touch temp.vcf.gz.tbi
    touch temp.visual_report.html
    """
}

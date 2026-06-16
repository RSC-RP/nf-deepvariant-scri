process MAKE_PEDIGREE {
    tag "$prefix"
    container 'https://depot.galaxyproject.org/singularity/python:3.14'

    input:
    path(bam_csv)
    path(sample_lookup)
    val(prefix)

    output:
    path(ped)

    script:
    ped = "${prefix}.somalier.ped"
    """
    make_pedigree.py $bam_csv $sample_lookup $ped
    """
}

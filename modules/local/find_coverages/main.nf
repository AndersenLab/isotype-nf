process FIND_COVERAGES {
    container null
    executor "local"

    input:
    path bams

    output:
    path "coverages.txt", emit: coverages
    path "versions.yml",  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    for I in ${bams}/*.coverage; do
        tail -n +2 \${I} | cut -f 1,2 >> coverages.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """

    stub:
    """
    touch coverages.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """
}

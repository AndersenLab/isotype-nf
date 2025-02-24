process GTCHECK {
    label 'gtcheck'
    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 24 * task.attempt }
    memory = { 80.GB * task.attempt }

    input:
    path "vcf.npy"
    path "gtcheck.py"

    output:
    path "gtcheck.txt" , emit: gtcheck
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python gtcheck.py vcf.npy gtcheck.txt ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
    END_VERSIONS
    """

    stub:
    """
    touch gtcheck.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
    END_VERSIONS
    """
}

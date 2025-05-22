process CALL_ISOTYPES {
    label 'call_isotypes'
    errorStrategy 'retry'
    time { 1.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    path "gtcheck.txt"
    path "call_isotypes.py"
    path "compare_isotype_calls.py"
    path "previous_isotypes.tsv"
    path coverages
    val concordance_cutoff

    output:
    path "isotype_groups.tsv",                                            emit: groups
    tuple path("isotype_comparison.pdf"), path("isotype_comparison.txt"), emit: comparison
    path "isotype_summary.txt",                                           emit: summary
    path "versions.yml",                                                  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR="."
    python call_isotypes.py gtcheck.txt previous_isotypes.tsv ${coverages} ${concordance_cutoff}
    python compare_isotype_calls.py gtcheck.txt previous_isotypes.tsv isotype_comparison.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
        scipy: \$( python -c "import scipy; print(scipy.version.version)" )
    END_VERSIONS
    """

    stub:
    """
    touch isotype_groups.txt
    touch isotype_comparison.pdf
    touch isotype_comparison.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
        scipy: \$( python -c "import scipy; print(scipy.version.version)" )
    END_VERSIONS
    """
}

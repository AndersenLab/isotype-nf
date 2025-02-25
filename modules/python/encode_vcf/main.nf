process ENCODE_VCF {
    label 'encode_vcf'
    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 1500000.MB * task.attempt }

    input:
    path vcf
    path "encode_vcf.py"

    output:
    path "vcf.npy"      , emit: encoded
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python encode_vcf.py ${vcf} vcf.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
    END_VERSIONS
    """

    stub:
    """
    touch vcf.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
        numpy: \$( python -c "import numpy; print(numpy.version.version)" )
    END_VERSIONS
    """
}

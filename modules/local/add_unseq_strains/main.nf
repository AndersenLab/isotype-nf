process ADD_UNSEQ_STRAINS {
    container null
    executor "local"

    input:
    path groups
    path unseq_strains

    output:
    path "isotype_groups_all.tsv", emit: groups
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo -e "group\tstrain\tisotype\tisotype_ref_strain" > isotype_groups_all.tsv
    awk 'BEGIN{OFS="\t"}{
          if ( NR == FNR ){
            if ( \$1 !~ /^group/ ){
              if ( \$3 in ISOTYPES ){
                ISOTYPES[\$3]=ISOTYPES[\$3]","\$2;
              } else {
                ISOTYPES[\$3]=\$2;
                REFSTRAINS[\$3]=\$4;
                GROUPNUM[\$3]=\$1;
              }
            }
          }
          else {
            if ( \$1 !~ /^strain/ ){
              if ( \$2 in ISOTYPES ){
                ISOTYPES[\$2]=ISOTYPES[\$2]","\$1;
              } else {
                ISOTYPES[\$2]=\$1;
                REFSTRAINS[\$2]=\$2;
                GROUPNUM[\$2]=length(ISOTYPES)+1;
              }
            }
          }
        }END{
          for (ISOTYPE in ISOTYPES){
            NUM=GROUPNUM[ISOTYPE];
            REFSTRAIN=REFSTRAINS[ISOTYPE];
            split(ISOTYPES[ISOTYPE],STRAINS,",");
            for (I=1;I<=length(STRAINS);I++){
              STRAIN=STRAINS[I];
              printf "%i\t%s\t%s\t%s\n", NUM, STRAIN, ISOTYPE, REFSTRAIN;
            }
          }
        }' ${groups} ${unseq_strains} | sort -k1,1n >> isotype_groups_all.tsv
    """

    stub:
    """
    touch isotype_groups_all.tsv
    """
}

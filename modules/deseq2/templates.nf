#!/usr/bin/env nextflow

/*
 * Example design and contrasts files for RNA-Seq pipeline
 * Creates template files for differential expression analysis
 */

process CREATE_DESIGN_TEMPLATE {
    label 'process_low'
    
    publishDir "${params.outdir}/templates", mode: 'copy'
    
    output:
    path "design_template.csv", emit: design
    
    script:
    """
    cat > design_template.csv << EOL
sample_id,condition,batch
sample1,control,1
sample2,control,1
sample3,control,2
sample4,control,2
sample5,treatment,1
sample6,treatment,1
sample7,treatment,2
sample8,treatment,2
EOL
    """
}

process CREATE_CONTRASTS_TEMPLATE {
    label 'process_low'
    
    publishDir "${params.outdir}/templates", mode: 'copy'
    
    output:
    path "contrasts_template.csv", emit: contrasts
    
    script:
    """
    cat > contrasts_template.csv << EOL
name,control,treatment
treatment_vs_control,control,treatment
EOL
    """
}

workflow CREATE_DE_TEMPLATES {
    main:
    CREATE_DESIGN_TEMPLATE()
    CREATE_CONTRASTS_TEMPLATE()
    
    emit:
    design = CREATE_DESIGN_TEMPLATE.out.design
    contrasts = CREATE_CONTRASTS_TEMPLATE.out.contrasts
}

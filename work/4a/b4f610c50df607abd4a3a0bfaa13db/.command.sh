#!/bin/bash -ue
# Run MultiQC with custom config
multiqc -f \
    -c multiqc_config.yaml \
    -n multiqc_report.html \
    .

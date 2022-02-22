#!/bin/bash
awk '{if($3=="transcript") {OFS="\t"; print $14, $10} }' genes.gtf | sed 's/[;\"]//g' > t2g.tsv

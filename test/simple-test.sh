#!/bin/bash
set -euo pipefail

m4=$1
d=$2
mkdir -p $d


ref_vcf=$d/$(basename $3).gz
ref_msav=$d/$(basename $3 .vcf).msav
tar_vcf=$d/$(basename $4).gz
imputed_vcf=$d/imputed.vcf.gz

which bcftools || (>&2 echo "Error: bcftools is required to run tests. On debian run 'apt install bcftools'"; exit 1)

bcftools view $3 -Oz -o $ref_vcf
bcftools index $ref_vcf
bcftools view $4 -Oz -o $tar_vcf
bcftools index $tar_vcf
$m4 --compress-reference $ref_vcf > $ref_msav
$m4 $ref_msav $tar_vcf -f GT -O vcf.gz --temp-buffer 2 > $imputed_vcf

gzip -cd $imputed_vcf | grep -v "^#" | cut -f9- > $d/imputed_gt_matrix.tsv
gzip -cd $ref_vcf | grep -v "^#" | cut -f9- > $d/ref_gt_matrix.tsv
diff $d/ref_gt_matrix.tsv $d/imputed_gt_matrix.tsv
#!/bin/bash
set -euo pipefail

d=m4_simple_test_output
mkdir -p $d

m4=$1
ref_vcf=$d/$(basename $2).gz
ref_msav=$d/$(basename $2 .vcf).msav
tar_vcf=$d/$(basename $3).gz
imputed_vcf=$d/imputed.vcf.gz

which bcftools || (>&2 echo "Error: bcftools is required to run tests. On debian run 'apt install bcftools'"; exit 1)

bcftools view $2 -Oz -o $ref_vcf
bcftools index $ref_vcf
bcftools view $3 -Oz -o $tar_vcf
bcftools index $tar_vcf
$m4 --compress-reference $ref_vcf > $ref_msav
$m4 $ref_msav $tar_vcf -f GT -O vcf.gz > $imputed_vcf

gzip -cd $imputed_vcf | grep -v "^#" | cut -f9- > $d/imputed_gt_matrix.tsv
gzip -cd $ref_vcf | grep -v "^#" | cut -f9- > $d/ref_gt_matrix.tsv
diff $d/ref_gt_matrix.tsv $d/imputed_gt_matrix.tsv
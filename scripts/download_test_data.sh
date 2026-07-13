#!/bin/bash
# Download chr22-subset reference data for testing vartriage.
# Total download: ~3-4 GB. Run from the project root.

set -e
cd "$(dirname "$0")/.."
mkdir -p data/references

echo "=== 1/4: GENCODE GTF (gene annotation) ==="
echo "    ~50 MB compressed"
wget -q --show-progress -O data/references/gencode.v46.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
gunzip -f data/references/gencode.v46.annotation.gtf.gz
echo "    Done: data/references/gencode.v46.annotation.gtf"

echo ""
echo "=== 2/4: ClinVar VCF (clinical significance) ==="
echo "    ~30 MB"
wget -q --show-progress -O data/references/clinvar.vcf.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -q --show-progress -O data/references/clinvar.vcf.gz.tbi \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

echo "    Converting to TSV format..."
# ClinVar GRCh38 VCF uses bare chromosome names (1, 2, ..., 22, X, Y)
# but vartriage VCFs typically use 'chr' prefix. Normalize to chr-prefixed.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\n' \
  data/references/clinvar.vcf.gz \
  | awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' \
  | sed '1i\chrom\tpos\tref\talt\tclinical_significance' \
  > data/references/clinvar.tsv
echo "    Done: data/references/clinvar.tsv (chr-prefixed)"

echo ""
echo "=== 3/4: gnomAD chr22 exomes (population frequency) ==="
echo "    ~4.7 GB - this will take a while"
wget -q --show-progress -O data/references/gnomad.chr22.vcf.bgz \
  https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1.1/vcf/exomes/gnomad.exomes.v4.1.1.sites.chr22.vcf.bgz
wget -q --show-progress -O data/references/gnomad.chr22.vcf.bgz.tbi \
  https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1.1/vcf/exomes/gnomad.exomes.v4.1.1.sites.chr22.vcf.bgz.tbi

echo "    Converting to TSV format (chrom, pos, ref, alt, af)..."
# gnomAD GRCh38 VCF uses chr-prefixed names already
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
  data/references/gnomad.chr22.vcf.bgz \
  | sed '1i\chrom\tpos\tref\talt\taf' \
  > data/references/gnomad_chr22.tsv
echo "    Done: data/references/gnomad_chr22.tsv"

echo ""
echo "=== 4/4: CADD scores for chr22 ==="
echo "    NOTE: CADD whole-genome is 350GB+. For testing, we create a minimal file."
echo "    For real analysis, download from: https://cadd.gs.washington.edu/download"
echo "    The sample file below covers a few known positions. After running the"
echo "    pipeline once with your VCF, use scripts/prepare_references.sh to generate"
echo "    a CADD file that overlaps with your actual variants."

# Create a minimal CADD test file (header uses # prefix which vartriage skips)
cat > data/references/cadd_chr22_sample.tsv << 'EOF'
#chrom	pos	ref	alt	score
chr22	16364843	G	A	23.4
chr22	16389542	C	T	28.1
chr22	17565849	G	A	34.0
chr22	19710952	C	T	15.7
chr22	20916938	A	G	22.3
chr22	24130192	G	A	32.5
chr22	29065370	G	A	26.8
chr22	36661906	A	G	19.2
chr22	42128277	C	T	35.0
chr22	50356620	G	A	24.9
EOF
echo "    Done: data/references/cadd_chr22_sample.tsv (sample only)"

echo ""
echo "=== All done ==="
echo "Files in data/references/:"
ls -lh data/references/

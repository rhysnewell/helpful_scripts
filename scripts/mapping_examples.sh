
#With Freebayes and naive variant calling
freebayes -v 08-freebayes/r2.parent.d77.vcf -f r2.parent.d77.assembly.fna --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic -c < <(minimap2 -ax sr r2.parent.d77.assembly.fna r2.parent.d77/r2.parent.d77.1.fq.gz.trimmed.1.fq.gz r2.parent.d77/r2.parent.d77.1.fq.gz.trimmed.2.fq.gz | samtools view -Sub -F4 | samtools sort)

freebayes -v 08-freebayes/r2.parent.d182.vcf -f r2.parent.d182.assembly.fna --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic -c < <(minimap2 -ax sr r2.parent.d182.assembly.fna r2.parent.d182/r2.parent.d182.1.fq.gz.trimmed.1.fq.gz r2.parent.d182/r2.parent.d182.1.fq.gz.trimmed.2.fq.gz | samtools view -Sub -F4 | samtools sort)

#BWA
bwa mem -t 20 r2.parent.d77.assembly.fna r2.parent.d77/r2.parent.d77.1.fq.gz.trimmed.1.fq.gz r2.parent.d77/r2.parent.d77.1.fq.gz.trimmed.2.fq.gz | samtools view -Sub -F4 | samtools sort

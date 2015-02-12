# geno_huffman

Huffman tree based compression algorithm for SNP data. <br>
<b>Input</b>: .bam file <br>
<b>Usage</b>:  -in (input .bam file name) -snp (number of snps) -ind (number of individuals)  <br>
<b>Options</b>: <br>-maf_min  -maf_max:  minimum and maximum threshold of Minor Allele Frequency for SNPs, 
that will be included into analysis<br>
-cor Usage of correlation between adjacent snps (0 - do not use correlation, 1 - use correlation)<br>
<b>Output</b>: <br>
"Huffman_result.txt" file will include information of compressed snps: number of correlated snps, frequencies of genotypes, used
in Huffman algorithm, optimal number of genotypes in block, final length of encoded string. <br>
"Huffman_ind_result.txt" file will include information about individuals: number of individual, bite cost of individual <br>





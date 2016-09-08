# filter.gz file

Quote from an email discussion with David Reich (2015/12/17) about which set of sites from the archaic admixture array to use for the "Ice Age" paper (the positions of these sites are in "filter.gz" file):

#### 954,456 “All"

This is the total number of SNPs in the file I sent, originally prepared by Qiaomei. This corresponds to the total number of SNPs in the archaic admixture capture reagent (1,749,385 SNPs) that have the characteristic of having at least one allele in a Neanderthal (Altai, Mezmaiskaya, or Vindija) that is different from the allele in the great majority of Yoruba from the 1000 Genomes Project. Qiaomei says that this list is chosen *BLINDED* to the allele type in Yoruba.

#### 783,974 “Set 1"
This is my (new) filtering of QIaomei’s file. I took her random alleles, and then looked for fixed difference between Altai on the one hand, and the B Panel Yoruba on the other hand (I have now dropped the Mbuti individual; this adds ~10,000 SNPs). I now realize that there was a difference relative to your filtering, though, which is that I should have used the diploid calls. I made a mistake here, as my description of the set I used was that it was based on diploid calls and this was incorrect.

#### 664,625 “Set 2"
This is a new filter of Qiaomei’s diploid call file that I just made (I didn’t send you the diploid file before). I took her random alleles, and then looked for a fixed difference between Altai on the one hand, and the B-panel Yoruba (on individual) on the other hand. That is, they are homozygous different from each other. This is still ~190K more SNPs than you have.

I think we *definitely* need to filter down from the 954,456 SNPs for our analysis, since that ascertainment uses Mezmaiskaya and/or Vindija, which we use as setting a baseline for the match rate to Altai expected for a true Neanderthal.

My preference is for “Set 1”, since:
1. It gives substantially more SNPs and thus power
2. Random allele calling is cleaner from a population genetic point of view (it means that the ascertainment is not be dependent on the demographic history of Neanderthals including the Altai Neanderthal since separation from the modern human lineage).

If you agree, I would suggest redoing your analyses (all analyses including individual ancestry estimates and regressions) with “Set 1”. I attach a gzipped “filter.gz” file with the 954,456 SNPs, in the same order as the random allele call file I sent you before. Columns are:
- Chromosome
- hg19 position	
- Set 1 filter (1/0)
- Set 2 filter (1/0)


# "individuals" file

This is a list of the 48 individuals in the same order that they appear in the genotype string in the “genosnp” file. Columns are:
1. Index (order in genotype string)
2. Sample ID
3. Name of sample in paper used for example in Table S3.1
4. Date of sample that is the mean of the range of Extended Data Table 1 and should be used in the regression
5. Flag indicated if the sample is listed in Table S3.1 (1=listed, 0=not listed). Please don’t use for analysis samples that don’t have a 1 in this column. For example, this flag indicates which version of Kostenki14 and Oase1 should be used (Oase1_d is Mateja’s processing).
6. Flag indicating if the sample should be used in the main least squares fit. A total of 29 samples are indicated, including 26 ancient modern humans along with Han, Dai and Karitiana. I suggest that you also do a second regression dropping Han, Dai, Karitiana, and Stuttgart (and possibly a third also dropping Ust’-Ishim) to ensure that the results are consistent.


# "genosnp_*.gz" files

This contains 956,456 rows corresponding to the SNPs from the archaic admixture array that are informative about Neanderthal ancestry (the SNPs you picked are a subset of these). The columns are:
1. Chromosome
2. hg19 position
3. Chromosome_hg19position
4. A silly conversion to genetic position that you should ignore
5. SNP allele 0
6. SNP allele 1
7. B-statistic decile (Shop’s version)
8. Genotype string with randomly selected alleles for the 48 individuals. Codes are: (0=allele 0), (1=allele 1), (9=missing data)
9. Absolute value of B-statistic x 1000
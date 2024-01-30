
# piphillinR

The package performs the functional prediction from 16S rRNA sequencing data by the algorithm in [Piphillin software](https://doi.org/10.1186/s12864-019-6427-1). The other paper is [here](https://doi.org/10.1371/journal.pone.0166104). This implementation referenced the Python implementation, pyphillin.

## Usage

We first need to compile data we need to use Piphillin algorithm. The pre-computed files are available at [TBU URL](URL).
After downloading, run `alignSequences` function to search the representative sequences in the reference 16S sequences. The function needs `vsearch` executable in `PATH` to work. Installation instruction of `vsearch` can be found at the [official repository](https://github.com/torognes/vsearch). After the global alignment is finished, `profileMetagenome` function can be used with the following input.

- ASV abundance table (typically obtained from DADA2)
- BLAST results from vsearch (the second column should match the copy number table)
- Gene copy number table (typically KEGG ORTHOLOGY)
- 16S rRNA copy number table

```r
res <- profileMetagenome(taxTable, copyNumTable, KOTable, blastRes)
```

## Constructing the customized database

If you have KEGG license, it would be relatively easy to construct the database (See the section "Reference databases" in the [original paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6427-1#Sec9)). The 16S sequences are extracted by the KO assignment of `K01977`.

## Comparison with PICRUSt2 (correlation analysis)

Although the database is different between PICRUSt2 (IMG) and Piphillin (KEGG), we conducted the correlation analysis between the KO abundance table produced by two software.

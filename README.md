# Generates genomes around ani values based on gene RBMs along ANI gradient

    - uses a gamma distribution to model RBM gene distribution around ANI
    - Genome evolution is neutral random based on RBM distriubtion.
    - Average gene length and stdev (default 1000 bp, 250 stdev)
    - * no deletions, insertions or rearrangements * (future?)
    - Gene synteny is conserved
    - Genes start with ATG and stop with TAG
    - 10 bps between every gene (see future additions)
    - Genome starts and ends with 10 additional basepairs before/after genes
    - Pangenome Gene classes randomly distributed Core or Accessory.
    - Creates draft genomes choosing n contigs from gaussian mu=10,sigma=5.
    - Contigs number stdev = average contigs / 2
    - creates 15% * source genes additional auxiliary genes.
    - shuffles auxillary genome choosing shared fraction from gaussion
      where mu=y=0.02(ANI)-1.05 and sigma=0.02

No input files are needed.

Writes a genome.fna, a cds.fnn for each genome.
Writes a gene cluster file with pangenome classifications.
Writes histogram PDF of RBM distribution for each simulated genome to source genome.

Needs Python version 3.9+ for random.sample counts parameter.

The random.seed is set at 42 for reproducibility

```bash
# print the help file and parameter list
python Simulate_population_genomes.py -h

# run with default settings
python Simulate_population_genomes.py -o my_simulated_genomes
```

# future additions
    - operons - sizes are randomly sampled from gaussian mu=10genes, std=5
    - intergenic regions are randomly sampled from gaussian mu=5bp, std=2
    - interoperon regions are randomly sampled from gaussian mu=100bp, std=25
    - recombination rate default 10%? simulate recombination event in 1/10
      of the genomes. How to keep ANI the same distance with increased F100?
      if some genes are more similar some genes need to be less similar.
      * also what fraction of genes are recombined?
    - accept input genome and gene fasta as seed genome
    - add mutation rate per generation per number of generations feature
        - no ANI model but keep gene and pangenome function
    - use daughter genomes as new source genomes for future generations
    - user input source genome

Currently mutations are added at random across each gene to match a gamma distribution with a mean approximate around a specified ANI value. There are no recombination events, insertions, or deletions. All new genomes are altered versions of the initial source genome that is simulated first.
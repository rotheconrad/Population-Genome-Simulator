# Generates genomes around ani values based on gene RBMs along ANI gradient

Step 1: randomly generate a source genome with a user specified number of genes of length mu with a standard deviation. Default is 3000 genes of length mu=1000 and stdev=250. So gene lengths are randomly chosen from a normal distribution. Then add ATG and TAG to the beginning and end of each gene, and add 10bps (random sequence) between each gene plus 10bps to the beginning and end. Designate 85% of the genes as core genes and 15% as accessory. Randomly generate additional pool of accessory genes that can be randomly distributed across the population.

Step 2: use the input ANI range (default 95-100) with a step size (default 0.1) and generate gamma distributions with a mean matching the ANI target. For each new genome simulation, generate new random RBM gene distribution from the ANI matched gamma distribution then add random point mutations to the RBM gene in the source genome to create the new genome. Lastly, shuffle random fraction of accessory genes.

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
# SVCollector: Optimized sample selection for cost-efficient long-read population sequencing

Structural Variations (SVs) are increasingly recognized for their importance in genomics. Short-read sequencing is the most widely-used approach for genotyping large numbers of samples for SVs but suffers from relatively poor accuracy. Here we present SVCollector, an open-source method that optimally selects samples to maximize variant discovery and validation using long read resequencing or PCR-based validation. SVCollector has two modes: selecting those samples that are individually the most diverse or those that collectively capture the largest number of variations.

If you experience problems or have suggestions please post an issue here or contact: fritz.sedlazeck@gmail.com


# How to build SVCollector

<pre>$ wget https://github.com/fritzsedlazeck/SVCollector/archive/master.zip -O SVCollector.tar.gz
$ tar xzvf SVCollector.tar.gz
$ cd SVCollector-master/Debug
$ make

$ ./SVCollector
</pre>

# Running SVCollector

SVCollector can be run in a few different modes (greedy, topN, or random) with an input multi-sample VCF file and outputs a ranked list of the samples to select. The top level command is run like this:


```
$ ./Debug/SVCollector
./SVCollector <option> my_svs_vcf_file output_ranked
<option>: greedy, topN or random
my_svs_vcf_file: A valid uncompressed multisample VCF file.
num_samples: The number of samples that should be ranked.
output_ranked: The file to write out the ranked list with additional information.
```

### Greedy Analysis

This is the recommended mode for all users as it will best optimize the selection of samples. The parameters for this mode are shown below.

```
$ ./Debug/SVCollector greedy
Input VCF file
Min allele count (-1 to disable)
Number of samples to select
Take AF into account (1) or not (0) per allele
Optionally: File of names to select anyways (NA to disable)
Optionally: Text File of names and weights (NA to disable)
Output file
```

### topN Analysis

This mode is provided for comparison purposes to evaluate how the greedy mode compares to this simplier selection mode. The parameters for this mode are shown below.

```
$ ./Debug/SVCollector topN
Input VCF file
Number of samples to select
Take AF into account (1) or not (0) per allele
Output file
```

### random Analysis

This is the most naive approach that just picks N samples at random from the entire input VCF file. The parameters for this mode are shown below.

```
$ ./Debug/SVCollector random
Input VCF file
Number of samples to select
Take AF into account (1) or not (0) per allele
Output file
```

### Running all modes

We also provide a helper script (SVCollector.sh) that will run all 3 modes (greedy, topN, and 10 trials of a random selection) and make a simple plot comparing the results over the first `numtoplot` samples from the input VCF file. If you have [GNU Parallel](https://www.gnu.org/software/parallel/) installed you should edit the script to replace the for loop with the much faster parallel version.

```
$ ./SVCollector.sh
USAGE: SVCollector.sh samples.vcf numtoplot workdir
```


# Demo

For evaluation, we include a simulation script that generates a multi-sample VCF file with an arbitrary population structure. Briefly, the simulator simulates `F` founder genotypes, that each contain on average `Normal(N,M)` variants placed at random along the genome (the initial genome size is fixed at 100,000,000 bp). Then for each founder population, a collection of `Normal(S,T)` individual samples are generated at random that contain the original founder variants plus an additional `Normal(X,Y)` variants. Consequently, the expected total number of variants in the collection is `F * N + F * S * X` variants. If `N > S`, then most of the variants will be shared within the population group, and if `S > N`, most variants will be unique to that sample. We emphasize this is not designed to simulate realistic pedigrees, but to examine the extremes of high or low levels of sharing among the individuals.


After compiling the code, you can run the demo like this:

<pre>
$ cd SVCollector/simul
$ ./demo.sh
</pre>

Note you will need to have `perl` and the `Math/Random` CPAN package installed. This can be installed with conda as:

```
$ conda install perl-math-random
```

The demo script will generate 2 simulated populations (simple10.vcf and complex10.vcf) and 2 working directories with the results of a greedy selection, a topN selection, and 10 trials of a random selection from these populations.

The first simulated population (simple10.vcf) has 10 founder genomes with exactly 1000 variants located at random. From each founder genome, 10 samples are simulated that contain the 1000 founder variants plus an additional 100 variants. The sharp inflection point for the greedy curve at N=10 illustrates how the code realizes there are 10 founder populations. After these 10 populations have been sampled, the rate at which additional SVs are identified reduces to a much slower rate as these variants are only contained in individual samples. The plot will be generated to `simul/simple10/simple10.vcf.png`.

The second simulated population (complex10.vcf) also has 10 founder genomes that each randomly contain `Normal(500,250)` variants. From each founder population, `Normal(10, 5)` samples are simulated that contain the founder variants, plus an additional `Normal(500,250)` variants unique to this sample. Notice that despite having the same number of founder genomes, the curves are substantially different, and lack the inflection point at `N=10`. This highlights how sample specific variants contribute at a similar level to the founder genotypes. The plot will be generated to `simul/complex10/complex10.vcf.png`.


# Citation

Please cite our preprint:
https://www.biorxiv.org/content/10.1101/2020.08.06.240390v1
# SVCollector

Structural Variations (SVs) are increasingly recognized for their importance in genomics. Short-read sequencing is the most widely-used approach for genotyping large numbers of samples for SVs but suffers from relatively poor accuracy. Here we present SVCollector, an open-source method that optimally selects samples to maximize variant discovery and validation using long read resequenc- ing or PCR-based validation. SVCollector has two modes: selecting those samples that are individually the most diverse or those that collectively capture the largest number of variations.
If you experience problems or have suggestions please post an issue here or contact: fritz.sedlazeck@gmail.com


# How to build SVCollector

<pre>wget https://github.com/fritzsedlazeck/SVCollector/archive/master.zip -O SVCollector.tar.gz
tar xzvf SVCollector.tar.gz
cd SVCollector-master/Debug
make

./SVCollector
</pre>


# Citation

Please cite our preprint:
https://www.biorxiv.org/content/early/2018/06/08/342386

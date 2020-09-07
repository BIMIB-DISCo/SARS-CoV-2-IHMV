# VERSO Variant Calling

The two scripts in this folder perform *variant calling* for *single*
and *paired ends* using standard tools.

Both scripts run on a UN\*X/Linux platform and require the following tools,
which must be installed separatedly and be available in the `PATH`.

- [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic "trimmomatic") version 0.39.
- [`bwa`](http://bio-bwa.sourceforge.net/ "bwa") version 07.17.
- [`samtools`](http://www.htslib.org/ "Samtools") version 1.10.
- [`picard`](https://broadinstitute.github.io/picard/ "Picard") version 2.23.4.
- [`varscan`](http://dkoboldt.github.io/varscan/ "VarScan") version 2.4.2.

The two scripts implement a traditional variant calling pipeline
following, e.g., [**Best Practices for Variant Calling using
GATK**](https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1).


### License

Please check the LICENSE file in the top directory for licensing
information.

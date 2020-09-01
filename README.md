# TensorSignature Data

A docker container environment to convert VCF data to a mutation count tensor for the TensorSignatures software.

## Installation

To use the image pull it from docker. Note that the image is approximately 5 GB large.

```
$ docker pull sagar87/tensorsignatures-data:latest
```

## Usage

To use the image switch into the folder containing your VCF data. Then run image using the following command and supply the VCF files as arguments as well as the name of the output file (must be the last argument) as arguments.

```
docker run -v $PWD:/usr/src/app/mount sagar87/tensorsignatures-data <vcf1.vcf> <vcf2.vcf> ... <vcfn.vcf> <output.h5>
```

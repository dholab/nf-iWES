# NextFlow pipeline for genotyping immunological genes from macaque whole exome sequence data

My latest invocation for testing purposes, as an example:
```
nextflow run iwes.nf \
--create_reference_files false \
--fastq_dir "/Volumes/T7 1/27734/Baylor_33/rt24330" \
--run_name Baylor_33 \
--animal mamu
```

This workflow is still under development and is not yet ready for replication by other users. Still, if you'd like to set up NextFlow and Docker and run it for yourself, we recommend you consult the setup instructions from [one of our previous workflows.](https://github.com/dholab/MHC-II-allele-discovery#detailed-instructions)
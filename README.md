Analysis workflow for community-level analysis:
```mermaid
flowchart TD;
    A[Create phyloseq object]-->B[Access omnicrobe database];
    A-->C[Raw dataset plots];
    B-->D[Assess and remove contaminants];
    D-->E[Filtered dataset plots];
    D-->F[Sample summary plots];
    D-->G[Hypothesis testing];
    D-->H[Core microbiome];
```

Analysis workflow for metagenomes-assembled genomes analysis:
```mermaids
flowchart TD;
    A[Prep MAG tree and metadata]-->B[Plot MAG tree]
    A-->C[Access omnicrobe database]
    A-->D[Codiversification tests]
```

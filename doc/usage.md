# GFFUtils Usages
Rui Yang
2024-01-30

- [<span class="toc-section-number">1</span>
  <span class="header-section-number">1</span> APIs](#apis)
- [<span class="toc-section-number">2</span>
  <span class="header-section-number">2</span> Example
  usages](#example-usages)
  - [<span class="toc-section-number">2.1</span>
    <span class="header-section-number">2.1</span> `using` related
    packages](#using-related-packages)
  - [<span class="toc-section-number">2.2</span>
    <span class="header-section-number">2.2</span> Read in GFF3
    file](#read-in-gff3-file)
    - [<span class="toc-section-number">2.2.1</span>
      <span class="header-section-number">2.2.1</span> Extract extended
      TSS regions](#extract-extended-tss-regions)
    - [<span class="toc-section-number">2.2.2</span>
      <span class="header-section-number">2.2.2</span> Extract extended
      TES regions](#extract-extended-tes-regions)
    - [<span class="toc-section-number">2.2.3</span>
      <span class="header-section-number">2.2.3</span> Extract extended
      feature body regions](#extract-extended-feature-body-regions)
  - [<span class="toc-section-number">2.3</span>
    <span class="header-section-number">2.3</span> Sort BED
    file](#sort-bed-file)

## <span class="header-section-number">1</span> APIs

- `GFFBox`

- `gff_read`

- `bed_write`

- `show`

- `tss`

- `tes`

- `region`

- `bed_clip`

- `bed_sort`

- `faidx_read`

- `bed_read`

For detailed usages, type `?<name>` in the Julia REPL.

## <span class="header-section-number">2</span> Example usages

### <span class="header-section-number">2.1</span> `using` related packages

``` julia
using GFFUtils
using DataFrames
```

    [ Info: Precompiling GFFUtils [9b58c212-aae8-41f9-86da-aeb22b670f37]

### <span class="header-section-number">2.2</span> Read in GFF3 file

``` julia
filename = "E:/BaiduSyncdisk/AutoBackup/Packages/Julia/GFFUtils/exampledata/example.gff3"

gffbox = gff_read(filename)
## gffbox is an object of type GFFBox with four fields:
# features: a DataFrame object storing features;
# sequences: a DataFrame object storing FASTA sequences;
# seqinfo: a DataFrame object storing chromosome ranges;
# buildinfo: a Dict object storing various metadata.
```

    Features: 
    20 features in total.
    24 fields in total: seqid, source, type, start, end, score, strand, phase, ID, Name, biotype, description, gene_id, logic_name, version, Parent, tag, transcript_id, constitutive, ensembl_end_phase, ensembl_phase, exon_id, rank, protein_id.

    Sequences: 
    2 sequences in total.

    Sequence info: 
    5 sequences in total.

    Build info: 

    genebuild-last-updated: 2022-03
    genome-build-accession: GCA_015227675.2
    genome-build: mRatBN7.2
    genome-date: 2020-11
    gff-version: 3
    genome-version: mRatBN7.2

``` julia
show(gffbox.features)
```

    20×24 DataFrame
     Row │ seqid    source   type             start   end     score     strand  ph ⋯
         │ String?  String?  String?          Int64?  Int64?  Float64?  Char?   In ⋯
    ─────┼──────────────────────────────────────────────────────────────────────────
       1 │ 1        ensembl  gene             455309  474795   missing  +       mi ⋯
       2 │ 1        ensembl  mRNA             455309  474795   missing  +       mi
       3 │ 1        ensembl  five_prime_UTR   455309  455393   missing  +       mi
       4 │ 1        ensembl  exon             455309  455396   missing  +       mi
       5 │ 1        ensembl  CDS              455394  455396   missing  +          ⋯
       6 │ 1        ensembl  exon             470250  470376   missing  +       mi
       7 │ 1        ensembl  CDS              470250  470376   missing  +
       8 │ 1        ensembl  exon             470563  470623   missing  +       mi
       9 │ 1        ensembl  CDS              470563  470623   missing  +          ⋯
      10 │ 1        ensembl  CDS              471754  473716   missing  +
      11 │ 1        ensembl  exon             471754  474795   missing  +       mi
      12 │ 1        ensembl  three_prime_UTR  473717  474795   missing  +       mi
      13 │ 10       ensembl  gene             610774  613518   missing  -       mi ⋯
      14 │ 10       ensembl  mRNA             610774  613518   missing  -       mi
      15 │ 10       ensembl  exon             610774  612004   missing  -       mi
      16 │ 10       ensembl  CDS              610774  612004   missing  -
      17 │ 10       ensembl  exon             613233  613293   missing  -       mi ⋯
      18 │ 10       ensembl  CDS              613233  613293   missing  -
      19 │ 10       ensembl  exon             613485  613518   missing  -       mi
      20 │ 10       ensembl  CDS              613485  613518   missing  -
                                                                  17 columns omitted

``` julia
show(gffbox.sequences)
```

    2×4 DataFrame
     Row │ seqid    description  seqsize  sequence                          
         │ String?  String?      UInt64?  String?                           
    ─────┼──────────────────────────────────────────────────────────────────
       1 │ seq1     seq1             100  cttctgggcgtacccgattctcggagaacttg…
       2 │ seq2     seq2             100  ttcaagtgctcagtcaatgtgattcacagtat…

``` julia
show(gffbox.seqinfo)
```

    5×3 DataFrame
     Row │ seqid    start   end       
         │ String?  Int64?  Int64?    
    ─────┼────────────────────────────
       1 │ 1             1  260522016
       2 │ 2             1  249053267
       3 │ 10            1    1072111
       4 │ X             1  152453651
       5 │ Y             1   18315841

``` julia
show(gffbox.buildinfo)
```

    Dict("genebuild-last-updated" => "2022-03", "genome-build-accession" => "GCA_015227675.2", "genome-build" => "mRatBN7.2", "genome-date" => "2020-11", "gff-version" => "3", "genome-version" => "mRatBN7.2")

#### <span class="header-section-number">2.2.1</span> Extract extended TSS regions

``` julia
tsss = tss(gffbox, ["gene"], 100000000, 100)
show(tsss)
```

    2×6 DataFrame
     Row │ chrom    chromStart  chromEnd  name                score     strand 
         │ String?  Int64       Int64     String              Float64?  Char?  
    ─────┼─────────────────────────────────────────────────────────────────────
       1 │ 1                 1    455409  ENSRNOG00000065285   missing  +
       2 │ 10           613418   1072111  ENSRNOG00000069636   missing  -

#### <span class="header-section-number">2.2.2</span> Extract extended TES regions

``` julia
tess = tes(gffbox, ["gene"], 100000000, 100)
show(tess)
```

    2×6 DataFrame
     Row │ chrom    chromStart  chromEnd  name                score     strand 
         │ String?  Int64       Int64     String              Float64?  Char?  
    ─────┼─────────────────────────────────────────────────────────────────────
       1 │ 1                 1    474895  ENSRNOG00000065285   missing  +
       2 │ 10           610674   1072111  ENSRNOG00000069636   missing  -

#### <span class="header-section-number">2.2.3</span> Extract extended feature body regions

``` julia
regions = region(gffbox, ["gene"], 100000000, 100)
show(regions)
```

    2×6 DataFrame
     Row │ chrom    chromStart  chromEnd  name                score     strand 
         │ String?  Int64       Int64     String              Float64?  Char?  
    ─────┼─────────────────────────────────────────────────────────────────────
       1 │ 1                 1    474895  ENSRNOG00000065285   missing  +
       2 │ 10           610674   1072111  ENSRNOG00000069636   missing  -

### <span class="header-section-number">2.3</span> Sort BED file

``` julia
df = DataFrame(
    "seqid" => ["1", "10", "2", "contig10", "contig6", "X", "Y", "2"],
    "start" => [100, 100, 100, 100, 100, 100, 100, 50],
    "end" => [890, 666, 500, 1000, 300, 1000, 250, 2000])

srt_df = bed_sort(df, sort_cols=["seqid", "start", "end"])  # Based on the first three columns in a BED file
show(srt_df)
```

    8×3 DataFrame
     Row │ seqid     start  end   
         │ String    Int64  Int64 
    ─────┼────────────────────────
       1 │ 1           100    890
       2 │ 2            50   2000
       3 │ 2           100    500
       4 │ 10          100    666
       5 │ X           100   1000
       6 │ Y           100    250
       7 │ contig6     100    300
       8 │ contig10    100   1000

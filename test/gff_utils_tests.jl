@testset "Utils tests" begin
    test_input = """
##gff-version 3
##sequence-region 1 1 260522016
##sequence-region 2 1 249053267
##sequence-region 10 1 1072111
##sequence-region X 1 152453651
##sequence-region Y 1 18315841
#!genome-build mRatBN7.2
#!genome-version mRatBN7.2
#!genome-date 2020-11
#!genome-build-accession GCA_015227675.2
#!genebuild-last-updated 2022-03
1	ensembl	gene	455309	474795	.	+	.	ID=gene:ENSRNOG00000065285;Name=LOC102548633;biotype=protein_coding;description=zinc finger protein 729-like [Source:RGD Symbol%3BAcc:7713719];gene_id=ENSRNOG00000065285;logic_name=ensembl;version=1
1	ensembl	mRNA	455309	474795	.	+	.	ID=transcript:ENSRNOT00000117015;Parent=gene:ENSRNOG00000065285;Name=LOC102548633-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=ENSRNOT00000117015;version=1
1	ensembl	five_prime_UTR	455309	455393	.	+	.	Parent=transcript:ENSRNOT00000117015
1	ensembl	exon	455309	455396	.	+	.	Parent=transcript:ENSRNOT00000117015;Name=ENSRNOE00000625275;constitutive=1;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSRNOE00000625275;rank=1;version=1
1	ensembl	CDS	455394	455396	.	+	0	ID=CDS:ENSRNOP00000090915;Parent=transcript:ENSRNOT00000117015;protein_id=ENSRNOP00000090915
1	ensembl	exon	470250	470376	.	+	.	Parent=transcript:ENSRNOT00000117015;Name=ENSRNOE00000661823;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSRNOE00000661823;rank=2;version=1
1	ensembl	CDS	470250	470376	.	+	0	ID=CDS:ENSRNOP00000090915;Parent=transcript:ENSRNOT00000117015;protein_id=ENSRNOP00000090915
1	ensembl	exon	470563	470623	.	+	.	Parent=transcript:ENSRNOT00000117015;Name=ENSRNOE00000583675;constitutive=1;ensembl_end_phase=2;ensembl_phase=1;exon_id=ENSRNOE00000583675;rank=3;version=1
1	ensembl	CDS	470563	470623	.	+	2	ID=CDS:ENSRNOP00000090915;Parent=transcript:ENSRNOT00000117015;protein_id=ENSRNOP00000090915
1	ensembl	CDS	471754	473716	.	+	1	ID=CDS:ENSRNOP00000090915;Parent=transcript:ENSRNOT00000117015;protein_id=ENSRNOP00000090915
1	ensembl	exon	471754	474795	.	+	.	Parent=transcript:ENSRNOT00000117015;Name=ENSRNOE00000641146;constitutive=1;ensembl_end_phase=-1;ensembl_phase=2;exon_id=ENSRNOE00000641146;rank=4;version=1
1	ensembl	three_prime_UTR	473717	474795	.	+	.	Parent=transcript:ENSRNOT00000117015
###
10	ensembl	gene	610774	613518	.	-	.	ID=gene:ENSRNOG00000069636;Name=Zfp748-ps3;biotype=protein_coding;description=zinc finger protein 748%2C pseudogene 3 [Source:RGD Symbol%3BAcc:40952346];gene_id=ENSRNOG00000069636;logic_name=ensembl;version=1
10	ensembl	mRNA	610774	613518	.	-	.	ID=transcript:ENSRNOT00000106239;Parent=gene:ENSRNOG00000069636;Name=Zfp748-ps3-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=ENSRNOT00000106239;version=1
10	ensembl	exon	610774	612004	.	-	.	Parent=transcript:ENSRNOT00000106239;Name=ENSRNOE00000625142;constitutive=1;ensembl_end_phase=0;ensembl_phase=2;exon_id=ENSRNOE00000625142;rank=3;version=1
10	ensembl	CDS	610774	612004	.	-	1	ID=CDS:ENSRNOP00000079178;Parent=transcript:ENSRNOT00000106239;protein_id=ENSRNOP00000079178
10	ensembl	exon	613233	613293	.	-	.	Parent=transcript:ENSRNOT00000106239;Name=ENSRNOE00000598372;constitutive=1;ensembl_end_phase=2;ensembl_phase=1;exon_id=ENSRNOE00000598372;rank=2;version=1
10	ensembl	CDS	613233	613293	.	-	2	ID=CDS:ENSRNOP00000079178;Parent=transcript:ENSRNOT00000106239;protein_id=ENSRNOP00000079178
10	ensembl	exon	613485	613518	.	-	.	Parent=transcript:ENSRNOT00000106239;Name=ENSRNOE00000486605;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSRNOE00000486605;rank=1;version=2
10	ensembl	CDS	613485	613518	.	-	0	ID=CDS:ENSRNOP00000079178;Parent=transcript:ENSRNOT00000106239;protein_id=ENSRNOP00000079178
##FASTA
>seq1
cttctgggcgtacccgattctcggagaacttgccgcaccattccgccttg
tgttcattgctgcctgcatgttcattgtctacctcggctacgtgtggcta
>seq2
ttcaagtgctcagtcaatgtgattcacagtatgtcaccaaatattttggc
agctttctcaagggatcaaaattatggatcattatggaatacctcggtgg
"""

    gffbox = gff_read(IOBuffer(test_input))

    ## Test tss
    tsss = tss(gffbox, "gene", 100000000, 100)

    tsss_column_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    @test all(occursin.(names(tsss), tsss_column_names))

    tsss_content = ["1" "1" "455409" "ENSRNOG00000065285" "missing" "+"; "10" "613418" "1072111" "ENSRNOG00000069636" "missing" "-"]
    @test all(occursin.(string.(Array(tsss)), tsss_content))

    ## Test tes
    tess = tes(gffbox, "gene", 100000000, 100)

    tess_column_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    @test all(occursin.(names(tess), tess_column_names))

    tess_content = ["1" "1" "474895" "ENSRNOG00000065285" "missing" "+"; "10" "610674" "1072111" "ENSRNOG00000069636" "missing" "-"]
    @test all(occursin.(string.(Array(tess)), tess_content))

    ## Test region
    regions = region(gffbox, "gene", 100000000, 100)

    regions_column_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    @test all(occursin.(names(regions), regions_column_names))

    regions_content = ["1" "1" "474895" "ENSRNOG00000065285" "missing" "+"; "10" "610674" "1072111" "ENSRNOG00000069636" "missing" "-"]
    @test all(occursin.(string.(Array(regions)), regions_content))

    ## Test bed_clip
    # Test for bed_clip has been contained in the above three tests

    ## Test bed_sort
    df = DataFrame(
        "seqid" => ["1", "10", "2", "contig10", "contig6", "X", "Y", "2"],
        "start" => [100, 100, 100, 100, 100, 100, 100, 50],
        "end" => [890, 666, 500, 1000, 300, 1000, 250, 2000])

    srt_df = bed_sort(df, sort_cols=["seqid", "start", "end"])

    srt_df_column_names = ["seqid", "start", "end"]
    @test all(occursin.(names(srt_df), srt_df_column_names))

    srt_df_content = ["1" "100" "890"; "2" "50" "2000"; "2" "100" "500"; "10" "100" "666"; "X" "100" "1000"; "Y" "100" "250"; "contig6" "100" "300"; "contig10" "100" "1000"]
    @test all(occursin.(string.(Array(srt_df)), srt_df_content))
end
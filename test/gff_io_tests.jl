@testset "IO tests" begin
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

    ## Test gff_read
    # Test features
    features_column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "ID", "Name", "biotype", "description", "gene_id", "logic_name", "version", "Parent", "tag", "transcript_id", "constitutive", "ensembl_end_phase", "ensembl_phase", "exon_id", "rank", "protein_id"]
    @test all(occursin.(names(gffbox.features), features_column_names))

    features_content = ["1" "ensembl" "gene" "455309" "474795" "missing" "+" "missing" "gene:ENSRNOG00000065285" "LOC102548633" "protein_coding" "zinc finger protein 729-like [Source:RGD Symbol;Acc:7713719]" "ENSRNOG00000065285" "ensembl" "1" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing"; "1" "ensembl" "mRNA" "455309" "474795" "missing" "+" "missing" "transcript:ENSRNOT00000117015" "LOC102548633-201" "protein_coding" "missing" "missing" "missing" "1" "gene:ENSRNOG00000065285" "Ensembl_canonical" "ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing"; "1" "ensembl" "five_prime_UTR" "455309" "455393" "missing" "+" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing"; "1" "ensembl" "exon" "455309" "455396" "missing" "+" "missing" "missing" "ENSRNOE00000625275" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000117015" "missing" "missing" "1" "0" "-1" "ENSRNOE00000625275" "1" "missing"; "1" "ensembl" "CDS" "455394" "455396" "missing" "+" "0" "CDS:ENSRNOP00000090915" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000090915"; "1" "ensembl" "exon" "470250" "470376" "missing" "+" "missing" "missing" "ENSRNOE00000661823" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000117015" "missing" "missing" "1" "1" "0" "ENSRNOE00000661823" "2" "missing"; "1" "ensembl" "CDS" "470250" "470376" "missing" "+" "0" "CDS:ENSRNOP00000090915" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000090915"; "1" "ensembl" "exon" "470563" "470623" "missing" "+" "missing" "missing" "ENSRNOE00000583675" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000117015" "missing" "missing" "1" "2" "1" "ENSRNOE00000583675" "3" "missing"; "1" "ensembl" "CDS" "470563" "470623" "missing" "+" "2" "CDS:ENSRNOP00000090915" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000090915"; "1" "ensembl" "CDS" "471754" "473716" "missing" "+" "1" "CDS:ENSRNOP00000090915" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000090915"; "1" "ensembl" "exon" "471754" "474795" "missing" "+" "missing" "missing" "ENSRNOE00000641146" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000117015" "missing" "missing" "1" "-1" "2" "ENSRNOE00000641146" "4" "missing"; "1" "ensembl" "three_prime_UTR" "473717" "474795" "missing" "+" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000117015" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing"; "10" "ensembl" "gene" "610774" "613518" "missing" "-" "missing" "gene:ENSRNOG00000069636" "Zfp748-ps3" "protein_coding" "zinc finger protein 748, pseudogene 3 [Source:RGD Symbol;Acc:40952346]" "ENSRNOG00000069636" "ensembl" "1" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "missing"; "10" "ensembl" "mRNA" "610774" "613518" "missing" "-" "missing" "transcript:ENSRNOT00000106239" "Zfp748-ps3-201" "protein_coding" "missing" "missing" "missing" "1" "gene:ENSRNOG00000069636" "Ensembl_canonical" "ENSRNOT00000106239" "missing" "missing" "missing" "missing" "missing" "missing"; "10" "ensembl" "exon" "610774" "612004" "missing" "-" "missing" "missing" "ENSRNOE00000625142" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000106239" "missing" "missing" "1" "0" "2" "ENSRNOE00000625142" "3" "missing"; "10" "ensembl" "CDS" "610774" "612004" "missing" "-" "1" "CDS:ENSRNOP00000079178" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000106239" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000079178"; "10" "ensembl" "exon" "613233" "613293" "missing" "-" "missing" "missing" "ENSRNOE00000598372" "missing" "missing" "missing" "missing" "1" "transcript:ENSRNOT00000106239" "missing" "missing" "1" "2" "1" "ENSRNOE00000598372" "2" "missing"; "10" "ensembl" "CDS" "613233" "613293" "missing" "-" "2" "CDS:ENSRNOP00000079178" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000106239" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000079178"; "10" "ensembl" "exon" "613485" "613518" "missing" "-" "missing" "missing" "ENSRNOE00000486605" "missing" "missing" "missing" "missing" "2" "transcript:ENSRNOT00000106239" "missing" "missing" "1" "1" "0" "ENSRNOE00000486605" "1" "missing"; "10" "ensembl" "CDS" "613485" "613518" "missing" "-" "0" "CDS:ENSRNOP00000079178" "missing" "missing" "missing" "missing" "missing" "missing" "transcript:ENSRNOT00000106239" "missing" "missing" "missing" "missing" "missing" "missing" "missing" "ENSRNOP00000079178"]
    @test all(occursin.(string.(Array(gffbox.features)), features_content))

    # Test sequences
    sequences_column_names = ["seqid", "description", "seqsize", "sequence"]
    @test all(occursin.(names(gffbox.sequences), sequences_column_names))

    sequences_content = ["seq1" "seq1" "100" "cttctgggcgtacccgattctcggagaacttgccgcaccattccgccttgtgttcattgctgcctgcatgttcattgtctacctcggctacgtgtggcta"; "seq2" "seq2" "100" "ttcaagtgctcagtcaatgtgattcacagtatgtcaccaaatattttggcagctttctcaagggatcaaaattatggatcattatggaatacctcggtgg"]
    @test all(occursin.(string.(Array(gffbox.sequences)), sequences_content))

    # Test seqinfo
    seqinfo_column_names = ["seqid", "start", "end"]
    @test all(occursin.(names(gffbox.seqinfo), seqinfo_column_names))

    seqinfo_content = ["1" "1" "260522016"; "2" "1" "249053267"; "10" "1" "1072111"; "X" "1" "152453651"; "Y" "1" "18315841"]
    @test all(occursin.(string.(Array(gffbox.seqinfo)), seqinfo_content))

    # Test buildinfo
    buildinfo = Dict("genebuild-last-updated" => "2022-03", "genome-build-accession" => "GCA_015227675.2", "genome-build" => "mRatBN7.2", "genome-date" => "2020-11", "gff-version" => "3", "genome-version" => "mRatBN7.2")
    @test gffbox.buildinfo == buildinfo
end
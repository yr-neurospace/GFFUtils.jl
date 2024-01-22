"""
    GFFBox(features::Union{DataFrame,Missing}, sequences::Union{DataFrame,Missing}, seqinfo::Union{DataFrame,Missing}, buildinfo::Union{Dict,Missing}) -> GFFBox

Type `GFFBox`.

You can construct the object of type `GFFBox` by passing either one, two, three, or four arguments to `GFFBox()`.

The easiest way to do this is to use the `gff_read(filename::AbstractString)` function.

See also `gff_read`.
"""
mutable struct GFFBox
    features::Union{DataFrame,Missing}
    sequences::Union{DataFrame,Missing}
    seqinfo::Union{DataFrame,Missing}
    buildinfo::Union{Dict,Missing}

    function GFFBox()
        new(missing, missing, missing, missing)
    end

    function GFFBox(features::Union{DataFrame,Missing})
        new(features, missing, missing, missing)
    end

    function GFFBox(features::Union{DataFrame,Missing}, sequences::Union{DataFrame,Missing})
        new(features, sequences, missing, missing)
    end

    function GFFBox(features::Union{DataFrame,Missing}, sequences::Union{DataFrame,Missing}, seqinfo::Union{DataFrame,Missing})
        new(features, sequences, seqinfo, missing)
    end

    function GFFBox(features::Union{DataFrame,Missing}, sequences::Union{DataFrame,Missing}, seqinfo::Union{DataFrame,Missing}, buildinfo::Union{Dict,Missing})
        new(features, sequences, seqinfo, buildinfo)
    end
end

"""
    gff_read(filename::AbstractString) -> GFFBox

Read a **GFF3** file and return an object of type `GFFBox`.
"""
function gff_read(filename::AbstractString)
    gff_read(open(filename))
end

"""
    gff_read(input::IO) -> GFFBox

Read a **GFF3** file and return an object of type `GFFBox`.
"""
function gff_read(input::IO)
    reader = GFF3.Reader(input, save_directives=true, skip_directives=false, skip_comments=false)

    gffbox = GFFBox(
        DataFrame(
            "seqid" => Union{String,Missing}[],
            "source" => Union{String,Missing}[],
            "type" => Union{String,Missing}[],
            "start" => Union{Int64,Missing}[],
            "end" => Union{Int64,Missing}[],
            "score" => Union{Float64,Missing}[],
            "strand" => Union{Char,Missing}[],
            "phase" => Union{UInt8,Missing}[]),
        DataFrame(
            "seqid" => Union{String,Missing}[],
            "description" => Union{String,Missing}[],
            "seqsize" => Union{UInt64,Missing}[],
            "sequence" => Union{String,Missing}[]),
        DataFrame(
            "seqid" => Union{String,Missing}[],
            "start" => Union{Int64,Missing}[],
            "end" => Union{Int64,Missing}[]),
        Dict{String,String}()
    )

    for record in reader
        if GFF3.isfeature(record)  # Parse features
            record_fields = [
                "seqid" => GFF3.hasseqid(record) ? GFF3.seqid(record) : missing,
                "source" => GFF3.hassource(record) ? GFF3.source(record) : missing,
                "type" => GFF3.hasfeaturetype(record) ? GFF3.featuretype(record) : missing,
                "start" => GFF3.hasseqstart(record) ? GFF3.seqstart(record) : missing,
                "end" => GFF3.hasseqend(record) ? GFF3.seqend(record) : missing,
                "score" => GFF3.hasscore(record) ? GFF3.score(record) : missing,
                "strand" => GFF3.hasstrand(record) ? convert(Char, GFF3.strand(record)) : missing,
                "phase" => GFF3.hasphase(record) ? GFF3.phase(record) : missing
            ]
            record_attributes = GFF3.attributes(record)

            record_attributes_max_length = max(map(x -> length(x.second), record_attributes)...)
            if record_attributes_max_length > 1
                record_fields = map(record_fields) do x
                    x.first => repeat([x.second], record_attributes_max_length)
                end
                record_attributes = map(record_attributes) do x
                    x.first => convert(Vector{Union{String,Missing}}, x.second)
                end
                record_attributes = map(record_attributes) do x
                    x.first => [x.second; repeat([missing], record_attributes_max_length - length(x.second))]
                end
            end
            append!(gffbox.features, DataFrame([record_fields; record_attributes]), cols=:union)
        elseif GFF3.isdirective(record)  # Parse directives, starting with double hashes "##"
            seq_region_m = match(r"^##sequence-region .+", string(record))
            if !isnothing(seq_region_m)
                seq_region_strs = split(seq_region_m.match, r" +")
                append!(gffbox.seqinfo,
                    DataFrame(
                        "seqid" => String(seq_region_strs[2]),
                        "start" => parse(Int64, seq_region_strs[3]),
                        "end" => parse(Int64, seq_region_strs[4])),
                    cols=:setequal)
                continue
            end

            gff_version_m = match(r"^##gff-version .+", string(record))
            if !isnothing(gff_version_m)
                gff_version_strs = split(gff_version_m.match, r" +")
                gffbox.buildinfo[replace(gff_version_strs[1], r"^##" => "")] = gff_version_strs[2]
            end
        elseif GFF3.iscomment(record)  # Parse comments, starting with a single "#"
            c_genome_build_m = match(r"#!genome-build .+", string(record))
            if !isnothing(c_genome_build_m)
                c_genome_build_strs = split(c_genome_build_m.match, r" +")
                gffbox.buildinfo[replace(c_genome_build_strs[1], r"^#!" => "")] = c_genome_build_strs[2]
            end

            c_genome_version_m = match(r"#!genome-version .+", string(record))
            if !isnothing(c_genome_version_m)
                c_genome_version_strs = split(c_genome_version_m.match, r" +")
                gffbox.buildinfo[replace(c_genome_version_strs[1], r"^#!" => "")] = c_genome_version_strs[2]
            end

            c_genome_date_m = match(r"#!genome-date .+", string(record))
            if !isnothing(c_genome_date_m)
                c_genome_date_strs = split(c_genome_date_m.match, r" +")
                gffbox.buildinfo[replace(c_genome_date_strs[1], r"^#!" => "")] = c_genome_date_strs[2]
            end

            c_genome_build_accession_m = match(r"#!genome-build-accession .+", string(record))
            if !isnothing(c_genome_build_accession_m)
                c_genome_build_accession_strs = split(c_genome_build_accession_m.match, r" +")
                gffbox.buildinfo[replace(c_genome_build_accession_strs[1], r"^#!" => "")] = c_genome_build_accession_strs[2]
            end

            c_genebuild_last_updated_m = match(r"#!genebuild-last-updated .+", string(record))
            if !isnothing(c_genebuild_last_updated_m)
                c_genebuild_last_updated_strs = split(c_genebuild_last_updated_m.match, r" +")
                gffbox.buildinfo[replace(c_genebuild_last_updated_strs[1], r"^#!" => "")] = c_genebuild_last_updated_strs[2]
            end
        end
    end

    # Parse FASTA sequences, after the line of "###FASTA"
    if GFF3.hasfasta(reader)
        fa_reader = GFF3.getfasta(reader)
        for fa_record in fa_reader
            append!(gffbox.sequences,
                DataFrame(
                    "seqid" => FASTX.identifier(fa_record),
                    "description" => replace(FASTX.description(fa_record), Regex(FASTX.identifier(fa_record) * " +") => ""),
                    "seqsize" => FASTX.seqsize(fa_record),
                    "sequence" => FASTX.sequence(fa_record)),
                cols=:setequal)
        end
    end

    close(reader)

    gffbox
end

"""
    bed_write(filename::AbstractString, bed::DataFrame; need_sort::Bool=true, sort_cols::Vector{String}=["chrom", "chromStart", "chromEnd"]) -> file

Write `bed` `DataFrame` into a file.

For more details about sorting `bed` `DataFrame`, see `bed_sort`.
"""
function bed_write(filename::AbstractString, bed::DataFrame; need_sort::Bool=true, sort_cols::Vector{String}=["chrom", "chromStart", "chromEnd"])
    if need_sort
        bed = bed_sort(bed, sort_cols=sort_cols)
    end

    CSV.write(filename, bed, delim="\t", append=false, writeheader=false)
end

"""
    show(io::IO, gffbox::GFFBox)

Show the object of type `GFFBox`.
"""
function show(io::IO, gffbox::GFFBox)
    println(io, "Features: ")
    println(io, nrow(gffbox.features), " features in total.")
    println(io, length(names(gffbox.features)), " fields in total: ", join(names(gffbox.features), ", "), ".")
    print(io, "\n")

    println(io, "Sequences: ")
    println(io, nrow(gffbox.sequences), " sequences in total.")
    print(io, "\n")

    println(io, "Sequence info: ")
    println(io, nrow(gffbox.seqinfo), " sequences in total.")
    print(io, "\n")

    println(io, "Build info: ")
    for (k, v) in gffbox.buildinfo
        println(k, ": ", v)
    end
end

"""
    faidx_read(filename::AbstractString) -> DataFrame

Read FASTA index file to extract the chromosome regions.

**Note:** the first **two** columns are mandatory while the others are ignored.
"""
function faidx_read(filename::AbstractString)
    faidx_read(open(filename))
end

"""
    faidx_read(input::IO) -> DataFrame

Read FASTA index file to extract the chromosome regions.

**Note:** the first **two** columns are mandatory while the others are ignored.
"""
function faidx_read(input::IO)
    df = CSV.read(input, DataFrame, header=false)
    df_bed = select(df,
        :Column1 => :chrom,
        :Column1 => (x -> repeat([1], length(x))) => :chromStart,
        :Column2 => :chromEnd)
    srt_df = bed_sort(df_bed, sort_cols=["chrom", "chromStart", "chromEnd"])

    srt_df
end

"""
    bed_read(filenname::AbstractString, need_sort::Bool=true) -> DataFrame

Read BED file with optional sorting option based on the first three columns.
"""
function bed_read(filenname::AbstractString; need_sort::Bool=true)
    bed_read(open(filenname), need_sort=need_sort)
end

"""
    bed_read(input::IO; need_sort::Bool=true) -> DataFrame

Read BED file with optional sorting option based on the first three columns.
"""
function bed_read(input::IO; need_sort::Bool=true)
    df = CSV.read(input, DataFrame, header=false)
    if need_sort
        srt_df = bed_sort(df, sort_cols=["Column1", "Column2", "Column3"])
    else
        srt_df = df
    end

    srt_df
end
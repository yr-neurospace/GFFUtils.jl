"""
    tss(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true) -> DataFrame

Extract transcription start sites (TSSs) of `level` ("gene" by default) `upstream`bp (0bp by default) upstream and `downstream`bp (0bp by default) downstream.

If `clip=true` (default), then clip extended TSS regions based on chromosome lengthes, which should have been stored in the field `seqinfo` of the object of type `GFFBox`.
"""
function tss(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true)
    tsss = @chain gffbox.features begin
        @subset :type .== level
        @select begin
            :chrom = :seqid
            $([:start, :end, :strand] => ByRow((x, y, z) -> z == '-' ? [y - downstream, y + upstream] : [x - upstream, x + downstream]) => [:chromStart, :chromEnd])
            :name = replace.(:ID, r"^\w+:" => "")
            :score
            :strand
        end
        @transform $([:chromStart, :chromEnd] => ByRow((x, y) -> Interval{Int64,Closed,Closed}(x, y)) => :interval)
    end

    if clip
        clipped_tsss = bed_clip(tsss, gffbox.seqinfo, bed_cols=["chrom", "chromStart", "chromEnd"])
    else
        clipped_tsss = tsss
    end

    clipped_tsss
end

"""
    tes(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true) -> DataFrame

Extract transcription end sites (TESs) of `level` ("gene" by default) `upstream`bp (0bp by default) upstream and `downstream`bp (0bp by default) downstream.

If `clip=true` (default), then clip extended TES regions based on chromosome lengthes, which should have been stored in the field `seqinfo` of the object of type `GFFBox`.
"""
function tes(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true)
    tess = @chain gffbox.features begin
        @subset :type .== level
        @select begin
            :chrom = :seqid
            $([:start, :end, :strand] => ByRow((x, y, z) -> z == '-' ? [x - downstream, x + upstream] : [y - upstream, y + downstream]) => [:chromStart, :chromEnd])
            :name = replace.(:ID, r"^\w+:" => "")
            :score
            :strand
        end
        @transform $([:chromStart, :chromEnd] => ByRow((x, y) -> Interval{Int64,Closed,Closed}(x, y)) => :interval)
    end

    if clip
        clipped_tess = bed_clip(tess, gffbox.seqinfo, bed_cols=["chrom", "chromStart", "chromEnd"])
    else
        clipped_tess = tess
    end

    clipped_tess
end

"""
    region(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true) -> DataFrame

    Extract regions of `level` ("gene" by default) `upstream`bp (0bp by default) upstream and `downstream`bp (0bp by default) downstream.

    If `clip=true` (default), then clip extended regions based on chromosome lengthes, which should have been stored in the field `seqinfo` of the object of type `GFFBox`.
"""
function region(gffbox::GFFBox, level::AbstractString="gene", upstream::Int=0, downstream::Int=0; clip::Bool=true)
    regions = @chain gffbox.features begin
        @subset :type .== level
        @select begin
            :chrom = :seqid
            $([:start, :end, :strand] => ByRow((x, y, z) -> z == '-' ? [x - downstream, y + upstream] : [x - upstream, y + downstream]) => [:chromStart, :chromEnd])
            :name = replace.(:ID, r"^\w+:" => "")
            :score
            :strand
        end
        @transform $([:chromStart, :chromEnd] => ByRow((x, y) -> Interval{Int64,Closed,Closed}(x, y)) => :interval)
    end

    if clip
        clipped_regions = bed_clip(regions, gffbox.seqinfo, bed_cols=["chrom", "chromStart", "chromEnd"])
    else
        clipped_regions = regions
    end

    clipped_regions
end

"""
    bed_clip(bed::DataFrame, seqinfo::DataFrame; bed_cols::Vector{String}=["chrom", "chromStart", "chromEnd"]) -> DataFrame

Clip regions in the `bed` based on regions in the `seqinfo`.

Both must have at least these three columns: "seqid", "start", and "end", but the names of those three columns don't need to be exactly the same as "seqid", "start", and "end".

For columns in `bed`, you **must** give those three column names you used.

For columns in `seqinfo`, the **order** of the columns must be exactly the same as "seqid", "start", and "end" regardless of what the names are.
"""
function bed_clip(bed::DataFrame, seqinfo::DataFrame; bed_cols::Vector{String}=["chrom", "chromStart", "chromEnd"])
    sym_bed_cols = Symbol.(bed_cols)
    bed_interval = transform(bed, [sym_bed_cols[2], sym_bed_cols[3]] => ByRow((x, y) -> Interval{Int64,Closed,Closed}(x, y)) => :interval)

    sym_seqinfo_cols = propertynames(seqinfo)
    seqinfo_dict = @chain seqinfo begin
        @select $([sym_seqinfo_cols[1], sym_seqinfo_cols[2], sym_seqinfo_cols[3]] => ByRow((id, x, y) -> id => Interval{Int64,Closed,Closed}(x, y)) => :interval)
        Dict(_.interval)
    end

    clipped_bed = @chain bed_interval begin
        @transform $([sym_bed_cols[1], :interval] => ByRow((id, i) -> intersect(i, seqinfo_dict[id])) => :interval)
        @transform $(:interval => ByRow(x -> [x.first, x.last]) => [sym_bed_cols[2], sym_bed_cols[3]])
        @select $(Not(:interval))
    end

    clipped_bed
end

"""
    bed_sort(bed::DataFrame, sort_cols::Vector{String}=["chrom", "chromStart", "chromEnd"]) -> DataFrame

Sort `bed` based on `chrom`, `chromStart`, and `chromEnd` by the function `NaturalSort.natural`.

**Note:** columns given in `sort_cols` must be able to be parsed into `Int` except for the 1st one.
"""
function bed_sort(bed::DataFrame; sort_cols::Vector{String}=["chrom", "chromStart", "chromEnd"])
    sym_sort_cols = Symbol.(sort_cols)

    for sym_col in sym_sort_cols
        bed[!, sym_col] = string.(bed[!, sym_col])
    end

    srt_bed = sort(bed, sym_sort_cols, lt=NaturalSort.natural)

    for sym_col in sym_sort_cols[2:end]
        srt_bed[!, sym_col] = parse.(Int64, srt_bed[!, sym_col])
    end

    srt_bed
end
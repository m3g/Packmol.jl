#
# Renders images of packed systems with VMD (offscreen, with the Tachyon renderer), for
# the tutorials of the docs (docs/src/tutorials/images). Requires `vmd` and ImageMagick's
# `convert` in the PATH.
#
# Each representation is a NamedTuple with the fields of a VMD representation:
#
#     (selection = "water", style = "Lines", color = "Name", material = "AOChalky")
#
# `rotate` is a list of rotations applied to the view, e.g. [(:x, -90), (:y, 30)], and
# `cell` draws the periodic cell of the CRYST1 record (pbctools), centered at the origin
# (`cell = :origin`) or with a corner at the origin (`cell = :unitcell`).
#
const default_material = "AOChalky"

function render(pdbfile::AbstractString, imagefile::AbstractString;
    reps::AbstractVector,
    rotate = Tuple{Symbol,Real}[],
    zoom::Real = 1.0,
    size::Tuple{Int,Int} = (1600, 1200),
    cell::Union{Nothing,Symbol} = nothing,
    clip::Union{Nothing,Tuple{Real,Real}} = nothing,
    width::Int = 1000,
)
    tcl = IOBuffer()
    println(tcl, "mol new {$(abspath(pdbfile))} waitfor all")
    println(tcl, "mol delrep 0 top")
    for rep in reps
        style = get(rep, :style, "Licorice 0.3 12 12")
        color = get(rep, :color, "Name")
        material = get(rep, :material, default_material)
        println(tcl, "mol representation $style")
        println(tcl, "mol selection {$(rep.selection)}")
        println(tcl, "mol color $color")
        println(tcl, "mol material $material")
        println(tcl, "mol addrep top")
    end
    print(tcl, """
    color Display Background white
    color change rgb tan 0.85 0.75 0.55
    display projection Orthographic
    display depthcue off
    display ambientocclusion on
    display aoambient 0.8
    display aodirect 0.4
    display shadows on
    axes location off
    display resize $(size[1]) $(size[2])
    display resetview
    """)
    for (axis, angle) in rotate
        println(tcl, "rotate $axis by $angle")
    end
    println(tcl, "scale by $zoom")
    if !isnothing(clip)
        println(tcl, "display nearclip set $(clip[1])")
        println(tcl, "display farclip set $(clip[2])")
    end
    if !isnothing(cell)
        println(tcl, "package require pbctools")
        println(tcl, "pbc box -color gray -width 2 -center $cell")
    end
    tga = tempname() * ".tga"
    println(tcl, "render TachyonInternal {$tga}")
    println(tcl, "quit")
    script = tempname() * ".tcl"
    write(script, take!(tcl))
    log = read(pipeline(`vmd -dispdev text`; stdin=script, stderr=stdout), String)
    isfile(tga) || error("VMD did not render $imagefile (script: $script):\n$log")
    run(`convert $tga -trim +repage -bordercolor white -border 20 -resize $(width)x\> -quality 85 $imagefile`)
    rm(tga); rm(script)
    return imagefile
end

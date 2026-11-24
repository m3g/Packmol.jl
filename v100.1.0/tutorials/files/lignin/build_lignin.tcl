#
# The lignin oligomers of the tutorial (lignin_*.pdb), built with LigninBuilder
# (https://github.com/jvermaas/LigninBuilder; Vermaas et al., ACS Sustainable Chem. Eng.
# 2019, 7, 3443) from its library of spruce lignin, and minimized with NAMD:
#
#   vmd -dispdev text -e build_lignin.tcl -args <LigninBuilder directory> <namd executable>
#
# 100 oligomers are built (in the directory Spruce); of the 68 of 120 to 400 atoms (sorted
# by size), one in seven is minimized, by LigninBuilder's protocol, and written with its
# bonds (CONECT records), as lignin_L<n>.pdb.
#
set ligninbuilder [lindex $argv 0]
set namd [lindex $argv 1]
lappend auto_path [file join $ligninbuilder LigninBuilderPlugin]
# the CHARMM carbohydrate parameters, distributed with VMD
set env(CHARMMPARDIR) [file join $env(VMDDIR) plugins noarch tcl readcharmmpar1.5]
package require ligninbuilder
package require topotools

::ligninbuilder::buildfromlibrary [file join $ligninbuilder libraries "Spruce Library.txt"] Spruce
::ligninbuilder::makelignincoordinates Spruce .

# the oligomers of 120 to 400 atoms, sorted by size: one in seven
set sizes {}
foreach psf [glob Spruce/L*.psf] {
    set mid [mol new $psf]
    set n [molinfo $mid get numatoms]
    mol delete $mid
    if { $n >= 120 && $n <= 400 } { lappend sizes [list [file rootname [file tail $psf]] $n] }
}
set sizes [lsort -integer -index 1 [lsort -dictionary -index 0 $sizes]]
file mkdir selected
set selected {}
for { set i 3 } { $i < [llength $sizes] } { incr i 7 } {
    set name [lindex $sizes $i 0]
    lappend selected $name
    file copy -force Spruce/$name.psf Spruce/$name.pdb selected
}

# minimization (bonded terms on the CPU: the CUDA builds of NAMD 3 fail on them)
::ligninbuilder::minimizestructures selected $namd "+p8" "bondedCUDA 0"

# the minimized oligomers, with their bonds
foreach name $selected {
    set mid [mol new selected/$name.psf]
    mol addfile selected/$name.pdb waitfor all molid $mid
    set all [atomselect $mid all]
    $all set element [lmap n [$all get name] {string index $n 0}]
    $all writepdb tmp.pdb
    set in [open tmp.pdb r]
    set out [open lignin_$name.pdb w]
    foreach line [split [read -nonewline $in] "\n"] {
        if { [string range $line 0 2] != "END" && [string range $line 0 5] != "CRYST1" } { puts $out $line }
    }
    close $in
    foreach i [$all get index] bonded [$all getbonds] {
        set line [format "CONECT%5d" [expr {$i + 1}]]
        foreach j $bonded { append line [format "%5d" [expr {$j + 1}]] }
        puts $out $line
    }
    puts $out "END"
    close $out
    file delete tmp.pdb
    mol delete $mid
}
quit

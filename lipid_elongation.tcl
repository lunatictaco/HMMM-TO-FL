#######################################################################################
#### This script elongates the truncated HMMM lipids to their full-length (FL) form, ##
#### This script only works for normal phospholipids and includes two steps:       ####
####      1. Manually set the coordinates of the 7th and 8th carbon of each tail   ####
####          to avoid undesired orientation of newly added tails;                 ####
####      2. Elongate the remainging tails by PSFGEN.                              ####
#### The script was originally developed by Dr. Defne Gorgun Ozgulbas and modified ####
#### by Yupeng Li, both from Dr. Emad Tajkhorshid's group.                         ####
#######################################################################################

#######################################################################################
####                              Function Definition                              ####
#######################################################################################
# Procedure to add a standard line for an atom to pdb file
proc writePDBformat {in OutFile} {
    # Extract each column from the input line of the pdb file
    set elem1 [lindex $in 0]
    set elem2 [lindex $in 1]
    set elem3 [lindex $in 2]
    set elem4 [lindex $in 3]
    set elem5 [lindex $in 4]
    set elem6 [lindex $in 5]
    set elem7 [lindex $in 6]
    set elem8 [lindex $in 7]
    set elem9 [lindex $in 8]
    set elem10 [lindex $in 9]
    set elem11 [lindex $in 10]

    # Write formatted line to output file (the number before each element defines the width)
    puts $OutFile [format "%-*4s  %*5d %*4s %-*5s%*4d    %*.3f%*.3f%*.3f%*.2f%*.2f      %-*4s" 4 $elem1 5 $elem2 4 $elem3 5 $elem4 4 $elem5 8 $elem6 8 $elem7 8 $elem8 6 $elem9 6 $elem10 4 $elem11]
}
#######################################################################################
####    Manually set the coordinates of C[23]7 and C[23]8 by editing pdb file      ####
#######################################################################################
set WDIR hmmm2fl_building
cd $WDIR; # Enter the work directory

# Open the input and output pdb files
set InFile [open "./memb.pdb" r]
set OutFile [open "./memb2.pdb" w]

gets $InFile line; # Read the first line from the input pdb file
set atmNum 1; # Initialize the atom number for the output pdb file

# Process each line in the input pdb file
while {[gets $InFile line] >= 0} {
    # Locate to C26 atom
    if {[lindex $line 3] != "CHL1M" && [lindex $line 2] == "C26"} {
        # Keep C26 atom in the output pdb
        set newLine $line
        lset newLine 2 C26
        writePDBformat $newLine $OutFile
        incr atmNum

        # Get the coordinates of C26 atom
        set xValue [lindex $line 5]
        set yValue [lindex $line 6]
        set zValue [lindex $line 7]

        # Add H6R atom
        gets $InFile line
        set newLine $line
        lset newLine 2 H6R
        lset newLine 1 $atmNum
        writePDBformat $newLine $OutFile
        incr atmNum

        # Add H6S atom
        gets $InFile line
        set newLine $line
        lset newLine 2 H6S
        lset newLine 1 $atmNum
        writePDBformat $newLine $OutFile
        incr atmNum
        
        # The H6T atom in the truncated HMMM lipid is ignnored
        # Add C27 and C28 atoms with adjusted z-coordinates referring to the z-position of C26
        gets $InFile line
        # If lipid is in the upper leaflet, put C27 and C28 atoms below C26 atom; otherwise, put them above C26 atom
        set c27zValue [expr {$zValue > 0 ? $zValue - 1.5 : $zValue + 1.5}]
        set c28zValue [expr {$zValue > 0 ? $zValue - 3.0 : $zValue + 3.0}]

        lset newLine 1 $atmNum
        lset newLine 7 $c27zValue
        lset newLine 2 C27
        writePDBformat $newLine $OutFile
        incr atmNum

        lset newLine 1 $atmNum
        lset newLine 7 $c28zValue
        lset newLine 2 C28
        writePDBformat $newLine $OutFile

      # Locate to C36 atom
    } elseif {[lindex $line 3] != "CHL1M" && [lindex $line 2] == "C36"} {
        # Keep C36 atom in the output pdb
        set newLine $line
        lset newLine 2 C36
        writePDBformat $newLine $OutFile
        incr atmNum

        # Get the coordinates of C36 atom
        set xValue [lindex $line 5]
        set yValue [lindex $line 6]
        set zValue [lindex $line 7]

        # Add H6X atom
        gets $InFile line
        set newLine $line
        lset newLine 2 H6X
        lset newLine 1 $atmNum
        writePDBformat $newLine $OutFile
        incr atmNum

        # Add H6Y atom
        gets $InFile line
        set newLine $line
        lset newLine 2 H6Y
        lset newLine 1 $atmNum
        writePDBformat $newLine $OutFile
        incr atmNum

        # The H6Z atom in the truncated HMMM lipid is ignnored
        # Add C37 and C38 atoms with adjusted z-coordinates referring to the z-position of C36
        gets $InFile line
        # If lipid is in the upper leaflet, put C37 and C38 atoms below C36 atom; otherwise, put them above C36 atom
        set c37zValue [expr {$zValue > 0 ? $zValue - 1.5 : $zValue + 1.5}]
        set c38zValue [expr {$zValue > 0 ? $zValue - 3.0 : $zValue + 3.0}]

        lset newLine 1 $atmNum
        lset newLine 7 $c37zValue
        lset newLine 2 C37
        writePDBformat $newLine $OutFile
        incr atmNum

        lset newLine 1 $atmNum
        lset newLine 7 $c38zValue
        lset newLine 2 C38
        writePDBformat $newLine $OutFile

    } else {
        # Write unmodified line to the output pdb file with updated atom number
        set newLine $line
        lset newLine 1 $atmNum
        writePDBformat $newLine $OutFile
    }
    incr atmNum
}

# Close the input and output files
close $OutFile
close $InFile
puts "The addition of C27, C28, C37, and C38 atoms is complete."
#######################################################################################
####                   Elongate the remaining tails by PSFGEN                      ####
#######################################################################################
# Get all the lipid types
mol new memb2.pdb
set all [atomselect top all]
set liptypes [lsort -unique [$all get resname]]

# Find the corresponding topology file for each lipid type
set requiredtopfiles {../toppar/top_all36_prot.rtf ../toppar/top_all36_lipid.rtf ../toppar/top_all36_carb.rtf ../toppar/top_all36_cgenff.rtf ../toppar/toppar_water_ions_namd.str ../toppar/toppar_all36_lipid_cholesterol.str ../toppar/toppar_all36_lipid_inositol.str ../toppar/toppar_all36_lipid_sphingo.str}
set extratopfiles [open "/home/bgworek2/Documents/Fulltail_grow/hmmm2fl_building/extratopfiles.dat" w]; # Save and export the extra topology files for the ring piercing
set alltopfiles [glob ../toppar/*]

foreach liptype $liptypes {
  foreach file $alltopfiles {
    if {[catch {exec grep "RESI $liptype" $file} result]} {
      continue
    } else {
      if {[lsearch $requiredtopfiles $file] == -1} {
        lappend requiredtopfiles $file
        puts $extratopfiles "parameters                          $file"
      }
    }
  }
}
close $extratopfile

mol delete all; # Clear the existing loaded molecule

# Get water segnames
mol new wat.pdb
set watsel [atomselect top water]
set watsegnames [lsort -unique [$watsel get segname]]
mol delete all; # Clear the existing loaded molecule

package require psfgen; # Load PSFGEN plugin
resetpsf; # Delete all segments in the structure

psfgen_logfile "./load_topology.log"; # Log file for loading topology files
# Load each of the required topology files
foreach topfile $requiredtopfiles {
  topology $topfile
}
psfgen_logfile close; # Close the log file for loading topology files

psfgen_logfile "./structure_preparation.log"; #Log the structure preparation process to a new log file

# Define the protein segment
segment P {
  pdb protein.pdb; # Load protein pdb file
  first ACE; # N-terminus cap. May vary depending on the system
  last CTER; # C-terminus cap. May vary depending on the system
}
coordpdb protein.pdb P; # Assign coordinates to the protein segment from the protein pdb file

# Define the water segname
foreach watseg $watsegnames {
  segment $watseg {
    auto none
    pdb $watseg.pdb
  }
  coordpdb $watseg.pdb $watseg
}

# Define the ions segment
segment ION {
  pdb ions.pdb
}
coordpdb ions.pdb ION

# Define the membrane segment
segment MEMB {
  first none
  pdb memb2.pdb; # Load membrane pdb file
  last none
}
coordpdb memb2.pdb MEMB; # Assign coordinates to the membrane segment
guesscoord; # Generate missing coordinates based on bond connectivity in topology files

# Define the membrane segment
segment NMEMB {
  first none
  pdb notmemb.pdb; # Load membrane pdb file
  last none
}
coordpdb notmemb.pdb NMEMB; # Assign coordinates to the membrane segment
guesscoord; # Generate missing coordinates based on bond connectivity in topology files

# Write the strcuture and coordinates of the system containing membran-bound protein and newly extended FL membrane to a psf and a pdb file
writepsf protein_flmemb.psf
writepdb protein_flmemb.pdb

puts "Elongation is complete."
psfgen_logfile close; # Close the log file for structure preparation
cd ..

quit; # Quit VMD

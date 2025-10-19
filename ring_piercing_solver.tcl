#######################################################################################
#### This script is my modification to the minimzation for solving ring piercings  ####
#### The script was originally developed by Dr. Defne Gorgun Ozgulbas and modified ####
#### by Yupeng Li, both from Dr. Emad Tajkhorshid's group.                         ####
####      Then I came in, and messed it up. I mean, make it better. -Bryan         ####
#######################################################################################

#######################################################################################
####                              Function Definition                              ####
#######################################################################################

# Load the longbondeliminator package
lappend auto_path ./LongBondElminiator-main/vmdplugin/
package require longbondeliminator

# Prompt user for a number (with default)
proc prompt_number {prompt default} {
    while 1 {
        puts -nonewline "${prompt} (default: $default): "
        flush stdout
        if {[gets stdin line] < 0 && [eof stdin]} {
            return -code error "EOF detected"
        }
        if {$line eq ""} {
            return $default
        } elseif {[string is double -strict $line]} {
            return $line
        } else {
            puts "Invalid number. Please enter a numeric value."
        }
    }
}


# Default values
set num_cores [prompt_number "Enter number of cores" 16]
set tolerance [prompt_number "Enter bond tolerance" 0.4]

#######################################################################################
####                            Script Execution                                   ####
#######################################################################################

# Create and enter minimization directory
set workdir "minimization"
if {![file exists $workdir]} {
    file mkdir $workdir
}
cd $workdir

# Define source paths (relative to script’s original path)
set src_pdb "../hmmm2fl_building/protein_flmemb.pdb"
set src_psf "../hmmm2fl_building/protein_flmemb.psf"
set dst_pdb "./protein_flmemb.pdb"
set dst_psf "./protein_flmemb.psf"
set src_namd "../minimize.namd"

# Check input files exist
foreach file [list $src_pdb $src_psf] {
    if {![file exists $file]} {
        puts "Error: Required input file '$file' does not exist."
        exit 1
    }
}

# Copy input files
file copy -force $src_pdb $dst_pdb
file copy -force $src_psf $dst_psf
file copy -force $src_namd "minimize.namd"

# Run the minimization
puts "Running minimization with $num_cores cores and tolerance $tolerance..."
::longbondeliminator::minimize ./minimize.namd $num_cores "/usr/local/NAMD_3.0.1_CUDA/namd3 +p$num_cores" $tolerance 0 1

# Get the list of loaded molecules
set mol_list [molinfo list]

if {[llength $mol_list] == 0} {
    puts "No molecule is loaded into VMD. Cannot save output."
} else {
    # Assume last loaded molecule is the one to save
    set molid [lindex $mol_list end]
    puts "Detected loaded molecule (molid $molid). Saving outputs..."
    set sel [atomselect top "all"]
    $sel writepdb ../fulltail_regrow.pdb
    file copy -force protein_flmemb.psf ../fulltail_regrow.psf

    puts "Saved fulltail_regrow.pdb and fulltail_regrow.psf in [pwd]"
}

quit; #Quit VMD



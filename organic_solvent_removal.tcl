#######################################################################################
#### This script removes the organic solvent molecules in the mebrane core.        ####
#### The script was originally developed by Dr. Defne Gorgun Ozgulbas and modified ####
#### by Yupeng Li, both from Dr. Emad Tajkhorshid's group.                         ####
#######################################################################################

#######################################################################################
####                              Function Definition                              ####
#######################################################################################
# Prompt for psf file
proc psfinput {prompt} {
  # Start an infinite loop to repeatedly prompt for valid input
  while 1 {
     # Print the prompt without a newline to keep user input on the same line
    puts -nonewline "${prompt}: "
    flush stdout; # Flush the output buffer to ensure the prompt appears immediately
        
    # Read a line of input from stdin
    if {[gets stdin line] < 0 && [eof stdin]} {
      return -code error "end of file detected"; # If end-of-file (EOF) is detected, return an error with a message
    } elseif {$line eq ""} {
      puts "Input cannot be empty. Please try again."; # If the input is empty, prompt the user to re-enter a valid value
    } else {
      if {[file exists $line]} {
        return $line; # If the input is valid and exists, return the input
      } else {
        puts "Error: File '$line' does not exist. Please enter a valid file path."; # If the file does not exist, display an error and re-prompt the user
      }
    }
  }
}

# Prompt for pdb file
proc pdbinput {prompt} {
  # Start an infinite loop to repeatedly prompt for valid input
  while 1 {
     # Print the prompt without a newline to keep user input on the same line
    puts -nonewline "${prompt}: "
    flush stdout; # Flush the output buffer to ensure the prompt appears immediately

    # Read a line of input from stdin
    if {[gets stdin line] < 0 && [eof stdin]} {
      return -code error "end of file detected"; # If end-of-file (EOF) is detected, return an error with a message
    } elseif {$line eq ""} {
      puts "Input cannot be empty. Please try again."; # If the input is empty, prompt the user to re-enter a valid value
    } else {
      if {[file exists $line]} {
        return $line; # If the input is valid and exists, return the input
      } else {
        puts "Error: File '$line' does not exist. Please enter a valid file path."; # If the file does not exist, display an error and re-prompt the user
      }
    }
  }
}

#######################################################################################
####                            Script Execution                                   ####
#######################################################################################
mkdir hmmm2fl_building; # Create work folder
set WDIR hmmm2fl_building

set psffile [psfinput "Enter your psf file"]; # Enter psf file
set pdbfile [pdbinput "Enther your pdb file"]; # Enter pdb file

mol new $psffile; # Load psf file into VMD
mol addfile $pdbfile; # Load pdb file on top of the psf file

set prot [atomselect top protein]; # Select every atom of protein
set notmemb [atomselect top "resname SAPI13 BSM CHL1"];
set memb [atomselect top "resname POPC POPE POPS"]; # Select atoms of membrane
set ions [atomselect top ions]; # Select ions
set all [atomselect top all]; # Select all the atoms

# Move the z-center of membrane to 0
set zcen [measure center $memb]; # Get z-center of membrane
$all moveby [list 0 0 [expr [lindex $zcen 2] * -1]]

$prot writepdb $WDIR/protein.pdb; # Store protein in a single pdb file
$notmemb writepdb $WDIR/notmemb.pdb; # Store non-truncated membrane components
$memb writepdb $WDIR/memb.pdb; # Store membrane in a single pdb file
$ions writepdb $WDIR/ions.pdb; # Store ions in a single pdb file

# Select and store water molecules by their segnames
set wat [atomselect top water]
$wat writepdb $WDIR/wat.pdb
set watsegnames [lsort -unique [$wat get segname]]
foreach watseg $watsegnames {
  set watsel [atomselect top "segname $watseg"]
  $watsel writepdb $WDIR/$watseg.pdb
  $watsel delete
}

mol delete all; # Clear the existing loaded molecule

quit; # Quit VMD



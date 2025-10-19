package require topotools
#######################################################################################
#### This script resolves potential ring piercing in the system.                   ####
#### The script was originally developed by Dr. Defne Gorgun Ozgulbas and Dr. Josh ####
#### Vermass, modified by Yupeng Li.                                               ####
#### It has been integrated into Ring Piercing Solver plugin in VMD. More info     ####
#### about its installation and usage can be found at the following github page:   ####
#### https://github.com/dgozgulbas/RPplugin?tab=readme-ov-file                     ####
#######################################################################################

#######################################################################################
####                              Function Definition                              ####
#######################################################################################
proc resolve_piercing {namdbin namdargs bondcutoff} {
  # Set path for working directory, structure directory, and parameter file directory
  set WDIR hmmm2fl_ringpiercing
  if {! [file exists $WDIR]} { 
    mkdir $WDIR
  }
  cd $WDIR

  set STRUCTDIR ../hmmm2fl_building

  # Set the input psf and pdb file
  # DO NOT change the names of psf and pdb files here  - they are generated from lipid_elongation.tcl
  set psf $STRUCTDIR/protein_flmemb.psf
  set pdb $STRUCTDIR/protein_flmemb.pdb
    
  # Load psf and pdb file. The molecule id is saved in variable moleid
  set moleid [mol new $psf] 
  mol addfile $pdb

  # The beta value for all atoms in the system is set to 0.
  set asel [atomselect $moleid "all"]
  $asel set beta 0

  # Export the psf/pdb file in which the beta for all the atoms is zero
  animate write psf "system.psf" sel $asel
  animate write pdb "system.pdb" sel $asel

  # A molecular dynamics flexible fitting (MDFF) density map is initialized and saved
  mdffi sim $asel -o "grid_rp.dx" -res 10 -spacing 1

  # The "finished" variable is initialized to 0, representing the state of the completeness of ring piercing. 0 means ring piericng hasn't been resolved yet.
  set finished 0
  # An "iteration" variable is set to 0 for counting 
  set iteration 0
  
  # Select the atoms with abnormal bonds
  set badbeta [atomselect $moleid "beta > 0"]

  ##################################################################################
  ### In each iteration of this loop, a minimization will run first, followed by ###
  ### determining if there's any ring piercing happening in the system.          ###
  ### If there is, the beta value for these atoms will be set to 1, and a new    ###
  ### MDFF grid force will be generated for these atoms, which will be used      ###
  ### in the minimization step of the next iteration.                            ###
  ### This loop continues until no ring piercing exists.                         ###
  ##################################################################################
  # If the variable "finished" is 0, which means ring piercing hasn't been resolved, the loop will continue executing
  while { ! $finished } {
    puts "##################################################################"
    puts "#############           Iteration [expr $iteration + 1] starts           #############"
    # Set NAMD configuration file and log file
    set conffile minimize.conf
    set logfile minimize.log

    # Execute minimization in NAMD
    puts "#############          Minimizing......              #############"
    # ::ExecTool::exec $namdbin $namdargs $conffile > $logfile

    exec $namdbin $namdargs $conffile > $logfile

    animate delete all $moleid; # Delete all frames from memory

    incr iteration; 
    incr finished; # Assume ring piercing is resolved

    # Load the atomic coordinates after minimization  
    mol addfile out.coor type namdbin waitfor all $moleid;

    #set unfinished [vecsum [$asel get beta]]; # Sum up the number of atoms involved in abnormal bonds

    set unfinished [vecsum [$asel get beta]]
    $asel set beta 0; # Reset beta value for all atoms to 0
    $asel set occupancy 0

    # If the bond length is greater than the user-defined cutoff, we will consider it as ring piercing
    # Test all the bonds in the system and find the abnormal bonds
    foreach bond [topo getbondlist -molid $moleid] { 
      if { [measure bond $bond molid $moleid] > $bondcutoff } {
        set ssel [atomselect $moleid "same residue as index $bond"]; # Select atoms with abnormal bonds
        $ssel set beta 1
        $ssel set occupancy 1
        set finished 0; # Variable "finished" is reset to 0 because ring piercing still exists
        $ssel delete
      }
    }

    # Terminate if iterates for more than 100 times
    if { ! $finished && $iteration > 1000 } {
      error "Minimization was not successful, even after 100 iterations."
    }

    # If variable "finished" is 0, which means ring piercing hasn't been resolved,
    # new grid force will be generated for atoms with abnormal bonds
    if { ! $finished } {
      $badbeta update
      mdffi sim $badbeta -o "grid_rp.dx" -res 6 -spacing 0.8
      puts "#############           [$badbeta num] atoms need to be fixed            #############"
    }
    
    if { $unfinished } {
      set finished 0
    }

    # Write final PDB file for the current iteration
    $asel writepdb "system.pdb"
    puts "#############          Iteration $iteration finished          #############"
    puts "##################################################################"
  }

  # Write the final PDB file
  puts "#############         Ring piercing solved         #############"
  $asel writepdb "system.pdb"

  # Clean up selections and loaded molecules
  $asel delete
  $badbeta delete
  mol delete $moleid

  cd ..
  # Copy the final files to the main folder
  puts "#############     Copying the final PSF/PDB     #############"
  exec cp "$WDIR/system.psf" "./PROT_FLMEMB.psf"
  exec cp "$WDIR/system.pdb" "./PROT_FLMEMB.pdb"
}
#######################################################################################
####                            Script Execution                                   ####
#######################################################################################
# Generate .str file which contains essential information about the system
mkdir hmmm2fl_ringpiercing

# Load psf and pdb
mol new hmmm2fl_building/protein_flmemb.psf
mol addfile hmmm2fl_building/protein_flmemb.pdb

# Get the number of lipids in each leaflet
set uppermemb [atomselect top "name P P1 O3 and z > 0"]; # Select lipids in upper leaflet
set lowermemb [atomselect top "name P P1 O3 and z < 0"]; # Select lipids in lower leaflet
set uppernum [$uppermemb num]
set lowernum [$lowermemb num]

# Estimate the system dimension by water
set wat [atomselect top water]; # Select water
set watnum [expr [$wat num] / 3]
set watminmax [measure minmax $wat]
set a [expr [lindex [lindex $watminmax 1] 0] - [lindex [lindex $watminmax 0] 0]]
set b [expr [lindex [lindex $watminmax 1] 1] - [lindex [lindex $watminmax 0] 1]]
set c [expr [lindex [lindex $watminmax 1] 2] - [lindex [lindex $watminmax 0] 2]]

# Get the ion types and the corresponding number of each ion
set ion [atomselect top ions]
set iontypes [lsort -unique [$ion get resname]]
foreach iontype $iontypes {
  set sel [atomselect top "resname $iontype"]
  if {[$sel get charge] > 0} {
    set iontype_pos $iontype
    set ionnum_pos [$sel num]
  } else {
    set iontype_neg $iontype
    set ionnum_neg [$sel num]
  }
  $sel delete
}

# Save the information above to .str file
set strfp [open "hmmm2fl_ringpiercing/system.str" w]

puts $strfp "set boxtype	rect";	# do not change
puts $strfp "set xtltype	CUBIC"; # change to CUBIC in the case of A = B = C, and ORTHorhombic if A != B != C
puts $strfp "set alpha	90.0";	# do not change
puts $strfp "set beta	90.0";	# do not change
puts $strfp "set gamma 90.0"; #do not change
puts $strfp "set zcen 0.0"; # do not change
puts $strfp "set a $a";	# the X dimension
puts $strfp "set b $b";	# the Y dimension
puts $strfp "set c $c"; # the Z dimension
puts $strfp "set nliptop $uppernum";	# the total number of lipids at the top leaflet
puts $strfp "set nlipbot $lowernum";	# the total number of lipids at the bottom leaflet
puts $strfp "set posid $iontype_pos";	# segment ID of positive ions
puts $strfp "set negid $iontype_neg";	# segment ID of negative ions
puts $strfp "set npos $ionnum_pos"; # the number of positive ions
puts $strfp "set nneg $ionnum_neg";	# the number of negative ions
puts $strfp "set nwater $watnum"; # the number of water molecules

close $strfp

# Copy minimize.conf to working folder
exec cp /home/bgworek2/Documents/Scripts/0Scripts/HMMM-TO-FL-SCRIPT/minimize_template.conf hmmm2fl_ringpiercing/minimize.conf

# Resolve ring piercing
resolve_piercing /usr/local/NAMD_3.0.1_CUDA/namd3 "+p30" 2.2

quit; #Quit VMD

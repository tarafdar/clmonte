# (C) 1992-2014 Altera Corporation. All rights reserved.                         
# Your use of Altera Corporation's design tools, logic functions and other       
# software and tools, and its AMPP partner logic functions, and any output       
# files any of the foregoing (including device programming or simulation         
# files), and any associated documentation or information are expressly subject  
# to the terms and conditions of the Altera Program License Subscription         
# Agreement, Altera MegaCore Function License Agreement, or other applicable     
# license agreement, including, without limitation, that your use is for the     
# sole purpose of programming logic devices manufactured by Altera and sold by   
# Altera or its authorized distributors.  Please refer to the applicable         
# agreement for further details.                                                 
    

# This is so that re-compiling within the same directory produces the same
# result from Quartus.
post_message "Deleting incremental_db to ensure imported partition is only used"
file delete -force incremental_db


# Make sure OpenCL SDK installation exists
post_message "Checking for OpenCL SDK installation, environment should have ALTERAOCLSDKROOT defined"
if {[catch {set sdk_root $::env(ALTERAOCLSDKROOT)} result]} {
  post_message -type warning "OpenCL SDK installation not found.  Make sure ALTERAOCLSDKROOT is correctly set"
} else {
  post_message "ALTERAOCLSDKROOT=$::env(ALTERAOCLSDKROOT)"
}

# If for some reason system.qsys wasn't generated, generate it!
if { ! [file exists system/synthesis/system.v] } {
  post_message "Generating system.qsys system since system.v does not exist:"
  post_message "    ip-generate --component-file=system.qsys --file-set=QUARTUS_SYNTH --output-directory=system/synthesis --report-file=qip:system/synthesis/system.qip"
  qexec "ip-generate --component-file=system.qsys --file-set=QUARTUS_SYNTH --output-directory=system/synthesis --report-file=qip:system/synthesis/system.qip"
}

# Replace the clock crosser component (allows for longer paths between clock domains)
file copy -force "ip/clock_crosser/altera_avalon_st_clock_crosser.v" "system/synthesis/submodules/altera_avalon_st_clock_crosser.v"

# copy over verilog files into Qsys generated system.qsys for CvP revision flow fixes
file copy -force scripts/cvpupdatefix/altpcie_hip_256_pipen1b.v system/synthesis/submodules/altpcie_hip_256_pipen1b.v
file copy -force scripts/cvpupdatefix/altpcie_sv_hip_ast_hwtcl.v system/synthesis/submodules/altpcie_sv_hip_ast_hwtcl.v
file copy -force scripts/cvpupdatefix/cvp_update_reset.v system/synthesis/submodules/cvp_update_reset.v
file copy -force scripts/cvpupdatefix/cvp_update_reset_zero.v system/synthesis/submodules/cvp_update_reset_zero.v

# Cvp revision or base revision
if { [string match [lindex $quartus(args) 2] "base" ] == 0 } {
  post_message "Compiling CvP!"
} else {
  post_message "Compiling Base!"
}

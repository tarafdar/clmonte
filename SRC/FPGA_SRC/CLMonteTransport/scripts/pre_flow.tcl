# (C) 1992-2013 Altera Corporation. All rights reserved.                         
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

# Remove top_clean.sdc before Quartus flow
post_message "Removing top_clean.sdc"
file delete -force top_clean.sdc
set fp [open "top_clean.sdc" w]
close $fp

# If for some reason system.qsys wasn't generated, generate it!
if { ! [file exists system/synthesis/system.v] } {
  post_message "Generating system.qsys system since system.v does not exist:"
  post_message "    ip-generate --component-file=system.qsys --file-set=QUARTUS_SYNTH --output-directory=system/synthesis --report-file=qip:system/synthesis/system.qip"
  qexec "ip-generate --component-file=system.qsys --file-set=QUARTUS_SYNTH --output-directory=system/synthesis --report-file=qip:system/synthesis/system.qip"
}

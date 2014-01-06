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
    


package require ::quartus::flow
package require ::quartus::project

post_message "Running post-flow script"

# Dynamically determine the project name by finding the qpf in the directory.
set qpf_files [glob -nocomplain *.qpf]

if {[llength $qpf_files] == 0} {
    error "No QSF detected"
} elseif {[llength $qpf_files] > 1} {
    post_message "Warning: More than one QSF detected. Picking the first one."
}

set qpf_file [lindex $qpf_files 0]
set project_name [string range $qpf_file 0 [expr [string first . $qpf_file] - 1]]

post_message "Project name: $project_name"

project_open $project_name

# Run PLL adjust script
post_message "Running PLL adjustment script"
source iface/ip/bsp/adjust_plls.tcl

# Generate CvP files
post_message "Generating CvP files"
if {[catch {execute_module -tool cpf -args "-c --cvp top.sof top.rbf"} result]} {
  post_message -type error "Error generating CvP files! $result"
  exit 2
}

# Generate fpga.bin used for reprogramming
post_message "Generating fpga.bin"
set argv [list top.sof top.core.rbf]
if {[catch {source scripts/create_fpga_bin.tcl} res]} {
  post_message -type error "Error in create_fpga_bin.tcl! $res"
  exit 2
}

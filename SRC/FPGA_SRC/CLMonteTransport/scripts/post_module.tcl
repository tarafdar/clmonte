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
    

# load required packages
load_package flow

#Get args
set module [lindex $quartus(args) 0]
set top [lindex $quartus(args) 1]

post_message "Running post_module.tcl script"

if [string match "quartus_map" $module] {
  if { ![file exists acl_iface_partition.qxp] }  {
    post_message "Skipping import since acl_iface_partition.qxp not found"
  } else {
    load_package incremental_compilation
    load_package project

    post_message "Running QXP import on project $top"

    # open project
    project_open $top

    # command to execute import
    execute_flow -incremental_compilation_import
  }
}

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
    

#  qsys-script --cmd="set instance instance_name" --script=scripts/remove_kernel.tcl --system-file=system.qsys
package require -exact qsys 12.0

if {[expr ! [info exists instance]]} {
    error "Must use --cmd=\"set instance instance_name \""
}

#set qsys [lindex $args 0]
#set instance [lindex $args 0]

remove_instance $instance

save_system

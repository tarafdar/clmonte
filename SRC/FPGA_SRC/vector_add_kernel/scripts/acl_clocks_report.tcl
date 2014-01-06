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
    

#Run using: quartus_sta -t <script.tcl> <project>
set report_file acl_clocks_report.txt

if { $argc != 1} {error "Error: Usage: quartus_sta -t <script.tcl> <project_name>" }

set proj [lindex $argv 0]

project_open $proj

create_timing_netlist

read_sdc
update_timing_netlist

report_timing -setup -from_clock [get_clocks { kernel_clk }] -to_clock [get_clocks { kernel_clk }] -npaths 10 -detail full_path -panel_name {Kernel 1x Clock Setup} -file $report_file

report_timing -setup -from_clock [get_clocks { kernel_clk2x }] -to_clock [get_clocks { kernel_clk2x }] -npaths 10 -detail full_path -panel_name {Kernel 2x Clock Setup} -file $report_file -append

project_close

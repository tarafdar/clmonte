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
    

#**************************************************************
# This .sbc file is created by Terasic Tool.
# Users are recommended to modify this file to match users logic.
#**************************************************************

#**************************************************************
# Create Clock
#**************************************************************
create_clock -period 100MHz [get_ports config_clk]
create_clock -period 166MHz [get_ports Mem0_RefClk]
create_clock -period 166MHz [get_ports Mem1_RefClk]
create_clock -period 100MHz [get_ports pcie_refclk]

# Override the default 10 MHz JTAG TCK:
create_clock -name altera_reserved_tck -period 30.00  -waveform {0.000 15.0} {altera_reserved_tck}
set_input_delay -clock altera_reserved_tck 8 [get_ports altera_reserved_tdi]
set_input_delay -clock altera_reserved_tck 8 [get_ports altera_reserved_tms]
set_output_delay -clock altera_reserved_tck -clock_fall  -fall -max 10 [get_ports altera_reserved_tdo]
set_output_delay -clock altera_reserved_tck -clock_fall  -rise -max 10 [get_ports altera_reserved_tdo]
set_output_delay -clock altera_reserved_tck -clock_fall  -fall -min .2 [get_ports altera_reserved_tdo]
set_output_delay -clock altera_reserved_tck -clock_fall  -rise -min .2 [get_ports altera_reserved_tdo]

#**************************************************************
# Set Clock Latency
#**************************************************************

#**************************************************************
# Set Input Delay
#**************************************************************

#**************************************************************
# Set Output Delay
#**************************************************************

#**************************************************************
# Set Multicycle Path
#**************************************************************

#**************************************************************
# Set Maximum Delay
#**************************************************************

#**************************************************************
# Set Minimum Delay
#**************************************************************

#**************************************************************
# Set Input Transition
#**************************************************************

#**************************************************************
# Set Load
#**************************************************************

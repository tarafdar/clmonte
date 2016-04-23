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

set_input_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } -max 2.5 [get_ports {UFM_SCL UFM_SDA UCD_SCL UCD_SDA TMP431C_SCL TMP431C_SDA}]

set_input_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_WAIT[*]} ]
set_input_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_DQ[*]} ]

#**************************************************************
# Set Output Delay
#**************************************************************

set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } -add_delay 3 [get_ports {UFM_SCL UFM_SDA UCD_SCL UCD_SDA TMP431C_SCL TMP431C_SDA}]

set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_A[*]} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_DQ[*]} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_CE_n[*]} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_OE_n} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_RST_n} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_WE_n} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_WP_n} ]
set_output_delay -clock { system_inst|board|pcie|altera_s5_a2p|altpcie_hip_256_pipen1b|stratixv_hssi_gen3_pcie_hip|coreclkout } 1  [ get_ports {FLASH_ADV_n} ]

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

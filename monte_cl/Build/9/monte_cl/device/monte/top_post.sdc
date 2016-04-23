#**************************************************************
# Create Generated Clock
#**************************************************************
derive_pll_clocks 

#**************************************************************
# Set Clock Uncertainty
#**************************************************************
derive_clock_uncertainty

#**************************************************************
# Set Clock Groups
#**************************************************************
set_clock_groups -asynchronous \
-group { \
   config_clk \
} -group { \
   Mem0_RefClk \
} -group { \
   Mem1_RefClk \
} -group { \
   pcie_refclk \
   [get_clocks {system_inst|board|pcie|*}] \
} -group { \
   [get_clocks {system_inst|board|acl_kernel_clk|kernel_pll|*}] \
} -group { \
   altera_reserved_tck \
} -group { \
    Mem0_Ck \
    Mem0_Ck_n \
    Mem0_Dqs[0]_IN \
    Mem0_Dqs[0]_OUT \
    Mem0_Dqs[1]_IN \
    Mem0_Dqs[1]_OUT \
    Mem0_Dqs[2]_IN \
    Mem0_Dqs[2]_OUT \
    Mem0_Dqs[3]_IN \
    Mem0_Dqs[3]_OUT \
    Mem0_Dqs[4]_IN \
    Mem0_Dqs[4]_OUT \
    Mem0_Dqs[5]_IN \
    Mem0_Dqs[5]_OUT \
    Mem0_Dqs[6]_IN \
    Mem0_Dqs[6]_OUT \
    Mem0_Dqs[7]_IN \
    Mem0_Dqs[7]_OUT \
    Mem0_Dqs_n[0]_OUT \
    Mem0_Dqs_n[1]_OUT \
    Mem0_Dqs_n[2]_OUT \
    Mem0_Dqs_n[3]_OUT \
    Mem0_Dqs_n[4]_OUT \
    Mem0_Dqs_n[5]_OUT \
    Mem0_Dqs_n[6]_OUT \
    Mem0_Dqs_n[7]_OUT \
    [get_clocks {system_inst|board|ddr3a|*}] \
} -group { \
    Mem1_Ck \
    Mem1_Ck_n \
    Mem1_Dqs[0]_IN \
    Mem1_Dqs[0]_OUT \
    Mem1_Dqs[1]_IN \
    Mem1_Dqs[1]_OUT \
    Mem1_Dqs[2]_IN \
    Mem1_Dqs[2]_OUT \
    Mem1_Dqs[3]_IN \
    Mem1_Dqs[3]_OUT \
    Mem1_Dqs[4]_IN \
    Mem1_Dqs[4]_OUT \
    Mem1_Dqs[5]_IN \
    Mem1_Dqs[5]_OUT \
    Mem1_Dqs[6]_IN \
    Mem1_Dqs[6]_OUT \
    Mem1_Dqs[7]_IN \
    Mem1_Dqs[7]_OUT \
    Mem1_Dqs_n[0]_OUT \
    Mem1_Dqs_n[1]_OUT \
    Mem1_Dqs_n[2]_OUT \
    Mem1_Dqs_n[3]_OUT \
    Mem1_Dqs_n[4]_OUT \
    Mem1_Dqs_n[5]_OUT \
    Mem1_Dqs_n[6]_OUT \
    Mem1_Dqs_n[7]_OUT \
    [get_clocks {system_inst|board|ddr3b|*}] \
}

#**************************************************************
# Set False Path
#**************************************************************

# Cut path to uniphy reset, this reset is asyncronous
set_false_path -from system:system_inst|system_board:board|altera_reset_controller:reset_controller_global|altera_reset_synchronizer:alt_rst_sync_uq1|altera_reset_synchronizer_int_chain_out -to [get_registers *umemphy*ureset*reset_reg*]

# Cut path to pcie npor - this signal is asynchronous
set_false_path -from system:system_inst|system_board:board|sw_reset:por_reset_counter|sw_reset_n_out -to system:system_inst|system_board:board|altpcie_sv_hip_avmm_hwtcl:pcie*

# Make the kernel reset multicycle
set_multicycle_path -to * -setup 2 -from {system:system_inst|system_board:board|system_board_kernel_interface:kernel_interface|altera_reset_controller:reset_controller_sw|altera_reset_synchronizer:alt_rst_sync_uq1|altera_reset_synchronizer_int_chain_out}
set_multicycle_path -to * -hold 1 -from {system:system_inst|system_board:board|system_board_kernel_interface:kernel_interface|altera_reset_controller:reset_controller_sw|altera_reset_synchronizer:alt_rst_sync_uq1|altera_reset_synchronizer_int_chain_out}

# FLASH Interface
set_false_path -to [get_ports {FLASH_RST_n}]

# Optionally overconstrain the kernel clock in the iface or relax the kernel
#if { $::TimeQuestInfo(nameofexecutable) == "quartus_fit" } {
#set_max_delay -to [get_clocks system_inst|board|acl_kernel_clk|kernel_pll|altera_pll_i|stratixv_pll|counter[0].output_counter|divclk] -through system_inst|board* 2.5
#set_max_delay -to [get_clocks system_inst|board|acl_kernel_clk|kernel_pll|altera_pll_i|stratixv_pll|counter[1].output_counter|divclk] -through system_inst|board* 1.25
#set non_iface_keepers [remove_from_collection [get_keepers *] [get_keepers system_inst\|board\|*]]
#set_max_delay 3.5 -from $non_iface_keepers -to $non_iface_keepers
#}


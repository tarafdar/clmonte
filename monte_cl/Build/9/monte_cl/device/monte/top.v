//------------------------------------------------------------------------------
//      Nallatech is providing this design, code, or information "as is".       
//      solely for use on Nallatech systems and equipment.                     
//      By providing this design, code, or information                         
//      as one possible implementation of this feature, application            
//      or standard, NALLATECH IS MAKING NO REPRESENTATION THAT THIS           
//      IMPLEMENTATION IS FREE FROM ANY CLAIMS OF INFRINGEMENT,                
//      AND YOU ARE RESPONSIBLE FOR OBTAINING ANY RIGHTS YOU MAY REQUIRE       
//      FOR YOUR IMPLEMENTATION.  NALLATECH EXPRESSLY DISCLAIMS ANY            
//      WARRANTY WHATSOEVER WITH RESPECT TO THE ADEQUACY OF THE                
//      IMPLEMENTATION, INCLUDING BUT NOT LIMITED TO ANY WARRANTIES OR         
//      REPRESENTATIONS THAT THIS IMPLEMENTATION IS FREE FROM CLAIMS OF        
//      INFRINGEMENT, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS        
//      FOR A PARTICULAR PURPOSE.                                              
//                                                                             
//      USE OF SOFTWARE. This software contains elements of software code      
//      which are the property of Nallatech Limited (Nallatech Software).      
//      Use of the Nallatech Software by you is permitted only if you hold a   
//      valid license from Nallatech Limited or a valid sub-license from a     
//      licensee of Nallatech Limited. Use of such software shall be governed  
//      by the terms of such license or sub-license agreement.                 
//      The Nallatech software is for use solely on Nallatech hardware         
//      unless you hold a license permitting use on other hardware.            
//                                                                             
//      This Nallatech Software is protected by copyright law and              
//      international treaties. Unauthorized reproduction or distribution of   
//      this software, or any portion of it, may result in severe civil and    
//      criminal penalties, and will be prosecuted to the maximum extent       
//      possible under law. Nallatech products are covered by one or more      
//      patents. Other US and international patents pending.                   
//      Please see www.nallatech.com for more information                      
//                                                                             
//      Nallatech products are not intended for use in life support            
//      appliances, devices, or systems. Use in such applications is           
//      expressly prohibited.                                                  
//                                                                             
//      Copyright © 1998-2014 Nallatech Limited. All rights reserved.          
//------------------------------------------------------------------------------
// $Id$
//------------------------------------------------------------------------------
//
//                         N
//                        NNN
//                       NNNNN
//                      NNNNNNN
//                     NNNN-NNNN          Nallatech
//                    NNNN---NNNN         (a subsidiary of
//                   NNNN-----NNNN        Interconnect Systems Inc.)
//                  NNNN-------NNNN
//                 NNNN---------NNNN
//                NNNNNNNN---NNNNNNNN
//               NNNNNNNNN---NNNNNNNNN
//                -------------------
//               ---------------------
//
//------------------------------------------------------------------------------
// Title       : Top Level HDL of OpenCL BSP
// Project     : PCIe385 HPC
//------------------------------------------------------------------------------
// Description : 
//
//
//------------------------------------------------------------------------------
// Known Issues and Omissions:
//
//
//------------------------------------------------------------------------------

module top(

   //////// CLOCK //////////
   input config_clk,  // 100MHz clock 
   input Mem0_RefClk, // 166MHz clock
   input Mem1_RefClk, // 166MHz clock

   //////// Power //////////
   output USER_1P35V_EN,

   //////// PCIe //////////
   input        pcie_refclk,
   input        perstl0_n,  // Reset to embedded PCIe
   input        hip_serial_rx_in0,
   input        hip_serial_rx_in1,
   input        hip_serial_rx_in2,
   input        hip_serial_rx_in3,
   input        hip_serial_rx_in4,
   input        hip_serial_rx_in5,
   input        hip_serial_rx_in6,
   input        hip_serial_rx_in7,
   output       hip_serial_tx_out0,
   output       hip_serial_tx_out1,
   output       hip_serial_tx_out2,
   output       hip_serial_tx_out3,
   output       hip_serial_tx_out4,
   output       hip_serial_tx_out5,
   output       hip_serial_tx_out6,
   output       hip_serial_tx_out7,

   //////// DDR3 //////////
   output        Mem0_Reset_n,
   output [15:0] Mem0_Addr,
   output [2:0]  Mem0_Bank,
   output        Mem0_Cas_n,
   output        Mem0_Cke,
   output        Mem0_Ck,
   output        Mem0_Ck_n,
   output        Mem0_Cs_n,
   output [8:0]  Mem0_Dm,
   inout [71:0]  Mem0_Dq,
   inout  [8:0]  Mem0_Dqs,
   inout  [8:0]  Mem0_Dqs_n,
   output        Mem0_Odt,
   output        Mem0_Ras_n,
   output        Mem0_We_n,
   input         Mem0_Rzq,

   //////// DDR3 //////////
   output        Mem1_Reset_n,
   output [15:0] Mem1_Addr,
   output [2:0]  Mem1_Bank,
   output        Mem1_Cas_n,
   output        Mem1_Cke,
   output        Mem1_Ck,
   output        Mem1_Ck_n,
   output        Mem1_Cs_n,
   output [8:0]  Mem1_Dm,
   inout [71:0]  Mem1_Dq,
   inout  [8:0]  Mem1_Dqs,
   inout  [8:0]  Mem1_Dqs_n,
   output        Mem1_Odt,
   output        Mem1_Ras_n,
   output        Mem1_We_n,
   input         Mem1_Rzq,

   //////// 10G Ethernet //////////
   output        TenGbE_RefClkSel,
   output [1:0]  TenGbE_TxDisable,

   //////// UFM ////////
   inout         UFM_SCL,
   inout         UFM_SDA,

   //////// UCD ////////
   inout         UCD_SCL,
   inout         UCD_SDA,
   
   //////// TMP431C ////////
   inout         TMP431C_SCL,
   inout         TMP431C_SDA,

   //////// FLASH Interface ////////
   output [26:2] FLASH_A,
   inout  [31:0] FLASH_DQ,
   output        FLASH_OE_n,
   output        FLASH_WE_n,
   output [1:0]  FLASH_CE_n,

   //////// FLASH Options ////////
   input [1:0]   FLASH_WAIT,
   output        FLASH_WP_n,
   output        FLASH_ADV_n,
   output        FLASH_CLK,
   output        FLASH_RST_n,

   //////// PFL Interface ////////
   output        PFL_FLASH_REQ,
   input         PFL_FLASH_GRNT,
   output        FPGA_SBOOT_EN,
   output        FPGA_SBOOT_PAGE,

   //////// LED //////////
   output [7:0]  leds
);

//=======================================================
//  PARAMETER declarations
//=======================================================

//=======================================================
//  REG/WIRE declarations
//=======================================================
wire          resetn;
wire          npor;
wire [26:0]   flash0_a;
wire [31:0]   flash_data;
wire          flash0_oe_n;
wire          flash1_oe_n;
wire          flash0_we_n;
wire          flash1_we_n;
wire [1:0]    flash_ce;
wire [7:0]    host_control;

//=======================================================
//  Board-specific 
//=======================================================

assign USER_1P35V_EN = 1'b1;       // Force 1.35V supply on

assign TenGbE_RefClkSel = 1'b1;    // Select 644.53125 MHz clock input
assign TenGbE_TxDisable = 2'b11;   // TX disable for SFP

assign FLASH_CE_n[0]   = ( (PFL_FLASH_GRNT==1'b1) ? flash_ce[0] : 1'bZ);
assign FLASH_CE_n[1]   = ( (PFL_FLASH_GRNT==1'b1) ? flash_ce[0] : 1'bZ);
assign FLASH_A[26:2]   = ( (PFL_FLASH_GRNT==1'b1) ? {flash0_a[26:2]} : 25'bZ);
assign FLASH_OE_n      = ( (PFL_FLASH_GRNT==1'b1) ? flash0_oe_n : 1'bZ);
assign FLASH_WE_n      = ( (PFL_FLASH_GRNT==1'b1) ? flash0_we_n : 1'bZ);
assign PFL_FLASH_REQ   = ~host_control[0];
assign FPGA_SBOOT_EN   = host_control[1];
assign FPGA_SBOOT_PAGE = host_control[2];
assign FLASH_WP_n      = ( (PFL_FLASH_GRNT==1'b1) ? 1'b1 : 1'bZ);
assign FLASH_ADV_n     = ( (PFL_FLASH_GRNT==1'b1) ? 1'b0 : 1'bZ);
assign FLASH_CLK       = ( (PFL_FLASH_GRNT==1'b1) ? 1'b0 : 1'bZ);
assign FLASH_RST_n     = ( (PFL_FLASH_GRNT==1'b1) ? 1'b1 : 1'bZ);

//=======================================================
//  Reset logic 
//=======================================================
assign resetn = perstl0_n;

//=======================================================
//  System instantiation
//=======================================================

system system_inst 
(
   // Global signals
   .global_reset_reset_n( resetn ),  // No hard reset !!!
   .config_clk_clk( config_clk ),
   .ddr3a_pll_ref_clk( Mem0_RefClk ),
   .ddr3b_pll_ref_clk( Mem1_RefClk ),

   // PCIe pins
   .pcie_npor_pin_perst(perstl0_n),
   .pcie_npor_npor(npor),
   .pcie_npor_out_reset_n(npor),
   .pcie_refclk_clk( pcie_refclk ),
   .pcie_rx_in0( hip_serial_rx_in0 ),
   .pcie_rx_in1( hip_serial_rx_in1 ),
   .pcie_rx_in2( hip_serial_rx_in2 ),
   .pcie_rx_in3( hip_serial_rx_in3 ),
   .pcie_rx_in4( hip_serial_rx_in4 ),
   .pcie_rx_in5( hip_serial_rx_in5 ),
   .pcie_rx_in6( hip_serial_rx_in6 ),
   .pcie_rx_in7( hip_serial_rx_in7 ),
   .pcie_tx_out0( hip_serial_tx_out0 ),
   .pcie_tx_out1( hip_serial_tx_out1 ),
   .pcie_tx_out2( hip_serial_tx_out2 ),
   .pcie_tx_out3( hip_serial_tx_out3 ),
   .pcie_tx_out4( hip_serial_tx_out4 ),
   .pcie_tx_out5( hip_serial_tx_out5 ),
   .pcie_tx_out6( hip_serial_tx_out6 ),
   .pcie_tx_out7( hip_serial_tx_out7 ),

   // DDR3 pins
   .ddr3a_mem_reset_n( Mem0_Reset_n ),
   .ddr3a_mem_a( Mem0_Addr ),
   .ddr3a_mem_ba( Mem0_Bank ),
   .ddr3a_mem_cas_n( Mem0_Cas_n ),
   .ddr3a_mem_ck( Mem0_Ck ),
   .ddr3a_mem_ck_n( Mem0_Ck_n ),
   .ddr3a_mem_cke( Mem0_Cke ),
   .ddr3a_mem_cs_n( Mem0_Cs_n ),
   .ddr3a_mem_dm( Mem0_Dm ),
   .ddr3a_mem_dq( Mem0_Dq ),
   .ddr3a_mem_dqs_n( Mem0_Dqs_n ),
   .ddr3a_mem_dqs( Mem0_Dqs ),
   .ddr3a_mem_oct_rzqin( Mem0_Rzq ),
   .ddr3a_mem_odt( Mem0_Odt ),
   .ddr3a_mem_ras_n( Mem0_Ras_n ),
   .ddr3a_mem_we_n( Mem0_We_n ),

   // DDR3 pins
   .ddr3b_mem_reset_n( Mem1_Reset_n ),
   .ddr3b_mem_a( Mem1_Addr ),
   .ddr3b_mem_ba( Mem1_Bank ),
   .ddr3b_mem_cas_n( Mem1_Cas_n ),
   .ddr3b_mem_ck( Mem1_Ck ),
   .ddr3b_mem_ck_n( Mem1_Ck_n ),
   .ddr3b_mem_cke( Mem1_Cke ),
   .ddr3b_mem_cs_n( Mem1_Cs_n ),
   .ddr3b_mem_dm( Mem1_Dm ),
   .ddr3b_mem_dq( Mem1_Dq ),
   .ddr3b_mem_dqs_n( Mem1_Dqs_n ),
   .ddr3b_mem_dqs( Mem1_Dqs ),
   .ddr3b_mem_oct_rzqin( Mem1_Rzq ),
   .ddr3b_mem_odt( Mem1_Odt ),
   .ddr3b_mem_ras_n( Mem1_Ras_n ),
   .ddr3b_mem_we_n( Mem1_We_n ),

   .reconfig_to_xcvr_reconfig_to_xcvr({10 {24'h0, 2'b11, 44'h0}}),
   .reconfig_from_xcvr_reconfig_from_xcvr(),

   .kernel_pll_refclk_clk( Mem1_RefClk ),

   // I2C Interfaces
   .i2c_ucd_scl_pad_io( UCD_SCL ),
   .i2c_ucd_sda_pad_io( UCD_SDA ),
   .i2c_ufm_scl_pad_io( UFM_SCL ),
   .i2c_ufm_sda_pad_io( UFM_SDA ),
   .i2c_tmp431c_scl_pad_io( TMP431C_SCL ),
   .i2c_tmp431c_sda_pad_io( TMP431C_SDA ),

   // Flash Interface
   .pfl_flash_req_export( host_control ),
   .pfl_flash_grnt_export( PFL_FLASH_GRNT ),
   .flash_if_tcm_address_out( flash0_a ),
   .flash_if_tcm_read_n_out( flash0_oe_n ),
   .flash_if_tcm_write_n_out( flash0_we_n ),
   .flash_if_tcm_data_out( FLASH_DQ[31:0] ),
   .flash_if_tcm_chipselect_n_out( flash_ce[0] )
);

assign leds[7:0] = 8'b0101000;

endmodule

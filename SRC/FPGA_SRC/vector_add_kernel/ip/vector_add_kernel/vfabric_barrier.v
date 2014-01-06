// (C) 1992-2013 Altera Corporation. All rights reserved.                         
// Your use of Altera Corporation's design tools, logic functions and other       
// software and tools, and its AMPP partner logic functions, and any output       
// files any of the foregoing (including device programming or simulation         
// files), and any associated documentation or information are expressly subject  
// to the terms and conditions of the Altera Program License Subscription         
// Agreement, Altera MegaCore Function License Agreement, or other applicable     
// license agreement, including, without limitation, that your use is for the     
// sole purpose of programming logic devices manufactured by Altera and sold by   
// Altera or its authorized distributors.  Please refer to the applicable         
// agreement for further details.                                                 
    


module vfabric_barrier(clock, resetn, i_dataa, i_dataa_valid, o_dataa_stall, 
    i_datab, i_datab_valid, o_datab_stall,
    i_datac, i_datac_valid, o_datac_stall, i_datad, i_datad_valid, o_datad_stall,
    o_dataout, o_dataout_valid, i_dataout_stall, i_start, i_counter_value, i_settings);

parameter DATA_WIDTH = 1;
parameter CONFIG_WIDTH = 20;
parameter FIFO_DEPTH = 64;
parameter BARRIER_MAX_SIMULTANEOUS_WORKGROUPS = 4;
parameter BARRIER_WORKGROUP_SIZE_BITS = 10;
  
  input clock, resetn;
  input [DATA_WIDTH-1:0] i_dataa;
  input [DATA_WIDTH-1:0] i_datab;
  input [DATA_WIDTH-1:0] i_datac;
  input [DATA_WIDTH-1:0] i_datad;
  input i_dataa_valid, i_datab_valid, i_datac_valid, i_datad_valid;
  output o_dataa_stall, o_datab_stall, o_datac_stall, o_datad_stall;
  output [DATA_WIDTH-1:0] o_dataout;
  output o_dataout_valid;
  input i_dataout_stall;
  
  input i_start;
  input [31:0] i_counter_value;
  input [CONFIG_WIDTH-1:0] i_settings;

  wire data_to_count;
  wire data_to_count_valid;
  
  vfabric_pla pla0(.clock(clock), .resetn(resetn), .i_dataa(i_dataa), .i_dataa_valid(i_dataa_valid), 
      .o_dataa_stall(o_dataa_stall), .i_datab(i_datab), .i_datab_valid(i_datab_valid), 
      .o_datab_stall(o_datab_stall), .i_datac(i_datac), .i_datac_valid(i_datac_valid),
      .o_datac_stall(o_datac_stall), .i_datad(i_datad), .i_datad_valid(i_datad_valid), 
      .o_datad_stall(o_datad_stall), .o_dataout(data_to_count), 
      .o_dataout_valid(data_to_count_valid), .i_stall(1'b0), .i_settings(i_settings),
      .i_start(i_start));
  defparam pla0.CONFIG_WIDTH = CONFIG_WIDTH;
  defparam pla0.FIFO_DEPTH = FIFO_DEPTH;
    
  acl_barrier_simple barrier0(
  	.clock(clock),
  	.resetn(resetn),
  	.valid_entry(data_to_count_valid),
  	.data_entry(),
  	.stall_entry(1'b0),
  	.valid_exit(o_dataout_valid),
  	.data_exit(),
  	.stall_exit(i_dataout_stall),
  	.num_live_workitems(),
    .workgroup_size(i_counter_value));
  defparam barrier0.DATA_WIDTH = DATA_WIDTH;
  defparam barrier0.MAX_SIMULTANEOUS_WORKGROUPS = BARRIER_MAX_SIMULTANEOUS_WORKGROUPS;
  defparam barrier0.WORKGROUP_SIZE_BITS = BARRIER_WORKGROUP_SIZE_BITS;
 
endmodule

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
    


module acl_embedded_workgroup_issuer_fifo #(
  parameter integer MAX_SIMULTANEOUS_WORKGROUPS = 2,    // >0
  parameter integer WORKGROUP_SIZE_ARG_WIDTH = 32       // >0
)
(
  input logic clock, 
  input logic resetn, 

  // Handshake for item entry into function.
  input logic valid_in, 
  output logic stall_out, 
  
  // Handshake with entry basic block
  output logic valid_entry, 
  input logic stall_entry,

  // Observe threads exiting the function .
  // This is asserted when items are ready to be retired from the workgroup.
  input logic valid_exit, 
  // This is asserted when downstream is not ready to retire an item from
  // the workgroup.
  input logic stall_exit, 
  
  // Observe which work item is exiting the function.
  // Really only care about hw work group id, which is stored
  // in the most significant bits.
  input logic [31:0] done_local_id_3, 
  
  // Need workgroup_size to know when one workgroup ends
  // and another begins.
  input logic [WORKGROUP_SIZE_ARG_WIDTH-1:0] workgroup_size, 
  
  // Tell the downstream block the what hw group the issuing
  // The most significant bits of local_id_3 encode the hw group
  // to which the issuing item belongs.
  output logic [31:0] local_id_3, 
  
  // Tell downstream blocks how many work items are still alive
  // in this workgroup.
  output logic [WORKGROUP_SIZE_ARG_WIDTH*MAX_SIMULTANEOUS_WORKGROUPS-1:0] live_thread_count,

  // Pass though global_id, local_id and group_id.
  input logic [31:0] global_id_in[2:0],
  input logic [31:0] local_id_in[2:0],
  input logic [31:0] group_id_in[2:0],
  output logic [31:0] global_id_out[2:0],
  output logic [31:0] local_id_out[2:0],
  output logic [31:0] group_id_out[2:0]
);
  localparam integer MAX_WG_SIZE = {WORKGROUP_SIZE_ARG_WIDTH{1'b1}};
  localparam integer WG_LIMIT_BITS = MAX_SIMULTANEOUS_WORKGROUPS > 1 ? $clog2(MAX_SIMULTANEOUS_WORKGROUPS) : 1;
  localparam integer SCALAR_LOCAL_ID_W = (WORKGROUP_SIZE_ARG_WIDTH > (32-$clog2(MAX_SIMULTANEOUS_WORKGROUPS))) ? 32-$clog2(MAX_SIMULTANEOUS_WORKGROUPS) : WORKGROUP_SIZE_ARG_WIDTH;

  logic [WG_LIMIT_BITS-1:0] entry_l_wgid;

  // Entry: 1 cycle latency
  // Exit: 1 cycle latency
  acl_work_group_limiter #(
    .WG_LIMIT(MAX_SIMULTANEOUS_WORKGROUPS),
    .KERNEL_WG_LIMIT(MAX_SIMULTANEOUS_WORKGROUPS),
    .MAX_WG_SIZE(MAX_WG_SIZE),
    .WG_FIFO_ORDER(1),
    .IMPL("kernel")   // this parameter is very important to get the right implementation
  )
  limiter(
    .clock(clock),
    .resetn(resetn),
    .wg_size(workgroup_size),

    .entry_valid_in(valid_in),
    .entry_k_wgid(),
    .entry_stall_out(stall_out),
    .entry_valid_out(valid_entry),
    .entry_l_wgid(entry_l_wgid),
    .entry_stall_in(stall_entry),

    .exit_valid_in(valid_exit & ~stall_exit),
    .exit_l_wgid(done_local_id_3[31:32-WG_LIMIT_BITS]),
    .exit_stall_out(),
    .exit_valid_out(),
    .exit_stall_in(1'b0)
  );

  // Pass through ids (global, local, group).
  // Match the latency of the work-group limiter, which is one cycle.
  always @(posedge clock)
    if( ~stall_entry ) 
    begin
      global_id_out <= global_id_in;
      local_id_out <= local_id_in;
      group_id_out <= group_id_in;
    end

  // scalarized local id generator
  logic [SCALAR_LOCAL_ID_W-1:0] scalar_local_id;
  always @(posedge clock or negedge resetn)
    if( ~resetn )
      scalar_local_id <= '0;
    else if( valid_entry & ~stall_entry )
    begin
      if( scalar_local_id == workgroup_size - 'd1 )
        scalar_local_id <= '0;
      else
        scalar_local_id <= scalar_local_id + 'd1;
    end

  // complete local_id_3
  assign local_id_3[SCALAR_LOCAL_ID_W-1:0] = scalar_local_id;
  generate
    if( MAX_SIMULTANEOUS_WORKGROUPS > 1 )
      assign local_id_3[31:32-WG_LIMIT_BITS] = entry_l_wgid;
  endgenerate

endmodule
// vim:set filetype=verilog_systemverilog:


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
    


module acl_embedded_workgroup_issuer #(
  parameter integer MAX_SIMULTANEOUS_WORKGROUPS = 2,    // >0
  parameter integer WORKGROUP_SIZE_ARG_WIDTH = 32,      // >0
  parameter string WORKGROUP_EXIT_ORDER = "fifo"     // fifo|noninterleaved|unknown
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
generate
if( WORKGROUP_EXIT_ORDER == "fifo" )
begin
  acl_embedded_workgroup_issuer_fifo #(
    .MAX_SIMULTANEOUS_WORKGROUPS(MAX_SIMULTANEOUS_WORKGROUPS),
    .WORKGROUP_SIZE_ARG_WIDTH(WORKGROUP_SIZE_ARG_WIDTH)
  )
  issuer(.*);
end
else if( WORKGROUP_EXIT_ORDER == "noninterleaved" || WORKGROUP_EXIT_ORDER == "unknown" )
begin
  acl_embedded_workgroup_issuer_complex #(
    .MAX_SIMULTANEOUS_WORKGROUPS(MAX_SIMULTANEOUS_WORKGROUPS),
    .WORKGROUP_SIZE_ARG_WIDTH(WORKGROUP_SIZE_ARG_WIDTH)
  )
  issuer(.*);
end
else
begin
    // synthesis translate off
    initial
      $fatal("%m: unsupported configuration (WORKGROUP_EXIT_ORDER=%s)", WORKGROUP_EXIT_ORDER);
    // synthesis translate on
end
endgenerate

endmodule


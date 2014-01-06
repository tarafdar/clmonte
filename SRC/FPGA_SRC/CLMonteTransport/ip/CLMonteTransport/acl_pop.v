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
    


//===----------------------------------------------------------------------===//
//
// C backend 'pop' primitive
//
//===----------------------------------------------------------------------===//
module acl_pop (
	clock,
	resetn,

	// input stream from kernel pipeline
	dir,
    valid_in,
	data_in,
    stall_out,

	// downstream, to kernel pipeline
    valid_out,
    stall_in,
    data_out,

	// feedback downstream, from feedback acl_push
    feedback_in,
    feedback_valid_in,
    feedback_stall_out
);

    parameter DATA_WIDTH = 32;

input clock, resetn, stall_in, valid_in, feedback_valid_in;
output stall_out, valid_out, feedback_stall_out;
input [DATA_WIDTH-1:0] data_in;
input dir;
output [DATA_WIDTH-1:0] data_out;
input [DATA_WIDTH-1:0] feedback_in;

wire feedback_downstream, data_downstream;

assign feedback_downstream = valid_in & ~dir & feedback_valid_in;
assign data_downstream = valid_in & dir;

assign valid_out = feedback_downstream | data_downstream;
assign data_out = ~dir ? feedback_in : data_in;

//assign stall_out = stall_in;
//assign stall_out = valid_in & ~((feedback_downstream | data_downstream) & ~stall_in);
// assign stall_out = ~((feedback_downstream | data_downstream) & ~stall_in);
// stall upstream if
//   downstream is stalling (stall_in)
//   I'm waiting for data from feedback (valid_in&~dir&~feedback_valiid_in)
 assign stall_out = ( valid_in & ~dir & ~feedback_valid_in ) | stall_in;

// don't accept data if:
//  downstream cannot accept data (stall_in)
//  data from upstream is selected (data_downstream)
//  no thread exists to read data (~valid_in)
assign feedback_stall_out = stall_in  | data_downstream | ~valid_in; 

endmodule


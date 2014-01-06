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
// C backend 'push' primitive
//
// Upstream are signals that go to the feedback (snk node is a acl_pop),
// downstream are signals that continue into our "normal" pipeline.
//
// dir indicates if you want to push it to the feedback
//   1 - push to feedback
//   0 - bypass, just push out to downstream
//===----------------------------------------------------------------------===//
module acl_push (
	clock,
	resetn,

	// interface from kernel pipeline, input stream
	dir,
	data_in,
    valid_in,
    stall_out,

	// interface to kernel pipeline, downstream
    valid_out,
    stall_in,
    data_out,

	// interface to pipeline feedback, upstream
    feedback_out,
    feedback_valid_out,
    feedback_stall_in
);

    parameter DATA_WIDTH = 32;
    parameter FIFO_DEPTH = 1;
    parameter MIN_FIFO_LATENCY = 0;
    // style can be "REGULAR", for a regular push
    // or "TOKEN" for a special fifo that hands out tokens 
    parameter string STYLE = "REGULAR";     // "REGULAR"/"TOKEN"

input clock, resetn, stall_in, valid_in, feedback_stall_in;
output stall_out, valid_out, feedback_valid_out;
input [DATA_WIDTH-1:0] data_in;
input dir;
output [DATA_WIDTH-1:0] data_out, feedback_out;

wire [DATA_WIDTH-1:0] feedback;
wire data_downstream, data_upstream;

assign data_upstream = valid_in & dir;
assign data_downstream = valid_in;

wire feedback_stall, feedback_valid;

reg consumed_downstream, consumed_upstream;

assign valid_out = data_downstream & !consumed_downstream;
assign feedback_valid = data_upstream & !consumed_upstream;
assign data_out = data_in;
assign feedback = data_in;

//assign stall_out = valid_in & ( ~(data_downstream & ~stall_in) & ~(data_upstream & ~feedback_stall));
// assign stall_out = valid_in & ( ~(data_downstream & ~stall_in) | ~(data_upstream & ~feedback_stall));
assign stall_out = stall_in | (feedback_stall & dir & valid_in) ; //valid_in & ( ~(data_downstream & ~stall_in) | ~(data_upstream & ~feedback_stall));

always @(posedge clock or negedge resetn) begin
   if (!resetn) begin
      consumed_downstream <= 1'b0;
      consumed_upstream <= 1'b0;
   end else begin
      if (consumed_downstream)
        consumed_downstream <= stall_out;
      else  
        consumed_downstream <= stall_out & (data_downstream & ~stall_in);

      if (consumed_upstream)
        consumed_upstream <= stall_out;
      else  
        consumed_upstream <= stall_out & (data_upstream & ~feedback_stall);
   end
end

localparam TYPE = MIN_FIFO_LATENCY < 1 ? (FIFO_DEPTH < 8 ? "zl_reg" : "zl_ram") : (MIN_FIFO_LATENCY < 3 ? (FIFO_DEPTH < 8 ? "ll_reg" : "ll_ram") : (FIFO_DEPTH < 8 ? "ll_reg" : "ram"));

  generate
    if ( STYLE == "TOKEN" )
    begin
      acl_token_fifo_counter 
      #(
        .DEPTH(FIFO_DEPTH)
       )
      fifo (
        .clock(clock),
        .resetn(resetn),
        .data_out(feedback_out),
        .valid_in(feedback_valid),
        .valid_out(feedback_valid_out),
        .stall_in(feedback_stall_in),
        .stall_out(feedback_stall)
      );
    end
    else
    begin
      acl_data_fifo #(
       .DATA_WIDTH(DATA_WIDTH),
       .DEPTH(((TYPE == "ram")  || (TYPE == "ll_ram") || (TYPE == "zl_ram")) ? FIFO_DEPTH + 1 : FIFO_DEPTH),
       .IMPL(TYPE),
       .ALLOW_FULL_WRITE(1)
       )
      fifo (
      .clock(clock),
      .resetn(resetn),
      .data_in(feedback),
      .data_out(feedback_out),
      .valid_in(feedback_valid),
      .valid_out(feedback_valid_out),
      .stall_in(feedback_stall_in),
      .stall_out(feedback_stall)
      );
    end
  endgenerate

endmodule


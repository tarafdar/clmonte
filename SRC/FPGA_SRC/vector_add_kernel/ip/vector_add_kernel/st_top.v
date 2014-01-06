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
    


/*
 * This is a variant of the staging reg - which allows you to add an
 * pre-initialized value; useful for feedback channels
 */
module init_reg
(
    clk, reset, i_init, i_data, i_valid, o_stall, o_data, o_valid, i_stall
);

/*************
* Parameters *
*************/
parameter WIDTH    = 32;
parameter INIT     = 0;
parameter INIT_VAL = 32'h0000000000000000;

/********
* Ports *
********/
// Standard global signals
input clk;
input reset;
input i_init;

// Upstream interface
input [WIDTH-1:0] i_data;
input i_valid;
output o_stall;

// Downstream interface
output [WIDTH-1:0] o_data;
output o_valid;
input i_stall;

/***************
* Architecture *
***************/
reg [WIDTH-1:0] r_data;
reg r_valid;

// Upstream
assign o_stall = r_valid;

// Downstream
assign o_data = (r_valid) ? r_data : i_data;
assign o_valid = (r_valid) ? r_valid : i_valid;

// Storage reg
always@(posedge clk or posedge reset)
begin
    if(reset == 1'b1)
    begin
        r_valid <= INIT;
        r_data  <= INIT_VAL;
    end
    else if (i_init) 
    begin
        r_valid <= INIT;
        r_data  <= INIT_VAL;
    end
    else
    begin
        if (~r_valid) r_data <= i_data;
        r_valid <= i_stall && (r_valid || i_valid);
    end
end

endmodule

//===----------------------------------------------------------------------===//
//
// Avalon Streaming Read Unit
//
//===----------------------------------------------------------------------===//
module st_read (
        clock,
        resetn,
        i_init,

        // input stream from kernel pipeline
        // this triggers the read request from the fifo
        i_predicate,
        i_valid,
        o_stall,

        // downstream (ouput), to kernel pipeline
        o_valid,
        i_stall,
        o_data,

        // input data from inter kernel pipeline
        i_fifodata,
        i_fifovalid,
        o_fifoready
        );

parameter DATA_WIDTH = 32;
parameter INIT = 0;
parameter INIT_VAL = 64'h0000000000000000;

input clock, resetn, i_stall, i_valid, i_fifovalid;
// init reinitializes the init fifo
input i_init;
output o_stall, o_valid, o_fifoready;
input i_predicate;
output [DATA_WIDTH-1:0] o_data;
input [DATA_WIDTH-1:0] i_fifodata;

wire feedback_downstream, data_downstream;
wire nop = i_predicate;

wire initvalid;
wire initready;

assign feedback_downstream = i_valid & ~nop & initvalid;
assign data_downstream = i_valid & nop;

wire init_reset;
wire r_o_stall;
wire init_val;

generate
if ( INIT ) begin
assign init_reset = ~resetn;
assign init_val   = i_init;
end
else begin
assign init_reset = ~resetn;
assign init_val   = 1'b0;
end
endgenerate

init_reg
  #( .WIDTH    ( DATA_WIDTH   ),
     .INIT     ( INIT     ),
     .INIT_VAL ( INIT_VAL ) )
reg_data ( 
      .clk     ( clock ),
      .reset   ( init_reset ),
      .i_init  ( init_val   ),
      .i_data  ( i_fifodata ),
      .i_valid ( i_fifovalid ), 
      .o_valid ( initvalid ),
      .o_data  ( o_data ),
      .o_stall ( r_o_stall ),
      .i_stall ( ~initready ) 
      );

assign o_fifoready = ~r_o_stall;

assign o_valid = feedback_downstream | data_downstream;
// assign o_data = i_fifodata ;

// stall upstream if
//   downstream is stalling (i_stall)
//   I'm waiting for data from fifo, don't stall if this read is
//   predicated
assign o_stall = ( i_valid & ~nop & ~initvalid ) | i_stall;

// don't accept data if:
//  downstream cannot accept data (i_stall)
//  data from upstream is selected (data_downstream)
//  no thread exists to read data (~i_valid)
// TODO: I should never set o_fifoready is this is
//       a fifo peek operation
assign initready = ~(i_stall  | data_downstream | ~i_valid); 

endmodule

//===----------------------------------------------------------------------===//
//
// Avalon Streaming Write Unit
// downstream are signals that continue into our "normal" pipeline.
//
//===----------------------------------------------------------------------===//
module st_write (
        clock,
        resetn,

        // interface from kernel pipeline, input stream
        i_predicate,
        i_data,
        i_valid,
        o_stall,

        // interface to kernel pipeline, downstream
        o_valid,
        i_stall,
        // data_out,

        // interface to kernel channel fifo, avalon master
        o_fifodata,
        o_fifovalid,
        i_fifoready
        );

parameter DATA_WIDTH = 32;

input clock, resetn, i_stall, i_valid, i_fifoready;
output o_stall, o_valid, o_fifovalid;
input [DATA_WIDTH-1:0] i_data;
input i_predicate;
output [DATA_WIDTH-1:0] o_fifodata;

wire nop;
assign nop = i_predicate;

assign o_fifovalid   = i_valid & ~nop & ~i_stall;

wire fifo_stall;

assign o_valid      = i_valid & (i_fifoready | nop);

assign o_stall = i_stall | (fifo_stall & (~nop) & i_valid) ;

assign fifo_stall = ~i_fifoready;
assign o_fifodata = i_data;

endmodule


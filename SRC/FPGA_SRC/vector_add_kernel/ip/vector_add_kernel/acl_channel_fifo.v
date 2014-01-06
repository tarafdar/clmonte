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
    



module acl_channel_fifo
#(
    parameter integer FIFO_DEPTH = 16,
    parameter integer DATA_W = 64              // > 0
)
(
    input logic clock,
    input logic resetn,

    input logic                 avst_in_valid,
    input logic    [DATA_W-1:0] avst_in_data,
    output logic                avst_in_ready,

    input logic                 avst_out_ready,
    output logic   [DATA_W-1:0] avst_out_data,
    output logic                avst_out_valid
);
    
    wire r_o_stall;
    wire r_o_valid;
    generate
    if (FIFO_DEPTH == 0)
    begin
        assign avst_out_data = avst_in_data;
        assign avst_out_valid = avst_in_valid;
        assign avst_in_ready = avst_out_ready;
    end
    else
    begin
        logic write;
        logic read;
        logic stall_out;
        logic valid_out;
        logic fifo_full;
        logic fifo_empty;

        assign write = avst_in_valid;
        assign read  = avst_out_ready;

        assign avst_in_ready  = ~r_o_stall;
        assign avst_out_valid = valid_out;

        wire  [DATA_W-1:0] r_out_data;
        acl_data_fifo
        #(
            .DATA_WIDTH(DATA_W),
            .DEPTH(FIFO_DEPTH),
            .IMPL("ram")
        )
        fifo
        (
            .clock     (clock),
            .resetn    (resetn), 
            .data_in   ( r_out_data ),
            .valid_in  ( r_o_valid ),
            .stall_out (stall_out),
            .data_out  (avst_out_data),
            .stall_in  ( ~read),
            .valid_out (valid_out),
            .empty     (fifo_empty),
            .full      (fifo_full)
        );
        assign r_out_data = avst_in_data;
        assign r_o_valid = write;
        assign r_o_stall = stall_out;
    end
    endgenerate
endmodule

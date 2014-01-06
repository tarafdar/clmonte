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
    


module vfabric_register(clock, resetn, i_settings, i_datain, i_datain_valid, o_datain_stall, o_dataout, o_dataout_valid, i_dataout_stall);

parameter DATA_WIDTH = 32;
parameter STALL_FREE = 0;

 input clock, resetn;
 input [DATA_WIDTH-1:0] i_settings;
 input [DATA_WIDTH-1:0] i_datain;
 input i_datain_valid;
 output o_datain_stall;
 
 output [DATA_WIDTH-1:0] o_dataout;
 output o_dataout_valid;
 input i_dataout_stall;

generate
  if( STALL_FREE == 1)
  begin
    reg [1:0] valid_counter;
    
    always @( posedge clock or negedge resetn )
    begin
      if( !resetn )
        valid_counter <= 2'b0;
      else
        valid_counter <= (i_datain_valid) ? valid_counter + 2'b01 : valid_counter;
    end

    assign o_dataout = {24'b0, i_settings[8*valid_counter +: 8]};
    assign o_dataout_valid = i_datain_valid;
    assign o_datain_stall = i_dataout_stall;
    
  end
  else
  begin
    assign o_dataout = i_settings;
    assign o_dataout_valid = i_datain_valid;
    assign o_datain_stall = i_dataout_stall;
  end
endgenerate
 
endmodule

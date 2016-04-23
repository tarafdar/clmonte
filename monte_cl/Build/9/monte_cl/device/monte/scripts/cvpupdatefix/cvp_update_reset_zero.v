// (C) 1992-2014 Altera Corporation. All rights reserved.                         
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
    

module cvp_update_reset # (
   parameter PIPE32_SIM_ONLY = 0
   ) (
   input pipe8_sim_only,
   input npor,
   input perst,
   output cvp_tx_analogreset
   );




//   assign serdes_tx_analogreset      = ((PIPE32_SIM_ONLY==1)||(pipe8_sim_only==1'b1))?1'b0           :(low_str(hip_hard_reset)=="disable")?~npor_int:1'b0;
//
   assign cvp_tx_analogreset      = 1'b0;



endmodule

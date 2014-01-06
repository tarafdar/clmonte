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
    


/***************************
* This module issues the workitems within a workgroup while supporting multiple
* workgroups in flight
****************************/

//
// It is designed to be embedded inside the function implementation,
// a function's implementation, i.e. along side the basic blocks.
//
// It has three states:
//
//    NEW:   Prepared to issue items from the next work group
//
//    ISSUE: Issuing items from a workgroup.
//          In this state, a new item will be issued in the same cycle
//          that both the following conditions hold:
//             The downstream block can accept it, and
//             The upstream function wrapper is ready to issue it.
//
//    WAIT:  Waiting for items to retire, to free up the next hw work group.
//
// All the work item intrinsic values are ready directly 
// from the wrapper's register set.
// We have to take care to ensure that the local_id_3 and hw_group_sel
// are consistent with the same-cycle values for work item values
// provided by the wrapper.
//
// TODO
// Out-of-order work-group exit may result in a sub-optimal number of
// work-groups being issued
//

module acl_embedded_workgroup_issuer_complex (clock, resetn, valid_in, stall_out, valid_entry, stall_entry,
  valid_exit, stall_exit, done_local_id_3, workgroup_size, local_id_3, live_thread_count,
  global_id_in, local_id_in, group_id_in, global_id_out, local_id_out, group_id_out);

  // Define module parameters
  parameter MAX_SIMULTANEOUS_WORKGROUPS = 2;
  parameter WORKGROUP_SIZE_ARG_WIDTH = 32; 
  
  input clock;
  input resetn;

  // Handshake for item entry into function.
  input valid_in;
  output stall_out;

  // Handshake with entry basic block
  output valid_entry;
  input stall_entry;

  // Observe threads exiting the function .
  // This is asserted when items are ready to be retired from the workgroup.
  input valid_exit;
  // This is asserted when downstream is not ready to retire an item from
  // the workgroup.
  input stall_exit;

  // Observe which work item is exiting the function.
  // Really only care about hw work group id, which is stored
  // in the most significant bits.
  input       [31:0] done_local_id_3;

  // Need workgroup_size to know when one workgroup ends
  // and another begins.
  input       [WORKGROUP_SIZE_ARG_WIDTH-1:0] workgroup_size;

  // Tell the downstream block the what hw group the issuing
  // The most significant bits of local_id_3 encode the hw group
  // to which the issuing item belongs.
  output      [31:0] local_id_3;


  // Tell downstream blocks how many work items are still alive
  // in this workgroup.
  output      [WORKGROUP_SIZE_ARG_WIDTH*MAX_SIMULTANEOUS_WORKGROUPS-1:0] live_thread_count;

  // Pass through global_id, local_id and group_id.
  input [31:0] global_id_in[2:0];
  input [31:0] local_id_in[2:0];
  input [31:0] group_id_in[2:0];
  output [31:0] global_id_out[2:0];
  output [31:0] local_id_out[2:0];
  output [31:0] group_id_out[2:0];

localparam LOG2_MAX_SIMULTANEOUS_WORKGROUPS=$clog2(MAX_SIMULTANEOUS_WORKGROUPS);

// States for state machine
//     NEW - Prepare to issue new workgroup
//     ISSUE - issue workitems in a workgroup (set valid high)
//     WAIT - Check if ready to accept new workgroup or if all done
localparam [2:0] wg_STATE_NEW=3'd0;
localparam [2:0] wg_STATE_ISSUE=3'd1;
localparam [2:0] wg_STATE_WAIT=3'd2;

reg  [2:0] present_state;
reg  [2:0] next_state;

////// Issuing info
// Which hw group id should we assign to incoming workgroups?
reg  [((LOG2_MAX_SIMULTANEOUS_WORKGROUPS>0) ? LOG2_MAX_SIMULTANEOUS_WORKGROUPS : 1)-1:0] hw_group_sel;
reg  [WORKGROUP_SIZE_ARG_WIDTH-1:0] num_issued_in_wg; // 0..workgroup_size-1

////// Retiring info
// Will an item be retired from the accelerator?
wire retire = valid_exit & ~stall_exit;
// When a thread exits the accelerator, determine which hw group id
// it was associated with.  Threads can exit in any order.
wire [((LOG2_MAX_SIMULTANEOUS_WORKGROUPS>0) ? LOG2_MAX_SIMULTANEOUS_WORKGROUPS : 1)-1:0] done_hw_group_sel;
reg  [WORKGROUP_SIZE_ARG_WIDTH-1:0] num_items_not_done[MAX_SIMULTANEOUS_WORKGROUPS-1:0];

generate
if (LOG2_MAX_SIMULTANEOUS_WORKGROUPS>0)
begin

   assign done_hw_group_sel= done_local_id_3[31:32-LOG2_MAX_SIMULTANEOUS_WORKGROUPS];

   // Update when transitioning out of issue state.  The scheduling is
   // simple, rotate through available hw groups and stall if the next one isn't
   // free (no bypassing).
   always @(posedge clock or negedge resetn)
   begin
    if (~(resetn)) 
      hw_group_sel<={LOG2_MAX_SIMULTANEOUS_WORKGROUPS{1'b0}};
    else if (present_state==wg_STATE_ISSUE && next_state==wg_STATE_WAIT)
    begin
      // Wrap-around if at the max. simultaneous workgroups value.
      if (hw_group_sel == MAX_SIMULTANEOUS_WORKGROUPS - 1)
        hw_group_sel <= '0;
      else
        hw_group_sel <= hw_group_sel + 'd1;
    end
   end

end
else
begin
   assign done_hw_group_sel=1'b0;
   always @(posedge clock) hw_group_sel=1'b0;
end
endgenerate

// Should we begin a new group?  We will not issue a new item at least until the next cycle.
wire begin_group = (present_state == wg_STATE_NEW)    & valid_in & (~stall_entry);
// Will we issue a new work item?
wire issue       = (present_state == wg_STATE_ISSUE)  & valid_in & (~stall_entry);

// Can consume items only when issuing state, and only when the downstream block
// is ready to consume.
assign stall_out = ~(present_state == wg_STATE_ISSUE) | stall_entry;
// Ready to send a new item along any time we're in issue state and we have
// a valid in.
assign valid_entry = (present_state == wg_STATE_ISSUE) & valid_in;

wire last_item_in_workgroup = (num_issued_in_wg == (workgroup_size - 2'b01));

// Have we just finished sending an entire workgroup?
wire workgroup_sent = (present_state==wg_STATE_ISSUE) && issue && last_item_in_workgroup;

// This must be consistent with the work item intrinsics from the wrapper.
// And it must be stable until the next item is issued.
wire [31:0] local_idx_of_issuing_item = num_issued_in_wg;

// Most significant bits encode hardware group id.  
// Least significant bits encode space index.
assign local_id_3 = { hw_group_sel, local_idx_of_issuing_item[31-LOG2_MAX_SIMULTANEOUS_WORKGROUPS:0] };

// The ids (global, local and group) are simply passed through.
assign global_id_out = global_id_in;
assign local_id_out = local_id_in;
assign group_id_out = group_id_in;

// Update num_issued_in_wg.
// Its value is kept stable until the next hw workgroup starts
always @(posedge clock or negedge resetn) begin
   if ( ~resetn )
      num_issued_in_wg <= {WORKGROUP_SIZE_ARG_WIDTH{1'b0}};
   else if ( begin_group )
      num_issued_in_wg <= {WORKGROUP_SIZE_ARG_WIDTH{1'b0}};
   else if ( issue )
      num_issued_in_wg <= num_issued_in_wg + 1'b1;
   else
      num_issued_in_wg <= num_issued_in_wg;
end



// Check if next hw group is free to accept a new workgroup.  WARNING: This
// signal is only valid in the WAIT state!
wire next_hw_group_free = ~(|num_items_not_done[hw_group_sel]);

// State machine
//     NEW - Prepare to issue new workgroup
//     ISSUE - issue workitems in a workgroup (set valid high)
//     WAIT - Check if ready to accept new workgroup or if all done
always@*
begin
  next_state = wg_STATE_NEW;
  case (present_state)
    wg_STATE_NEW:
      next_state = (begin_group) ? wg_STATE_ISSUE : wg_STATE_NEW;
    wg_STATE_ISSUE:
      next_state = (workgroup_sent) ? wg_STATE_WAIT : wg_STATE_ISSUE;
    wg_STATE_WAIT:
      next_state = (next_hw_group_free) ? wg_STATE_NEW : wg_STATE_WAIT;
  endcase
end

always @(posedge clock or negedge resetn)
begin
  if (~(resetn))
    present_state <= wg_STATE_NEW;
  else
    present_state <= next_state;
end


// Compute num_items_not_done for each HW group
generate
genvar i;
  for (i=0; i<MAX_SIMULTANEOUS_WORKGROUPS; i=i+1)
  begin : numdone_gen
    always @(posedge clock or negedge resetn)
    begin
      if (~(resetn))
        num_items_not_done[i] <= {WORKGROUP_SIZE_ARG_WIDTH{1'b0}};
      else 
        case (present_state)
          wg_STATE_NEW:
            if ( retire && (done_hw_group_sel==i) ) // retiring from hw group i
              num_items_not_done[i] <= (num_items_not_done[i] - 2'b01);

          wg_STATE_ISSUE:
            // The trick is we might be issuing items into the same hw group that might be retiring an item
            if (       ( issue  && (hw_group_sel==i) )          // issuing into hw group i
                  &&  ~( retire && (done_hw_group_sel==i) )  // not retiring from hw group i
               )
               num_items_not_done[i] <= (num_items_not_done[i] + 2'b01);
            else if ( ~( issue  && (hw_group_sel==i) )          // not issuing into hw group i
                  &&   ( retire && (done_hw_group_sel==i) )  // retiring from hw group i
               )
               num_items_not_done[i] <= (num_items_not_done[i] - 2'b01);

          wg_STATE_WAIT: // same as NEW
            if ( retire && done_hw_group_sel==i ) // retiring from hw group i
              num_items_not_done[i] <= (num_items_not_done[i] - 2'b01);
        endcase
    end
  end
endgenerate

// We need to track how many threads are alive for each HW workgroup.  This is
// because barriers apply only to live threads (ie. some threads can exit the
// loop while the remaining go into a barrier).
generate
genvar j;
  for (j=0; j<MAX_SIMULTANEOUS_WORKGROUPS; j=j+1)
  begin : livethread_gen
    assign live_thread_count[j*WORKGROUP_SIZE_ARG_WIDTH +: WORKGROUP_SIZE_ARG_WIDTH] = num_items_not_done[j];
  end
endgenerate

endmodule
// vim:set filetype=verilog:

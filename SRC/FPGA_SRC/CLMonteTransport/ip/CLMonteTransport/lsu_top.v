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
    


//
// Top level load/store unit
//
// Attributes of load/store units
//  Coalesced:  Accesses to neighbouring memory locations are grouped together
//              to improve efficiency and efficiently utilize memory bandwidth.
//  Hazard-Safe:The LSU is not susceptable to data hazards.
//  Ordered:    The LSU requires accesses to be in-order to properly coalesce.
//  Pipeline:   The LSU can handle multiple requests at a time without
//              stalling.  Improves throughput.
//
// Supports the following memory access patterns:
//  Simple    - STYLE="SIMPLE"
//              Coalesced: No, Ordered: N/A, Hazard-Safe: Yes, Pipelined, No
//              Simple un-pipelined memory access.  Low throughput.
//  Pipelined - STYLE="PIPELINED"
//              Coalesced: No, Ordered: N/A, Hazard-Safe: Yes, Pipelined: Yes
//              Requests are submitted as soon as they are received.
//              Pipelined access to memory so multiple requests can be 
//              in flight at a time.
//  Coalesced - STYLE="BASIC-COALESCED"
//   "basic"    Coalesced: Yes, Ordered: Yes, Hazard-Safe: Yes, Pipelined: Yes
//              Requests are submitted as soon as possible to memory, stalled
//              requests are coalesced with neighbouring requests if they
//              access the same page of memory.
//  Coalesced - STYLE="BURST-COALESCED"
//   "burst"    Coalesced: Yes, Ordered: Yes, Hazard-Safe: Yes, Pipelined: Yes
//              Requests are buffered until the biggest possible burst can
//              be made.
//  Streaming - STYLE="STREAMING"
//              Coalesced: Yes, Ordered: Yes, Hazard-Safe: No, Pipelined: ?
//              A FIFO is instantiated which burst reads large blocks from 
//              memory to keep the FIFO full of valid data.  This block can 
//              only be used if accesses are in-order, and addresses can be
//              simply calculated from (base_address + n * word_width).  The
//              block has no built-in hazard protection.
//  Atomic - STYLE="ATOMIC-PIPELINED"
//"pipelined"
//              Coalesced: No, Ordered: N/A, Hazard-Safe: Yes, Pipelined: Yes
//              Atomic: Yes
//              Requests are submitted as soon as they are received.
//              Pipelined access to memory so multiple requests can be 
//              in flight at a time.
//              Response is returned as soon as read is complete, 
//              write is issued subsequently by the atomic module at the end
//              of arbitration.

module lsu_top
(
    clock, clock2x, resetn, stream_base_addr, stream_size, stream_reset, i_atomic_op, o_stall, 
    i_valid, i_address, i_writedata, i_cmpdata, i_predicate, i_bitwiseor, i_stall, o_valid, o_readdata, avm_address, 
    avm_read, avm_readdata, avm_write, avm_writeack, avm_writedata, avm_byteenable, 
    avm_waitrequest, avm_readdatavalid, avm_burstcount,
    o_active,
    o_input_fifo_depth,
    o_writeack,
    flush
);

/*************
* Parameters *
*************/
parameter STYLE="PIPELINED"; // The LSU style to use (see style list above)
parameter AWIDTH=32;         // Address width (32-bits for Avalon)
parameter ATOMIC_WIDTH=6;    // Width of operation operation indices
parameter WIDTH_BYTES=4;     // Width of the request (bytes)
parameter MWIDTH_BYTES=32;   // Width of the global memory bus (bytes)
parameter WRITEDATAWIDTH_BYTES=32;  // Width of the readdata/writedata signals, 
                                    // may be larger than MWIDTH_BYTES for atomics
parameter ALIGNMENT_BYTES=2; // Request address alignment (bytes)
parameter READ=1;            // Read or write?
parameter ATOMIC=0;          // Atomic?
parameter BURSTCOUNT_WIDTH=6;// Determines max burst size
// Why two latencies? E.g. A streaming unit prefetches data, its latency to
// the kernel is very low because data is for the most part ready and waiting.
// But the lsu needs to know how much data to buffer to hide the latency to
// memory, hence the memory side latency.
parameter KERNEL_SIDE_MEM_LATENCY=1;  // Effective Latency in cycles as seen by the kernel pipeline
parameter MEMORY_SIDE_MEM_LATENCY=1;  // Latency in cycles between LSU and memory
parameter USE_WRITE_ACK=0;   // Enable the write-acknowledge signal
parameter USECACHING=0;
parameter CACHESIZE=1024;
parameter PROFILE_ADDR_TOGGLE=0;
parameter USEINPUTFIFO=1;        // FIXME specific to lsu_pipelined
parameter USEOUTPUTFIFO=1;       // FIXME specific to lsu_pipelined
parameter FORCE_NOP_SUPPORT=0;   // Stall free pipeline doesn't want the NOP fifo
parameter HIGH_FMAX=1;       // Enable optimizations for high Fmax

// Profiling
parameter ACL_PROFILE=0;      // Set to 1 to enable stall/valid profiling
parameter ACL_PROFILE_ID=1;   // Each LSU needs a unique ID
parameter ACL_PROFILE_COUNTER_WIDTH=64;

// Verilog readability and parsing only - no functional purpose
parameter ADDRSPACE=0;

// Local memory parameters
parameter ENABLE_BANKED_MEMORY=0;// Flag enables address permutation for banked local memory config
parameter ABITS_PER_LMEM_BANK=0; // Used when permuting lmem address bits to stride across banks
parameter NUMBER_BANKS=1;        // Number of memory banks - used in address permutation (1-disable)
parameter LMEM_ADDR_PERMUTATION_STYLE=0; // Type of address permutation (currently unused)
// The following localparams have if conditions, and the second is named
// "HACKED..." because address bit permutations are controlled by the
// ENABLE_BANKED_MEMORY parameter.  The issue is that this forms the select
// input of a MUX (if statement), and synthesis evaluates both inputs.
// When not using banked memory, the bit select ranges don't make sense on
// the input that isn't used, so we need to hack them in the non-banked case
// to get through ModelSim and Quartus.
localparam BANK_SELECT_BITS = (ENABLE_BANKED_MEMORY==1) ? $clog2(NUMBER_BANKS) : 1; // Bank select bits in address permutation
localparam HACKED_ABITS_PER_LMEM_BANK = (ENABLE_BANKED_MEMORY==1) ? ABITS_PER_LMEM_BANK : $clog2(MWIDTH_BYTES)+1;

// Parameter limitations:
//    AWIDTH: Only tested with 32-bit addresses
//    WIDTH_BYTES: Must be a power of two
//    MWIDTH_BYTES: Must be a power of 2 >= WIDTH_BYTES
//    ALIGNMENT_BYTES: Must be a power of 2 satisfying,
//                     WIDTH_BYTES <= ALIGNMENT_BYTES <= MWIDTH_BYTES
//
//    The width and alignment restrictions ensure we never try to read a word
//    that strides across two "pages" (MWIDTH sized words)

// TODO: Convert these back into localparams when the back-end supports it
parameter WIDTH=8*WIDTH_BYTES;                      // Width in bits
parameter MWIDTH=8*MWIDTH_BYTES;                    // Width in bits
parameter WRITEDATAWIDTH=8*WRITEDATAWIDTH_BYTES;              // Width in bits
localparam ALIGNMENT_ABITS=$clog2(ALIGNMENT_BYTES); // Address bits to ignore
localparam LSU_CAPACITY=256;   // Maximum number of 'in-flight' load/store operations

// Performance monitor signals
parameter INPUTFIFO_USEDW_MAXBITS=8;

// LSU unit properties
localparam ATOMIC_PIPELINED_LSU=(STYLE=="ATOMIC-PIPELINED");
localparam PIPELINED_LSU=( (STYLE=="PIPELINED") || (STYLE=="BASIC-COALESCED") || (STYLE=="BURST-COALESCED") || (STYLE=="BURST-NON-ALIGNED") );
localparam SUPPORTS_NOP=( (STYLE=="STREAMING") || (STYLE=="SEMI-STREAMING") || (STYLE=="BURST-NON-ALIGNED") || FORCE_NOP_SUPPORT==1);
localparam SUPPORTS_BURSTS=( (STYLE=="STREAMING") || (STYLE=="BURST-COALESCED") || (STYLE=="SEMI-STREAMING") || (STYLE=="BURST-NON-ALIGNED") );

/********
* Ports *
********/
// Standard global signals
input clock;
input clock2x;
input resetn;
input flush;

// Streaming interface signals
input [AWIDTH-1:0] stream_base_addr;
input [31:0] stream_size;
input stream_reset;

// Atomic interface
input [WIDTH-1:0] i_cmpdata; // only used by atomic_cmpxchg
input [ATOMIC_WIDTH-1:0] i_atomic_op;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;
input [WIDTH-1:0] i_writedata;
input i_predicate;
input [AWIDTH-1:0] i_bitwiseor;

// Downstream interface
input i_stall;
output o_valid;
output [WIDTH-1:0] o_readdata;

// Avalon interface
output [AWIDTH-1:0] avm_address;
output avm_read;
input [WRITEDATAWIDTH-1:0] avm_readdata;
output avm_write;
input avm_writeack;
output o_writeack;
output [WRITEDATAWIDTH-1:0] avm_writedata;
output [WRITEDATAWIDTH_BYTES-1:0] avm_byteenable;
input avm_waitrequest;
input avm_readdatavalid;
output [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

output reg o_active;

// For profiling/performance monitor
output [INPUTFIFO_USEDW_MAXBITS-1:0] o_input_fifo_depth;

wire lsu_active;

// For handling dependents of this lsu
assign o_writeack = avm_writeack;

// If this is a banked local memory LSU, then permute address bits so that
// consective words in memory are in different banks.  Do this by
// taking the k lowest bits of the word-address and shifting them to the top
// of the aggregate local memory address width.  The number of bits k
// corresponds to the number of banks parameter.

localparam MWIDTH_BYTES_CLIP = (MWIDTH_BYTES==1) ? 2 : MWIDTH_BYTES; //to get around modelsim looking at addr[-1:0] if MWIDTH_BYTES==1
function [AWIDTH-1:0] permute_addr ( input [AWIDTH-1:0] addr);
  if (ENABLE_BANKED_MEMORY==1)
  begin
    if (MWIDTH_BYTES==1) begin
      permute_addr= {
        addr[(AWIDTH-1) : (HACKED_ABITS_PER_LMEM_BANK+BANK_SELECT_BITS)], // High order bits unchanged
        addr[($clog2(MWIDTH_BYTES)+BANK_SELECT_BITS-1) : $clog2(MWIDTH_BYTES)], // Bank select from lsbits
        addr[(HACKED_ABITS_PER_LMEM_BANK + BANK_SELECT_BITS-1) : ($clog2(MWIDTH_BYTES) + BANK_SELECT_BITS)]
        };
    end
    else begin
      permute_addr= {
        addr[(AWIDTH-1) : (HACKED_ABITS_PER_LMEM_BANK+BANK_SELECT_BITS)], // High order bits unchanged
        addr[($clog2(MWIDTH_BYTES)+BANK_SELECT_BITS-1) : $clog2(MWIDTH_BYTES)], // Bank select from lsbits
        addr[(HACKED_ABITS_PER_LMEM_BANK + BANK_SELECT_BITS-1) : ($clog2(MWIDTH_BYTES) + BANK_SELECT_BITS)],
        addr[($clog2(MWIDTH_BYTES_CLIP)-1) : 0]         // Byte address within a word
        };
    end
   end
   else
   begin
     permute_addr= addr;
   end
endfunction

wire [AWIDTH-1:0] avm_address_raw;
assign avm_address=permute_addr(avm_address_raw);

/***************
* Architecture *
***************/

// Tie off the unused read/write signals
generate
// atomics dont have unused signals
if(ATOMIC==0) begin
  if(READ==1)
  begin
     assign avm_write = 1'b0;
     //assign avm_writedata = {MWIDTH{1'bx}};
     assign avm_writedata = {MWIDTH{1'b0}}; // make writedata 0 because it is used by atomics
  end
  else // WRITE
  begin
    assign avm_read = 1'b0;
  end
end
else begin //ATOMIC
  assign avm_write = 1'b0;
end
endgenerate

// Write acknowledge support: If WRITEACK is not to be supported, than assume
// that a write is fully completed as soon as it is accepted by the fabric.
// Otherwise, wait for the writeack signal to return.
wire lsu_writeack;
generate
if(USE_WRITE_ACK==1)
begin
   assign lsu_writeack = avm_writeack;
end
else
begin
   assign lsu_writeack = avm_write && !avm_waitrequest;
end
endgenerate

// NOP support: The least-significant address bit indicates if this is a NOP
// instruction (i.e. we do not wish a read/write to be performed).  
// Appropriately adjust the valid and stall inputs to the core LSU block to
// ensure NOP instructions are not executed and preserve their ordering with
// other threads.
wire lsu_i_valid;
wire lsu_o_valid;
wire lsu_i_stall;
wire lsu_o_stall;
wire [AWIDTH-1:0] address;
wire nop;

generate
if(SUPPORTS_NOP)
begin
   // Module intrinsicly supports NOP operations, just pass them on through
   assign lsu_i_valid = i_valid;
   assign lsu_i_stall = i_stall;
   assign o_valid = lsu_o_valid;
   assign o_stall = lsu_o_stall;
   assign address = i_address | i_bitwiseor;
end
else if(PIPELINED_LSU || ATOMIC_PIPELINED_LSU)
begin
   // No built-in NOP support.  Pipelined LSUs without NOP support need us to 
   // build a fifo along side the core LSU to track NOP instructions
   wire nop_fifo_empty;
   wire nop_fifo_full;
   wire nop_next;

   assign nop = i_predicate;
   assign address = i_address | i_bitwiseor;

   // Store the NOP status flags along side the core LSU
   // Assume (TODO eliminate this assumption?) that KERNEL_SIDE_MEM_LATENCY is the max 
   // number of simultaneous requests in flight for the LSU. The NOP FIFO will
   // will be sized to KERNEL_SIDE_MEM_LATENCY+1 to prevent stalls when the LSU is
   // full.
   //
   // For smaller latency values, use registers to implement the FIFO.
   if(KERNEL_SIDE_MEM_LATENCY <= 64)
   begin
      acl_ll_fifo #(
         .WIDTH(1),
         .DEPTH(KERNEL_SIDE_MEM_LATENCY+1)
      ) nop_fifo (
         .clk(clock),
         .reset(~resetn),
         .data_in(nop),
         .write(i_valid && !o_stall),
         .data_out(nop_next),
         .read(o_valid && !i_stall),
         .full(nop_fifo_full),
         .empty(nop_fifo_empty)
      );
   end
   else
   begin
      scfifo #(
         .add_ram_output_register( "OFF" ),
         .intended_device_family( "Stratix IV" ),
         .lpm_numwords( KERNEL_SIDE_MEM_LATENCY+1 ),
         .lpm_showahead( "ON" ),
         .lpm_type( "scfifo" ),
         .lpm_width( 1 ),
         .lpm_widthu( $clog2(KERNEL_SIDE_MEM_LATENCY+1) ),
         .overflow_checking( "OFF" ),
         .underflow_checking( "OFF" )
      ) nop_fifo (
         .clock(clock),
         .data(nop),
         .rdreq(o_valid && !i_stall),
         .wrreq(i_valid && !o_stall),
         .empty(nop_fifo_empty),
         .full(nop_fifo_full),
         .q(nop_next),
         .aclr(!resetn),
         .almost_full(),
         .almost_empty(),
         .usedw(),
         .sclr()
      );
   end

   // Logic to prevent NOP instructions from entering the core
   assign lsu_i_valid = !nop && i_valid && !nop_fifo_full;
   assign lsu_i_stall = nop_fifo_empty || nop_next || i_stall;

   // Logic to generate the valid bit for NOP instructions that have bypassed
   // the LSU.  The instructions must be kept in order so they are consistant
   // with data propagating through pipelines outside of the LSU.
   assign o_valid = (lsu_o_valid || nop_next) && !nop_fifo_empty;
   assign o_stall = nop_fifo_full || lsu_o_stall;
end
else
begin
   // An unpipelined LSU will only have one active request at a time.  We just
   // need to track whether there is a pending request in the LSU core and
   // appropriately bypass the core with NOP requests while preserving the
   // thread ordering.  (A NOP passes straight through to the downstream 
   // block, unless there is a pending request in the block, in which case
   // we stall until the request is complete).
   reg pending;
   always@(posedge clock or negedge resetn)
   begin
      if(resetn == 1'b0)
         pending <= 1'b0;
      else
         pending <= pending ? ((lsu_i_valid && !lsu_o_stall) || !(lsu_o_valid && !lsu_i_stall)) :
                              ((lsu_i_valid && !lsu_o_stall) && !(lsu_o_valid && !lsu_i_stall));
   end

   assign nop = i_predicate;
   assign address = i_address | i_bitwiseor;

   assign lsu_i_valid = i_valid && !nop;
   assign lsu_i_stall = i_stall;
   assign o_valid = lsu_o_valid || (!pending && i_valid && nop);
   assign o_stall = lsu_o_stall || (pending && nop);
end
endgenerate

// Styles with no burst support require burstcount=1
generate
if(!SUPPORTS_BURSTS)
begin
   assign avm_burstcount = 1;
end
endgenerate

// Profiling signals.
logic [63:0] req_cache_hit_count;
logic [63:0] extra_unaligned_reqs;

// Generate different architectures based on the STYLE parameter
generate

////////////////
// Simple LSU //
////////////////
if(STYLE=="SIMPLE")
begin
    if(READ == 1)
    begin
        lsu_simple_read #(
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .HIGH_FMAX(HIGH_FMAX)
        ) simple_read (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_active(lsu_active),
            .o_readdata(o_readdata),
            .avm_address(avm_address_raw),
            .avm_read(avm_read),
            .avm_readdata(avm_readdata),
            .avm_waitrequest(avm_waitrequest),
            .avm_byteenable(avm_byteenable),
            .avm_readdatavalid(avm_readdatavalid)
        );
    end
    else
    begin
        lsu_simple_write #(
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS)
        ) simple_write (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_writedata(i_writedata),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_write(avm_write),
            .avm_writeack(lsu_writeack),
            .avm_writedata(avm_writedata),
            .avm_byteenable(avm_byteenable),
            .avm_waitrequest(avm_waitrequest)
        );
    end
end

///////////////
// Pipelined //
///////////////
else if(STYLE=="PIPELINED")
begin
    if(READ == 1)
    begin
        lsu_pipelined_read #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .USEINPUTFIFO(USEINPUTFIFO),
            .USEOUTPUTFIFO(USEOUTPUTFIFO)
        ) pipelined_read (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_readdata(o_readdata),
            .o_input_fifo_depth(o_input_fifo_depth),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_read(avm_read),
            .avm_readdata(avm_readdata),
            .avm_waitrequest(avm_waitrequest),
            .avm_byteenable(avm_byteenable),
            .avm_readdatavalid(avm_readdatavalid)
        );
    end
    else
    begin
        lsu_pipelined_write #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .USEINPUTFIFO(USEINPUTFIFO)
        ) pipelined_write (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_writedata(i_writedata),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_input_fifo_depth(o_input_fifo_depth),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_write(avm_write),
            .avm_writeack(lsu_writeack),
            .avm_writedata(avm_writedata),
            .avm_byteenable(avm_byteenable),
            .avm_waitrequest(avm_waitrequest)
        );
    end
end

//////////////////////
// Atomic Pipelined //
//////////////////////
else if(STYLE=="ATOMIC-PIPELINED")
begin
    lsu_atomic_pipelined #(
           .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
           .AWIDTH(AWIDTH),
           .WIDTH_BYTES(WIDTH_BYTES),
           .MWIDTH_BYTES(MWIDTH_BYTES),
           .WRITEDATAWIDTH_BYTES(WRITEDATAWIDTH_BYTES),
           .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
           .USEINPUTFIFO(USEINPUTFIFO),
           .USEOUTPUTFIFO(USEOUTPUTFIFO),
           .ATOMIC_WIDTH(ATOMIC_WIDTH)
    ) atomic_pipelined (
           .clk(clock),
           .reset(!resetn),
           .o_stall(lsu_o_stall),
           .i_valid(lsu_i_valid),
           .i_address(address),
           .i_stall(lsu_i_stall),
           .o_valid(lsu_o_valid),
           .o_readdata(o_readdata),
           .o_input_fifo_depth(o_input_fifo_depth),
           .o_active(lsu_active),
           .avm_address(avm_address_raw),
           .avm_read(avm_read),
           .avm_readdata(avm_readdata),
           .avm_waitrequest(avm_waitrequest),
           .avm_byteenable(avm_byteenable),
           .avm_readdatavalid(avm_readdatavalid),
           .i_atomic_op(i_atomic_op),
           .i_writedata(i_writedata),
           .i_cmpdata(i_cmpdata),
           .avm_writeack(lsu_writeack),
           .avm_writedata(avm_writedata)
    );
end

/////////////////////
// Basic Coalesced //
/////////////////////
else if(STYLE=="BASIC-COALESCED")
begin
    if(READ == 1)
    begin
        lsu_basic_coalesced_read #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS)
        ) basic_coalesced_read (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_readdata(o_readdata),
            .avm_address(avm_address_raw),
            .avm_read(avm_read),
            .avm_readdata(avm_readdata),
            .avm_waitrequest(avm_waitrequest),
            .avm_byteenable(avm_byteenable),
            .avm_readdatavalid(avm_readdatavalid)
        );
    end
    else
    begin
        lsu_basic_coalesced_write #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS)
        ) basic_coalesced_write (
            .clk(clock),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_writedata(i_writedata),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_write(avm_write),
            .avm_writeack(lsu_writeack),
            .avm_writedata(avm_writedata),
            .avm_byteenable(avm_byteenable),
            .avm_waitrequest(avm_waitrequest)
        );
    end
end

/////////////////////
// Burst Coalesced //
/////////////////////
else if(STYLE=="BURST-COALESCED")
begin
    if(READ == 1)
    begin
        lsu_bursting_read #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
            .USECACHING(USECACHING),
            .HIGH_FMAX(HIGH_FMAX),
            .ACL_PROFILE(ACL_PROFILE),
            .CACHE_SIZE_N(CACHESIZE)
        ) bursting_read (
            .clk(clock),
            .clk2x(clock2x),
            .reset(!resetn),
            .flush(flush),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_readdata(o_readdata),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_read(avm_read),
            .avm_readdata(avm_readdata),
            .avm_waitrequest(avm_waitrequest),
            .avm_byteenable(avm_byteenable),
            .avm_burstcount(avm_burstcount),
            .avm_readdatavalid(avm_readdatavalid),
            .req_cache_hit_count(req_cache_hit_count)
        );
    end
    else
    begin
        // Non-writeack stores are similar to streaming, where the pipeline
        // needs only few threads which just drop off data, and internally the
        // LSU must account for arbitration contention and other delays.
        lsu_bursting_write #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
            .USE_WRITE_ACK(USE_WRITE_ACK),
            .HIGH_FMAX(HIGH_FMAX)
        ) bursting_write (
            .clk(clock),
            .clk2x(clock2x),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_writedata(i_writedata),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_active(lsu_active),
            .word_byte_enable({WIDTH_BYTES{1'b1}}),
            .word_bit_enable({(WIDTH_BYTES*8){1'b1}}),
            .avm_address(avm_address_raw),
            .avm_write(avm_write),
            .avm_writeack(lsu_writeack),
            .avm_writedata(avm_writedata),
            .avm_byteenable(avm_byteenable),
            .avm_burstcount(avm_burstcount),
            .avm_waitrequest(avm_waitrequest)
        );
    end
end


/////////////////////////////////
// Burst Coalesced Non Aligned //
/////////////////////////////////
else if(STYLE=="BURST-NON-ALIGNED")
begin
    if(READ == 1)
    begin
        lsu_non_aligned_read #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
            .USECACHING(USECACHING),
            .HIGH_FMAX(HIGH_FMAX),
            .ACL_PROFILE(ACL_PROFILE)
        ) bursting_non_aligned_read (
            .clk(clock),
            .clk2x(clock2x),
            .reset(!resetn),
            .flush(flush),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_nop(i_predicate),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_readdata(o_readdata),
            .o_active(lsu_active),
            .avm_address(avm_address_raw),
            .avm_read(avm_read),
            .avm_readdata(avm_readdata),
            .avm_waitrequest(avm_waitrequest),
            .avm_byteenable(avm_byteenable),
            .avm_burstcount(avm_burstcount),
            .avm_readdatavalid(avm_readdatavalid),
            .extra_unaligned_reqs(extra_unaligned_reqs),
            .req_cache_hit_count(req_cache_hit_count)
        );
    end
    else
    begin
        lsu_non_aligned_write #(
            .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
            .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
            .AWIDTH(AWIDTH),
            .WIDTH_BYTES(WIDTH_BYTES),
            .MWIDTH_BYTES(MWIDTH_BYTES),
            .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
            .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
            .USE_WRITE_ACK(USE_WRITE_ACK),
            .HIGH_FMAX(HIGH_FMAX)
        ) bursting_non_aligned_write (
            .clk(clock),
            .clk2x(clock2x),
            .reset(!resetn),
            .o_stall(lsu_o_stall),
            .i_valid(lsu_i_valid),
            .i_address(address),
            .i_nop(i_predicate),
            .i_writedata(i_writedata),
            .i_stall(lsu_i_stall),
            .o_valid(lsu_o_valid),
            .o_active(lsu_active),
            //.word_byte_enable({WIDTH_BYTES{1'b1}}),
            //.word_bit_enable({(WIDTH_BYTES*8){1'b1}}),
            .avm_address(avm_address_raw),
            .avm_write(avm_write),
            .avm_writeack(lsu_writeack),
            .avm_writedata(avm_writedata),
            .avm_byteenable(avm_byteenable),
            .avm_burstcount(avm_burstcount),
            .avm_waitrequest(avm_waitrequest)
        );
    end
end
///////////////
// Streaming //
///////////////
else if(STYLE=="STREAMING")
begin
   if(READ==1)
   begin
      lsu_streaming_read #(
         .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
         .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
         .AWIDTH(AWIDTH),
         .WIDTH_BYTES(WIDTH_BYTES),
         .MWIDTH_BYTES(MWIDTH_BYTES),
         .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
         .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH)
      ) streaming_read (
         .clk(clock),
         .reset(!resetn),
         .o_stall(lsu_o_stall),
         .i_valid(lsu_i_valid),
         .i_stall(lsu_i_stall),
         .o_valid(lsu_o_valid),
         .o_readdata(o_readdata),
         .o_active(lsu_active),
         .i_nop(i_predicate),
         .base_address(stream_base_addr),
         .size(stream_size),
         .avm_address(avm_address_raw),
         .avm_burstcount(avm_burstcount),
         .avm_read(avm_read),
         .avm_readdata(avm_readdata),
         .avm_waitrequest(avm_waitrequest),
         .avm_byteenable(avm_byteenable),
         .avm_readdatavalid(avm_readdatavalid)
      );
   end
   else
   begin
     lsu_streaming_write #(
         .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
         .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
         .AWIDTH(AWIDTH),
         .WIDTH_BYTES(WIDTH_BYTES),
         .MWIDTH_BYTES(MWIDTH_BYTES),
         .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
         .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH)
     ) streaming_write (
         .clk(clock),
         .reset(!resetn),
         .o_stall(lsu_o_stall),
         .i_valid(lsu_i_valid),
         .i_stall(lsu_i_stall),
         .o_valid(lsu_o_valid),
         .o_active(lsu_active),
         .i_writedata(i_writedata),
         .i_nop(i_predicate),
         .base_address(stream_base_addr),
         .size(stream_size),
         .avm_address(avm_address_raw),
         .avm_burstcount(avm_burstcount),
         .avm_write(avm_write),
         .avm_writeack(lsu_writeack),
         .avm_writedata(avm_writedata),
         .avm_byteenable(avm_byteenable),
         .avm_waitrequest(avm_waitrequest)
     );
   end
end
////////////////////
// SEMI-Streaming //
////////////////////
else if(STYLE=="SEMI-STREAMING")
begin
   if(READ==1)
   begin
      lsu_read_cache #(
         .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
         .AWIDTH(AWIDTH),
         .WIDTH_BYTES(WIDTH_BYTES),
         .MWIDTH_BYTES(MWIDTH_BYTES),
         .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
         .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
         .REQUESTED_SIZE(CACHESIZE)
      ) read_cache (
         .clk(clock),
         .reset(!resetn),
         .flush(flush),
         .o_stall(lsu_o_stall),
         .i_valid(lsu_i_valid),
         .i_address(address),
         .i_stall(lsu_i_stall),
         .o_valid(lsu_o_valid),
         .o_readdata(o_readdata),
         .o_active(lsu_active),
         .i_nop(i_predicate),
         .avm_address(avm_address_raw),
         .avm_burstcount(avm_burstcount),
         .avm_read(avm_read),
         .avm_readdata(avm_readdata),
         .avm_waitrequest(avm_waitrequest),
         .avm_byteenable(avm_byteenable),
         .avm_readdatavalid(avm_readdatavalid)
      );
   end
end
endgenerate

always@(posedge clock or negedge resetn)
   if (!resetn)
      o_active <= 1'b0;
    else
      o_active <= lsu_active;

// Profile the valids and stalls of the LSU
generate
if(ACL_PROFILE==1)
begin
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] bw;
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] total_ivalid;
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] total_req;
 
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] clock_cycles_total;
   logic [15:0] clock_cycles_no_valid;

   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] i_stall_count;
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] o_stall_count;

   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] avm_read_count;
   logic [ACL_PROFILE_COUNTER_WIDTH-1:0] avm_burstcount_total;

   profile_sink #(.SINK_ID(ACL_PROFILE_ID),.PROBE_WIDTH(4*ACL_PROFILE_COUNTER_WIDTH)) ps (
      .probe( {bw,total_ivalid,total_req,clock_cycles_total} ));
   
   // TODO generalize profile id assignment
   profile_sink #(
      .SINK_ID(ACL_PROFILE_ID + 64),
      .PROBE_WIDTH($bits(o_stall_count) + $bits(i_stall_count) + $bits(extra_unaligned_reqs) + $bits(req_cache_hit_count))
      ) ps2 (
      .probe( {o_stall_count, i_stall_count, extra_unaligned_reqs, req_cache_hit_count} ));

   profile_sink #(
      .SINK_ID(ACL_PROFILE_ID + 96),
      .PROBE_WIDTH($bits(avm_read_count) + $bits(avm_burstcount_total))
      ) ps3 (
      .probe( {avm_read_count, avm_burstcount_total} ));

   always@(posedge clock or negedge resetn)
      if (!resetn)
      begin
         bw <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
         clock_cycles_total <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
         total_ivalid <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
         total_req <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
         clock_cycles_no_valid <= 16'b0;
         i_stall_count <= '0;
         o_stall_count <= '0;
         avm_read_count <= '0;
         avm_burstcount_total <= '0;
      end
      else
      begin
         if (flush)
         begin
            bw <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
            total_ivalid <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
            clock_cycles_total <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
            total_req <= {ACL_PROFILE_COUNTER_WIDTH{1'b0}};
            clock_cycles_no_valid <= 16'b0;
            i_stall_count <= '0;
            o_stall_count <= '0;
            avm_read_count <= '0;
            avm_burstcount_total <= '0;
         end
         else
         begin
            bw <= bw + ((READ==1) ? (avm_readdatavalid ? MWIDTH_BYTES : 0) : (avm_write & ~avm_waitrequest ? MWIDTH_BYTES : 0));

            if (~i_valid & ~clock_cycles_no_valid[15]) 
               clock_cycles_no_valid <= clock_cycles_no_valid + 1;
            else if (i_valid) 
               clock_cycles_no_valid <= 16'b0;

            if (~clock_cycles_no_valid[15]) clock_cycles_total <= clock_cycles_total + 1;

            total_ivalid <= total_ivalid + (i_valid & ~o_stall);
            total_req <= total_req + (i_valid & ~i_predicate & ~o_stall);

            if(i_stall & o_valid)
              i_stall_count <= i_stall_count + 'd1;

            if(o_stall & i_valid)
              o_stall_count <= o_stall_count + 'd1;

            if(avm_read & ~avm_waitrequest)
            begin
              avm_read_count <= avm_read_count + 'd1;
              avm_burstcount_total <= avm_burstcount_total + avm_burstcount;
            end
         end
      end
end
endgenerate

// synthesis translate_off
// Profiling data - for simulation only
reg  [31:0] bw_kernel;
reg  [31:0] bw_avalon;

// Measure Bandwidth on Avalon signals
always@(posedge clock or negedge resetn)
begin
   if (!resetn)
     bw_avalon <= 0;
   else 
     if (READ==1 && avm_readdatavalid)
       bw_avalon <= bw_avalon + MWIDTH_BYTES;
     else if (READ==0 && avm_write && ~avm_waitrequest)
       bw_avalon <= bw_avalon + MWIDTH_BYTES;
end

// Measure Bandwidth on kernel signals
always@(posedge clock or negedge resetn)
begin
   if (!resetn)
     bw_kernel <= 0;
   else if (i_valid && !o_stall && ~nop)
     bw_kernel <= bw_kernel + WIDTH_BYTES;
end
// synthesis translate_on

generate
if(PROFILE_ADDR_TOGGLE==1 && STYLE!="SIMPLE")
begin
  localparam COUNTERWIDTH=12;
  // We currently assume AWIDTH is always 32, but we really need to set this to
  // a tight lower bound to avoid wasting area here.
  logic [COUNTERWIDTH-1:0] togglerate[AWIDTH-ALIGNMENT_ABITS+1];

  acl_toggle_detect 
    #(.WIDTH(AWIDTH-ALIGNMENT_ABITS), .COUNTERWIDTH(COUNTERWIDTH)) atd (
      .clk(clock),
      .resetn(resetn),
      .valid(i_valid && ~o_stall && ~nop),
      .value({i_address >> ALIGNMENT_ABITS,{ALIGNMENT_ABITS{1'b0}}}),
      .count(togglerate));

  acl_debug_mem #(.WIDTH(COUNTERWIDTH), .SIZE(AWIDTH-ALIGNMENT_ABITS+1)) dbg_mem (
      .clk(clock),
      .resetn(resetn),
      .write(i_valid && ~o_stall && ~nop),
      .data(togglerate));

end
endgenerate

endmodule

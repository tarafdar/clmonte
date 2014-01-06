module lsu_non_aligned_read
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_stall, o_valid, o_readdata, 
    o_active, //Debugging signal
    avm_address, avm_read, avm_readdata, avm_waitrequest, avm_byteenable,
    avm_readdatavalid,
    avm_burstcount,
    flush, i_nop,

    // Profiling
    extra_unaligned_reqs,
    req_cache_hit_count
);

parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=160;   // The max number of live threads
parameter MEMORY_SIDE_MEM_LATENCY=32;    
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USECACHING=0;
parameter CACHE_SIZE_N=1024;
parameter HIGH_FMAX=1;
parameter ACL_PROFILE=0;

localparam WIDTH=WIDTH_BYTES*8;
localparam MWIDTH=MWIDTH_BYTES*8;
// -------- Interface Declarations ------------

// Standard global signals
input clk;
input clk2x;
input reset;
input flush;
input i_nop;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;

// Downstream interface
input i_stall;
output reg o_valid;
output reg [WIDTH-1:0] o_readdata;
output o_active;

// Avalon interface
output [AWIDTH-1:0] avm_address;
output avm_read;
input [MWIDTH-1:0] avm_readdata;
input avm_waitrequest;
output [MWIDTH_BYTES-1:0] avm_byteenable;
input avm_readdatavalid;
output [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

// Profiling
output logic [63:0] extra_unaligned_reqs;
output logic [63:0] req_cache_hit_count;

reg reg_lsu_i_valid;
reg [AWIDTH-1:0] reg_lsu_i_address;
reg reg_nop;
reg reg_consecutive;

wire stall_signal_directly_from_lsu;
assign o_stall = reg_lsu_i_valid & stall_signal_directly_from_lsu;

// --------------- Pipeline stage : Consecutive Address Checking  --------------------
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_lsu_i_valid <= 1'b0;
      reg_lsu_i_address <= 'x;
      reg_consecutive <= 'x;
      reg_nop <= 'x;
   end
   else
   begin
      if (~o_stall & i_valid & ~i_nop)
      begin
         reg_lsu_i_address <= i_address;
      end

      if (~o_stall)
      begin
         reg_lsu_i_valid <= i_valid;
         reg_nop <= i_nop;
         reg_consecutive <= ((reg_lsu_i_address+WIDTH_BYTES) == i_address);
      end
   end
end
// -------------------------------------------------------------------

lsu_non_aligned_read_internal #(
   .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
   .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
   .AWIDTH(AWIDTH),
   .WIDTH_BYTES(WIDTH_BYTES),
   .MWIDTH_BYTES(MWIDTH_BYTES),
   .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
   .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
   .USECACHING(USECACHING),
   .CACHE_SIZE_N(CACHE_SIZE_N),
   .HIGH_FMAX(HIGH_FMAX),
   .ACL_PROFILE(ACL_PROFILE)
) non_aligned_read (
   .clk(clk),
   .clk2x(clk2x),
   .reset(reset),
   .flush(flush),
   .o_stall(stall_signal_directly_from_lsu),
   .i_valid(reg_lsu_i_valid),
   .i_address(reg_lsu_i_address),
   .i_stall(i_stall),
   .o_valid(o_valid),
   .o_readdata(o_readdata),
   .o_active(o_active),
   .avm_address(avm_address),
   .avm_read(avm_read),
   .avm_readdata(avm_readdata),
   .avm_waitrequest(avm_waitrequest),
   .avm_byteenable(avm_byteenable),
   .avm_burstcount(avm_burstcount),
   .avm_readdatavalid(avm_readdatavalid),
   .i_nop(reg_nop),
   .consecutive(reg_consecutive),
   .extra_unaligned_reqs(extra_unaligned_reqs),
   .req_cache_hit_count(req_cache_hit_count)
);
endmodule

// Non-aligned read wrapper for LSUs
//
module lsu_non_aligned_read_internal
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_stall, o_valid, o_readdata, 
    o_active, //Debugging signal
    avm_address, avm_read, avm_readdata, avm_waitrequest, avm_byteenable,
    avm_readdatavalid,
    avm_burstcount,
    flush, i_nop,
    consecutive,

    // Profiling
    extra_unaligned_reqs,
    req_cache_hit_count
);


// Paramaters to pass down to lsu_top
//
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=160;   // The max number of live threads
parameter MEMORY_SIDE_MEM_LATENCY=32;    
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USECACHING=0;
parameter CACHE_SIZE_N=1024;
parameter HIGH_FMAX=1;

parameter TIMEOUT=8;

parameter ACL_PROFILE=0;

localparam WIDTH=WIDTH_BYTES*8;
localparam MWIDTH=MWIDTH_BYTES*8;
localparam TRACKING_FIFO_DEPTH=KERNEL_SIDE_MEM_LATENCY+1;
localparam WIDTH_ABITS=$clog2(WIDTH_BYTES);
localparam TIMEOUTBITS=$clog2(TIMEOUT);

// 
// Suppose that we vectorize 4 ways and are accessing a float4 but are only guaranteed float alignment
//
// WIDTH_BYTES=16           --> $clog2(WIDTH_BYTES) = 4
// ALIGNMENT_ABITS          --> 2
// UNALIGNED_SELECTION_BITS --> 2
//
// +----+----+----+----+----+----+
// | X  | Y  |  Z | W  | A  | B  |
// +----+----+----+----+----+----+
//  0000 0100 1000 1100  ...
//
// float4 access at 1000
//   requires two aligned access
//        0000 -> mux out Z , W
//       10000 -> mux out A , B
//
localparam UNALIGNED_SELECTION_BITS=$clog2(WIDTH_BYTES)-ALIGNMENT_ABITS;

// How much alignment are we guaranteed in terms of bits
// float -> ALIGNMENT_ABITS=2 -> 4 bytes -> 32 bits
localparam ALIGNMENT_DBITS=8*(2**ALIGNMENT_ABITS);

// -------- Interface Declarations ------------

// Standard global signals
input clk;
input clk2x;
input reset;
input flush;
input i_nop;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;

// Downstream interface
input i_stall;
output reg o_valid;
output reg [WIDTH-1:0] o_readdata;
output o_active;

// Avalon interface
output [AWIDTH-1:0] avm_address;
output avm_read;
input [MWIDTH-1:0] avm_readdata;
input avm_waitrequest;
output [MWIDTH_BYTES-1:0] avm_byteenable;
input avm_readdatavalid;
output [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

// external help for tracking addresses
input consecutive;

// Profiling
output logic [63:0] extra_unaligned_reqs;
output logic [63:0] req_cache_hit_count;

// ------- Bursting LSU instantiation ---------

wire lsu_o_stall;
wire lsu_o_valid;
wire [WIDTH-1:0] lsu_o_readdata;
wire lsu_i_valid;
wire lsu_i_stall;
wire [AWIDTH-1:0] lsu_i_address;

reg reg_lsu_i_valid;
reg [AWIDTH-1:0] reg_lsu_i_address;

wire stall_signal_directly_from_lsu;
assign lsu_o_stall = reg_lsu_i_valid & stall_signal_directly_from_lsu;

// --- Pipeline before going into the LSU ---
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_lsu_i_valid <= 1'b0;
      reg_lsu_i_address <= 'x;
   end
   else
   begin
      if (~lsu_o_stall)
      begin
         reg_lsu_i_valid <= lsu_i_valid;
         reg_lsu_i_address <= lsu_i_address;
      end
   end
end
// --- ---------------------------------- ---
lsu_bursting_read #(
   .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
   .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
   .AWIDTH(AWIDTH),
   .WIDTH_BYTES(WIDTH_BYTES),
   .MWIDTH_BYTES(MWIDTH_BYTES),
   .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
   .ALIGNMENT_ABITS(WIDTH_ABITS),
   .USECACHING(USECACHING),
   .CACHE_SIZE_N(CACHE_SIZE_N),
   .HIGH_FMAX(HIGH_FMAX),
   .FORCE_OUTPUT_REGISTER(HIGH_FMAX),
   .ACL_PROFILE(ACL_PROFILE)
) bursting_read (
   .clk(clk),
   .clk2x(clk2x),
   .reset(reset),
   .flush(flush),
   .o_stall(stall_signal_directly_from_lsu),
   .i_valid(reg_lsu_i_valid),
   .i_address(reg_lsu_i_address),
   .i_stall(lsu_i_stall),
   .o_valid(lsu_o_valid),
   .o_readdata(lsu_o_readdata),
   .o_active(o_active),
   .avm_address(avm_address),
   .avm_read(avm_read),
   .avm_readdata(avm_readdata),
   .avm_waitrequest(avm_waitrequest),
   .avm_byteenable(avm_byteenable),
   .avm_burstcount(avm_burstcount),
   .avm_readdatavalid(avm_readdatavalid),
   .req_cache_hit_count(req_cache_hit_count)
);

reg [WIDTH-1:0] reg_o_readdata;
reg reg_o_readdata_valid;
wire reg_stalled;

wire i_stall_to_oreg;
assign reg_stalled = i_stall_to_oreg & reg_o_readdata_valid;
assign lsu_i_stall = reg_stalled;
 
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_o_readdata_valid <= 1'b0;
      reg_o_readdata <= 'x;         
   end
   else
   begin
      if (~reg_stalled)
      begin
         reg_o_readdata_valid <= lsu_o_valid;
         reg_o_readdata <= lsu_o_readdata;            
      end
   end
end

// --------- Module Internal State -------------

// Timeout counter --> issue whatever stuff we have lying around
reg [TIMEOUTBITS:0] timeout;
wire timedout;
assign timedout=timeout[TIMEOUTBITS];
// The signals to the tracking fifo
wire tracking_fifo_valid_out;
wire tracking_fifo_stall_in;
wire tracking_fifo_stall_out;
wire tracking_fifo_valid_in;

// Do we want to make a request
wire request;

// Need to have space in the tracking FIFO and the LSU has to be able to accept
wire can_accept_request;
assign can_accept_request = ~(tracking_fifo_stall_out | lsu_o_stall); 

// Is the request accepted
wire request_accepted;
assign request_accepted = request & can_accept_request;

reg [AWIDTH-1:0] last_address;
reg still_need_to_issue_2nd_word;

// When do we need to issue the 2nd word ?
//  1. The previous address needed a 2nd word and the current requested address isn't 
//     consecutive with the previous
//  2. timeout
wire issue_2nd_word;
assign issue_2nd_word = (i_valid & (~consecutive | i_nop) | timedout) & still_need_to_issue_2nd_word;

// The actual requested address going into the LSU
assign lsu_i_address = issue_2nd_word ? {last_address[AWIDTH-1:WIDTH_ABITS]+1,{WIDTH_ABITS{1'b0}}} : 
                                            {i_address[AWIDTH-1:WIDTH_ABITS],{WIDTH_ABITS{1'b0}}};

// Stall out if we 
//  1. can't accept the request right now because of fifo fullness or lsu stalls
//  2. we need to issue the 2nd word from previous requests before proceeding to this one
assign o_stall = i_valid & ~issue_2nd_word & ~can_accept_request |
                 issue_2nd_word;

// Look at the lower bits to determine how much we need to shift the two aligned accesses
wire [UNALIGNED_SELECTION_BITS-1:0] unaligned_part_select;
assign unaligned_part_select = i_address[UNALIGNED_SELECTION_BITS-1+ALIGNMENT_ABITS:ALIGNMENT_ABITS];

// Is this request access already aligned .. then no need to do anything special
wire is_access_aligned;
assign is_access_aligned = ~(|unaligned_part_select);

assign request = issue_2nd_word | i_valid;
 
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      timeout <= {(TIMEOUTBITS+1){1'b0}};
   end
   else
   begin
      if (~still_need_to_issue_2nd_word | request_accepted | i_valid)
      begin
         timeout <= {(TIMEOUTBITS+1){1'b0}};
      end
      else if (~i_valid & ~timedout)
      begin
         timeout <= timeout + 1;
      end
   end
end

// Create the requests to the tracking fifo and the lsu itself
assign tracking_fifo_valid_in = request & ~lsu_o_stall;
assign lsu_i_valid = (i_valid & ~i_nop | issue_2nd_word) & ~tracking_fifo_stall_out;

always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      still_need_to_issue_2nd_word <= 1'b0;
      last_address <= 'x;
   end
   else
   begin
      still_need_to_issue_2nd_word <= issue_2nd_word ? ~request_accepted : 
                           (i_valid & ~i_nop & ~o_stall ? ~is_access_aligned : still_need_to_issue_2nd_word); 
      last_address <= (i_valid & ~i_nop & ~o_stall) ? i_address : last_address;
   end
end

wire [UNALIGNED_SELECTION_BITS-1:0] shift_to_align;
wire flush_word;
wire return_nop_val;

acl_data_fifo #(
    .DATA_WIDTH(2+UNALIGNED_SELECTION_BITS),
    .DEPTH(TRACKING_FIFO_DEPTH),
    .IMPL(HIGH_FMAX ? "sandwich" : "ram")
) tracking_fifo (
    .clock(clk),
    .resetn(~reset),
    .data_in( {i_nop & ~issue_2nd_word,issue_2nd_word,unaligned_part_select} ),
    .valid_in( tracking_fifo_valid_in ),
    .data_out( {return_nop_val,flush_word,shift_to_align} ),
    .valid_out( tracking_fifo_valid_out ),
    .stall_in( tracking_fifo_stall_in ),
    .stall_out( tracking_fifo_stall_out )
);

wire i_stall_from_pipeline_reg;
wire o_valid_int;

// the tracking fifo stalls until we output useful data
assign tracking_fifo_stall_in = ~((o_valid_int & ~i_stall_from_pipeline_reg) | (tracking_fifo_valid_out & flush_word & reg_o_readdata_valid)); 
assign i_stall_to_oreg = ~((o_valid_int & ~i_stall_from_pipeline_reg & ~return_nop_val) | (tracking_fifo_valid_out & flush_word & reg_o_readdata_valid)); 

// The output is valid whenever
//   1. there is a valid entry in the tracking fifo
//   2. the flush_word bit is not 0
//   3. the final output register has valid data 
//   4.   in the case where the access is unaligned, we need the lsu itself to also have valid data (ie. the 2nd word)
//
assign o_valid_int = tracking_fifo_valid_out & ~flush_word & ~return_nop_val & reg_o_readdata_valid & ( ~(|shift_to_align) | lsu_o_valid ) | tracking_fifo_valid_out & return_nop_val;

// Multiplexing logic
wire [2*WIDTH-1:0] output_data;

wire [WIDTH-1:0] o_readdata_int;
assign output_data = { lsu_o_readdata, reg_o_readdata };
assign o_readdata_int = output_data[ (shift_to_align*ALIGNMENT_DBITS) +: WIDTH ]; 

// --- Pipeline after the LSU ---
assign i_stall_from_pipeline_reg = o_valid & i_stall;

always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      o_valid <= 1'b0; 
      o_readdata <= 'x;
   end
   else
   begin
      if (~i_stall_from_pipeline_reg)
      begin
         o_valid <= o_valid_int;
         o_readdata <= o_readdata_int;
      end
   end
end

// Profiling
generate
if(ACL_PROFILE == 1)
begin
  always @(posedge clk or posedge reset)
  begin
    if(reset)
    begin
      extra_unaligned_reqs <= '0;
    end
    else if(flush)
    begin
      extra_unaligned_reqs <= '0;
    end
    else if(tracking_fifo_valid_in & ~tracking_fifo_stall_out & issue_2nd_word & ~i_nop)  // TODO verify
      extra_unaligned_reqs <= extra_unaligned_reqs + 'd1;
    begin
    end
  end
end
endgenerate

endmodule

module lsu_non_aligned_write
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_writedata, i_stall, o_valid, 
    o_active, //Debugging signal
    avm_address, avm_write, avm_writeack, avm_writedata, avm_byteenable, avm_waitrequest,
    avm_burstcount, i_nop
);
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=32;    // Memory latency in threads
parameter MEMORY_SIDE_MEM_LATENCY=32;    
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USE_WRITE_ACK=0;      // Wait till the write has actually made it to global memory
parameter HIGH_FMAX=1;

localparam WIDTH=8*WIDTH_BYTES;
localparam MWIDTH=8*MWIDTH_BYTES;
localparam BYTE_SELECT_BITS=$clog2(MWIDTH_BYTES);
/********
* Ports *
********/
// Standard global signals
input clk;
input clk2x;
input reset;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;
input [WIDTH-1:0] i_writedata;

// Downstream interface
input i_stall;
output o_valid;
output reg o_active;

// Avalon interface
output [AWIDTH-1:0] avm_address;
output avm_write;
input avm_writeack;
output [MWIDTH-1:0] avm_writedata;
output [MWIDTH_BYTES-1:0] avm_byteenable;
input avm_waitrequest;
output [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

input i_nop;

reg reg_lsu_i_valid;
reg [AWIDTH-1:0] reg_lsu_i_address;
reg [WIDTH-1:0] reg_lsu_i_writedata;
reg reg_nop;
reg reg_consecutive;

wire stall_signal_directly_from_lsu;
assign o_stall = reg_lsu_i_valid & stall_signal_directly_from_lsu;

// --------------- Pipeline stage : Consecutive Address Checking  --------------------
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_lsu_i_valid <= 1'b0;
      reg_lsu_i_address <= 'x;
      reg_lsu_i_writedata <= 'x;
      reg_consecutive <= 'x;
      reg_nop <= 'x;
   end
   else
   begin
      if (~o_stall & i_valid & ~i_nop)
      begin
         reg_lsu_i_address <= i_address;
      end

      if (~o_stall)
      begin
         reg_lsu_i_valid <= i_valid;
         reg_lsu_i_writedata <= i_writedata;
         reg_nop <= i_nop;
         reg_consecutive <= ((reg_lsu_i_address+WIDTH_BYTES) == i_address);
      end
   end
end
// -------------------------------------------------------------------

lsu_non_aligned_write_internal #(
   .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
   .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
   .AWIDTH(AWIDTH),
   .WIDTH_BYTES(WIDTH_BYTES),
   .MWIDTH_BYTES(MWIDTH_BYTES),
   .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
   .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
   .USE_WRITE_ACK(USE_WRITE_ACK),
   .HIGH_FMAX(HIGH_FMAX)
) non_aligned_write (
   .clk(clk),
   .clk2x(clk2x),
   .reset(reset),
   .o_stall(stall_signal_directly_from_lsu),
   .i_valid(reg_lsu_i_valid),
   .i_address(reg_lsu_i_address),
   .i_writedata(reg_lsu_i_writedata),
   .i_stall(i_stall),
   .o_valid(o_valid),
   .o_active(o_active),
   .avm_address(avm_address),
   .avm_write(avm_write),
   .avm_writeack(avm_writeack),
   .avm_writedata(avm_writedata),
   .avm_byteenable(avm_byteenable),
   .avm_burstcount(avm_burstcount),
   .avm_waitrequest(avm_waitrequest),
   .i_nop(reg_nop),
   .consecutive(reg_consecutive)
);
endmodule

//
// Non-aligned write wrapper for LSUs
//
module lsu_non_aligned_write_internal
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_writedata, i_stall, o_valid, 
    o_active, //Debugging signal
    avm_address, avm_write, avm_writeack, avm_writedata, avm_byteenable, avm_waitrequest,
    avm_burstcount,
    i_nop,
    consecutive
);

// Paramaters to pass down to lsu_top
//
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=160;   // Determines the max number of live requests.
parameter MEMORY_SIDE_MEM_LATENCY=0;   // Determines the max number of live requests.
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USECACHING=0;
parameter USE_WRITE_ACK=0;
parameter TIMEOUT=8;
parameter HIGH_FMAX=1;

localparam WIDTH=WIDTH_BYTES*8;
localparam MWIDTH=MWIDTH_BYTES*8;
localparam TRACKING_FIFO_DEPTH=KERNEL_SIDE_MEM_LATENCY+1;
localparam WIDTH_ABITS=$clog2(WIDTH_BYTES);
localparam TIMEOUTBITS=$clog2(TIMEOUT);

// 
// Suppose that we vectorize 4 ways and are accessing a float4 but are only guaranteed float alignment
//
// WIDTH_BYTES=16           --> $clog2(WIDTH_BYTES) = 4
// ALIGNMENT_ABITS          --> 2
// UNALIGNED_SELECTION_BITS --> 2
//
// +----+----+----+----+----+----+
// | X  | Y  |  Z | W  | A  | B  |
// +----+----+----+----+----+----+
//  0000 0100 1000 1100  ...
//
// float4 access at 1000
//   requires two aligned access
//        0000 -> mux out Z , W
//       10000 -> mux out A , B
//
localparam UNALIGNED_SELECTION_BITS=$clog2(WIDTH_BYTES)-ALIGNMENT_ABITS;

// How much alignment are we guaranteed in terms of bits
// float -> ALIGNMENT_ABITS=2 -> 4 bytes -> 32 bits
localparam ALIGNMENT_DBYTES=2**ALIGNMENT_ABITS;
localparam ALIGNMENT_DBITS=8*ALIGNMENT_DBYTES;

// -------- Interface Declarations ------------

// Standard global signals
input clk;
input clk2x;
input reset;
input i_nop;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;
input [WIDTH-1:0] i_writedata;

// Downstream interface
input i_stall;
output o_valid;
output o_active;

// Avalon interface
output [AWIDTH-1:0] avm_address;
output avm_write;
input avm_writeack;
output [MWIDTH-1:0] avm_writedata;
output [MWIDTH_BYTES-1:0] avm_byteenable;
input avm_waitrequest;
output [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

// help from outside to track addresses
input consecutive;

// ------- Bursting LSU instantiation ---------

wire lsu_o_stall;
wire lsu_o_valid;
wire lsu_i_valid;
wire lsu_i_stall;
wire [AWIDTH-1:0] lsu_i_address;
wire [2*WIDTH-1:0] lsu_i_writedata;
wire [2*WIDTH_BYTES-1:0] lsu_i_byte_enable;
wire [2*WIDTH-1:0] lsu_i_bit_enable;

reg [WIDTH-1:0] last_writedata;
reg [WIDTH-1:0] last_bit_enable;
reg [WIDTH_BYTES-1:0] last_byte_enable;

reg reg_lsu_i_valid;
reg [AWIDTH-1:0] reg_lsu_i_address;
reg [WIDTH-1:0] reg_lsu_i_writedata;
reg [WIDTH-1:0] reg_lsu_i_bit_enable;
reg [WIDTH_BYTES-1:0] reg_lsu_i_byte_enable;

wire stall_signal_directly_from_lsu;
assign lsu_o_stall = reg_lsu_i_valid & stall_signal_directly_from_lsu;

// --- Pipeline before going into the LSU ---
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_lsu_i_valid <= 1'b0;
      reg_lsu_i_address <= 'x;
      reg_lsu_i_writedata <= 'x;
      reg_lsu_i_bit_enable <= 'x;
      reg_lsu_i_byte_enable <= 'x;
   end
   else
   begin
      if (~lsu_o_stall)
      begin
         reg_lsu_i_valid <= lsu_i_valid;
         reg_lsu_i_address <= lsu_i_address;
         reg_lsu_i_writedata <= lsu_i_writedata[WIDTH-1:0];
         reg_lsu_i_bit_enable <= lsu_i_bit_enable[WIDTH-1:0];
         reg_lsu_i_byte_enable <= lsu_i_byte_enable[WIDTH_BYTES-1:0];
      end
   end
end
// --- ---------------------------------- ---

lsu_bursting_write #(
   .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
   .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
   .AWIDTH(AWIDTH),
   .WIDTH_BYTES(WIDTH_BYTES),
   .MWIDTH_BYTES(MWIDTH_BYTES),
   .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
   .ALIGNMENT_ABITS(WIDTH_ABITS),
   .USE_WRITE_ACK(USE_WRITE_ACK),
   .HIGH_FMAX(HIGH_FMAX)
) bursting_write (
   .clk(clk),
   .clk2x(clk2x),
   .reset(reset),
   .o_stall(stall_signal_directly_from_lsu),
   .i_valid(reg_lsu_i_valid),
   .i_address(reg_lsu_i_address),
   .i_writedata(reg_lsu_i_writedata),
   .i_stall(lsu_i_stall),
   .o_valid(lsu_o_valid),
   .o_active(o_active),
   .word_byte_enable(reg_lsu_i_byte_enable),
   .word_bit_enable(reg_lsu_i_bit_enable),
   .avm_address(avm_address),
   .avm_write(avm_write),
   .avm_writeack(avm_writeack),
   .avm_writedata(avm_writedata),
   .avm_byteenable(avm_byteenable),
   .avm_burstcount(avm_burstcount),
   .avm_waitrequest(avm_waitrequest)
);

reg reg_o_writedata;
wire reg_stalled;

wire i_stall_to_oreg;
assign reg_stalled = i_stall_to_oreg & reg_o_writedata;
assign lsu_i_stall = reg_stalled;
 
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_o_writedata <= 1'b0;
   end
   else
   begin
      if (~reg_stalled)
      begin
         reg_o_writedata <= lsu_o_valid;
      end
   end
end

// --------- Module Internal State -------------

// Timeout counter --> issue whatever stuff we have lying around
reg [TIMEOUTBITS:0] timeout;
wire timedout;
assign timedout=timeout[TIMEOUTBITS];
// The signals to the tracking fifo
wire tracking_fifo_valid_out;
wire tracking_fifo_stall_in;
wire tracking_fifo_stall_out;
wire tracking_fifo_valid_in;

// Do we want to make a request
wire request;

// Need to have space in the tracking FIFO and the LSU has to be able to accept
wire can_accept_request;
assign can_accept_request = ~(tracking_fifo_stall_out | lsu_o_stall); 

// Is the request accepted
wire request_accepted;
assign request_accepted = request & can_accept_request;

reg [AWIDTH-1:0] last_address;
reg still_need_to_issue_2nd_word;

// When do we need to issue the 2nd word ?
//  1. The previous address needed a 2nd word and the current requested address isn't 
//     consecutive with the previous
//  2. timeout
wire issue_2nd_word;
assign issue_2nd_word = (i_valid & (~consecutive | i_nop) | timedout) & still_need_to_issue_2nd_word;

// The actual requested address going into the LSU
assign lsu_i_address = issue_2nd_word ? {last_address[AWIDTH-1:WIDTH_ABITS]+1,{WIDTH_ABITS{1'b0}}} : 
                                            {i_address[AWIDTH-1:WIDTH_ABITS],{WIDTH_ABITS{1'b0}}};

// Look at the lower bits to determine how much we need to shift the two aligned accesses
wire [UNALIGNED_SELECTION_BITS-1:0] unaligned_part_select;
assign unaligned_part_select = i_address[UNALIGNED_SELECTION_BITS-1+ALIGNMENT_ABITS:ALIGNMENT_ABITS];

// The actual data to be written and corresponding byte/bit enables
assign lsu_i_bit_enable   = last_bit_enable & {WIDTH{still_need_to_issue_2nd_word}} | 
                            ({{WIDTH{1'b0}},{WIDTH{~issue_2nd_word}}} << (unaligned_part_select * ALIGNMENT_DBITS));

assign lsu_i_byte_enable  = last_byte_enable & {WIDTH_BYTES{still_need_to_issue_2nd_word}} | 
                            ({{WIDTH_BYTES{1'b0}},{WIDTH_BYTES{~issue_2nd_word}}} << (unaligned_part_select * ALIGNMENT_DBYTES));

assign lsu_i_writedata    = (({{WIDTH{1'b0}}, i_writedata & {WIDTH{~issue_2nd_word}}}) << (unaligned_part_select * ALIGNMENT_DBITS)) | 
                            last_writedata & {WIDTH{still_need_to_issue_2nd_word}};

// Stall out if we 
//  1. can't accept the request right now because of fifo fullness or lsu stalls
//  2. we need to issue the 2nd word from previous requests before proceeding to this one
assign o_stall = i_valid & ~issue_2nd_word & ~can_accept_request |
                 issue_2nd_word;

// Is this request access already aligned .. then no need to do anything special
wire is_access_aligned;
assign is_access_aligned = ~(|unaligned_part_select);

assign request = issue_2nd_word | i_valid;
 
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      timeout <= {(TIMEOUTBITS+1){1'b0}};
   end
   else
   begin
      if (~still_need_to_issue_2nd_word | request_accepted | i_valid)
      begin
         timeout <= {(TIMEOUTBITS+1){1'b0}};
      end
      else if (~i_valid & ~timedout)
      begin
         timeout <= timeout + 1;
      end
   end
end

// Create the requests to the tracking fifo and the lsu itself
assign tracking_fifo_valid_in = request & ~lsu_o_stall;
assign lsu_i_valid = (i_valid & ~i_nop | issue_2nd_word) & ~tracking_fifo_stall_out;

always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      still_need_to_issue_2nd_word <= 1'b0;
      last_address <= 'x;
      last_writedata <= 'x;
      last_bit_enable <= {WIDTH{1'b0}};
      last_byte_enable <= {WIDTH_BYTES{1'b0}};
   end
   else
   begin
      still_need_to_issue_2nd_word <= issue_2nd_word ? ~request_accepted : 
                           (i_valid & ~i_nop & ~o_stall ? ~is_access_aligned : still_need_to_issue_2nd_word); 
      last_address <= (i_valid & ~i_nop & ~o_stall) ? i_address : last_address;
      last_writedata <= (i_valid & ~i_nop & ~o_stall) ? lsu_i_writedata[2*WIDTH-1:WIDTH] : last_writedata;
      last_bit_enable <= (i_valid & ~i_nop & ~o_stall) ? lsu_i_bit_enable[2*WIDTH-1:WIDTH] : last_bit_enable;
      last_byte_enable <= (i_valid & ~i_nop & ~o_stall) ? lsu_i_byte_enable[2*WIDTH_BYTES-1:WIDTH_BYTES] : last_byte_enable;
   end
end

wire no_shift_needed;
wire flush_word;
wire return_nop_val;

// Use the LL fifo for small depths
generate
if(TRACKING_FIFO_DEPTH <= 8)
begin
   wire tracking_fifo_full;  
   wire tracking_fifo_empty;
   acl_ll_fifo #(
    .WIDTH(3),
    .DEPTH(TRACKING_FIFO_DEPTH)
   ) tracking_fifo (
      .clk(clk),
      .reset(reset),
      .data_in( {i_nop & ~issue_2nd_word,issue_2nd_word,~(|unaligned_part_select)} ),
      .write( tracking_fifo_valid_in & ~tracking_fifo_stall_out ),
      .data_out( {return_nop_val,flush_word,no_shift_needed} ),
      .read(tracking_fifo_valid_out & ~tracking_fifo_stall_in),
      .full(tracking_fifo_full),
      .empty(tracking_fifo_empty)
   );
   assign tracking_fifo_valid_out = ~tracking_fifo_empty;
   assign tracking_fifo_stall_out = tracking_fifo_full;
end
else
begin
   acl_data_fifo #(
      .DATA_WIDTH(3),
      .DEPTH(TRACKING_FIFO_DEPTH),
      .IMPL(HIGH_FMAX ? "sandwich" : "ram")
   ) tracking_fifo (
      .clock(clk),
      .resetn(~reset),
      .data_in( {i_nop & ~issue_2nd_word,issue_2nd_word,~(|unaligned_part_select)} ),
      .valid_in( tracking_fifo_valid_in ),
      .data_out( {return_nop_val,flush_word,no_shift_needed} ),
      .valid_out( tracking_fifo_valid_out ),
      .stall_in( tracking_fifo_stall_in ),
      .stall_out( tracking_fifo_stall_out )
   );
end
endgenerate

// the tracking fifo stalls until we output useful data
assign tracking_fifo_stall_in = ~((o_valid & ~i_stall) | (tracking_fifo_valid_out & flush_word & reg_o_writedata)); 
assign i_stall_to_oreg = ~((o_valid & ~i_stall & ~return_nop_val) | (tracking_fifo_valid_out & flush_word & reg_o_writedata)); 

// The output is valid whenever
//   1. there is a valid entry in the tracking fifo
//   2. the flush_word bit is not 0
//   3. the final output register has valid data 
//   4.   in the case where the access is unaligned, we need the lsu itself to also have valid data (ie. the 2nd word)
//
assign o_valid = tracking_fifo_valid_out & ~flush_word & ~return_nop_val & reg_o_writedata & ( no_shift_needed | lsu_o_valid ) | tracking_fifo_valid_out & return_nop_val;

endmodule


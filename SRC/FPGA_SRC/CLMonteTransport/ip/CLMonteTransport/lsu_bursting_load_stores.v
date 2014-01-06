// Bursting read LSU
// I promise I'll add more top level comments later in the week. Right now, I need to check this in
//
module lsu_bursting_read
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_stall, o_valid, o_readdata, 
    o_active, //Debugging signal
    avm_address, avm_read, avm_readdata, avm_waitrequest, avm_byteenable,
    avm_readdatavalid,
    avm_burstcount,
    flush,

    // Profiling
    req_cache_hit_count
);

/*************
* Parameters *
*************/
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=160;   // The max number of live threads
parameter MEMORY_SIDE_MEM_LATENCY=32;    // Memory latency in cycles
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USECACHING=0;
parameter MAX_THREADS=64;       // Maximum # of threads to group into a burst request
parameter HIGH_FMAX=1;
parameter FORCE_OUTPUT_REGISTER=0;

parameter CACHE_SIZE_N=1024;

parameter ACL_PROFILE = 0;

// Derived cache parameters
localparam CACHE_SIZE_LOG2N=$clog2(CACHE_SIZE_N);

// Derived parameters
localparam MAX_BURST=2**(BURSTCOUNT_WIDTH-1);

localparam WIDTH=8*WIDTH_BYTES;
localparam MWIDTH=8*MWIDTH_BYTES;
localparam BYTE_SELECT_BITS=$clog2(MWIDTH_BYTES);
localparam SEGMENT_SELECT_BITS=BYTE_SELECT_BITS-ALIGNMENT_ABITS;
localparam PAGE_SELECT_BITS=AWIDTH-BYTE_SELECT_BITS;

// Constants
localparam REQUEST_FIFO_DEPTH=2*KERNEL_SIDE_MEM_LATENCY;

/********
* Ports *
********/
// Standard global signals
input clk;
input clk2x;
input reset;
input flush;

// Upstream interface
output o_stall;
input i_valid;
input [AWIDTH-1:0] i_address;

// Downstream interface
input i_stall;
output o_valid;
output [WIDTH-1:0] o_readdata;
output o_active;

// Avalon interface
output reg [AWIDTH-1:0] avm_address;
output reg avm_read;
input [MWIDTH-1:0] avm_readdata;
input avm_waitrequest;
output [MWIDTH_BYTES-1:0] avm_byteenable;
input avm_readdatavalid;
output reg [BURSTCOUNT_WIDTH-1:0] avm_burstcount;

// Profiling
output logic [63:0] req_cache_hit_count;

/***************
* Architecture *
***************/
wire [PAGE_SELECT_BITS-1:0] page_addr;
wire [SEGMENT_SELECT_BITS:0] rf_data_in;
wire [BYTE_SELECT_BITS-1:0] byte_addr;
wire next_new_page;
wire c_stall;
wire c_new_page;
wire [PAGE_SELECT_BITS-1:0] c_req_addr;
wire c_req_valid;
wire rf_full;
wire rf_valid;
wire [SEGMENT_SELECT_BITS:0] rf_data;
wire rf_next_valid;
wire [SEGMENT_SELECT_BITS:0] rf_next_data;
wire rf_stall_in;
wire rm_stall;
wire rm_valid ;
wire [MWIDTH-1:0] rm_data;
wire rm_stall_in;
wire [BURSTCOUNT_WIDTH-1:0] burstcount;

wire [AWIDTH-1:0] i_address_from_cache;
wire o_stall_to_cache;
wire i_valid_from_cache;

reg in_my_cache;

wire i_stall_from_oreg;

generate if(USECACHING)
begin
   reg [AWIDTH-1:0] i_address_reg;
   reg i_address_valid;

   reg [AWIDTH-1:0] i_address_reg_d2;
   reg i_address_valid_d2;

   wire o_stall_reg2;

   reg [AWIDTH-1-BYTE_SELECT_BITS-CACHE_SIZE_LOG2N:0] cached_base_address [0:CACHE_SIZE_N-1];
   reg [0:CACHE_SIZE_N-1] cached_base_address_valid;

   assign i_address_from_cache = i_address_reg_d2;
   assign i_valid_from_cache = i_address_valid_d2;
   assign o_stall = i_address_valid && o_stall_reg2;
   assign o_stall_reg2 = i_address_valid_d2 && o_stall_to_cache;

   wire [CACHE_SIZE_LOG2N-1:0] cache_index;
   assign cache_index = i_address_reg[CACHE_SIZE_LOG2N-1+BYTE_SELECT_BITS:BYTE_SELECT_BITS];

   always@(posedge clk or posedge reset)
   begin
      if(reset)
      begin
         i_address_reg <= 'x;       // don't need to reset
         i_address_reg_d2 <= 'x;    // don't need to reset
         in_my_cache <= 'x;         // don't need to reset
         // note: cached_base_address is not reset because it stops the RAM
         // from being inferred and that leads the creation of a huge number
         // of registers

         i_address_valid <= 1'b0;
         i_address_valid_d2 <= 1'b0;
         cached_base_address_valid <= {CACHE_SIZE_N{1'b0}};
      end
      else
      begin
         // Input Register for addresses
         if (!o_stall) 
         begin
             i_address_reg <= i_address;
             i_address_valid <= i_valid;
         end
          
         if (!o_stall_reg2)
         begin
             i_address_reg_d2 <= i_address_reg;
             i_address_valid_d2 <= i_address_valid;
             // Should infer a RAM read port here
             in_my_cache <= cached_base_address_valid[cache_index] && (cached_base_address[cache_index] == i_address_reg[AWIDTH-1:BYTE_SELECT_BITS+CACHE_SIZE_LOG2N]);
         end

         if (flush)
         begin
             cached_base_address_valid <= {CACHE_SIZE_N{1'b0}};
         end
         else if (!o_stall && i_address_valid)
         begin
             // Should infer a RAM write port here
             cached_base_address[cache_index] <= i_address_reg[AWIDTH-1:BYTE_SELECT_BITS+CACHE_SIZE_LOG2N];
             cached_base_address_valid[cache_index] <= 1'b1;
         end
      end
   end

end
else
begin
   assign i_address_from_cache = i_address;
   assign o_stall = o_stall_to_cache;
   assign i_valid_from_cache = i_valid;
end
endgenerate

// --------------- Pipeline stage : Burst Checking -------------------
reg i_valid_from_burst_check_int;
reg [AWIDTH-1:0] i_address_from_burst_check_int;
reg reg_common_burst_int;
reg pipe_in_my_cache_int;

wire i_valid_from_burst_check;
wire [AWIDTH-1:0] i_address_from_burst_check;
wire reg_common_burst;
wire pipe_in_my_cache;

reg [AWIDTH-1:0] last_non_cached_address;

wire o_stall_to_burst_check;
wire staging_reg_stall;
assign o_stall_to_cache = i_valid_from_burst_check_int & staging_reg_stall;

// --- Pipeline before going into the LSU ---
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      i_valid_from_burst_check_int <= 1'b0;
      i_address_from_burst_check_int <= 'x;
      reg_common_burst_int <= 'x;
      pipe_in_my_cache_int <= 1'b0;
      last_non_cached_address <= 'x;
   end
   else
   begin
      if (USECACHING ? (~o_stall_to_cache & i_valid_from_cache & ~in_my_cache) : (~o_stall_to_cache & i_valid_from_cache))
      begin
         last_non_cached_address <= i_address_from_cache;
      end 

      if (~o_stall_to_cache)
      begin
         i_address_from_burst_check_int <= i_address_from_cache;
         i_valid_from_burst_check_int <= i_valid_from_cache;
         reg_common_burst_int <= (last_non_cached_address[AWIDTH-1:BYTE_SELECT_BITS+BURSTCOUNT_WIDTH-1] == i_address_from_cache[AWIDTH-1:BYTE_SELECT_BITS+BURSTCOUNT_WIDTH-1]);
         pipe_in_my_cache_int <= in_my_cache;
      end
   end
end

acl_staging_reg #(
  .WIDTH(AWIDTH+2)
) staging_reg (
  .clk(clk), 
  .reset(reset), 
  .i_data({i_address_from_burst_check_int,reg_common_burst_int,pipe_in_my_cache_int}), 
  .i_valid(i_valid_from_burst_check_int), 
  .o_stall(staging_reg_stall), 
  .o_data({i_address_from_burst_check,reg_common_burst,pipe_in_my_cache}), 
  .o_valid(i_valid_from_burst_check), 
  .i_stall(o_stall_to_burst_check)
);

generate
if(USECACHING == 1 && ACL_PROFILE == 1)
begin
  always @(posedge clk or posedge reset)
  begin
    if(reset)
    begin
      req_cache_hit_count <= '0;
    end
    else if(flush)
    begin
      req_cache_hit_count <= '0;
    end
    else if(i_valid_from_burst_check & ~o_stall_to_burst_check)
    begin
      req_cache_hit_count <= req_cache_hit_count + pipe_in_my_cache;
    end
  end
end
else
begin
  assign req_cache_hit_count = 'x;
end
endgenerate

// -------------------------------------------------------------------

// Coalescer - Groups subsequent requests together if they are compatible and 
// the avalon bus is stalled.
assign page_addr = i_address_from_burst_check[AWIDTH-1:BYTE_SELECT_BITS];

wire burst_coalescer_valid_in;
generate if(USECACHING)
begin
assign burst_coalescer_valid_in = !pipe_in_my_cache && i_valid_from_burst_check && !rf_full;
end
else
begin
assign burst_coalescer_valid_in = i_valid_from_burst_check && !rf_full;
end
endgenerate

bursting_coalescer #(
    .PAGE_ADDR_WIDTH(PAGE_SELECT_BITS),
    .TIMEOUT(16),
    .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
    .MAXBURSTCOUNT(MAX_BURST),
    .MAX_THREADS(MAX_THREADS),
    .USECACHING(USECACHING)
) coalescer (
    .clk(clk),
    .reset(reset),
    .i_page_addr(page_addr),
    .i_valid(burst_coalescer_valid_in),
    .o_stall(c_stall),
    .o_new_page(c_new_page),
    .o_req_addr(c_req_addr),
    .o_req_valid(c_req_valid),
    .i_stall(rm_stall),
    .o_burstcount(burstcount),
    .common_burst(reg_common_burst),
    .i_input_accepted_from_wrapper_lsu(i_valid_from_burst_check && !o_stall_to_burst_check),
    .i_reset_timeout(i_valid_from_burst_check)
    //.i_reset_timeout(i_valid_from_burst_check & !rf_full) -- Consider guarding with rf_full, but I don't think it's necessary
);

// Response FIFO - Buffers the requests so they can be extracted from the 
// wider memory bus.  Stores the segment address to extract the requested
// word from the response data, and a bit indicating if the request comes
// from a new page.
wire [ SEGMENT_SELECT_BITS + 1 + 1 + CACHE_SIZE_LOG2N - 1 : 0 ] cache_i_data;
 
generate if(SEGMENT_SELECT_BITS > 0)
begin
   wire [SEGMENT_SELECT_BITS-1:0] segment_addr;
   assign segment_addr = i_address_from_burst_check[BYTE_SELECT_BITS-1:ALIGNMENT_ABITS];
   assign rf_data_in = {c_new_page, segment_addr};
   assign byte_addr = (rf_data[SEGMENT_SELECT_BITS-1:0] << ALIGNMENT_ABITS);
   assign cache_i_data = {pipe_in_my_cache,i_address_from_burst_check[CACHE_SIZE_LOG2N+BYTE_SELECT_BITS-1:BYTE_SELECT_BITS],rf_data_in[SEGMENT_SELECT_BITS]&!pipe_in_my_cache,rf_data_in[SEGMENT_SELECT_BITS-1:0]}; 
end
else
begin
   assign rf_data_in = c_new_page;
   assign byte_addr = {BYTE_SELECT_BITS{1'b0}};
   assign cache_i_data = {pipe_in_my_cache,i_address_from_burst_check[CACHE_SIZE_LOG2N+BYTE_SELECT_BITS-1:BYTE_SELECT_BITS],rf_data_in[SEGMENT_SELECT_BITS]&!pipe_in_my_cache}; 
end
endgenerate

wire use_data_from_cache;
wire [CACHE_SIZE_LOG2N-1:0] data_from_cache_index;
wire use_data_from_cache_next;

generate if(USECACHING)
begin
wire [CACHE_SIZE_LOG2N-1:0] dummy_next_data;
lookahead_fifo #(
    .WIDTH( SEGMENT_SELECT_BITS + 1 + 1 + CACHE_SIZE_LOG2N ),
    .DEPTH( REQUEST_FIFO_DEPTH )
) request_fifo (
    .clk(clk), 
    .reset(reset),
    .i_data(cache_i_data),
    .i_valid( i_valid_from_burst_check && !c_stall ),
    .o_full(rf_full),
    .i_stall(rf_stall_in),
    .o_valid(rf_valid),
    .o_data({use_data_from_cache,data_from_cache_index,rf_data}),
    .o_next_valid(rf_next_valid),
    .o_next_data({use_data_from_cache_next,dummy_next_data,rf_next_data})
);
end
else
begin
lookahead_fifo #(
    .WIDTH( SEGMENT_SELECT_BITS + 1 ),
    .DEPTH( REQUEST_FIFO_DEPTH )
) request_fifo (
    .clk(clk), 
    .reset(reset),
    .i_data( rf_data_in ),
    .i_valid( i_valid_from_burst_check && !c_stall ),
    .o_full(rf_full),
    .i_stall(rf_stall_in),
    .o_valid(rf_valid),
    .o_data(rf_data),
    .o_next_valid(rf_next_valid),
    .o_next_data(rf_next_data)
);
end
endgenerate

// Read master - Handles pipelined read transactions through MM-Avalon.
lsu_pipelined_read #(
    .AWIDTH( AWIDTH ),
    .WIDTH_BYTES( MWIDTH_BYTES ),
    .MWIDTH_BYTES( MWIDTH_BYTES ),
    .ALIGNMENT_ABITS( BYTE_SELECT_BITS ),
    .KERNEL_SIDE_MEM_LATENCY( MEMORY_SIDE_MEM_LATENCY ),
    .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
    .USEBURST( 1 ),
    .USEINPUTFIFO( 1 ),
    .INPUTFIFOSIZE( KERNEL_SIDE_MEM_LATENCY ), // Make sure that we have enough capacity
    .PIPELINE_INPUT( 1 ),
    .SUPERPIPELINE( 1 ),
    .HIGH_FMAX(HIGH_FMAX)
) read_master (
    .clk(clk),
    .reset(reset),
    .o_stall(rm_stall),
    .i_valid(c_req_valid),
    .i_address({c_req_addr, {BYTE_SELECT_BITS{1'b0}}}),
    .i_burstcount(burstcount),
    .i_stall(rm_stall_in),
    .o_valid(rm_valid),
    .o_active(o_active),
    .o_readdata(rm_data),
    .avm_address(avm_address),
    .avm_read(avm_read),
    .avm_readdata(avm_readdata),
    .avm_waitrequest(avm_waitrequest),
    .avm_byteenable(avm_byteenable),
    .avm_readdatavalid(avm_readdatavalid),
    .avm_burstcount(avm_burstcount)
);

generate if(USECACHING)
begin
   assign o_stall_to_burst_check = !pipe_in_my_cache && c_stall || rf_full;

   reg [MWIDTH-1:0] data_cache[0:CACHE_SIZE_N-1];
   always@(posedge clk)
   begin
      if(rm_valid && rf_valid && !use_data_from_cache)
      begin
         data_cache[data_from_cache_index] <= rm_data;
      end
   end

   wire [MWIDTH-1:0] reused_data_value;
   assign reused_data_value = data_cache[data_from_cache_index];


   // Stalling the request FIFO is a bit different now since
   // we don't care about rm data in the cace when we're reading
   // from the cache
   assign rf_stall_in = i_stall_from_oreg || (!rm_valid && !(rf_valid && use_data_from_cache));


   // This is crazy complicated now -- yikes
   assign rm_stall_in = !((!rf_valid || !rf_stall_in) && rf_next_valid && next_new_page && rm_valid); 

   //
   // Put in an output register since there is alot of data from the cache and muxing
   acl_data_fifo #(
     .DATA_WIDTH(WIDTH),
     .DEPTH(2),
     .IMPL(HIGH_FMAX && !FORCE_OUTPUT_REGISTER ? "passthrough" : "ll_reg")
   )
   lsu_oreg_fifo (
     .clock(clk),
     .resetn(~reset),
     .data_in(use_data_from_cache ? reused_data_value[8*byte_addr +: WIDTH] : rm_data[8*byte_addr +: WIDTH]),
     .data_out(o_readdata),
     .valid_in(rf_valid && (use_data_from_cache || rm_valid)),
     .valid_out(o_valid),
     .stall_in(i_stall),
     .stall_out(i_stall_from_oreg)
   );
end
else
begin
   assign rf_stall_in = i_stall_from_oreg || !rm_valid;
   assign o_stall_to_burst_check = c_stall || rf_full;
   assign rm_stall_in = (!next_new_page && rf_next_valid) || rf_stall_in;

   // Output register and cut stall path
   acl_data_fifo #(
     .DATA_WIDTH(WIDTH),
     .DEPTH(2),
     .IMPL(HIGH_FMAX ? "passthrough" : "ll_reg")
   )
   lsu_oreg_fifo (
     .clock(clk),
     .resetn(~reset),
     .data_in(rm_data[8*byte_addr +: WIDTH]),
     .data_out(o_readdata),
     .valid_in(rm_valid && rf_valid),
     .valid_out(o_valid),
     .stall_in(i_stall),
     .stall_out(i_stall_from_oreg)
   );
end
endgenerate

// Control logic
// Highest bit of rf_next_data indicates whether this is a new avalon request
// (new page) or was coalesced into the previous request.
assign next_new_page = rf_next_data[SEGMENT_SELECT_BITS];

endmodule

module lsu_bursting_write
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_writedata, i_stall, o_valid, 
    o_active, //Debugging signal
    avm_address, avm_write, avm_writeack, avm_writedata, avm_byteenable, avm_waitrequest,
    avm_burstcount,
    word_byte_enable,
    word_bit_enable
);
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=32;    // The max number of live threads
parameter MEMORY_SIDE_MEM_LATENCY=32;    // The latency to get to the iface (no response needed from DDR, we generate writeack right before the iface).
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

// Byte enable control 
input [WIDTH_BYTES-1:0] word_byte_enable;
input [WIDTH-1:0] word_bit_enable;

reg reg_lsu_i_valid;
reg [AWIDTH-1:0] reg_lsu_i_address;
reg [WIDTH-1:0] reg_lsu_i_writedata;
reg [WIDTH-1:0] reg_lsu_i_bit_enable;
reg [WIDTH_BYTES-1:0] reg_lsu_i_byte_enable;
reg reg_common_burst;

wire stall_signal_directly_from_lsu;
assign o_stall = reg_lsu_i_valid & stall_signal_directly_from_lsu;

// --------------- Pipeline stage : Burst Checking  --------------------
always@(posedge clk or posedge reset)
begin
   if (reset)
   begin
      reg_lsu_i_valid <= 1'b0;
      reg_lsu_i_address <= 'x;
      reg_lsu_i_writedata <= 'x;
      reg_lsu_i_bit_enable <= 'x;
      reg_lsu_i_byte_enable <= 'x;
      reg_common_burst <= 'x;
   end
   else
   begin
      if (~o_stall & i_valid)
      begin
         reg_lsu_i_address <= i_address;
      end

      if (~o_stall)
      begin
         reg_lsu_i_valid <= i_valid;
         reg_lsu_i_writedata <= i_writedata;
         reg_lsu_i_bit_enable <= word_bit_enable;
         reg_lsu_i_byte_enable <= word_byte_enable;
         reg_common_burst <= (reg_lsu_i_address[AWIDTH-1:BYTE_SELECT_BITS+BURSTCOUNT_WIDTH-1] == i_address[AWIDTH-1:BYTE_SELECT_BITS+BURSTCOUNT_WIDTH-1]);
      end
   end
end
// -------------------------------------------------------------------

lsu_bursting_write_internal #(
   .KERNEL_SIDE_MEM_LATENCY(KERNEL_SIDE_MEM_LATENCY),
   .MEMORY_SIDE_MEM_LATENCY(MEMORY_SIDE_MEM_LATENCY),
   .AWIDTH(AWIDTH),
   .WIDTH_BYTES(WIDTH_BYTES),
   .MWIDTH_BYTES(MWIDTH_BYTES),
   .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
   .ALIGNMENT_ABITS(ALIGNMENT_ABITS),
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
   .i_stall(i_stall),
   .o_valid(o_valid),
   .o_active(o_active),
   .word_byte_enable(reg_lsu_i_byte_enable),
   .word_bit_enable(reg_lsu_i_bit_enable),
   .avm_address(avm_address),
   .avm_write(avm_write),
   .avm_writeack(avm_writeack),
   .avm_writedata(avm_writedata),
   .avm_byteenable(avm_byteenable),
   .avm_burstcount(avm_burstcount),
   .avm_waitrequest(avm_waitrequest),
   .common_burst(reg_common_burst)
);
endmodule

//
// Burst coalesced write module
// Again, top level comments later
//
module lsu_bursting_write_internal
(
    clk, clk2x, reset, o_stall, i_valid, i_address, i_writedata, i_stall, o_valid, 
    o_active, //Debugging signal
    avm_address, avm_write, avm_writeack, avm_writedata, avm_byteenable, avm_waitrequest,
    avm_burstcount,
    word_byte_enable,
    word_bit_enable,
    common_burst
);

/*************
* Parameters *
*************/
parameter AWIDTH=32;            // Address width (32-bits for Avalon)
parameter WIDTH_BYTES=4;        // Width of the memory access (bytes)
parameter MWIDTH_BYTES=32;      // Width of the global memory bus (bytes)
parameter ALIGNMENT_ABITS=2;    // Request address alignment (address bits)
parameter KERNEL_SIDE_MEM_LATENCY=32;    // Memory latency in cycles
parameter MEMORY_SIDE_MEM_LATENCY=32;    // Memory latency in cycles
parameter BURSTCOUNT_WIDTH=6;   // Size of Avalon burst count port
parameter USE_WRITE_ACK=0;      // Wait till the write has actually made it to global memory
parameter HIGH_FMAX=1;

// WARNING: Kernels will hang if InstrDataDep claims that a store
// has more capacity than this number
//
parameter MAX_CAPACITY=128;     // Must be a power-of-2 to keep things simple

// Derived parameters
localparam MAX_BURST=2**(BURSTCOUNT_WIDTH-1);
//
// Notice that in the non write ack case, the number of threads seems to be twice the sensible number
// This is because MAX_THREADS is usually the limiter on the counter width. Once one request is assemembled,
// we want to be able to start piplining another burst. Thus the factor of 2.
// The MEMORY_SIDE_MEM_LATENCY will further increase this depth if the compiler
// thinks the lsu will see a lot of contention on the Avalon side.
//
localparam __WRITE_FIFO_DEPTH = (WIDTH_BYTES==MWIDTH_BYTES) ? 3*MAX_BURST : 2*MAX_BURST;
// No reason this should need more than max MLAB depth
localparam _WRITE_FIFO_DEPTH = ( __WRITE_FIFO_DEPTH < 64 ) ? __WRITE_FIFO_DEPTH : 64;
// Need at least 4 to account for fifo push-to-pop latency
localparam WRITE_FIFO_DEPTH = ( _WRITE_FIFO_DEPTH > 8 ) ? _WRITE_FIFO_DEPTH : 8;

// If writeack, make this equal to 
localparam MAX_THREADS=(USE_WRITE_ACK ? KERNEL_SIDE_MEM_LATENCY : (2*MWIDTH_BYTES/WIDTH_BYTES*MAX_BURST)); // Maximum # of threads to group into a burst request
//
localparam WIDTH=8*WIDTH_BYTES;
localparam MWIDTH=8*MWIDTH_BYTES;
localparam BYTE_SELECT_BITS=$clog2(MWIDTH_BYTES);
localparam SEGMENT_SELECT_BITS=BYTE_SELECT_BITS-ALIGNMENT_ABITS;
localparam PAGE_SELECT_BITS=AWIDTH-BYTE_SELECT_BITS;
localparam NUM_SEGMENTS=2**SEGMENT_SELECT_BITS;
localparam SEGMENT_WIDTH=8*(2**ALIGNMENT_ABITS);
localparam SEGMENT_WIDTH_BYTES=(2**ALIGNMENT_ABITS);

// Constants
localparam COUNTER_WIDTH=(($clog2(MAX_THREADS)+1 < $clog2(MAX_CAPACITY+1)) ? 
                          $clog2(MAX_CAPACITY+1) : ($clog2(MAX_THREADS)+1)); // Determines the max writes 'in-flight'

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

// Byte enable control 
input [WIDTH_BYTES-1:0] word_byte_enable;
input [WIDTH-1:0] word_bit_enable;

// Help from outside
input common_burst;

/***************
* Architecture *
***************/
wire input_accepted;
wire output_acknowledged;
wire write_accepted;
wire [PAGE_SELECT_BITS-1:0] page_addr;
wire c_new_page;
wire c_page_done;
wire [PAGE_SELECT_BITS-1:0] c_req_addr;
wire c_req_valid;
wire c_stall;
reg [COUNTER_WIDTH-1:0] occ_counter;
reg [COUNTER_WIDTH-1:0] inflight_stores_counter;

// Replicated version of the occ and stores counters that decrement instead of increment
// This allows me to check the topmost bit to determine if the counter is non-empty
reg [COUNTER_WIDTH-1:0] occ_counter_neg;
reg [COUNTER_WIDTH-1:0] inflight_stores_counter_neg;

reg [COUNTER_WIDTH-1:0] ack_counter;
reg [COUNTER_WIDTH-1:0] next_counter;
reg [MWIDTH-1:0] wm_writedata;
reg [MWIDTH_BYTES-1:0] wm_byteenable;
reg [MWIDTH-1:0] wm_wide_wdata;
reg [MWIDTH_BYTES-1:0] wm_wide_be;
reg [MWIDTH-1:0] wm_wide_bite;

wire w_fifo_full;
wire [BURSTCOUNT_WIDTH-1:0] c_burstcount;
// Track the current item in the write burst since we issue c_burstcount burst reqs
reg [BURSTCOUNT_WIDTH-1:0] burstcounter;

// The address components
assign page_addr = i_address[AWIDTH-1:BYTE_SELECT_BITS];

wire oc_full;

// Coalescer - Groups subsequent requests together if they are compatible
// and the output register stage is stalled
bursting_coalescer #(
    .PAGE_ADDR_WIDTH(PAGE_SELECT_BITS),
    .TIMEOUT(16),
    .BURSTCOUNT_WIDTH(BURSTCOUNT_WIDTH),
    .MAXBURSTCOUNT(MAX_BURST),
    .MAX_THREADS(MAX_THREADS)
) coalescer (
    .clk(clk),
    .reset(reset),
    .i_page_addr(page_addr),
    .i_valid(i_valid && !oc_full && !w_fifo_full),
    .o_stall(c_stall),
    .o_new_page(c_new_page),
    .o_page_done(c_page_done),
    .o_req_addr(c_req_addr),
    .o_req_valid(c_req_valid),
    .i_stall(w_fifo_full),
    .o_burstcount(c_burstcount),
    .common_burst(common_burst)
);

// Writedata MUX
generate if( SEGMENT_SELECT_BITS > 0 )
begin
   wire [SEGMENT_SELECT_BITS-1:0] segment_select;
   assign segment_select = i_address[BYTE_SELECT_BITS-1:ALIGNMENT_ABITS];
   always@(*)
   begin
      wm_wide_wdata = {MWIDTH{1'bx}};
      wm_wide_wdata[segment_select*SEGMENT_WIDTH +: WIDTH] = i_writedata;

      wm_wide_be = {MWIDTH_BYTES{1'b0}};
      wm_wide_be[segment_select*SEGMENT_WIDTH_BYTES +: WIDTH_BYTES] = word_byte_enable;

      wm_wide_bite = {MWIDTH{1'b0}};
      wm_wide_bite[segment_select*SEGMENT_WIDTH +: WIDTH] = word_bit_enable;
   end
end
else
begin
   always@(*)
   begin
      wm_wide_wdata = {MWIDTH{1'bx}};
      wm_wide_wdata[0 +: WIDTH] = i_writedata;

      wm_wide_be = {MWIDTH_BYTES{1'b0}};
      wm_wide_be[0 +: WIDTH_BYTES] = word_byte_enable;

      wm_wide_bite = {MWIDTH{1'b0}};
      wm_wide_bite[0 +: WIDTH] = word_bit_enable;
   end
end
endgenerate

// Track the current write burst data - coalesce writes together until the 
// output registers are ready for a new request.
always@(posedge clk or posedge reset)
begin
    if(reset)
    begin
        wm_writedata <= {MWIDTH{1'b0}};
        wm_byteenable <= {MWIDTH_BYTES{1'b0}};
    end
    else
    begin
        if(c_new_page)
        begin
            wm_writedata <= wm_wide_wdata;
            wm_byteenable <= wm_wide_be;
        end
        else if(input_accepted)
        begin
            wm_writedata <= (wm_wide_wdata & wm_wide_bite) | (wm_writedata & ~wm_wide_bite);
            wm_byteenable <= wm_wide_be | wm_byteenable;
        end
    end
end

wire [COUNTER_WIDTH-1:0] num_threads_written;

// This FIFO stores the actual data to be written
//
//
wire w_data_fifo_full;
acl_data_fifo #(
    .DATA_WIDTH(COUNTER_WIDTH+MWIDTH+MWIDTH_BYTES),
    .DEPTH(WRITE_FIFO_DEPTH),
    .IMPL(HIGH_FMAX ? "ram_plus_reg" : "ram")
) req_fifo (
    .clock(clk),
    .resetn(~reset),
    .data_in( {wm_writedata,wm_byteenable} ),
    .valid_in( c_page_done & !w_fifo_full ),
    .data_out( {avm_writedata,avm_byteenable} ),
    .stall_in( ~write_accepted ),
    .stall_out( w_data_fifo_full )
);

// This FIFO stores the number of valid's to release with each writeack 
//
wire w_ack_fifo_full;
acl_data_fifo #(
    .DATA_WIDTH(COUNTER_WIDTH),
    .DEPTH(2*WRITE_FIFO_DEPTH),
    .IMPL(HIGH_FMAX ? "ram_plus_reg" : "ram")
) ack_fifo (
    .clock(clk),
    .resetn(~reset),
    .data_in( next_counter ),
    .valid_in( c_page_done & !w_fifo_full ),
    .data_out( num_threads_written ),
    .stall_in( ~avm_writeack ),
    .stall_out( w_ack_fifo_full )
);

// This FIFO hold the request information { address & burstcount }
//
wire w_fifo_stall_in;
assign w_fifo_stall_in = !(write_accepted && (burstcounter == avm_burstcount));

wire w_request_fifo_full;
acl_data_fifo #(
    .DATA_WIDTH(AWIDTH+BURSTCOUNT_WIDTH),
    .DEPTH(WRITE_FIFO_DEPTH),
    .IMPL(HIGH_FMAX ? "ram_plus_reg" : "ram")
) req_fifo2 (
    .clock(clk),
    .resetn(~reset),
    .data_in( {c_req_addr,{BYTE_SELECT_BITS{1'b0}},c_burstcount} ),
    .valid_in( c_req_valid & !w_fifo_full ), // The basical coalescer stalls on w_fifo_full holding c_req_valid high
    .data_out( {avm_address,avm_burstcount} ),
    .valid_out( avm_write ),
    .stall_in( w_fifo_stall_in ),
    .stall_out( w_request_fifo_full )
);

// The w_fifo_full is the OR of the data or request fifo's being full
assign w_fifo_full = w_data_fifo_full | w_request_fifo_full | w_ack_fifo_full;

// Occupancy counter - track the number of successfully transmitted writes
// and the number of writes pending in the next request.
//    occ_counter - the total occupancy (in threads) of the unit
//    next_counter - the number of threads coalesced into the next transfer
//    ack_counter - the number of pending threads with write completion acknowledged
assign input_accepted = i_valid && !o_stall;
assign write_accepted = avm_write && !avm_waitrequest;
assign output_acknowledged = o_valid && !i_stall;

always@(posedge clk or posedge reset)
begin
   if(reset == 1'b1)
   begin
      occ_counter <= {COUNTER_WIDTH{1'b0}};
      occ_counter_neg <= {COUNTER_WIDTH{1'b0}};
      ack_counter <= {COUNTER_WIDTH{1'b0}};
      next_counter <= {COUNTER_WIDTH{1'b0}};
      inflight_stores_counter <= {COUNTER_WIDTH{1'b0}};
      inflight_stores_counter_neg <= {COUNTER_WIDTH{1'b0}};
      burstcounter <= 6'b000001;
      o_active <= 1'b0;
   end
   else
   begin
      occ_counter <= occ_counter + input_accepted - output_acknowledged;
      occ_counter_neg <= occ_counter_neg - input_accepted + output_acknowledged;

      inflight_stores_counter <= inflight_stores_counter + input_accepted -
              (avm_writeack ? num_threads_written : {COUNTER_WIDTH{1'b0}});
      inflight_stores_counter_neg <= inflight_stores_counter_neg - input_accepted + 
              (avm_writeack ? num_threads_written : {COUNTER_WIDTH{1'b0}});

      next_counter <= input_accepted + (c_page_done ? {COUNTER_WIDTH{1'b0}} : next_counter);
      ack_counter <= ack_counter 
              - (USE_WRITE_ACK ? (avm_writeack ? num_threads_written : {COUNTER_WIDTH{1'b0}}) : input_accepted)
              + output_acknowledged;
      burstcounter <= write_accepted ? ((burstcounter == avm_burstcount) ? 6'b000001 : burstcounter+1) : burstcounter;
      o_active <= occ_counter_neg[COUNTER_WIDTH-1] | inflight_stores_counter_neg[COUNTER_WIDTH-1];
   end
end
assign oc_full = occ_counter[COUNTER_WIDTH-1] | inflight_stores_counter[COUNTER_WIDTH-1];

// Pipeline control signals
assign o_stall = oc_full || c_stall || w_fifo_full;
assign o_valid = ack_counter[COUNTER_WIDTH-1];

endmodule

// BURST COALESCING MODULE 
//
// Similar to the basic coalescer but supports checking if accesses are in consecutive DRAM "pages"
// Supports the ad-hocly discovered protocols for bursting efficiently with avalaon
// - Don't burst from an ODD address
// - If not on a burst boundary, then just burst up to the next burst bondary
//
// Yes, I know, this could be incorporated into the basic coalescer. But that's really not my "thing"
//
module bursting_coalescer
( 
    clk, reset, i_page_addr, i_valid, o_stall, o_new_page, o_page_done, o_req_addr, o_burstcount,
    o_req_valid, i_stall, 
    common_burst,

    // For the purposes of maintaining latency correctly, we need to know if total # of threads 
    // accepted by the caching LSU
    i_input_accepted_from_wrapper_lsu,
    i_reset_timeout
);

parameter PAGE_ADDR_WIDTH=32;
parameter TIMEOUT=8; 
parameter BURSTCOUNT_WIDTH=6; // Size of Avalon burst count port
parameter MAXBURSTCOUNT=32;   // This isn't the max supported by Avalon, but the max that the instantiating module needs
parameter MAX_THREADS=64;     // Must be a power of 2
parameter USECACHING=0;

localparam SLAVE_MAX_BURST=2**(BURSTCOUNT_WIDTH-1);
localparam THREAD_COUNTER_WIDTH=$clog2(MAX_THREADS+1);

input clk;
input reset;

input [PAGE_ADDR_WIDTH-1:0] i_page_addr;
input i_valid;
output o_stall;
output o_new_page;
output o_page_done;

output [PAGE_ADDR_WIDTH-1:0] o_req_addr;
output o_req_valid;
output [BURSTCOUNT_WIDTH-1:0] o_burstcount;
input i_stall;

input common_burst;

input i_input_accepted_from_wrapper_lsu;
input i_reset_timeout;

reg [PAGE_ADDR_WIDTH-1:0] page_addr;
reg [PAGE_ADDR_WIDTH-1:0] last_page_addr;
reg [PAGE_ADDR_WIDTH-1:0] last_page_addr_p1;
reg [BURSTCOUNT_WIDTH-1:0] burstcount;
reg valid;
wire ready;
wire waiting;
wire match;


wire timeout;
reg [$clog2(TIMEOUT):0] timeout_counter;

reg [THREAD_COUNTER_WIDTH-1:0] thread_counter;

generate if(USECACHING)
begin
   assign timeout = timeout_counter[$clog2(TIMEOUT)] | thread_counter[THREAD_COUNTER_WIDTH-1];
end
else
begin
   assign timeout = timeout_counter[$clog2(TIMEOUT)];
end
endgenerate

// Internal signal logic
wire match_burst_address;
wire match_next_page;
wire match_current_page;

generate
if ( BURSTCOUNT_WIDTH > 1 )
begin
  assign match_next_page     = (i_page_addr[BURSTCOUNT_WIDTH-2:0] == last_page_addr_p1[BURSTCOUNT_WIDTH-2:0]) && (|last_page_addr_p1[BURSTCOUNT_WIDTH-2:0]);
  assign match_current_page  = (i_page_addr[BURSTCOUNT_WIDTH-2:0] == last_page_addr[BURSTCOUNT_WIDTH-2:0]); 
end
else
begin
  assign match_next_page     = 1'b0;
  assign match_current_page  = 1'b1;
end
endgenerate

assign match_burst_address = common_burst;//(i_page_addr[PAGE_ADDR_WIDTH-1:BURSTCOUNT_WIDTH-1] == last_page_addr[PAGE_ADDR_WIDTH-1:BURSTCOUNT_WIDTH-1]);

assign match = (match_burst_address && (match_current_page || match_next_page)) && !thread_counter[THREAD_COUNTER_WIDTH-1];

assign ready = !valid || !(i_stall || waiting);
assign waiting = !timeout && (!i_valid || match);

wire input_accepted = i_valid && !o_stall;

always@(posedge clk or posedge reset)
begin
    if(reset)
    begin
        page_addr <= {PAGE_ADDR_WIDTH{1'b0}};
        last_page_addr <= {PAGE_ADDR_WIDTH{1'b0}};
        last_page_addr_p1 <= {PAGE_ADDR_WIDTH{1'b0}};
        burstcount <= 1;
        valid <= 1'b0;
        timeout_counter <= 0;
        thread_counter <= {THREAD_COUNTER_WIDTH{1'b0}};
    end
    else
    begin
        page_addr <= ready ? i_page_addr : page_addr;
        last_page_addr <= ready ? i_page_addr : (input_accepted && match_next_page ? i_page_addr : last_page_addr );
        last_page_addr_p1 <= ready ? i_page_addr+1 : (input_accepted && match_next_page ? i_page_addr+1 : last_page_addr_p1 );
        valid <= ready ? i_valid : valid;
        burstcount <= ready ? 6'b000001 : (input_accepted &&  match_next_page ? burstcount+1 : burstcount );
        thread_counter <= ready ? 1 : (USECACHING ? 
                                        (i_input_accepted_from_wrapper_lsu && !thread_counter[THREAD_COUNTER_WIDTH-1] ? thread_counter+1 : thread_counter ) : 
                                        (input_accepted ? thread_counter+1 : thread_counter));

        if( USECACHING && i_reset_timeout || !USECACHING && i_valid )
            timeout_counter <= 'd1;
        else if( valid && !timeout )
            timeout_counter <= timeout_counter + 'd1;
    end
end

// Outputs
assign o_stall = !match && !ready && i_valid;
// We're starting a new page (used by loads)
assign o_new_page = ready || i_valid && match_burst_address && !thread_counter[THREAD_COUNTER_WIDTH-1] && match_next_page;
assign o_req_addr = page_addr;
assign o_burstcount = burstcount;
assign o_req_valid = valid && !waiting;
// We're just finished with a page (used by stores)
assign o_page_done = valid && !waiting && !i_stall || !ready && i_valid && match_burst_address && !thread_counter[THREAD_COUNTER_WIDTH-1] && match_next_page; 
endmodule



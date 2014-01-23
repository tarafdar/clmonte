package require -exact qsys 12.1

set_validation_property AUTOMATIC_VALIDATION false

add_instance vector_add_kernel_system vector_add_kernel_system

add_connection acl_iface.kernel_clk vector_add_kernel_system.clock_reset
add_connection acl_iface.kernel_clk2x vector_add_kernel_system.clock_reset2x
add_connection acl_iface.kernel_reset vector_add_kernel_system.clock_reset_reset

add_connection vector_add_kernel_system.avm_memgmem0_port_0_0_rw acl_iface.kernel_mem0
add_connection vector_add_kernel_system.avm_memgmem0_port_1_0_rw acl_iface.kernel_mem1


add_connection acl_iface.kernel_irq vector_add_kernel_system.kernel_irq

add_connection acl_iface.kernel_cra vector_add_kernel_system.avs_vector_add_cra
set_connection_parameter_value acl_iface.kernel_cra/vector_add_kernel_system.avs_vector_add_cra baseAddress "0x0"


save_system

# (C) 1992-2014 Altera Corporation. All rights reserved.                         
# Your use of Altera Corporation's design tools, logic functions and other       
# software and tools, and its AMPP partner logic functions, and any output       
# files any of the foregoing (including device programming or simulation         
# files), and any associated documentation or information are expressly subject  
# to the terms and conditions of the Altera Program License Subscription         
# Agreement, Altera MegaCore Function License Agreement, or other applicable     
# license agreement, including, without limitation, that your use is for the     
# sole purpose of programming logic devices manufactured by Altera and sold by   
# Altera or its authorized distributors.  Please refer to the applicable         
# agreement for further details.                                                 
    

# Wrapper to call a script as if it was a function.  We can't make the script
# into a function because then we can't call it independently from command line.
# This wrapper basically just sets argc and argv appropriately, to avoid manually
# setting them before sourcing a child script.

proc call_script_as_function { script_name args } {
  post_message "Script_as_function wrapper: Executing script: $script_name"

  # Save original global argv/argc
  set orig_argv $::argv
  set orig_argc $::argc

  # Set global argv/argc
  set ::argv $args
  set ::argc [llength $args]

  # Execute script in new stack frame
  set rc [catch {uplevel [list source $script_name]} ret]

  # Restore original argv/argc
  set ::argv $orig_argv
  set ::argc $orig_argc
  return -code $rc $ret
}


